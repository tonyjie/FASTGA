#!/bin/bash
# Three-way measurement harness: baseline (ddeea32, no -C) vs Opt C (agent-optimization,
# re-scans genome per chunk) vs fused-B (fused-scan-once, -R scan-once + tmpfs ktab).
#
# Usage: run_fused.sh <dataset> <T>
#   dataset = EXAMPLE | human   (EXAMPLE = docs deliverable for Task 4; human is Task 5)
#   T       = thread count (8 or 32 per the task brief)
#
# For each config emits stage_data/<dataset>/T<T>/<cfg>/{run.Llog,run.time,run.stderr,
# timeline.tsv,meta.txt,md5.txt,count.txt}.
#
# Storage model (see CLAUDE.md "Storage Problem" + progress ledger note):
#   GIXmake writes the .ktab + .gcnt sidecar into TPATH (the genome's own directory), NOT
#   into -P (SORT_PATH is only pos-list/seed-temp scratch). So the "real disk vs tmpfs" axis
#   is about *where the genome directory itself lives*, not just -P:
#     baseline / optC : genome dir + -P both on real disk (local /scratch)
#     fusedB          : genome dir AND -P both staged under /dev/shm -> ~zero real-disk bytes
#
# Sampling technique mirrors docs/agent_optimization/human_stages/run_stages_human.sh:
#   - `du -sb` on the genome/work directory for named (visible) files: GDB .bps, GIX .ktab,
#     per-chunk ktab sawtooth, output .1aln.
#   - unlinked-but-open temp (seed-pair _pair./_uniq./_algn. files, invisible to `du`) summed
#     from /proc/<pid>/fd `st_size` for FastGA + its direct children (system()-spawned
#     GIXmake keeps the same PID across execve, so `pgrep -P` catches it).
#   - RSS timeline sampled via `ps -o rss` summed over the same PID set (approximate; the
#     authoritative peak is /usr/bin/time -v "Maximum resident set size" in run.time).
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
DATASET="${1:?usage: run_fused.sh <EXAMPLE|human> <T>}"
T="${2:?usage: run_fused.sh <EXAMPLE|human> <T>}"

case "$DATASET" in
  EXAMPLE)
    G1="${G1:-$REPO/EXAMPLE/HAP1.fasta.gz}"
    G2="${G2:-$REPO/EXAMPLE/HAP2.fasta.gz}"
    ;;
  human)
    # Task 5 concern: ddeea32's FAtoGDB segfaults on CHM13 (see human_stages/); reusing this
    # harness for human requires substituting a working FAtoGDB into the baseline worktree,
    # exactly as human_stages/run_stages_human.sh does. Not wired here (out of Task 4 scope).
    G1="${G1:-/scratch/jl4257/seq_align/fastga_datasets/GRCh38/GCF_000001405.40_GRCh38.p14_genomic.fna}"
    G2="${G2:-/scratch/jl4257/seq_align/fastga_datasets/CHM13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna}"
    ;;
  *)
    echo "unknown dataset '$DATASET' (expected EXAMPLE or human)"; exit 1 ;;
esac

SCRATCH="${SCRATCH:-/scratch/jl4257/fused_stages/$DATASET}"   # real-disk work root
SD="$HERE/stage_data/$DATASET/T$T"
mkdir -p "$SD" "$SCRATCH"

WT_DDEEA32="${WT_DDEEA32:-$REPO/../FASTGA-fused-wt-ddeea32}"
WT_OPTC="${WT_OPTC:-$REPO/../FASTGA-fused-wt-optC}"

# For human, every config's FAtoGDB (all ddeea32-descended) segfaults building CHM13's GDB
# (ANO regression in the 11 commits after 5671357). Build a working 5671357 FAtoGDB once and
# shadow each config's binary with it (the GDB is binary-independent 2-bit sequence, so this
# only fixes the crash; it does not affect the GIX/merge/align under test).
WORKING_FATOGDB=""
if [ "$DATASET" = "human" ]; then
  WORKING_FATOGDB="$SCRATCH/FAtoGDB_working"
  if [ ! -x "$WORKING_FATOGDB" ]; then
    echo "[$(date +%H:%M:%S)] building working 5671357 FAtoGDB (ddeea32's segfaults on CHM13)"
    WT57="$REPO/../FASTGA-fused-wt-fatogdb57"
    git -C "$REPO" worktree remove --force "$WT57" 2>/dev/null; rm -rf "$WT57"
    git -C "$REPO" worktree add --detach "$WT57" 5671357 >/dev/null 2>&1
    ( cd "$WT57" && make FAtoGDB >/dev/null 2>&1 )
    cp "$WT57/FAtoGDB" "$WORKING_FATOGDB"
    git -C "$REPO" worktree remove --force "$WT57" 2>/dev/null; rm -rf "$WT57"
  fi
fi
install_fatogdb() { [ -n "$WORKING_FATOGDB" ] && cp "$WORKING_FATOGDB" "$1/FAtoGDB"; }

echo "[$(date +%H:%M:%S)] building fused-B (current repo checkout, $(git -C "$REPO" rev-parse --short HEAD))"
( cd "$REPO" && make FastGA GIXmake ONEview >/dev/null 2>&1 )
install_fatogdb "$REPO"

build_worktree() {
  wt=$1; ref=$2
  git -C "$REPO" worktree remove --force "$wt" 2>/dev/null; rm -rf "$wt"
  git -C "$REPO" worktree add --detach "$wt" "$ref" >/dev/null 2>&1
  ( cd "$wt" && make FastGA GIXmake ONEview >/dev/null 2>&1 )
  install_fatogdb "$wt"
}

# ---- one monitored run --------------------------------------------------------------------
# run_one <name> <bindir> <flags> <loc: real|tmpfs>
run_one() {
  name=$1; bindir=$2; flags=$3; loc=$4
  OUT="$SD/$name"; rm -rf "$OUT"; mkdir -p "$OUT"

  if [ "$loc" = "tmpfs" ]; then
    D="/dev/shm/fused_stages_${DATASET}_T${T}_${name}"
  else
    D="$SCRATCH/T${T}_${name}"
  fi
  rm -rf "$D"; mkdir -p "$D/work" "$D/tmp"
  B1=$(basename "$G1"); B2=$(basename "$G2")
  ln -sf "$G1" "$D/work/$B1"; ln -sf "$G2" "$D/work/$B2"

  echo "[$(date +%H:%M:%S)] $name (bin=$bindir flags='${flags:-none}' loc=$loc T=$T) start"
  echo "bindir=$bindir flags=$flags loc=$loc T=$T G1=$G1 G2=$G2" > "$OUT/meta.txt"

  MON="$OUT/timeline.tsv"
  echo -e "elapsed_s\trealdisk_mb\ttmpfs_mb\trss_mb" > "$MON"

  export PATH="$bindir:$PATH"
  cd "$D/work"
  /usr/bin/time -v -o "$OUT/run.time" \
    "$bindir/FastGA" -v -T"$T" $flags -P"$D/tmp" -L:"$OUT/run.Llog" -1:"$D/work/out" \
      "$B1" "$B2" > /dev/null 2> "$OUT/run.stderr" &
  FPID=$!; START=$(date +%s.%N)
  while kill -0 "$FPID" 2>/dev/null; do
    PIDS="$FPID $(pgrep -P "$FPID" 2>/dev/null||true)"
    if [ "$loc" = "tmpfs" ]; then
      REAL=0; TMPFS=$(du -sb "$D/work" 2>/dev/null|cut -f1)
    else
      REAL=$(du -sb "$D/work" 2>/dev/null|cut -f1); TMPFS=0
    fi
    Tfd=$(for p in $PIDS; do
            [ -d /proc/$p/fd ] || continue
            for fd in /proc/$p/fd/*; do [ -e "$fd" ] || continue
              case "$(readlink "$fd" 2>/dev/null||true)" in "$D/tmp"*) stat -L -c '%i %s' "$fd" 2>/dev/null;; esac
            done
          done | sort -u -k1,1 | awk '{s+=$2}END{print s+0}')
    if [ "$loc" = "tmpfs" ]; then TMPFS=$((${TMPFS:-0}+${Tfd:-0})); else REAL=$((${REAL:-0}+${Tfd:-0})); fi
    RSS=$(for p in $PIDS; do ps --no-headers -o rss -p "$p" 2>/dev/null; done | awk '{s+=$1}END{print s+0}')
    echo -e "$(echo "$(date +%s.%N)-$START"|bc)\t$(echo "scale=1;${REAL:-0}/1048576"|bc)\t$(echo "scale=1;${TMPFS:-0}/1048576"|bc)\t$(echo "scale=1;${RSS:-0}/1024"|bc)" >> "$MON"
    sleep 1
  done
  wait "$FPID"; RC=$?
  export PATH="${PATH#$bindir:}"

  if [ $RC -ne 0 ]; then
    echo "[$(date +%H:%M:%S)] $name FAILED (rc=$RC) -- see $OUT/run.stderr"
  else
    # payload md5: exclude '!' (provenance: invocation cmdline + timestamp) and '<' (input
    # file path references, e.g. "/dev/shm/mdchk 3") -- both vary with the run's working
    # directory/timestamp even when the alignment payload is byte-identical. Verified this
    # filter is stable across two independent same-content runs in different directories.
    MD5=$("$bindir/ONEview" "$D/work/out.1aln" 2>/dev/null | grep -v '^!' | grep -v '^<' | md5sum | cut -d' ' -f1)
    CNT=$("$bindir/ONEview" "$D/work/out.1aln" 2>/dev/null | grep -c '^A ')
    echo "$MD5" > "$OUT/md5.txt"; echo "$CNT" > "$OUT/count.txt"
    # leftover-scratch check: confirm the -R cleanup actually removed pos-lists/.gcnt
    LEFTOVER=$(find "$D" -maxdepth 2 \( -name '.post.*' -o -name '*.gcnt' \) 2>/dev/null | wc -l)
    echo "$LEFTOVER" > "$OUT/leftover_scratch.txt"
    echo "[$(date +%H:%M:%S)] $name done: md5=$MD5 count=$CNT leftover_scratch=$LEFTOVER"
  fi
  cd "$HERE"; rm -rf "$D"
}

# ---- configs --------------------------------------------------------------------------------
echo "[$(date +%H:%M:%S)] building baseline worktree (ddeea32)"
build_worktree "$WT_DDEEA32" ddeea32
run_one baseline "$WT_DDEEA32" "" real

echo "[$(date +%H:%M:%S)] building Opt C worktree (agent-optimization)"
build_worktree "$WT_OPTC" agent-optimization
run_one optC_C4 "$WT_OPTC" "-C4" real
run_one optC_C8 "$WT_OPTC" "-C8" real

run_one fusedB_C4 "$REPO" "-C4" tmpfs
run_one fusedB_C8 "$REPO" "-C8" tmpfs

git -C "$REPO" worktree remove --force "$WT_DDEEA32" 2>/dev/null; rm -rf "$WT_DDEEA32"
git -C "$REPO" worktree remove --force "$WT_OPTC" 2>/dev/null; rm -rf "$WT_OPTC"

# ---- gating summary -------------------------------------------------------------------------
echo
echo "=== gating checks (T=$T) ==="
REF_MD5=$(cat "$SD/baseline/md5.txt" 2>/dev/null || echo MISSING)
REF_CNT=$(cat "$SD/baseline/count.txt" 2>/dev/null || echo MISSING)
echo "baseline: md5=$REF_MD5 count=$REF_CNT"
for cfg in optC_C4 optC_C8 fusedB_C4 fusedB_C8; do
  m=$(cat "$SD/$cfg/md5.txt" 2>/dev/null || echo MISSING)
  c=$(cat "$SD/$cfg/count.txt" 2>/dev/null || echo MISSING)
  ok="OK"; [ "$m" = "$REF_MD5" ] && [ "$c" = "$REF_CNT" ] || ok="MISMATCH"
  echo "$cfg: md5=$m count=$c  [$ok]"
done
echo "ALL_DONE"
