#!/bin/bash
# Incremental optimization profiling on the HUMAN pair (GRCh38 x CHM13, T=32).
# For each cumulative optimization stage it builds that stage's binary, then runs
# one monitored FastGA -T32 that captures BOTH:
#   - storage footprint over time  -> stage_data/<stage>/timeline.tsv
#   - per-stage runtime + CPU%      -> stage_data/<stage>/run.Llog  (FastGA -L)
#   - end-to-end                    -> stage_data/<stage>/run.time
#
# Stages (cumulative, on the agent-optimization branch):
#   baseline  ddeea32           stock upstream
#   opt1      50e4a16           + early GIX deletion
#   opt3      6f90a69           + drop mask byte
#   opt4      a946173           + drop LCP byte (on-the-fly)
#   optC_C16  1b6014e  -C16     + bilateral chunked build/merge
#
# GDB stage uses a working 5671357 FAtoGDB (ddeea32's segfaults on CHM13); it is
# identical across stages. Storage runs on local /scratch, so temp (unlinked-open,
# invisible to du on local FS) is summed from /proc/PID/fd via st_size.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
G1="${G1:-/scratch/jl4257/seq_align/fastga_datasets/GRCh38/GCF_000001405.40_GRCh38.p14_genomic.fna}"
G2="${G2:-/scratch/jl4257/seq_align/fastga_datasets/CHM13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna}"
SCRATCH="${SCRATCH:-/scratch/jl4257/human_stages}"
WT="${WT:-$REPO/../FASTGA-stage-wt}"
T=32
FATOGDB="$SCRATCH/FAtoGDB_working"

# --- build a working 5671357 FAtoGDB once ---
mkdir -p "$SCRATCH"
FG57="$REPO/../FASTGA-fatogdb57"
git -C "$REPO" worktree remove --force "$FG57" 2>/dev/null; rm -rf "$FG57"
git -C "$REPO" worktree add --detach "$FG57" 5671357 >/dev/null 2>&1
( cd "$FG57" && make FAtoGDB >/dev/null 2>&1 )
cp "$FG57/FAtoGDB" "$FATOGDB"
git -C "$REPO" worktree remove --force "$FG57" 2>/dev/null

run_stage() {
  name=$1; commit=$2; flags=$3
  # build this stage's binary
  git -C "$REPO" worktree remove --force "$WT" 2>/dev/null; rm -rf "$WT"
  git -C "$REPO" worktree add --detach "$WT" "$commit" >/dev/null 2>&1
  ( cd "$WT" && make FastGA GIXmake >/dev/null 2>&1 )
  cp "$FATOGDB" "$WT/FAtoGDB"                      # working GDB builder
  export PATH="$WT:$PATH"

  D="$SCRATCH/run"; rm -rf "$D"; mkdir -p "$D/work" "$D/tmp"
  ln -sf "$G1" "$D/work/g1.fna"; ln -sf "$G2" "$D/work/g2.fna"
  OUT="$HERE/stage_data/$name"; mkdir -p "$OUT"
  MON="$OUT/timeline.tsv"; echo -e "elapsed_s\tpersistent_mb\ttemp_mb" > "$MON"
  echo "[$(date +%H:%M:%S)] $name ($commit ${flags:-none}) start"
  cd "$D/work"
  /usr/bin/time -v -o "$OUT/run.time" \
    "$WT/FastGA" -v -T$T $flags -P"$D/tmp" -L:"$OUT/run.Llog" g1.fna g2.fna > /dev/null 2> "$OUT/run.stderr" &
  FPID=$!; START=$(date +%s.%N)
  while kill -0 "$FPID" 2>/dev/null; do
    W=$(du -sb "$D/work" 2>/dev/null|cut -f1)
    Tfd=$(for p in "$FPID" $(pgrep -P "$FPID" 2>/dev/null||true); do
            [ -d /proc/$p/fd ] || continue
            for fd in /proc/$p/fd/*; do [ -e "$fd" ] || continue
              case "$(readlink "$fd" 2>/dev/null||true)" in "$D/tmp"*) stat -L -c '%i %s' "$fd" 2>/dev/null;; esac
            done
          done | sort -u -k1,1 | awk '{s+=$2}END{print s+0}')
    echo -e "$(echo "$(date +%s.%N)-$START"|bc)\t$(echo "scale=1;${W:-0}/1048576"|bc)\t$(echo "scale=1;${Tfd:-0}/1048576"|bc)" >> "$MON"
    sleep 2
  done
  wait "$FPID"
  export PATH="${PATH#$WT:}"
  echo "[$(date +%H:%M:%S)] $name done"
  cd "$HERE"; rm -rf "$D"
}

run_stage baseline ddeea32 ""
run_stage opt1     50e4a16 ""
run_stage opt3     6f90a69 ""
run_stage opt4     a946173 ""
run_stage optC_C16 1b6014e "-C16"
git -C "$REPO" worktree remove --force "$WT" 2>/dev/null; rm -rf "$SCRATCH"
echo "ALL_DONE"
