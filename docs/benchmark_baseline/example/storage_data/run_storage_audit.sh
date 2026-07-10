#!/bin/bash
# Storage audit for STOCK UPSTREAM FastGA on the EXAMPLE dataset. Records the
# combined scratch footprint (persistent GDB/GIX + temp seed-pairs) over time
# per thread count, plus per-stage boundaries (FastGA -L). Writes to audit/:
#   summary.csv                 peak persistent / temp / total per thread
#   T{NN}_timeline.tsv          elapsed_s, persistent_mb, temp_mb  (work+tmp du)
#   T{NN}.Llog                  FastGA -L per-stage resource log (phase boundaries)
#
# Baseline binary must be stock-upstream (NOT optimize-memory / agent-optimization):
#   FASTGA=/path/to/upstream/FastGA bash run_storage_audit.sh
#
# IMPORTANT: the scratch dir MUST be on a filesystem where `du` still counts
# open()-then-unlink()ed files. NFS does (via silly-rename .nfsXXXX entries);
# local ext4/xfs and tmpfs do NOT (the directory entry is gone, so `du` on the
# dir skips those blocks) and would hide the seed-pair temp peak. The default
# ($HERE/.scratch, on the repo's NFS mount) works; override with BENCH_SCRATCH.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
FGA="${FASTGA:-$REPO/FastGA}"
export PATH="$(dirname "$FGA"):$PATH"          # FastGA shells out to FAtoGDB/GIXmake
MONITOR="$HERE/monitor_storage.sh"
HAP1="$REPO/EXAMPLE/HAP1.fasta.gz"; HAP2="$REPO/EXAMPLE/HAP2.fasta.gz"
SCRATCH="${BENCH_SCRATCH:-$HERE/.scratch}"     # disk-backed; override if $HOME is tmpfs/networked-slow
THREADS="${THREADS:-1 2 4 8 16 32}"

AUDIT="$HERE/audit"; mkdir -p "$AUDIT"; mkdir -p "$SCRATCH"
CSV="$AUDIT/summary.csv"
echo "threads,peak_persistent_mb,peak_temp_mb,peak_total_mb" > "$CSV"
echo "Binary : $FGA"
echo "Scratch: $SCRATCH  (must be disk-backed, not tmpfs)"

for T in $THREADS; do
  TDIR=$(printf "T%02d" "$T")
  WORK="$SCRATCH/work"; TMP="$SCRATCH/tmp"
  rm -rf "$WORK" "$TMP"; mkdir -p "$WORK" "$TMP"
  cp "$HAP1" "$HAP2" "$WORK/"; cd "$WORK"
  IN=$(( $(stat -c '%s' HAP1.fasta.gz) + $(stat -c '%s' HAP2.fasta.gz) ))
  MON="$AUDIT/${TDIR}_timeline.tsv"; Llog="$AUDIT/${TDIR}.Llog"; : > "$Llog"
  "$FGA" -v -k -T"$T" -P"$TMP" -L:"$Llog" HAP1.fasta.gz HAP2.fasta.gz > /dev/null 2> "$AUDIT/${TDIR}.log" &
  FPID=$!
  bash "$MONITOR" "$FPID" "$WORK" "$TMP" "$IN" "$MON" & MPID=$!
  wait $FPID || true; sleep 0.3; kill $MPID 2>/dev/null || true; wait $MPID 2>/dev/null || true

  read pk_p pk_t pk_tot < <(awk -F'\t' 'NR>1{if($2+0>p)p=$2; if($3+0>t)t=$3; if($2+$3>tot)tot=$2+$3}
                                        END{printf "%.1f %.1f %.1f",p,t,tot}' "$MON")
  echo "  $TDIR: persistent=${pk_p} MB  temp=${pk_t} MB  total=${pk_tot} MB"
  echo "$T,$pk_p,$pk_t,$pk_tot" >> "$CSV"
  cd "$HERE"; rm -rf "$WORK" "$TMP"
done
rm -rf "$SCRATCH"
echo "DONE -> $CSV + audit/T*_timeline.tsv + audit/T*.Llog"
