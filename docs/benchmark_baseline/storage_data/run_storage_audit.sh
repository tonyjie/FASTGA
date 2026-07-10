#!/bin/bash
# Storage audit for STOCK UPSTREAM FastGA on the EXAMPLE dataset. Writes
# audit/summary.csv + audit/T*_tmpdir_monitor.tsv next to this script.
#
# Baseline binary must be a stock-upstream build (NOT optimize-memory /
# agent-optimization). Point FASTGA at it — see ../README.md "Reproduce":
#   FASTGA=/path/to/upstream/FastGA bash run_storage_audit.sh
#
# Method: GDB/GIX go to a workdir; temp files go to a dedicated tmpdir; a
# background monitor polls `du -sb` on the tmpdir every 50 ms to capture the
# open()-then-unlink()ed temp files that are invisible to `ls`. Persistent
# files are measured with du after each run.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
FGA="${FASTGA:-$REPO/FastGA}"
export PATH="$(dirname "$FGA"):$PATH"      # FastGA shells out to GIXmake
MONITOR="$HERE/monitor_tmpdir.sh"
HAP1="$REPO/EXAMPLE/HAP1.fasta.gz"; HAP2="$REPO/EXAMPLE/HAP2.fasta.gz"
THREADS="${THREADS:-1 2 4 8 16 32}"

AUDIT="$HERE/audit"; mkdir -p "$AUDIT"
CSV="$AUDIT/summary.csv"
echo "threads,peak_workdir_mb,peak_tmpdir_mb,peak_total_mb,peak_tmpdir_bytes" > "$CSV"
echo "Binary: $FGA"

for T in $THREADS; do
  TDIR=$(printf "T%02d" "$T")
  WORK="$(mktemp -d)"; TMP="$(mktemp -d)"
  cp "$HAP1" "$HAP2" "$WORK/"; cd "$WORK"
  MON="$AUDIT/${TDIR}_tmpdir_monitor.tsv"
  # -k keeps GIX so we can measure the persistent footprint
  "$FGA" -v -k -T"$T" -P"$TMP" HAP1.fasta.gz HAP2.fasta.gz > /dev/null 2> "$AUDIT/${TDIR}.log" &
  FPID=$!
  bash "$MONITOR" "$FPID" "$TMP" "$MON" & MPID=$!
  wait $FPID || true; sleep 0.5; kill $MPID 2>/dev/null || true; wait $MPID 2>/dev/null || true

  WB=$(du -sb "$WORK" | cut -f1)
  IN=$(( $(stat -c '%s' "$WORK/HAP1.fasta.gz") + $(stat -c '%s' "$WORK/HAP2.fasta.gz") ))
  WB=$(( WB - IN ))
  PT=$(awk -F'\t' 'NR>1{if($2+0>m)m=$2+0}END{print m+0}' "$MON")
  printf "  %s: persistent=%.1f MB  temp=%.1f MB  total=%.1f MB\n" \
    "$TDIR" "$(bc<<<"scale=1;$WB/1048576")" "$(bc<<<"scale=1;$PT/1048576")" "$(bc<<<"scale=1;($WB+$PT)/1048576")"
  echo "$T,$(bc<<<"scale=1;$WB/1048576"),$(bc<<<"scale=1;$PT/1048576"),$(bc<<<"scale=1;($WB+$PT)/1048576"),$PT" >> "$CSV"
  cd "$HERE"; rm -rf "$WORK" "$TMP"
done
echo "DONE -> $CSV"
