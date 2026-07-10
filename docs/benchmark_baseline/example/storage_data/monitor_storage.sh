#!/bin/bash
# monitor_storage.sh <PID> <workdir> <tmpdir> <input_bytes> <outfile>
# Polls `du -sb` on BOTH the workdir (persistent GDB/GIX) and the tmpdir (temp
# seed-pair / alignment files) every 50 ms while PID is alive, and records the
# combined storage footprint over time. Columns: elapsed_s, persistent_mb, temp_mb.
#
# NOTE: workdir/tmpdir MUST be on a filesystem where `du` counts open-but-
# unlinked files. NFS does (silly-rename); local ext4/xfs and tmpfs do NOT, so
# the seed-pair temp (the bulk of the peak) would be invisible there.
PID=$1; WORK=$2; TMP=$3; IN=$4; OUT=$5
echo -e "elapsed_s\tpersistent_mb\ttemp_mb" > "$OUT"
START=$(date +%s.%N)
while kill -0 "$PID" 2>/dev/null; do
  NOW=$(date +%s.%N); EL=$(echo "$NOW - $START" | bc)
  W=$(du -sb "$WORK" 2>/dev/null | cut -f1); T=$(du -sb "$TMP" 2>/dev/null | cut -f1)
  P=$(( ${W:-0} - IN )); [ "$P" -lt 0 ] && P=0
  echo -e "${EL}\t$(echo "scale=2;$P/1048576"|bc)\t$(echo "scale=2;${T:-0}/1048576"|bc)" >> "$OUT"
  sleep 0.05
done
