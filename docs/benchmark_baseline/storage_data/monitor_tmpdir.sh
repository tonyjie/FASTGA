#!/bin/bash
# monitor_tmpdir.sh <PID> <tmpdir> <output_file>
# Monitors actual disk blocks consumed in a directory using du.
# This catches deleted-but-open files because the filesystem still
# accounts for their blocks until the FDs are closed.
# We use a dedicated empty tmpdir so du only counts FastGA's temp files.

PID=$1
TMPDIR_PATH=$2
OUTFILE=$3
INTERVAL=0.05

if [ -z "$PID" ] || [ -z "$TMPDIR_PATH" ] || [ -z "$OUTFILE" ]; then
    echo "Usage: $0 <PID> <tmpdir> <output_file>"
    exit 1
fi

echo -e "elapsed_s\tdu_bytes\tls_bytes\tls_files" > "$OUTFILE"
START=$(date +%s.%N)

while kill -0 "$PID" 2>/dev/null; do
    NOW=$(date +%s.%N)
    ELAPSED=$(echo "$NOW - $START" | bc)

    # du -sb counts actual disk blocks used (includes deleted-but-open files on some FS)
    DU_BYTES=$(du -sb "$TMPDIR_PATH" 2>/dev/null | cut -f1)

    # ls visible files
    LS_BYTES=0
    LS_COUNT=0
    for f in "$TMPDIR_PATH"/*; do
        if [ -e "$f" ] && [ -f "$f" ]; then
            SZ=$(stat -c '%s' "$f" 2>/dev/null) || continue
            LS_BYTES=$((LS_BYTES + SZ))
            LS_COUNT=$((LS_COUNT + 1))
        fi
    done

    echo -e "${ELAPSED}\t${DU_BYTES}\t${LS_BYTES}\t${LS_COUNT}" >> "$OUTFILE"
    sleep $INTERVAL
done

# Final measurement after process exits
NOW=$(date +%s.%N)
ELAPSED=$(echo "$NOW - $START" | bc)
DU_BYTES=$(du -sb "$TMPDIR_PATH" 2>/dev/null | cut -f1)
echo -e "${ELAPSED}\t${DU_BYTES}\t0\t0" >> "$OUTFILE"
