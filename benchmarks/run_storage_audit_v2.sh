#!/bin/bash
set -euo pipefail

FASTGA_DIR="/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA"
FASTGA="$FASTGA_DIR/FastGA"
EXAMPLE_DIR="$FASTGA_DIR/EXAMPLE"
HAP1="$EXAMPLE_DIR/HAP1.fasta.gz"
HAP2="$EXAMPLE_DIR/HAP2.fasta.gz"
BENCH_DIR="$FASTGA_DIR/benchmarks"
MONITOR="$BENCH_DIR/monitor_tmpdir.sh"
THREADS=(1 2 4 8 16 32)

export PATH="$FASTGA_DIR:$PATH"

AUDIT_DIR="$BENCH_DIR/storage_audit_v2"
mkdir -p "$AUDIT_DIR"

CSV="$AUDIT_DIR/summary.csv"
echo "threads,peak_workdir_mb,peak_tmpdir_mb,peak_total_mb,peak_tmpdir_bytes" > "$CSV"

echo "=== Storage Audit V2: Separate workdir and tmpdir ==="
echo "Strategy: GDB/GIX go to workdir, temp files go to dedicated tmpdir"
echo "du -sb on tmpdir captures deleted-but-open files"
echo ""

for T in "${THREADS[@]}"; do
    TDIR=$(printf "T%02d" "$T")
    WORKDIR="$AUDIT_DIR/$TDIR/work"
    TMPDIR_PATH="$AUDIT_DIR/$TDIR/tmp"
    mkdir -p "$WORKDIR" "$TMPDIR_PATH"
    cd "$WORKDIR"

    # Clean
    rm -rf "$WORKDIR"/* "$WORKDIR"/.* "$TMPDIR_PATH"/* 2>/dev/null || true

    cp "$HAP1" "$HAP2" .

    MONFILE="$AUDIT_DIR/${TDIR}_tmpdir_monitor.tsv"
    LOGFILE="$AUDIT_DIR/${TDIR}.log"

    echo "--- T=$T ---"

    # Run FastGA with -P pointing to separate tmpdir
    # Use -k to keep GIX files so we can measure workdir size
    "$FASTGA" -v -k -T"${T}" -P"$TMPDIR_PATH" HAP1.fasta.gz HAP2.fasta.gz \
        > /dev/null 2> "$LOGFILE" &
    FPID=$!

    # Monitor the tmpdir
    bash "$MONITOR" $FPID "$TMPDIR_PATH" "$MONFILE" &
    MPID=$!

    wait $FPID || true
    sleep 0.5
    kill $MPID 2>/dev/null || true
    wait $MPID 2>/dev/null || true

    # Measure workdir (persistent files)
    WORKDIR_BYTES=$(du -sb "$WORKDIR" 2>/dev/null | cut -f1)
    # Subtract input fasta.gz files
    INPUT_BYTES=0
    for f in "$WORKDIR"/HAP1.fasta.gz "$WORKDIR"/HAP2.fasta.gz; do
        [ -e "$f" ] && INPUT_BYTES=$((INPUT_BYTES + $(stat -c '%s' "$f")))
    done
    WORKDIR_BYTES=$((WORKDIR_BYTES - INPUT_BYTES))
    WORKDIR_MB=$(echo "scale=1; $WORKDIR_BYTES / 1048576" | bc)

    # Peak tmpdir from monitor
    PEAK_TMPDIR=$(awk -F'\t' 'NR>1 {if($2+0 > max) max=$2+0} END {print max+0}' "$MONFILE")
    PEAK_TMPDIR_MB=$(echo "scale=1; $PEAK_TMPDIR / 1048576" | bc)

    PEAK_TOTAL_MB=$(echo "scale=1; ($WORKDIR_BYTES + $PEAK_TMPDIR) / 1048576" | bc)

    echo "  Persistent files (workdir): ${WORKDIR_MB} MB"
    echo "  Peak temp files (tmpdir):   ${PEAK_TMPDIR_MB} MB"
    echo "  Peak total:                 ${PEAK_TOTAL_MB} MB"

    echo "$T,$WORKDIR_MB,$PEAK_TMPDIR_MB,$PEAK_TOTAL_MB,$PEAK_TMPDIR" >> "$CSV"

    # Show what's in tmpdir after run (should be empty if all temps were cleaned)
    REMAINING=$(ls -la "$TMPDIR_PATH" 2>/dev/null | wc -l)
    echo "  Files remaining in tmpdir:  $((REMAINING - 3))"

    # Clean up
    rm -rf "$WORKDIR"/* "$WORKDIR"/.* "$TMPDIR_PATH"/* 2>/dev/null || true
done

echo ""
echo "=== Summary ==="
column -t -s',' "$CSV"
echo ""
echo "Detailed monitor logs: $AUDIT_DIR/"
