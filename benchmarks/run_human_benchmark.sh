#!/bin/bash
set -euo pipefail

# ============================================================
# FastGA Human Genome Benchmark: GRCh38 vs CHM13, T=32
# ============================================================

FASTGA_DIR="/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA"
FASTGA="$FASTGA_DIR/FastGA"
MONITOR="$FASTGA_DIR/benchmarks/monitor_tmpdir.sh"

SCRATCH="/scratch/jl4257/seq_align/fastga_datasets"
CHM13_DIR="$SCRATCH/CHM13"
GRCH38_DIR="$SCRATCH/GRCh38"
RESULT_DIR="$SCRATCH/GRCh38_vs_CHM13"

CHM13_FA="$CHM13_DIR/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
GRCH38_FA="$GRCH38_DIR/GCF_000001405.40_GRCh38.p14_genomic.fna"

THREADS=32

export PATH="$FASTGA_DIR:$PATH"

echo "============================================================"
echo "FastGA Human Genome Benchmark: GRCh38 vs CHM13"
echo "============================================================"
echo "Date:    $(date)"
echo "Host:    $(hostname)"
echo "CPU:     $(lscpu | grep 'Model name' | sed 's/.*: *//')"
echo "Threads: $THREADS"
echo ""

# ---- Check inputs exist ----
for f in "$CHM13_FA" "$GRCH38_FA"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Input not found: $f"
        exit 1
    fi
done
echo "CHM13:  $CHM13_FA ($(du -h "$CHM13_FA" | cut -f1))"
echo "GRCh38: $GRCH38_FA ($(du -h "$GRCH38_FA" | cut -f1))"
echo ""

# ---- Check disk space ----
FREE_KB=$(df -k "$SCRATCH" | tail -1 | awk '{print $4}')
FREE_GB=$((FREE_KB / 1048576))
echo "Free disk space on scratch: ${FREE_GB} GB"
if [ "$FREE_GB" -lt 150 ]; then
    echo "ERROR: Less than 150 GB free. Need ~86 GB peak + safety margin."
    echo "Free up space or reduce dataset size."
    exit 1
fi
echo ""

# ---- Clean old cached files ----
echo "=== Cleaning old cached files ==="
for DIR in "$CHM13_DIR" "$GRCH38_DIR"; do
    echo "Cleaning $DIR ..."
    rm -f "$DIR"/*.1gdb "$DIR"/.*.bps "$DIR"/*.gix "$DIR"/.*.ktab.* \
          "$DIR"/.*.post.* "$DIR"/*.1ano "$DIR"/*.1aln 2>/dev/null || true
done

echo "Cleaning $RESULT_DIR ..."
rm -rf "$RESULT_DIR"/thread_* "$RESULT_DIR"/_tmp_sort \
       "$RESULT_DIR"/tmp "$RESULT_DIR"/logs "$RESULT_DIR"/timing \
       "$RESULT_DIR"/monitor "$RESULT_DIR"/*.1aln 2>/dev/null || true
echo "Done cleaning."
echo ""

# ---- Create directory structure ----
TMPDIR_PATH="$RESULT_DIR/tmp"
LOGDIR="$RESULT_DIR/logs"
TIMINGDIR="$RESULT_DIR/timing"
MONDIR="$RESULT_DIR/monitor"
mkdir -p "$TMPDIR_PATH" "$LOGDIR" "$TIMINGDIR" "$MONDIR"

# ---- Record disk space before ----
FREE_BEFORE=$(df -k "$SCRATCH" | tail -1 | awk '{print $4}')
echo "Disk free before run: $((FREE_BEFORE / 1048576)) GB"
echo ""

# ---- Start monitor ----
MONFILE="$MONDIR/human_T${THREADS}_tmpdir.tsv"
echo "=== Starting tmpdir monitor ==="

# ---- Run FastGA ----
echo "=== Running FastGA -T${THREADS} ==="
echo "Start time: $(date)"
echo ""

/usr/bin/time -v -o "$TIMINGDIR/human_T${THREADS}.time" \
    "$FASTGA" -v -k -T"${THREADS}" \
    -P"$TMPDIR_PATH" \
    -L:"$LOGDIR/fastga_internal.log" \
    -1:"$RESULT_DIR/output.1aln" \
    "$GRCH38_FA" \
    "$CHM13_FA" \
    2> "$LOGDIR/verbose.log" &
FPID=$!

# Start tmpdir monitor
bash "$MONITOR" $FPID "$TMPDIR_PATH" "$MONFILE" &
MPID=$!

echo "FastGA PID: $FPID"
echo "Monitor PID: $MPID"
echo "Waiting for completion..."
echo ""

# Also periodically report disk usage to stderr
while kill -0 "$FPID" 2>/dev/null; do
    FREE_NOW=$(df -k "$SCRATCH" | tail -1 | awk '{print $4}')
    FREE_NOW_GB=$((FREE_NOW / 1048576))
    TMPDU=$(du -sb "$TMPDIR_PATH" 2>/dev/null | cut -f1)
    TMPDU_MB=$((TMPDU / 1048576))
    echo "[$(date +%H:%M:%S)] Disk free: ${FREE_NOW_GB} GB | Tmpdir: ${TMPDU_MB} MB"

    # Abort if disk gets too low
    if [ "$FREE_NOW_GB" -lt 50 ]; then
        echo "CRITICAL: Disk free dropped below 50 GB! Killing FastGA."
        kill "$FPID" 2>/dev/null || true
        kill "$MPID" 2>/dev/null || true
        exit 1
    fi

    sleep 30
done

# Wait for FastGA to finish
wait $FPID
FASTGA_EXIT=$?

sleep 1
kill $MPID 2>/dev/null || true
wait $MPID 2>/dev/null || true

echo ""
echo "FastGA exit code: $FASTGA_EXIT"
echo "End time: $(date)"
echo ""

# ---- Record disk space after ----
FREE_AFTER=$(df -k "$SCRATCH" | tail -1 | awk '{print $4}')
echo "Disk free after run: $((FREE_AFTER / 1048576)) GB"
echo "Disk consumed by run: $(( (FREE_BEFORE - FREE_AFTER) / 1048576 )) GB"
echo ""

# ---- Measure persistent file sizes ----
echo "=== Persistent File Sizes ==="
SIZES_FILE="$RESULT_DIR/file_sizes_report.txt"
echo "=== Persistent File Sizes ===" > "$SIZES_FILE"
echo "Date: $(date)" >> "$SIZES_FILE"
echo "" >> "$SIZES_FILE"

for GENOME_DIR in "$GRCH38_DIR" "$GRCH38_DIR/../CHM13"; do
    GENOME_DIR=$(realpath "$GENOME_DIR")
    GENOME_NAME=$(basename "$GENOME_DIR")
    echo "--- $GENOME_NAME ---" >> "$SIZES_FILE"

    GDB_SIZE=0; BPS_SIZE=0; GIX_SIZE=0; KTAB_TOTAL=0

    for f in "$GENOME_DIR"/*.1gdb; do
        [ -e "$f" ] && GDB_SIZE=$(stat -c '%s' "$f")
    done
    for f in "$GENOME_DIR"/.*.bps; do
        [ -e "$f" ] && BPS_SIZE=$((BPS_SIZE + $(stat -c '%s' "$f")))
    done
    for f in "$GENOME_DIR"/*.gix; do
        [ -e "$f" ] && GIX_SIZE=$(stat -c '%s' "$f")
    done
    for f in "$GENOME_DIR"/.*.ktab.*; do
        [ -e "$f" ] && KTAB_TOTAL=$((KTAB_TOTAL + $(stat -c '%s' "$f")))
    done

    GIX_ALL=$((GIX_SIZE + KTAB_TOTAL))
    TOTAL=$((GDB_SIZE + BPS_SIZE + GIX_ALL))

    GIX_ALL_GB=$(echo "scale=2; $GIX_ALL / 1073741824" | bc)
    TOTAL_GB=$(echo "scale=2; $TOTAL / 1073741824" | bc)
    BPS_MB=$(echo "scale=1; $BPS_SIZE / 1048576" | bc)

    printf "  GDB (.1gdb):     %15s bytes\n" "$GDB_SIZE" >> "$SIZES_FILE"
    printf "  BPS (.bps):      %15s bytes  (%s MB)\n" "$BPS_SIZE" "$BPS_MB" >> "$SIZES_FILE"
    printf "  GIX stub (.gix): %15s bytes\n" "$GIX_SIZE" >> "$SIZES_FILE"
    printf "  KTAB partitions: %15s bytes\n" "$KTAB_TOTAL" >> "$SIZES_FILE"
    printf "  GIX total:       %15s bytes  (%s GB)\n" "$GIX_ALL" "$GIX_ALL_GB" >> "$SIZES_FILE"
    printf "  TOTAL:           %15s bytes  (%s GB)\n" "$TOTAL" "$TOTAL_GB" >> "$SIZES_FILE"
    echo "" >> "$SIZES_FILE"

    echo "$GENOME_NAME: GIX=$GIX_ALL_GB GB, BPS=$BPS_MB MB, Total=$TOTAL_GB GB"
done

# Output .1aln size
if [ -f "$RESULT_DIR/output.1aln" ]; then
    ALN_SIZE=$(stat -c '%s' "$RESULT_DIR/output.1aln")
    ALN_MB=$(echo "scale=1; $ALN_SIZE / 1048576" | bc)
    echo "output.1aln: ${ALN_MB} MB" >> "$SIZES_FILE"
    echo "Output .1aln: ${ALN_MB} MB"
fi

cat "$SIZES_FILE"
echo ""

# ---- Parse timing ----
echo "=== Timing Results ==="
TIMEFILE="$TIMINGDIR/human_T${THREADS}.time"
if [ -f "$TIMEFILE" ]; then
    WALL=$(grep "Elapsed (wall clock)" "$TIMEFILE" | sed 's/.*): //')
    USER_T=$(grep "User time" "$TIMEFILE" | sed 's/.*): //')
    SYS_T=$(grep "System time" "$TIMEFILE" | sed 's/.*): //')
    RSS=$(grep "Maximum resident" "$TIMEFILE" | sed 's/.*): //')
    CPUPCT=$(grep "Percent of CPU" "$TIMEFILE" | sed 's/.*): //' | tr -d '%')

    echo "Wall clock:    $WALL"
    echo "User CPU:      ${USER_T}s"
    echo "System CPU:    ${SYS_T}s"
    echo "Peak RSS:      ${RSS} KB ($(echo "scale=1; $RSS / 1048576" | bc) GB)"
    echo "CPU%:          ${CPUPCT}%"
fi
echo ""

# ---- Parse peak temp from monitor ----
echo "=== Peak Temp Storage ==="
if [ -f "$MONFILE" ]; then
    PEAK_TMP=$(awk -F'\t' 'NR>1 {if($2+0 > max) max=$2+0} END {print max+0}' "$MONFILE")
    PEAK_TMP_GB=$(echo "scale=2; $PEAK_TMP / 1073741824" | bc)
    echo "Peak temp dir usage: ${PEAK_TMP_GB} GB (${PEAK_TMP} bytes)"
fi
echo ""

# ---- Show verbose log tail (per-phase timing) ----
echo "=== Per-Phase Timing (from verbose log) ==="
if [ -f "$LOGDIR/verbose.log" ]; then
    grep -E "(Resources for phase|Total Resources|Creating|Starting|Sorting|Converting|Total seeds|Total hits)" \
        "$LOGDIR/verbose.log" | head -30
fi
echo ""

echo "============================================================"
echo "Human genome benchmark complete!"
echo "============================================================"
echo "Results dir:   $RESULT_DIR"
echo "File sizes:    $SIZES_FILE"
echo "Timing:        $TIMEFILE"
echo "Verbose log:   $LOGDIR/verbose.log"
echo "Monitor data:  $MONFILE"
echo "Output .1aln:  $RESULT_DIR/output.1aln"
