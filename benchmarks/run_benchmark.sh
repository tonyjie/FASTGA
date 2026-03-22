#!/bin/bash
set -euo pipefail

# ============================================================
# FastGA Benchmark Suite
# Thread scaling study + intermediate file size measurement
# ============================================================

FASTGA_DIR="/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA"
FASTGA="$FASTGA_DIR/FastGA"
BENCH_DIR="$FASTGA_DIR/benchmarks"
EXAMPLE_DIR="$FASTGA_DIR/EXAMPLE"
HAP1="$EXAMPLE_DIR/HAP1.fasta.gz"
HAP2="$EXAMPLE_DIR/HAP2.fasta.gz"
MONITOR="$BENCH_DIR/monitor_fd_sizes.sh"
THREADS=(1 2 4 8 16 32 64)
REPEATS=3

# Ensure PATH includes FastGA binaries
export PATH="$FASTGA_DIR:$PATH"

# ---- Setup ----
mkdir -p "$BENCH_DIR"/{results/{T01,T02,T04,T08,T16,T32,T64,file_sizes},logs,timing,fd_monitor}

# Verify binaries
for bin in FastGA FAtoGDB GIXmake; do
    if [ ! -x "$FASTGA_DIR/$bin" ]; then
        echo "ERROR: $FASTGA_DIR/$bin not found or not executable"
        exit 1
    fi
done

echo "============================================================"
echo "FastGA Benchmark Suite"
echo "============================================================"
echo "Date:    $(date)"
echo "Host:    $(hostname)"
echo "CPU:     $(lscpu | grep 'Model name' | sed 's/.*: *//')"
echo "Cores:   $(nproc)"
echo "Threads: ${THREADS[*]}"
echo "Repeats: $REPEATS"
echo ""

# ============================================================
# Phase 1: Persistent file sizes
# ============================================================
echo "=== Phase 1: Measuring persistent file sizes ==="
FSIZE_DIR="$BENCH_DIR/results/file_sizes"
cd "$FSIZE_DIR"

# Clean any previous run
rm -f *.1gdb .*.bps *.gix .*.ktab.* .*.post.* *.1aln HAP1.fasta.gz HAP2.fasta.gz 2>/dev/null || true

# Copy inputs (not symlink — FastGA writes outputs alongside the input path)
cp "$HAP1" "$HAP2" .

echo "Running FastGA -v -k -T8 to measure persistent files..."
/usr/bin/time -v -o "$BENCH_DIR/timing/file_sizes.time" \
    "$FASTGA" -v -k -T8 -P. HAP1.fasta.gz HAP2.fasta.gz 2>&1 | tee "$BENCH_DIR/logs/file_sizes.log"

# Record persistent file sizes
REPORT="$BENCH_DIR/results/file_sizes_report.txt"
echo "=== Persistent File Sizes ===" > "$REPORT"
echo "Date: $(date)" >> "$REPORT"
echo "" >> "$REPORT"

echo "--- Individual Files ---" >> "$REPORT"
for f in HAP1.1gdb .HAP1.bps HAP1.gix .HAP1.ktab.* .HAP1.post.* \
         HAP2.1gdb .HAP2.bps HAP2.gix .HAP2.ktab.* .HAP2.post.*; do
    if [ -e "$f" ]; then
        SIZE=$(stat -c '%s' "$f")
        SIZE_MB=$(echo "scale=2; $SIZE / 1048576" | bc)
        printf "%-30s %15s bytes  (%s MB)\n" "$f" "$SIZE" "$SIZE_MB" >> "$REPORT"
    fi
done

# Check for output .1aln
for f in *.1aln; do
    if [ -e "$f" ]; then
        SIZE=$(stat -c '%s' "$f")
        SIZE_MB=$(echo "scale=2; $SIZE / 1048576" | bc)
        printf "%-30s %15s bytes  (%s MB)\n" "$f" "$SIZE" "$SIZE_MB" >> "$REPORT"
    fi
done

echo "" >> "$REPORT"
echo "--- Summary by Category ---" >> "$REPORT"
for genome in HAP1 HAP2; do
    GDB_SIZE=$(stat -c '%s' "${genome}.1gdb" 2>/dev/null || echo 0)
    BPS_SIZE=$(stat -c '%s' ".${genome}.bps" 2>/dev/null || echo 0)
    GIX_SIZE=$(stat -c '%s' "${genome}.gix" 2>/dev/null || echo 0)

    KTAB_TOTAL=0
    for f in ".${genome}.ktab."*; do
        [ -e "$f" ] && KTAB_TOTAL=$((KTAB_TOTAL + $(stat -c '%s' "$f")))
    done

    POST_TOTAL=0
    for f in ".${genome}.post."*; do
        [ -e "$f" ] && POST_TOTAL=$((POST_TOTAL + $(stat -c '%s' "$f")))
    done

    TOTAL=$((GDB_SIZE + BPS_SIZE + GIX_SIZE + KTAB_TOTAL + POST_TOTAL))
    TOTAL_MB=$(echo "scale=2; $TOTAL / 1048576" | bc)
    GIX_ALL=$((GIX_SIZE + KTAB_TOTAL + POST_TOTAL))
    GIX_ALL_MB=$(echo "scale=2; $GIX_ALL / 1048576" | bc)

    echo "" >> "$REPORT"
    echo "$genome:" >> "$REPORT"
    printf "  GDB (.1gdb):     %12s bytes\n" "$GDB_SIZE" >> "$REPORT"
    printf "  BPS (.bps):      %12s bytes\n" "$BPS_SIZE" >> "$REPORT"
    printf "  GIX stub (.gix): %12s bytes\n" "$GIX_SIZE" >> "$REPORT"
    printf "  KTAB partitions: %12s bytes\n" "$KTAB_TOTAL" >> "$REPORT"
    printf "  POST partitions: %12s bytes\n" "$POST_TOTAL" >> "$REPORT"
    printf "  GIX total:       %12s bytes  (%s MB)\n" "$GIX_ALL" "$GIX_ALL_MB" >> "$REPORT"
    printf "  TOTAL:           %12s bytes  (%s MB)\n" "$TOTAL" "$TOTAL_MB" >> "$REPORT"
done

echo ""
echo "Persistent file sizes recorded in: $REPORT"
cat "$REPORT"
echo ""

# ============================================================
# Phase 2: Thread scaling
# ============================================================
echo "=== Phase 2: Thread scaling study ==="

CSV="$BENCH_DIR/results/benchmark_summary.csv"
echo "threads,repeat,wall_clock_s,user_time_s,system_time_s,max_rss_kb,peak_fd_bytes,cpu_pct" > "$CSV"

for T in "${THREADS[@]}"; do
    TDIR=$(printf "T%02d" "$T")

    for REP in $(seq 1 $REPEATS); do
        echo "--- T=$T, repeat=$REP ---"

        WORKDIR="$BENCH_DIR/results/$TDIR/rep${REP}"
        mkdir -p "$WORKDIR"
        cd "$WORKDIR"

        # Clean previous run artifacts
        rm -f *.1gdb .*.bps *.gix .*.ktab.* .*.post.* *.1aln HAP1.fasta.gz HAP2.fasta.gz 2>/dev/null || true

        # Copy inputs (safe approach — avoids symlink path resolution issues)
        cp "$HAP1" "$HAP2" .

        TIMEFILE="$BENCH_DIR/timing/${TDIR}_rep${REP}.time"
        LOGFILE="$BENCH_DIR/logs/${TDIR}_rep${REP}.log"
        FDFILE="$BENCH_DIR/fd_monitor/${TDIR}_rep${REP}.tsv"
        FASTGA_LOG="$BENCH_DIR/logs/${TDIR}_rep${REP}.fastga_log"

        # Launch FastGA in background so we can capture its PID
        /usr/bin/time -v -o "$TIMEFILE" \
            "$FASTGA" -v -T"${T}" -P. -L:"$FASTGA_LOG" HAP1.fasta.gz HAP2.fasta.gz \
            > /dev/null 2> "$LOGFILE" &
        FASTGA_PID=$!

        # Start FD monitor in background
        bash "$MONITOR" $FASTGA_PID "$FDFILE" &
        MONITOR_PID=$!

        # Wait for FastGA to finish
        wait $FASTGA_PID || true

        # Give monitor time to detect process exit
        sleep 1
        kill $MONITOR_PID 2>/dev/null || true
        wait $MONITOR_PID 2>/dev/null || true

        # Parse /usr/bin/time output
        WALL=$(grep "Elapsed (wall clock)" "$TIMEFILE" | sed 's/.*): //')
        USER_T=$(grep "User time" "$TIMEFILE" | sed 's/.*): //')
        SYS_T=$(grep "System time" "$TIMEFILE" | sed 's/.*): //')
        RSS=$(grep "Maximum resident" "$TIMEFILE" | sed 's/.*): //')
        CPUPCT=$(grep "Percent of CPU" "$TIMEFILE" | sed 's/.*): //' | tr -d '%')

        # Convert wall clock (h:mm:ss or m:ss.ss) to seconds
        WALL_S=$(echo "$WALL" | awk -F: '{
            if (NF==3) printf "%.2f", $1*3600+$2*60+$3;
            else if (NF==2) printf "%.2f", $1*60+$2;
            else printf "%.2f", $1
        }')

        # Parse peak FD size from monitor
        PEAK_FD=$(awk -F'\t' 'NR>1 {if($2+0 > max) max=$2+0} END {print max+0}' "$FDFILE" 2>/dev/null || echo 0)

        echo "$T,$REP,$WALL_S,$USER_T,$SYS_T,$RSS,$PEAK_FD,$CPUPCT" >> "$CSV"

        # Clean up generated files
        rm -f *.1gdb .*.bps *.gix .*.ktab.* .*.post.* *.1aln HAP1.fasta.gz HAP2.fasta.gz 2>/dev/null || true

        echo "  Wall: ${WALL_S}s  User: ${USER_T}s  Sys: ${SYS_T}s  RSS: ${RSS}KB  Peak FD: ${PEAK_FD}  CPU: ${CPUPCT}%"
    done
done

echo ""
echo "=== Phase 2 complete ==="
echo "Raw data: $CSV"

# ============================================================
# Phase 3: Generate report
# ============================================================
echo ""
echo "=== Phase 3: Generating report ==="

REPORT_MD="$FASTGA_DIR/docs/benchmark_results.md"
mkdir -p "$FASTGA_DIR/docs"

cat > "$REPORT_MD" << 'HEADER'
# FastGA Benchmark Results: Thread Scaling & Intermediate File Sizes

## Test Configuration

| Parameter | Value |
|---|---|
HEADER

echo "| Dataset | EXAMPLE/HAP1.fasta.gz vs EXAMPLE/HAP2.fasta.gz (~86 Mbp each) |" >> "$REPORT_MD"
echo "| Server | $(hostname) |" >> "$REPORT_MD"
echo "| CPU | $(lscpu | grep 'Model name' | sed 's/.*: *//') |" >> "$REPORT_MD"
echo "| Total threads | $(nproc) |" >> "$REPORT_MD"
echo "| Thread counts tested | 1, 2, 4, 8, 16, 32, 64 |" >> "$REPORT_MD"
echo "| Repeats per config | 3 |" >> "$REPORT_MD"
echo "| Date | $(date) |" >> "$REPORT_MD"
echo "| Measurement tools | /usr/bin/time -v, FastGA -L log, /proc/PID/fd monitor |" >> "$REPORT_MD"

echo "" >> "$REPORT_MD"
echo "## Persistent File Sizes (with -k flag, T=8)" >> "$REPORT_MD"
echo "" >> "$REPORT_MD"
echo '```' >> "$REPORT_MD"
cat "$BENCH_DIR/results/file_sizes_report.txt" >> "$REPORT_MD"
echo '```' >> "$REPORT_MD"

echo "" >> "$REPORT_MD"
echo "## Thread Scaling Results (averaged over $REPEATS repeats)" >> "$REPORT_MD"
echo "" >> "$REPORT_MD"
echo '| Threads | Wall Clock (s) | User CPU (s) | Sys CPU (s) | Peak RSS (MB) | Peak Temp FD (MB) | CPU% | Speedup | Efficiency |' >> "$REPORT_MD"
echo '|--------:|---------------:|-------------:|------------:|--------------:|------------------:|-----:|--------:|-----------:|' >> "$REPORT_MD"

# Compute averages and speedup from CSV
awk -F, '
NR>1 {
    t=$1
    wall[t]+=$3; user[t]+=$4; sys[t]+=$5; rss[t]+=$6; fd[t]+=$7; cpu[t]+=$8; n[t]++
}
END {
    # Sort by thread count
    split("1 2 4 8 16 32 64", order, " ")
    base = wall[1]/n[1]
    for (i in order) {
        t = order[i]
        if (n[t] > 0) {
            w = wall[t]/n[t]
            u = user[t]/n[t]
            s = sys[t]/n[t]
            r = rss[t]/n[t]/1024
            f = fd[t]/n[t]/1048576
            c = cpu[t]/n[t]
            speedup = base/w
            eff = speedup/t * 100
            printf "| %d | %.2f | %.2f | %.2f | %.0f | %.1f | %.0f | %.2fx | %.1f%% |\n", t, w, u, s, r, f, c, speedup, eff
        }
    }
}' "$CSV" >> "$REPORT_MD"

echo "" >> "$REPORT_MD"
echo "## Per-Run Raw Data" >> "$REPORT_MD"
echo "" >> "$REPORT_MD"
echo '```csv' >> "$REPORT_MD"
cat "$CSV" >> "$REPORT_MD"
echo '```' >> "$REPORT_MD"

echo "" >> "$REPORT_MD"
echo "## Notes" >> "$REPORT_MD"
echo "" >> "$REPORT_MD"
cat >> "$REPORT_MD" << 'NOTES'
- **Small dataset caveat**: The EXAMPLE dataset (~86 Mbp per genome) is small. At high thread counts, process spawning and I/O overhead may dominate, limiting observed scaling.
- **Full pipeline measured**: Each run includes FAtoGDB (single-threaded) + GIXmake (multi-threaded) + FastGA alignment (multi-threaded). The FAtoGDB phase is a fixed overhead regardless of thread count.
- **Temp file monitoring**: FD monitor polls every 0.5s. Very short-lived temp files may be missed between polls.
- **Peak Temp FD**: Measures the maximum total size of all open file descriptors related to FastGA temp files (including unlinked files that are still held open).
- **NUMA**: This server has 2 NUMA nodes. At T>16, cross-socket memory accesses may affect performance.
NOTES

echo ""
echo "============================================================"
echo "Benchmark complete!"
echo "============================================================"
echo "Report:     $REPORT_MD"
echo "CSV data:   $CSV"
echo "File sizes: $BENCH_DIR/results/file_sizes_report.txt"
echo "Timing:     $BENCH_DIR/timing/"
echo "FD logs:    $BENCH_DIR/fd_monitor/"
echo "FastGA logs: $BENCH_DIR/logs/"
