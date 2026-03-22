# Benchmark Instructions

How to reproduce the performance and storage benchmarks documented in `benchmark_performance.md` and `benchmark_storage.md`.

## Prerequisites

- FastGA built in repo root (`make`)
- EXAMPLE dataset present: `EXAMPLE/HAP1.fasta.gz` and `EXAMPLE/HAP2.fasta.gz`
- Python 3 with matplotlib (for plotting)

## Scripts Overview

All scripts are in `benchmarks/`. Only these 4 are current:

| Script | Purpose | Output |
|---|---|---|
| `run_benchmark.sh` | Thread scaling study (T=1,2,4,8,16,32 x 3 repeats). Measures wall clock, CPU time, peak RSS per run. Also measures persistent file sizes with `-k`. | `benchmarks/results/benchmark_summary.csv`, `benchmarks/results/file_sizes_report.txt`, `docs/benchmark_results.md` (auto-generated) |
| `run_storage_audit_v2.sh` | Measures peak temp file disk usage at each thread count using a dedicated tmpdir. | `benchmarks/storage_audit_v2/summary.csv`, per-thread `*_tmpdir_monitor.tsv` files |
| `monitor_tmpdir.sh` | Background monitor used by `run_storage_audit_v2.sh`. Polls `du -sb` on a directory at 50ms intervals to capture disk blocks consumed by unlinked temp files. | TSV file with columns: `elapsed_s`, `du_bytes`, `ls_bytes`, `ls_files` |
| `plot_storage_timeline.py` | Generates storage timeline figures from `storage_audit_v2` monitor data. | `docs/storage_timeline_T08.png`, `docs/storage_timeline_comparison.png` |

## How to Run

### 1. Thread Scaling Benchmark

```bash
cd /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA
bash benchmarks/run_benchmark.sh
```

- Runs FastGA at T=1,2,4,8,16,32 with 3 repeats each (21 runs + 1 file size run)
- Takes ~15-20 minutes for the EXAMPLE dataset
- Results in `benchmarks/results/benchmark_summary.csv` and `benchmarks/timing/`
- Per-phase timing is in verbose logs at `benchmarks/logs/T*_rep*.log`
- Persistent file sizes from the `-k` run are in `benchmarks/results/file_sizes_report.txt`

### 2. Storage Audit

```bash
cd /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA
bash benchmarks/run_storage_audit_v2.sh
```

- Runs FastGA at T=1,2,4,8,16,32 with `-P <dedicated_tmpdir>` to isolate temp files
- `monitor_tmpdir.sh` runs in background polling `du -sb` at 50ms intervals
- This captures unlinked-but-open temp files (`_pair.*`, `_uniq.*`, `_algn.*`) that are invisible to `ls`
- Takes ~10-15 minutes for the EXAMPLE dataset
- Summary in `benchmarks/storage_audit_v2/summary.csv`
- Per-thread monitor logs in `benchmarks/storage_audit_v2/T*_tmpdir_monitor.tsv`

### 3. Generate Figures

```bash
cd /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA
python3 benchmarks/plot_storage_timeline.py
```

- Reads `benchmarks/storage_audit_v2/T*_tmpdir_monitor.tsv`
- Generates `docs/storage_timeline_T08.png` (detailed T=8 timeline) and `docs/storage_timeline_comparison.png` (T=1,4,8,32 comparison)
- Phase boundaries in the plot script are hardcoded for the EXAMPLE dataset timing — update `PHASES_T08` and `phase_timings` dict if re-running on different data

## Key Methodological Notes

### Why a dedicated tmpdir?

FastGA's temp files (`_pair.*`, `_uniq.*`, `_algn.*`) are `open()`-ed then immediately `unlink()`-ed. They consume disk blocks but are invisible to `ls` or `/proc/pid/fd` monitoring. By pointing `-P` to a dedicated empty directory, `du -sb` on that directory captures the actual block usage including unlinked files — because the filesystem tracks blocks until all FDs are closed.

### What the monitor captures vs. misses

- **Captures**: All files written to the `-P` tmpdir — seed pair temps, alignment temps, GIXmake distribution temps
- **Does not capture**: Persistent files in the working directory (GDB, GIX) — these are measured separately via `ls -la`
- **Polling interval**: 50ms. Very short-lived files (created and closed within 50ms) may be missed, but the seed pair files persist for seconds to minutes

### Per-phase timing extraction

The verbose logs (`benchmarks/logs/T*_rep*.log`) contain `Resources for phase:` lines with format `Xu Ys Zw P%` where X=user CPU, Y=sys CPU, Z=wall clock, P=CPU utilization. Times may be in `min:sec.ms` format (e.g., `1:22.414u` = 82.414s). The phases appear in order: GDB HAP1, GIX HAP1, GDB HAP2, GIX HAP2, Seed Merge, Sort+Align, PAF Conversion.

## Output Files Summary

After running all benchmarks:

```
benchmarks/
  results/
    benchmark_summary.csv          # Thread scaling: threads, wall, user, sys, rss, peak_fd, cpu%
    file_sizes_report.txt          # Persistent file sizes (GDB, GIX, ktab breakdown)
  storage_audit_v2/
    summary.csv                    # Peak storage per thread count
    T01_tmpdir_monitor.tsv         # Temp disk usage timeline, T=1
    T02_tmpdir_monitor.tsv         # ...
    T04_tmpdir_monitor.tsv
    T08_tmpdir_monitor.tsv
    T16_tmpdir_monitor.tsv
    T32_tmpdir_monitor.tsv
  logs/
    T*_rep*.log                    # Verbose stderr with per-phase timing

docs/
  benchmark_performance.md         # Analysis document (manually written)
  benchmark_storage.md             # Analysis document (manually written)
  storage_timeline_T08.png         # Generated figure
  storage_timeline_comparison.png  # Generated figure
  benchmark_instructions.md        # This file
```
