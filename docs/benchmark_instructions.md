# Benchmark Instructions

How to reproduce the performance and storage benchmarks documented in `benchmark_performance.md` and `benchmark_storage.md`.

## Prerequisites

- FastGA built in repo root (`make`)
- EXAMPLE dataset present: `EXAMPLE/HAP1.fasta.gz` and `EXAMPLE/HAP2.fasta.gz`
- For human genome benchmark: FASTA files on scratch (see Section 4 below)
- Python 3 with matplotlib (for plotting)

## Scripts Overview

All scripts are in `benchmarks/`:

| Script | Purpose | Output |
|---|---|---|
| `run_benchmark.sh` | EXAMPLE dataset thread scaling study (T=1,2,4,8,16,32 x 3 repeats). Measures wall clock, CPU time, peak RSS per run. Also measures persistent file sizes with `-k`. | `benchmarks/results/benchmark_summary.csv`, `benchmarks/results/file_sizes_report.txt` |
| `run_storage_audit_v2.sh` | EXAMPLE dataset: measures peak temp file disk usage at each thread count using a dedicated tmpdir. | `benchmarks/storage_audit_v2/summary.csv`, per-thread `*_tmpdir_monitor.tsv` files |
| `run_human_benchmark.sh` | Human genome benchmark: GRCh38 vs CHM13 at T=32. Cleans old cached files, monitors disk usage, measures persistent and temp file sizes. | Results under `/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/` |
| `monitor_tmpdir.sh` | Background monitor used by audit/benchmark scripts. Polls `du -sb` on a directory at 50ms intervals to capture disk blocks consumed by unlinked temp files. | TSV file with columns: `elapsed_s`, `du_bytes`, `ls_bytes`, `ls_files` |
| `plot_storage_timeline.py` | Generates EXAMPLE storage timeline figures from `storage_audit_v2` monitor data. | `docs/storage_timeline_T08.png`, `docs/storage_timeline_comparison.png` |
| `plot_human_storage_timeline.py` | Generates human genome storage timeline figure from disk free space polling data. | `docs/storage_timeline_human_T32.png` |

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

### 4. Human Genome Benchmark (GRCh38 vs CHM13)

```bash
cd /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA
bash benchmarks/run_human_benchmark.sh
```

- Runs FastGA at T=32 comparing GRCh38 (~3.2 GB, 705 contigs) vs CHM13 (~3.0 GB, 24 contigs)
- **Data location**: Input FASTA and all generated files on `/scratch/jl4257/seq_align/fastga_datasets/`
  - GRCh38 FASTA: `/scratch/.../GRCh38/GCF_000001405.40_GRCh38.p14_genomic.fna`
  - CHM13 FASTA: `/scratch/.../CHM13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna`
  - GDB/GIX files created alongside input FASTA in their respective dirs
  - Temp files, logs, output → `/scratch/.../GRCh38_vs_CHM13/`
- **Cleans old cached files** before starting (GDB, GIX, ktab, bps, 1ano, partial results)
- **Safety checks**: Aborts if < 150 GB free before start; kills FastGA if free space drops below 50 GB during run
- **Monitors disk usage** via `df` every 30s (since tmpdir monitor doesn't capture temp files when FastGA writes them relative to the input FASTA path rather than the `-P` dir)
- Takes ~10 minutes at T=32
- Results: file sizes report, `/usr/bin/time -v` output, verbose per-phase log, disk usage polling

### 5. Generate Human Genome Storage Figure

```bash
cd /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA
python3 benchmarks/plot_human_storage_timeline.py
```

- Reconstructs the storage timeline from per-phase wall clock times and disk free space polling data
- Phase boundaries and disk polling data are hardcoded in the script from the benchmark run
- Generates `docs/storage_timeline_human_T32.png`

## Key Methodological Notes

### Why a dedicated tmpdir?

FastGA's temp files (`_pair.*`, `_uniq.*`, `_algn.*`) are `open()`-ed then immediately `unlink()`-ed. They consume disk blocks but are invisible to `ls` or `/proc/pid/fd` monitoring. By pointing `-P` to a dedicated empty directory, `du -sb` on that directory captures the actual block usage including unlinked files — because the filesystem tracks blocks until all FDs are closed.

### What the monitor captures vs. misses

- **Captures**: All files written to the `-P` tmpdir — seed pair temps, alignment temps, GIXmake distribution temps
- **Does not capture**: Persistent files in the working directory (GDB, GIX) — these are measured separately via `ls -la`
- **Polling interval**: 50ms. Very short-lived files (created and closed within 50ms) may be missed, but the seed pair files persist for seconds to minutes

### Per-phase timing extraction

The verbose logs contain `Resources for phase:` lines with format `Xu Ys Zw P%` where X=user CPU, Y=sys CPU, Z=wall clock, P=CPU utilization. Times may be in `min:sec.ms` format (e.g., `1:22.414u` = 82.414s). The phases appear in order: GDB genome1, GIX genome1, GDB genome2, GIX genome2, Seed Merge, Sort+Align, PAF Conversion.

- EXAMPLE logs: `benchmarks/logs/T*_rep*.log`
- Human logs: `/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/logs/verbose.log`

### Human genome benchmark: disk monitoring approach

For the human genome run, the tmpdir monitor (`monitor_tmpdir.sh`) did not capture temp files because FastGA writes `_pair.*` files relative to the working directory rather than the `-P` path in some cases. Instead, we monitor total disk consumption via `df` polling at 30s intervals. The temp file overhead is computed as: `(BASELINE_FREE - current_free) - persistent_file_sizes`.

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
  benchmark_performance.md         # Performance analysis (EXAMPLE + human genome)
  benchmark_storage.md             # Storage analysis (EXAMPLE + human genome)
  storage_timeline_T08.png         # EXAMPLE dataset timeline (T=8)
  storage_timeline_comparison.png  # EXAMPLE dataset thread comparison (T=1,4,8,32)
  storage_timeline_human_T32.png   # Human genome timeline (T=32)
  benchmark_instructions.md        # This file

# Human genome benchmark results (on scratch, not in repo)
/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/
  output.1aln                      # Alignment output
  file_sizes_report.txt            # Persistent file sizes
  logs/verbose.log                 # Per-phase timing
  logs/fastga_internal.log         # FastGA -L metrics
  timing/human_T32.time            # /usr/bin/time output
  monitor/human_T32_tmpdir.tsv     # Tmpdir monitor (limited data for human run)
```
