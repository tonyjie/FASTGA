# FastGA Thread Scaling — New Upstream Build (2026-07-09)

Thread scaling study of the **current upstream FastGA** on the repo EXAMPLE
dataset. This is a fresh run of the same methodology used in
[`old_archive/benchmark_performance.md`](old_archive/benchmark_performance.md), which recorded the
`optimize-memory` build back on 2026-03-21. Use this doc for the new-build
numbers; the older doc remains the record of the previous build.

> **Companion storage study:** [`benchmark_storage_upstream.md`](benchmark_storage_upstream.md)
> compares the disk footprint of the old vs new upstream (verdict: byte-identical, thread-independent).

## Build Under Test

| Item | Value |
|---|---|
| Branch / commit | `main` @ `10ebff7` = upstream `ddeea32` ("Fixed a bug for the uninitialised field ANO_PAIR.parse") + local `.gitignore` |
| Source | `thegenemyers/FASTGA` upstream/main, freshly `make clean && make` |
| Note | This is the **stock upstream** aligner — it does NOT contain the `optimize-memory` Opt1/Opt3 storage changes. |

## Test Configuration

| Parameter | Value |
|---|---|
| Dataset | EXAMPLE/HAP1.fasta.gz vs EXAMPLE/HAP2.fasta.gz (~86 Mbp each) |
| Server | en-ec-zhang-x4 (same machine as the 2026-03-21 study) |
| CPU | AMD EPYC 9124 16-Core (2 sockets x 16 cores = 64 logical threads) |
| Thread counts tested | 1, 2, 4, 8, 16, 32 (and 64, which fails — see below) |
| Repeats per config | 3 (median reported) |
| Date | 2026-07-09 |
| Measurement | `/usr/bin/time -v` for wall/CPU/RSS; FastGA `-L:` log for phases |
| Command | `FastGA -v -T<N> -P. -L:<log> HAP1.fasta.gz HAP2.fasta.gz > out.paf` |
| Method | Generated GDB/GIX/index files removed before every run, so each measurement is a full end-to-end build+align. |

**T=64 fails**: `GIXmake: # of threads can be at most 32, more doesn't help.`
FastGA hard-caps threads at 32. Results cover T=1..32.

## Overall Thread Scaling (median of 3 repeats)

| Threads | Wall (s) | User (s) | Sys (s) | CPU% | Peak RSS (MB) | Speedup | Efficiency | Non-redundant aln |
|--------:|---------:|---------:|--------:|-----:|--------------:|--------:|-----------:|------------------:|
| 1  | 128.4 | 101.4 | 27.2 | 100% | 513 | 1.00x | 100.0% | 323,569 |
| 2  | 69.5  | 101.2 | 28.0 | 185% | 523 | 1.85x | 92.5%  | 323,569 |
| 4  | 45.6  | 101.8 | 27.9 | 284% | 552 | 2.82x | 70.4%  | 323,569 |
| 8  | 28.5  | 102.6 | 27.8 | 457% | 653 | 4.50x | 56.3%  | 323,569 |
| 16 | 20.2  | 105.4 | 28.2 | 662% | 682 | 6.37x | 39.8%  | 323,569 |
| 32 | 16.2  | 122.5 | 31.1 | 946% | 740 | 7.94x | 24.8%  | 323,569 |

![Thread scaling: new upstream build](thread_scaling_new_upstream.png)

*Left: end-to-end wall clock (log-log). Right: speedup vs T=1 against ideal
linear; the 32-thread hard cap bounds the x-axis.*

## Key Observations

**Correctness invariant to thread count.** All runs at every thread count
produced exactly **323,569 non-redundant alignments** (ave len 1953 bp). Thread
count does not change results — a good multi-threading correctness check.

**Sub-linear scaling, Amdahl-bounded.** Near-linear through T=4 (2.82x, 70%
efficient), then diminishing: T=32 reaches 7.94x (25% efficient). At T=32 the
average CPU utilisation is only ~946% (~9.5 of 32 cores busy on average),
because serial-heavy phases dominate on this small dataset:
- `FAtoGDB` (GDB creation) is single-threaded (~0.4s x 2 genomes).
- `ALNtoPAF` (PAF conversion) is single-threaded.
- The seed-merge phase has limited parallelism.

On the human genome (3.1 Gbp) the parallel alignment phase dominates, so
scaling there is materially better (see the human T=32 section of
`old_archive/benchmark_performance.md`).

**Memory grows modestly.** Peak RSS 513 MB (T=1) → 740 MB (T=32), 1.44x for 32x
threads — per-thread buffers are small relative to shared structures.

## Comparison vs the 2026-03-21 optimize-memory build

Same machine, same dataset, same method — but **different sessions/dates**, so
system load was not controlled. Treat this as indicative, not a rigorous
back-to-back comparison.

| Threads | Wall new (s) | Wall old (s) | Speedup new | Speedup old |
|--------:|-------------:|-------------:|------------:|------------:|
| 1  | 128.4 | 133.7 | 1.00x | 1.00x |
| 2  | 69.5  | 73.1  | 1.85x | 1.83x |
| 4  | 45.6  | 48.8  | 2.82x | 2.74x |
| 8  | 28.5  | 32.6  | 4.50x | 4.10x |
| 16 | 20.2  | 24.6  | 6.37x | 5.43x |
| 32 | 16.2  | 22.0  | 7.94x | 6.07x |

The new upstream build is consistently faster and scales somewhat better,
most visibly at high thread counts (T=32: 16.2s vs 22.0s, ~26% faster). This is
suggestive of upstream improvements between the two builds, but a controlled
alternating re-run on an idle node would be needed to confirm it rather than
attribute it to load differences.

## Reproduce

Raw data and scripts are in
[`benchmark_thread_scaling_upstream/`](benchmark_thread_scaling_upstream/):
- `run_thread_scaling.sh` — the driver (edit `THREADS`/`REPS`)
- `results.tsv` — parsed per-run metrics
- `analyze.py` — builds the median table and the figure

```bash
bash docs/benchmark_thread_scaling_upstream/run_thread_scaling.sh
python3 docs/benchmark_thread_scaling_upstream/analyze.py
```
