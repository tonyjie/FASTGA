# FastGA Benchmark: Performance & Thread Scaling

## Test Configuration

| Parameter | Value |
|---|---|
| Dataset | EXAMPLE/HAP1.fasta.gz vs EXAMPLE/HAP2.fasta.gz (~86 Mbp each) |
| Server | en-ec-zhang-x4 |
| CPU | AMD EPYC 9124 16-Core Processor (2 sockets x 16 cores = 64 threads) |
| Thread counts tested | 1, 2, 4, 8, 16, 32 |
| Repeats per config | 3 |
| Date | 2026-03-21 |
| Measurement tools | `/usr/bin/time -v`, FastGA `-L` log |

**Note**: T=64 failed because GIXmake enforces a maximum of 32 threads (`GIXmake: # of threads can be at most 32`). Results below cover T=1 through T=32.

## Pipeline Phases Explained

When you run `FastGA genome1.fa genome2.fa`, it orchestrates the full pipeline internally by spawning subprocesses. Here are the phases, the underlying tool for each, and the files they generate:

```
FastGA genome1.fa genome2.fa
  │
  ├─ Phase 1: GDB Creation     FAtoGDB genome1.fa  →  genome1.1gdb + .genome1.bps
  ├─ Phase 2: GIX Build        GIXmake genome1     →  genome1.gix + .genome1.ktab.*
  ├─ Phase 3: GDB Creation     FAtoGDB genome2.fa  →  genome2.1gdb + .genome2.bps
  ├─ Phase 4: GIX Build        GIXmake genome2     →  genome2.gix + .genome2.ktab.*
  │
  ├─ Phase 5: Seed Merge       Merge two GIX indices  →  _pair.* temp files
  ├─ Phase 6: Sort + Align     Sort seeds → chain → align  →  _uniq.*/_algn.* temp → .1aln
  │
  └─ Phase 7: PAF Conversion   ALNtoPAF  →  PAF text to stdout
```

| Phase | Tool | Threading | What It Does |
|---|---|---|---|
| **GDB Creation** | `FAtoGDB` (spawned via `system()`) | Single-threaded | Converts FASTA → `.1gdb` (ONEcode metadata) + `.bps` (2-bit compressed DNA). Also creates `.1ano` mask file if soft masking is detected. |
| **GIX Build** | `GIXmake` (spawned via `system()`) | Multi-threaded (max 32) | Builds sorted k-mer index. Sub-phases: (1) scan genome for syncmer-filtered k-mers, partition into distribution files; (2) MSD radix sort each partition, compress, output. |
| **Seed Merge** | Inside `FastGA.c` (Phase 1) | Multi-threaded | Linearly merges two GIX indices to find "adaptamers" — adaptive seeds where k-mers match between genomes. For each k-mer in genome A, finds the longest prefix match in genome B. Seeds with frequency <= `-f` (default 10) in genome B are emitted as position pairs. |
| **Seed Sort + Alignment** | Inside `FastGA.c` (Phases 2-4) | Multi-threaded | Three sub-steps bundled together: (1) **Reimport + radix sort**: reads seed pairs back, transforms to (anti-diagonal, diagonal-bucket) representation, in-place reverse MSD radix sort; (2) **Chaining + alignment**: linear scan to build chains of nearby seeds, runs wave-based local aligner (`Local_Alignment()`) on each chain; (3) **Sort + merge**: sorts per-thread overlap files, k-way merge into final output. |
| **PAF Conversion** | `ALNtoPAF` (spawned via `system()`) | Single-threaded | Converts compact trace-point encoded `.1aln` to PAF text format. Reconstructs alignment details from trace points on the fly. |

**Shortcut**: If you pre-build GDB and GIX files (via `-k` flag or manual `FAtoGDB` + `GIXmake` calls), then `FastGA genome1.gix genome2.gix` skips Phases 1-4 and jumps straight to the seed merge. This is recommended when comparing one genome against many others.

## Overall Thread Scaling (averaged over 3 repeats)

| Threads | Wall Clock (s) | User CPU (s) | Sys CPU (s) | Peak RSS (MB) | CPU% | Speedup | Efficiency |
|--------:|---------------:|-------------:|------------:|--------------:|-----:|--------:|-----------:|
| 1 | 133.5 | 101.2 | 28.9 | 513 | 97% | 1.00x | 100.0% |
| 2 | 73.3 | 101.2 | 28.9 | 522 | 177% | 1.82x | 91.1% |
| 4 | 49.2 | 101.1 | 29.2 | 552 | 265% | 2.71x | 67.9% |
| 8 | 32.6 | 103.2 | 29.3 | 653 | 406% | 4.09x | 51.1% |
| 16 | 25.1 | 106.6 | 31.7 | 682 | 552% | 5.32x | 33.3% |
| 32 | 22.0 | 121.4 | 37.1 | 740 | 720% | 6.07x | 19.0% |

### Key Observations

**Runtime Scaling:**
- Good scaling from T=1 to T=2 (91% efficient) and T=1 to T=4 (68% efficient)
- Diminishing returns beyond T=8: going from 8 to 32 threads only improves wall clock from 32.6s to 22.0s (1.48x for 4x more threads)
- The single-threaded FAtoGDB phase (~0.5s per genome = ~1s total) is a small fixed overhead
- At T=32, user CPU time increases to 121.4s (vs 101.2s at T=1), indicating thread coordination overhead
- System CPU time also increases significantly at T=32 (37.1s vs 28.9s at T=1)

**Memory Scaling:**
- Peak RSS grows modestly with threads: 513 MB (T=1) to 740 MB (T=32), a 1.44x increase for 32x threads
- This suggests per-thread buffers are relatively small compared to shared data structures

**Scaling Bottleneck:**
- The pipeline has a serial fraction: FAtoGDB is single-threaded, and the seed merge phase in FastGA has limited parallelism
- For this small dataset, Amdahl's law limits speedup: even with infinite threads, the serial portions prevent further improvement
- Larger genomes would likely show better scaling since the parallelizable alignment phase dominates

## Per-Phase Runtime Breakdown (rep1 data)

Wall clock time per phase at each thread count:

| Threads | GDB Total (s) | GIX Build (s) | Seed Merge (s) | Sort + Align (s) | PAF Conv (s) | Total (s) |
|--------:|---------------:|---------------:|----------------:|------------------:|--------------:|----------:|
| 1 | 1.0 | 18.9 | 30.0 | 83.1 | 0.5 | 133.5 |
| 2 | 1.1 | 11.1 | 16.5 | 43.9 | 0.3 | 72.8 |
| 4 | 1.0 | 7.5 | 11.2 | 28.7 | 0.2 | 48.5 |
| 8 | 1.0 | 6.2 | 8.3 | 16.6 | 0.2 | 32.3 |
| 16 | 1.0 | 6.2 | 7.1 | 9.8 | 0.2 | 24.3 |
| 32 | 1.0 | 7.7 | 5.7 | 7.1 | 0.2 | 21.7 |

Percentage of total wall time at each thread count:

| Threads | GDB | GIX Build | Seed Merge | Sort + Align | PAF |
|--------:|----:|----------:|-----------:|-------------:|----:|
| 1 | 0.7% | 14.1% | 22.5% | **62.3%** | 0.4% |
| 2 | 1.5% | 15.3% | 22.6% | **60.3%** | 0.4% |
| 4 | 2.0% | 15.5% | 23.0% | **59.1%** | 0.4% |
| 8 | 3.0% | 19.3% | 25.6% | **51.4%** | 0.6% |
| 16 | 4.1% | 25.4% | 29.2% | **40.4%** | 0.8% |
| 32 | 4.6% | **35.4%** | **26.4%** | **32.6%** | 1.0% |

## Per-Phase Detailed Resources (rep1 data)

### GDB Creation (HAP1 + HAP2 combined) — `FAtoGDB`, single-threaded

| Threads | User CPU (s) | Sys CPU (s) | Wall (s) | CPU% | Speedup | Efficiency |
|--------:|-------------:|------------:|----------:|-----:|--------:|-----------:|
| 1 | 0.87 | 0.06 | 0.99 | 93% | 1.00x | 100.0% |
| 2 | 0.83 | 0.09 | 1.07 | 86% | 0.93x | 46.4% |
| 4 | 0.85 | 0.08 | 0.98 | 95% | 1.01x | 25.3% |
| 8 | 0.82 | 0.10 | 0.98 | 94% | 1.01x | 12.6% |
| 16 | 0.82 | 0.10 | 0.98 | 94% | 1.01x | 6.3% |
| 32 | 0.83 | 0.10 | 0.99 | 94% | 1.00x | 3.1% |

Entirely single-threaded. Wall clock is constant ~1.0s regardless of thread count.

### GIX Build (HAP1 + HAP2 combined) — `GIXmake`, multi-threaded (max 32)

| Threads | User CPU (s) | Sys CPU (s) | Wall (s) | CPU% | Speedup | Efficiency |
|--------:|-------------:|------------:|----------:|-----:|--------:|-----------:|
| 1 | 14.13 | 2.60 | 18.89 | 89% | 1.00x | 100.0% |
| 2 | 13.89 | 2.47 | 11.14 | 147% | 1.70x | 84.8% |
| 4 | 13.90 | 2.47 | 7.48 | 219% | 2.53x | 63.1% |
| 8 | 15.57 | 2.40 | 6.23 | 289% | 3.04x | 37.9% |
| 16 | 18.34 | 2.58 | 6.17 | 339% | 3.06x | 19.1% |
| 32 | 21.46 | 3.04 | 7.68 | 319% | 2.46x | 7.7% |

Scales well up to T=4 (2.53x). Plateaus at T=8-16 (~6.2s). **Regresses at T=32** (7.68s) — user CPU inflates from 14.1s to 21.5s due to thread coordination overhead.

### Seed Merge — inside `FastGA.c`, multi-threaded

| Threads | User CPU (s) | Sys CPU (s) | Wall (s) | CPU% | Speedup | Efficiency |
|--------:|-------------:|------------:|----------:|-----:|--------:|-----------:|
| 1 | 3.81 | 25.22 | 30.05 | 97% | 1.00x | 100.0% |
| 2 | 3.83 | 25.17 | 16.45 | 176% | 1.83x | 91.3% |
| 4 | 3.87 | 25.37 | 11.15 | 262% | 2.69x | 67.4% |
| 8 | 4.16 | 25.45 | 8.27 | 358% | 3.64x | 45.5% |
| 16 | 4.42 | 27.06 | 7.11 | 443% | 4.23x | 26.4% |
| 32 | 4.91 | 31.14 | 5.73 | 630% | 5.25x | 16.4% |

Steady scaling but well below linear. **System CPU dominates** (~25s, largely from memory-mapping and I/O of the GIX index files), which limits parallelism. User CPU is tiny (~4s) and barely changes. This phase is I/O-bound.

### Seed Sort + Alignment — inside `FastGA.c`, multi-threaded

| Threads | User CPU (s) | Sys CPU (s) | Wall (s) | CPU% | Speedup | Efficiency |
|--------:|-------------:|------------:|----------:|-----:|--------:|-----------:|
| 1 | 82.42 | 0.68 | 83.14 | 100% | 1.00x | 100.0% |
| 2 | 81.99 | 0.91 | 43.87 | 189% | 1.90x | 94.8% |
| 4 | 82.14 | 0.89 | 28.66 | 290% | 2.90x | 72.5% |
| 8 | 82.41 | 0.84 | 16.59 | 502% | 5.01x | 62.7% |
| 16 | 81.96 | 1.11 | 9.82 | 846% | 8.46x | 52.9% |
| 32 | 92.64 | 2.29 | 7.05 | 1346% | 11.79x | 36.8% |

**The most parallelizable phase** and the dominant runtime contributor. Near-linear scaling up to T=2 (94.8% efficient). At T=8, 5.01x speedup (62.7% efficiency). At T=32, 11.79x but user CPU inflates from 82.4s to 92.6s. This is almost entirely CPU-bound (sys CPU < 1s at low thread counts), making it the primary beneficiary of threading.

### PAF Conversion — `ALNtoPAF`, single-threaded

| Threads | User CPU (s) | Sys CPU (s) | Wall (s) | CPU% | Speedup | Efficiency |
|--------:|-------------:|------------:|----------:|-----:|--------:|-----------:|
| 1 | 0.00 | 0.00 | 0.47 | 1% | 1.00x | 100.0% |
| 2 | 0.00 | 0.00 | 0.31 | 1% | 1.50x | 75.1% |
| 4 | 0.00 | 0.00 | 0.22 | 3% | 2.13x | 53.2% |
| 8 | 0.00 | 0.01 | 0.19 | 7% | 2.44x | 30.4% |
| 16 | 0.00 | 0.01 | 0.21 | 9% | 2.23x | 13.9% |
| 32 | 0.00 | 0.02 | 0.22 | 10% | 2.11x | 6.6% |

Trivially fast (<0.5s). Near-zero CPU — wall time is I/O latency. Not a scaling target.

## Phase Scaling Summary

- **Seed sort + alignment** is the dominant phase (62% of T=1 runtime) and scales the best: **11.8x speedup at T=32**. It's CPU-bound and highly parallelizable.
- **Seed merge** is the second most expensive (22% of T=1 runtime) but is **I/O-bound** (sys CPU ~25s vs user CPU ~4s), capping its speedup at **5.3x at T=32**.
- **GIX build** scales up to T=8 (3.0x speedup) then **regresses at T=32** due to thread coordination overhead.
- **GDB creation** and **PAF conversion** are negligible and single-threaded.
- At T=32, time is evenly split across GIX build (35%), seed merge (26%), and sort+align (33%) — no single phase dominates anymore, which is why overall scaling flattens.
- **Best overall sweet spot**: T=8 to T=16, where the two heavy phases both show good speedup without excessive overhead.

## Per-Run Raw Data

| Threads | Rep | Wall (s) | User (s) | Sys (s) | RSS (KB) | CPU% |
|--------:|----:|---------:|---------:|--------:|---------:|-----:|
| 1 | 1 | 133.84 | 101.60 | 28.84 | 525524 | 97% |
| 1 | 2 | 133.68 | 101.29 | 28.90 | 525520 | 97% |
| 1 | 3 | 133.02 | 100.69 | 28.90 | 525480 | 97% |
| 2 | 1 | 73.13 | 101.01 | 28.93 | 535148 | 177% |
| 2 | 2 | 73.59 | 101.25 | 28.97 | 534964 | 176% |
| 2 | 3 | 73.09 | 101.34 | 28.65 | 534900 | 177% |
| 4 | 1 | 48.80 | 101.23 | 29.11 | 565664 | 267% |
| 4 | 2 | 48.83 | 101.06 | 29.11 | 565516 | 266% |
| 4 | 3 | 49.85 | 101.15 | 29.49 | 565672 | 262% |
| 8 | 1 | 32.54 | 103.49 | 29.09 | 668912 | 407% |
| 8 | 2 | 32.72 | 103.20 | 29.36 | 668952 | 405% |
| 8 | 3 | 32.63 | 102.97 | 29.32 | 668904 | 405% |
| 16 | 1 | 24.60 | 106.14 | 31.18 | 698644 | 558% |
| 16 | 2 | 26.41 | 107.92 | 33.17 | 698644 | 534% |
| 16 | 3 | 24.19 | 105.81 | 30.75 | 698552 | 564% |
| 32 | 1 | 22.04 | 120.43 | 37.02 | 757752 | 714% |
| 32 | 2 | 21.90 | 121.62 | 37.45 | 757856 | 726% |
| 32 | 3 | 22.03 | 122.28 | 36.76 | 757520 | 721% |
