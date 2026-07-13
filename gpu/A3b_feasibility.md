# A3b feasibility: can a batched GPU beat FastGA's CPU aligner? (2026-07-13)

Before committing to the multi-day `search_seeds` batching redesign, measure whether the
batched GPU kernels can actually beat the CPU align phase. Benchmarks over **all 322,531
alignments** of the EXAMPLE part (`gpu/trace_bench.cu`, `gpu/disc_bench.cu`; resident
contigs, one batched call, best of 3).

## Trace kernel (warp-parallel) — FAST, occupancy-bound

| concurrent warps (TR_CHUNK) | wave scratch | time | aln/s | vs CPU align (T=8) |
|---:|---:|---:|---:|---:|
| 256  | 1 GB  | 13.15 s | 24,527  | 1.3× |
| 512  | 2 GB  | 7.98 s  | 40,423  | 2.1× |
| 1024 | 4 GB  | 4.74 s  | 68,012  | 3.5× |
| 2048 | 8 GB  | 2.82 s  | 114,236 | 5.9× |
| 4096 | 16 GB | 1.74 s  | 185,393 | 9.6× |
| 8192 | 32 GB | **1.06 s** | **305,787** | **15.8×** |

The trace kernel is **occupancy-limited**, not compute-limited: throughput scales ~linearly
with the number of resident warps, which the per-warp wave scratch (4 MB) caps. The A100's
80 GB easily holds 8192 warps (32 GB) alongside the ~1.2 GB resident contigs. Trace is **not
the bottleneck** — the earlier A3a slowness was pure per-tube launch granularity.

## Discovery kernel (1 thread per task) — the current bottleneck

100,000 discoveries in 5.58 s = **17,918 disc/s** — about the CPU align rate. This kernel is
*not* warp-parallel (`disc_batch` runs one thread per task, serial x-drop). Warp-parallelizing
it with the same band-across-lanes pattern the trace kernel already uses is a **bounded**
change and should give a similar ~10–20× → ~200–350k disc/s.

## Verdict

**The align PHASE can be won on the GPU.** Reference: CPU align phase = 16.7 s for 323 k
alignments (T=8), ~19,400 aln/s.
- Trace: **15.8×** already (1.06 s), no algorithmic work needed — just scratch.
- Discovery: at parity now, but warp-parallelizable to match (bounded, proven pattern).
- Combined GPU wave (both warp-parallel, occupancy-tuned): **~2–4 s** vs the CPU wave's ~14 s.

**But the END-TO-END speedup is Amdahl- and dataset-limited:**
- EXAMPLE total = seed-merge 41 s + align 16.7 s. Align is only **29%** of runtime, so even a
  zero-cost align caps the total speedup at **1.4×**. EXAMPLE is *not* representative.
- On align-heavy mammalian genomes (M1: sort+align is 80–98% of runtime), the same kernels
  would translate to a **~3–5× end-to-end** win. That is the regime the acceleration targets.

## What A3b actually requires (now de-risked)
1. **Warp-parallelize the discovery kernel** (bounded; copy the trace kernel's lane pattern).
2. **Large trace scratch** for occupancy (done: `TR_CHUNK` tunable; 8 GB → 5.9×, 32 GB → 15.8×).
3. **`search_seeds` batching redesign** — collect first-tubes across a whole part, one big
   discover+trace batch, iterate continuations (`alow=eant`), preserve inline redundancy
   removal. **This is the multi-day cost and the real risk** (orchestration overhead, load
   balance, keeping the GPU saturated across contig-pairs).

**Recommendation:** the kernels clear the bar; pursue A3b **if** the target is align-heavy
genomes (mammalian WGA), where the end-to-end payoff justifies the redesign. For EXAMPLE-scale
or seed-merge-dominated workloads, it will not move the total materially.
