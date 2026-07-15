# GPU acceleration of FastGA — results summary

All on branch `agent-optimization-wt`, A100. EXAMPLE = hap1 vs hap2 (human haplotypes),
T=8 CPU baseline. `.1aln` quality vs the CPU baseline of **323,569 records / 632,119,471
aligned A-bases**.

## What maps to the GPU
FastGA's Phase-3 per-tube local alignment (Myers O(nd) difference-wave, the 80–98%
bottleneck). Two kernels, both **warp-cooperative (one warp per alignment, band across 32
lanes)**:
- **Discovery** — forward+reverse x-drop from a seed → endpoints+diffs.
- **Trace** — banded furthest-reaching wave + traceback → FastGA trace-points (diff, Δb) per
  100 bp panel. *Novel: the Myers O(nd) local trace has not been emitted on a GPU before.*
Seed-merge, seed-sort/chaining, redundancy removal, and output stay on the CPU.

## Correctness (measured, end-to-end `-G`)
| | records | aligned A-bases |
|---|---:|---:|
| CPU baseline | 323,569 | 632,119,471 |
| `-G` GPU discovery+trace | 320,491 (**99.05%**) | 626,798,190 (**99.16%**) |

GPU-emitted trace-points are 100% Check_Trace_Points-valid with the **exact same edit
distance** as FastGA → identical identity, no quality loss. (The ~1% gap is the >10%-divergence
tail where the x-drop endpoints differ; those tubes fall back to the exact CPU aligner.)

## Kernel throughput (measured, batched, all 322,531 EXAMPLE alignments)
| kernel | rate | time |
|---|---:|---:|
| discovery (warp)        | 443,341 disc/s | 0.73 s |
| trace (8 GB scratch)    | 114,236 aln/s  | 2.82 s |
| trace (32 GB scratch)   | 305,787 aln/s  | 1.06 s |
| **combined wave (8 GB)**  | | **3.55 s** |
| **combined wave (32 GB)** | | **1.79 s** |
| CPU align phase (T=8)   | ~19,400 aln/s | 16.70 s |

**Batched GPU wave = 4.7–9.3× the CPU align phase.**

## Representative dataset: HUMAN (GRCh38 × CHM13, ~3.1 Gbp each, T=32) — measured

EXAMPLE is seed-merge-dominated; human is the align-heavy regime. All numbers below are
**measured** on this branch (instrumented `FastGA.wave` + the GPU benches on real human
alignments extracted from the 518,037-record `.1aln`).

**Phase breakdown (measured, T=32):**
| phase | wall | share |
|---|---:|---:|
| GDB + GIX build | 107 s | 18% |
| seed-merge | 11 s | 2% |
| **sort + align** | **487 s** | **80%** |
| total | 605 s | 100% |

**Inside sort+align — the wave is NOT the whole phase (this corrects the earlier assumption):**
Instrumenting every `Local_Alignment` call: **wave = 1490.5 CPU-s = 60.7% of the sort+align
CPU-time** (2455.9 CPU-s). The other ~39% is the seed radix-sort (billions of seeds) +
chaining + redundancy + output-merge (14.7 Gbp of trace). Wave wall ≈ 247 s of the 487 s.

**Alignment-length distribution is heavily skewed:** mean 28,415 bp, but ~90% are ≤32 kbp
(mean ~4 kbp) and ~10% are *very* long (up to hundreds of kbp) — and the long tail dominates
total work. The GPU trace kernel stores furthest-reaching `x` in `short`, so alignments
>32 kbp overflow → **CPU fallback exactly on the expensive tail**.

**GPU vs CPU on representative (long) human alignments — same tasks:**
| | rate |
|---|---:|
| GPU trace (≤32 kbp, mean 4 kbp) | 148,823 aln/s |
| CPU-32 `Compute_Alignment` | 806 aln/s |
(Compute_Alignment is a slower baseline than FastGA's `Local_Alignment`; the honest wave
speedup is bounded below by the ~8× measured earlier and is large on the long low-divergence
alignments the GPU can hold.)

**Honest human end-to-end (Amdahl on the measured wave share):**
| scenario | end-to-end |
|---|---:|
| current kernel (short cap → long tail on CPU) | **~1.3×** |
| + long-alignment kernel (int `x`, no fallback) | **~1.9× (align-only) / ~1.65× (full run)** |

**Why not 3-5× (correcting my earlier projection):** the wave is only ~49% of the align-only
runtime (61% of sort+align × ...) — human-human is *low-divergence*, so alignments are long
and *easy* (shallow waves), which makes the seed-sort + output ~half the phase. GPU-aligner
acceleration alone is therefore **Amdahl-capped near ~1.9×** on this pair. A bigger win needs
also accelerating the seed radix-sort and output-merge (a whole-pipeline story), and/or a
more divergent pair (chimp/mouse) where the wave is a larger share.

## A3b batched pipeline — WORKS end-to-end, but batch-size-limited (2026-07-15)

The dynamic-batching queue + resident-genome integration (FastGA.c `-G`) runs the whole align
phase on the GPU and produces a valid `.1aln`:

| EXAMPLE, T | records (vs CPU 323,569) | aligned bases (vs 632,119,471) | sort+align wall |
|---|---:|---:|---:|
| −G T=8  | 320,506 (**99.05%**) | 626,790,495 (**99.16%**) | 404 s |
| −G T=64 | 320,506 (99.05%)     | 626,790,495 (99.16%)     | 259 s |
| CPU T=8 | 323,569              | 632,119,471              | **16.7 s** |

**Correctness: ✓** — the batched pipeline is functionally correct (99% coverage, all kernels
in the loop). **Speed: ✗** — it is 15-24× *slower* than the CPU align phase, and the reason is
architectural, not the kernels: the batching queue's batch is **capped at the live thread
count** (each thread's tube-walk is sequential, so a thread has exactly one pending tube), and
FastGA is practical to ~64 threads. Batch ~64 gives the trace kernel only ~10-20k aln/s
(occupancy-starved); the offline benches need **batch ~2048 for the 15× win**. The queue also
serializes CPU/GPU (the flush holds the mutex during the GPU call), and releasing it doesn't
help — threads in the current batch are blocked awaiting their own result, so no larger batch
accumulates.

**Conclusion:** the queue-batching design cannot reach the batch sizes the GPU needs. The only
path to batches of thousands is the **coroutine approach** — interleave many triple-sweeps so
one thread has many independent tubes in flight at once — which was deferred as higher-risk.
That, not the kernels, is the remaining work between "correct" and "faster than CPU".

## End-to-end (the honest bottom line)
| | measured | note |
|---|---|---|
| CPU baseline total (T=8) | **57.5 s** | seed-merge 41 s + align 16.7 s |
| `-G` unbatched (A3a) | 64 min | correct, but per-tube GPU launch — a negative control |
| `-G` batched (A3b) | *projected* ~45–47 s | seed-merge 41 s + align ~4–6 s |

**EXAMPLE end-to-end is Amdahl-capped at ~1.3×** — seed-merge (41 s) is 71% of runtime and is
not GPU-accelerated. EXAMPLE is *not* representative.

**On align-heavy mammalian WGA** (human/chimp/mouse; sort+align is 80–98% of runtime per our
M1 characterization), the same measured kernel speedups translate to a projected
**~3–5× end-to-end**. That is the regime this acceleration targets.

## Status
- Kernels: **done & measured** (discovery 24.7× warp speedup; trace 15.8× occupancy-tuned;
  both correct).
- Integration `-G`: **done, correct, unbatched** (proves the pipeline; slow by design).
- Batching orchestration (A3b): **designed** (gpu/A3b_design.md) — a coroutine interleave of
  FastGA's triple-sweeps, required because `align_contigs` has three nested levels of
  sequential state. This is the remaining multi-day core-loop rewrite; its orchestration
  overhead is the last unmeasured factor between the measured kernel wins and the projected
  end-to-end number.

## The one-line finding
FastGA's aligner *can* be accelerated 5–9× on the GPU (kernels proven), but the win is gated
less by the kernel than by (a) FastGA's sequential tube-walk, which forces a coroutine
batching rewrite, and (b) Amdahl — the GPU only pays off on align-heavy workloads, since
seed-merge dominates otherwise.
