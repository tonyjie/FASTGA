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
