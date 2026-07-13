# Phase A: GPU trace-point emission (the actual speedup)

**Goal:** `-G` end-to-end *faster* than CPU T=32 by having the GPU emit FastGA
trace-points directly, removing the CPU `Compute_Alignment` that makes the current
correctness-first `-G` ~160× slower.

## Why this is the lever
The current `-G` aligns every tube **twice**: GPU x-drop for endpoints, then CPU
`Compute_Alignment` (`split_nd`) to build the trace. If the GPU emits the trace-points,
the CPU only *packs* them (`Compress_TraceTo8`) — the alignment work happens once, on the
GPU, batched.

## The trace-point contract (confirmed, align.c/alncode.c)
`path->trace` = `uint16` pairs, one per A-panel k covering A-coords `[k*100,(k+1)*100)`:
even = diffs in panel (subs+ins+del, unit cost), odd = Δb (B-bases consumed).
`tlen = 2*(ceil(aepos/100) − floor(abpos/100))`; `sum(Δb)=bepos−bbpos`,
`sum(diffs)=path->diffs`. Build by walking the base-level path with an A-cursor, binning
`(diff,Δb)` by `floor(a/100)`. Must pass `Check_Trace_Points` (align.c:4006).

## Algorithm (per alignment, one warp)
1. Endpoints `(abpos,bbpos)→(aepos,bepos)` come from the existing x-drop discovery.
2. **Banded O(ND) furthest-reaching wave**, forward from `(abpos,bbpos)`, band width
   `2*KBAND+|drift|` (< ~768, structurally bounded — M2), storing the furthest x per
   `(d,k)` in **global memory** (short; ~D_max×band ≤ ~1.5 MB/alignment).
3. Stop when the wave reaches `(aepos,bepos)` at `d = path->diffs`.
4. **Traceback** through the stored wave (warp lane 0, or lane-parallel segments):
   recover match/sub/ins/del ops, maintaining the A-cursor, and accumulate
   `(diff,Δb)` into panel bins `floor(a/100)`. Emit the uint16 pair vector.
5. Handle the boundary-fold (align.c:5557) so inserts on a TP boundary land left.

Memory bound: batch B alignments → B×1.5 MB. A100 80 GB fits thousands per launch;
cap batch size and loop. Overflow tubes (band or D too large) → CPU fallback.

## Milestones (each independently testable)

### A1 — Offline trace kernel + validation harness  ← START HERE (highest risk, self-contained)
- New kernel `trace_kernel_warp` (extend `gpu_align.cu`): banded wave with global-mem
  wave store + traceback → uint16 trace-point vector.
- Harness: extend `extract_tasks`/a new `extract_trace` to dump, per real FastGA
  alignment, the endpoints + FastGA's own trace-point vector (from `Compute_Alignment`).
- Validate: run the GPU kernel on those tasks, compare the emitted vector to FastGA's
  pair-by-pair. Metrics: % panels exact `(diff,Δb)`; `sum(Δb)` and `sum(diffs)` match;
  `Check_Trace_Points` passes. Accept near-exact (path ties allowed) as long as the
  decoded alignment identity matches.
- **Exit criterion:** ≥99% of tasks produce a `Check_Trace_Points`-valid vector whose
  decoded identity equals FastGA's within rounding.

### A1 — DONE (2026-07-13) ✓
`gpu/trace_validate.cu` (`make trace_validate`) + `gpu/extract_trace.c` (`make extract_trace`,
dumps FastGA's `Compute_Alignment(DIFF_TRACE)` reference per alignment). The trace kernel:
one warp per task, banded furthest-reaching wave storing every d-layer in global scratch,
lane-0 traceback binning (diff, Δb) into global-`tspace`-phased panels. EXAMPLE, 20,000
real alignment tasks:

| metric | result |
|---|---|
| Check_Trace_Points-VALID (panel count + Σdiffs=edit-dist + ΣΔb=bw) | **100.00%** (19,999/19,999; 1 band-overflow) |
| total edit distance == FastGA | **100.00%** |
| bit-exact per-panel vs FastGA | 54.9% |

**Interpretation:** every GPU trace is a valid FastGA trace-point vector with the *exact same
edit distance* → identical alignment identity, zero quality loss. The 45% that aren't
bit-exact are equally-optimal alternative paths (the furthest-reaching wave places indels
rightmost within a run; FastGA left-aligns). This is cosmetic for the goal: `.1aln`
coverage/#alignments/identity/dedup depend only on endpoints + diffs (both 100% correct),
not on internal gap placement. Tie-break reordering doesn't help (54.6%→54.9%); true
left-alignment would need an indel-shift pass — deferred unless A4 shows dedup sensitivity.

### A2 — Batched library API
- `gpu_trace_batch(g, n, ab,ae,bb,be, out_trace, out_tlen, out_off)` in
  `fastga_gpu.{h,cu}`: one warp/alignment, variable-length outputs via offset array.
- Self-test mirrors A1 but through the library.

### A2 — DONE (2026-07-13) ✓
`gpu_trace_batch(g, n, ab,ae,bb,be, tspace, out_trace, out_tlen)` in `fastga_gpu.{h,cu}`:
one warp per alignment (endpoints index the resident contigs — no per-alignment memcpy),
banded wave in a chunked global scratch (`TR_CHUNK=256` warps, ~1 GB), fixed
`FGA_TRACE_MAX_PAIRS=512`-uint16 output slots. `gpu/trace_lib_test.cu`
(`make trace_lib_test`) packs all tasks into one resident buffer (phase-aligning each
A-offset to its `abpos mod tspace`) and validates the batched path: EXAMPLE 20k tasks →
100% tlen match, 54.88% bit-exact — identical to the A1 kernel, confirming a faithful port.

### A3 — Integrate into `gpu_align_tube`
- Replace `Compute_Alignment` with: GPU endpoints + GPU trace-points; CPU only
  `Compress_TraceTo8` + write `.las`. Keep the CPU `Local_Alignment` fallback for
  guard-rejected / overflow tubes.
- **Two-pass batching** in `align_contigs`: pass 1 collects all tube seeds for the
  contig-pair; one `gpu_discover_batch` + one `gpu_trace_batch`; pass 2 packs results.
  (Removes the one-tube-per-launch stall.)

### A3a — DONE (2026-07-13): correct, but confirms batching is the whole ballgame ✓/✗
`gpu_align_tube` now emits the trace on the GPU (`gpu_trace_batch`, n=1 per tube) instead of
CPU `Compute_Alignment`; a per-thread `pair->gpu_tracebuf` holds the vector, `path->trace`
points at it, `Compress_TraceTo8` packs it. CPU `Local_Alignment` fallback on overflow.

EXAMPLE, T=8:

| version | records | aligned A-bases | sort+align wall |
|---|---:|---:|---:|
| CPU baseline            | 323,569 | 632,119,471 | **16.7 s** |
| −G `Compute_Alignment`  | 320,433 | 626,012,667 | 44.5 min |
| −G A3a (GPU trace)      | 320,491 (**99.05%**) | 626,798,190 (**99.16%**) | **63 min** |

**Correctness: ✓.** GPU-emitted trace-points produce a valid `.1aln` with 99% of alignments
and bases — the trace kernel works in the real pipeline, no quality loss.

**Speed: ✗ (as predicted).** A3a is *slower* than the Compute_Alignment version: it does TWO
per-tube GPU launches (discover n=1 + trace n=1) with per-launch memcpy+sync overhead, and
the trace kernel stores a multi-MB wave per tube. Un-batched GPU per tube loses to the CPU
by ~230×. This is not a kernel problem — it is a **launch-granularity** problem.

### The batching blocker (the real finding)
Batching is blocked by a **sequential data dependency in FastGA's tube walk**: inside a chain,
each alignment's end anti-diagonal `eant` sets the next tube's start (`alow = eant`, FastGA.c
~3494). So tubes within a chain cannot be launched together — tube *i+1* needs tube *i*'s
result. Only the FIRST tube of each chain, and different chains / contig-pairs, are
independent. Real batching therefore needs a **redesign of `search_seeds`**, not a local edit:
  1. Collect the first tube of every chain across a whole seed-sort *part* (many contig-pairs).
  2. `gpu_discover_batch` + `gpu_trace_batch` them all at once (thousands of warps — the
     regime where A1/A2 are fast and the 16.7× design applies).
  3. Iterate: chains whose alignment didn't reach `ahgh` re-enter a second batched round
     (most resolve in one). Preserve the inline redundancy-removal in the write-back pass.
This is a multi-day restructure and is where the actual end-to-end speedup lives.

### A4 — End-to-end timing + quality (pending the batching redesign)
- `-G` vs CPU T=32 wall on EXAMPLE (and a larger pair). Target: align phase faster than
  CPU. Report honestly if PCIe/H2D or fallback rate caps the speedup (Amdahl).
- `.1aln` coverage/identity vs CPU (expect ≥99%, as the correctness-first version).

## Risks
- **Trace semantics fidelity:** path ties vs FastGA's GREEDIEST tie-break may shift a diff
  between adjacent panels; acceptable if `Check_Trace_Points` passes and identity matches.
- **Memory / batch granularity:** wave store per alignment × batch; cap and loop.
- **Giant contigs:** H2D per contig-pair amortized across its tubes (already done in A-load).
- **Fallback rate on divergent tubes:** the >10% divergence tail (Validation B) still
  falls back to CPU; that's correct but caps the speedup on divergent data.
