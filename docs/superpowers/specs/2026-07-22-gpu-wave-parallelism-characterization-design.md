# GPU Wave-Parallelism Characterization — Design Spec

**Status:** draft (brainstormed 2026-07-22)
**Author:** jl4257 + agent (worktree `agent-opt`, branch `agent-optimization-wt`)
**Node:** en-ec-zhang-x4 — 2× AMD EPYC 9124 (32 physical cores / 64 SMT), 2× NVIDIA A100 80GB
**Supersedes framing of:** the prior `gpu/` kernels (discover + trace two-kernel design). Those are
kept only as a *coding reference*; their conclusions are not carried forward (see §9).

---

## 1. Goal (one sentence)

Measure, honestly and on the real workload, **how much of the parallelism that already exists in
FastGA's `Local_Alignment` a GPU can actually exploit, and what speedup that yields** versus the
improved CPU baseline — via a standalone, faithful warp-cooperative port and a benchmark harness.
This is a **characterization study**, not a production integration.

## 2. First-principles grounding (the ground truth this study rests on)

### 2.1 What the CPU baseline does (known, verified from `align.c` this session)
`Local_Alignment(align, work, spec, low, hgh, anti, lbord, hbord)` is FastGA's real per-alignment
aligner (the align-phase hot path; **not** `Compute_Alignment`, which FastGA never calls). For one
alignment it runs, from a seed anti-diagonal/diagonal:
- `reverse_wave` then `forward_wave` — each **one** x-drop-terminated furthest-reaching
  (Myers O(nd), unit-cost) wavefront;
- it records **sparse checkpoints** along the sweep and **never stores the whole DP matrix**;
- trace-points come from a **cheap backward pointer-walk** over those checkpoints;
- x-drop / trim (`TRIM_MLAG` accounting) is the stop condition. Memory is O(sparse), independent of
  depth.

### 2.2 The parallelism that already exists inside it
- **Intra-alignment (band):** at a fixed depth `d`, the ~128 active band diagonals are mutually
  independent and can be computed in parallel. **Depth `d` itself is sequential** — `d` depends on
  `d−1`. We do **not** break this.
- **Inter-alignment:** the ~812K alignments (`forward_wave` calls in the human profiling) are
  mutually independent.

### 2.3 How the CPU harvests it (known)
Only **inter-alignment**, across 28–32 cores (workstealing-v3c: a contig-pair work-pool, 5.33×
align-phase, bit-exact). The CPU does **not** touch the intra-alignment band parallelism — one core
sweeps a band scalarly.

### 2.4 The measured workload distribution (from `assets/fastga-profiling/results/C2_bandhist.txt`)
- Band ≈ 128 (96% in [64,256]) — near-constant, so per-wave work ≈ band × depth ≈ depth-driven.
- Depth heavy-tailed: 812,662 waves; **88% of work in the 5.3% with depth ≥ 8K**, 58% in the 1.4%
  with depth ≥ 64K. **No single monster wave** (single-wave util ceiling 59–237×) — the tail is
  many *independent* deep waves. Favorable for inter-alignment parallelism.

### 2.5 What we want
Map **both** existing parallelism levels onto the GPU — band → warp lanes, alignments → warps —
**without changing the algorithm or breaking depth-sequentiality** — and measure the real speedup,
understanding *where* (in the distribution) the GPU wins or loses.

## 3. Non-goals / integrity guardrails (binding)

1. **No FastGA integration.** No `-G` path, no host-loop/phase-split changes. That is a separate,
   later spec, justified only if this study's measured fit warrants it.
2. **No synthetic or GPU-flattering datasets.** The **full real distribution** is the headline
   number. A GPU-friendly subset (e.g. the deep tail) may appear **only** as a labelled diagnostic,
   never as the headline. Any sampling of the full set must be **uniform / distribution-preserving**
   and explicitly distinguished from a biased subset.
3. **Do not break depth-sequentiality or invent new parallelism.** Only harvest the band +
   inter-alignment parallelism that already exists.
4. **Faithful algorithm.** The kernel replicates `align.c`'s actual wave (its real x-drop/trim), not
   a reinvented score criterion. (A reinvented `matches−2·diffs`, x-drop-40 criterion produced only
   ~60% endpoint fidelity in a throwaway bench — that is a warning, not the target.)

## 4. Approach — the mapping (decided)

**Warp-per-alignment (Approach 1).** One warp = one alignment; its 32 lanes tile the ~128 band
diagonals at each depth step; the depth loop is sequential within the warp; thousands of alignments
run on thousands of warps. This is the literal hardware image of §2.2.

Rejected: **thread-per-alignment** (discards band parallelism, huge per-thread band state, each
thread very slow); **block-per-alignment** (only needed for band ≫ 32×; band~128 fits a warp at
4 diagonals/lane; heavier sync).

**Honest mechanism this study will quantify** (from a throwaway measurement this session, to be
re-established faithfully): one warp is *slower* than one CPU core on a single deep wave (~40×), and
the GPU's case rests entirely on having ~100× more concurrent warps than the CPU has cores. The
study must report this decomposition (per-warp rate × concurrent warps vs per-core rate × cores),
not just the aggregate ratio.

## 5. Architecture — three parts

### Part A — Faithful kernel (`gpu/wave_kernel.cu`, new)
A warp-cooperative CUDA port of `align.c`'s **actual** `forward_wave` / `reverse_wave`:
- **Stage 1 (sweep):** forward + reverse x-drop wave → endpoints + diff-count, replicating align.c's
  real trim/x-drop so endpoints match the CPU (target: ≫60% — as high as the parallel order allows,
  see §6). Ping-pong band state (2 layers), int positions (deep waves reach x > 2^15), **no depth
  cap**, band across 32 lanes with a warp argmax where align.c takes a max.
- **Stage 2 (trace):** sparse-checkpoint save during the sweep + backward-walk → FastGA trace-points
  (the same mechanism align.c uses; **not** the old "store all d-layers" trace kernel that overflowed
  at depth 2048). Memory O(checkpoints × band), no depth cap.

**Staging rule (fit-driven):** build + validate Stage 1 first. Proceed to Stage 2 only as far as the
measured Stage-1 fit justifies; the ultimate target is end-to-end output matching the CPU (§6).

**Tail scheduling (sub-design, measured, not assumed):** the depth heavy-tail imbalances warps.
Evaluate naive order, depth-sorted order, and (if needed) a persistent-kernel work-stealing variant.
Report which is used and the imbalance penalty.

### Part B — Benchmark harness
**Genome-resident, seed-driven** (this both scales to the full distribution and is the primary
timing basis — no per-alignment window extraction, no hundreds of GB):
- Load both genomes resident: CPU = NUMERIC contig buffers (as FastGA does); GPU = resident 2-bit or
  NUMERIC genome (~750 MB each).
- Extract, from the human `.1aln`, the **full list of per-alignment seeds** = (contig pair, strand,
  seed anti/diag at the alignment midpoint, contig bounds, reference endpoints+diffs for validation).
  This is small (~812K/518K records × a few ints) — the whole real distribution, streamed, **no
  biased subset**.
  - *Methodology note (documented limitation):* midpoint-of-real-alignment seeds are guaranteed to
    sit inside a real alignment; they characterize the cost of the waves that produce real output
    (what the align phase spends time on), not FastGA's failed/short discoveries (which are cheap).
- **CPU baseline:** isolated `Local_Alignment` on the resident contigs — **1 core** (per-core rate)
  and **32 cores** (real T=32, `OMP schedule(dynamic)`, pinned). The real function, no seed-sort /
  decompress / output confounds.
- **GPU:** Part-A kernel, warp-per-seed, reading the resident genome.

**Reporting:**
- **Headline:** GPU vs CPU on the **full distribution** (expected *below* the deep-tail figure —
  shallow waves idle lanes; that is the honest number).
- **Diagnostic:** stratified by (band, depth) bucket → GPU-vs-CPU per regime → a *where it fits* map
  and a routing statement (e.g. deep→GPU, shallow→CPU).
- **Two timing bases (honest range):** (i) wave-engine — genome resident, kernel compute vs isolated
  CPU; (ii) with per-batch H2D/D2H + launch overhead counted. Report both.

### Part C — Correctness / divergence analyzer
**Tiered, optimality-based bar** (decided): a GPU result is *correct* if it has the **same edit
distance (score)** as the CPU and its trace-points are **`Check_Trace_Points`-valid**. Report:
- % byte-exact endpoints, % same-score, % valid;
- **classify every divergence** (same-score-different-path vs genuinely-worse) and **attribute** it to
  the GPU's parallel band-order / argmax tie-breaking vs the CPU's sequential scan schedule.
Respect the GPU's parallel nature (do not force CPU scan order just to match bytes); the deliverable
is a clear account of *where and why* GPU output differs from the CPU baseline.

## 6. Correctness methodology (precise)

- **Reference:** the human `.1aln` (518,037 non-redundant alignments) and the CPU `Local_Alignment`
  re-run in the harness (endpoints+diffs, and trace-points in Stage 2).
- **Stage 1 gate:** GPU endpoints+diffs vs CPU — report exact-match %, same-score % (|Δdiff|=0),
  and depth-work ratio (GPU diffs / CPU diffs) so a throughput comparison is judged on matched work.
- **Stage 2 gate:** GPU trace-points `Check_Trace_Points`-valid at 100%, same edit distance; report
  byte-exact vs equally-optimal split.
- Divergence is **analyzed, not eliminated** — but every divergence class must be named and
  attributed.

## 7. Measurement methodology (what is counted)

- **Primary (wave-engine):** genome pre-resident on GPU (one-time, uncounted); kernel compute timed
  with CUDA events, best-of-N. CPU = isolated `Local_Alignment` wall (best-of-N, pinned).
- **Secondary (realistic):** add per-batch H2D of seeds + D2H of results + launch. (Genome H2D is a
  one-time residency cost, reported separately, not per-batch.)
- **Decomposition:** per-warp effective rate × concurrent warps (from achieved occupancy) vs
  per-core rate × cores — so the *mechanism* of any win/loss is explicit, not just the ratio.
- Full-distribution headline + per-(band,depth)-bucket diagnostic, both timing bases.

## 8. Success criteria (what the study delivers)

1. **A full-distribution GPU-vs-CPU throughput result (+ honest range)** on the real workload — no
   cherry-picking. This is the headline answer to "can the GPU exploit this parallelism, and how much."
2. **A fit map:** GPU-vs-CPU per (band, depth) regime → a clear verdict on which parts of the real
   workload are a good/poor GPU fit, with the per-warp-vs-per-core mechanism decomposition.
3. **A correctness report:** tiered pass rates + a named, attributed divergence analysis vs the CPU
   schedule.
4. A written conclusion: is this workload a good GPU fit, where, why, and what the measured speedup is.

## 9. What is explicitly dropped from the prior GPU work (and why)

The prior two-kernel design (`disc_batch` + `trace_batch`) computes a **different algorithm** than
`Local_Alignment`: it sweeps once to find endpoints, then **re-runs a full banded DP storing all
d-layers** for traceback (O(band×depth) memory, depth-capped at 2048 → overflow on the deep waves
holding 88% of the work), i.e. **two passes + whole-matrix storage** vs the CPU's **one pass +
sparse checkpoints**. Its throughput numbers (8.1×/306k aln/s) were measured on shallow EXAMPLE data
and on `Compute_Alignment` (not FastGA's hot path). All of that is superseded. This study starts from
the faithful port.

## 10. Risks & open questions

- **Fidelity ceiling under parallel order.** Matching align.c's exact endpoints may be limited by the
  warp's band-order/argmax tie-breaking. Mitigation: port align.c's *actual* trim logic (not a
  reinvented criterion); measure how close the parallel order gets; treat the residual as attributed
  divergence (§6), acceptable per the tiered bar.
- **Tail imbalance on the GPU.** The depth heavy-tail imbalances warps; the fix is scheduling (Part A
  tail sub-design), measured not assumed.
- **Full-distribution scale.** Handled by genome-residency (seeds only, not windows); if the seed set
  is still large, process in streamed batches — never a biased subset.
- **Seed choice.** Midpoint seeds characterize real-alignment wave cost, not failed discoveries;
  documented in §5-B.

## 11. Environment & existing assets

- **Dataset:** `/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/output.1aln` (518,037 alns)
  + the two `.1gdb`/`.bps`. Human run used stock defaults (ALIGN_RATE .3 → `New_Align_Spec(0.7,100,
  freq,0)`, ALIGN_MIN 100).
- **Reuse (coding reference only):** `gpu/extract_deep.c` (seed/window extraction from `.1aln`),
  `gpu/deep_cpu_bench.c` (isolated `Local_Alignment` timing pattern, sentinel-correct buffers),
  `gpu/fastga_gpu.cu` (warp-cooperative writing style). The kernel itself is **rewritten** as a
  faithful align.c port.
- **Source of truth for the wave:** `align.c` `forward_wave` / `reverse_wave` / trim logic — to be
  read and replicated in Part A.
