# FastGA acceleration: a measurement-grounded thesis

## Whole-pipeline streaming dataflow, not a lone kernel port

- **Date:** 2026-07-12
- **Branch:** `agent-optimization-wt`
- **Status:** design/vision doc, grounded in kernel-level measurements (M1–M4 below)
- **Companion:** related-work & publishability survey in [`gpu_fpga_opportunity.md`](gpu_fpga_opportunity.md);
  workload characterization (divergence axis) in [`../benchmark_matrix/`](../benchmark_matrix/)
- **Raw data:** `docs/acceleration/kernel_profile/{human,chimp,mouse}.log` (from a
  compile-guarded profiling build, `-DPROFILE_STAGES` in `FastGA.c` + `-DPROFILE_KERNEL` in `align.c`;
  baseline binary unaffected).

---

## 0. Thesis in one paragraph

FastGA is **not** "a fast tool with one hot alignment kernel." It is a **cache-coherent streaming
sort-merge pipeline**: adaptamer-merge → radix-sort seeds → wave-align → sort+merge output. Prior
GPU/FPGA accelerators each grabbed one compute kernel out of a fundamentally *random-access* tool
(LASTZ SW-extend, minimap2 chaining). We measured FastGA's actual kernels and found: (a) the
**accelerable fraction shifts with divergence** — the wave is 98% of sort+align for similar genomes
but only 56% for divergent ones, where the memory-bandwidth-bound seed sort becomes co-dominant
(M1); (b) the wave's **band width is structurally bounded (<256, never >512) across all divergence
levels** — a fixed systolic array / warp covers ~100% of alignments (M2); (c) alignment length and
wave depth are **heavy-tailed** (M2, M3) — the load-imbalance source. Together these say the
differentiated opportunity is **a resident streaming dataflow accelerator that pipelines the whole
sort→merge→align→sort→merge chain on-device (HBM), exploiting FastGA's cache-coherence in
hardware** — accelerating the *workload's architecture*, not one kernel.

---

## 1. The measurements

All on the divergence axis (CHM13 × query, uppercased, T=32; queries by increasing distance).
Wave metrics are recorded per primary forward-wave call in `Local_Alignment` (`align.c`).

### M1 — the accelerable fraction shifts with divergence

`sort+align` (the 80–90% phase) split into Phase 2 (`reimport`+`rmsd_sort`, the seed radix sort,
memory-bandwidth-bound), Phase 3 (`search_seeds` = chain + Myers wave, compute-bound), Phase 4
(`la_sort`+`la_merge`, output):

| query (≈divergence) | sort+align total | **Phase 3 wave** | Phase 2 seed-sort | Phase 4 output | # wave calls |
|---|--:|--:|--:|--:|--:|
| human (~0.1%) | 606 s | **592 s — 97.7%** | 11 s — 1.8% | 2.9 s — 0.5% | 1.42 M |
| chimp (~1.2%) | 1062 s | **1045 s — 98.4%** | 14 s — 1.3% | 2.8 s — 0.3% | 1.66 M |
| mouse (~15%) | 34 s | **19 s — 56.1%** | **15 s — 43.5%** | 0.1 s — 0.4% | 0.38 M |

**Reading:** for align-heavy pairs the wave is ~98% of sort+align → near-ideal Amdahl for a wave
accelerator. For divergent mouse the wave is only 56%; the **seed radix sort (Phase 2) jumps to 43%**
because there is little to align but all seeds must still be sorted. Output (Phase 4) is always
negligible (<0.5%). *The bottleneck within sort+align moves* — so a lone wave accelerator is
sufficient only in the align-heavy corner.

### M2a — band width is bounded AND divergence-invariant (the feasibility unlock)

Distribution of each alignment's max band width `hgh−low` (log2 bins):

| band | human | chimp | mouse |
|---|--:|--:|--:|
| [64,128) | 70.2% | 68.3% | 53.2% |
| [128,256) | 24.1% | 31.5% | 45.9% |
| [256,512) | 0.07% | 0.07% | 0.03% |
| **≥512** | **0%** | **0%** | **0%** |

**~99.9% of alignments have band < 256, and none exceed 512, at every divergence level.** The band
does *not* grow with divergence. Root cause is architectural: FastGA chains inside a **128-wide
diagonal band** (`BUCK_ANTI` in `FastGA.c`), and x-drop pruning (`TRIM_MLAG` in `align.c`) caps
growth. **Consequence:** a fixed **256-wide** (512 for safety) systolic PE array (FPGA) or a
warp/small-block (GPU) covers essentially every alignment in all regimes — the "adaptive band"
concern is empirically dead.

### M2b — wave depth (differences) is what grows with divergence

Distribution of differences `dif` per alignment (the wave's cycle count / work):

- **human, mouse:** peak at 128–512 differences; mouse tail is light (essentially none >2048).
- **chimp:** much heavier tail — **~12% of alignments have >4096 differences, out to ~10^6.**

This **quantitatively explains chimp's 1062 s outlier** (from the divergence axis): same band width
as human, but far **deeper** waves (more cycles per alignment). Hardware picture: **fixed PE width
(~256), variable cycle count** — and the depth variance is precisely the load-imbalance source.

### M3 — alignment length is heavy-tailed (the size skew)

Sample (first thousands of alignments per dataset; **start-of-file / chr1-biased — indicative only**):

| | median A-len | p99 | max | median div% |
|---|--:|--:|--:|--:|
| human | 1.9 kb | 64 kb | 446 kb | ~7% |
| chimp | 1.8 kb | 22 kb | 42 kb | ~10% |
| mouse | 0.53 kb | 3.3 kb | 11 kb | ~21% |

Most alignments are **short (hundreds–thousands bp)** with a **tail to hundreds of kb** (chromosome-
scale segments). This is the tube-size skew: the bulk is many small tasks; a few giants dominate
wall time — the same imbalance the CPU shows (chimp's align uses ~3/32 cores).

### M4 — arithmetic intensity / roofline (analytical; `perf` unavailable on this node)

- **Wave (Phase 3):** per wave-step, O(band) integer max/add + the slide's 2-bit comparisons;
  working set = `V/M/T` arrays (band × ~16 B ≈ 4 KB at band 256) + a small sequence window — **hot
  in registers/L1**. It is **not** DRAM-bandwidth-bound; it is latency/dependency-bound (the wave is
  sequential in `d`) and **throughput-bound across many independent alignments**. It is O(nd)
  **sparse** — *not* a dense-SW GCUPS workload, so GCUPS comparisons to SW accelerators are
  apples-to-oranges. → wins come from **parallel task throughput** (GPU many warps; FPGA many lanes),
  not FLOPs.
- **Seed sort (Phase 2) & GIX-build:** MSD radix sort — **streaming, memory-bandwidth-bound**
  (multiple passes over the seed/k-mer array). On CPU it is bandwidth-starved (caps at ~7–8 cores,
  ~200–400 GB/s). On an **HBM device (~3 TB/s) it has ~10–15× more bandwidth** → the memory-bound
  stages are accelerable *when resident on-device* (mature GPU radix sort, e.g. CUB, runs near HBM
  bandwidth; a 3 Gbp index ~33 GB fits in 80 GB HBM).

---

## 2. What the measurements imply

### For the kernel (both platforms)
- **Band bounded → fixed hardware.** A 256-wide PE array (FPGA) / warp (GPU) is enough; no dynamic
  band-window machinery. This is the single biggest de-risking result.
- **Depth variance → the load-balance problem.** Alignment cost = (fixed width) × (variable cycles),
  and cycles span 128 → 10^6. The few deep/long alignments (chimp tail; the chr-scale tubes in M3)
  are the load-imbalance case that must be handled explicitly.

### FPGA
- **Fixed 256-wide systolic wavefront array** (PE_k holds `V[k]`; recurrence `V'[k]=max(V[k+1]+1,
  V[k-1]+1, V[k]+2)` is nearest-neighbor; slide = 2-bit comparator; `BVEC` trace = shift register).
- **Trace-point tiling** at the 100-bp panel boundaries bounds on-chip state to O(band×100) —
  breaking GeneTEK's ≤1000-bp wall; FastGA's own panel structure gives the tile cuts.
- Deep tail: giants stream through more tiles (natural); no divergence penalty (each lane
  independent). **Headline: energy / throughput-per-watt** (GeneTEK precedent ~100× energy).

### GPU
- **Warp-per-alignment** (band <256 fits) + **block-cooperation for the deep/long tail** (chimp) +
  **persistent-kernel work-queue** to absorb the skew. 2-bit genomes (~3 GB) resident; trace-point
  recompute for traceback. **Headline: throughput** (esp. all-vs-all); single-run ~5× (Amdahl).

### Whole pipeline (the differentiator)
M1 proves a lone wave accelerator covers 98% of sort+align for similar genomes but only 56% for
divergent ones. M4 shows the memory-bound sort is accelerable **on HBM**. Therefore the strong design
is a **resident streaming dataflow accelerator**: keep the 2-bit genomes + indices in HBM and pipeline
**seed-sort → adaptamer-merge → wave → output-sort/merge** on-device, data never returning to host
memory — the hardware realization of FastGA's cache-coherence. This:
1. covers **all regimes** (wave for align-heavy, HBM-bandwidth sort for divergent);
2. accelerates FastGA-**specific** primitives no one has touched (**adaptamer-merge**: sorted-index
   linear merge with LCP — a merge-network/comparator design);
3. composes with the **storage co-design** (early GIX deletion, chunked build) → many concurrent
   comparisons per node → **node-level all-vs-all throughput**, which **dodges the single-run Amdahl
   ceiling** and rides the pangenome/T2T trend.

---

## 3. Staged roadmap (each stage is a standalone result)

| Stage | Deliverable | De-risked by |
|---|---|---|
| **0 — Characterization** *(done)* | Divergence-axis regime map + M1–M4 kernel profile | this doc + `benchmark_matrix/` |
| **1 — Wave accelerator** | GPU warp/block feasibility **or** FPGA 256-wide systolic + trace-point tiling | M2 band<256; M2b depth model |
| **2 — On-device seed sort** | Radix sort + adaptamer-merge resident on HBM | M1 (43% for mouse); M4 (HBM 10–15×) |
| **3 — Resident all-vs-all** | Indices resident, reused across N² comparisons + storage co-design | storage work; 2-bit compactness |

Recommended: **Stage 1 on GPU** for a fast feasibility/upper-bound number, **Stage 1 on FPGA** as the
novel contribution (trace-point tiling), then **Stage 2–3** for the differentiated whole-pipeline
throughput story.

---

## 4. Honest risk register

- **Amdahl is regime-dependent (M1).** Wave-only acceleration: ~5–6× ceiling for align-heavy, far
  less for divergent (where sort dominates). Mitigation = Stage 2 (accelerate the sort too) + the
  throughput/all-vs-all framing.
- **The deep/long tail (M2b, M3)** is the load-imbalance risk on *both* platforms; it is the core
  engineering problem, not an afterthought.
- **PCIe / host↔device** for a single comparison; amortized only in all-vs-all with resident data.
- **Band-bounded is measured on primate+mouse.** Confirm on a **highly repetitive / large plant or
  amphibian genome** (maize, newt) where chaining density and band behavior may differ before
  hard-claiming "band <256 universally."
- **M3 identity/length is start-of-file sampled** (chr1-biased); re-measure with a uniform pass
  (keep the `.1aln`, `ALNshow` all records) before publishing the length distribution.
- **Engineering cost**: the full resident dataflow is a large system (Darwin/ASIC-scale); the staged
  roadmap keeps each step independently publishable.

---

## 5. Venues & framing

- **Characterization (Stage 0)** → **IISWC** (precedents: PangenomicsBench IISWC'25, GenomicsBench):
  "FastGA is a new cache-coherent WGA workload; here is where it bottlenecks and how that moves with
  divergence." M1's shifting fraction + M2's bounded band are the headline figures.
- **Co-design (Stages 1–3)** → **SC / IPDPS / ICS**; **ISCA/MICRO/HPCA** if the cache-coherence↔
  streaming-dataflow correspondence is argued as an architectural insight; **FPGA/FCCM/FPL / DAC**
  for the FPGA aligner (energy headline).

## 6. Next steps

- [ ] Re-run one dataset keeping the `.1aln` and dumping **all** alignment records → unbiased M3
      (length distribution) + a slide-length sample.
- [ ] Add a **repeat-rich / large genome** (maize or newt) to M2 to test the band bound universally.
- [ ] Prototype Stage 1 GPU wave kernel (batched, band≤256) to get a measured upper-bound speedup.
- [ ] Model the HBM roofline for Stage 2 (radix sort at 3 TB/s) to quantify the divergent-regime win.
