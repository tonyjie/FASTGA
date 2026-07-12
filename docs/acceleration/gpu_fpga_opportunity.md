# GPU/FPGA acceleration of FastGA — related work & publishability

- **Date:** 2026-07-11
- **Branch:** `agent-optimization-wt`
- **Provenance:** deep-research survey (22 sources fetched, 25 claims 3-vote adversarially
  verified → 22 confirmed / 3 refuted). Confidence flagged per claim below. This is a
  literature/strategy scan, not a substitute for a full related-work section — verify the two
  weakly-sourced figures and run the targeted novelty search (§7) before writing hard claims.

## 1. TL;DR verdict

There **is** a credible publishable opportunity. Every existing GPU/FPGA genome-alignment
accelerator targets a **different kernel/algorithm** than FastGA's aligner. FastGA's actual
aligner — **Myers' 2014 O(nd) *local trace-point* difference-wave (the DALIGNER LA-finder)** at
whole-genome scale — appears to be an **open, unaddressed acceleration target** (medium
confidence; §3). A bare kernel port would be Amdahl-limited and incremental; the **strong** angle
is a **workload-characterization + full-system co-design** paper that (a) characterizes FastGA as
a new important cache-coherent WGA workload, (b) fixes its sort+align **load imbalance** (only
~3–6 of 32 cores), (c) accelerates the Myers local-wave kernel, and (d) folds in our
**storage/memory co-design** for concurrent large-genome runs.

## 2. Related work (verified)

Every row accelerates a kernel/algorithm **distinct** from FastGA's unit-cost O(nd) local
trace-point wave.

| Accelerator | Accelerates | Platform | Reported speedup | Venue / year |
|---|---|---|---|---|
| **SegAlign** | LASTZ seed-filter-extend (ungapped X-drop dominates) | 8×V100 GPU | up to **14×** vs 96-vCPU parallel LASTZ; 180× seed-filter vs 1 core | SC '20 |
| **FastZ** | gapped LASTZ seed-and-extend SW-DP (y-drop pruning) | GPU (Titan X/V100/RTX3080) | **43× / 93× / 111×** vs sequential LASTZ | SC '21 |
| **Darwin / GACT** | D-SOFT seed filtration + **tiled Smith-Waterman** DP (systolic anti-diagonal) | FPGA/ASIC | ~**15,000×** (assembly, full-system); 15.6× over DALIGNER | ASPLOS '18 (Best Paper) |
| **GPU-Darwin (GACT)** | confirms GACT = tiled SW on GPU | GPU | — | BMC Bioinformatics '20 |
| **KSW2-GPU** (Zeni et al.) | Suzuki-Kasahara **gap-affine** (Gotoh) score-difference DP | H100 GPU | **1145 GCUPS**; up to 8.51× in minimap2 | PACT '24 |
| **WFA-GPU** | Marco-Sola **gap-affine** Wavefront Alignment | GPU | mature; gap-affine only | Bioinformatics '23 |
| **eWFA-GPU** | **global** unit-cost edit-distance WFA (CIGAR) | GPU | **19–100×** over an O(ND) CPU impl; 24–102× over BPM | IEEE Access '22 |
| **Guo et al.** | minimap2 **chaining** (1-D DP, ~70% of runtime) | FPGA + GPU | first to target chaining (as of '19) | FCCM '19 |
| **mm2-ax** (Sadasivan) | minimap2 **chaining** only (heterogeneous CPU/GPU) | A100 GPU | **2.57–5.41×** vs mm2-fast (30 AVX-512 cores) | 2023 |
| **GeneTEK** | Myers **exact unit-cost edit distance**, **≤1000 bp** reads (HLS FIFO-capped) | Zynq UltraScale+ FPGA | ~**2.1×** speed; up to **111× energy** | 2025 (arXiv 2509.01020 / CBM) |
| **FPGA Myers bit-vector** | Myers' **1999** bit-parallel edit distance (different algorithm) | FPGA | — | 2016 (MEDICON) |

**Refuted during verification** (do not cite as-is): (a) that GeneTEK is "the closest precedent
for FastGA's aligner" — 0-3, it is read-scale (≤1000 bp), not genome-scale trace-point;
(b) that WFA-GPU is "long-read pairwise only, not WGA" — 1-2, scope overstated.

## 3. Novelty gap — is FastGA's aligner an open target?

**Verdict: appears open (medium confidence).** The survey exhaustively maps existing
accelerators onto *other* algorithms: LASTZ seed-filter-extend (SegAlign, FastZ), tiled/gap-affine
SW DP (Darwin/GACT, KSW2-GPU), gap-affine or global WFA (WFA-GPU, eWFA-GPU), Myers **bit-vector**
1999 and **read-scale** unit-cost edit distance (GeneTEK), and chaining (Guo, mm2-ax). **None**
accelerate Myers' 2014 O(nd) *local trace-point* wave at genome scale.

Two data points shape feasibility:
- **eWFA-GPU** proves the unit-cost / O(ND) kernel family **is GPU-tractable** (19–100× over CPU
  O(ND)) — de-risks the GPU path, while leaving the *local trace-point* variant open.
- **GeneTEK** shows FPGA unit-cost edit distance yields large **energy** wins but only ~2× raw
  speed, and is capped at ≤1000 bp — a caution for the FPGA path at genome scale.

⚠️ **This is an inductive "absence of evidence" claim**, not a sourced negative. Before any hard
novelty statement in a paper, run §7's targeted search.

## 4. Publishability & venues

- **Bare kernel port = incremental.** Amdahl ceiling: seed-merge ~1% + index build ~15% ⇒ even
  infinite align speedup ≈ **5–6× end-to-end**. Systems reviewers scrutinize exactly this.
- **What elevates it to strong:** (1) characterize FastGA as a new, important cache-coherent WGA
  workload; (2) attack the true **80–90% sort+align** bottleneck end-to-end; (3) **fix the load
  imbalance** so the hardware is actually used; (4) fold in the **storage/memory co-design**
  (lower peak disk → concurrent large-genome / all-vs-all runs); (5) report real whole-genome
  speedup **and energy/cost**, not microbenchmark GCUPS.
- **Venue fit:**
  - **IISWC** — the workload-characterization piece. Direct precedents: *PangenomicsBench*
    (IISWC '25, Cornell/Batten) and *GenomicsBench* (UMich). Our divergence-axis data + the
    load-imbalance finding are already ~60% of such a paper.
  - **SC / IPDPS / ICS** — full-system co-design (parallelization + acceleration + storage).
  - **ISCA / MICRO / HPCA** — only with a genuine architectural insight (the cache-coherence vs
    random-access memory tradeoff; the memory-BW-bound index build).
  - **FPGA / FCCM / FPL / DAC** — a kernel-focused FPGA aligner.
  - **Bioinformatics / BIBM / TCBB** — tool/impact framing.

## 5. Precedents (accelerating an already-fast CPU tool)

The **minimap2** line is the closest precedent and validates the pattern: **profile first,
accelerate the real hotspot, split heterogeneously mindful of Amdahl/PCIe.** Guo et al. (FCCM '19)
and **mm2-ax** (2023) offload only **chaining** (~60–70% of runtime; ~68% for ONT >100 kb) to
FPGA/GPU, leaving seeding/base-alignment on CPU — 2.6–5.4× end-to-end with cost-awareness. This is
almost exactly FastGA's situation with sort+align. Pitfalls they navigated and reviewers expect
addressed: **Amdahl** (offload must target the dominant phase), **PCIe transfer**, and
**load balance / occupancy**.

## 6. Trends (2023–2025) that help the framing

- **Pangenome / T2T all-vs-all** is exploding (94-haplotype draft human pangenome; Nature '23),
  driving demand for fast, storage-frugal WGA at scale.
- **wfmash** is the standard all-vs-all engine for pangenome pipelines (nf-core/pangenome) — and
  it uses **gap-affine WFA**, a *different* algorithm; FastGA's cache-coherent adaptamer approach
  is a distinct point in the design space.
- **NVIDIA Parabricks v4.4** GPU-accelerates **Giraffe** (pangenome *graph*, short-read-to-graph)
  — **not** whole-genome assembly alignment. It does **not** preempt FastGA's niche.
- **Mumemto** (2025) — a CPU pangenome MUM tool (prefix-free parsing), algorithmic not hardware.

## 7. The pivotal open question (decide before committing)

**Can the atomic-contig load imbalance be fixed on the CPU alone** — finer-grained work
decomposition (per-diagonal-band / per-chain tasks) — without breaking FastGA's cache-coherence,
recovering much of the 32 cores **before any accelerator?**
- If **yes**: a strong, cheap first contribution — but it *weakens* the "you need a GPU/FPGA"
  argument, so it must be measured honestly.
- Either way, the honest sequencing is the paper's backbone: **fix load balance on CPU → measure
  the new ceiling → quantify remaining headroom for GPU/FPGA → decide speedup vs energy/cost vs
  concurrency as the headline.**

Relevant code to assess this: `align_contigs()` and `search_seeds()` in `FastGA.c`, the reverse
MSD sort partition boundaries in `RSDsort.c`, and the per-thread contig assignment
(`IDBsplit`/`Select`).

## 8. Recommended angle

A **two-paper arc**: (1) **IISWC** — workload characterization of FastGA (regime map from our
profiling matrix + the load-imbalance analysis); (2) **SC / IPDPS** — the co-design system
(CPU load-balance fix → Myers local-wave GPU/FPGA offload → storage co-design), reporting
whole-genome speedup and energy/cost. If aiming for a single paper, target SC/IPDPS with the
characterization as motivation.

## 9. Caveats on this survey

- Novelty verdict is **medium** (inductive), not sourced-negative.
- **GACT "80,000×" could not be confirmed** — canonical Darwin numbers are 15,000× full-system,
  15.6× over DALIGNER; use those.
- **mm2-ax** figures come via a low-tier publisher (Fortune Journals), though mirrored on
  bioRxiv/PMC.
- **GeneTEK "113%" = ~2.13×**, not 113×.
- The "no prior chaining acceleration" precedent claim is **time-bound to ~2019**; chaining is now
  an active target (mm2-gb 2024, minimap2-fpga 2023) — don't present as current.
- All venue/publishability guidance is analyst judgment extrapolated from precedent.

## 10. Next steps

- [ ] **Targeted novelty search** to harden §3: `daligner GPU/FPGA`, `trace-point alignment
      hardware`, `LA-finder acceleration`, and forward-citations of Myers WABI 2014.
- [ ] **CPU load-imbalance study** (§7) — the pivotal question.
- [ ] Confirm the chimp align-time outlier with ≥3 reps; extend the profiling matrix (size axis,
      pig) to strengthen the characterization.

## Sources

- SegAlign (SC'20): https://public.gi.ucsc.edu/~yatisht/files/segalign-sc20.pdf
- FastZ (SC'21): https://ieeexplore.ieee.org/document/9910115/
- Darwin/GACT (ASPLOS'18): http://bejerano.stanford.edu/papers/p199-turakhia.pdf
- GPU-Darwin (BMC Bioinf'20): https://link.springer.com/article/10.1186/s12859-020-03685-1
- KSW2-GPU / difference-recurrence (PACT'24): https://dl.acm.org/doi/10.1145/3656019.3676894
- WFA-GPU (Bioinformatics'23): https://academic.oup.com/bioinformatics/article/39/12/btad701/7425447 · https://pmc.ncbi.nlm.nih.gov/articles/PMC10697739/
- eWFA-GPU (IEEE Access'22): https://www.researchgate.net/publication/361280209_Accelerating_Edit-Distance_Sequence_Alignment_on_GPU_Using_the_Wavefront_Algorithm
- GeneTEK (FPGA, 2025): https://arxiv.org/abs/2509.01020
- FPGA Myers bit-vector (2016): https://link.springer.com/chapter/10.1007/978-3-319-32703-7_104
- minimap2 chaining accel — Guo (FCCM'19): https://vast.cs.ucla.edu/sites/default/files/publications/minimap2-acc-approved.pdf
- mm2-ax (2023): https://www.fortunejournals.com/articles/accelerating-minimap2-for-accurate-long-read-alignment-on-gpus.html
- PangenomicsBench (IISWC'25): https://www.csl.cornell.edu/~cbatten/pdfs/kaplan-panobench-iiswc2025.pdf
- GenomicsBench: https://genomicsbench.eecs.umich.edu/
- NVIDIA Parabricks pangenome: https://developer.nvidia.com/blog/discover-new-biological-insights-with-accelerated-pangenome-alignment-in-nvidia-parabricks/
- wfmash: https://github.com/waveygang/wfmash
- Human pangenome reference (Nature'23): https://www.nature.com/articles/s41586-023-05896-x
