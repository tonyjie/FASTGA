# GPU Wave-Parallelism Characterization — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Measure how much of `Local_Alignment`'s existing band + inter-alignment parallelism a GPU
can exploit, and the resulting speedup vs the improved CPU baseline, on the *full real* human
workload — via a faithful warp-cooperative port of `align.c`'s wave and a genome-resident benchmark
harness.

**Architecture:** A standalone benchmark (NO FastGA integration). Both engines read from a resident
genome and a compact per-alignment **seed list** extracted from the real `.1aln` (so the full
distribution costs seeds, not hundreds of GB of windows). CPU engine = FastGA's real
`Local_Alignment`. GPU engine = a warp-per-alignment port of `align.c`'s `forward_wave`/`reverse_wave`
(band → 32 lanes, alignments → warps, depth sequential). Built in two stages: Stage 1 sweep
(endpoints+diffs), Stage 2 sparse-checkpoint trace-points. Correctness is tiered/optimality-based
with divergence attributed to parallel order.

**Tech Stack:** C (C11) + CUDA (nvcc, `-arch=sm_80`, A100 80GB); OpenMP for the CPU baseline; FastGA's
own `align.c` / `GDB.c` / `alncode.c` / `gene_core.c` / `ONElib.c` linked in. Node: 2× EPYC 9124
(32 phys cores), 2× A100 80GB.

## Global Constraints

- **Faithful algorithm.** The GPU wave replicates `align.c`'s actual `forward_wave`/`reverse_wave`
  (its real x-drop/trim: `TRIM_MLAG=250`, score/reach from the `Align_Spec`), NOT a reinvented
  criterion. Do not break depth-sequentiality; parallelize only the band (across 32 lanes) and
  across alignments (across warps).
- **CPU baseline = isolated `Local_Alignment`** (the real hot-path function), never `Compute_Alignment`.
  Spec: `New_Align_Spec(1.-ALIGN_RATE, 100, gdb->freq, 0)` with stock `ALIGN_RATE=.3`, `ALIGN_MIN=100`.
- **Sentinel convention:** every sequence buffer has `seq[-1]==4` and `seq[len]==4` (NUMERIC 0–3
  between), exactly as `New_Contig_Buffer`+`Get_Contig` produce; the wave's inner match loop relies
  on the `4` sentinel to stop.
- **Full real distribution is the headline.** Any subset (e.g. deep tail) is a labelled diagnostic
  only. Any down-sampling is uniform/distribution-preserving, never biased toward GPU-friendly waves.
- **Two timing bases, always both:** wave-engine (genome resident, kernel compute vs isolated CPU)
  AND with per-batch H2D/D2H+launch. Report the honest range.
- **Correctness = tiered:** same edit distance (score) + `Check_Trace_Points`-valid; classify &
  attribute every divergence (same-score-different-path vs worse) to parallel order. Do not force CPU
  scan order to chase bytes.
- **Dataset:** `/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/output.1aln` (518,037 alns)
  + the two `.1gdb`/`.bps`. Reference paper spec: `docs/superpowers/specs/2026-07-22-gpu-wave-parallelism-characterization-design.md`.
- **No commits without the human's approval.** Build artifacts (binaries) are not committed.

## File Structure

| File | Responsibility |
|---|---|
| `gpu/WAVE_PORT_NOTES.md` (new) | Porting reference: the exact `align.c` wave contract to replicate (params, x-drop, checkpoints, trace emission) + the warp mapping. |
| `gpu/wave_harness.h` (new) | Shared formats: `SeedRec` (per-alignment seed+bounds+ref), file header, small helpers. |
| `gpu/extract_seeds.c` (new) | Extract the full-distribution seed list from `.1aln` (no windows). |
| `gpu/wave_bench_cpu.c` (new) | Genome-resident CPU baseline: real `Local_Alignment` per seed, 1+32 core, stratified. |
| `gpu/wave_kernel.cu` (new) | Faithful warp-cooperative port: Stage-1 sweep, then Stage-2 trace. Device kernels + C-callable launchers. |
| `gpu/wave_bench_gpu.cu` (new) | Genome-resident GPU driver: loads genome + seeds, runs the kernel, times (both bases), validates vs CPU-produced reference, stratifies. |
| `gpu/wave_validate.c` (new) | Correctness/divergence analyzer: compares GPU vs CPU per-alignment (score/endpoints/trace), classifies + attributes divergence. |
| `Makefile` (modify) | Add targets for the six new binaries. |

Reuse (coding reference; do not import verbatim): `gpu/extract_deep.c`, `gpu/deep_cpu_bench.c`,
`gpu/fastga_gpu.cu`.

---

## Task 1: Porting reference — document align.c's wave contract

**Files:**
- Create: `gpu/WAVE_PORT_NOTES.md`
- Read: `align.c:385-916` (`forward_wave`), `align.c:919-1450` (`reverse_wave`), `align.h:205-236`.

**Interfaces:**
- Produces: a written contract later tasks port against — the exact inputs/outputs, the x-drop rule,
  the band iteration, the `Pebble cells[]` checkpoint layout, and where trace-points are emitted.

- [ ] **Step 1:** Read `forward_wave`/`reverse_wave` and write `WAVE_PORT_NOTES.md` capturing, with
  `align.c` line cites: (a) inputs — `aseq/bseq` (sentinel `4`-terminated), band `[*mind,maxd]`, seed
  anti-diagonal `mida`, `aoff`, `Align_Spec` (`trace_space=100`, `ave_path`, `reach`, `score`,
  `table`); (b) the per-diagonal furthest-reach inner loop (`while c==d` slide, `c==4`/`d==4` clip);
  (c) the x-drop/trim rule (`besta/trima`, `TRIM_MLAG=250`, `while more && lasta >= besta-TRIM_MLAG`);
  (d) `Pebble{ptr,diag,diff,mark}` = the sparse checkpoint at each `tspace` boundary (`mark=na`), and
  the backward walk (`for h=cells[h].ptr; h>=0; ...`) that emits trace-points; (e) outputs written to
  `apath` (`abpos/aepos/bbpos/bepos/diffs/trace/tlen`).
- [ ] **Step 2:** In the same doc, state the warp mapping decision: 32 lanes tile the band diagonals
  `[low,hgh]` (≤ ~256 wide → ≤ 8 diagonals/lane) at each `dif`; the `dif` loop is sequential; a warp
  argmax replaces the scalar best/trim tracking; `cells[]` become a per-warp checkpoint buffer.
  List the 3 known GPU deviations to expect and attribute later: parallel band-order tie-breaking,
  argmax vs scalar first-best, and int-vs-native position width.
- [ ] **Step 3:** Commit (`docs: WAVE_PORT_NOTES — align.c wave contract for the GPU port`).

*No automated test; the deliverable is the reference doc. Acceptance: every claim cites an `align.c`
line and a reviewer can follow the port from it.*

---

## Task 2: Full-distribution seed extractor

**Files:**
- Create: `gpu/wave_harness.h`, `gpu/extract_seeds.c`
- Modify: `Makefile`
- Reference: `gpu/extract_deep.c` (orientation/reading logic), `ALNtoPAF.c` (record sweep).

**Interfaces:**
- Produces: `SeedRec` + a `.seeds` file consumed by Tasks 3 and 6.

```c
// wave_harness.h
#ifndef WAVE_HARNESS_H
#define WAVE_HARNESS_H
#include <stdint.h>
#define WAVE_SEEDS_MAGIC 0x53564157u   // "WAVS"
typedef struct { uint32_t magic, nseeds, tspace, reserved; } WaveSeedHeader;
typedef struct {
  int32_t aread, bread;      // contig indices into gdb1 / gdb2
  int32_t flags;             // COMP bit = B reverse-complemented
  int32_t alen, blen;        // contig lengths
  int32_t seed_anti, seed_diag;   // midpoint seed: anti=ab+ae... , diag=diag(mid); see extractor
  int32_t ref_ab, ref_ae, ref_bb, ref_be, ref_diffs;  // .1aln reference, for validation
} SeedRec;
#endif
```

- [ ] **Step 1 (test first):** Write `gpu/test_extract_seeds.sh` that will run the extractor on the
  human `.1aln` and assert: header `nseeds == 518037`, `tspace == 100`, and that the first record's
  `ref_ab/ref_ae` match `ALNtoPAF` row 1 (extract A-start/end via `ALNtoPAF | head -1`). Run it →
  expect FAIL (no binary yet).
- [ ] **Step 2:** Implement `extract_seeds.c` by adapting `extract_deep.c`: open the `.1aln` + both
  GDBs, sweep every overlap, and for each write a `SeedRec` — `aread/bread/flags/alen/blen` from the
  overlap+GDB, `ref_*` from `path`, and the seed via the midpoint convention verified in
  `deep_cpu_bench`: `sa=(ab+ae)/2`, `sb=(bb+be)/2`, then `seed_diag = sa-sb`, `seed_anti = sa+sb`.
  Emit **no sequence windows**. Patch the header at the end. NO
  `dmin`/`wmax` filtering — the full distribution.
- [ ] **Step 3:** Add `Makefile` target `extract_seeds` (mirror the `extract_disc` target: link
  `align.c GDB.c alncode.c gene_core.c ONElib.c -lpthread -lm -lz`).
- [ ] **Step 4:** Run `test_extract_seeds.sh` → expect PASS (`nseeds==518037`, tspace 100, row-1 match).
- [ ] **Step 5:** Commit (`feat: full-distribution seed extractor for the wave harness`).

---

## Task 3: CPU baseline harness (genome-resident, real Local_Alignment)

**Files:**
- Create: `gpu/wave_bench_cpu.c`
- Modify: `Makefile`
- Reference: `gpu/deep_cpu_bench.c` (Local_Alignment call + sentinel buffers), `extract_deep.c`
  (B-contig orientation via `Get_Contig_Piece`+`Complement_Seq`).

**Interfaces:**
- Consumes: `.seeds` (Task 2), the two GDBs.
- Produces: a reference dump `.cpuref` (per-seed `abpos/aepos/bbpos/bepos/diffs`) consumed by Task 6/7
  validation, plus stratified throughput to stdout.

- [ ] **Step 1 (test first):** Write `gpu/test_wave_cpu.sh`: run the CPU bench on the full `.seeds`;
  assert it processes 518,037 seeds and that endpoint agreement vs `ref_*` within 50 bp is ≥ 80%
  (matches the deep-set 82% we measured; a regression below this means the seed/orientation setup is
  wrong). Run → FAIL (no binary).
- [ ] **Step 2:** Implement `wave_bench_cpu.c`: open both GDBs, load **all contigs resident** (loop
  `Get_Contig` into per-contig NUMERIC buffers with `[-1]`/`[len]` sentinels; complement handled per
  seed via a forward B cache + `Complement_Seq` for COMP seeds, as `extract_deep.c` does). For each
  seed set `aln.aseq/bseq/alen/blen/flags`, `dg=seed_diag`, `anti=seed_anti`, call
  `Local_Alignment(&aln,work,spec,dg,dg,anti,-1,-1)` with `spec=New_Align_Spec(0.7,100,gdb1->freq,0)`.
  Write `.cpuref`. Time single-thread (per-core, `SKIP_SINGLE=1` to bypass) and 32-thread
  (`OMP schedule(dynamic,16)`, pinned).
- [ ] **Step 3:** Add stratification: bucket each seed by (band proxy, depth=`path.diffs`) and print a
  per-bucket aln/s table. (Band proxy: `|(ae-ab)-(be-bb)|` net drift — the value we can compute
  cheaply; note in-code it is a proxy for the true active band.)
- [ ] **Step 4:** Add `Makefile` target `wave_bench_cpu` (like `deep_cpu_bench`, add `-fopenmp`).
- [ ] **Step 5:** Run `test_wave_cpu.sh` → PASS. Record 1-core and 32-core aln/s in the ledger.
- [ ] **Step 6:** Commit (`feat: genome-resident CPU Local_Alignment baseline harness`).

---

## Task 4: GPU genome residency + seed upload (no kernel yet)

**Files:**
- Create: `gpu/wave_kernel.cu` (context + loaders only in this task), `gpu/wave_bench_gpu.cu` (driver
  skeleton).
- Modify: `Makefile`
- Reference: `gpu/fastga_gpu.cu` (`gpu_open`/`gpu_load_seqs`/`gpu_load_seqs_2bit` residency pattern).

**Interfaces:**
- Produces: `wave_ctx*` holding both genomes resident (NUMERIC, contig-base offset table) + seed
  buffers on device. Consumed by Tasks 5,6.

```c
// in wave_kernel.cu — C-callable
typedef struct wave_ctx wave_ctx;
wave_ctx *wave_open(void);
void      wave_close(wave_ctx*);
// upload whole genomes resident (NUMERIC, sentinel-padded per contig), + contig base offsets
void      wave_load_genomes(wave_ctx*, const unsigned char*A, const long*Abase, int nA, long Alen,
                                        const unsigned char*B, const long*Bbase, int nB, long Blen);
```

- [ ] **Step 1 (test first):** In `wave_bench_gpu.cu` add a `--selftest` path that uploads the two
  genomes and copies back a few contig bytes, asserting they match the host NUMERIC. Build target
  `wave_bench_gpu`. Run → FAIL (loaders unimplemented).
- [ ] **Step 2:** Implement `wave_open/close/wave_load_genomes` in `wave_kernel.cu` (cudaMalloc the two
  genome byte arrays + base-offset int64 arrays; H2D copy; store lengths). Report genome H2D time
  separately (one-time residency cost).
- [ ] **Step 3:** Implement the driver skeleton: read `.seeds`, build device seed arrays
  (contig→genome absolute coord via base offsets), run `--selftest`.
- [ ] **Step 4:** Add `Makefile` targets (`nvcc -O3 -arch=sm_80 -Igpu`). Run `--selftest` → PASS.
- [ ] **Step 5:** Commit (`feat: GPU genome residency + seed upload (loaders, no wave yet)`).

---

## Task 5: Stage-1 forward sweep kernel (faithful, warp-cooperative)

**Files:**
- Modify: `gpu/wave_kernel.cu` (add `forward_sweep_warp` + launcher), `gpu/wave_bench_gpu.cu`.
- Reference: `WAVE_PORT_NOTES.md`, `align.c:385-916`.

**Interfaces:**
- Produces: `int wave_forward_batch(wave_ctx*, int n, seeds…, int*ae,int*be,int*fdiff)` — forward
  endpoints + forward diff count, genome-relative.

- [ ] **Step 1 (test first):** In the driver add `--fwd-validate`: for a fixed 5,000-seed uniform
  sample, compare GPU forward endpoint `(ae,be)` + `fdiff` to a **host reference computed by calling
  `align.c`'s wave forward-only** (reuse `wave_bench_cpu` logic restricted to the forward half via a
  helper that records `apath->aepos/bepos` and the forward diff). Assert **same-score (|Δdiff|==0)
  ≥ 95%** and exact-endpoint ≥ (record, no hard gate — analyzed later). Run → FAIL.
- [ ] **Step 2:** Implement `forward_sweep_warp`: one warp per seed; band `[low,hgh]` across lanes
  (≤ 8 diagonals/lane, `TR`-style ping-pong of the furthest-reach `V` and match `M` state, **int**
  positions, **no depth cap**); each `dif` step = warp sweeps the active `[max(low,-dif),…]` range,
  slides on matches until sentinel `4`, warp-argmax for best score; x-drop stop when
  `lasta < besta - TRIM_MLAG` (port the real rule, not a `−40` score gap). Bound reads to the seed's
  contig `[base,base+len)`.
- [ ] **Step 3:** Wire `wave_forward_batch`; run `--fwd-validate` on the 5k sample.
- [ ] **Step 4:** If same-score < 95%, debug against `WAVE_PORT_NOTES` (usual causes: band range,
  trim rule, sentinel handling) — iterate until ≥ 95% same-score. Record exact-endpoint %.
- [ ] **Step 5:** Commit (`feat: Stage-1 faithful forward sweep (warp-cooperative)`).

---

## Task 6: Stage-1 reverse sweep + full discovery, validated on the full distribution

**Files:**
- Modify: `gpu/wave_kernel.cu` (add `reverse_sweep_warp`, combine), `gpu/wave_bench_gpu.cu`.
- Reference: `align.c:919-1450`.

**Interfaces:**
- Produces: `wave_discover_batch(...) → ab,ae,bb,be,diffs` (full Local_Alignment-equivalent endpoints,
  genome-relative), matching the CPU `.cpuref` fields.

- [ ] **Step 1 (test first):** `--discover-validate` over the **full** `.seeds`: compare GPU
  `ab/ae/bb/be/diffs` to `.cpuref` (Task 3). Assert: **same-score ≥ 95%**, and report exact-endpoint %,
  depth-work ratio `mean(GPU diffs / CPU diffs)` (must be in [0.97,1.05] — neither cutting nor
  inflating work). Run → FAIL.
- [ ] **Step 2:** Implement `reverse_sweep_warp` (mirror of forward per `align.c:919-1450`; opposite
  extension direction, `lasta <= besta + TRIM_MLAG`). Combine fwd+rev into `wave_discover_batch`
  (endpoints = seed ± reverse/forward reach; diffs = fdiff+rdiff).
- [ ] **Step 3:** Run `--discover-validate` on the full 518,037. Iterate to ≥ 95% same-score and
  depth-ratio in range. Record exact-endpoint % and the same-score % (these feed the correctness
  report).
- [ ] **Step 4:** Commit (`feat: Stage-1 full discovery (fwd+rev), validated vs CPU full distribution`).

---

## Task 7: Stage-1 benchmark — full-distribution headline + fit map

**Files:**
- Modify: `gpu/wave_bench_gpu.cu` (timing + stratification + occupancy readout).
- Create: `gpu/results/stage1_fit.md` (generated results, committed).

**Interfaces:**
- Consumes: resident genome + full `.seeds`; CPU numbers from Task 3.

- [ ] **Step 1:** Add CUDA-event timing (best-of-N) for `wave_discover_batch` over the full seed set,
  **two bases**: (i) resident/kernel-only; (ii) + per-batch seed H2D and result D2H. Print both.
- [ ] **Step 2:** Add the per-(band,depth)-bucket GPU aln/s table (same buckets as Task 3) → the fit
  map; and the mechanism decomposition: achieved occupancy (via
  `cudaOccupancyMaxActiveBlocksPerMultiprocessor`) → concurrent warps; per-warp rate = throughput /
  concurrent-warps; print vs CPU per-core rate.
- [ ] **Step 3:** Run on the full distribution. Write `gpu/results/stage1_fit.md`: headline GPU-vs-CPU
  (full dist, both bases), the fit map, the per-warp×count vs per-core×cores decomposition, and the
  Stage-1 correctness numbers (same-score %, exact %, depth-ratio). **Label the deep-tail row as a
  diagnostic, not the headline.**
- [ ] **Step 4:** Commit (`results: Stage-1 full-distribution GPU-vs-CPU fit map`).

**Fit gate for Stage 2:** review `stage1_fit.md` with the human. Proceed to Stage 2 (trace) only if the
sweep fit justifies it; otherwise the study may conclude at Stage 1 (still a valid characterization).

---

## Task 8: Stage-2 trace-points (sparse-checkpoint) + divergence analysis

**Files:**
- Modify: `gpu/wave_kernel.cu` (checkpoints during sweep + backward-walk emit), `gpu/wave_bench_gpu.cu`.
- Create: `gpu/wave_validate.c`
- Modify: `Makefile`
- Reference: `align.c:860-895` (backward walk / trace emit), `align.h` `Check_Trace_Points`.

**Interfaces:**
- Produces: per-alignment GPU trace-points (diff, Δb per 100-bp panel), same encoding FastGA writes.

- [ ] **Step 1 (test first):** `wave_validate.c`: for a uniform 20k sample, run GPU trace + CPU
  `Local_Alignment` trace; assert GPU trace **100% `Check_Trace_Points`-valid** and **same edit
  distance**; compute byte-exact-vs-equally-optimal split; classify divergences (same-score-diff-path
  vs worse) and attribute (parallel band-order / argmax / width). Run → FAIL (no trace yet).
- [ ] **Step 2:** Add **sparse checkpoints** to the warp sweep: at each `tspace` boundary store a
  per-warp `Pebble`-equivalent (ptr,diag,diff,mark) in a bounded per-warp buffer (O(panels×band), NOT
  all d-layers — this is the key difference from the discarded trace kernel). After the sweep, lane-0
  walks back the checkpoints and bins (diff, Δb) per panel — port `align.c:860-895`.
- [ ] **Step 3:** Run `wave_validate` on 20k → iterate to 100% valid + same-distance. Produce the
  divergence classification+attribution table.
- [ ] **Step 4:** Add `Makefile` target `wave_validate`. Commit
  (`feat: Stage-2 sparse-checkpoint trace-points + divergence analyzer`).

---

## Task 9: End-to-end correctness + final benchmark + write-up

**Files:**
- Modify: `gpu/wave_bench_gpu.cu`
- Create: `gpu/results/final.md`

- [ ] **Step 1:** Full-distribution run WITH trace: GPU vs CPU throughput (both timing bases,
  trace included), and the full correctness report (same-score %, valid %, exact %, divergence
  attribution) from `wave_validate` scaled to the full set (or a large uniform sample, labelled).
- [ ] **Step 2:** Write `gpu/results/final.md`: the headline full-distribution speedup + range, the
  fit map, the mechanism decomposition, the correctness/divergence report, and a plain conclusion —
  *is this workload a good GPU fit, where, why, and what is the measured speedup*. Deep tail appears
  only as a labelled diagnostic.
- [ ] **Step 3:** Update the hub task doc + `project_gpu_acceleration` memory with the final numbers
  (draft only; do not commit/push without the human's approval).
- [ ] **Step 4:** Commit results (`results: final GPU wave-parallelism characterization`).

---

## Review / Self-check (writing-plans)

- **Spec coverage:** every spec §maps to a task — study(all)/faithful-port(T1,5,6,8)/full-dist-first
  (T2,3,7)/tiered-correctness(T6,8)/timing-range(T7,9)/no-integration(scope)/divergence-attrib(T8).
- **Type consistency:** `SeedRec`/`WaveSeedHeader` defined in Task 2, consumed unchanged in 3/4/6;
  `wave_ctx` + launcher signatures defined in Task 4, extended in 5/6/8.
- **Known risk carried in-plan:** Task 5/6 fidelity may plateau < 95% same-score under parallel order;
  the plan's response is debug-against-`WAVE_PORT_NOTES` then, if it is a genuine parallel-order
  limit, record it as attributed divergence (per the tiered bar) rather than force CPU scan order —
  this is the one place a task may "fail" its 95% gate and instead convert to a documented result;
  **surface it to the human at the Task 7 fit gate.**
