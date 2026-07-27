# Stage-1 benchmark: full-distribution GPU-vs-CPU headline + fit map

**Study:** GPU wave-parallelism characterization, Stage 1 (discovery: `wave_discover_batch` =
forward + reverse sweep → endpoints + edit distance, NO trace-points).
**Dataset:** GRCh38 × CHM13 (~3.1 Gbp each), full real `.1aln`, **518,037 seeds** — the complete
distribution extracted by `extract_seeds` (Task 2), no subsetting.
**Hardware:** 1× NVIDIA A100 80GB PCIe (`CUDA_VISIBLE_DEVICES=0`, 108 SMs), 2× EPYC 9124 (32
physical cores) for the CPU baseline (pinned `OMP_PROC_BIND=close OMP_PLACES=cores`).
**GPU kernel:** `wave_discover_batch` / `wave_discover_batch_timed`, committed at `0f4704f`
(Task 6), timing/fit-map instrumentation added this task (Task 7, NOT YET committed — see
milestone-commit policy at the end).

> **The full 518,037-seed distribution below is the HEADLINE.** The deep-tail row in the fit map
> (depth > 1000 edits) is a labelled **diagnostic** only — see the explicit tag in that table —
> never treat it as the result.

---

## Honest caveat (read this before the numbers)

The GPU `wave_discover_batch` computes the wave **sweep only** (forward+reverse endpoints and
edit-distance/`diffs`) — it does **not** emit trace-points. The CPU baseline, `Local_Alignment`
(FastGA's real hot-path aligner), **also emits trace-points in the same pass** as part of its
normal operation — it cannot be made to skip that work without patching FastGA itself. So this
Stage-1 comparison is **discovery-vs-(discovery+trace)**, which **slightly favors the GPU**. The
trace-inclusive, apples-to-apples comparison is Stage 2 (Tasks 8–9); this Stage-1 result is a
valid but partial characterization on its own terms, not a final verdict on GPU fit.

---

## Headline: full distribution, both timing bases

CUDA-event timed, best-of-5, `wave_discover_batch_timed`. The one-time genome H2D upload
(9.37 GB: A + Bfwd + Brev) took **0.64–1.0 s** across runs and is reported separately — it is
**not** included in either basis below (a resident engine pays it once, not per batch).

| Timing basis | Definition | Time (518,037 seeds) | Throughput |
|---|---|---:|---:|
| **(i) wave-engine / kernel-only** | genome resident; only the fwd+rev sweep kernels | 7.86 s | **65,884 seeds/s** |
| **(ii) realistic** | (i) + per-batch seed H2D + result D2H (excludes one-time genome H2D and per-call device-buffer alloc/free) | 7.87 s | **65,835 seeds/s** |

Bases (i) and (ii) are within 0.08% of each other here: at 518,037 seeds the per-batch seed
upload (~8.29 MB of `wave_seed` structs) and the 5 result arrays (~10 MB total) are negligible
next to the ~7.9 s of kernel compute — measured H2D 0.7–5.1 ms, D2H 0.8–6.0 ms across the 5 reps,
i.e. the "realistic" tax is **< 0.1%** of wall time at this batch size. The two-basis split
matters more at the per-bucket scale (below), where individual buckets are much smaller batches
and transfer overhead is a larger fraction of the (much smaller) kernel time.

Device-buffer `cudaMalloc`/`cudaFree` (unavoidable in this one-shot-call harness design) is
excluded from both bases — a long-lived engine would allocate such scratch once, so charging it
per-batch would be dishonest in the pessimistic direction, same as the genome upload is excluded
in the optimistic direction.

### CPU baseline (re-confirmed this task: pinned, best-of-3 for the 32-core number)

| | Throughput | Notes |
|---|---:|---|
| **1-core** | **1,200.7 aln/s** | `OMP_PROC_BIND=close OMP_PLACES=cores`; single measurement (T=1 has no thread-placement ambiguity); Task 3's original unpinned number was 1,212.0 aln/s — 0.9% apart, pinning made no material difference at T=1 |
| **32-core** | **38,713.8 aln/s** | best of 3 pinned reruns: 38,600.3 / 38,482.9 / **38,713.8**; `OMP schedule(dynamic,16)`, `SKIP_SINGLE=1` for reruns 2–3; Task 3's original unpinned number was 37,782.9 aln/s — pinning gave a genuine +2.5% here |

### The full-distribution ratio (the headline number)

| | basis (i) kernel-only | basis (ii) realistic |
|---|---:|---:|
| GPU vs CPU 1-core  | **54.87×** | **54.83×** |
| GPU vs CPU 32-core | **1.70×** | **1.70×** |

**GPU wave-discovery is ~55× a single CPU core, and ~1.70× all 32 pinned cores together** — the
honest headline, with the no-trace caveat above applying to both numbers.

---

## Mechanism decomposition: why this ratio, not some other number

`cudaOccupancyMaxActiveBlocksPerMultiprocessor` on the two Stage-1 kernels (block = 32 threads =
1 warp; `ptxas -v` reports 56/40 registers/thread for forward/reverse, 56 bytes shared mem/block
— small enough that occupancy is **not** register- or shared-mem-limited):

| kernel | registers/thread | max blocks/SM (measured) | limiting factor |
|---|---:|---:|---|
| `forward_sweep_warp` | 56 | 32 | **architectural max-resident-blocks/SM (32 on sm_80)** — not registers (56×32=1,792 regs/block leaves room for 36 by the 65,536-register/SM budget) or threads (32-thread blocks would allow 64 by the 2,048-thread/SM budget) |
| `reverse_sweep_warp` | 40 | 32 | same — also block-count-capped, not register-capped |

**Concurrent warps = min(32,32) blocks/SM × 108 SMs = 3,456 warps resident at once**, out of
518,037 seeds total (the launch grid is capped at `WV_POOL=8192` blocks, so each of the 3,456
concurrently-resident warps grid-strides through further seeds as earlier ones finish — roughly
150 seeds/warp over the whole run).

| | rate | unit count | per-unit rate |
|---|---:|---:|---:|
| GPU (basis i) | 65,884.3 aln/s | 3,456 concurrent warps | **19.06 aln/s per warp** |
| CPU (1-thread) | 1,200.7 aln/s | 32 cores | **1,200.7 aln/s per core** |

**Per-core / per-warp = 1,200.7 / 19.06 ≈ 63.0×.** This is the mechanistic explanation for the
aggregate ratio: a CPU core executing the real scalar `Local_Alignment` (full-width ALUs, deep
out-of-order pipelines, large per-core caches) is **~63× faster per lane** than one GPU warp
running the same recurrence with 32 SIMT lanes forced into lockstep on a problem that most of the
time doesn't need all 32 lanes — but the GPU has **3,456 concurrent
warps vs 32 CPU cores (108× more parallel units)**, so the aggregate ratio is the product of a
large concurrency advantage mostly (not entirely) cancelled by a much larger per-lane
efficiency disadvantage:

```
aggregate ratio ≈ (concurrent warps / cores) × (per-warp rate / per-core rate)
                ≈ (3456/32) × (19.06/1200.7)
                ≈ 108.0 × 0.01588
                ≈ 1.71×
```

This reconstructs the measured 1.70× (basis i vs 32-core) to within rounding — the decomposition
is internally consistent and is the "how" behind the headline ratio: **many slow warps, not few
fast ones.**

> **⚠ Read the decomposition as accounting, not as two measurements (noted 2026-07-27).** Only the
> two aggregate throughputs are measured: GPU 65,884 seeds/s and CPU-32-core 38,714 aln/s, whose
> ratio *is* 1.70×. The per-core figure 1,200.7 aln/s is also a direct stopwatch number
> (`wave_bench_cpu.c`, single-threaded loop over all 518,037 seeds calling the real
> `Local_Alignment`, genome preloaded, 431.5 s). But **19.06 aln/s per warp is not measured — it is
> 65,884 ÷ 3,456**, where 3,456 is the CUDA occupancy calculator's *analytic* resident-warp ceiling.
> Nsight Compute measures **achieved** occupancy at 29.26%, i.e. ~2,023 warps actually busy. Redo the
> split with that: (2023/32) × (32.6/1200.7) = 63.2 × 0.02715 ≈ **1.72×** — the same answer, because
> the warp count cancels between the two factors. So "108× more workers, each 63× weaker" and "63×
> more workers, each 37× weaker" are the same fact stated two ways; neither pair is independently
> observed.
>
> **What *is* independently measured (`ncu_profile.log`, `ncu_sass_fwd.csv`):** per-SASS
> `Avg. Threads Executed` shows **76.5% of instruction issue runs on exactly one lane** (the
> `if (lane == 0)` phase-2/3 code) and returns only 19.2% of the useful work, while phase-1 band
> work at 17–32 live lanes is 12.1% of issue and 64.8% of the work — whole-kernel average 4.52 of 32
> live lanes. That, not memory bandwidth, is the per-lane efficiency disadvantage this section
> describes. Nothing in either kernel is bandwidth-bound (DRAM 0.01–0.44% of peak, L1 98–99.8%).

**Unexploited headroom (observation, not implemented/measured further here):** the limiting
factor is the block-count cap (32 blocks/SM), not registers — a 32-thread block occupies only
1,024 of each SM's 2,048 threads at the max resident block count, so half the SM's thread
capacity is architecturally unusable by this one-warp-per-block design regardless of how light
the kernel is. Packing multiple independent warps per block (e.g. 4–8 warps/block, unchanged
per-warp logic) could in principle raise concurrent warps toward ~6,912–13,824 and roughly
double-to-quadruple the aggregate ratio. This is a real design lever for a future task; it is
flagged here, not exercised, since Task 7's scope is characterization, not further optimization.

---

## Fit map: GPU aln/s by (band, depth) bucket

Same bucket boundaries as `wave_bench_cpu.c`'s stratified table (`BAND_HI = {8,32,128,512,∞}`,
`DEPTH_HI = {10,50,200,1000,∞}`; band = `|(ae-ab)-(be-bb)|`, depth = `diffs`), so GPU and CPU
per-bucket throughput are directly comparable bucket-for-bucket. Bucketing uses the GPU's **own**
discovered endpoints/diffs (not the CPU reference), so bucket populations differ from the CPU
table by only a handful of seeds (99.2% same-score means the two engines almost always agree on
which bucket a seed lands in). Per-bucket GPU rate = kernel-only time re-measured on just that
bucket's seed subset (best-of-3), matching the CPU table's pure-compute (no-I/O) methodology.

### GPU (this task)

```
fit map: GPU aln/s (kernel-only basis), same (band,depth) buckets as wave_bench_cpu.c
    depth\band  <=8         <=32        <=128       <=512       <=999999999
            <=10    719844.5       722.6           -           -           -   (n=23236)
            <=50    456583.2    120709.0      7141.2           -           -   (n=42698)
            <=200    637851.7    476760.4    134464.1      1944.1           -   (n=161816)
            <=1000    258341.5    291746.6    213906.3      9448.0           -   (n=236164)
            <=999999999      4314.5      6277.3      4541.4      1773.7       274.0   (n=54123)   <-- DIAGNOSTIC (deep tail), NOT the headline
```

### CPU (this task's re-confirmed pinned run, 1-core, per-seed timers — same run that produced 1,200.7 aln/s above)

```
stratified aln/s  (rows=depth(path.diffs) bucket, cols=band-proxy |Δa-Δb| bucket)
    depth\band  <=8         <=32        <=128       <=512       <=999999999
            <=10      8188.2      7539.0           -           -           -   (n=23313)
            <=50      7852.1      7613.8      6892.9           -           -   (n=42693)
            <=200      7169.2      6565.0      5660.8      8117.2           -   (n=161626)
            <=1000      4827.4      4476.6      3655.1      2232.5           -   (n=235771)
            <=999999999       695.6       543.6       236.3        42.7        19.1   (n=54634)
```

### Where the GPU wins / loses (bucket-by-bucket ratio, GPU/CPU)

| depth \ band | ≤8 | ≤32 | ≤128 | ≤512 | ≤999999999 |
|---|---:|---:|---:|---:|---:|
| ≤10   | **87.9×** | 0.10× | — | — | — |
| ≤50   | **58.2×** | **15.9×** | 1.04× | — | — |
| ≤200  | **89.0×** | **72.6×** | **23.8×** | 0.24× | — |
| ≤1000 | **53.5×** | **65.2×** | **58.5×** | 4.23× | — |
| >1000 (**diagnostic**) | 6.20× | **11.6×** | **19.2×** | **41.5×** | **14.3×** |

**The fit is exactly the shape the study expected: the GPU wins big (24–90×) wherever the active
band is narrow (≤128), across every depth including the shallowest/most common seeds** — one
warp's 32 lanes cover the whole band in a single pass there, so the GPU pays roughly the same
tiny fixed cost per seed regardless of depth, while the CPU's per-seed cost scales with the
actual diff count. **The GPU loses once drift exceeds ~128–512** — but the causal chain is more
indirect than "wider band → more tiles." **drift is bounded by depth** (net indels can never
exceed the total edit distance, so `drift ≤ diffs` for every seed; the per-cell N matrix below is
consequently lower-triangular). So the high-drift cells are populated *almost entirely by the
deepest, longest alignments*: of the **8,059** seeds with drift in (128, 512], **6,789** sit in the
deep-tail (>1000-edit) row, and **all 966** seeds with drift > 512 are in it (corrected counts —
an earlier draft here misread the per-cell *rate* 9,448.0 aln/s as a count; the validated per-cell
populations are tabulated below). A larger sustained band does force the 32 lanes to serialize a
few more tiles per step (a 2–5× effect), but the dominant factor is that per-warp work explodes on
these deep seeds exactly where per-warp efficiency is already worst — and, as the per-cell
time-share below shows, that thin deep-tail minority (10.5% of seeds) is where **~80% of CPU wall
time** actually goes. The huge 24–90× narrow-band wins cost almost no time, so they cannot move the
aggregate; the headline lands at 1.70× because it is set by the deep tail, where both machines
spend their clock.

The `≤10 depth × ≤32 band` cell (722.6 aln/s, ratio 0.10×; its seeds are part of the depth≤10
row, n=23,236) is a below-CPU outlier at this granularity. It is most likely a **small-N
measurement artifact**: per-bucket GPU rates isolate each bucket into a standalone re-run, so a
small, trivially-shallow bucket is dominated by fixed launch overhead rather than true per-seed
cost — which *understates* the GPU, i.e. errs against it, not for it. It is at worst weakly
consistent with the "wide band hurts" trend; we do not lean on this single cell.

**Deep-tail diagnostic (labelled, not the headline):** the last row (depth > 1000 edits,
n=54,123, ~10.4% of all seeds) is where the CPU is slowest per-seed by far (19.1–695.6 aln/s) —
exactly the regime a GPU should help most in principle, and indeed the GPU/CPU ratio there climbs
to 41.5× at the widest resolvable band — but this row is a small, biased slice of the
distribution (the hardest seeds only) and **must never stand in for the full-distribution
headline above.**

---

## Per-cell population and time-share (added post-hoc; not part of the original Task-7 run)

Bucketed directly from `wave.cpuref` (`Local_Alignment`'s own endpoints/diffs) with the same
`{8,32,128,512,∞}` × `{10,50,200,1000,∞}` boundaries. **Row totals reproduce the CPU stratified
table above bit-for-bit** (23,313 / 42,693 / 161,626 / 235,771 / 54,634 → 518,037), which validates
the bucketing. Script: `percell.py` (reads cpuref, histograms, cross-checks row totals).

### Per-cell N (count of alignments)

```
 depth\drift    <=8      <=32     <=128     <=512      >512       ROW
     <=10     23311        2         0         0         0     23313
     <=50     38565     4108        20         0         0     42693
    <=200     82392    67163     12065         6         0    161626
   <=1000     54498   102539     77470      1264         0    235771
    >1000      9766    15575     21538      6789       966     54634
      COL    208532   189387    111093      8059       966    518037
```

Two structural facts fall out: (1) the matrix is **lower-triangular** — `drift ≤ depth` always, so
drift and depth are correlated axes, not independent ones; (2) the distribution is **extremely
narrow-skewed** — drift ≤ 32 covers **77%** of all seeds, while drift > 512 is **966 seeds =
0.19%**. The `depth≤10 × drift≤32` cell that read 0.10× in the ratio table holds exactly **N = 2**
alignments — confirming that outlier is pure small-N launch-overhead noise, not signal.

### Time-share % (share of total wall time each cell consumes)

`time_cell = N_cell / rate_cell`, normalized. **CPU is the rigorous one** (per-seed timers are
additive; the cell times sum to 431.5 s = the measured 1-core wall, 518,037 / 1200.7). GPU per-cell
rates come from *isolated per-bucket re-runs*, so they are **not** additive (they sum to 18.4 s, not
the real full-batch 7.86 s) and the small deep-tail batches carry inflated fixed launch overhead —
so the GPU shares below are **directional, not exact**, and overstate the deep-tail concentration.

```
CPU 1-core time-share %          GPU kernel-only time-share % (directional)
 depth\drift  <=8   <=32  <=128 <=512  >512  ROW    depth\drift  <=8   <=32  <=128 <=512  >512  ROW
     <=10    0.66  0.00     .     .     .   0.7%        <=10    0.18  0.02     .     .     .   0.2%
     <=50    1.14  0.13  0.00     .     .   1.3%        <=50    0.46  0.18  0.02     .     .   0.7%
    <=200    2.66  2.37  0.49  0.00     .   5.5%       <=200    0.70  0.76  0.49  0.02     .   2.0%
   <=1000    2.62  5.31  4.91  0.13     .  13.0%      <=1000    1.15  1.91  1.97  0.73     .   5.7%
    >1000    3.25  6.64 21.12 36.84 11.72  79.6%       >1000   12.29 13.47 25.75 20.78 19.14  91.4%
```

**The deep-tail row (depth > 1000) is 10.5% of the alignments but 79.6% of CPU wall time**; a single
cell — `depth>1000 × drift≤512`, N = 6,789 (1.3% of all seeds) — is **36.8%** of CPU time. This is
the exact same workload-imbalance shape the CPU sort+align phase suffers: a ~1% minority of monster
alignments sets the clock. It is also *why* discovery only reaches 1.70× rather than the 24–90× of
the narrow-band majority — the majority is nearly time-free, and the aggregate is decided by the
deep tail where GPU/32-core is close to parity. The direct implication for further GPU work: a
uniform per-warp optimization (e.g. shared-memory checkpoints) chips at the wrong thing; the payoff
is in giving the deep-tail monsters a different execution strategy (multi-warp cooperation, or
heterogeneous CPU/GPU split), the same lesson as the CPU work-stealing fix.

---

## Stage-1 correctness (from Task 6 / `--discover-validate`, full 518,037 seeds; re-confirmed bit-for-bit this task)

| Metric | Result | Gate | Verdict |
|---|---|---|---|
| **same-score** (`\|Δdiffs\|==0`) | **513,897/518,037 = 99.201%** | ≥ 95% | **PASS** |
| **depth-work ratio** mean(GPU diffs / CPU diffs) | **1.0055** | ∈ [0.97, 1.05] | **PASS** |
| exact endpoint (all 4 of ab/ae/bb/be, tol 0) | 510,580/518,037 = **98.561%** | (report) | — |
| endpoint within 2bp (all 4) | 513,273/518,037 = 99.080% | (report) | — |
| band overflow | 0/518,037 | — | — |

This task's `--stage1-fit` (a fully independent code path from `--discover-validate`, sharing
only the underlying kernels) reproduced **same-score 99.201%, depth-ratio 1.0055, 0 overflow**
bit-for-bit — confirming the new `wave_discover_batch_timed` timing instrumentation is
functionally identical to the existing, already-validated `wave_discover_batch` (same kernels,
same launch parameters, only cudaEvent markers and a restructured allocation path added).

The exact-endpoint 98.56% is what proves the reverse sweep is load-bearing and correct — `ab`/
`bb` are produced only by `reverse_sweep_warp`, so matching all four endpoints byte-exact on
98.56% of seeds means the reverse wave reproduces `Local_Alignment`'s `abpos`/`bbpos`, not merely
the forward half (see Task 6's report for the full divergence attribution of the residual 0.8%).

---

## Fit gate for Stage 2

**Recommendation: proceed to Stage 2 (trace-points, Tasks 8–9).** The Stage-1 sweep is faithful
(99.2% same-score, a correct and load-bearing reverse sweep) and shows a real,
mechanistically-understood win in its natural regime — narrow diagonal band, which covers the
large majority of the distribution (every bucket with band ≤128 beats the CPU by 24–90×,
regardless of depth). The 1.70× full-distribution aggregate undersells that: it is dragged down
by a shrinking-but-nonzero wide-band minority whose per-warp cost compounds with depth. Because
Stage 1 omits trace-point emission — work the CPU baseline unavoidably includes — this aggregate
also structurally favors the GPU a bit; Stage 2 is needed before a genuine end-to-end verdict can
be drawn, and this document explicitly declines to draw one yet (see the caveat near the top).

---

## Files changed this task (NOT committed — milestone-commit policy; controller commits after human review)

- `gpu/wave_kernel.h` — added `wave_discover_batch_timed` (CUDA-event-split H2D/kernel/D2H
  timing) and `wave_query_occupancy` declarations.
- `gpu/wave_kernel.cu` — implemented both; `wave_discover_batch` (Task 6) is untouched.
- `gpu/wave_bench_gpu.cu` — added `--stage1-fit` (self-contained: best-of-5 full-distribution
  timing on both bases, occupancy-based mechanism decomposition, per-(band,depth) fit map via
  best-of-3 per-bucket subset re-timing, correctness re-check) plus `--cpu-1core=`/
  `--cpu-32core=` flags so the freshly re-measured pinned CPU numbers feed the printed ratios
  instead of stale constants.
- `gpu/results/stage1_fit.md` — this document (new).

## Reproduce

```bash
# CPU: pinned re-confirmation (1-core once + .cpuref, then 32-core best-of-3)
OMP_PROC_BIND=close OMP_PLACES=cores ./gpu/wave_bench_cpu \
  /scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/output.1aln \
  gpu/wave.seeds gpu/wave.cpuref 32                        # run 1: 1-core + 32-core + .cpuref
SKIP_SINGLE=1 OMP_PROC_BIND=close OMP_PLACES=cores ./gpu/wave_bench_cpu \
  /scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/output.1aln \
  gpu/wave.seeds gpu/wave.cpuref 32                        # runs 2,3: 32-core only, best-of-3

# GPU: full-distribution headline + fit map + occupancy decomposition
make wave_bench_gpu
CUDA_VISIBLE_DEVICES=0 ./gpu/wave_bench_gpu \
  /scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/output.1aln \
  gpu/wave.seeds --stage1-fit --cpu-1core=1200.7 --cpu-32core=38713.8
```
