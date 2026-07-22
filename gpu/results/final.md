# Final benchmark: trace-inclusive, apples-to-apples GPU-vs-CPU end-to-end (Stage 2)

**Study:** GPU wave-parallelism characterization, Stage 2 (**trace-inclusive**: `wave_trace_batch` =
forward + reverse sweep **plus** sparse-checkpoint trace-point emission, in one pass — the same
discover+trace work FastGA's CPU `Local_Alignment` does).
**Dataset:** GRCh38 × CHM13 (~3.1 Gbp each), full real `.1aln`, **518,037 seeds** — the complete
distribution from `extract_seeds` (Task 2), no subsetting.
**Hardware:** 1× NVIDIA A100 80GB PCIe (`CUDA_VISIBLE_DEVICES=0`, 108 SMs), 2× EPYC 9124 (32
physical cores) for the CPU baseline (pinned `OMP_PROC_BIND=close OMP_PLACES=cores`).
**GPU kernel:** `wave_trace_warp` / `wave_trace_batch` (Task 8, in the working tree), timed via the
new `wave_trace_batch_timed` (Task 9). CPU baseline reused verbatim from Stage-1
(`gpu/results/stage1_fit.md`) — the CPU side is unchanged, it always emitted trace.

---

## Why this document exists: the apples-to-apples correction

Stage-1's headline **1.70× vs 32 CPU cores** was **discovery-only** and **structurally favored the
GPU**: the GPU computed only endpoints + edit distance, while the CPU baseline (`Local_Alignment`)
*also* emitted trace-points in the same pass and could not be made to skip that work. Task 8 gave the
GPU a trace kernel too, so Stage 2 measures the **trace-inclusive** GPU throughput — GPU
discover+**trace** over all 518,037 seeds — against the **same** CPU discover+trace baseline. **This
is the honest end-to-end number.** As anticipated in the Stage-1 caveat, it is markedly lower than
1.70×.

---

## (1) HEADLINE — trace-inclusive, full distribution, both timing bases

CUDA-event timed, best-of-3, `wave_trace_batch_timed`, chunked at 20,000 seeds/batch (bounds the
`O(n·trace_stride)` device/host trace buffer and matches a real batched engine). The one-time genome
H2D upload (9.37 GB) is reported separately, as in Stage-1, and is **not** in either basis. Kernel
time was stable to ±1 ms across the 3 reps (58572.7 / 58573.6 / 58573.4 ms).

| Timing basis | Definition | Time (518,037 seeds) | Throughput |
|---|---|---:|---:|
| **(i) wave-engine / kernel-only** | genome resident; only the `wave_trace_warp` launch (discover+trace) | 58.57 s | **8,844.3 seeds/s** |
| **(ii) realistic** | (i) + per-batch seed H2D + result-meta D2H + **compacted** trace D2H | 58.63 s | **8,835.8 seeds/s** |

Bases (i) and (ii) are within **0.1%** — the per-batch transfers are negligible next to the ~58.6 s
kernel: seed H2D 25.5 ms (8.3 MB), meta D2H 3.7 ms (6 arrays × 2.07 MB), and the **compacted** trace
D2H is **378 MB total** (Σ tlen·2 bytes over all seeds) → a few ms.

> **Transfer-honesty note.** The harness's `trace` buffer is strided at a worst-case 32,768 uint16
> per seed, so a literal D2H of it moves **33.95 GB** of mostly-zero padding (2.46 s) — a **harness
> artifact**, not an engine cost. A real engine compacts traces device-side and moves only the
> 378 MB of real trace data. Basis (ii) charges the **compacted** cost (the honest figure); charging
> the full strided copy instead would still only drop throughput to 8,483 seeds/s (0.22×), so this
> choice does not change the verdict.

### CPU baseline (pinned, reused verbatim from Stage-1 — CPU side is unchanged)

| | Throughput | Notes |
|---|---:|---|
| **1-core** | **1,200.7 aln/s** | `OMP_PROC_BIND=close OMP_PLACES=cores`; discover+trace (always) |
| **32-core** | **38,713.8 aln/s** | best of 3 pinned reruns; `schedule(dynamic,16)` |

### The trace-inclusive ratio — THE headline

| | basis (i) kernel-only | basis (ii) realistic |
|---|---:|---:|
| GPU vs CPU 1-core  | **7.37×** | **7.36×** |
| GPU vs CPU 32-core (committed 512-cap) | **0.23×** | **0.23×** |

**Honest headline: as a drop-in warp-per-seed trace port, the trace-inclusive GPU wave engine runs
at 0.23–0.32× of 32 pinned CPU cores — i.e. below 1×, the number this study set out to produce.**
The range's two ends are both real, measured, and committed:

| Launch grid (concurrency) | resident warps | seeds/s (basis i) | vs 32-core | log |
|---|---:|---:|---:|---|
| **committed default (TR_POOL=512)** | 512 | **8,844.9** | **0.23×** | `stage2_fit.log` |
| TR_POOL=2048 | 2,048 | 11,921.0 | 0.31× | `stage2_fit_pool2048.log` |
| **kernel's own occupancy ceiling (TR_POOL=2160)** | 2,160 | **12,317.4** | **0.32×** | `stage2_fit_pool2160_ceiling.log` |

**The 512 cap is a self-imposed memory budget, NOT a hardware limit.** `wave_query_trace_occupancy`
reports the kernel's *own* occupancy ceiling as **2,160 warps** (20 blocks/SM × 108 SMs,
register-limited); the committed launch grid is capped at 512 only to hold the checkpoint scratch to
~2 GB. Lifting that arbitrary cap to the kernel's own ceiling — **no redesign, only a larger launch
grid** — already lifts the ratio to **0.32×** while using just **~19.2 GB of the A100's 80 GB**
(measured peak at 2,160 warps). Both timing bases agree to 0.1% at every pool size. The realistic
"best this drop-in kernel does" is **0.32×**; the pessimistic committed-cap edge is **0.23×**.

### Discovery-only 1.70× — clearly labelled SUB-result (NOT the headline)

For reference, Stage-1's discovery-only measurement (`gpu/results/stage1_fit.md`) on the same seeds:
GPU 65,884 seeds/s = **54.87× vs 1-core, 1.70× vs 32-core**. That number **excluded** trace on the
GPU while including it on the CPU; it is a valid partial characterization of the *sweep* alone, not
an end-to-end verdict. Adding trace to the GPU drops the aggregate from 1.70× to 0.23× — a **7.4×**
fall, decomposed next.

---

## (2) Fit map — reference Stage-1's; how trace changes it

The per-(band,depth) fit map is carried from Stage-1 (`gpu/results/stage1_fit.md` §"Fit map"); it is
**not** independently re-measured per-bucket for the trace kernel here. Its **shape is preserved**,
and the reason is the mechanism finding below: at the committed concurrency (512 warps) the trace
kernel's **per-warp** rate (17.27 aln/s) is within ~10% of discovery's (19.06 aln/s), and the trace
work a warp adds is proportional to alignment length (one checkpoint per `tspace`=100 bp panel), not
to band width. So the same qualitative fit holds — narrow diagonal band (≤128) is where a warp's 32
lanes cover the band in one pass and the GPU is most efficient; wide band still serializes. What
trace changes is **not the per-bucket shape but the global scale**: every bucket is pulled down by
the ~6.75× drop in *resident warps* (below), a uniform aggregate factor, so the fit map's relative
structure carries over while its absolute throughput scales down ~7.4×.

*(Limitation: this is an argued, not a re-measured, per-bucket claim — the per-warp near-equality and
the length-proportional trace cost justify it, but a bucket-by-bucket trace fit map was not run.)*

---

## (3) Mechanism decomposition — why 0.23×, and why trace costs 7.4×

`wave_query_trace_occupancy`: `wave_trace_warp` allows **20 blocks/SM × 108 SMs = 2,160** by
occupancy, but the committed launch grid is capped at **TR_POOL = 512** blocks. **That cap is a
self-imposed memory budget, not a hardware limit** — it exists only because each warp needs an
`O(checkpoints)` scratch (`TR_CELLS_MAX = 262,144` × 4 int arrays), so 512 warps already consume
~2.15 GB of checkpoint buffers, and the default keeps that conservative. At the committed default,
**resident warps = min(2160, 512) = 512** — the number the 0.23× headline is measured at. Lifting the
cap to the 2,160-warp occupancy ceiling (measured in §1, ~19.2 GB of 80 GB) raises the ratio to
0.32× with no algorithm change; the concurrency sweep below shows why that is the concurrency knob's
ceiling — and why a *traffic*-reduction redesign, not concurrency, is the remaining lever.

| | rate | unit count | per-unit rate |
|---|---:|---:|---:|
| GPU trace (basis i) | 8,844.3 aln/s | 512 resident warps | **17.27 aln/s per warp** |
| GPU discovery (Stage-1) | 65,884.3 aln/s | 3,456 resident warps | 19.06 aln/s per warp |
| CPU (1-thread) | 1,200.7 aln/s | 32 cores | 1,200.7 aln/s per core |

**Aggregate reconstruction (trace-inclusive):**
```
ratio vs 32-core ≈ (resident warps / cores) × (per-warp rate / per-core rate)
                 ≈ (512 / 32) × (17.27 / 1200.7)
                 ≈ 16.0 × 0.01438
                 ≈ 0.23×      (matches the measured 0.23×)
```

**Where the 7.4× discovery→trace fall comes from — it is concurrency, not per-warp slowdown.**
Discovery ran **3,456** resident warps at 19.06 aln/s/warp; trace runs **512** at 17.27 aln/s/warp.
The **per-warp rate barely changed** (19.06 → 17.27, ×0.91) — a warp emits trace at almost the same
rate it discovers. The whole aggregate loss is the **3,456 → 512 = 6.75× collapse in resident
warps**, forced by the trace kernel's per-warp checkpoint memory. 6.75 × (1/0.91) ≈ 7.4×, matching
65,884 → 8,844.

**Two distinct levers — separate them carefully.**

**(a) The CONCURRENCY knob (raise resident warps) — measured, and it is nearly closed by bandwidth.**
Raising `TR_POOL` from the committed 512 up to the kernel's own occupancy ceiling was measured
directly and committed (`stage2_fit_pool2048.log`, `stage2_fit_pool2160_ceiling.log`):

| launch grid | resident warps | seeds/s | per-warp aln/s | vs 32-core |
|---|---:|---:|---:|---:|
| 512 (committed) | 512 | 8,844.9 | 17.27 | 0.23× |
| 2,048 | 2,048 | 11,921.0 | 5.82 | 0.31× |
| 2,160 (occupancy ceiling) | 2,160 | 12,317.4 | 5.70 | 0.32× |

4.2× the warps (512 → 2,160) returned only **1.39× throughput** (8,844 → 12,317), and per-warp rate
**cratered 17.27 → 5.70 aln/s** (×0.33): the added warps drove proportionally more global-memory
checkpoint traffic and saturated bandwidth, so ~3/4 of the added concurrency was cancelled. So
**this lever is robustly closed** — the concurrency knob alone tops out at 0.32×. (This differs from
Stage-1's discovery kernel, whose limiter was the block-count cap with bandwidth to spare; trace's
limiter is bandwidth itself, because of the checkpoint scratch it must stream.)

**(b) The TRAFFIC-REDUCTION lever (shrink per-warp checkpoint traffic) — UNTESTED, and it attacks
exactly the bandwidth this study blames.** The concurrency sweep saturated bandwidth *at the current
per-warp traffic*; it says nothing about lowering that traffic. The per-warp checkpoint buffer
`TR_CELLS_MAX = 262,144` is over-provisioned **~1,400×** vs the typical **~183-panel** wave — almost
every warp's live checkpoints would fit in a right-sized **shared-memory** buffer, with a global
overflow fallback for the rare deep wave. That redesign would move the ≥99% common-case checkpoint
stream off the global-memory bus entirely — directly cutting the traffic that caps this kernel. It
is **unquantified headroom**: by the per-warp decomposition it could **plausibly approach parity**.
It is *not* a tuning knob (it is a real kernel redesign, out of scope here), and it is **not** free
of ceilings — trace is register-capped at **2,160 resident warps** vs discovery's 3,456, so even a
bandwidth-freed trace kernel's realistic ceiling is **likely ≲1×**, not the 1.70× of discovery.

**Bottom line: the sub-1× verdict is robust to the concurrency knob (0.23–0.32×); the
traffic-reduction redesign is unquantified headroom, plausibly near-parity, still likely ≲1×.**

---

## (4) Correctness / divergence — trace-point validity (100k uniform sample)

`gpu/wave_validate` on a **100,000-seed uniform sample** (5× the Task-8 20k gate; the 20k result at
`gpu/results/wave_validate_20k.log` reproduced identically at this scale):

| Metric | Result | Gate | Verdict |
|---|---|---|---|
| **Check_Trace_Points-valid** (align.c's own validator) | **99,866 / 99,866 = 100.000%** | 100% | **PASS** |
| GPU trace INVALID | 0 | — | — |
| same edit distance vs CPU (\|Δdiffs\|==0) | 99,123 / 99,866 = 99.256% | — | — |
| endpoint-exact (all 4, tol 0) | 98,465 / 99,866 = 98.597% | — | — |
| **in-scope byte-exact** (endpoints+tlen+every (d,b) identical) | **97,904 / 97,904 = 100.000%** | — | — |
| **in-scope same-score** | 97,904 / 97,904 = 100.000% | ≥ 99% | **PASS** |
| in-scope genuinely-worse path | 0 | — | — |
| out-of-scope (CPU DUB_TRIM short-end re-center) | 1,962 (byte-exact 5) | — | (excluded) |
| trace-buffer overflow (tlen > 32,768 ⇒ > ~1.6 Mb aln) | 134 / 100,000 = 0.134% | — | (excluded) |

**Reading:** every GPU trace is a structurally legal FastGA trace (100% `Check_Trace_Points`), and on
the **97,904 in-scope seeds** (CPU alignment long at both ends, so no DUB_TRIM re-center) the GPU
trace is **byte-for-byte identical** to `Local_Alignment`'s — same endpoints, same `tlen`, every
`(diff, Δb)` panel. All divergence lives in the **out-of-scope** set (1,962): those are cases where
the CPU's `Local_Alignment` applied its **DUB_TRIM short-end re-center/re-extend** (align.c:1535,
1551–1576) — a full-routine post-process the wave port deliberately omits, the **same exclusion
Stage-1 made**. It is not a GPU trace bug. The 0.134% overflow is the deep tail of >1.6 Mb
alignments exceeding the per-seed trace scratch; they are excluded, not miscomputed.

### LIMITATION (stated honestly, per the Task-8 review)

The in-scope / out-of-scope split classifies "short end" using the **DUB_TRIM predicate on the GPU's
own pre-recenter endpoints** (the principled detector, since `Local_Alignment` tests fshort/rshort on
its *pre*-recenter wave endpoints — exactly what the GPU produces). The out-of-scope
"genuinely-worse" seeds were **not independently re-checked against the CPU's own pre-recenter
shortness** — doing so would require intercepting `forward_wave`/`reverse_wave` before the DUB_TRIM
branch inside align.c (an invasive patch), which was judged not cheap and not done. The residual is
tightly bounded by the whole-sample facts — **99.256% global same-distance, 100% Check_Trace_Points,
100% in-scope byte-exact, 0 in-scope genuinely-worse** — but the out-of-scope worse-path set is not
separately audited against CPU pre-recenter endpoints. This is a known, bounded limit, not a silent
gap.

---

## (5) CONCLUSION — is this workload a good GPU fit?

**Plain answer: as a drop-in warp-per-seed trace port, the full end-to-end aligner (discover + trace)
runs at 0.23–0.32× of 32 pinned CPU cores — below 1×, on both timing bases, at proven correctness.
0.23× is the committed 512-warp memory budget; 0.32× is the same kernel at its own 2,160-warp
occupancy ceiling (no redesign). This is where the GPU stops being competitive as a drop-in port —
NOT a proof the GPU fundamentally cannot.**

- The GPU wave engine is **faithful**: 100% valid trace, 100% byte-exact on all 97,904 in-scope
  seeds, all divergence attributed to an out-of-scope CPU post-process. Correctness is not the issue.
- The GPU **beats a single CPU core ~7.4×**, but a single A100 does **not** beat 32 EPYC cores once
  trace is included. The mechanism is clear and honest: trace-point emission does not slow a warp
  much (per-warp 17.27 vs discovery's 19.06 aln/s at 512 warps), but at the committed cap its
  `O(checkpoint)` per-warp memory holds resident warps at 512 (vs discovery's 3,456) — the concurrency
  collapse behind the 1.70× → 0.23× fall.
- **The CONCURRENCY knob is robustly closed, but it is not the only lever.** Raising the launch grid
  to the kernel's own occupancy ceiling was **measured and committed**: 512 → 2,160 warps returned
  only 1.39× throughput (0.23× → 0.32× vs 32-core, using just ~19.2 GB of 80 GB) because the kernel
  is **global-memory-bandwidth-bound** on the checkpoint stream. So concurrency alone tops out at
  0.32×. But the untested **traffic-reduction** redesign — right-sized **shared-memory** checkpoints
  (the per-warp `TR_CELLS_MAX`=262,144 buffer is over-provisioned ~1,400× vs the typical ~183-panel
  wave) with a global overflow fallback for the rare deep wave — attacks exactly that bandwidth, and
  by the per-warp decomposition could **plausibly approach parity**. It is unquantified headroom, not
  a tuning knob; and trace is register-capped at 2,160 warps (vs discovery's 3,456), so its realistic
  ceiling is **likely ≲1×**.
- **Where a GPU could still fit:** the narrow-band majority of the *discovery* sweep is genuinely
  GPU-friendly (Stage-1: 24–90× per-warp-efficient buckets), and this cross-checks the memory note's
  independent A3b finding that the trace kernel — not discovery — is what starves (~10–20k aln/s
  there, ~12k here). A productive GPU design would keep discovery on the GPU and either keep trace on
  the CPU or redesign trace to move its checkpoint stream into shared memory. As a **drop-in**
  warp-per-seed port, the trace stage is where the GPU stops being competitive.

**Full-distribution headline is the verdict; the deep tail (>1000 edits / >1.6 Mb) appears only as
the 0.134% overflow diagnostic and never stands in for it. No cherry-picking: the honest range is
0.23× (committed 512-cap) to 0.32× (occupancy ceiling) — a drop-in trace port is sub-1×, with a
shared-memory-checkpoint redesign as unquantified, plausibly-near-parity headroom.**

---

## Files changed this task (NOT committed — milestone-commit policy)

- `gpu/wave_kernel.h` — declared `wave_trace_batch_timed` (4-way CUDA-event split: seed H2D /
  kernel / meta D2H / strided-trace D2H) and `wave_query_trace_occupancy`.
- `gpu/wave_kernel.cu` — implemented both; the Task-8 `wave_trace_warp` **kernel algorithm is
  untouched**. The only change is that the launch-grid cap `TR_POOL` (committed default **512**) is now
  run-time overridable via env `TR_POOL` through a `tr_pool_cap()` helper, so the §1/§3 concurrency
  sweep (2,048 and the 2,160 occupancy ceiling) is a pure launch-size measurement, not an algorithm
  change. The committed source default remains 512 (verified: env-off reproduces 8,844.9 seeds/s,
  17.275 aln/s/warp, 0.23× — bit-for-bit the prior committed number).
- `gpu/wave_bench_gpu.cu` — added `--stage2-fit` (chunked best-of-N trace-inclusive full-distribution
  timing on both bases, occupancy/concurrency readout, compacted-vs-strided trace-transfer honesty,
  same-score cross-check vs `.cpuref`) and the `trace_time_chunked` helper.
- `gpu/results/final.md` — this document (new).
- `gpu/results/stage2_fit.log` (committed 512-cap), **`gpu/results/stage2_fit_pool2048.log`**,
  **`gpu/results/stage2_fit_pool2160_ceiling.log`** (concurrency sweep — the committed evidence behind
  the 0.31×/0.32× upper end), `gpu/results/wave_validate_100k.log` — raw logs (new).

## Reproduce

```bash
make wave_bench_gpu wave_validate

# trace-inclusive headline + mechanism (best-of-3), committed 512-warp cap -> 0.23x
STAGE2_NREPS=3 CUDA_VISIBLE_DEVICES=0 ./gpu/wave_bench_gpu \
  /scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/output.1aln \
  gpu/wave.seeds --stage2-fit --cpu-1core=1200.7 --cpu-32core=38713.8

# concurrency sweep (launch-size only; env TR_POOL) -> pool2048 = 0.31x, ceiling 2160 = 0.32x
TR_POOL=2048 STAGE2_NREPS=3 CUDA_VISIBLE_DEVICES=0 ./gpu/wave_bench_gpu ... --stage2-fit ...
TR_POOL=2160 STAGE2_NREPS=3 CUDA_VISIBLE_DEVICES=0 ./gpu/wave_bench_gpu ... --stage2-fit ...

# trace-point correctness, 100k uniform sample
CUDA_VISIBLE_DEVICES=0 ./gpu/wave_validate \
  /scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/output.1aln \
  gpu/wave.seeds --sample=100000
```
