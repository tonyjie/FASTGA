# Comprehensive hotspot profile — FastGA align run (human, T=32)

`perf record -F 499` over the full FastGA align run (GRCh38 × CHM13, pre-built indices, so
this covers seed-merge + sort+align = the 498 s run, not GDB/GIX build). Flat self-time:

| % self | function | what it is | GPU fit |
|---:|---|---|---|
| 31.2% | `reverse_wave` | align wave (reverse) | ✓ mapped (A1/A2) |
| 28.8% | `forward_wave`  | align wave (forward) | ✓ mapped |
| **29.4%** | **`Uncompress_Read`** | **2-bit DNA → 1-byte/base unpack** (`Get_Contig`) | **excellent (unmapped)** |
| 4.0% | `new_merge_thread` | adaptamer seed merge | I/O + merge; poor |
| 2.4% | `radix_sort` | reverse-MSD seed sort | good, but small |
| 1.5% | `align_contigs` | chain scan / tube walk | control; stays CPU |
| 0.9% | `reimport_thread` | seed reimport/transform | small |
| ~0%  | `Compress_TraceTo8`, `fwrite`, `oneWriteLine` | output-merge / write | negligible |

Wave (forward+reverse) = **60.0%**, cross-validating the WAVE_TIMING instrumentation (61%).

## The finding: the #2 kernel is sequence DECOMPRESSION, not sort/output
`Uncompress_Read` is **29%** — bigger than everything except the wave. My earlier "sort+output
is the other ~40%" was wrong: radix sort is only 2.4% and output is ~0%. The non-wave time is
dominated by **decompressing the 2-bit genome into the aligner's byte buffer**.

Why so large: `align_contigs` calls `Get_Contig` whenever the A- or B-contig changes; the
B-contig is re-decompressed once per A-contig it pairs with. FastGA streams (low memory) rather
than caching decompressed contigs, so the same chromosome is unpacked many times.

## GPU-mapping assessment
- **Wave (60%)** — mapped; warp-cooperative kernels (A1/A2), measured fast.
- **`Uncompress_Read` (29%)** — *embarrassingly parallel* (independent 2-bit unpacks,
  bandwidth-bound). Two ways it disappears:
  1. **Keep the genome 2-bit-resident on the GPU (780 MB each) and unpack on-device** (or read
     2-bit directly in the wave kernel). The GPU aligner needs the sequences resident anyway, so
     this removes BOTH the CPU decompression AND the repeated H2D of decompressed bytes.
  2. (CPU-only quick win, orthogonal) cache decompressed contigs so each is unpacked once — trades
     ~3-12 GB RAM for eliminating the redundant re-decompression.
- **radix_sort (2.4%), merge (4%)** — small; radix sort is GPU-friendly (cub::DeviceRadixSort) but
  low priority given its size. Output-merge is I/O-bound, poor GPU fit, but ~0% anyway.

## What this does to the Amdahl ceiling
Wave-only offload addresses 60% → end-to-end ceiling ~1/(1-0.49)=~1.9× (wave is ~49% of the
run after the ~18% GDB/GIX + non-wave). **Wave + on-device decompression addresses ~89%** of the
align run → the ceiling jumps to ~1/(1-0.89) ≈ **~9×** on the align-only run (before GDB/GIX,
which is a separate one-time cost). This is the real reason to map the second kernel: it is the
difference between a ~2× and a ~5-9× system.

## On-device decompression — IMPLEMENTED + measured (2026-07-13)

`gpu/decomp_bench.cu` (`make decomp_bench`) + `unpack2bit` / `gpu_load_seqs_2bit` in the
library. 3 Gbp (human-scale), output validated byte-identical to the `Uncompress_Read` formula:

| | rate |
|---|---:|
| CPU 1-thread | 1.14 Gbase/s |
| CPU 32-thread (single pass) | 16.1 Gbase/s (0.186 s) |
| **GPU kernel** | **257 Gbase/s (0.012 s) — 16× CPU-32** |
| GPU incl. 2-bit H2D (750 MB) | 68 Gbase/s |

The deeper point: a *single* 3 Gbp CPU pass is only 0.186 s, yet `Uncompress_Read` is 29%
(~140 s) of the run — so the cost is **redundant re-decompression** (the same contig unpacked
per contig-pair). Keeping each genome 2-bit-resident on the GPU (750 MB) and unpacking once
(12 ms) — via `gpu_load_seqs_2bit`, which uploads the packed `.bps` and unpacks on-device
instead of the host copying pre-decompressed NUMERIC — makes the 29% vanish. The library now
supports this; wiring it into the `-G` path is a one-line swap (`Get_Contig(...,COMPRESSED)` +
`gpu_load_seqs_2bit` instead of `Get_Contig(...,NUMERIC)` + `gpu_load_seqs`), best done as part
of A3b so the resident 2-bit genome is shared across the whole batched pipeline.

## CPU decompression cache (`-DSEQ_CACHE`) — implemented + measured

`FastGA.c` under `#ifdef SEQ_CACHE`: decompress every B-contig once into a shared read-only
array (`build_bcache`), and point `align->bseq` at it instead of re-`Get_Contig`-ing per
contig-pair. (A is already loaded once per A-contig; B is the redundant one; B is never
complemented in the main path, so one forward cache serves both strands.)

Human GRCh38 × CHM13, T=32 — **per-alignment output byte-identical (md5 match, 518,037 records)**:
| | wall | CPU-time (user) | RSS |
|---|---:|---:|---:|
| stock | 494.7 s | 2584.7 CPU-s | 19.9 GB |
| cache | 448.2 s | 1761.2 CPU-s | 15.3 GB |
| **ratio** | **1.10×** | **1.47×** | — |

**Key nuance:** removing the 29% decompression cuts CPU-time 1.47× but wall only **1.10×** —
because decompression is *parallel* work spread across cores, while the wall is bound by the
wave's ~6-core critical path (the wave is the least-parallel phase). So the decompression cache
is an energy/CPU-efficiency win (and a real one, worth upstreaming), but **the wave is the wall
lever** — which is exactly what the GPU offload targets. (A GPU pipeline still eliminates the
decompression by keeping the genome 2-bit-resident, per `gpu_load_seqs_2bit`; it just isn't the
main wall win.)

## Recommended next step
Map decompression to the GPU together with the wave: keep each genome's `.bps` (2-bit) resident,
add an on-device unpack (or 2-bit-aware sequence access in the wave/discovery kernels), and drop
`gpu_load_seqs`' host-side NUMERIC copy. Re-profile to confirm the 29% is gone. Then the batching
redesign (A3b) operates on a pipeline whose two dominant kernels are both on the GPU.
