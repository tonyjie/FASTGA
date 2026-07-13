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

## Recommended next step
Map decompression to the GPU together with the wave: keep each genome's `.bps` (2-bit) resident,
add an on-device unpack (or 2-bit-aware sequence access in the wave/discovery kernels), and drop
`gpu_load_seqs`' host-side NUMERIC copy. Re-profile to confirm the 29% is gone. Then the batching
redesign (A3b) operates on a pipeline whose two dominant kernels are both on the GPU.
