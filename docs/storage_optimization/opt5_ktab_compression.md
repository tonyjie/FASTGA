# Optimization 5: Block-Compressed ktab (In Progress)

## Summary

Compress ktab partition files using zstd block compression for streaming decompression during the seed merge. Prototype validated ~14% compression ratio (1.16x).

| Property | Value |
|---|---|
| **Type** | Zero-cost (intended) |
| **Quality tier** | Tier 1 (Bit-exact, intended) |
| **Target** | Reduce ktab file size by ~14% |
| **Status** | **Blocked** — reader works, writer has unresolved bug |

## Compression Analysis

Tested on EXAMPLE dataset HAP1 ktab (740 MB, 12 bytes/entry after Opt3):

| Component | Size | zstd -9 | Ratio |
|---|---:|---:|---|
| K-mer suffix (7B/entry) | 56.6 MB | 44.0 MB | 1.3x (78%) |
| LCP (1B/entry) | 8.1 MB | 4.0 MB | 2.0x (50%) |
| Position payload (4B/entry) | 32.3 MB | 30.7 MB | **1.1x (95%)** — nearly incompressible |
| **Total per part** | **96.9 MB** | **83.5 MB** | **1.16x (86%)** |

The position payload (genomic coordinates) is essentially random data, limiting overall compression to ~14%.

## What Works

### Reader (libfastk.c) — Verified Correct
- Detects compressed format via negative `nels` in ktab header
- `More_Kmer_Stream`: decompresses one block per call, matching STREAM_BLOCK size
- `GoTo_Kmer_Index`: computes target block, decompresses, positions within block
- `Clone_Kmer_Stream`: allocates separate compression buffers per clone
- `Free_Kmer_Stream`: frees compression buffers
- **Verified**: produces 51,082,720 seeds (identical to baseline) when reading uncompressed ktab

### Compression Script (benchmarks/compress_ktab.py) — Verified Correct
- Post-processing tool that compresses existing uncompressed ktab files
- Produces identical decompressed bytes (verified byte-for-byte)
- Block format: negative nels header + block index + zstd compressed blocks

### Compressed Format
```
Header:
  int    KMER              (same as before)
  int64  nels              (NEGATIVE = -actual_nels, signals compressed)
  int    block_entries     (entries per block = 131072 = STREAM_BLOCK)
  int    num_blocks
  int64  block_offsets[num_blocks+1]

Body:
  [zstd compressed block 0]
  [zstd compressed block 1]
  ...
```

## What's Blocked

### Writer (GIXmake.c) — Has Bug
When the compressed output code is added to GIXmake.c, the ktab content changes even though:
- The memcpy reads from the same sarray addresses as the original big_write
- The byte counts are identical
- The sort/compress phases are not modified

**Symptoms**: Seed count changes from 51.1M to 53.6M (5% increase). Average seed length decreases from 28.5 to 25.9. Non-redundant alignments decrease from 323,569 to 309,494.

**Extensive debugging performed**:
1. Reader verified correct: gives 51.1M on uncompressed ktab (both with and without Opt5 libfastk.c changes)
2. Linking with `-lzstd` alone doesn't change GIXmake output
3. GIXmake is deterministic (two runs of same binary produce identical ktab)
4. `memcpy` collection (no compression) produces identical data — confirmed 51.1M
5. `ZSTD_compress` entire part_data at once and discard — no corruption, 51.1M
6. `ZSTD_compress` in blocks and discard — no corruption, 51.1M
7. Dual output test: uncompressed "verify" files written alongside compressed files from same `part_data` → verify files give 51.1M, compressed files give 53.6M
8. This means the C reader's decompression produces different bytes from the original data, even though Python's zstd decompression of the same compressed files matches

**Suspected root cause**: Mismatch between how the C reader (libfastk.c) decompresses blocks and how GIXmake wrote them. The block boundaries, entry counts per block, or the decompressed buffer size calculation may be subtly off, causing entries to be misaligned after decompression. The next debugging step: add byte-level verification inside the C reader to compare decompressed data against known-good data.

**Possible workaround**: Use `compress_ktab.py` as a post-processing step (build uncompressed, then compress externally, then read with the new libfastk.c reader). This approach was verified to produce correct decompressed data via Python.

## Estimated Impact (if completed)

| Genome Size | Current ktab (Opt3) | With Compression | Savings |
|---|---:|---:|---:|
| EXAMPLE (86 Mbp, per genome) | 740 MB | 638 MB | 102 MB (14%) |
| Human (3.1 Gbp, per genome) | ~30 GB | ~25.8 GB | ~4.2 GB |
| Human (both genomes) | ~58 GB | ~50 GB | ~8 GB |

## Git Info

Not committed — changes were reverted to keep codebase at Opt1+3. The reader code (libfastk.c) and compression script (compress_ktab.py) are available but not staged.

## Next Steps

1. Debug the GIXmake writer — likely a subtle variable scope or compiler optimization issue
2. Alternative: implement as a post-processing step using compress_ktab.py
3. If writer bug is resolved, verify bit-exact output and measure performance

## Checklist

- [x] Compression ratio validated (14% with zstd)
- [x] Reader implemented and verified correct
- [x] Compression script (compress_ktab.py) working
- [ ] Writer bug resolved
- [ ] Bit-exact output verified end-to-end
- [ ] Performance measured
- [ ] Documentation complete
- [ ] Git commits recorded
