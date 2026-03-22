# Storage Optimization Plan for FastGA

## Context

Benchmarking on human genomes (GRCh38 vs CHM13) confirmed that FastGA's peak disk usage is **71 GB** for a single comparison, with GIX indices at **62.7 GB (88%)** being the dominant cost. Running multiple comparisons on shared HPC nodes is often infeasible. We want to reduce storage requirements while preserving runtime performance and alignment quality.

## Evaluation Metrics

### 1. Storage (primary optimization target)
- Peak disk usage during run
- Persistent file sizes (GIX, GDB)
- Temp file peak (seed pairs, alignment temps)

### 2. Performance (must not significantly degrade)
- Wall clock time (measured via `/usr/bin/time -v`)
- Per-phase breakdown (from verbose log)
- Compare against baseline: 9m48s at T=32 for human genomes

### 3. Alignment Quality (two tiers)
- **Tier 1 — Bit-exact**: Output `.1aln` is identical to baseline. No quality evaluation needed. This is the standard for "zero-cost" optimizations.
- **Tier 2 — Comparable quality**: Output differs but alignment quality is acceptable. Evaluate using:
  - Genome coverage (% of bases covered by alignments)
  - Number of alignments and average length
  - Comparison against baseline using ALNshow statistics
  - For rigorous evaluation: simulated genomes with known-truth regions (paper's Section 5.1 methodology)

## Optimization Opportunities (ranked by impact and feasibility)

### Opportunity 1: Early GIX Deletion (Zero-cost, Bit-exact)

**What**: Delete GIX files immediately after seed merge completes, instead of at program exit.

**Why it works**: The `Kmer_Stream` (which reads `.ktab.*`) is already freed at the end of `adaptamer_merge()` (line 2483). The sort+align phase never touches GIX files — it only reads `.bps` and `_pair.*` temp files. Currently, GIX sits on disk through the entire 8-minute sort+align phase unnecessarily.

**Impact**: Does NOT reduce peak disk (peak occurs during seed merge when GIX is still needed), but frees ~63 GB during sort+align. Critical for concurrent runs: two human comparisons overlapping in sort+align would need 71+9=80 GB instead of 71+71=142 GB.

**Implementation** (in `FastGA.c`):
1. After `adaptamer_merge()` returns (~line 5154), free P1/P2 Post_Lists
2. Call `GIXrm -f` (not `-fg` — keep GDB since it's needed for sort+align) if `!KEEP && TYPE <= IS_GDB`
3. Set `TYPE1 = TYPE2 = IS_GDB + 1` so `Clean_Exit()` doesn't try to delete again
4. Move `Free_Post_List(P1/P2)` from end of main (~line 5283) to before GIXrm

**Risk**: Low. Need to verify `Close_GDB()` doesn't access `.1gdb` files (it likely only frees memory). Must handle SELF mode (P2==P1). Use `-f` not `-fg` to preserve GDB for safety.

**Verification**: Run on EXAMPLE + human datasets, confirm bit-exact `.1aln` output, measure storage timeline.

---

### Opportunity 2: Sparse Prefix Index (Bit-exact, Moderate effort)

**What**: Replace the 128 MB fixed `.gix` stub (16M int64 prefix index for all 12-bp prefixes) with a two-level index: coarse 8-bp level (64K entries, 512 KB) + fine-grain binary search.

**Impact**: Saves ~127 MB per genome (~254 MB total). Small fraction of total GIX but a constant overhead regardless of genome size.

**Risk**: Medium. The prefix index is used for O(1) lookups in the merge inner loop. Binary search adds latency. Need to benchmark whether the seed merge phase (currently 5.3s for human) slows down noticeably.

**Verification**: Bit-exact `.1aln` output (same seeds found, same alignments).

---

### Opportunity 3: Eliminate Mask Byte When Unused (Bit-exact, Low effort)

**What**: When no soft masks are applied (the common case — `-M` flag is off by default), skip the mask byte in ktab entries, reducing `MBYTES` from `KBYTES+1` to `KBYTES`.

**Impact**: Saves 1 byte per k-mer entry = ~7% of ktab size. For human genomes: ~2.2 GB per genome, ~4.4 GB total.

**Risk**: Low-medium. Requires changes in both GIXmake.c (write) and libfastk.c (read). The mask byte is checked during seed merge even when masking is off. Need to handle format versioning (old vs new ktab files).

**Verification**: Bit-exact `.1aln` output when no masks are used.

---

### Opportunity 4: On-the-fly LCP Computation (Bit-exact, Medium effort)

**What**: Remove the stored LCP byte from ktab entries. Recompute LCP during sequential scan in the seed merge by comparing consecutive k-mer suffixes.

**Impact**: Saves 1 byte per k-mer entry = ~7% of ktab. For human genomes: ~2.2 GB per genome. Combined with Opportunity 3: ~14% savings.

**Risk**: Medium. LCP computation in the inner loop adds CPU work. The merge already computes some LCP-like comparisons. Need to benchmark carefully — the seed merge is currently only 5.3s but is I/O-bound, so adding CPU work may not slow it down if I/O bandwidth is the bottleneck.

**Verification**: Bit-exact `.1aln` output (same LCP values computed, same seeds found).

---

### Opportunity 5: Streaming Ktab with Delta Compression (Bit-exact, High effort)

**What**: Delta-encode consecutive k-mer suffixes in ktab (they share long prefixes since they're sorted), then compress with a fast codec (LZ4/zstd). Decompress on-the-fly during the sequential merge scan.

**Impact**: Potentially 30-50% reduction in ktab I/O and disk size. For human genomes: ~10-16 GB per genome.

**Risk**: High complexity. Requires changes to GIXmake (write compressed), libfastk (read compressed), and careful buffering. Decompression CPU must not bottleneck the merge.

**Verification**: Bit-exact `.1aln` output.

---

### Opportunity 6: More Aggressive Syncmer Filtering (Trade-off, Low code effort)

**What**: Change syncmer parameters from (TMER=12, SMER=8) to e.g. (TMER=16, SMER=8), reducing the sampling rate from ~75% to ~44% of positions.

**Impact**: ~40% reduction in ktab + post-list size. For human genomes: ~13 GB per genome. Also reduces seed count, speeding up merge and sort+align.

**Risk**: **Alignment quality trade-off**. Fewer sampled positions means longer gaps between seeds, potentially missing short or divergent alignments. Requires Tier 2 evaluation (genome coverage, alignment count, sensitivity on simulated data).

**Verification**: Compare alignment coverage and sensitivity against baseline. NOT bit-exact.

---

## Recommended Implementation Order

### Phase A: Zero-cost win (Opportunity 1)
Implement early GIX deletion. No quality risk, no performance risk, immediate benefit for concurrent runs. Validates our code modification workflow on a simple change.

### Phase B: Format optimizations (Opportunities 3 + 4)
Eliminate mask byte and stored LCP when appropriate. ~14% ktab reduction with bit-exact output. Requires GIXmake + libfastk changes but the merge algorithm itself is unchanged.

### Phase C: Evaluate trade-off space (Opportunity 6)
Benchmark more aggressive syncmer filtering on EXAMPLE + human datasets. This requires building an evaluation framework first (alignment comparison scripts). If sensitivity loss is acceptable for the user's use case, this gives the biggest single reduction.

### Phase D: Advanced compression (Opportunity 5)
Only if Phase B savings are insufficient. High engineering effort.

## Immediate Next Steps

1. **Implement Opportunity 1** (early GIX deletion) — simplest, zero-risk
2. **Set up evaluation framework** — script to compare two `.1aln` files (alignment count, coverage, identical check)
3. **Benchmark Opportunity 1** on EXAMPLE and human datasets
4. **Then decide** whether to pursue Phase B or Phase C next based on results and priorities

## Files to Modify

- `FastGA.c` — Early GIX deletion (Opportunity 1)
- `GIXmake.c` — Format changes (Opportunities 3, 4)
- `libfastk.c` / `libfastk.h` — Reader changes (Opportunities 3, 4, 5)
- `benchmarks/` — New evaluation/comparison scripts

## Verification Strategy

For each optimization:
1. Run on EXAMPLE dataset, compare `.1aln` output (bit-exact for Tier 1)
2. Run on human genomes, compare `.1aln` output
3. Measure storage timeline (reuse `monitor_tmpdir.sh` and `df` polling)
4. Measure runtime (reuse `/usr/bin/time -v` + verbose log parsing)
5. Document results in `docs/benchmark_storage.md`
