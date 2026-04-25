# Storage Optimization Efforts

This directory tracks our efforts to reduce FastGA's disk storage requirements.

## Baseline

Measured on EXAMPLE dataset (~86 Mbp per genome) and human genomes (GRCh38 vs CHM13, ~3.1 Gbp each):

| Metric | EXAMPLE (T=32) | Human (T=32) |
|---|---:|---:|
| GIX (both genomes) | 1,864 MB | 62.7 GB |
| Seed pair temps (peak) | 438 MB | ~7 GB |
| Peak total disk | 2,344 MB | 71 GB |
| Sort+align disk (baseline) | ~1,905 MB | ~64 GB |

GIX indices are 88% of peak storage. See `../benchmark_storage.md` for full analysis.

## Evaluation Metrics

- **Storage**: Peak disk, persistent files, temp file peak
- **Performance**: Wall clock time, per-phase breakdown
- **Alignment Quality**:
  - **Tier 1 (Bit-exact)**: Output `.1aln` identical to baseline — no quality evaluation needed
  - **Tier 2 (Comparable)**: Evaluate genome coverage, alignment count/length vs baseline

## Optimization Plan

See [optimization_plan.md](optimization_plan.md) for the full plan with 6 identified opportunities.

## Completed Optimizations

| # | Optimization | Storage Impact | Quality | Status | Commits |
|---|---|---|---|---|---|
| 1 | [Early GIX Deletion](opt1_early_gix_deletion.md) | Frees ~1,860 MB during sort+align (EXAMPLE) / ~63 GB (human) | Bit-exact | Implemented, verified | Docs: `049e311`, Code: `77ac29d` |
| 3 | [Eliminate Mask Byte](opt3_eliminate_mask_byte.md) | -7.7% ktab size = -124 MB (EXAMPLE) / -4.6 GB (human) | Bit-exact | Implemented, verified | Docs: `6460e0a`, Code: `ea31b73` |
| 7 | [Chunk-wise GIX (single-sided)](opt7_chunkwise_gix.md) | -36% peak (EXAMPLE) / projected -37% (human) | Bit-exact (after Opt 8 buck-bzero fix) | Implemented; alignment-level verification only after Opt 8 fix | Code: `685a071` |
| 8 | [Bilateral Chunking](opt8_bilateral_chunking.md) | **-60.2% peak (1,810 → 721 MB)** on EXAMPLE; projected **-69%** (human, ~20 GB peak) | Bit-exact (51,082,720 seeds + 323,569 alignments verified) | Implemented, verified end-to-end | Uncommitted (working tree) |

## Cumulative Effect

With Opt 1 + Opt 3 combined:

| Metric | Baseline | With Opt 1+3 | Savings |
|---|---:|---:|---|
| GIX size (EXAMPLE, both) | 1,864 MB | 1,740 MB | -124 MB (-6.6%) |
| Peak total (EXAMPLE) | 2,344 MB | ~2,220 MB | -124 MB (-5.3%) |
| Sort+align disk (EXAMPLE) | ~1,905 MB | ~438 MB | -1,467 MB (-77%) |
| GIX size (human, both) | 62.7 GB | ~58.1 GB | -4.6 GB (-7.3%) |
| Peak total (human) | 71 GB | ~66.4 GB | -4.6 GB (-6.5%) |
| Sort+align disk (human) | ~64 GB | ~7 GB | -57 GB (-89%) |

## Planned Optimizations

| # | Optimization | Estimated Impact | Quality | Status |
|---|---|---|---|---|
| 2 | Sparse prefix index | ~254 MB savings | Bit-exact | Not started |
| 4 | On-the-fly LCP | -8.3% ktab = -127 MB (EXAMPLE) / -4.5 GB (human) | Bit-exact | **Verified** — full pipeline exact match (seeds + alignments) |
| 5 | [Ktab Compression](opt5_ktab_compression.md) | ~14% ktab reduction (not 30-50% as projected) | Bit-exact (intended) | **FAILED** — decompression round-trip bug, reverted |
| 6 | Aggressive syncmer filtering | ~40% ktab reduction | Trade-off | Not started |

## Comprehensive Report

See [optimization_report.md](optimization_report.md) for a detailed analysis of all attempted optimizations, rationale, results, and recommendations.
