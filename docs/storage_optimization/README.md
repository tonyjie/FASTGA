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

| # | Optimization | Storage Impact | Quality | Status |
|---|---|---|---|---|
| 1 | [Early GIX Deletion](opt1_early_gix_deletion.md) | Frees ~1,860 MB during sort+align (EXAMPLE) / ~63 GB (human) | Bit-exact | Implemented, verified |

## Planned Optimizations

| # | Optimization | Estimated Impact | Quality | Status |
|---|---|---|---|---|
| 2 | Sparse prefix index | ~254 MB savings | Bit-exact | Not started |
| 3 | Eliminate mask byte | ~7% ktab reduction | Bit-exact | Not started |
| 4 | On-the-fly LCP | ~7% ktab reduction | Bit-exact | Not started |
| 5 | Streaming compression | 30-50% ktab reduction | Bit-exact | Not started |
| 6 | Aggressive syncmer filtering | ~40% ktab reduction | Trade-off | Not started |
