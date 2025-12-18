# Iteration 0 Report

## Performance
- **Baseline Wall Time**: 13.084s
- **Current Wall Time**: 10.338s
- **Speedup**: 21% (Phase 2 speedup: 35%)
- **Correctness**: Verified (PASS)

## Changes
- Implemented a **Worker Pool** (Producer-Consumer model) for `search_seeds`.
- Previously, `search_seeds` was parallelized by partitioning the contigs using `rmsd_sort`. Due to the small number of contigs in the input, `rmsd_sort` produced few partitions (likely 2-6), leaving most of the 32 threads idle.
- The new implementation uses the existing partitions to define "Producer" threads. These threads identify seed groups (which must be processed sequentially per group) and push "Alignment Jobs" to a shared queue.
- All 32 threads act as "Consumers" (Workers), popping jobs from the queue and executing `align_contigs`.
- This ensures that the heavy lifting (`align_contigs`) is distributed across all available threads, regardless of the number of contigs/partitions.

## Analysis
- Phase 2 time dropped from 8.3s to 5.3s.
- CPU utilization in Phase 2 increased from ~1300% to ~2300%.
- Phase 1 (Adaptive seed merge) remains a bottleneck (4.9s wall, 30s system). It has high system time overhead.

## Next Steps
- Optimize Phase 1 to reduce system time overhead (likely due to frequent IO or small reads).
