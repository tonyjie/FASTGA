# Optimization Summary

## Goal
Optimize FastGA algorithm on multi-threaded CPUs.

## Baseline
- Time: **13.084s**
- CPU Utilization: ~1100% (of 3200% possible).
- Bottleneck: Phase 2 "Seed sort and alignment search" (8.3s) had low utilization. Phase 1 "Adaptive seed merge" (4.7s) had high system time.

## Iteration 0: Worker Pool
- **Strategy**: Implemented a Producer-Consumer (Worker Pool) model for `search_seeds`.
- **Problem**: Original code partitioned work by contigs (`rmsd_sort`). For small numbers of contigs (e.g. Human genome), this resulted in few partitions (<32), limiting thread usage.
- **Solution**: Producers (assigned to partitions) identify alignment candidates (groups) and push them to a shared Queue. All 32 threads act as Consumers, processing alignments in parallel.
- **Result**: Time **10.338s**. Speedup 21%. Phase 2 time dropped to 5.3s.

## Iteration 1: Batching (Failed)
- **Strategy**: Tried to batch jobs in the Worker Queue to reduce lock contention.
- **Result**: Time **13.2s** (Regression).
- **Reason**: Batching likely introduced load imbalance or overhead that outweighed locking benefits.

## Iteration 2: Tuned Worker Pool (Best)
- **Strategy**: Reverted to Single Item Queue (Iter 0 logic). Increased IO Buffer size in Phase 1 (8MB) to reduce system calls. Increased Queue capacity.
- **Result**: Time **9.688s**. Speedup 26%.
- **Details**: Phase 2 time **4.796s** (vs 8.3s baseline). Utilization ~2300%. Phase 1 time remained ~4.8s.

## Iteration 3: Large Read Buffers (Failed)
- **Strategy**: Increased `STREAM_BLOCK` in `libfastk.c` to 1MB items (10MB+).
- **Result**: Time **92s** (Massive Regression). System time exploded.
- **Reason**: Memory thrashing or inefficiency with very large reads.

## Final Result
- **Time**: 9.688s
- **Speedup**: 26%
- **Code Changes**: Modified `FastGA.c` to implement Worker Pool for `search_seeds` and increased IO buffers.
