# Iteration 2 Report

## Performance
- **Baseline Wall Time**: 13.084s
- **Iter 0 Wall Time**: 10.338s
- **Current Wall Time**: 9.688s
- **Speedup**: 26% vs Baseline.
- **Correctness**: Verified (PASS)

## Changes
- **Worker Pool**: Reverted to Single Item Queue (from Iter 0), rejecting Batching (Iter 1).
- **Buffer Size**: Increased IO Buffers in Phase 1 (`new_merge_thread`) from 1MB to 8MB.
- **Queue Size**: Increased Global Queue size to 1,000,000 items.

## Analysis
- Phase 2 time improved to 4.796s (vs 5.333s in Iter 0 and 8.4s in Iter 1).
- CPU Utilization in Phase 2: ~2275% (23 threads).
- Phase 1 time remained ~4.8s. System time is high (30s). Larger buffers didn't significantly reduce system time (maybe frequent `read` calls inside `Kmer_Stream` library, which uses `STREAM_BLOCK` 128KB, are the bottleneck, and I didn't change `STREAM_BLOCK`).

## Next Steps
- Try to improve Phase 2 utilization further (to >30 threads) by reducing lock contention in the Worker Queue using **Work Stealing** with multiple queues.
