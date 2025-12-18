# Iteration 5 Report

## Performance
- **Wall Time**: 16.190s (Regression).
- **Phase 1 Time**: 10.716s (Regression vs 4.8s in Iter 2).
- **Correctness**: Verified (PASS)

## Changes
- Reduced Phase 1 (Adaptive seed merge) threads from 32 to 4.
- Hypothesis was that 32 threads caused system contention.

## Analysis
- Phase 1 system time dropped slightly (30s -> 25.6s).
- But Wall time INCREASED significantly (4.8s -> 10.7s).
- This proves that parallelism WAS helping in Phase 1, despite the high system time overhead. The system time is likely per-thread overhead (syscalls), not contention on a single lock.
- Reducing threads reduced total throughput.

## Conclusion
- Phase 1 benefits from high thread count.
- Revert this change.
- Iteration 2 (9.688s) remains the best.
