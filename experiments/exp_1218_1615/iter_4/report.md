# Iteration 4 Report

## Performance
- **Wall Time**: 15.336s (Regression).
- **Phase 1 Time**: 10.420s (Regression).
- **Correctness**: Verified (PASS)

## Changes
- Increased `STREAM_BLOCK` to 0x40000 (256KB items, 4MB reads) in `libfastk.c`.
- Added `posix_fadvise(..., POSIX_FADV_SEQUENTIAL)` in `libfastk.c`.

## Analysis
- Phase 1 system time increased from ~30s (Iter 2) to ~90s.
- Large read buffers (4MB) and/or `posix_fadvise` caused performance degradation, likely due to memory pressure or cache thrashing.
- Reverted `libfastk.c`.

## Conclusion
- Default `STREAM_BLOCK` (128KB items, 2MB reads) seems optimal for this system.
