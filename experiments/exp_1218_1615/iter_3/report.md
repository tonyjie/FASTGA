# Iteration 3 Report

## Performance
- **Wall Time**: 1m 32s (92s).
- **Result**: FAILED (Huge regression in Phase 1).
- **Correctness**: Verified (PASS)

## Changes
- Tried to increase `STREAM_BLOCK` in `libfastk.c` from 128KB (items) to 1MB (items).
- This aimed to reduce `read` syscalls in Phase 1.

## Analysis
- Phase 1 system time exploded to 11 minutes (cpu time) / 86s (wall time).
- This suggests massive memory thrashing or inefficient large I/O handling.
- Reverted this change.

## Conclusion
- Iteration 2 remains the best configuration.
