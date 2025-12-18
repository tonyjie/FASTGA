# Iteration 3 Report

## Summary
- **Baseline wall time**: 13.570s
- **Iter_1 wall time (2MB)**: 13.344s (best)
- **Iter_3 wall time (4MB)**: 13.623s
- **Regression from iter_1**: 0.279s (2.1% slower)
- **Verification**: PASSED

## Code Change
Modified `FastGA.c` - increased IO_BUFFER_SIZE from 2MB to 4MB:
```c
// Before:
#define    IO_BUFFER_SIZE  2000000  // 2MB buffer

// After:
#define    IO_BUFFER_SIZE  4000000  // 4MB buffer
```

## Performance Breakdown

### Phase 1: Adaptive seed merge
- Iter_1: 4.877s wall, 701.0% CPU
- Iter_3: 4.757s wall, 723.4% CPU
- **Improvement: 2.5% faster** (minor improvement)

### Phase 2: Seed sort and alignment search
- Iter_1: 8.441s wall, 1318.7% CPU
- Iter_3: 8.835s wall, 1334.6% CPU
- **Regression: 4.7% slower**

### Total
- Iter_1: 13.344s wall (best)
- Iter_3: 13.623s wall
- **Net regression: 2.1% slower**

## Analysis
**FAILED OPTIMIZATION**: Increasing IO_BUFFER_SIZE from 2MB to 4MB made performance worse.

The 4MB buffer:
1. Slightly improved Phase 1 (fewer I/O syscalls)
2. **Significantly worsened Phase 2** - likely due to:
   - Increased memory pressure (8GB+ total buffer space with 32 threads)
   - L3 cache pollution from larger buffers
   - Diminishing returns from fewer syscalls vs. cache efficiency

**Finding**: 2MB buffer appears to be the optimal size, balancing syscall reduction with cache efficiency.

## Decision
**REVERT**: This change will be reverted. The 2MB buffer from iter_1 remains optimal.

## Summary of Buffer Size Experiments
| Buffer Size | Wall Time | Delta vs Baseline |
|-------------|-----------|-------------------|
| 1MB (baseline) | 13.570s | - |
| 2MB (iter_1) | 13.344s | -1.7% ✓ |
| 4MB (iter_3) | 13.623s | +0.4% ✗ |

## Next Steps
- Revert to 2MB buffer
- Try a different optimization approach (not buffer-related)
