# Iteration 2 Report

## Summary
- **Baseline wall time**: 13.570s
- **Iter_1 wall time**: 13.344s (best)
- **Iter_2 wall time**: 13.556s
- **Regression from iter_1**: 0.212s (1.6% slower)
- **Verification**: PASSED

## Code Change
Modified `FastGA.c` - increased POST_BLOCK size from 128KB to 256KB:
```c
// Before:
#define POST_BLOCK 0x20000

// After:
#define POST_BLOCK 0x40000  // 256KB (was 128KB) - larger reads from index files
```

## Performance Breakdown

### Phase 1: Adaptive seed merge
- Iter_1: 4.877s wall, 701.0% CPU
- Iter_2: 4.796s wall, 713.7% CPU
- **Improvement: 1.7% faster** (minor improvement)

### Phase 2: Seed sort and alignment search
- Iter_1: 8.441s wall, 1318.7% CPU
- Iter_2: 8.737s wall, 1332.7% CPU
- **Regression: 3.5% slower**

### Total
- Iter_1: 13.344s wall, 1090.6% CPU
- Iter_2: 13.556s wall, 1111.6% CPU
- **Net regression: 1.6% slower**

## Analysis
**FAILED OPTIMIZATION**: Increasing POST_BLOCK from 128KB to 256KB made performance worse.

The larger POST_BLOCK:
1. Slightly improved Phase 1 (fewer reads from index files)
2. **Significantly worsened Phase 2** - likely due to:
   - Increased memory pressure from larger buffers per thread
   - Cache pollution from reading more data than needed
   - The 128KB size was already well-tuned for the access patterns

## Decision
**REVERT**: This change will be reverted before the next iteration since it worsened overall performance.

## Next Steps
- Revert POST_BLOCK to 0x20000 (128KB)
- Try different optimization approaches:
  - Optimize hot loop memory access patterns
  - Use compiler hints like __builtin_prefetch
  - Try different IO_BUFFER_SIZE values
