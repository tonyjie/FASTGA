# Iteration 1 Report

## Summary
- **Baseline wall time**: 13.570s
- **Iter_0 wall time**: 13.485s
- **Iter_1 wall time**: 13.344s
- **Total improvement from baseline**: 0.226s (1.7% faster)
- **Improvement from iter_0**: 0.141s (1.0% faster)
- **Verification**: PASSED

## Code Change
Modified `FastGA.c` - increased I/O buffer size from 1MB to 2MB:
```c
// Added constant definition:
#define    IO_BUFFER_SIZE  2000000  // 2MB buffer for I/O operations (was 1MB)

// Replaced all occurrences of 1000000 (IO buffer related) with IO_BUFFER_SIZE
```

This reduces the number of write() and read() system calls by half.

## Performance Breakdown

### Phase 1: Adaptive seed merge
- Baseline: 5.002s wall, 692.3% CPU
- Iter_1: 4.877s wall, 701.0% CPU
- **Improvement: 2.5% faster** (0.125s saved)

### Phase 2: Seed sort and alignment search
- Baseline: 8.550s wall, 1394.6% CPU
- Iter_1: 8.441s wall, 1318.7% CPU
- **Improvement: 1.3% faster** (0.109s saved)

### Total
- Baseline: 13.570s wall, 1134.1% CPU
- Iter_1: 13.344s wall, 1090.6% CPU
- **Net improvement: 1.7% faster**

## Analysis
Increasing the I/O buffer size from 1MB to 2MB improved both phases:
1. **Phase 1 improved** - Fewer write() system calls during seed pair output
2. **Phase 2 improved** - Fewer read() calls during seed import for sorting

The 2x buffer size means approximately 2x fewer I/O system calls, which reduces kernel overhead. The memory cost is modest (doubling buffer allocation from ~2GB to ~4GB for 32 threads).

## Cumulative Changes
1. RSDsort.c thresholds (iter_0)
2. IO_BUFFER_SIZE 1MB -> 2MB (iter_1)

## Next Steps
Consider:
1. Further increase buffer sizes (4MB or 8MB)
2. Use O_DIRECT for bypassing kernel buffering
3. Optimize memory access patterns in hot loops
