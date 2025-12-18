# Iteration 0 Report

## Summary
- **Baseline wall time**: 13.570s
- **Iter_0 wall time**: 13.485s
- **Improvement**: 0.085s (0.6% faster)
- **Verification**: PASSED

## Code Change
Modified `RSDsort.c` - increased sorting thresholds to delay switch from radix sort to shell sort:
```c
// Before:
#define THR0 15
#define THR1 15
#define THR2  8
#define GAP1  9
#define GAP2  4

// After:
#define THR0 25
#define THR1 25
#define THR2 12
#define GAP1 13
#define GAP2  6
```

## Performance Breakdown

### Phase 1: Adaptive seed merge
- Baseline: 5.002s wall, 692.3% CPU
- Iter_0: 4.702s wall, 739.4% CPU
- **Improvement: 6% faster** (0.3s saved)

### Phase 2: Seed sort and alignment search
- Baseline: 8.550s wall, 1394.6% CPU
- Iter_0: 8.766s wall, 1270.4% CPU
- **Regression: 2.5% slower** (0.2s lost)

### Total
- Baseline: 13.570s wall, 1134.1% CPU
- Iter_0: 13.485s wall, 1083.8% CPU
- **Net improvement: 0.6% faster**

## Analysis
The increased sorting thresholds had a **mixed effect**:
1. **Phase 1 improved** - The higher THR0 threshold allows radix sort to continue longer before switching to shell sort, which benefited the seed merge phase.
2. **Phase 2 regressed slightly** - The larger gap sizes in shell sort (GAP1, GAP2) may have added overhead for the specific data patterns in the alignment search phase.

The overall improvement is modest (0.6%) but demonstrates that the sorting thresholds can impact performance. The next iteration should focus on more impactful optimizations.

## Next Steps
Consider:
1. Optimizing I/O buffer sizes in the merge threads
2. Improving thread load balancing
3. Reducing memory copies in hot paths
