# Iteration 4 Report

## Summary
- **Baseline wall time**: 13.570s
- **Iter_1 wall time**: 13.344s
- **Iter_4 wall time**: 13.341s (new best!)
- **Total improvement from baseline**: 0.229s (1.7% faster)
- **Improvement from iter_1**: 0.003s (0.02% faster)
- **Verification**: PASSED

## Code Changes

### 1. Reverted iter_2 and iter_3 regressions
- IO_BUFFER_SIZE restored to 2MB (optimal)
- POST_BLOCK restored to 0x20000 (128KB)

### 2. Modified Makefile compiler flags
```makefile
# Before:
CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

# After:
CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing -march=native -funroll-loops
```

**New flags:**
- `-march=native`: Optimize for the current CPU architecture
- `-funroll-loops`: Unroll loops for potentially better performance

## Performance Breakdown

### Phase 1: Adaptive seed merge
- Baseline: 5.002s wall, 692.3% CPU
- Iter_4: 4.836s wall, 702.5% CPU
- **Improvement: 3.3% faster**

### Phase 2: Seed sort and alignment search
- Baseline: 8.550s wall, 1394.6% CPU
- Iter_4: 8.481s wall, 1354.7% CPU
- **Improvement: 0.8% faster**

### Total
- Baseline: 13.570s wall, 1134.1% CPU
- Iter_4: 13.341s wall, 1116.1% CPU
- **Net improvement: 1.7% faster**

## Final Summary of All Iterations

| Iteration | Wall Time | Delta vs Baseline | Key Change |
|-----------|-----------|-------------------|------------|
| Baseline | 13.570s | - | Original code |
| Iter_0 | 13.485s | -0.6% | RSDsort thresholds |
| Iter_1 | 13.344s | -1.7% ✓ | IO_BUFFER_SIZE 2MB |
| Iter_2 | 13.556s | -0.1% ✗ | POST_BLOCK 256KB (reverted) |
| Iter_3 | 13.623s | +0.4% ✗ | IO_BUFFER_SIZE 4MB (reverted) |
| **Iter_4** | **13.341s** | **-1.7%** ✓ | Compiler flags |

## Cumulative Optimizations (Final Code)
1. **RSDsort.c**: THR0=25, THR1=25, THR2=12, GAP1=13, GAP2=6
2. **FastGA.c**: IO_BUFFER_SIZE = 2MB
3. **Makefile**: Added `-march=native -funroll-loops`

## Analysis
The compiler optimizations provided a small additional improvement:
- `-march=native` allows GCC to use CPU-specific instructions (e.g., AVX2, SSE4.2)
- `-funroll-loops` reduces loop overhead in hot paths

The combined optimizations achieved:
- **1.7% overall improvement** (13.570s → 13.341s)
- **3.3% improvement in Phase 1** (seed merge)
- **0.8% improvement in Phase 2** (alignment search)

## Recommendations for Further Optimization
1. Profile-guided optimization (PGO) could provide additional gains
2. Consider SIMD vectorization of key loops
3. Memory alignment of frequently accessed structures
4. Thread affinity binding for NUMA systems
