# FastGA Optimization 2: Sparse Prefix Index Design

## Overview

This design implements Optimization 2 from the storage reduction plan: replacing FastGA's fixed 128 MB prefix index with a sparse two-level index to achieve ~252 MB storage savings with bit-exact output.

## Problem Statement

FastGA's current GIX indices consume 128 MB per genome for a fixed prefix index, regardless of actual k-mer distribution density. This 128 MB consists of 16M int64 entries covering all possible 12-byte (3×4-bit packed) k-mer prefixes. Most of these entries are likely sparse since k-mers aren't uniformly distributed across sequence space.

**Current Impact**:
- 128 MB × 2 genomes = 256 MB wasted on mostly-empty prefix indices
- Fixed overhead regardless of genome size or k-mer density
- Represents ~0.4% of total GIX storage for human genomes, but eliminates constant overhead

## Design Approach

**Selected Strategy**: Fixed 4-Byte Sparse Index (prioritizing maximum storage reduction over performance)

Replace the dense 12-byte prefix system with:
1. **Coarse Level**: 4-byte prefix index (256 entries, 2 KB)
2. **Fine Level**: Binary search within coarse ranges using full 12-byte comparison

This provides maximum storage reduction with predictable implementation complexity.

## Architecture

### Current System
```
.gix file structure:
├── Header (kmer size, parts, etc.)
├── Dense prefix index: 16M × int64 = 128 MB
└── Metadata (PostBytes, ContBytes, etc.)

Lookup process:
prefix_value = (kmer[0] << 16) | (kmer[1] << 8) | kmer[2]  // 12-byte prefix
position = prefix[prefix_value]  // O(1) direct access
```

### New System
```
.gix file structure:  
├── Header (kmer size, parts, format_version, etc.)
├── Sparse prefix index: 256 × int64 = 2 KB  
└── Metadata (PostBytes, ContBytes, etc.)

Lookup process:
coarse_key = kmer[0]  // 4-byte prefix (first byte)
start = (coarse_key == 0) ? 0 : sparse_index[coarse_key-1]
end = sparse_index[coarse_key]  
position = binary_search(ktab_range[start:end], full_12byte_prefix)  // O(log n)
```

## Component Implementation

### 1. GIXmake.c Changes (Index Creation)

**Target Functions**: `compress_thread()`, main index writing section

**Current Code**:
```c
// Line ~1228: prefix accumulation  
prefix += ((beg>>2) << 16);  // 16M entry addressing
prefix[idx] += w;            // idx from 12-byte prefix

// Line ~1343: allocation
prefix = (int64 *) MyBlock;
bzero(prefix,sizeof(int64)*0x1000000);  // 16M entries

// Line ~1563: writing to .gix
write(tab,prefix,sizeof(int64)*0x1000000);
```

**New Code**:
```c
// Use only first byte for sparse index
sparse_prefix += beg>>2;  // 256 entry addressing
sparse_prefix[kmer_first_byte] += w;  

// Allocation for 256 entries only
sparse_prefix = (int64 *) MyBlock;  
bzero(sparse_prefix,sizeof(int64)*0x100);  // 256 entries

// Write sparse index to .gix
write(tab,sparse_prefix,sizeof(int64)*0x100);
```

**Index Build Logic**:
- Change prefix increment from `((sarray[x+1] << 8) | sarray[x+2])` to `sarray[x+1]` 
- Reduce prefix array size from 16M to 256 entries
- Maintain cumulative transformation: `prefix[x] += prefix[x-1]`

### 2. libfastk.c Changes (Index Lookup)

**Target Functions**: `Find_Kmer_Index()`, `GoTo_Kmer_Entry()`, `More_Kmer_Stream()`

**Current Code**:
```c
// Line ~687: direct lookup
l = index[m-1];  // m is 12-byte prefix value
r = index[m];
```

**New Code**:
```c
// Two-stage lookup
coarse_key = extract_first_byte(kmer);
coarse_start = (coarse_key == 0) ? 0 : sparse_index[coarse_key-1];
coarse_end = sparse_index[coarse_key];

// Binary search within coarse range
position = binary_search_ktab(ktab, coarse_start, coarse_end, full_kmer_prefix);
```

**Binary Search Implementation**:
```c
static int64 binary_search_ktab(uint8 *ktab, int64 start, int64 end, 
                                uint8 *target_prefix, int prefix_bytes) {
    int64 left = start, right = end;
    while (left < right) {
        int64 mid = (left + right) / 2;
        uint8 *entry = ktab + mid * entry_size;
        int cmp = memcmp(entry, target_prefix, prefix_bytes);
        if (cmp < 0) left = mid + 1;
        else right = mid;
    }
    return left < end ? left : -1;  // -1 if not found
}
```

### 3. Format Versioning

**Backwards Compatibility**: Add format version field to distinguish sparse vs dense indices.

**Header Changes**:
```c
// In GIXmake.c - writing header
int format_version = 2;  // Version 2 = sparse prefix
write(tab, &format_version, sizeof(int));

// In libfastk.c - reading header  
int format_version;
read(f, &format_version, sizeof(int));
if (format_version == 1) {
    // Legacy: read 16M entries
} else if (format_version == 2) {
    // New: read 256 entries  
}
```

## Data Flow Changes

### Index Creation Flow (GIXmake)
1. **K-mer Processing**: During MSD radix sort, extract k-mers and increment `sparse_prefix[kmer[0]]++`
2. **Accumulation**: Convert counts to cumulative offsets: `sparse_prefix[i] += sparse_prefix[i-1]`
3. **File Writing**: Write 256 int64 entries (2 KB) instead of 16M entries (128 MB)

### Index Lookup Flow (libfastk)
1. **Target Extraction**: Extract first byte of target k-mer as coarse key
2. **Range Lookup**: Use coarse key to get ktab range `[start, end)`  
3. **Binary Search**: Search within range for exact 12-byte prefix match
4. **Position Return**: Return ktab position or -1 if not found

## Error Handling

**Invalid Lookups**: Binary search returns -1 for non-existent k-mers, maintaining same API contract
**File Format Errors**: Version check catches incompatible .gix files, fails gracefully with error message
**Memory Allocation**: Sparse index requires 256× less memory, reduces allocation failure risk

## Testing Strategy

### Bit-Exact Verification
1. **EXAMPLE Dataset**: Run full pipeline, compare output `.1aln` byte-for-byte with baseline
2. **Human Genomes**: Verify identical alignment results on GRCh38 vs CHM13  
3. **Self-Comparison**: Test symmetric mode to ensure correct handling

### Storage Measurement
1. **Before/After Comparison**: Measure .gix file sizes with `du -sb`
2. **Expected Savings**: 
   - EXAMPLE: 128 MB → 2 KB per genome = ~254 MB total savings
   - Human: 128 MB → 2 KB per genome = ~254 MB total savings
   - Percentage reduction: ~0.4% of total storage but eliminates fixed overhead

### Performance Baseline
1. **Seed Merge Timing**: Measure time for seed merge phase specifically
2. **Acceptance Criteria**: Any performance degradation acceptable (storage priority)
3. **Memory Usage**: Verify 256× reduction in prefix index memory footprint

## Implementation Plan

### Phase 1: Index Creation (GIXmake.c)
- [ ] Modify prefix array allocation: `0x1000000` → `0x100` 
- [ ] Update prefix increment logic to use first byte only
- [ ] Add format version writing to .gix header
- [ ] Update file write to output 256 entries

### Phase 2: Index Lookup (libfastk.c)  
- [ ] Implement binary search helper function
- [ ] Modify `Find_Kmer_Index()` for two-stage lookup
- [ ] Update `GoTo_Kmer_Entry()` and related functions
- [ ] Add format version reading and compatibility handling

### Phase 3: Integration & Testing
- [ ] Test on EXAMPLE dataset, verify bit-exact output
- [ ] Test on human genomes, verify bit-exact output  
- [ ] Measure storage savings and performance impact
- [ ] Update documentation with results

## Risk Assessment

**Low Risk**: 
- Isolated to prefix index logic, core k-mer storage unchanged
- Binary search is well-understood algorithm
- Format versioning enables rollback capability

**Medium Risk**:
- Performance degradation (acceptable given storage priority)
- Additional complexity in lookup path (mitigated by comprehensive testing)

## Success Criteria

1. **Storage Reduction**: Achieve ~252 MB savings (128 MB × 2 genomes) 
2. **Correctness**: Bit-exact `.1aln` output on both test datasets
3. **Compatibility**: Graceful handling of old .gix files via versioning
4. **Integration**: No changes to public APIs or workflow

## Future Enhancements

If further optimization needed:
- **Adaptive Granularity**: Choose prefix level based on actual k-mer distribution
- **Compressed Sparse Index**: Apply lightweight compression to 256-entry index
- **Profile-Guided Optimization**: Analyze real workload patterns for optimal coarse level

---

**Note**: This spec was written while user was unavailable for review. User should review and approve before proceeding to implementation planning.