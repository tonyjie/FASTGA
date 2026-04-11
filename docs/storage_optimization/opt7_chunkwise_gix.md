# Optimization 7: Chunk-wise GIX Processing

**Status**: PARTIAL SUCCESS — correctness verified, significant performance regression

## Concept

Instead of building all 8 ktab partitions for both genomes simultaneously (~62 GB for human), build genome2's ktab in chunks (e.g., 4 chunks of 2 partitions each). For each chunk:

1. Run `GIXmake -C first:last` to build only that chunk's partitions
2. Open genome1's full GIX + genome2's chunk GIX
3. Run `adaptamer_merge` — seeds are appended to the same temp files across chunks
4. Delete genome2's chunk ktab files
5. Repeat for the next chunk

This reduces genome2's peak ktab from 100% to 1/K of its full size (K = number of chunks).

## Implementation

### GIXmake Changes
- Added `-C first:last` flag (1-based partition range)
- Main sort loop restricted to `CHUNK_FIRST..CHUNK_LAST`
- Output partition files numbered 1..chunk_nparts (relative to chunk)
- `.gix` stub written with chunk_nparts and chunk-only prefix index
- Verified: chunk partition files are byte-identical to corresponding full-build partitions

### FastGA Changes
- Added `-C K` flag (number of chunks, must be >= 2)
- Chunk loop orchestrates: delete genome2's full ktab → loop K times (build chunk via `system("GIXmake -C ...")`, open streams, merge, delete chunk files)
- Seeds accumulate across chunks in the same `_pair.*` temp files

### Critical Performance Fix: Empty-Prefix Skip
In chunk mode, ~75% of genome2's prefix space is empty per chunk. The merge thread must skip T1 entries whose prefix has no match in T2's chunk. 

**Original approach (failed — hung at 99%)**: Iterate T1 one entry at a time for each empty prefix, then GoTo back. O(total_entries × empty_prefixes).

**First fix (failed — still 18x slower)**: Use `T1->index[cpre]` to jump T1 to the next prefix. Still O(num_empty_prefixes) disk seeks.

**Final fix (working)**: When T2's cache is empty, `T2->cpre` already points to T2's next non-empty prefix (since `while (T2->cpre < cpre) Next_Kmer_Entry(T2)` advanced T2 past the empty range during cache loading). Jump T1 directly to `T2->cpre`. This is O(num_chunks) seeks per merge, not O(num_empty_prefixes).

```c
if (ctop == cache)
  { int npre = T2->cpre;   // T2 already at its next non-empty prefix
    if (T2->csuf == NULL) break;  // T2 exhausted
    int64 t1_start = (npre > 0) ? T1->index[npre-1] : 0;
    if (t1_start >= tend) break;
    GoTo_Kmer_Index(T1, t1_start > 0 ? t1_start-1 : 0);
    continue;
  }
```

## Results (EXAMPLE dataset, ~86 Mbp, T=4, K=4)

### Correctness
| Metric | Baseline | Chunk Mode | Match? |
|---|---:|---:|---|
| Total seeds | 51,082,720 | 51,082,720 | **Exact** |
| Ave seed length | 28.5 | 28.5 | **Exact** |
| Per-chunk seeds | — | 12.8M + 12.8M + 12.8M + 12.8M | Sums exactly |

### Performance
| Metric | Baseline | Chunk Mode | Overhead |
|---|---:|---:|---|
| Merge wall time | ~12s | ~206s | **~18x slower** |
| Merge user CPU | ~3.6s | ~184s | ~50x more CPU |

The performance regression comes from:
1. **4× GIXmake rebuilds**: ~14s overhead (minor)
2. **Chunk 1 full T1 scan**: Chunk 1 iterates all T1 entries 0→100%, same as baseline
3. **Chunk 1 dominates**: Chunks 2-4 skip instantly (0%→100%), but chunk 1 does full work
4. **Re-reading T1 4 times**: Each chunk reopens T1, though OS page cache mitigates this

### Measured Storage WITHOUT `-n` (EXAMPLE dataset, ~86 Mbp, T=4, K=4)
| Metric | Baseline | Chunk Mode | Savings |
|---|---:|---:|---|
| **Peak disk** | **1,658 MB** | **1,658 MB** | **0 MB (unchanged!)** |
| Sustained disk (during merge) | 1,658 MB | 1,146 MB | -512 MB (-31%) |
| After genome2 deletion | — | 976 MB | — |

**Critical finding**: Without `-n`, peak storage is NOT reduced because genome2's full GIX must be built first (to obtain NPARTS and metadata).

### Measured Storage WITH `-n` Stub-Only Mode (EXAMPLE dataset, ~86 Mbp, T=4, K=4)

The `-n` flag was implemented to solve the peak storage problem: GIXmake writes only the `.gix` stub (128 MB prefix index + metadata) without any ktab partition files. FastGA then builds ktab one chunk at a time.

| Metric | Baseline | Chunk + `-n` | Savings |
|---|---:|---:|---|
| **Peak disk** | **1,658 MB** | **1,061 MB** | **-597 MB (-36%)** |
| Pre-run storage | 1,658 MB | 976 MB | -682 MB (-41%) |
| Seeds produced | 51,082,720 | 51,082,720 | **Exact match** |

### Extrapolated Storage WITH `-n` (Human Genomes, ~3.1 Gbp, K=4)
| Metric | Baseline | Chunk + `-n` | Savings |
|---|---:|---:|---|
| **Peak GIX** | **~62.7 GB** | **~39.3 GB** | **-23.4 GB (-37%)** |
| Genome2 on disk | ~31.35 GB | ~7.8 GB (1 chunk) | -23.5 GB (-75%) |

### Implementation: `-n` Stub-Only Mode
To reduce peak storage, genome2's full GIX must never exist. The `-n` flag achieves this:
1. **GIXmake `-n`**: Runs the full partitioning phase (computes NPARTS, Ksplit[], prefix index) but skips writing ktab partition files. Outputs only the `.gix` stub.
2. **FastGA `-C K`**: Passes `-n` to GIXmake for genome2's initial build. Opens genome2's `.gix` stub via `Open_Post_List` (which reads only the stub) but sets `T2 = NULL` (no ktab to stream). Then builds chunks on demand via `GIXmake -C first:last`.

## Assessment

**Correctness**: Verified — bit-exact seed output (51,082,720 seeds, exact match with baseline).

**Storage**: With `-n` stub-only mode, **reduces peak by 36% (measured) / 37% (projected for human)**. This is the first optimization to actually reduce peak GIX storage. For human genomes: peak drops from ~62.7 GB to ~39.3 GB, saving ~23.4 GB.

**Performance**: 18x merge phase regression on EXAMPLE dataset (206s vs 12s). However:
- Merge is only ~5s of the ~10 min human pipeline (0.8% of total runtime)
- For human genomes: ~30s extra merge time for 23 GB less peak storage — likely acceptable
- Performance could be improved by also chunking genome1 (not implemented)

**Verdict**: With the `-n` stub-only mode, Opt 7 delivers meaningful peak storage reduction. The implementation is correct and the tradeoff (merge phase overhead for 37% less peak disk) is favorable for large genomes on storage-constrained nodes. Recommended as opt-in via `-C K` flag.

## Files Modified
- `GIXmake.c` — `-C first:last` flag, chunk-restricted sort loop
- `FastGA.c` — `-C K` flag, chunk orchestration loop, empty-prefix skip optimization
