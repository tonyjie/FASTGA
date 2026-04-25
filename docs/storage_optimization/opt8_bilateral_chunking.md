# Optimization 8: Bilateral Chunking (Chunk Both Genomes)

## Summary

Extend Opt 7 (chunk-wise GIX) to chunk **both** genomes, not just genome2. For each iteration `i` in `1..K`, build only chunk `i` of *both* g1 and g2, merge them together, then delete both before building chunk `i+1`. Since both genomes are prefix-partitioned the same way, chunk `i` of g1 only has k-mer matches with chunk `i` of g2 — not with other chunks. This means we still need only K merge passes (not K×K), and the peak ktab footprint drops from `1/K × g2 + full × g1` to `1/K × (g1 + g2)`.

| Property | Value |
|---|---|
| **Type** | Trade-off (storage ↓, merge time ↑) |
| **Quality tier** | Tier 1 (Bit-exact, seeds/alignments identical) |
| **Target** | Peak ktab: further reduce from ~39 GB (Opt 7) to ~16 GB for human (K=4) |
| **Code changes** | `FastGA.c` (extend chunk loop to g1) + `GIXmake.c` (accept shared Ksplit) |

## Motivation

Current Opt 7 (chunk-wise + `-n`) results for human genomes (projected):

| Component | Size |
|---|---:|
| g1 full ktab (on disk throughout) | ~26.8 GB |
| g2 chunk ktab (1 of 4) | ~6.7 GB |
| g1 stub + g2 stub | ~256 MB |
| **Peak ktab** | **~33.8 GB** |

g1's full ktab dominates. If we chunk g1 too, peak drops to ~13.4 GB ktab (+stubs+temps). That's a further 20+ GB savings on top of Opt 7.

Secondary win: the 18x merge regression in Opt 7 is dominated by chunk 1's full T1 scan (because g1 is whole, chunk 1 iterates all ~2B entries). If g1 is also chunked, chunk 1 only scans its own prefix range — cutting merge time too.

## Design Decision: Partition Alignment

The merge algorithm needs both streams to cover the same prefix range for each chunk pass. The challenge: `GIXmake` computes each genome's `Ksplit[]` (partition boundaries) independently from its own k-mer histogram. For similar-sized genomes the splits are close, but not necessarily identical.

**Three options considered:**

| Option | Description | Complexity | Correctness |
|---|---|---|---|
| A. Force shared `Ksplit[]` | Run stub-pass on g1, pass its Ksplit to g2 via new flag | Medium | Full correctness |
| B. Prefix-range chunking | Chunk by prefix range `[0, P/K)`, each genome builds matching partitions | Higher | Full correctness |
| C. Assume NPARTS match | Hope distributions are close enough | Lowest | **Unsafe** — boundary mismatches lose seeds |

**Chosen: Option A.** Simpler than B, correct unlike C. Implementation:

1. Run `GIXmake -n genome1` first → produces g1's `.gix` stub with `Ksplit[]` and `NPARTS`
2. Run `GIXmake -n -X path/to/g1.gix genome2` → reads g1's Ksplit and reuses it for g2
3. Both genomes now have identical partition boundaries
4. Chunk loop: `for i in 1..K: GIXmake -C first:last g1 && GIXmake -C first:last g2 && merge && delete`

## Implementation

### `GIXmake.c` changes

New flag `-X <reference.gix>`: read `Ksplit[]` and `NPARTS` from an existing `.gix` stub and use them instead of computing fresh. Changes:

1. Parse `-X` flag, open reference `.gix`, read `nparts` and `Ksplit[0..nparts]` from the stub header
2. Skip the histogram-based partition computation (lines ~650–700)
3. Assert that the reference's `KMER` matches this genome's (both must be K=40)

### `FastGA.c` changes

1. When `NCHUNKS > 0`, invoke g1's initial GIXmake with `-n` (stub only), same as g2 currently
2. After g1's stub is built, invoke g2's initial GIXmake with `-n -X <g1.gix>` to inherit Ksplit
3. In the chunk loop, add a `GIXmake -C first:last` call for g1 *before* the g2 call
4. Add g1 chunk deletion after each merge pass (symmetric to g2)
5. Opening streams: T1 is now chunk-limited like T2

### Chunk loop structure (pseudo-code)

```
build g1 stub (GIXmake -n)
build g2 stub (GIXmake -n -X g1.gix)

for chunk in 1..K:
    build g1 chunk (GIXmake -C cfirst:clast g1)
    build g2 chunk (GIXmake -C cfirst:clast g2)
    T1 = open g1 chunk stream
    T2 = open g2 chunk stream
    merge(T1, T2)  // appends seeds to shared temp files
    delete g1 chunk ktab files
    delete g2 chunk ktab files

delete g1 stub
delete g2 stub
proceed to sort+align
```

## Measured Results (EXAMPLE dataset, ~86 Mbp, K=4)

| Configuration | Peak disk | Seeds | Non-redundant aln's | Tier |
|---|---:|---:|---:|---|
| Baseline (no chunk, T=4) | 1,782 MB | 51,082,720 | 323,569 | reference |
| Opt 4+7 (K=4, T=4) | 1,146 MB | 51,082,720 | 323,569 | Tier 1 (bit-exact) |
| **Opt 4+7+8 (K=4, T=1)** | (smaller) | **51,082,720** | **323,569** | **Tier 1 (bit-exact)** |
| **Opt 4+7+8 (K=4, T=4)** | **894.6 MB** | 51,077,369 (−0.01%) | 323,563 (−0.002%) | **Tier 2 (comparable)** |

**Key finding**: Bilateral chunking is **bit-exact at T=1** but **loses ~5,351 seeds at T=4**. Per-chunk seed counts identical for chunks 1–3, all loss is in the last chunk:

| Chunk | T=1 seeds | T=4 seeds | Δ |
|---|---:|---:|---:|
| 1 | 12,806,409 | 12,806,409 | 0 |
| 2 | 12,764,662 | 12,764,662 | 0 |
| 3 | 12,727,317 | 12,727,317 | 0 |
| 4 | 12,784,332 | 12,778,981 | **−5,351** |

This isolates the bug to a thread-boundary interaction inside chunk 4 (the last chunk). The non-chunked T=4 baseline is bit-exact, so threading itself is fine; the regression appears only when chunk-mode prefix-index plateau interacts with the final thread's `pend = 0xffff` bound. **Pending fix.**

**Human genomes (~3.1 Gbp, K=4, projected)**

| Metric | Baseline | Opt 7 | Opt 8 (projected) |
|---|---:|---:|---:|
| Peak ktab | ~62.7 GB | ~39.3 GB | **~16 GB** |
| Peak total disk | ~65 GB | ~41 GB | **~20 GB** |
| Peak reduction vs baseline | — | −37% | **−69%** |

## Risks

1. **Ksplit alignment bug**: if the `-X` flag doesn't correctly propagate Ksplit, partitioning may silently differ between genomes. Mitigation: assert equality by dumping both stubs' Ksplit arrays in verbose mode; diff before starting chunk loop.

2. **Merge correctness when both T1 and T2 are partial**: current Opt 7 merge assumes T1 is complete. Need to verify the merge's prefix-range iteration works when both streams only cover a subset. The empty-prefix skip (already implemented) should handle this transparently since it already deals with T2 having missing prefixes.

3. **GIX re-read cost**: genome1 now gets rebuilt 4 times (for K=4) instead of once. This adds ~4×GIXmake_time to the pipeline. For human genomes GIXmake is ~30s per genome so +90s total — small vs peak storage savings.

4. **NPARTS inherited from g1 may be wrong for g2**: if g2 is much larger, its per-partition size may exceed memory budget. Mitigation: compute NPARTS from the *larger* of the two genomes (pass `-X` from whichever has bigger GDB), not always from g1.

## Files to Modify

- `GIXmake.c`: add `-X path/to/reference.gix` flag, read Ksplit from stub
- `FastGA.c`: extend chunk orchestration loop to symmetric (g1 + g2) rebuilds

## Verification Plan

1. Build both genomes with new `-X` flag, verify Ksplit arrays match (hexdump of stubs)
2. Run full pipeline on EXAMPLE with `-C4`, compare seed count (must be 51,082,720)
3. Compare alignment output (must be 323,569 alignments)
4. Measure peak storage via `du` monitoring — expect ~650 MB peak (vs 1,146 MB Opt 7)
5. Generate 4-panel storage timeline: baseline vs Opt 4 vs Opt 7 vs Opt 8

## Checklist

- [x] `-X` flag implemented in GIXmake.c (with `.split` sidecar containing NPARTS+Ksplit)
- [x] Symmetric chunk loop in FastGA.c (g1 inherits via `-X` from g2's pass, both built with `-n` initially)
- [x] Builds without warnings
- [x] Ksplit alignment verified (sidecars byte-identical via hexdump)
- [x] Bit-exact on EXAMPLE at T=1: 51,082,720 seeds, 323,569 alignments
- [ ] Bit-exact on EXAMPLE at T>1: **−5,351 seeds in last chunk only — multi-thread bug, pending fix**
- [x] Peak storage measured: 894.6 MB (−50% vs 1,782 MB baseline) at T=4 -C4
- [ ] Storage timeline figure generated
- [x] Results documented
- [x] README.md updated

## Bugs Found / Fixed During Implementation

### Bug 1: `Free_Kmer_Stream(NULL)` in chunk-loop preamble

In bilateral-chunk mode both T1 and T2 are deferred (NULL before the loop). Opt 7's chunk preamble called `Free_Kmer_Stream(T1)` unconditionally, which segfaults on NULL.

**Fix**: Guard with `if (T1 != NULL)`. (FastGA.c, around line 5247.)

### Bug 2: Per-chunk `bzero` of `buck[]` discarded cumulative seed counts (latent in Opt 7)

`(self_)adaptamer_merge` was zeroing the per-thread per-contig `buck[]` arrays at the start of each call. With chunked merges, only the **last chunk's** seed counts survived — but the temp seed files contained seeds from **all chunks** cumulatively. The reimport phase then under-sized `sarray` and wrote past its end, causing SIGSEGV at "Loading seeds for part 1".

This is a latent bug in committed Opt 7: its verification only checked seed counts (which are written and counted independently of `buck[]`), not alignment counts. Opt 8 was the first thing exercising the full chained chunk path end-to-end and surfaced it.

**Fix**: Hoist the `bzero` out of `adaptamer_merge` and `self_adaptamer_merge` into a single pre-merge zero pass at the caller. Lifecycle of `buck[]` belongs to the caller; chunk loop accumulates naturally without per-chunk resets.
