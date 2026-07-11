# Improvement Plan — chunked build/merge (Opt C)

Proposed, not-yet-implemented improvements to the bilateral chunked build/merge (Opt C). See
[`agent_optimization_report.md`](agent_optimization_report.md) for what Opt C does today and
[`human_stages/README.md`](human_stages/README.md) for its measured cost.

---

## 1. Fuse per-bin sort with per-bin merge — never write the ktab to disk

### The insight

GIX build and seed merge decompose into **independent per-bin work**: bin *i*'s k-mers of A can
only ever match bin *i*'s k-mers of B (a match needs ≥12 shared prefix bases ⊃ the same 5-base
bucket ⊃ the same bin; both genomes share the same `Ksplit` boundaries). So you **never need the
whole sorted index at once** — you can *sort one bin, merge it, throw it away, and move to the
next*. The sorted ktab is a throwaway intermediate, not something that has to be persisted.

Today FastGA persists it anyway: GIXmake sorts every bin and writes **all** of it to `.ktab`
(~62 GB for human), then the seed merge streams that same 62 GB back off disk. The preferred fix is
to **fuse the two so the ktab is never materialized**:

```
scan each genome ONCE  →  pos-lists for all bins (held in RAM ~4 GB, or small on disk)
for bin i in 0..NPARTS-1:
    sort bin i of A in RAM            (from A's pos-lists)
    sort bin i of B in RAM            (from B's pos-lists — same Ksplit boundaries)
    merge A-bin_i vs B-bin_i in RAM   →  emit seeds
    discard both sorted bins
```

Peak extra RAM = one A-bin + one B-bin held sorted at a time (~4 GB each, the existing ~4 GB/bin
sizing) plus the pos-lists. The **62 GB ktab is never written**.

### Why it is strictly better than both current paths

- **vs. baseline:** skips writing the ~62 GB ktab *and* reading it back for the merge (~124 GB of
  I/O avoided). Not slower — likely a touch faster.
- **vs. Opt C:** Opt C re-scans the whole genome **per chunk** (a separate GIXmake process per
  chunk; pos-lists are process-local temp), costing **+204 s / +33 % wall** on human. Fusing scans
  each genome **once** — no redundant computation at all.
- **Correctness:** the per-bin merge is complete because matches never cross bins — *provided both
  genomes use identical bin boundaries* (exactly what the shared `Ksplit` / `-X` guarantees). It is
  a reorder + fuse only, so **bit-exact**.

For reference, why Opt C's re-scan happens today (the thing this proposal removes): `-C` gates only
Phase B (`for (part = CHUNK_FIRST; ...; part++)`, `GIXmake.c:1453`); Phase A (`distribute` =
`sample_thread` + `scan_thread`, `GIXmake.c:669-780`) runs in full every invocation; and the chunk
loop calls GIXmake as one `system()` process per chunk (`FastGA.c:5350-5436`) whose pos-lists are
`unlink`ed at exit — so every chunk re-scans the whole genome. With `-C16` that is **32 full-genome
scans vs. 2**. Fusing does the scan once and keeps its output in-process.

### What it costs / caveats

- **The pos-lists are still needed.** The scan must finish before *any* bin is complete (a bin's
  members are scattered across the whole genome), so both genomes' pos-lists (~2 GB each, the
  delta-encoded position stream) must exist after the scan. Hold them in RAM (~4 GB → truly zero
  disk) or on disk (~4 GB, negligible next to 62 GB).
- **The ~9 GB seed-pair temp still materializes.** This removes the *GIX* from disk, not the seeds;
  the seeds are consumed later by sort+align. (Fusing that too is the harder Backlog item.)
- **You give up the reusable persistent GIX.** Ideal for a **one-shot pairwise A-vs-B** (or
  disk-constrained / RAM-rich runs); worse for **one-vs-many**, where a persistent index is built
  once and reused — here you would re-scan A on every comparison.
- **Higher peak RAM** — bins are held in RAM instead of streamed from disk. This is the intended
  trade when RAM is the abundant resource (see
  [`../basics/fastga_storage_and_memory.md`](../basics/fastga_storage_and_memory.md)).
- **Engineering.** Fuse GIXmake's `k_sort` (bin sort) with FastGA's `adaptamer_merge` (bin merge)
  into one in-process pipeline; the merge must read from an in-RAM sorted array instead of a
  `Kmer_Stream` (disk). Moderate refactor, clean logically — and it also drops the multi-process
  `system()` split that the chunked path uses today.

### Lesser variant (if not fusing)

If GIXmake and the merge stay separate steps, at minimum **scan once and retain the pos-lists**,
then loop per-chunk sort → merge → delete reading those retained pos-lists. This alone removes Opt
C's re-scan penalty (recovers most of the +204 s) while still briefly writing each chunk's ktab to
disk and reading it back — simpler, but keeps a small ktab-on-disk footprint.

### Expected outcome (to be measured)

| | baseline | Opt C `-C16` today | **fused (preferred)** |
|---|--:|--:|--:|
| GIX on disk | ~62 GB | ~5 GB sawtooth | **~0 (RAM)** |
| Redundant full-genome scans | none | 32× | **none** |
| ktab write + read-back | yes | yes (per chunk) | **skipped** |
| Total wall vs baseline | — | +33 % | **~0 to slightly faster** |
| Peak scratch (disk) | ~73 GB | ~18.7 GB | **~9 GB (just the seeds)** |
| Peak RAM | ~19 GB | ~19 GB | ~19 + ~8 GB (two bins) |
| Correctness | ref | bit-exact | bit-exact |

The ~9 GB seed-pair temp is the remaining disk floor in every mode; fusing removes essentially
everything else. Numbers are design estimates — validate on the human pair with the `human_stages/`
harness.

---

## Backlog (not yet analyzed)

- **Chunk the seed-pair temp too?** Chunking removes the *GIX* from the peak but not the ~9 GB
  seed-pair temp (it accumulates across all chunks until sort+align consumes it). Consuming seeds
  incrementally — align each chunk's seeds before building the next chunk — could cut that too, but
  it entangles the (currently clean) build/merge and sort/align phases and needs its own design +
  correctness study.
- **Parallelize the two genomes' index builds** (orthogonal to Opt C): today `FAtoGDB`/`GIXmake`
  run strictly serially per genome, each using full `-T`, leaving cores idle during the serial
  `FAtoGDB` and the ~7–8-core GIX build (see
  [`../basics/fastga_implementation.md`](../basics/fastga_implementation.md)). Overlapping A and B
  could recover wall time independent of storage.
