# Improvement Plan — chunked build/merge (Opt C)

Proposed, not-yet-implemented improvements to the bilateral chunked build/merge (Opt C). See
[`agent_optimization_report.md`](agent_optimization_report.md) for what Opt C does today and
[`human_stages/README.md`](human_stages/README.md) for its measured cost.

---

## 1. Scan the genome once, not once per chunk

### The problem — redundant full-genome scans

Opt C's wall-time penalty is almost entirely **redundant re-computation of GIXmake's Phase A**
(the syncmer scan / "distribute"). On the human pair at `-C16`, `Index + merge` jumps
**102 → 305 s (+204 s, +33 % total wall)** — and that extra time is the scan being repeated.

**Why it repeats (code):**
- `-C cfirst:clast` gates **only Phase B** — the sort/output loop is
  `for (part = CHUNK_FIRST; part <= CHUNK_LAST; part++)` (`GIXmake.c:1453`).
- **Phase A (`distribute` = `sample_thread` + `scan_thread`, `GIXmake.c:669-780`) is *not* gated**
  by `-C`; it runs in full every invocation — it scans the whole genome, finds all syncmers, and
  writes the pos-lists for **all** `NPARTS` partitions.
- The chunk loop drives GIXmake as **one separate process per chunk** (blocking `system()` calls,
  `FastGA.c:5350-5436`), and the pos-list intermediates are **process-local temp files** (opened
  then immediately `unlink`ed, gone at process exit). So each chunk's fresh GIXmake process must
  re-run the full Phase A scan just to regenerate pos-lists it already computed on the previous
  chunk.

**Net:** with `-C16`, each genome is fully scanned **16×** (once per chunk) → **32 full-genome
scans** vs. the baseline's **2**. The scan can't be narrowed to one prefix range (you must compute
every position's k-mer to know its prefix), so the whole-genome scan is paid in full each time.
Phase B (the sort) is *not* redundant — 64 partitions are still sorted once overall, just spread
across the 16 processes.

### The fix — do Phase A once, loop Phase B per chunk

Reorder so the expensive scan happens a single time:

```
scan each genome ONCE  →  keep all NPARTS pos-lists on disk
for chunk k in 0..K-1:
    sort + write only chunk-k partitions of A   (Phase B on the retained pos-lists)
    sort + write only chunk-k partitions of B
    merge chunk-k of A vs chunk-k of B  →  emit seeds
    delete chunk-k ktabs
delete the retained pos-lists
```

This is a pure **reordering of when Phase A runs** — it changes no k-mer, no partition boundary,
no merge result, so it stays **bit-exact**.

Two ways to implement it:
- **(A) Integrated single process.** Teach FastGA's chunked path (or a new GIXmake mode) to run
  `distribute` once, retain the pos-lists, then loop `k_sort` over one chunk at a time and hand
  each chunk to `adaptamer_merge` in-process. Cleanest, but touches the GIXmake/FastGA boundary.
- **(B) Cache Phase A across processes.** Add a "scan-only" GIXmake mode that writes **persistent**
  pos-lists once; each per-chunk GIXmake then reads those instead of re-scanning. Smaller code
  change (keeps the one-process-per-chunk structure), at the cost of a persistent pos-list dir.

### Trade-off — pos-lists must live through the loop

The catch: retaining all pos-lists for both genomes across the chunk loop costs some **persistent
temp storage** that today is freed between chunks. Rough sizing: the pos-lists are the
delta-encoded position stream (≈ the old `.post` intermediate, ~2 GB/genome for human) → **~4 GB**
for both genomes held through the loop. That is small next to the current ~18.7 GB real peak, so
the peak should rise only modestly while the +204 s scan penalty largely disappears.

### Expected outcome (to be measured)

| | today (`-C16`) | after "scan once" |
|---|--:|--:|
| Index + merge wall | 305 s | ~110–130 s (target: near baseline 102 s) |
| Total wall vs baseline | +33 % | ~+2–5 % |
| Real peak scratch | 18.7 GB | ~20–23 GB (pos-lists add ~4 GB) |
| Correctness | bit-exact | bit-exact (reorder only) |

So the improvement trades a **small storage bump** for **most of Opt C's time penalty back** —
turning "−74 % space for +33 % time" into roughly "−70 % space for ~free". Numbers above are
design estimates; validate on the human pair with the `human_stages/` harness.

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
