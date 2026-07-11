# FastGA Implementation Notes — Steps, Parallelism, and Cost

The companion to [`fastga_guide.md`](fastga_guide.md). The guide explains *what* FastGA does and
*why*. This document answers, for each stage, the three questions that actually matter for
understanding the performance profile:

1. **How is it parallelized?** — what does one worker (thread) own?
2. **How is the work split into chunks/bins?** — by what criterion, and is it balanced?
3. **What does it cost?** — each step's share of the wall time and its CPU% (≈ how many of the 32
   cores it actually uses), measured at T=32.

**How to read this doc:** each stage opens with an **"In plain terms"** paragraph (intuition, no
code). Read just those four for the big picture. Concrete function/line pointers are parked in a
one-line **Code pointers** note at the end of each stage, so they never get in the way.

> Numbers are measured at T=32 (EXAMPLE HAP1×HAP2 ~86 Mbp / human GRCh38×CHM13 ~3.1 Gbp), median of
> 3 reps, from the FastGA/GIXmake `-L` per-phase logs. Behavior described is stock upstream FastGA.

---

## The whole pipeline at a glance

| Stage | What it produces | Parallelism | Bound by | CPU% (ex / human) |
|---|---|---|---|---|
| 1. GDB | 2-bit packed genome | **1 worker (serial)** | serial parse | 100% / 100% |
| 2. GIX build | sorted k-mer dictionary | genome-slice scan + per-bin sort | mem-bandwidth + writes | 710% / 829% |
| 3. Seed merge | raw match "seeds" | genome-A-slice streaming | **disk I/O** | 730% / 2513% |
| 4. Sort + align | final `.1aln` alignments | per-A-chromosome, serial merge tail | **compute**, capped by load balance | 1337% / 598% |

Two rules of thumb that explain the whole profile:
- **Seed merge is parallel *I/O*** — its "cores" are the kernel streaming the index off disk, so it
  scales with data *volume* (way up on human).
- **Sort + align is parallel *compute*** — capped by how evenly the genome's chromosomes divide, so
  it scales with chromosome *uniformity* (down on human, which has a few giant chromosomes).

> **The two genomes are indexed one after another, not concurrently.** Stages 1–2 process genome A
> and genome B *independently*, so they could in principle run side by side — but when you hand
> FastGA two FASTA files it builds them **serially**, shelling out (blocking `system()` calls)
> `FAtoGDB(A) → GIXmake(A) → FAtoGDB(B) → GIXmake(B)` in sequence, each step using the **full** `-T`
> (not split 16+16). So during `FAtoGDB` only 1 of 32 cores is busy, and during GIX only ~7–8 — and
> that idle time happens **twice, back to back**. Overlapping the two genomes' builds (or hiding the
> serial `FAtoGDB` behind the other genome's GIX) is a real missed opportunity upstream leaves on
> the table; if you build indices yourself you can just run the two `GIXmake` commands concurrently.
> (`FastGA.c:4853-4927`.)

---

## Stage 1 — GDB build (`FAtoGDB`)

> **In plain terms.** The genome arrives as a giant text file of A/C/G/T letters. Stage 1 just
> repacks it into a compact binary (2 bits per base) that every later stage can read quickly and
> index into without re-parsing text — like zipping a book so you can jump to any page instantly.
> One worker, one linear pass.

- **Parallelism:** none. A single thread reads the FASTA and writes the packed `.bps`. Runs once
  per genome, and the two runs are sequential.
- **Cost:** flat **100% CPU** (1 core) no matter what `-T` is — ~0.9 s (EXAMPLE) / ~20 s (human).
  It is the fixed serial floor that caps overall speedup.

*Code pointers:* `FAtoGDB` (FastGA shells out to it via `system()`).

---

## Stage 2 — GIX build (`GIXmake`)

> **In plain terms.** To compare two genomes fast, FastGA first builds each one a **sorted
> dictionary** of its 40-letter "words" (k-mers), so later two genomes can be compared by sweeping
> their dictionaries in lock-step instead of doing billions of random lookups. Building that
> dictionary over ~3 billion positions is too big to sort in one shot, so GIXmake works in two
> moves:
> 1. **Scatter** (Phase A): walk the genome, keep only the "interesting" positions (*syncmers*,
>    ~⅓ of them), and toss each into one of many *bins* by its leading letters — like sorting mail
>    into zip-code bins so no single pile is too big to sort in memory.
> 2. **Sort each bin & write** (Phase B): sort each bin's words alphabetically, write it out as a
>    `.ktab` file, and build a small "thumb-tab" index (the `.gix` stub) so lookups can jump.

The two moves show up as two verbose phases with very different cost (per genome, T=32):

| Phase | share of GIXmake wall | CPU% (ex / human) | ≈ cores |
|---|--:|--:|--:|
| A — scatter / scan | ~8% | 1360% / 1950% | ~14 / ~19 |
| B — sort each bin + write | ~92% | 677% / 780% | ~7 / ~8 |

**The takeaway:** Phase A is very parallel but tiny; **Phase B owns ~92% of the wall and only
reaches ~7–8 cores** — that is what sets GIXmake's overall ~7–8 cores.

### Phase A — scatter / scan

- **How it's parallelized:** split the **genome** across the workers — each worker scans an
  equal-length slice (contig ranges summing to ~equal sequence) completely independently. Disjoint
  regions, no shared state, so it hits ~14–19 of 32 cores. This is the most parallel part of
  GIXmake.
- **How positions go into bins:** each kept position is dropped into one of *N* bins ("parts")
  chosen by the **first few letters of its k-mer**. The bin boundaries are **not** fixed — a quick
  first pass over the genome counts how many positions fall under each prefix, and the boundaries
  are placed so every bin ends up with roughly the **same number of positions**. That keeps the
  later per-bin sorts balanced. (This is why Phase A scans the genome twice: once to count for
  balanced bins, once to actually place.)
- **How many bins (*N*):** chosen so each bin is about **4 GB** when expanded for sorting, rounded
  up to a multiple of the thread count and clamped to [8, 64]. (EXAMPLE → 8; human → up to 64.)
- **Cost:** ~8% of GIXmake's wall, ~14 cores (ex) / ~19 cores (human).

### Phase B — sort each bin & write

- **How it's parallelized — two levels:**
  - The bins are processed **one at a time** (a *serial* outer loop). This is deliberate: all bins
    reuse **one large in-memory sort buffer** (sized to the biggest bin) to cap RAM.
  - *Within* one bin, the sort and the compaction **are** parallel — the work is split by the
    k-mer's **leading byte** into finer sub-buckets, and workers grab sub-buckets from a pool.
  - Writing each finished bin to its `.ktab` file is **serial**.
- **"bin" vs "sub-bucket":** a **bin (part)** is a ~4 GB slice of k-mer-prefix space — the unit the
  outer loop steps through one at a time. A **sub-bucket** is a finer leading-byte slice *inside* a
  bin — the unit the workers parallelize the sort over.
- **Why only ~7–8 cores** despite a parallel sort: three limits stack up — (a) the radix sort is
  **memory-bandwidth-bound** (many cores contend for the same memory bus), (b) the outer per-bin
  loop is **serial**, and (c) each bin's file **write is serial**. Together they cap it well below
  32.
- **Cost:** ~92% of GIXmake's wall, ~7–8 cores. After the last bin, a final **serial** step writes
  the 128 MB `.gix` thumb-tab index + trailer.

*Code pointers:* Phase A = `sample_thread` (count) → `scan_thread`/`push` (place) inside
`distribute()`; Phase B = `setup_thread_*` (rebuild k-mers) → `msd_sort` (`MSDsort.c`, the radix
sort) → `compress_thread` (compact + build prefix index) → ktab write, inside `k_sort()`; the
genome split (`DBsplit`) and bin count (`NPARTS`) are set up in `main()`.

---

## Stage 3 — Seed merge (`adaptamer_merge`)

> **In plain terms.** Now both genomes have sorted dictionaries. Finding matches is a **zipper**:
> walk both dictionaries from A to Z at the same time, and wherever the same word turns up in both,
> record "position X in genome A looks like position Y in genome B." One linear sweep, no lookups.
> These raw matches are called *seeds*; words that are too common (repeats) are skipped.

- **How it's parallelized:** *N* workers, each owning a **slice of genome A** (a contig range),
  sweeping its slice against genome B's index. Same "split the genome into contig ranges" idea as
  Phase A.
- **How work is divided:** by genome-A contig range, ~equal sequence per worker.
- **What it costs — and the twist:** it is **tiny** (5–6 s; ~1% of human runtime) yet shows a very
  high CPU% — **730% on EXAMPLE, 2513% (~25 cores) on human**. But almost all of that is **kernel
  I/O time**, not computation: the "parallelism" is many threads streaming/paging the index off
  disk at once. So it **scales with index size** — human's 62 GB index gives the parallel reads far
  more to do → ~25 cores; EXAMPLE's small index → ~7. Great scaling, but it barely moves the total.
- After this step the index is never read again — which is why deleting it here (our Opt 1) is safe.

*Code pointers:* `new_merge_thread` inside `adaptamer_merge`; streams the `.ktab` in ~1.5 MB blocks.

---

## Stage 4 — Sort + align (`pair_sort_search`)

> **In plain terms.** Stage 3 produced a huge, unordered pile of raw hints ("spot X in A resembles
> spot Y in B"). Stage 4 turns that pile into real alignments in four moves:
> 1. **Sort** the hints so that hints belonging to the same region/diagonal sit next to each other.
> 2. **Chain** nearby hints into runs that trace out one candidate alignment.
> 3. **Align** each promising chain for real — the base-by-base alignment, via a fast wavefront method.
> 4. **Clean up & write** — drop duplicate/overlapping alignments, then merge every worker's
>    results into the final `.1aln`.

This is the **dominant** phase — 45% of EXAMPLE, **81% of human**. Its overall cost:
**1337% (~13 cores) on EXAMPLE, 598% (~6 cores) on human.**

> **No finer breakdown exists yet.** The `-L` log times all four moves as a **single** phase
> ("Sorting and merging alignments" = reimport + sort + chain + align + per-worker sort + final
> merge). `-v` prints per-part progress ("Loading/Sorting/Searching seeds for part N") but with **no
> timings**. The only sub-signal we have is the phase's user/sys split — EXAMPLE is almost pure
> `user` (so chain+align dominates and the serial merge tail is small); human has notable `sys`
> (contig fetch, gather files, the final write). Getting real per-move numbers would require either
> instrumenting the five sub-steps with timers or a sampling profiler (`perf`) — not yet done. So
> the per-move split below is by *character*, not a stopwatch.

- **How the work is divided (moves 1–3):** each worker gets a contiguous set of **genome-A
  chromosomes/contigs**, chosen so each worker's set holds ~equal numbers of seeds — **but the split
  points can only fall *between* chromosomes, never inside one.** The sort, the chaining, and the
  alignment all use this same per-chromosome split.
- **Why ~13 cores on EXAMPLE but only ~6 on human** (the key result): on EXAMPLE the chromosomes are
  similar-sized, so the seed load divides evenly and ~13 cores stay busy. On human a few chromosomes
  are each **bigger than one worker's fair share and can't be split** — chr1 alone is ~8% of the
  genome (2.6× the T=32 target), chr2 ~7.8%, … so a handful of workers grind those giants long after
  the rest finish their tiny scaffolds and go **idle**. That indivisible-giant-chromosome imbalance
  is the ~6-core ceiling — *not* an I/O limit.
- **The serial tail (move 4):** dropping duplicates is done per-worker (parallel), but **merging all
  workers' alignments into the single `.1aln` is done by one thread** — it streams and writes
  everything. It doesn't parallelize and adds a fixed tail, heavier on human's 518 K alignments.
- **Compute vs I/O:** this phase is **compute-dominated** (the wavefront aligner). On EXAMPLE it is
  almost pure compute (~13 cores, nearly all user time); on human there is also substantial I/O
  (fetching contig sequence, the per-worker gather files, and the final merge write), but user time
  still dominates. It **never touches the index** — only the packed genome and the seed/alignment
  temp files.

*Code pointers:* driver `pair_sort_search()`; move 1 = `reimport_thread` + `rmsd_sort`
(`RSDsort.c`, and the per-chromosome split limiter is `RSDsort.c:320-345`); moves 2–3 =
`search_seeds` → `align_contigs` → `Local_Alignment` (the wavefront aligner in `align.c`); move 4 =
`la_sort` (parallel) then `la_merge` (serial).

---

## Parallelism cheat-sheet

For each step: what one worker owns, and whether the step is parallel or serial.

| Stage | Step | One worker owns… | Parallel? |
|---|---|---|---|
| 1 GDB | pack genome | the whole genome | **serial** |
| 2 GIX-A | scatter scan | a slice of the genome | parallel (~14–19 cores) |
| 2 GIX-A | choose bins | — (a reduction) | **serial** |
| 2 GIX-B | sort a bin | a leading-byte sub-bucket | parallel, but **bins done 1-at-a-time** |
| 2 GIX-B | write `.ktab` / `.gix` | — | **serial** |
| 3 merge | zip the two indices | a genome-A contig slice | parallel (I/O-bound) |
| 4 sort+align | sort / chain / align | a set of genome-A chromosomes | parallel, **capped by chromosome sizes** |
| 4 sort+align | final merge → `.1aln` | — | **serial tail** |

The recurring pattern: FastGA parallelizes by **splitting data into independent, lock-free chunks**
(no shared-state contention), then joining at a barrier. Its two ceilings are (1) **chunk
granularity** — a chromosome is atomic, so skewed genomes leave cores idle — and (2) the **serial
spans** — GDB, choosing bins, each file write, and the final merge.
