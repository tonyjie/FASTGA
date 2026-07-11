# FastGA — Storage & Memory: what lives where, and why disk at all?

A companion to [`fastga_implementation.md`](fastga_implementation.md) (which covers runtime and
parallelism). Our `agent-optimization` work is all about shrinking FastGA's **disk** footprint — so
the natural question is: *is disk fundamental, or just a workaround for limited memory? If we have
plenty of RAM, could we skip disk entirely?* This doc discusses what FastGA keeps in RAM vs on
disk, **why**, how it moves data between the two, and what "just use RAM" would actually mean.

---

## The headline: the RAM footprint is far smaller than the disk footprint

On the human pair (GRCh38 × CHM13, T=32):

| Data | Size | Lives… |
|---|--:|---|
| GIX index (both genomes) | ~62 GB | **disk** — streamed, never fully in RAM |
| GDB `.bps` (both genomes) | ~1.5 GB | disk — random-access small pieces |
| Seed-pair temp files | ~9 GB | disk — written, then read back |
| **Peak disk scratch** | **~73 GB** | disk |
| **Peak process RSS** | **~19 GB** | RAM |

The one fact to internalize: **the 62 GB GIX is never all in RAM at once.** Peak RSS is only 19 GB
because FastGA *streams* the GIX through small (~1.5 MB) sliding buffers rather than loading it. So
FastGA is already, by construction, a "small memory, big disk" design — it deliberately keeps most
of its data off the heap.

---

## Why disk, not RAM? ("RAM is too small" is only one of the reasons)

It's tempting to assume FastGA writes to disk because the data doesn't fit in RAM. That's part of
it, but not the main story:

1. **It's an external-memory, *streaming* design — and the streaming itself is the speed strategy.**
   FastGA's whole premise ("cache coherence") is that comparing two genomes by **sweeping sorted
   data sequentially** beats random lookups into an in-RAM hash table (which thrash the CPU cache;
   see the guide). The sequential sweep is what makes it fast. Crucially, **that would be true even
   if the bytes lived in RAM** — sequential beats random regardless of the medium. So putting the
   GIX on disk isn't fighting the algorithm; the algorithm is built to stream, and disk is just the
   cheapest place to stream a 62 GB object from.

2. **The GIX is a build-once, reuse-many artifact.** Like a BLAST database, you build a genome's
   index once and can align it against many others. A persistent, on-disk index is the point — you
   don't want to rebuild a 30 GB structure on every run. (`GIXmake` and `FastGA` are even separate
   programs for exactly this reason.)

3. **Cross-process communication.** FastGA orchestrates separate executables — `FAtoGDB`,
   `GIXmake`, `FastGA` — that hand data to each other. Separate processes can only share large state
   through the filesystem.

4. **Robustness and scale.** The disk-backed design runs unchanged on a laptop or a cluster node,
   on a small bacterium or a 3 Gbp mammal, and lets several jobs share a node. It never *assumes*
   the whole index fits in RAM — so it never crashes when it doesn't.

So disk is a mix of a *performance* choice (stream sorted data), an *artifact* choice (reusable
index), an *architecture* choice (multi-process), and a *robustness* choice — with "RAM might be
too small" as just the safety net underneath.

---

## "Could we keep everything in RAM if RAM is sufficient?" — yes, and you already can

Our nodes do have large RAM, so this is a fair question. There are two levels of answer.

### (a) The free way — no code change: put the files on RAM-backed storage

You don't need to rewrite FastGA to make it "RAM-only" — you redirect its files to a RAM-backed
filesystem:

- **Temp/scratch → tmpfs:** point the sort/seed temp dir at `/dev/shm` with `-P /dev/shm`. Every
  "disk" read/write of the seed-pair and sort temps is then actually RAM.
- **GIX/GDB → tmpfs:** build (or copy) the genome directory under `/dev/shm` too, so the `.ktab`
  streaming reads come from RAM.
- **Even without either, the OS page cache already helps:** when there's free RAM, recently
  read/written file pages stay cached, so re-reads are served from memory transparently. The large
  `sys` time we measured in the seed merge (~25 "cores" of kernel I/O — it's paging the GIX in) is
  exactly what shrinks when the GIX already sits in page cache or tmpfs.

**RAM budget, human, RAM-only:** ~73 GB (files on tmpfs) + ~15–19 GB (process working set) ≈
**~90 GB**. Comfortable on a big node — but remember it multiplies with concurrent jobs.

### (b) The rewrite way — hold the GIX in an in-process data structure, no files

- **Upside:** saves the write-then-read round-trip of the GIX and temps.
- **But:** it does *not* speed up the merge itself (already sequential, not I/O-stalled once the
  data is cached); it breaks build-once-reuse and the clean multi-process split; and it fails the
  moment RAM is tight (big genome, modest node, many concurrent jobs). Narrow upside, real downside.

**Bottom line:** if you're RAM-rich, the practical answer is **tmpfs / `-P /dev/shm`**, not a
rewrite. FastGA's streaming design means it doesn't even notice that "disk" is now RAM — you get
RAM-only behavior for free.

---

## How FastGA moves data between memory and disk

The recurring rule: **big, cold, sequential data (the index, the temps) lives on disk and is
streamed; hot working sets (the current sort, the contigs being aligned, the merge heap) live in
RAM.**

| Data | On disk as | How it's moved | In RAM as |
|---|---|---|---|
| Genome sequence | `.bps` (2-bit) | random-access **small pieces** on demand (`Get_Contig`) | only the contigs currently being aligned |
| GIX k-mer table | `.ktab.*` | **streamed** via `read()` in ~1.5 MB blocks (`Kmer_Stream`) | one sliding buffer per thread |
| GIX prefix index | `.gix` stub (128 MB) | loaded once at open | a small index array (+ its inverse) |
| Seed pairs | `_pair.*` temp (open-then-`unlink`) | written via ~1 MB `IOBuffer`s, then read back | buffers only |
| The seed sort | — (stays in RAM) | in-place radix sort | the **big sort array** (`nelmax × swide`, several GB) — the main RAM consumer |
| Alignments | `_algn`/`_uniq` `.las` temp | per-thread gather, read back to sort + dedup | per-thread buffers |
| Final output | `.1aln` | serial k-way merge (heap) → write | ~4 GB merge read-buffers |

Two consequences worth noting:
- **Peak RSS (19 GB) ≈ the sort arrays + stream/merge buffers + current contigs — *not* the index.**
  The index's size barely touches RSS because it is streamed, not held.
- The temps are **open-then-`unlink`ed**: they occupy disk blocks but have no directory entry
  (invisible to `ls`; `du` sees them only on NFS via silly-rename), and are auto-reclaimed if the
  process dies. That's a crash-safety choice, not a memory one.

---

## How this connects to `agent-optimization`

`agent-optimization` reduces peak **disk scratch**. But if you run FastGA RAM-only (on tmpfs), that
scratch *is* RAM — so the very same optimizations directly reduce the **RAM** you'd need to go
disk-free:

| Config | Peak scratch | RAM-only budget (scratch + ~19 GB working set) |
|---|--:|--:|
| baseline | ~73 GB | ~90 GB |
| + Opt3/Opt4 (smaller entries) | ~65 GB | ~80 GB |
| + Opt C `-C16` (chunked) | **~18.7 GB** | **~38 GB** |

So "just add RAM" and "shrink the footprint" are not competing answers — they're complementary.
A smaller footprint is smaller whether it lands on SSD or in `/dev/shm`, and Opt C's chunking is
exactly what turns a ~90 GB RAM-only requirement into ~38 GB — i.e. it makes the RAM-only route
feasible on a far smaller node (or lets many more jobs share a big one). agent-optimization isn't
made obsolete by abundant RAM; it's what makes abundant RAM go further.
