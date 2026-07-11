# Performance Profiling — upstream FastGA (overview)

Per-stage runtime profile of **stock upstream FastGA**, across two datasets. The pipeline runs
as four stages that scale very differently, so the mix shifts with genome size.

| Dataset | Detail | Threads |
|---|---|---|
| **EXAMPLE** (HAP1×HAP2, ~86 Mbp) | [`example/performance_profiling.md`](example/performance_profiling.md) — full T=1…32 sweep, per-stage thread scaling | 1…32 |
| **Human** (GRCh38×CHM13, ~3.1 Gbp) | [`human/performance_profiling.md`](human/performance_profiling.md) — per-stage breakdown + CPU% | 32 |

## Per-stage runtime at T=32 (both datasets)

| Stage | EXAMPLE wall | share | CPU% (EXAMPLE) | human wall | share | CPU% (human) |
|---|--:|--:|--:|--:|--:|--:|
| GDB (`FAtoGDB` ×2) | 0.9 s | 6% | 100% (serial) | 19.6 s | 3% | 100% (serial) |
| GIX (`GIXmake` ×2) | 3.4 s | 21% | 710% (~7) | 87.9 s | 15% | 829% (~8) |
| Seed merge | 4.6 s | 28% | 730% (~7) | 5.5 s | 1% | 2513% (~25) |
| Sort + align (+output) | 7.3 s | 45% | 1337% (~13) | 491.3 s | **81%** | 598% (~6) |
| **Total** | 16.2 s | 100% | — | 604 s | 100% | — |
| Peak RSS | 740 MB | | | 19.0 GB | | — |

*(CPU% = `(user+sys)/wall`, 100% = 1 core, from `-L` per-phase resources, median of 3 reps.
The node has 32 cores.)*

## The cross-dataset story

- **Sort + align dominance grows with scale**: 45% of runtime on EXAMPLE → **81%** on human.
  Alignment work grows faster than genome size, so at human scale it swamps everything else.
  On EXAMPLE the seed merge (28%) and GIX build (21%) are still comparable; on human they are
  1% and 15%. Any wall-clock optimization must target sort+align — which is exactly the phase
  the storage optimizations leave untouched.
- **CPU% exposes where parallelism is wasted** (100% = 1 core; node has 32), and the two
  dominant phases **flip** between datasets:
  - *Seed merge* is the **most parallel** on human (~25 cores, 2513%) but only ~7 cores on
    EXAMPLE — yet it is trivially small (1% of human runtime), so its great scaling barely
    matters.
  - *Sort + align*, the **dominant** stage, is the opposite: ~13 cores on EXAMPLE but only
    ~6 cores (598%) on human — most of the 32 cores sit idle in the phase that owns 81% of the
    run. That flip is the whole reason overall speedup is sub-linear, and it is **not** simply
    "I/O bound" — see the parallelism model below.
- **`FAtoGDB` (GDB) is a serial floor** (100% CPU) on both — a fixed cost that grows in
  relative weight only as everything else parallelizes down.
- **32-thread hard cap**: `-T` > 32 is rejected by `GIXmake` (`GIXmake.c:1819`) — a deliberate
  performance cap ("more doesn't help"), see the EXAMPLE doc.

## How each stage is parallelized (and why the profile looks like this)

Each stage uses a different scheme; the CPU% above is a direct readout of each scheme's ceiling.
Splitting `user` vs `sys` time separates parallel **compute** from parallel **I/O**.

**1. GDB — `FAtoGDB`, single-threaded.** FASTA parse + 2-bit packing is one serial pass. Flat
100% regardless of `-T`; the Amdahl floor (~0.9 s EXAMPLE / ~20 s human).

**2. GIX build — `GIXmake`, two multi-threaded sub-phases, run once per genome (back-to-back).**
The DB is split into `NTHREADS` contig-ranges (`DBsplit`): a *distribute* pass (syncmer scan,
one thread per range) runs very hot but short (~13 cores, 1360% on EXAMPLE), then a *sort+output*
pass (MSD radix sort in `MSDsort.c`, then compress + write `.ktab`) is memory-bandwidth- and
write-bound at ~7 cores (677%) and owns most of GIXmake's wall. Net ≈ 7–8 cores.

**3. Seed merge — `adaptamer_merge`, `NTHREADS` streaming threads (`FastGA.c:2470`).** A single
linear sweep of both sorted GIXs, partitioned by genome-A contig slices. It is **I/O-, not
compute-, bound**: on EXAMPLE its 730% is ~0.9 core `user` + ~6.5 cores `sys` — the "parallelism"
is the kernel streaming/paging the GIX files across threads. That is why it **scales up with data
volume**: human's 62 GB GIX gives the parallel reads far more to chew on → 2513% (~25 cores).
Excellent scaling, but only 1% of human runtime.

**4. Sort + align — dominant, with a hard load-balance ceiling.** Per `2×NPARTS` (strand ×
partition) slice, in a barrier'd sequence (`FastGA.c:4362`…): reimport (parallel) → `rmsd_sort`
(parallel, memory-bound) → `search_seeds` = chain + wave-align (parallel) → `la_sort` (parallel)
→ `la_merge` = **single-threaded** k-way merge + output write. The align work is real compute
(mostly `user`: 95.6 s user / 7.2 s wall on EXAMPLE ≈ 13 cores). **But its granularity is one
genome-A contig = one atomic unit that cannot be split across threads** — `rmsd_sort` places
thread boundaries only at contig boundaries (`RSDsort.c:328-345`); a single contig `part[x]` is
never divided. With the per-thread target `thr = total/32 ≈ 3.1%`, any contig bigger than that
becomes an indivisible oversized chunk: human chr1 is ~8% of the genome (2.6× the target), chr2
~7.8%, … so a few giant chromosomes each pin one thread long past the rest, which finish their
tiny scaffolds and idle → **~6 cores**. EXAMPLE's contigs are comparably sized, so the same code
balances cleanly → ~13 cores. The serial `la_merge` tail (heavier on human's 518 K alignments)
pulls the average down further.

**One line:** seed merge is parallel *I/O* (scales with data **volume** — up on human); sort+align
is parallel *compute* capped by contig-granular load balance (scales with contig **uniformity** —
down on human). Same binary, opposite outcome, driven entirely by data shape.

Both datasets reproduce Ashir Rao's report shape (human §4: FAtoGDB ~20 s, GIXmake ~85 s, seed
merge ~6 s, sort+align ~484 s).
