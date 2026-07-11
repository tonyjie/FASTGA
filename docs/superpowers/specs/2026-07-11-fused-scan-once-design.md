# Design — Opt C "scan-once + tmpfs ktab" (Approach B)

**Date:** 2026-07-11
**Branch:** `fused-scan-once` (off `agent-optimization`)
**Status:** approved design, pre-implementation

## Goal

Empirically validate the "fused" chunked-merge improvement from
[`../../agent_optimization/plan_improvement.md`](../../agent_optimization/plan_improvement.md),
using the lower-risk **Approach B**: make the bilateral chunked build/merge (Opt C) **scan each
genome once** instead of once per chunk, and run its scratch on **tmpfs** so the ktab never touches
real disk. Measure performance, memory, and storage on the EXAMPLE and human datasets, proving
bit-exact correctness.

Approach B captures the two measurable wins of the full in-RAM fusion (no redundant scans →
wall back near baseline; ~zero real-disk GIX) **without** the high-risk in-process merge refactor
(that was Approach C, explicitly out of scope here).

## Background — why Opt C re-scans today

Opt C (`FastGA.c:5350-5456`) drives `GIXmake` as **one blocking `system()` process per chunk**.
Every `GIXmake` invocation runs `distribute()` (the full-genome syncmer scan) in full, because:
- pos-lists are named with `getpid()` (`GIXmake.c:1912`) and `open()`-then-`unlink()`ed at creation
  (`GIXmake.c:2105-2111`) — so they are process-local and vanish at exit;
- `-C cfirst:clast` gates only `k_sort()` / Phase B (`GIXmake.c:1453`), never `distribute()`.

At `-C16` on human that is **34 full-genome scans** (2 initial `-n` stub builds + 16 chunks × 2
genomes) vs. the baseline's 2 — the source of Opt C's `Index + merge` 102 → 305 s (+204 s / +33 %
wall). `k_sort` reads pos-lists read-only (`GIXmake.c:1459`), and `.split` is written right after
`distribute()` (`GIXmake.c:2133`, before `k_sort`) — so pos-lists are safe to persist and reuse.

## Design

### GIXmake changes (gated on a new cooperative-chunk flag; normal path untouched)

1. **Persistent, deterministic pos-lists.** Enable this whenever **`-n` OR `-R`** is set — both are
   cooperative-chunk signals FastGA already emits only in chunked mode (`-n` on the initial scan
   build, `-R` on the per-chunk sort builds). In that mode, name pos-lists by genome root —
   `<SORT_PATH>/.post.<ROOT>.<k>.idx` — instead of `getpid()`, and **do not** `unlink` them at
   creation. They already survive `k_sort` (read-only), so the `-n` scan build's pos-lists remain on
   disk for the `-R` builds to reuse. The plain path (neither `-n` nor `-R`) keeps the pid-named,
   unlink-at-create behavior, so ordinary/concurrent GIXmake runs never collide or leak.
   (A user running `-n` standalone would leave pos-lists behind — a minor, documented wart, since
   `-n` is effectively FastGA-internal; the chunk loop deletes them after the run.)

2. **`-R` reuse flag (skip scan).** When `-R` is given:
   - **skip `distribute()`** entirely (no sample pass, no scan pass, no pos-list writes);
   - take `NPARTS` / `Ksplit` from `-X` (already wired at `GIXmake.c:2016`);
   - **open** the existing persistent pos-lists by their deterministic name (read-only, no
     `O_CREAT`/`O_TRUNC`);
   - run `k_sort()` for the `-C` partition range only, exactly as today.

   `-R` is always used together with `-C` and `-X`. If the expected pos-lists are missing, error out
   (a caller bug).

### FastGA changes (chunk path only)

- The two initial `-n` stub builds (`FastGA.c:4863`, `4917`) already scan once each. In
  cooperative-chunk mode they now **persist** their pos-lists (via the GIXmake change above). No new
  invocation is added.
- The per-chunk `GIXmake -C … -X …` calls (`FastGA.c:5365-5391`) **add `-R`** so they reuse the
  persisted pos-lists instead of re-scanning.
- After the chunk loop (`FastGA.c:5436`), **delete both genomes' persistent pos-lists**
  (`.post.<ROOT1>.*.idx`, `.post.<ROOT2>.*.idx`), alongside the existing stub deletions.

### tmpfs (run-time, no code)

Run with `-P /dev/shm`. All scratch — persistent pos-lists, per-chunk ktab, and the seed-pair
temps — then lives on tmpfs (RAM), so real-disk usage drops to ~0 (only the persistent `.1gdb`/
`.bps` inputs and `.1aln` output remain on disk). "Storage" is thereby converted to RAM, measured
as tmpfs occupancy + process RSS.

## Correctness

**Invariant:** the single scan's pos-lists are byte-identical to what each per-chunk scan would
regenerate — same genome, same `-T` (→ same `DBsplit`, `NPARTS`), same `Ksplit` (fixed by the
`.split` sidecar). Reuse is therefore a pure reordering of *when* the scan runs, changing no k-mer,
no partition, no merge input.

**Verification (gating, both datasets):** `ONEview <out.1aln> | <strip-provenance> | md5sum`, at
matched thread count, must equal the `ddeea32` baseline —
- EXAMPLE: **323,569** alignments, matched T;
- human: **518,037** alignments, T=32 (the count/md5 already reproduced in the report).

Any mismatch means the pos-list naming/reuse diverged from the scan; the bit-exact check catches it
immediately.

## Measurement plan

Extend the `docs/agent_optimization/human_stages/` harness with a `fused_B` stage. Compare three
configs per dataset:

| Config | build | scratch FS |
|---|---|---|
| baseline | upstream `ddeea32` | real disk |
| Opt C `-C16` | current `agent-optimization` | real disk |
| **fused-B `-C16`** | `fused-scan-once` | **tmpfs (`-P /dev/shm`)** |

Datasets/threads: **EXAMPLE** (T=8 and T=32), **human GRCh38×CHM13** (T=32).

Metrics per run:
- **Performance:** end-to-end wall + the 3-phase breakdown (GDB / Index+merge / Sort+align) from
  `-L`. Key check: does fused-B's `Index+merge` drop from Opt C's ~305 s back toward baseline ~102 s?
- **Storage / memory:** real-disk peak, tmpfs occupancy over time (the "storage that is now RAM"),
  and process peak RSS (`/usr/bin/time`). Total RAM ≈ tmpfs occupancy + RSS.
- **Correctness:** the bit-exact md5 above.

Expected (from the plan, to be confirmed): fused-B wall ≈ baseline (re-scan penalty gone); real
disk ≈ 0; total RAM ≈ ~38 GB on human (vs Opt C's ~18.7 GB real disk); all bit-exact.

## Scope / non-goals

- **In scope:** the `-R` reuse flag + root-named persistent pos-lists in GIXmake; the chunk-loop
  rewire + pos-list cleanup in FastGA; the tmpfs run mode; the three-way measurement.
- **Out of scope:** the in-process merge refactor / feeding `adaptamer_merge` from an in-RAM sorted
  array (Approach C); removing the ~9 GB seed-pair temp (still materializes, on tmpfs here);
  parallelizing the two genomes' builds.

## Risks

- **Pos-list name/lifecycle mismatch** between the scan (`-n`) writer and the `-R` reader → wrong or
  missing input. Contained: deterministic naming is shared code, and the bit-exact check fails loudly
  if reuse diverges.
- **Concurrent runs on the same genome+`-P`** could collide on the now-deterministic pos-list names.
  Mitigated by gating deterministic naming to cooperative-chunk mode only; document that concurrent
  chunked runs must use distinct `-P` dirs.
- **tmpfs sizing:** human needs ~38 GB free in `/dev/shm`; confirm node headroom before the human run.
