# Baseline Storage Profiling — current upstream FastGA (EXAMPLE)

Storage footprint of **stock upstream FastGA** (`ddeea32`) on the EXAMPLE dataset
(`HAP1.fasta.gz` × `HAP2.fasta.gz`, ~86 Mbp each), default (no-mask) run.
Companion to [`performance_profiling.md`](performance_profiling.md).

![Baseline storage profiling](storage_profiling.png)

## Peak storage

| Component | Size | Note |
|---|---:|---|
| Persistent GIX + GDB | 1905.3 MB | the `.gix`/`.ktab`/`.post` index + 2-bit `.bps` |
| Temp (seed-pair files) | 438.4 MB | `_pair.*`/`_uniq.*`/`_algn.*`, unlinked-but-open |
| **Peak total** | **2343.7 MB** | both coexist during the seed merge |

**Thread-independent:** peak is identical (2343.7 MB) from T=1 to T=32 — file sizes are set by
genome content + parameters (K=40, `-f`), not by parallelism. More threads only compress the
timeline in wall-clock; they do not change the peak (Panel B).

## Where the storage goes (Panel A, T=8)

The temp-dir timeline shows the pipeline shape:
- two small bumps (~60 MB) — GIX build for HAP1, then HAP2;
- a ramp to the **seed-merge peak (~438 MB)** as seed-pair temp files accumulate while both
  whole GIXs sit on disk;
- a stepped decline through **sort + chain + align** as temp files are consumed.

In stock upstream the GIX is **not** freed after the seed merge — it stays on disk (the ~1.9 GB
persistent block) through the long sort+align tail, even though that phase never re-reads it.
(That plateau is exactly what the `optimize-memory` Opt1 "early GIX deletion" targets.)

## Method

FastGA run with `-k` (keep GIX) and `-P <dedicated tmpdir>`, while a background monitor polls
`du -sb` on the tmpdir every 50 ms — this captures the `open()`-then-`unlink()`ed temp files
that are invisible to `ls`. Persistent files measured with `du` after each run. Threads swept
T = 1…32. All scripts + data live under [`storage_data/`](storage_data/).

## Reproduce

The baseline is **stock upstream** FastGA. This branch's own `make` bakes in the
optimizations, so build `main` (= upstream `ddeea32` + `.gitignore`) separately:

```bash
# 1. Build a stock-upstream baseline binary in a throwaway worktree
git worktree add /tmp/fastga-baseline main
make -C /tmp/fastga-baseline

# 2. Run the storage audit against it (writes storage_data/audit/)
FASTGA=/tmp/fastga-baseline/FastGA \
  bash docs/benchmark_baseline/storage_data/run_storage_audit.sh

# 3. Regenerate this figure from storage_data/audit/
python3 docs/benchmark_baseline/storage_data/plot.py

# cleanup
git worktree remove /tmp/fastga-baseline
```

`storage_data/`: `run_storage_audit.sh` (driver), `monitor_tmpdir.sh` (du poller),
`audit/` (committed `summary.csv` + per-thread monitor TSVs), `plot.py` (→ `storage_profiling.png`).

> Old upstream (`5671357`) was verified byte-identical to `ddeea32` on this profile (the 11
> intervening commits touch only ANO/mask and read paths, not the write path), so these numbers
> are representative of upstream generally. That old-vs-new comparison is archived under
> [`../old_archive/`](../old_archive/).
