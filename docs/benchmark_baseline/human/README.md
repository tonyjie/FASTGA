# Human-genome Profiling — upstream FastGA (T=32)

Profiling of **stock upstream FastGA** on the full human pair **GRCh38 × CHM13**
(~3.1 Gbp each), at T=32. Companion to the EXAMPLE studies one level up
([`../performance_profiling.md`](../performance_profiling.md),
[`../storage_profiling.md`](../storage_profiling.md)).

| Study | Doc | Result |
|---|---|---|
| **Performance** | [`performance.md`](performance.md) | per-stage runtime breakdown + CPU% (end-to-end ~604 s, peak RSS 19 GB) |
| **Storage** | [`storage.md`](storage.md) | footprint over time (peak **72.7 GB**), by phase |

Both reproduce Ashir Rao's report numbers (~594 s / ~71 GB).

## Setup (shared by both)

Binary: **ddeea32 FastGA/GIXmake + a 5671357 FAtoGDB**. `ddeea32`'s own FAtoGDB
**segfaults** building CHM13's GDB (the masked-sequence/`.ano` regression — see
`../../agent_optimization_report.md`), so the GDB stage uses the working
old-upstream FAtoGDB; FastGA/GIXmake are current upstream.

Measurement filesystems:
- **Performance** runs on local `/scratch` (fast, representative timing).
- **Storage** also runs on local `/scratch`, but there `du` can't see the
  `open()`-then-`unlink()`ed temp files (no directory entry; NFS silly-rename is
  what makes `du` work elsewhere). Temp is instead summed from the process's open
  **file descriptors** (`/proc/PID/fd`, inode-deduped, using `st_size` = bytes
  written — `st_blocks` would over-count ~2× because XFS speculatively
  preallocates blocks for the actively-written seed files). Persistent GDB/GIX are
  real files, so `du(work)` is exact.

## Reproduce

```bash
# baseline binary: ddeea32 + working (5671357) FAtoGDB
git worktree add /tmp/fastga-baseline main && make -C /tmp/fastga-baseline
git worktree add --detach /tmp/fga57 5671357 && make -C /tmp/fga57 FAtoGDB
cp /tmp/fga57/FAtoGDB /tmp/fastga-baseline/FAtoGDB

# 3 perf reps + 1 storage run (edit G1/G2 to your FASTA paths)
BL=/tmp/fastga-baseline bash docs/benchmark_baseline/human/run_human.sh
python3 docs/benchmark_baseline/human/analyze_perf.py    # -> perf_breakdown.png + table
python3 docs/benchmark_baseline/human/plot_storage.py     # -> storage_timeline.png
```

`perf_data/` (rep*.Llog + rep*.time), `storage_data/` (timeline.tsv + storage.Llog),
`run_human.sh`, `analyze_perf.py`, `plot_storage.py`.
