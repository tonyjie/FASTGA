# Storage Profiling — upstream FastGA (overview)

Scratch-storage profile of **stock upstream FastGA**, across two datasets. Peak disk is set by
the two whole GIX indices plus the seed-pair temp files coexisting during the seed merge.

| Dataset | Detail | Threads |
|---|---|---|
| **EXAMPLE** (HAP1×HAP2, ~86 Mbp) | [`example/storage_profiling.md`](example/storage_profiling.md) — timeline by phase + peak vs threads (T=1…32) | 1…32 |
| **Human** (GRCh38×CHM13, ~3.1 Gbp) | [`human/storage_profiling.md`](human/storage_profiling.md) — footprint over time | 32 |

## Peak storage (T=32, both datasets)

| Component | EXAMPLE | Human | scaling |
|---|--:|--:|--:|
| Persistent GIX + GDB | 1905 MB | 64.2 GB | ~34× |
| Temp (seed pairs, alignment) | 438 MB | 8.5 GB | ~20× |
| **Peak total** | **2344 MB** | **72.7 GB** | ~32× |
| Persistent share of peak | 81% | 88% | — |

Storage scales ~linearly with genome size (~36× larger genome → ~32× more peak disk), and the
persistent GIX dominates (81–88%). Peak is **thread-independent** — file sizes are set by
genome content + parameters, not by parallelism (see the EXAMPLE thread sweep).

## The shape over time (human, T=32)

![human storage footprint](human/storage_timeline.png)

The same three-part shape holds at both scales:
1. **GIX build** — the persistent block ramps up (to ~64 GB on human) over the first ~110 s.
2. **Seed merge** (a brief sliver) — temp seed-pairs pile up *on top of* both whole GIXs →
   the peak (72.7 GB human / 2344 MB EXAMPLE).
3. **Sort + align** (80%+ of the run) — temp drains, but the **whole GIX plateau sits on disk
   the entire time even though this phase never re-reads it**.

That idle GIX plateau (≈64 GB doing nothing for 80% of a human run) is the central waste the
optimization work attacks: **Opt1** frees the GIX right after the merge; **Opt3/Opt4** shrink
per-entry size; **Opt C** never materializes both whole GIXs at once, cutting the peak itself
(to ~5 GB at `-C16` on human). See [`../agent_optimization/agent_optimization_report.md`](../agent_optimization/agent_optimization_report.md).

## Why measuring this is tricky

The persistent GIX/GDB are ordinary files (`du`-exact anywhere). The **temp** seed-pair files
are `open()`-then-`unlink()`ed — no directory entry — so they are invisible to `ls`, and `du`
only counts them on a filesystem that keeps a deleted-but-open entry (**NFS** silly-rename;
local ext4/xfs and tmpfs do **not**). The EXAMPLE audit runs on NFS (`du`-correct); the human
audit runs on local `/scratch` and sums the temp from the process's open file descriptors
(`st_size`). See each dataset's doc for details.
