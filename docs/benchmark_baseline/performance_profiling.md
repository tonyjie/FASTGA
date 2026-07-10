# Performance Profiling — upstream FastGA (overview)

Per-stage runtime profile of **stock upstream FastGA**, across two datasets. The pipeline runs
as four stages that scale very differently, so the mix shifts with genome size.

| Dataset | Detail | Threads |
|---|---|---|
| **EXAMPLE** (HAP1×HAP2, ~86 Mbp) | [`example/performance_profiling.md`](example/performance_profiling.md) — full T=1…32 sweep, per-stage thread scaling | 1…32 |
| **Human** (GRCh38×CHM13, ~3.1 Gbp) | [`human/performance_profiling.md`](human/performance_profiling.md) — per-stage breakdown + CPU% | 32 |

## Per-stage runtime at T=32 (both datasets)

| Stage | EXAMPLE wall | share | human wall | share | CPU% (human) |
|---|--:|--:|--:|--:|--:|
| GDB (`FAtoGDB` ×2) | 0.9 s | 6% | 19.6 s | 3% | 100% (serial) |
| GIX (`GIXmake` ×2) | 3.4 s | 21% | 87.9 s | 15% | 829% |
| Seed merge | 4.6 s | 28% | 5.5 s | 1% | 2513% |
| Sort + align (+output) | 7.3 s | 45% | 491.3 s | **81%** | 598% |
| **Total** | 16.2 s | 100% | 604 s | 100% | — |
| Peak RSS | 740 MB | | 19.0 GB | | — |

## The cross-dataset story

- **Sort + align dominance grows with scale**: 45% of runtime on EXAMPLE → **81%** on human.
  Alignment work grows faster than genome size, so at human scale it swamps everything else.
  On EXAMPLE the seed merge (28%) and GIX build (21%) are still comparable; on human they are
  1% and 15%. Any wall-clock optimization must target sort+align — which is exactly the phase
  the storage optimizations leave untouched.
- **CPU% exposes where parallelism is wasted** (100% = 1 core; node has 32): the **most
  parallel** stage, seed merge (~25 cores, 2513%), is trivially small (1%), while the
  **dominant** stage, sort+align, averages only ~6 cores (598%) — it is sort/merge/I/O bound,
  leaving most cores idle. This is why overall speedup is sub-linear (see EXAMPLE's T=1…32
  sweep and per-stage scaling).
- **`FAtoGDB` (GDB) is a serial floor** (100% CPU) on both — a fixed cost that grows in
  relative weight only as everything else parallelizes down.
- **32-thread hard cap**: `-T` > 32 is rejected by `GIXmake` (`GIXmake.c:1819`) — a deliberate
  performance cap ("more doesn't help"), see the EXAMPLE doc.

Both datasets reproduce Ashir Rao's report shape (human §4: FAtoGDB ~20 s, GIXmake ~85 s, seed
merge ~6 s, sort+align ~484 s).
