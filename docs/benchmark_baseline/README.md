# Baseline Profiling — upstream FastGA

Profiling of **stock upstream FastGA** (no `optimize-memory` / `agentic-steps` changes) —
the reference numbers every optimization is compared against. Build under test:
`main` @ `10ebff7` = upstream `ddeea32` + local `.gitignore`.

## Layout

Two cross-dataset **overview** docs at this level, and one subdirectory per dataset with the
detailed studies, data, and scripts:

| | Overview (both datasets) | EXAMPLE detail | Human detail |
|---|---|---|---|
| **Performance** | [`performance_profiling.md`](performance_profiling.md) | [`example/performance_profiling.md`](example/performance_profiling.md) | [`human/performance_profiling.md`](human/performance_profiling.md) |
| **Storage** | [`storage_profiling.md`](storage_profiling.md) | [`example/storage_profiling.md`](example/storage_profiling.md) | [`human/storage_profiling.md`](human/storage_profiling.md) |

- [`example/`](example/) — HAP1×HAP2 (~86 Mbp). Full **T=1…32 sweep**: overall + per-stage
  thread scaling, storage timeline by phase + peak-vs-threads. Self-contained
  (`example/performance_data/`, `example/storage_data/` — scripts + committed data).
- [`human/`](human/) — GRCh38×CHM13 (~3.1 Gbp), **T=32**: per-stage runtime + CPU%, storage
  footprint over time. See [`human/README.md`](human/README.md) for its setup (FAtoGDB
  workaround, temp-measurement method).

## Reproduce

Each dataset's doc has a copy-paste **Reproduce** section. In all cases the baseline binary is
**stock upstream** — NOT `optimize-memory` / `agent-optimization` (whose `make` bakes in the
optimizations) — so build `main` in a throwaway worktree and point the drivers at it via
`FASTGA=…` (the human runs additionally need a working `FAtoGDB`; see `human/README.md`).
