# Baseline Profiling — upstream FastGA

Profiling of **stock upstream FastGA** (no `optimize-memory` / `agentic-steps` changes) on
the repo EXAMPLE dataset. These are the reference numbers every optimization is compared
against.

Build under test: `main` @ `10ebff7` = upstream `ddeea32` + local `.gitignore`
(freshly `make clean && make`).

## Contents

| Study | Doc | What it measures |
|---|---|---|
| **Performance profiling** | [`performance_profiling.md`](performance_profiling.md) | Overall thread scaling **plus a per-stage breakdown** (GDB / GIX / seed-merge / sort+align) — each stage scales differently — on EXAMPLE, T=1…32 |
| **Storage profiling** | [`storage_profiling.md`](storage_profiling.md) | Peak scratch disk (persistent GIX/GDB + temp) and thread-independence |

Everything needed to reproduce lives in this directory (self-contained):
- `performance_data/` — `run_performance.sh` (driver; `-L` per-stage logs), `results.tsv` + `logs/` (data), `analyze.py` (→ `overall_scaling.png`, `stage_profile.png`).
- `storage_data/` — `run_storage_audit.sh` (driver), `monitor_storage.sh` (polls work+tmp du), `audit/` (data), `plot.py` (→ `storage_profiling.png`).

## Reproduce (both studies)

The baseline is **stock upstream** FastGA — NOT `optimize-memory` / `agent-optimization`
(whose `make` bakes in the optimizations). Build `main` (= upstream `ddeea32` + `.gitignore`)
in a throwaway worktree, then point the drivers at it via `FASTGA=…`:

```bash
git worktree add /tmp/fastga-baseline main && make -C /tmp/fastga-baseline
FASTGA=/tmp/fastga-baseline/FastGA bash performance_data/run_performance.sh
FASTGA=/tmp/fastga-baseline/FastGA bash storage_data/run_storage_audit.sh
python3 performance_data/analyze.py && python3 storage_data/plot.py
git worktree remove /tmp/fastga-baseline
```

See each study's doc for the full step-by-step.
