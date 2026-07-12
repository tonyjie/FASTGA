# Branch registry

What each git branch in this repo is for — so work doesn't get lost over time. This is the
canonical copy; it lives on **`agent-optimization`** (the hub that feature branches derive from).
From another branch, read it with `git show agent-optimization:docs/BRANCHES.md`. Update it here
when a branch's purpose or status changes.

_Last updated: 2026-07-11._

## Active

| Branch | Purpose | Status |
|---|---|---|
| **main** | Upstream FastGA mirror + one `.gitignore` commit. **Never develop here** — sync via rebase + force-push. | Mirror of upstream `ddeea32`. |
| **agent-optimization** | The primary working base. Consolidated, verified storage optimizations on top of upstream: **Opt1** (early GIX deletion), **Opt3** (drop mask byte), **Opt4** (drop LCP byte, on-the-fly), **Opt C** (bilateral chunked build/merge, `-C`). All bit-exact; benchmarked on EXAMPLE + human. Holds the docs tree (`docs/benchmark_baseline`, `docs/agent_optimization`, `docs/basics`). | Active base. |
| **fused-scan-once** | Validation of the "fused" Opt-C improvement, **Approach B**: GIXmake `-R` scan-once (persist pos-lists + `.gcnt` counts sidecar via `-J<pid>`, skip the per-chunk re-scan) + tmpfs ktab. **Bit-exact**; measured **zero real disk**, **~70 s** scan-once win on human `-C8`, at the cost of **+21 GB RAM** (a disk→RAM relocation, not a net footprint reduction — that needs Approach C). Off `agent-optimization`. Docs: `docs/agent_optimization/fused_stages/`; design/plan under `docs/superpowers/`. | Pushed; not merged. |

## Superseded / historical

| Branch | Purpose | Status |
|---|---|---|
| **optimize-memory** | Earliest storage-opt effort: Opt1 + Opt3 + **Opt5 (ktab compression — FAILED)**. Carries the `old_archive` docs. | Superseded by `agent-optimization`. |
| **agentic-steps** | Ashir Rao's branch: Opt4 + Opt C **as originally shipped** — had a `KBYTE` off-by-one that crashed / silently dropped alignments; root-caused and fixed when folding into `agent-optimization`. | Superseded (fix lives in `agent-optimization`). |
| **agent-optimize-1218_1422**, **agent-optimize-1218_1615**, **0v** | Older dated snapshots (Dec 2025 / Mar 2026) — early optimization summaries / guide work. | Purpose superseded; verify before reusing. |
| **agent-optimization-wt** | Transient git-worktree branch, not a tracked line of work. | Ephemeral. |

## Conventions
- Feature/experiment work branches **off `agent-optimization`**, not `main`.
- Each substantial effort gets a spec + plan under `docs/superpowers/` and a results dir under `docs/`.
- Update this file (and the mirror note in Claude's memory, `project_branch_registry`) whenever a branch's status changes.
