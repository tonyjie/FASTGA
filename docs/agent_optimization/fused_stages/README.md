# fused-B: scan-once + tmpfs ktab — measured results

Empirical validation of the "fused" chunked-merge improvement from
[`../plan_improvement.md`](../plan_improvement.md), via the lower-risk **Approach B**
(design: [`../../superpowers/specs/2026-07-11-fused-scan-once-design.md`](../../superpowers/specs/2026-07-11-fused-scan-once-design.md)).

**What was built.** A GIXmake `-R` (reuse) mode: the initial `-n` scan build persists its pos-lists
(root-named) plus a tiny `Buckets` counts sidecar (`.gcnt`); per-chunk sort builds pass `-R` to
**skip the genome re-scan** and reuse that state. FastGA's chunk loop passes `-R` and deletes the
persisted scratch afterward. Running with the genome + `-P` under **tmpfs (`/dev/shm`)** keeps the
ktab off real disk. This is the **lesser variant** of `plan_improvement.md` item 1 — the ktab is
still written (to RAM), not eliminated in-process (that would be Approach C).

**Three configs compared:** `baseline` (upstream `ddeea32`, no chunking, real disk), `optC`
(`agent-optimization`, `-C`, re-scans per chunk, real disk), `fusedB` (`fused-scan-once`, `-C` +
`-R` scan-once, tmpfs). `real disk` = `du` on the genome/work dir; `tmpfs` = the same on the
`/dev/shm` staging dir; `RAM` = peak RSS (`/usr/bin/time`) + tmpfs occupancy.

## Results

### EXAMPLE (HAP1 × HAP2, ~86 Mbp)

| T | Config | index+merge (s) | wall (s) | real disk (GB) | RAM (GB) | aln | bit-exact |
|--:|---|--:|--:|--:|--:|--:|:--:|
| 32 | baseline | 9 | 17 | 2.29 | 0.72 | 323569 | ref |
| 32 | Opt C `-C8` | 44 | 52 | 0.83 | 0.71 | 323569 | ✅ |
| 32 | **fused-B `-C8`** | **41** | **49** | **0.00** | 1.55 | 323569 | ✅ |
| 8 | baseline | 11 | 29 | 2.26 | 0.64 | 323569 | ref |
| 8 | Opt C `-C8` | 70 | 87 | 0.88 | 0.39 | 323569 | ✅ |
| 8 | **fused-B `-C8`** | **65** | **83** | **0.00** | 1.27 | 323569 | ✅ |

*(RAM for fusedB = RSS + tmpfs; for baseline/optC = RSS only, their scratch is on disk.)*

### Human (GRCh38 × CHM13, ~3.1 Gbp each), T=32

| Config | index+merge (s) | wall (s) | real disk (GB) | RAM (GB) | aln | bit-exact |
|---|--:|--:|--:|--:|--:|:--:|
| baseline (no chunk) | 101 | 607 | **72.6** | 19.1 | 518037 | ref |
| Opt C `-C4` | 193 | 708 | 27.1 | 19.1 | 518037 | ✅ |
| Opt C `-C8` | 234 | 745 | 19.2 | 19.1 | 518037 | ✅ |
| **fused-B `-C4`** | **159** | **668** | **0.00** | 46.5 (19.1+27.4) | 518037 | ✅ |
| **fused-B `-C8`** | **164** | **678** | **0.00** | 40.3 (19.1+21.2) | 518037 | ✅ |

## Findings

1. **Bit-exact, everywhere.** All configs at every T produce the byte-identical alignment payload
   (ONEview, provenance stripped) — EXAMPLE **323,569** aln, human **518,037** aln — matching the
   `ddeea32` baseline at matched thread count. The scan-once reorder is lossless, as designed.

2. **`-R` genuinely skips the scan.** On EXAMPLE, the per-GIXmake "partition" (distribute/scan)
   phase is ~0.2 s for the 2 initial `-n` scans and **0.002 s** for every `-R` chunk build — the
   scan is removed, not merely fast.

3. **Zero real disk.** fused-B's real-disk peak is **0** (vs baseline 72.6 GB, Opt C 19.2 GB on
   human). The whole scratch — persisted pos-lists, per-chunk ktab, seed temps — lives in tmpfs.

4. **Scan-once pays off at scale.** Human `-C8` index+merge drops **234 → 164 s (−70 s)** and wall
   **745 → 678 s (−67 s)** vs Opt C; `-C4` recovers 34 s. More chunks ⇒ more re-scans avoided ⇒
   bigger recovery. On EXAMPLE the win is only ~3 s because scanning 86 Mbp is ~free.

## Honest caveats

- **Disk → RAM relocation, not a net reduction.** fused-B's tmpfs peak (human `-C8`: 21.2 GB) ≈ Opt
  C's real-disk peak (19.2 GB) — the same bytes, moved from disk to RAM. fused-B does **not** shrink
  the footprint below Opt C; the 72.6 → ~19 GB reduction is *chunking* (shared with Opt C). The win
  is *which resource* it uses: **0 disk** in exchange for **~21 GB more RAM** (human). Right trade
  on a disk-constrained, RAM-rich node; wrong trade if RAM is scarce.
- **RSS is unchanged** (~19 GB on human) — the ktab is not held in-process; tmpfs is kernel RAM, not
  process RSS. "Memory" therefore splits into flat RSS vs higher total-RAM-incl-tmpfs.
- **Residual chunking cost.** fused-B wall is still +71 s over baseline on human — that is the
  per-chunk `k_sort`/process overhead, which scan-once does not remove; only the ~70 s of redundant
  scanning is recovered.
- **A true net-reduction** (never materialize the ktab, even in tmpfs) needs Approach C (in-process
  merge from an in-RAM sorted array) — deliberately out of scope here.

## Reproduce

```bash
# EXAMPLE (quick): builds 3 binaries in throwaway worktrees, runs baseline/Opt C/fused-B at each T
bash docs/agent_optimization/fused_stages/run_fused.sh EXAMPLE 8
bash docs/agent_optimization/fused_stages/run_fused.sh EXAMPLE 32
python3 docs/agent_optimization/fused_stages/analyze.py EXAMPLE

# human (~1 h): builds a working 5671357 FAtoGDB (ddeea32's segfaults on CHM13) and shadows each
# config's; fused-B stages the genome + -P under /dev/shm (needs ~40 GB free there)
bash docs/agent_optimization/fused_stages/run_fused.sh human 32
python3 docs/agent_optimization/fused_stages/analyze.py human
```

`stage_data/<dataset>/T<T>/<cfg>/`: `run.Llog` (FastGA/GIXmake `-L`), `run.time`
(`/usr/bin/time -v`), `timeline.tsv` (elapsed, real-disk MB, tmpfs MB, RSS MB), `md5.txt`,
`count.txt`, `leftover_scratch.txt` (0 = the `-R` cleanup worked).
