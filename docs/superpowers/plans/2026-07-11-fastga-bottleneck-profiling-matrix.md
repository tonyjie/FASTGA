# FastGA Bottleneck Profiling Matrix (v1 divergence axis) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the reusable tooling + dataset inventory to profile FastGA's per-phase bottleneck across a divergence gradient (CHM13 × {GRCh38, chimp, siamang, pig, mouse}, T=32) and show the `sort+align` share falling from ~81% as divergence grows.

**Architecture:** Generalize the existing, proven human harness (`docs/benchmark_baseline/human/run_human.sh` + `analyze_perf.py`) into per-pair run scripts + a cross-dataset aggregator. Deliverable mode is **plan + scripts only**: the executable tasks build and unit-test the tooling (the parser is tested against committed `-L` logs; shell scripts are tested via `--dry-run` and the disk guard). The actual multi-GB genome downloads and ~80 GB/run executions are triggered by the user afterward.

**Tech Stack:** Bash (harness, downloader), Python 3 + matplotlib/numpy (aggregation + figure), pytest (parser unit test), FastGA `-L` per-phase logs + `/usr/bin/time -v`.

## Global Constraints

- **All genome I/O on `/scratch/jl4257`** (local `VG00-scratch`, ~877 GB free, no per-user quota). Never write genome data to `/home/jl4257` or the shared-users NFS.
- **Disk guard mandatory:** no run proceeds if free space on the scratch target `< 150 GB` (configurable via `MIN_FREE_GB`). Exit non-zero with a clear message.
- **Threads: T=32** for every run (matches the existing human datapoint).
- **Baseline binary dir** via `BL=<dir>`; known location `/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA-baseline` (upstream `ddeea32` FastGA/GIXmake + `5671357` FAtoGDB). Scripts must not hardcode a different build.
- **Reference genome:** CHM13 v2.0 at `/scratch/jl4257/seq_align/fastga_datasets/CHM13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna` (already local). GRCh38 already local at `…/GRCh38/GCF_000001405.40_GRCh38.p14_genomic.fna`.
- **Primary divergence x-axis = ordinal phylogenetic rank** (human 0 < chimp 1 < siamang 2 < pig 3 < mouse 4). Secondary continuous proxy = non-redundant aligned bp (`N_aln × ave_len`) ÷ query genome size, recorded but not required for the trend.
- Deliverables live under `docs/benchmark_matrix/`. Commit after each task.

---

## File Structure

```
docs/benchmark_matrix/
  datasets_inventory.md          # Task 1  — full paper-dataset reference
  aggregate_matrix.py            # Task 2+3 — parse -L logs, cross-dataset table + figure
  tests/
    fixtures/human/logs/         # Task 2  — copies of committed rep1-3.Llog for the unit test
    test_parse_llog.py           # Task 2  — pytest for the -L parser
  run_pair.sh                    # Task 4  — one genome pair: perf reps + disk guard
  fetch_genomes.sh               # Task 5  — download the 4 query genomes, disk precheck
  run_divergence_axis.sh         # Task 6  — orchestrate the 5 points
  README.md                      # Task 6  — how to run + results write-up
  divergence/                    # produced by runs (gitignored except results.tsv + png)
    <label>/{logs,time}/
    results.tsv
    divergence_phase_share.png
```

---

### Task 1: Dataset inventory document

**Files:**
- Create: `docs/benchmark_matrix/datasets_inventory.md`

**Interfaces:**
- Consumes: nothing.
- Produces: a reference doc; no code depends on it.

- [ ] **Step 1: Write the inventory doc**

Create `docs/benchmark_matrix/datasets_inventory.md` with exactly this content:

````markdown
# FastGA paper — real dataset inventory

Datasets used in Myers, Durbin & Zhou, *FastGA: Fast Genome Alignment*, bioRxiv 2025
(doi:10.1101/2025.06.15.659750), §5. Recorded here as the reference space for the
bottleneck profiling matrix. **v1 profiles the mammalian divergence axis only**; the DToL
and simulated rows are logged for future size/species/sensitivity axes.

## §5.2 Mammalian → CHM13 (divergence axis, ~3 Gbp fixed)

| Role | Species | Size (Mb) | Accession / source | Paper FastGA CPU (min) |
|---|---|--:|---|--:|
| reference | CHM13 v2.0 (T2T) | 3,117 | GCF_009914755.1 | — |
| query | human GRCh38.p14 | 3,298 | GCF_000001405.40 | 70.5 |
| query | chimpanzee | 3,178 | GCF_028858775.2 | 28.7 |
| query | siamang | 3,263 | GCF_028878055.3 | 26.8 |
| query | pig | 2,612 | http://gigadb.org/dataset/102692 | 17.4 |
| query | mouse | 2,731 | GCA_964188535.1 | 16.6 |

All mammalian runs < 20 GB peak RSS (paper). Runtime falls with evolutionary distance
because closely related genomes have more alignable sequence → more alignment work.

## §5.3 DToL — six genera, within- + between-species (future size axis)

Per-species accessions are in the paper's **Supplementary Table S3** (retrieve when the size
axis is scheduled). Sizes and paper CPU-time from Table 1:

| Class | Genus | A (Mb) | B (Mb) | within CPU (min) | between CPU (min) |
|---|---|--:|--:|--:|--:|
| Insect | Acronicta (moth) | 405 | 466 | 10.7 | 9.2 |
| Fish | Thunnus (tuna) | 792 | 782 | 22.9 | 21.2 |
| Bird | Ammospiza (sparrow) | 1,241 | 1,398 | 45.4 | 33.2 |
| Reptile | Vipera (snake) | 1,632 | 1,695 | 153.8 | 70.9 |
| Mammal | Molossus (bat) | 2,505 | 2,567 | 43.5 | 51.2 |
| Amphibian | Lissotriton (newt) | 24,226 | 23,170 | 4,611 | 2,539 |

Newt (~24 Gbp) is the storage-stress point: naive GIX peak ≈ 600 GB; use `agent-optimization`
`-C16` chunked build (persistent peak ~5 GB) when scheduled.

## §5.1 Simulated genomes (future sensitivity axis)

Pair of ~84 Mb genomes: 10 kb blocks, each a similarity region (length 100 bp–5 kb) at
divergence 1%–65% (SNV 80% / ins 10% / del 10% on B) then random sequence, blocks shuffled
(no long-range alignments). 100 replicates per (length, divergence) cell. Construction
parameters only — not a download.
````

- [ ] **Step 2: Validate required rows are present**

Run:
```bash
cd docs/benchmark_matrix
for acc in GCF_009914755.1 GCF_000001405.40 GCF_028858775.2 GCF_028878055.3 GCA_964188535.1 102692; do
  grep -q "$acc" datasets_inventory.md || { echo "MISSING $acc"; exit 1; }
done
echo "inventory OK"
```
Expected: `inventory OK`

- [ ] **Step 3: Commit**

```bash
git add docs/benchmark_matrix/datasets_inventory.md
git commit -m "docs(matrix): paper dataset inventory for profiling matrix"
```

---

### Task 2: `-L` log parser (TDD against committed human logs)

**Files:**
- Create: `docs/benchmark_matrix/aggregate_matrix.py`
- Create: `docs/benchmark_matrix/tests/test_parse_llog.py`
- Create (copy fixtures): `docs/benchmark_matrix/tests/fixtures/human/logs/rep1.Llog` (+rep2,rep3)

**Interfaces:**
- Consumes: FastGA `-L` logs in the format of `docs/benchmark_baseline/human/perf_data/logs/rep*.Llog`.
- Produces:
  - `sec(tok: str) -> float` — parse `"6.706"`, `"2:06.205"`, `"8:14.618"` (s, m:s, h:m:s) to seconds.
  - `res(line: str) -> tuple[float,float,float] | None` — parse `"…Xu Ys Zw P%"` → `(user, sys, wall)`.
  - `parse_llog(path: str) -> dict` — returns `{"GDB": wall_s, "GIX": wall_s, "Seed merge": wall_s, "Sort+align": wall_s, "cpu": {stage: percent}, "rss_mb": float, "n_aln": int, "ave_len": int}`. Walls are summed over the two genomes for GDB/GIX; FastGA phases come from the **last** FastGA invocation's phase lines.

- [ ] **Step 1: Copy the committed logs as fixtures**

Run:
```bash
mkdir -p docs/benchmark_matrix/tests/fixtures/human/logs
cp docs/benchmark_baseline/human/perf_data/logs/rep*.Llog docs/benchmark_matrix/tests/fixtures/human/logs/
ls docs/benchmark_matrix/tests/fixtures/human/logs/
```
Expected: `rep1.Llog  rep2.Llog  rep3.Llog`

- [ ] **Step 2: Write the failing test**

Create `docs/benchmark_matrix/tests/test_parse_llog.py`:
```python
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from aggregate_matrix import sec, res, parse_llog

FIX = os.path.join(HERE, "fixtures", "human", "logs")

def test_sec_formats():
    assert sec("6.706") == 6.706
    assert abs(sec("2:06.205") - 126.205) < 1e-6
    assert abs(sec("8:14.618") - 494.618) < 1e-6

def test_res_line():
    assert res("  Resources for phase:  40:05.580u  9:19.434s  8:14.618w  599.5%") == \
        (2405.58, 559.434, 494.618)
    assert res("no resources here") is None

def test_parse_human_shares():
    p = parse_llog(os.path.join(FIX, "rep1.Llog"))
    total = p["GDB"] + p["GIX"] + p["Seed merge"] + p["Sort+align"]
    # Known human breakdown: GDB ~3%, GIX ~15%, merge ~1%, sort+align ~81%
    assert 0.78 <= p["Sort+align"] / total <= 0.84
    assert p["Seed merge"] / total < 0.03
    assert 0.10 <= p["GIX"] / total <= 0.20
    assert p["n_aln"] == 518037
    assert p["rss_mb"] == 19.0
```

- [ ] **Step 3: Run test to verify it fails**

Run: `cd docs/benchmark_matrix && python -m pytest tests/test_parse_llog.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'aggregate_matrix'` (or ImportError).

- [ ] **Step 4: Write minimal implementation**

Create `docs/benchmark_matrix/aggregate_matrix.py` (parser portion only for now):
```python
#!/usr/bin/env python3
"""Aggregate FastGA -L per-phase logs across a dataset matrix.
Parser reproduces the human baseline breakdown (docs/benchmark_baseline/human)."""
import re

STAGES = ["GDB", "GIX", "Seed merge", "Sort+align"]

def sec(tok):
    """'6.706' | '2:06.205' | '8:14.618' -> seconds."""
    v = tok.split(":")
    if len(v) == 1:
        return float(v[0])
    if len(v) == 2:
        return float(v[0]) * 60 + float(v[1])
    return float(v[0]) * 3600 + float(v[1]) * 60 + float(v[2])

def res(line):
    """'... Xu Ys Zw P%' -> (user, sys, wall) or None."""
    m = re.search(r'([\d:.]+)u\s+([\d:.]+)s\s+([\d:.]+)w', line)
    return (sec(m.group(1)), sec(m.group(2)), sec(m.group(3))) if m else None

def _cpu(line):
    m = re.search(r'([\d.]+)%', line)
    return float(m.group(1)) if m else None

def parse_llog(path):
    """Per-phase wall seconds. GDB/GIX summed over the two genomes; FastGA phases from the
    last FastGA invocation. Also returns cpu%, peak RSS (MB), n_aln, ave_len."""
    out = {s: 0.0 for s in STAGES}
    cpu = {s: [] for s in STAGES}
    rss_mb = 0.0; n_aln = 0; ave_len = 0
    sect = None            # 'gdb' | 'gix' | 'fastga'
    fastga_phase = 0       # 0=seed merge next, 1=sort+align next
    lines = open(path).read().splitlines()
    for ln in lines:
        if re.search(r'FAtoGDB\b.* -', ln): sect = "gdb"
        elif re.search(r'GIXmake\b.* -', ln): sect = "gix"
        elif re.search(r'FastGA\b.* -T', ln): sect = "fastga"; fastga_phase = 0
        elif "Total Resources" in ln:
            r = res(ln); c = _cpu(ln)
            if sect == "gdb" and r: out["GDB"] += r[2]; cpu["GDB"].append(c)
            elif sect == "gix" and r: out["GIX"] += r[2]; cpu["GIX"].append(c)
            elif sect == "fastga":
                m = re.search(r'([\d.]+)MB', ln)
                if m: rss_mb = float(m.group(1))
        elif "Resources for phase" in ln and sect == "fastga":
            r = res(ln); c = _cpu(ln)
            if r is None: continue
            if fastga_phase == 0:
                out["Seed merge"] = r[2]; cpu["Seed merge"] = [c]; fastga_phase = 1
            elif fastga_phase == 1:
                out["Sort+align"] = r[2]; cpu["Sort+align"] = [c]; fastga_phase = 2
            # phase 2 = PAF conversion, ignored
        elif "non-redundant aln" in ln:
            m = re.search(r'([\d]+)\s+non-redundant aln.s of ave len\s+([\d]+)', ln)
            if m: n_aln = int(m.group(1)); ave_len = int(m.group(2))
    out["cpu"] = {s: (sum(v)/len(v) if v and all(x is not None for x in v) else None)
                  for s, v in cpu.items()}
    out["rss_mb"] = rss_mb; out["n_aln"] = n_aln; out["ave_len"] = ave_len
    return out
```

- [ ] **Step 5: Run test to verify it passes**

Run: `cd docs/benchmark_matrix && python -m pytest tests/test_parse_llog.py -q`
Expected: PASS (3 passed).

- [ ] **Step 6: Commit**

```bash
git add docs/benchmark_matrix/aggregate_matrix.py docs/benchmark_matrix/tests/
git commit -m "feat(matrix): -L log parser with human-baseline unit test"
```

---

### Task 3: Cross-dataset aggregation + stacked-share figure

**Files:**
- Modify: `docs/benchmark_matrix/aggregate_matrix.py` (add point aggregation, table, figure, CLI)
- Modify: `docs/benchmark_matrix/tests/test_parse_llog.py` (add an aggregation test)

**Interfaces:**
- Consumes: `parse_llog` (Task 2); a directory `divergence/<label>/logs/rep*.Llog`.
- Produces:
  - `point_median(logs_dir: str) -> dict` — median across reps of each stage wall + cpu + rss + n_aln.
  - `aggregate(base_dir: str, points: list[tuple[str,int]]) -> str` — writes `results.tsv` and `divergence_phase_share.png` under `base_dir`; returns the tsv path. `points` = `[(label, rank), …]`.
  - CLI: `python aggregate_matrix.py <divergence_dir> label:rank [label:rank ...]`.

- [ ] **Step 1: Write the failing test**

Append to `docs/benchmark_matrix/tests/test_parse_llog.py`:
```python
def test_aggregate_one_point(tmp_path):
    from aggregate_matrix import aggregate
    d = tmp_path / "divergence" / "human" / "logs"
    d.mkdir(parents=True)
    for r in (1, 2, 3):
        (d / f"rep{r}.Llog").write_text(open(os.path.join(FIX, f"rep{r}.Llog")).read())
    tsv = aggregate(str(tmp_path / "divergence"), [("human", 0)])
    body = open(tsv).read()
    assert "human" in body and "Sort+align" in body
    assert os.path.exists(os.path.join(str(tmp_path / "divergence"), "divergence_phase_share.png"))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd docs/benchmark_matrix && python -m pytest tests/test_parse_llog.py::test_aggregate_one_point -q`
Expected: FAIL — `ImportError: cannot import name 'aggregate'`.

- [ ] **Step 3: Implement aggregation + figure + CLI**

Append to `docs/benchmark_matrix/aggregate_matrix.py`:
```python
import os, glob, statistics as _st, sys

def point_median(logs_dir):
    reps = [parse_llog(p) for p in sorted(glob.glob(os.path.join(logs_dir, "rep*.Llog")))]
    if not reps:
        raise FileNotFoundError(f"no rep*.Llog in {logs_dir}")
    med = {}
    for s in STAGES:
        med[s] = _st.median(r[s] for r in reps)
    med["rss_mb"] = _st.median(r["rss_mb"] for r in reps)
    med["n_aln"]  = int(_st.median(r["n_aln"] for r in reps))
    med["ave_len"] = int(_st.median(r["ave_len"] for r in reps))
    med["total"]  = sum(med[s] for s in STAGES)
    return med

def aggregate(base_dir, points):
    rows = []
    for label, rank in sorted(points, key=lambda x: x[1]):
        m = point_median(os.path.join(base_dir, label, "logs"))
        rows.append((label, rank, m))
    tsv = os.path.join(base_dir, "results.tsv")
    with open(tsv, "w") as f:
        f.write("label\trank\t" + "\t".join(STAGES) +
                "\ttotal_s\tsort_align_share\trss_mb\tn_aln\tave_len\n")
        for label, rank, m in rows:
            f.write(f"{label}\t{rank}\t" + "\t".join(f"{m[s]:.1f}" for s in STAGES) +
                    f"\t{m['total']:.1f}\t{m['Sort+align']/m['total']:.3f}"
                    f"\t{m['rss_mb']:.0f}\t{m['n_aln']}\t{m['ave_len']}\n")
    _plot(base_dir, rows)
    return tsv

def _plot(base_dir, rows):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt, numpy as np
    C = {"GDB": "#9e9e9e", "GIX": "#ff9800", "Seed merge": "#1976d2", "Sort+align": "#2e7d32"}
    labels = [r[0] for r in rows]
    x = np.arange(len(rows)); bottom = np.zeros(len(rows))
    fig, ax = plt.subplots(figsize=(1.6 * len(rows) + 3, 4.5))
    for s in STAGES:
        share = np.array([r[2][s] / r[2]["total"] for r in rows])
        ax.bar(x, share, bottom=bottom, color=C[s], label=s, width=0.62)
        bottom += share
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("share of runtime"); ax.set_ylim(0, 1)
    ax.set_title("FastGA per-phase share vs divergence (most similar → most divergent)")
    ax.legend(loc="upper right", framealpha=0.95)
    fig.tight_layout()
    fig.savefig(os.path.join(base_dir, "divergence_phase_share.png"), dpi=140)

if __name__ == "__main__":
    base = sys.argv[1]
    pts = [(a.split(":")[0], int(a.split(":")[1])) for a in sys.argv[2:]]
    print("wrote", aggregate(base, pts))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd docs/benchmark_matrix && python -m pytest tests/test_parse_llog.py -q`
Expected: PASS (5 passed).

- [ ] **Step 5: Commit**

```bash
git add docs/benchmark_matrix/aggregate_matrix.py docs/benchmark_matrix/tests/test_parse_llog.py
git commit -m "feat(matrix): cross-dataset aggregation + stacked-share figure"
```

---

### Task 4: `run_pair.sh` — one genome pair with disk guard

**Files:**
- Create: `docs/benchmark_matrix/run_pair.sh`

**Interfaces:**
- Consumes: `BL` env (baseline bin dir); positional `<label> <G1_fasta> <G2_fasta>`; env `MIN_FREE_GB` (default 150), `REPS` (default 3), `SCRATCH` (default `/scratch/jl4257/matrix_run`), `OUT` (default `docs/benchmark_matrix/divergence`).
- Produces: `OUT/<label>/logs/rep*.Llog`, `OUT/<label>/time/rep*.time`, `OUT/<label>/meta.tsv`. Exit 3 if disk guard trips; exit 2 on bad args.

- [ ] **Step 1: Write the script**

Create `docs/benchmark_matrix/run_pair.sh`:
```bash
#!/bin/bash
# Profile one FastGA genome pair (T=32) -> per-phase -L logs + /usr/bin/time.
# Usage: BL=<bin dir> ./run_pair.sh <label> <G1.fna> <G2.fna>
#   --dry-run : print the plan and exit without running FastGA.
set -u
DRY=0; [ "${1:-}" = "--dry-run" ] && { DRY=1; shift; }
LABEL="${1:?usage: run_pair.sh <label> <G1> <G2>}"
G1="${2:?need G1 fasta}"; G2="${3:?need G2 fasta}"
BL="${BL:?set BL=<baseline bin dir>}"
T=32; REPS="${REPS:-3}"; MIN_FREE_GB="${MIN_FREE_GB:-150}"
SCRATCH="${SCRATCH:-/scratch/jl4257/matrix_run}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${OUT:-$HERE/divergence}/$LABEL"

for f in "$G1" "$G2"; do [ -r "$f" ] || { echo "missing genome: $f" >&2; exit 2; }; done
[ -x "$BL/FastGA" ] || { echo "no FastGA in BL=$BL" >&2; exit 2; }

disk_guard() {  # refuse if free space on $1 < MIN_FREE_GB
  mkdir -p "$1"
  local avail; avail=$(df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc '0-9')
  [ -n "$avail" ] || { echo "disk guard: cannot stat $1" >&2; exit 3; }
  if [ "$avail" -lt "$MIN_FREE_GB" ]; then
    echo "disk guard: $1 has ${avail}GB free < ${MIN_FREE_GB}GB required — refusing" >&2
    exit 3
  fi
  echo "disk guard: ${avail}GB free on $1 (need ${MIN_FREE_GB}) — OK"
}

echo "=== pair '$LABEL': $(basename "$G1") x $(basename "$G2"), T=$T, $REPS reps ==="
disk_guard "$SCRATCH"
if [ "$DRY" = 1 ]; then
  echo "[dry-run] would run $REPS reps of: $BL/FastGA -v -T$T -P<tmp> -L:<log> g1 g2"
  echo "[dry-run] outputs -> $OUT/{logs,time}/"; exit 0
fi

mkdir -p "$OUT/logs" "$OUT/time"
# record measured genome sizes (bp) for the secondary divergence proxy
size_bp() { grep -v '^>' "$1" | tr -d '\n' | wc -c; }
{ echo -e "genome\tsize_bp"; echo -e "G1\t$(size_bp "$G1")"; echo -e "G2\t$(size_bp "$G2")"; } > "$OUT/meta.tsv"

for rep in $(seq 1 "$REPS"); do
  W="$SCRATCH/work"; TMP="$SCRATCH/tmp"; rm -rf "$W" "$TMP"; mkdir -p "$W" "$TMP"
  ln -sf "$G1" "$W/g1.fna"; ln -sf "$G2" "$W/g2.fna"
  ( cd "$W" && /usr/bin/time -v -o "$OUT/time/rep${rep}.time" \
      "$BL/FastGA" -v -T$T -P"$TMP" -L:"$OUT/logs/rep${rep}.Llog" g1.fna g2.fna \
      > /dev/null 2> "$OUT/logs/rep${rep}.stderr" )
  echo "[$(date +%H:%M:%S)] $LABEL rep$rep done"
  rm -rf "$W" "$TMP"
done
echo "PAIR_DONE $LABEL"
```

- [ ] **Step 2: Make executable and test the disk guard trips**

Run:
```bash
chmod +x docs/benchmark_matrix/run_pair.sh
cd docs/benchmark_matrix
BL=/tmp MIN_FREE_GB=999999 SCRATCH=/scratch/jl4257/_guardtest \
  bash -c 'touch /tmp/FastGA; chmod +x /tmp/FastGA;
           ./run_pair.sh --dry-run demo /etc/hostname /etc/hostname'; echo "exit=$?"
```
Expected: prints `disk guard: … < 999999GB required — refusing` and `exit=3`.

- [ ] **Step 3: Test dry-run passes the guard with a low threshold**

Run:
```bash
cd docs/benchmark_matrix
BL=/tmp MIN_FREE_GB=1 SCRATCH=/scratch/jl4257/_guardtest \
  ./run_pair.sh --dry-run demo /etc/hostname /etc/hostname; echo "exit=$?"
rm -rf /scratch/jl4257/_guardtest
```
Expected: prints `disk guard: …GB free … OK`, `[dry-run] would run 3 reps…`, `exit=0`.

- [ ] **Step 4: Commit**

```bash
git add docs/benchmark_matrix/run_pair.sh
git commit -m "feat(matrix): generalized per-pair profiling harness with disk guard"
```

---

### Task 5: `fetch_genomes.sh` — download query genomes with disk precheck

**Files:**
- Create: `docs/benchmark_matrix/fetch_genomes.sh`

**Interfaces:**
- Consumes: env `DEST` (default `/scratch/jl4257/seq_align/fastga_datasets`), `MIN_FREE_GB` (default 150); flag `--dry-run`.
- Produces: downloaded + gunzipped FASTA under `DEST/<name>/`. Idempotent (skips if the target `.fna` already exists). Exit 3 if disk guard trips.

- [ ] **Step 1: Write the script**

Create `docs/benchmark_matrix/fetch_genomes.sh`:
```bash
#!/bin/bash
# Download the divergence-axis query genomes (chimp, siamang, pig, mouse) to /scratch.
# CHM13 + GRCh38 are already local. Usage: ./fetch_genomes.sh [--dry-run]
set -u
DRY=0; [ "${1:-}" = "--dry-run" ] && DRY=1
DEST="${DEST:-/scratch/jl4257/seq_align/fastga_datasets}"
MIN_FREE_GB="${MIN_FREE_GB:-150}"

# name | approx download GB | URL  (NCBI datasets accessions; pig from GigaDB)
GENOMES=(
"chimpanzee|1.0|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna.gz"
"siamang|1.0|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/878/055/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.fna.gz"
"mouse|0.9|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/188/535/GCA_964188535.1/GCA_964188535.1_genomic.fna.gz"
)
# NOTE: pig is on GigaDB (dataset 102692), not NCBI — its URL must be confirmed manually.
# Add once verified:  "pig|0.8|<gigadb direct fna.gz url>"

avail=$(df -BG --output=avail "$(dirname "$DEST")" 2>/dev/null | tail -1 | tr -dc '0-9')
echo "dest=$DEST  free=${avail:-?}GB  (need >= ${MIN_FREE_GB})"
[ -n "$avail" ] && [ "$avail" -lt "$MIN_FREE_GB" ] && { echo "disk guard: refusing" >&2; exit 3; }

for spec in "${GENOMES[@]}"; do
  IFS='|' read -r name gb url <<< "$spec"
  out="$DEST/$name"; fna="$out/$(basename "${url%.gz}")"
  if ls "$out"/*.fna >/dev/null 2>&1; then echo "skip $name (present)"; continue; fi
  echo ">>> $name  (~${gb}GB)  <- $url"
  [ "$DRY" = 1 ] && continue
  mkdir -p "$out"
  curl -fL --retry 3 -o "$fna.gz" "$url" || { echo "download failed: $name" >&2; exit 4; }
  gunzip -f "$fna.gz"
  echo "    -> $fna ($(du -h "$fna" | cut -f1))"
done
echo "FETCH_DONE"
```

- [ ] **Step 2: Test dry-run lists genomes without downloading**

Run:
```bash
chmod +x docs/benchmark_matrix/fetch_genomes.sh
cd docs/benchmark_matrix
DEST=/scratch/jl4257/seq_align/fastga_datasets MIN_FREE_GB=1 ./fetch_genomes.sh --dry-run
echo "exit=$?"
```
Expected: prints `dest=… free=…GB`, three `>>> chimpanzee/siamang/mouse … <- https://…` lines (GRCh38/CHM13 not listed), `FETCH_DONE`, `exit=0`. No files downloaded.

- [ ] **Step 3: Test disk guard refuses**

Run:
```bash
cd docs/benchmark_matrix
MIN_FREE_GB=999999 ./fetch_genomes.sh --dry-run; echo "exit=$?"
```
Expected: `disk guard: refusing` and `exit=3`.

- [ ] **Step 4: Commit**

```bash
git add docs/benchmark_matrix/fetch_genomes.sh
git commit -m "feat(matrix): query-genome downloader with disk precheck"
```

---

### Task 6: Orchestrator, README, self-review, handoff

**Files:**
- Create: `docs/benchmark_matrix/run_divergence_axis.sh`
- Create: `docs/benchmark_matrix/README.md`
- Modify: `.gitignore` (ignore heavy per-rep logs, keep results.tsv + png)

**Interfaces:**
- Consumes: `run_pair.sh`, `fetch_genomes.sh`, `aggregate_matrix.py`, the local CHM13/GRCh38 FASTA, env `BL`.
- Produces: the full axis run + `divergence/results.tsv` + `divergence/divergence_phase_share.png`.

- [ ] **Step 1: Write the orchestrator**

Create `docs/benchmark_matrix/run_divergence_axis.sh`:
```bash
#!/bin/bash
# Run the divergence axis: CHM13 x {GRCh38, chimp, siamang, pig, mouse}, T=32.
# GRCh38 point reuses the existing human baseline logs unless RERUN_HUMAN=1.
# Usage: BL=<bin dir> ./run_divergence_axis.sh [--dry-run]
set -u
DRY="${1:-}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DS="${DS:-/scratch/jl4257/seq_align/fastga_datasets}"
CHM13="$DS/CHM13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
declare -A Q=(
  [human]="$DS/GRCh38/GCF_000001405.40_GRCh38.p14_genomic.fna"
  [chimpanzee]="$DS/chimpanzee"/*.fna
  [siamang]="$DS/siamang"/*.fna
  [pig]="$DS/pig"/*.fna
  [mouse]="$DS/mouse"/*.fna
)
declare -A RANK=( [human]=0 [chimpanzee]=1 [siamang]=2 [pig]=3 [mouse]=4 )
ORDER=(human chimpanzee siamang pig mouse)

# reuse the committed human logs as the 'human' point
if [ "${RERUN_HUMAN:-0}" != 1 ]; then
  mkdir -p "$HERE/divergence/human/logs"
  cp -n "$HERE/../benchmark_baseline/human/perf_data/logs/rep"*.Llog \
        "$HERE/divergence/human/logs/" 2>/dev/null || true
fi

for lbl in "${ORDER[@]}"; do
  g2=$(ls ${Q[$lbl]} 2>/dev/null | head -1)
  if [ "$lbl" = human ] && [ "${RERUN_HUMAN:-0}" != 1 ]; then
    echo "== $lbl: reuse existing baseline logs =="; continue
  fi
  [ -r "$g2" ] || { echo "== $lbl: genome missing (run fetch_genomes.sh) — skip =="; continue; }
  echo "== $lbl (rank ${RANK[$lbl]}) =="
  OUT="$HERE/divergence" "$HERE/run_pair.sh" $DRY "$lbl" "$CHM13" "$g2"
done

pts=""; for lbl in "${ORDER[@]}"; do
  [ -d "$HERE/divergence/$lbl/logs" ] && pts="$pts $lbl:${RANK[$lbl]}"; done
echo "== aggregating:$pts =="
[ "$DRY" = "--dry-run" ] || python "$HERE/aggregate_matrix.py" "$HERE/divergence" $pts
```

- [ ] **Step 2: Test orchestrator dry-run**

Run:
```bash
chmod +x docs/benchmark_matrix/run_divergence_axis.sh
cd docs/benchmark_matrix
BL=/tmp MIN_FREE_GB=1 ./run_divergence_axis.sh --dry-run; echo "exit=$?"
```
Expected: `== human: reuse existing baseline logs ==`, then for chimp/siamang/pig/mouse either `genome missing … skip` (not yet downloaded) or a dry-run guard line; ends with `== aggregating: human:0 … ==`; `exit=0`.

- [ ] **Step 3: Smoke-test the full aggregation on the human point**

Run:
```bash
cd docs/benchmark_matrix
mkdir -p divergence/human/logs
cp -n ../benchmark_baseline/human/perf_data/logs/rep*.Llog divergence/human/logs/
python aggregate_matrix.py divergence human:0
column -t divergence/results.tsv
```
Expected: `results.tsv` shows the `human` row with `sort_align_share` ≈ `0.81`; `divergence/divergence_phase_share.png` exists.

- [ ] **Step 4: Write the README**

Create `docs/benchmark_matrix/README.md`:
````markdown
# FastGA bottleneck profiling matrix

Characterizes how FastGA's per-phase bottleneck shifts across datasets. **v1 = divergence
axis** (CHM13 × {GRCh38, chimp, siamang, pig, mouse}, ~3 Gbp fixed, T=32). Design spec:
`docs/superpowers/specs/2026-07-11-fastga-bottleneck-profiling-matrix-design.md`.
Dataset reference: `datasets_inventory.md`.

## Run it (all I/O on /scratch; disk guard at 150 GB free)

```bash
export BL=/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA-baseline   # baseline bin
./fetch_genomes.sh                 # download chimp, siamang, mouse (pig: add GigaDB URL first)
BL=$BL ./run_divergence_axis.sh    # profiles each pair, reuses committed human logs
column -t divergence/results.tsv
open divergence/divergence_phase_share.png
```

Per-pair peak ≈ 80 GB transient (built, then cleaned each rep); persistent downloads ≈ 10–25 GB.
`pig` is on GigaDB (dataset 102692), not NCBI — confirm its `.fna.gz` URL and add the row to
`fetch_genomes.sh` before running the pig point.

## Hypothesis

`sort+align` share falls from ~81% (human, most similar) toward the divergent end as less
sequence is alignable; GIX/seed shares rise. See the stacked figure.

## Results

_(fill in after the runs: paste the `results.tsv` table and the figure.)_
````

- [ ] **Step 5: Ignore heavy artifacts, keep results**

Run:
```bash
cd "$(git rev-parse --show-toplevel)"
cat >> .gitignore <<'EOF'

# profiling-matrix raw per-rep artifacts (keep aggregated results + figure)
docs/benchmark_matrix/divergence/*/logs/
docs/benchmark_matrix/divergence/*/time/
docs/benchmark_matrix/divergence/*/*.stderr
EOF
```

- [ ] **Step 6: Self-review vs spec, then commit**

Verify each spec section maps to a task: inventory (Task 1), aggregator+figure (Tasks 2–3), harness+disk-guard (Task 4), downloader (Task 5), orchestrator+README+disk-plan (Task 6). Confirm no `divergence/*/logs` were staged.

Run:
```bash
cd "$(git rev-parse --show-toplevel)"
git status --short
git add docs/benchmark_matrix/run_divergence_axis.sh docs/benchmark_matrix/README.md .gitignore \
        docs/benchmark_matrix/divergence/results.tsv docs/benchmark_matrix/divergence/divergence_phase_share.png
git commit -m "feat(matrix): divergence-axis orchestrator, README, results scaffold"
```
Expected: clean commit; no `*/logs/` files tracked.

---

## Self-Review (plan vs spec)

- **Spec §3 inventory** → Task 1 (with accession-presence validation). ✓
- **Spec §4 divergence-axis 5 points** → Task 6 orchestrator (`ORDER`/`RANK`), reusing human. ✓
- **Spec §5 harness reuse / metrics / T=32 / divergence proxy** → Tasks 2–4 (parser reproduces 81%; `run_pair.sh` T=32, records genome sizes; ordinal rank primary, aligned-bp proxy via `n_aln×ave_len`). ✓
- **Spec §6 disk plan / guard** → Task 4 `disk_guard` + Task 5 precheck, `MIN_FREE_GB=150`, all under `/scratch`. ✓
- **Spec §7 scripts** → fetch (T5), run_pair (T4), run_divergence_axis (T6), aggregate (T2–3), README (T6). ✓
- **Spec §8 layout** → matches created paths. ✓
- **Spec §9 success criteria** → parser unit test + dry-run tests + human-point aggregation smoke test. ✓
- **Type consistency:** `parse_llog`→`point_median`→`aggregate` keys (`STAGES`, `total`, `rss_mb`, `n_aln`) consistent across Tasks 2–3. ✓
- **Placeholder scan:** the only deferred items are the pig GigaDB URL and DToL Supp-S3 accessions — both explicitly flagged, out of v1's critical path (chimp/siamang/mouse suffice to show the trend). ✓
