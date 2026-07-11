# FastGA Bottleneck Profiling Matrix — v1 (divergence axis)

- **Date:** 2026-07-11
- **Branch:** `agent-optimization-wt` (based on `agent-optimization`)
- **Status:** Design approved; spec pending user review → writing-plans
- **Owner:** jl4257
- **Deliverable mode:** **Plan + scripts only.** This spec produces reusable scripts and a
  documented methodology; the actual runs are executed by the user (no SLURM on this node;
  runs go on the local `/scratch` disk). Disk budget is a first-class constraint.

## 1. Motivation & hypothesis

FastGA is a new, architecturally unusual, and increasingly important whole-genome aligner
(cache-coherent, sort/merge-based, memory-system-optimized). We are evaluating it as a target
for hardware (GPU/FPGA) acceleration and want a paper-quality **workload characterization** as
the foundation.

Existing profiling on this branch covers exactly one real dataset: GRCh38 × CHM13 (human,
~3.1 Gbp, T=32), where the runtime splits GDB 3% / GIX 15% / seed-merge 1% / **sort+align 81%**.
That 81% is the aligner (Myers' O(nd) local wavefront, the DALIGNER LA-finder — **not** the
Marco-Sola WFA, so off-the-shelf WFA-GPU/WFA-FPGA are not drop-ins).

**Hypothesis:** the per-phase bottleneck is **not universal** — it shifts structurally with
dataset characteristics. The paper's own runtimes already hint at it (mammal→CHM13 runtime falls
70.5→16.6 min as evolutionary distance grows, because closely related genomes have far more
alignable sequence and thus more alignment work).

**v1 tests one axis only — divergence** — at fixed genome size (~3 Gbp): as the query diverges
from CHM13, confirm the `sort+align` share falls monotonically from ~81% and the
seeding/indexing share rises. This is the first block of the "characterize the workload →
argue where acceleration matters" narrative.

## 2. Scope

### In scope (v1)
1. **Dataset inventory** of every real dataset the FastGA paper used (reference for future axes).
2. **Divergence-axis runs**: CHM13 × {GRCh38, chimpanzee, siamang, pig, mouse}, T=32, per-phase
   profiling, and a cross-dataset aggregation (table + the key figure).
3. Reusable, parameterized scripts (fetch / run / aggregate) with a disk guard.

### Out of scope (YAGNI — recorded for future phases only)
- Size axis (DToL species 405 Mb → 24 Gb), the newt 24 Gbp extreme, within- vs between-species
  contrasts, thread scaling, T=16 paper-parity re-runs, the simulated-genome sensitivity sweep.
- Any GPU/FPGA implementation. This spec is characterization only.

## 3. Deliverable 1 — dataset inventory (`docs/benchmark_matrix/datasets_inventory.md`)

A table recording, for every dataset in the FastGA paper (Myers, Durbin & Zhou, bioRxiv 2025,
§5), the fields: taxon/species, genus, accession/source URL, genome size (Mb), role
(reference / query / within-species pair / between-species pair), and the paper's reported
FastGA CPU-time / memory / coverage.

**Mammalian → CHM13 (§5.2), sizes and accessions from the paper:**

| Role | Species | Size (Mb) | Accession / source |
|---|---|--:|---|
| reference | CHM13 v2.0 (T2T) | 3,117 | GCF_009914755.1 |
| query | human GRCh38.p14 | 3,298 | GCF_000001405.40 |
| query | chimpanzee | 3,178 | GCF_028858775.2 |
| query | siamang | 3,263 | GCF_028878055.3 |
| query | pig | 2,612 | gigadb.org/dataset/102692 |
| query | mouse | 2,731 | GCA_964188535.1 |

**DToL 12 species (§5.3) — six genera, within- + between-species (future size axis):** Insect
*Acronicta* (A. psi ~405 Mb, A. aceris ~466 Mb), Fish *Thunnus* (T. albacares ~792 Mb,
T. maccoyii ~782 Mb), Bird *Ammospiza* (A. caudacuta ~1,241 Mb, A. maritima ~1,398 Mb), Reptile
*Vipera* (V. latastei ~1,632 Mb, V. berus ~1,695 Mb), Mammal *Molossus* (M. alvarezi ~2,505 Mb,
M. nigricans ~2,567 Mb), Amphibian *Lissotriton* (L. vulgaris ~24,226 Mb, L. helveticus
~23,170 Mb). **Per-species accessions live in the paper's Supplementary Table S3** and are
retrieved when the size axis is scheduled — not needed for v1.

**Simulated genomes (§5.1, future sensitivity axis):** pair of ~84 Mb genomes built from 10 kb
blocks, each a similarity region (length 100 bp–5 kb) at divergence 1%–65% (SNV 80% / ins 10% /
del 10% on B) followed by random sequence, blocks shuffled. Recorded as construction parameters,
not a download.

## 4. Deliverable 2 — divergence-axis runs

Reference = CHM13; queries by increasing evolutionary distance (paper's exact gradient):

| # | Query | Size | Status |
|--:|---|--:|---|
| 1 | GRCh38 | 3.3 Gb | **already profiled (T=32, 81% align) — reuse existing result** |
| 2 | chimpanzee | 3.2 Gb | download + run |
| 3 | siamang | 3.3 Gb | download + run (optional 5th point; smooths the curve) |
| 4 | pig | 2.6 Gb | download + run |
| 5 | mouse | 2.7 Gb | download + run |

## 5. Methodology

- **Harness reuse.** Generalize the existing `docs/benchmark_baseline/human/run_human.sh` +
  `analyze_perf.py` (they are already parameterized by `G1`/`G2`/`BL`) into:
  - `run_pair.sh <label> <G1_fasta> <G2_fasta>` — one genome pair: 3 perf reps capturing
    FastGA's `-L` per-phase log (user/sys/wall/CPU% per phase) + `/usr/bin/time -v` (peak RSS).
    Given FASTA inputs, FastGA auto-invokes FAtoGDB + GIXmake, so the `-L` log yields the same
    GDB / GIX / seed-merge / sort+align breakdown unchanged.
  - Storage timeline is **skipped in v1** (near-constant at fixed size; it belongs to the size
    axis).
- **Binary.** Same baseline build used for the human baseline (upstream `ddeea32` FastGA/GIXmake
  + the `5671357` FAtoGDB that fixes the CHM13-GDB segfault). Passed via `BL=<bin dir>`. Rationale:
  keep v1 comparable to the existing 81%-align human datapoint. (agent-optimization binaries are a
  separate future comparison, not v1.)
- **Threads:** **T=32** (matches the existing human datapoint and our own baseline; paper's T=16
  parity is a future add-on).
- **Metrics per point:** per-phase wall / share / CPU% (median of 3 reps); peak RSS; and — as the
  divergence quantifier for the x-axis — **fraction of query genome covered by alignments**
  (primary; the paper's own metric), plus #alignments and total alignment size, parsed from the
  FastGA stderr/`-L` summary. Ordinal phylogenetic rank (human<chimp<siamang<pig<mouse) recorded
  as a fallback x-axis.
- **Key figure:** x = divergence proxy (coverage fraction, most-similar → most-divergent),
  y = stacked per-phase share — showing the `sort+align` band contracting from ~81% and the
  GIX/seed bands expanding.
- **Aggregator** `aggregate_matrix.py`: reads each point's parsed results, emits the cross-dataset
  table (per-phase wall + share + CPU% + RSS + coverage) and the stacked-share figure into
  `docs/benchmark_matrix/divergence/`.

## 6. Disk plan (hard constraint)

- **All I/O on `/scratch/jl4257`** (local `VG00-scratch`, 877 GB free, **no per-user quota**).
  `/home/jl4257` (224 GB free, 90% full, NFS) and the shared-users NFS (74 GB headroom, quota)
  are **not** used for genome data.
- **Per-run peak (single pair, ~3 Gbp each):** ~2×GIX (~33 GB each) + temp ≈ **~80 GB transient**,
  built then cleaned per run. Runs are serial, one pair at a time.
- **Persistent footprint after all downloads:** 4 new query FASTA (gzipped ~1 GB each) + their
  GDBs (~1.5 GB each) ≈ **~10–25 GB**, on top of the existing ~71 GB human data. Far below 877 GB.
- **Disk guard:** `run_pair.sh` checks `df` free space on the scratch target before each run and
  **refuses (non-zero exit, clear message)** if free space < a configurable threshold
  (default 150 GB) — never risk filling `/scratch`. `fetch_genomes.sh` reports projected vs.
  available space before downloading.

## 7. Scripts to produce (interfaces)

Under `docs/benchmark_matrix/`:

1. `fetch_genomes.sh` — download the 4 query genomes (chimp, siamang, pig, mouse) by accession to
   `/scratch/jl4257/seq_align/fastga_datasets/<name>/`; verify size; report disk before/after.
   Idempotent (skip if present). CHM13/GRCh38 already local.
2. `run_pair.sh <label> <G1> <G2>` — generalized perf harness (3 reps, `-L` + `/usr/bin/time`),
   writes `divergence/<label>/{logs,time}/`. Includes the disk guard.
3. `run_divergence_axis.sh` — loops the 5 points calling `run_pair.sh`; the GRCh38 point reuses the
   existing human result instead of re-running (configurable).
4. `aggregate_matrix.py` — generalizes `analyze_perf.py`; per-point parse → cross-dataset table +
   stacked-share figure.
5. `docs/benchmark_matrix/README.md` — how to run, disk notes, and the results write-up.

## 8. Directory layout

```
docs/benchmark_matrix/
  datasets_inventory.md          # Deliverable 1
  README.md                      # method + results write-up
  fetch_genomes.sh
  run_pair.sh
  run_divergence_axis.sh
  aggregate_matrix.py
  divergence/
    <label>/{logs,time}/         # per-point raw
    results.tsv                  # aggregated
    divergence_phase_share.png   # the key figure
```

## 9. Success criteria

- Inventory doc records all paper datasets with accessions/sizes for the mammals (DToL flagged to
  Supp S3).
- Scripts run end-to-end on one pair with the disk guard active, producing a parsed per-phase
  breakdown, on a dry-run/smoke test the user can trigger.
- The aggregated figure/table can express the v1 hypothesis: `sort+align` share vs divergence,
  with the GRCh38 (81%) point anchored and ≥3 more-divergent points to show the trend.
- Zero risk to `/scratch`: no run proceeds under the free-space threshold.

## 10. Risks / open questions

- **Genome availability/format:** pig is from GigaDB (not NCBI); URL/format must be confirmed in
  `fetch_genomes.sh`. Mouse `GCA_964188535.1` is a specific assembly — pin the exact accession.
- **Coverage as divergence proxy:** coverage conflates divergence with assembly completeness and
  repeat content; acceptable for v1 (single size, closely controlled), revisit for other axes.
- **Reuse of the human point:** assumes the existing human run used the same baseline binary and
  T=32 — it did (per `docs/benchmark_baseline/human/`).
