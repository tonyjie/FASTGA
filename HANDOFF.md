# FastGA fork — handoff

This fork of [FastGA](https://github.com/thegenemyers/FASTGA) carries a **CPU load-balancing
optimization of the `sort+align` phase** (the aligner's dominant stage). On a human×human run it
takes that phase from **6 → 28 of 32 cores, 5.33× faster, byte-for-byte identical output**. This
document tells you what's in the repo, where the datasets live on the server, and how to build,
run, and reproduce the result so you can build on it.

Repo: **`git@github.com:tonyjie/FASTGA.git`** (fork of `thegenemyers/FASTGA`).

---

## 1. Branches — what's here and how to use them

### `main` — clean upstream mirror
`main` is a faithful mirror of upstream FastGA plus a single commit adding `.gitignore` (build
artifacts, indices, local docs). Use it as the **stock baseline** — the reference the optimization
is measured against, and the branch to rebase future upstream syncs onto. Tip commit: `10ebff7`.

### `workstealing-v3c` — the optimization (use this) ⭐
The load-balancing work, four incremental **bit-exact** commits on top of `main`, all in a single
file (`FastGA.c`, functions `pair_sort_search` / `search_seeds`). The actual alignment routine
(`align_contigs`) is **never touched** — only *which thread runs which contig-pair* changes.

```
e40e6c5  v3c  reimport+sort write directly into the resident buffer (drop double buffer)  ← tip
db45c59  v3a  global interleaved pool, both strand passes resident, no barrier
dd4f6f0  v2   one pool per strand-pass (barrier only between the 2 passes)
59a5db3  v1   contig-pair shared work-pool, within each part (barrier kept)
10ebff7  (= main, stock upstream baseline)
```

Measured (GRCh38 × CHM13, T=32, cached indices, on `en-ec-zhang-x4`):

| version | sort+align wall | cores | peak RAM | speedup |
|---|--:|--:|--:|--:|
| stock (`main`)   | 491.6 s | 6.0  | 19 GiB | 1.00× |
| v1               | 423.5 s | 7.8  | 21 GiB | 1.16× |
| v2               | 150.5 s | 17.1 | 33 GiB | 3.27× |
| v3a              | 101.9 s | 24.9 | 35 GiB | 4.82× |
| **v3c (tip)**    | **92.3 s** | **28.0** | **27 GiB** | **5.33×** |

Whole FastGA with cached indices: **5.06×**; from FASTA incl. index build: **2.95×** (Amdahl — the
fixed ~113 s GDB+GIX build then dominates). This is also **PR #2** (open draft):
<https://github.com/tonyjie/FASTGA/pull/2>.

**The idea in one paragraph.** Stock FastGA splits the genome into ~66 "parts" by *length* and hands
each thread a fixed range of contig-pairs, balancing by *seed count* — but the real cost is wavefront
*depth*, uncorrelated with seed count, so a few "monster" contig-pairs (up to ~58 s each) strand
threads while the rest idle. The fix replaces the static split with a **shared atomic work-pool**:
all contig-pairs from both strand passes go into one job list (`build_jobs()`), and 28 threads pull
the next job via one atomic counter — self-balancing, no per-part barrier. v3c additionally sorts
directly into the resident seed buffer, dropping a redundant ~9 GiB copy (faster *and* less RAM).

### Other branches (context, not needed to reproduce)
`workstealing-v1` = the PR #2 head (identical tip to `workstealing-v3c`). `workstealing-v3a`,
`workstealing-v3b` = checkpoints (v3b tried a decompressed-sequence cache; it gave **zero** speedup
and was dropped — kept for comparison). The `agent-*`, `optimize-memory`, `fused-scan-once`,
`agentic-steps` branches are older/separate exploration lines (GPU, memory, etc.) — ignore for this work.

---

## 2. The datasets on this server

**Server:** `en-ec-zhang-x4`. The data is on **node-local `/scratch`**, so you must be logged into
**this specific host** to see it. You do **not** need the `jl4257` account — the tree is
world-readable (`0755` dirs, `0644` files); read it directly with your own account.

**Location:** `/scratch/jl4257/seq_align/fastga_datasets/` (~80 GB total), two human genomes:

| genome | prefix |
|---|---|
| GRCh38 | `/scratch/jl4257/seq_align/fastga_datasets/GRCh38/GCF_000001405.40_GRCh38.p14_genomic` |
| CHM13 (T2T) | `/scratch/jl4257/seq_align/fastga_datasets/CHM13/GCF_009914755.1_T2T-CHM13v2.0_genomic` |

Each prefix already has the **FASTA and the prebuilt FastGA indices**, so you can skip index
building entirely and point FastGA straight at the `.gix`:

- `<prefix>.fna` — source FASTA (~3.3 GB)
- `<prefix>.1gdb`, `.1ano`, `.<prefix>.bps` — the **GDB** (FastGA's genome database)
- `<prefix>.gix`, `.<prefix>.ktab.1 … .ktab.32` — the **GIX** (40-mer k-mer index)

> You can **read** these but not write into `jl4257`'s directory. That's fine: FastGA only *reads*
> the `.gdb`/`.gix`; it writes temporaries to a `-P <dir>` you choose and the result to a `-1:<path>`
> you choose. Point both at your own directory (e.g. `/scratch/$USER/...`).
>
> If you ever need to rebuild the CHM13 index yourself, note `FAtoGDB` in recent upstream **segfaults
> on soft-masked genomes** (CHM13 is soft-masked); work around it with an older `FAtoGDB` or by
> upper-casing the FASTA first. You won't hit this if you reuse the prebuilt indices above.

---

## 3. Quickstart — build, run, verify

All from a checkout of this fork on `en-ec-zhang-x4`. Toolchain present: `gcc 11`, `make`, `zlib`.

```bash
# 0. clone + pick the optimized branch
git clone git@github.com:tonyjie/FASTGA.git
cd FASTGA
git checkout workstealing-v3c

# 1. build (the whole suite; FastGA + ALNtoPAF are what you need here)
make FastGA ALNtoPAF          # ~10 s; a few benign "may be used uninitialized" warnings are normal

# 2. set up paths (use YOUR OWN scratch dir for temp + output)
G1=/scratch/jl4257/seq_align/fastga_datasets/GRCh38/GCF_000001405.40_GRCh38.p14_genomic
G2=/scratch/jl4257/seq_align/fastga_datasets/CHM13/GCF_009914755.1_T2T-CHM13v2.0_genomic
T=/scratch/$USER/fastga_tmp; mkdir -p $T

# 3. run the aligner on the two human genomes, reusing the prebuilt indices
#    -T32 threads, -k keep indices, -P temp dir, -L: log, -1: output alignment file
./FastGA -T32 -k -P$T -L:$T/run.log -1:$T/out.1aln $G1.gix $G2.gix

# 4. inspect the phase timing (this is the number the optimization moves)
grep -A1 'Seed sort and alignment' $T/run.log   # "Resources for phase: ... w  <CPU%>"
```

Expect the `sort+align` phase around **~92 s wall / ~2800% CPU (28 cores)** on this branch.

### Verify correctness (the regression gate)
Output is **bit-exact** vs stock upstream. The gate hashes the *sorted* alignments (order-independent,
so thread scheduling doesn't matter):

```bash
./ALNtoPAF $T/out.1aln | sort | sha256sum
# MUST equal:  0ef422c2b9d54aae8103919059f8b8bb49a1f77781ed1f5dba8f745231212e10
# and the log must report:  518037 non-redundant alignments
```

### Reproduce the speedup vs stock
Build & run `main` the same way, compare the `sort+align` phase wall/CPU%:

```bash
git worktree add ../FASTGA-stock main && cd ../FASTGA-stock && make FastGA
./FastGA -T32 -k -P$T -L:$T/stock.log -1:$T/stock.1aln $G1.gix $G2.gix
grep -A1 'Seed sort and alignment' $T/stock.log   # ~490 s / ~600% (6 cores)
```
Both `out.1aln` must hash to the same `0ef422c2…` value.

---

## 4. Where the deeper analysis lives

The full design rationale, per-version measurements, memory analysis, and a GPU-parallelism study
are written up outside this repo (they were the working notes for this optimization):

- **Task write-up + progress log:** `project-hub/projects/FastGA/2026-07-18-workload-imbalance.md`
- **Design/spec/plan + run data + scripts:** `project-hub/projects/FastGA/assets/workload-imbalance/`
  (per-version `.log`/`.time`, memory traces, the exact run scripts `ws_*.sh`)
- **Visual explainer (Artifact):** private link in the task doc's `artifact:` field
- **GPU handoff (if you go that direction):** `project-hub/projects/FastGA/GPU-HANDOFF-from-workload-imbalance.md`

On this server the hub is at `/work/shared/users/phd/jl4257/Project/project-hub/` (world-readable).

---

## 5. Where you could take it next

- **Memory guard (the main loose end).** v3c holds all seeds resident (~27 GiB, scales with genome
  size × divergence). Add an **adaptive fallback** to the old streaming/per-part mode when host RAM is
  tight, so it's safe on smaller machines. This is the one thing blocking "ship it everywhere."
- **Land PR #2.** It's an open *draft*; a review pass + un-drafting would merge the 5.33× into the fork.
- **Bigger / more datasets.** Try other genome pairs (or self-vs-self) to confirm the win generalizes;
  the pool self-balances, so it should, but the memory footprint grows with total seed count.
- **GPU.** The measured verdict (see the GPU handoff) is that a faithful drop-in port of the aligner is
  *below 1×* against this now-optimized CPU baseline; the real opportunity, if any, is batching at the
  **wavefront** level (≈812K uniform-band waves) rather than the contig-pair level. Start from the GPU
  handoff doc, and re-baseline against `workstealing-v3c`, not stock.

**One-line orientation:** build `workstealing-v3c`, point it at the prebuilt `.gix` indices under
`/scratch/jl4257/seq_align/fastga_datasets/`, and check the sorted-PAF hash is `0ef422c2…` — that's a
working, verified 5.33× aligner you can extend.
