# FastGA bottleneck profiling matrix

Characterizes how FastGA's per-phase bottleneck shifts across datasets. **v1 = divergence
axis** (CHM13 × {GRCh38, chimp, siamang, pig, mouse}, ~3 Gbp fixed, T=32). Design spec:
`docs/superpowers/specs/2026-07-11-fastga-bottleneck-profiling-matrix-design.md`.
Dataset reference: `datasets_inventory.md`.

## Run it (all I/O on /scratch; disk guard at 150 GB free)

```bash
export BL=/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA-baseline   # baseline bin
./fetch_genomes.sh                 # download chimp, siamang, mouse (pig: add GigaDB URL first)
BL=$BL ./run_divergence_axis.sh    # profiles each pair (CHM13 as genome-1)
column -t divergence/results.tsv
open divergence/divergence_phase_share.png
```

All 5 points use **CHM13 as genome-1 (reference)**: CHM13 x {GRCh38, chimp, siamang, pig, mouse}.
The human point is run as CHM13 x GRCh38 like every other point — no reused/reversed-orientation
logs. `results.tsv` and the figure are produced by actually running `run_divergence_axis.sh`;
nothing is committed as a sample, so run it before reading `divergence/`.

`rss_mb` in `results.tsv` is the **true peak RSS** (MB), read from `/usr/bin/time -v` output
(`Maximum resident set size`) captured per rep in `divergence/<label>/time/`, not FastGA's own
`-L` log (which underreports peak memory).

Per-pair peak ≈ 80 GB transient (built, then cleaned each rep); persistent downloads ≈ 10–25 GB.
`pig` is on GigaDB (dataset 102692), not NCBI — confirm its `.fna.gz` URL and add the row to
`fetch_genomes.sh` before running the pig point.

## Hypothesis

`sort+align` share falls from ~81% (human, most similar) toward the divergent end as less
sequence is alignable; GIX/seed shares rise. See the stacked figure.

## Results

_(fill in after the runs: paste the `results.tsv` table and the figure.)_
