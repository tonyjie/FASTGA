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
