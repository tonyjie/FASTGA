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
