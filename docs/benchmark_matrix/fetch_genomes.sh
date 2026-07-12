#!/bin/bash
# Download the divergence-axis query genomes (chimp, siamang, pig, mouse) to /scratch.
# CHM13 + GRCh38 are already local. Usage: ./fetch_genomes.sh [--dry-run]
set -u
DRY=0; [ "${1:-}" = "--dry-run" ] && DRY=1
DEST="${DEST:-/scratch/jl4257/seq_align/fastga_datasets}"
MIN_FREE_GB="${MIN_FREE_GB:-150}"

# name | approx download GB | URL  (NCBI datasets accessions; pig from GigaDB)
GENOMES=(
"chimpanzee|3.0|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.1_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz"
"siamang|2.9|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/878/055/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.fna.gz"
"mouse|2.6|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/188/535/GCA_964188535.1_C57BL_6J_T2T_v1/GCA_964188535.1_C57BL_6J_T2T_v1_genomic.fna.gz"
)
# NOTE: pig is on GigaDB (dataset 102692), not NCBI — its URL must be confirmed manually.
# Add once verified:  "pig|0.8|<gigadb direct fna.gz url>"

avail=$(df -BG --output=avail "$(dirname "$DEST")" 2>/dev/null | tail -1 | tr -dc '0-9')
echo "dest=$DEST  free=${avail:-?}GB  (need >= ${MIN_FREE_GB})"
[ -n "$avail" ] || { echo "disk guard: cannot determine free space for $(dirname "$DEST")" >&2; exit 3; }
[ "$avail" -lt "$MIN_FREE_GB" ] && { echo "disk guard: refusing" >&2; exit 3; }

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
