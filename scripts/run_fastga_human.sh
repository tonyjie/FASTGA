#!/bin/bash
## Usage: bash run_fastga.sh <output_directory>

# Check if output directory is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <output_directory>"
    echo "  <output_directory>: Directory where output files will be saved"
    exit 1
fi

OUTPUT_DIR="$1"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "Running FastGA with 32 threads"
echo "Output directory: ${OUTPUT_DIR}"

ROOT_DIR="/scratch/jl4257/seq_align/fastga_datasets"
GENOME1_FASTA="${ROOT_DIR}/GRCh38/GCF_000001405.40_GRCh38.p14_genomic.fna"
GENOME2_FASTA="${ROOT_DIR}/CHM13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"

# remove the log file if it exists
rm -f ${OUTPUT_DIR}/fastga_human.log

/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/FastGA \
    -vk \
    -v \
    -T32 \
    -L:${OUTPUT_DIR}/fastga_human.log \
    -P${OUTPUT_DIR} \
    -1:${OUTPUT_DIR}/GRCh38_vs_CHM13.1aln \
    ${GENOME1_FASTA} \
    ${GENOME2_FASTA} 
