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

ROOT_DIR="/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/EXAMPLE"
HAP1_FASTA="${ROOT_DIR}/HAP1.fasta.gz"
HAP2_FASTA="${ROOT_DIR}/HAP2.fasta.gz"

# remove the log file if it exists
rm -f ${OUTPUT_DIR}/fastga.log

/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/FastGA \
    -vk \
    -v \
    -T32 \
    -L:${OUTPUT_DIR}/fastga.log \
    -P${OUTPUT_DIR} \
    -1:${OUTPUT_DIR}/H1vH2.1aln \
    ${HAP1_FASTA} \
    ${HAP2_FASTA} 
