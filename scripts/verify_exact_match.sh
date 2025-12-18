#!/usr/bin/env bash
set -euo pipefail

#######################################
# Alignment Output Correctness Verification
# Compares two .1aln files by hashing their sorted PAF outputs
#######################################

# Default values
ROOT_DIR="/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA"
THREADS=32
WITH_CIGAR=false

print_usage() {
    cat <<EOF
Usage: $0 --baseline <file> --candidate <file> [options]

Required:
  --baseline <file>      Baseline .1aln file
  --candidate <file>     Candidate .1aln file to compare

Optional:
  --threads <N>          Threads for ALNtoPAF (default: 32)
  --tmpdir <dir>         Temp directory for sorting
  --with-cigar           Also compare with CIGAR strings (-x flag)
  -h, --help             Show this help

Exit codes: 0=match, 1=mismatch, 2=error
EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --baseline)   BASELINE="$2"; shift 2 ;;
        --candidate)  CANDIDATE="$2"; shift 2 ;;
        --threads)    THREADS="$2"; shift 2 ;;
        --tmpdir)     TMPDIR_SORT="$2"; shift 2 ;;
        --with-cigar) WITH_CIGAR=true; shift ;;
        -h|--help)    print_usage; exit 0 ;;
        *)            echo "Unknown option: $1"; print_usage; exit 2 ;;
    esac
done

# Validate required arguments
if [[ -z "${BASELINE:-}" || -z "${CANDIDATE:-}" ]]; then
    echo "Error: Missing required arguments"
    print_usage
    exit 2
fi

# Validate ALNtoPAF exists
ALNtoPAF="$ROOT_DIR/ALNtoPAF"
if [[ ! -x "$ALNtoPAF" ]]; then
    echo "Error: ALNtoPAF not found at $ALNtoPAF"
    exit 2
fi

# Setup temp directory
if [[ -z "${TMPDIR_SORT:-}" ]]; then
    TMPDIR_SORT="$(mktemp -d)"
    trap 'rm -rf "$TMPDIR_SORT"' EXIT
else
    mkdir -p "$TMPDIR_SORT"
fi

# Compute hash of sorted PAF output
# Args: $1 = aln file, $2 = extra flags (e.g., "-x" for CIGAR)
compute_hash() {
    "$ALNtoPAF" $2 -T"$THREADS" "$1" \
        | LC_ALL=C sort -T "$TMPDIR_SORT" \
        | sha256sum \
        | cut -d' ' -f1
}

# Compare hashes and report result
# Args: $1 = description, $2 = baseline hash, $3 = candidate hash
compare_hashes() {
    local desc="$1" base_hash="$2" cand_hash="$3"
    
    if [[ "$base_hash" != "$cand_hash" ]]; then
        echo "MISMATCH: $desc"
        echo "  baseline:  $base_hash"
        echo "  candidate: $cand_hash"
        exit 1
    fi
    echo "OK: $desc - $base_hash"
}

# Run comparisons
echo "Comparing: $(basename "$BASELINE") vs $(basename "$CANDIDATE")"

compare_hashes "PAF (no CIGAR)" \
    "$(compute_hash "$BASELINE" "")" \
    "$(compute_hash "$CANDIDATE" "")"

if $WITH_CIGAR; then
    compare_hashes "PAF with CIGAR (-x)" \
        "$(compute_hash "$BASELINE" "-x")" \
        "$(compute_hash "$CANDIDATE" "-x")"
fi

echo "PASS: All checks passed"
