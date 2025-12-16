#!/usr/bin/env bash
set -euo pipefail

# verify_scheme1_exact_match.sh
#
# Scheme-1 correctness gate for FastGA: ensure alignment output is identical to a baseline
# by comparing order-independent hashes of ALNtoPAF output (optionally with CIGAR).
#
# Requires: ALNtoPAF (built in FASTGA repo), coreutils, sort, sha256sum
# Optional: lots of temp space for sorting large PAF streams.

usage() {
  cat <<'USAGE'
Usage:
  verify_scheme1_exact_match.sh \
    --fastga-bin <FASTGA_repo_dir> \
    --baseline <baseline.1aln> \
    --candidate <candidate.1aln> \
    [--threads <T_for_ALNtoPAF>] \
    [--tmpdir <dir_for_sort_temp>] \
    [--with-cigar]

Meaning:
  - Compares baseline vs candidate using:
      1) ALNtoPAF (no CIGAR) : sha256 of sorted PAF lines
      2) (optional) ALNtoPAF -x (CIGAR X/=) : sha256 of sorted PAF lines

Exit codes:
  0  match
  1  mismatch
  2  invalid usage / missing tools

Example:
  ./verify_scheme1_exact_match.sh \
    --fastga-bin /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA \
    --baseline  /scratch/.../thread_8/GRCh38_vs_CHM13.1aln \
    --candidate /scratch/.../thread_32/GRCh38_vs_CHM13.1aln \
    --threads 32 \
    --tmpdir /scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/_tmp_sort \
    --with-cigar
USAGE
}

FASTGA_BIN=""
BASELINE=""
CANDIDATE=""
THREADS=32
TMPDIR_SORT=""
WITH_CIGAR=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastga-bin) FASTGA_BIN="$2"; shift 2;;
    --baseline) BASELINE="$2"; shift 2;;
    --candidate) CANDIDATE="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --tmpdir) TMPDIR_SORT="$2"; shift 2;;
    --with-cigar) WITH_CIGAR=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$FASTGA_BIN" || -z "$BASELINE" || -z "$CANDIDATE" ]]; then
  usage
  exit 2
fi

ALNtoPAF_BIN="$FASTGA_BIN/ALNtoPAF"
if [[ ! -x "$ALNtoPAF_BIN" ]]; then
  echo "ERROR: cannot execute ALNtoPAF at: $ALNtoPAF_BIN" >&2
  exit 2
fi

command -v sha256sum >/dev/null 2>&1 || { echo "ERROR: sha256sum not found" >&2; exit 2; }
command -v sort >/dev/null 2>&1 || { echo "ERROR: sort not found" >&2; exit 2; }

if [[ -z "$TMPDIR_SORT" ]]; then
  TMPDIR_SORT="$(mktemp -d)"
  trap 'rm -rf "$TMPDIR_SORT"' EXIT
else
  mkdir -p "$TMPDIR_SORT"
fi

hash_sorted_paf() {
  local aln="$1"
  local extra="$2"  # e.g. "" or "-x"

  # Order-independent hash: sort full lines (LC_ALL=C for determinism), then sha256.
  "$ALNtoPAF_BIN" $extra -T"$THREADS" "$aln" \
    | LC_ALL=C sort -T "$TMPDIR_SORT" \
    | sha256sum \
    | awk '{print $1}'
}

# 1) Plain PAF
base_plain=$(hash_sorted_paf "$BASELINE" "")
cand_plain=$(hash_sorted_paf "$CANDIDATE" "")

if [[ "$base_plain" != "$cand_plain" ]]; then
  echo "MISMATCH: sorted PAF (no CIGAR)" >&2
  echo "  baseline : $base_plain" >&2
  echo "  candidate: $cand_plain" >&2
  exit 1
fi

echo "OK: sorted PAF (no CIGAR) hash matches: $base_plain"

# 2) Optional CIGAR PAF
if [[ "$WITH_CIGAR" -eq 1 ]]; then
  base_x=$(hash_sorted_paf "$BASELINE" "-x")
  cand_x=$(hash_sorted_paf "$CANDIDATE" "-x")

  if [[ "$base_x" != "$cand_x" ]]; then
    echo "MISMATCH: sorted PAF with CIGAR (-x)" >&2
    echo "  baseline : $base_x" >&2
    echo "  candidate: $cand_x" >&2
    exit 1
  fi

  echo "OK: sorted PAF with CIGAR (-x) hash matches: $base_x"
fi

echo "PASS: Scheme-1 exact-match gate satisfied."
