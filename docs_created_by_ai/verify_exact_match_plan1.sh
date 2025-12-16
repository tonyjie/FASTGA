#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
verify_exact_match_plan1.sh

Plan 1 回归验证：比较两个 FastGA 的 .1aln 输出是否“结果等价”。

核心判据（推荐）：
  ALNtoPAF -x | sort | sha256sum  必须一致

为什么要 sort？
  多线程/不同 run 可能改变输出顺序，但 alignment 集合不变。
  排序后 hash 是 order-independent 的强判据。

用法：
  verify_exact_match_plan1.sh <baseline.1aln> <candidate.1aln> [--threads N] [--tmpdir DIR]

参数：
  --threads N   用于 ALNtoPAF 的线程数（默认 32）
  --tmpdir DIR  sort 的临时目录（默认 /tmp）

返回码：
  0  通过（hash 相同）
  2  失败（hash 不同）

USAGE
}

if [[ $# -lt 2 ]]; then
  usage
  exit 1
fi

BASELINE="$1"
CANDIDATE="$2"
shift 2

THREADS=32
TMPDIR_SORT="/tmp"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads)
      THREADS="$2"; shift 2;;
    --tmpdir)
      TMPDIR_SORT="$2"; shift 2;;
    -h|--help)
      usage; exit 0;;
    *)
      echo "Unknown arg: $1" >&2
      usage
      exit 1;;
  esac
done

if [[ ! -f "$BASELINE" ]]; then
  echo "Baseline not found: $BASELINE" >&2
  exit 1
fi
if [[ ! -f "$CANDIDATE" ]]; then
  echo "Candidate not found: $CANDIDATE" >&2
  exit 1
fi

# Resolve tool paths relative to this script location (expects to live under FASTGA/docs_created_by_ai/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FASTGA_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ALNtoPAF_BIN="$FASTGA_DIR/ALNtoPAF"
ALNshow_BIN="$FASTGA_DIR/ALNshow"

if [[ ! -x "$ALNtoPAF_BIN" ]]; then
  echo "ALNtoPAF not executable at: $ALNtoPAF_BIN" >&2
  echo "(Hint) build FastGA tools first (make) or adjust script paths." >&2
  exit 1
fi
if [[ ! -x "$ALNshow_BIN" ]]; then
  echo "ALNshow not executable at: $ALNshow_BIN" >&2
  exit 1
fi

mkdir -p "$TMPDIR_SORT" || true

echo "=== Quick sanity: record counts (ALNshow) ==="
BASE_HDR="$($ALNshow_BIN "$BASELINE" 2>/dev/null | head -1 || true)"
CAND_HDR="$($ALNshow_BIN "$CANDIDATE" 2>/dev/null | head -1 || true)"
echo "baseline:  $BASE_HDR"
echo "candidate: $CAND_HDR"
echo

hash_sorted_paf_x() {
  local aln="$1"
  "$ALNtoPAF_BIN" -x -T"$THREADS" "$aln" \
    | LC_ALL=C sort -T "$TMPDIR_SORT" \
    | sha256sum \
    | awk '{print $1}'
}

hash_sorted_paf() {
  local aln="$1"
  "$ALNtoPAF_BIN" -T"$THREADS" "$aln" \
    | LC_ALL=C sort -T "$TMPDIR_SORT" \
    | sha256sum \
    | awk '{print $1}'
}

echo "=== Strong gate: sorted PAF(+CIGAR X/=) sha256 ==="
BASE_X_HASH="$(hash_sorted_paf_x "$BASELINE")"
CAND_X_HASH="$(hash_sorted_paf_x "$CANDIDATE")"
echo "baseline(sorted PAF -x):  $BASE_X_HASH"
echo "candidate(sorted PAF -x): $CAND_X_HASH"

# Optional: also report sorted PAF without cigar (cheaper)
echo
echo "=== Info: sorted PAF (no CIGAR) sha256 ==="
BASE_HASH="$(hash_sorted_paf "$BASELINE")"
CAND_HASH="$(hash_sorted_paf "$CANDIDATE")"
echo "baseline(sorted PAF):  $BASE_HASH"
echo "candidate(sorted PAF): $CAND_HASH"

if [[ "$BASE_X_HASH" == "$CAND_X_HASH" ]]; then
  echo
echo "PASS: exact-match (order-independent) under ALNtoPAF -x | sort" 
  exit 0
else
  echo >&2
  echo "FAIL: mismatch under ALNtoPAF -x | sort" >&2
  exit 2
fi
