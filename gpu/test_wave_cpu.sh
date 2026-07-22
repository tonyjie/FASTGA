#!/bin/bash
# test_wave_cpu.sh -- TDD acceptance test for gpu/wave_bench_cpu.c (Task 3: CPU baseline
# harness, genome-resident, real Local_Alignment).
#
# Asserts:
#   1. it processes the full seed distribution: 518037 seeds.
#   2. endpoint agreement vs the .1aln reference (ref_ab/ae/bb/be) within 50bp is >= 80%
#      (matches the deep-set 82% measured earlier; a regression below this means the
#      seed/orientation setup is wrong -- the signal this test guards).
#   3. a .cpuref file is written.
#
# Usage: gpu/test_wave_cpu.sh [in.1aln]
set -u

ALN="${1:-/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/output.1aln}"
ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
EXTRACT_SEEDS="$ROOT_DIR/extract_seeds"
BIN="$ROOT_DIR/gpu/wave_bench_cpu"
SEEDS="$ROOT_DIR/gpu/.test_wave.seeds"
CPUREF="$ROOT_DIR/gpu/.test_wave.cpuref"

fail() { echo "FAIL: $1"; exit 1; }

[ -x "$EXTRACT_SEEDS" ] || fail "extract_seeds binary not found: $EXTRACT_SEEDS (build with 'make extract_seeds')"
[ -x "$BIN" ]           || fail "wave_bench_cpu binary not found: $BIN (build with 'make wave_bench_cpu')"
[ -f "$ALN" ]           || fail "input .1aln not found: $ALN"

echo "Extracting seeds from $ALN ..."
"$EXTRACT_SEEDS" "$ALN" "$SEEDS" || fail "extract_seeds exited non-zero"
[ -f "$SEEDS" ] || fail "seeds file was not created: $SEEDS"

echo "Running CPU wave bench (this loads both genomes resident, then 1-core + 32-core passes)..."
OUT=$("$BIN" "$ALN" "$SEEDS" "$CPUREF" 32 2>&1) || fail "wave_bench_cpu exited non-zero:
$OUT"
echo "$OUT"

NSEEDS=$(echo "$OUT" | grep -oE '^loaded [0-9]+ seeds' | grep -oE '[0-9]+')
[ -n "$NSEEDS" ] || fail "could not parse seed count from output"
[ "$NSEEDS" -eq 518037 ] || fail "seed count mismatch: got $NSEEDS want 518037"

AGREE=$(echo "$OUT" | grep -oE 'within 50bp: [0-9]+/[0-9]+ \([0-9.]+%\)' | grep -oE '\([0-9.]+%\)' | tr -d '(%)')
[ -n "$AGREE" ] || fail "could not parse endpoint-agreement percentage from output"

PASS_THRESH=80
AWK_RESULT=$(awk -v a="$AGREE" -v t="$PASS_THRESH" 'BEGIN{ print (a>=t) ? "ok" : "no" }')
[ "$AWK_RESULT" = "ok" ] || fail "endpoint agreement $AGREE% < ${PASS_THRESH}% threshold (seed/orientation setup likely wrong)"

[ -f "$CPUREF" ] || fail ".cpuref file was not written: $CPUREF"

echo "PASS: nseeds=$NSEEDS, endpoint agreement=${AGREE}% (>= ${PASS_THRESH}% threshold), .cpuref written to $CPUREF"
rm -f "$SEEDS" "$CPUREF"
exit 0
