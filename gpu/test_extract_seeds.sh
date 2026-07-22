#!/bin/bash
# test_extract_seeds.sh -- TDD test for gpu/extract_seeds.c (Task 2: full-distribution
# seed extractor).  Asserts:
#   1. header.nseeds == 518037
#   2. header.tspace == 100
#   3. record 0's ref_ab/ref_ae match ALNtoPAF row 1's A-start/A-end (cols 3/4, 0-based)
#      via the KNOWN contig sbeg offset (10000), checked independently at both ends --
#      not merely that the two ends agree with each other (which would pass even for
#      a constant wrong offset applied to both).
#   4. record 0's flags (COMP bit), blen, and ref_diffs cross-check against PAF row 1's
#      strand column, target-length column, and df:i: tag.
#
# Usage: gpu/test_extract_seeds.sh [in.1aln]
set -u

ALN="${1:-/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/output.1aln}"
ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
EXTRACTOR="$ROOT_DIR/extract_seeds"
ALNTOPAF="$ROOT_DIR/ALNtoPAF"
OUT="$ROOT_DIR/gpu/.test_output.seeds"

fail() { echo "FAIL: $1"; exit 1; }

[ -x "$EXTRACTOR" ] || fail "extractor binary not found: $EXTRACTOR (build with 'make extract_seeds')"
[ -x "$ALNTOPAF" ]  || fail "ALNtoPAF binary not found: $ALNTOPAF"
[ -f "$ALN" ]        || fail "input .1aln not found: $ALN"

echo "Running extractor on $ALN ..."
"$EXTRACTOR" "$ALN" "$OUT" || fail "extractor exited non-zero"

[ -f "$OUT" ] || fail "output .seeds file was not created: $OUT"

# --- Parse header (WaveSeedHeader: uint32 magic, nseeds, tspace, reserved = 16 bytes) ---
read -r MAGIC NSEEDS TSPACE RESERVED < <(python3 - "$OUT" <<'EOF'
import struct, sys
with open(sys.argv[1], "rb") as f:
    magic, nseeds, tspace, reserved = struct.unpack("<IIII", f.read(16))
print(magic, nseeds, tspace, reserved)
EOF
)

echo "header: magic=0x$(printf '%x' "$MAGIC") nseeds=$NSEEDS tspace=$TSPACE reserved=$RESERVED"

EXPECT_MAGIC=1398161751   # 0x53564157 "WAVS"
[ "$MAGIC" -eq "$EXPECT_MAGIC" ] || fail "magic mismatch: got $MAGIC want $EXPECT_MAGIC"
[ "$NSEEDS" -eq 518037 ]  || fail "nseeds mismatch: got $NSEEDS want 518037"
[ "$TSPACE" -eq 100 ]     || fail "tspace mismatch: got $TSPACE want 100"

# --- Row-1 check: ALNtoPAF cols 3/4 (0-based) are A-start/A-end.
#
# NOTE: ALNtoPAF reports coordinates *scaffold*-relative (it adds
# contigs1[acontig].sbeg, the contig's start offset within its scaffold --
# see ALNtoPAF.c's `aoff = contigs1[acontig].sbeg` / `aoff + path->abpos`).
# SeedRec's ref_ab/ref_ae are the raw `path->abpos/aepos`, i.e. *contig*-relative
# (per the brief: "ref_* from path"). For GRCh38 row 1 specifically, the first
# aligned contig of chr1 starts after a 10,000bp leading telomeric N-gap, so
# aoff=10000 != 0 here -- a direct equality would fail even though the extractor
# is correct.
#
# A weaker check (PAF_ab - ref_ab == PAF_ae - ref_ae) only proves the offset is
# *consistent* between start and end -- it cancels out and would still pass if
# the extractor applied any constant wrong offset (e.g. 0, or some other contig's
# sbeg) to both ref_ab and ref_ae. To actually catch that bug class we assert the
# offset equals the KNOWN value for this dataset/row: 10000, the sbeg of row 1's
# target contig (GRCh38 chr1's leading telomeric N-gap), checked independently at
# both ends.
PAF_ROW1=$("$ALNTOPAF" "$ALN" | head -1)
[ -n "$PAF_ROW1" ] || fail "ALNtoPAF produced no output"
# PAF cols (0-based): 2=A-start 3=A-end 4=strand 5=target-name 6=target-len
# (awk fields are 1-based, so $3/$4/$5/$6/$7 below); df:i: is a trailing tag.
PAF_AB=$(echo "$PAF_ROW1" | awk '{print $3}')
PAF_AE=$(echo "$PAF_ROW1" | awk '{print $4}')
PAF_STRAND=$(echo "$PAF_ROW1" | awk '{print $5}')
PAF_TLEN=$(echo "$PAF_ROW1" | awk '{print $7}')
PAF_DF=$(echo "$PAF_ROW1" | grep -oE 'df:i:[0-9]+' | cut -d: -f3)
echo "ALNtoPAF row1: A-start=$PAF_AB A-end=$PAF_AE strand=$PAF_STRAND target-len=$PAF_TLEN df:i:=$PAF_DF"
[ -n "$PAF_DF" ] || fail "could not parse df:i: tag from PAF row 1"

# --- Parse record 0 (SeedRec: 12 x int32 = 48 bytes) ---
read -r AREAD BREAD FLAGS ALEN BLEN SANTI SDIAG REF_AB REF_AE REF_BB REF_BE REF_DIFFS < <(python3 - "$OUT" <<'EOF'
import struct, sys
with open(sys.argv[1], "rb") as f:
    f.seek(16)  # skip header
    fields = struct.unpack("<12i", f.read(48))
print(*fields)
EOF
)

echo "record0: aread=$AREAD bread=$BREAD flags=$FLAGS alen=$ALEN blen=$BLEN " \
     "seed_anti=$SANTI seed_diag=$SDIAG ref_ab=$REF_AB ref_ae=$REF_AE ref_bb=$REF_BB ref_be=$REF_BE ref_diffs=$REF_DIFFS"

KNOWN_SBEG=10000   # GRCh38 chr1 leading telomeric N-gap; row 1's target contig sbeg
OFF_B=$((PAF_AB - REF_AB))
OFF_E=$((PAF_AE - REF_AE))
echo "contig->scaffold offset implied: from ab=$OFF_B, from ae=$OFF_E"
[ "$OFF_B" -eq "$KNOWN_SBEG" ] || fail "A-start offset wrong: PAF_ab-ref_ab=$OFF_B want $KNOWN_SBEG"
[ "$OFF_E" -eq "$KNOWN_SBEG" ] || fail "A-end offset wrong: PAF_ae-ref_ae=$OFF_E want $KNOWN_SBEG"

# --- Cross-check record 0 against the SAME PAF row 1 on fields with no
# coordinate-offset ambiguity: strand/COMP flag, target length (blen), and
# diff count (df:i:). These do not involve sbeg, so they independently
# corroborate that record 0 really corresponds to PAF row 1.
if [ "$PAF_STRAND" = "-" ]; then
  EXPECT_FLAGS_COMP=1
else
  EXPECT_FLAGS_COMP=0
fi
ACTUAL_FLAGS_COMP=$((FLAGS & 1))
[ "$ACTUAL_FLAGS_COMP" -eq "$EXPECT_FLAGS_COMP" ] || \
  fail "COMP flag mismatch: record0 flags=$FLAGS (comp bit=$ACTUAL_FLAGS_COMP) vs PAF strand=$PAF_STRAND"

[ "$BLEN" -eq "$PAF_TLEN" ] || fail "blen mismatch: record0 blen=$BLEN vs PAF target-len=$PAF_TLEN"

[ "$REF_DIFFS" -eq "$PAF_DF" ] || fail "ref_diffs mismatch: record0 ref_diffs=$REF_DIFFS vs PAF df:i:=$PAF_DF"

echo "PASS: nseeds=$NSEEDS tspace=$TSPACE row-1 ref_ab/ref_ae match ALNtoPAF ($PAF_AB/$PAF_AE)" \
     "via known offset $KNOWN_SBEG; flags/blen/ref_diffs cross-checked against PAF (strand=$PAF_STRAND blen=$PAF_TLEN df=$PAF_DF)"
rm -f "$OUT"
exit 0
