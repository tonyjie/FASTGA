#!/bin/bash
# Profile one FastGA genome pair (T=32) -> per-phase -L logs + /usr/bin/time.
# Usage: BL=<bin dir> ./run_pair.sh <label> <G1.fna> <G2.fna>
#   --dry-run : print the plan and exit without running FastGA.
set -u
DRY=0; [ "${1:-}" = "--dry-run" ] && { DRY=1; shift; }
LABEL="${1:?usage: run_pair.sh <label> <G1> <G2>}"
G1="${2:?need G1 fasta}"; G2="${3:?need G2 fasta}"
BL="${BL:?set BL=<baseline bin dir>}"
T=32; REPS="${REPS:-3}"; MIN_FREE_GB="${MIN_FREE_GB:-150}"
SCRATCH="${SCRATCH:-/scratch/jl4257/matrix_run}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${OUT:-$HERE/divergence}/$LABEL"

for f in "$G1" "$G2"; do [ -r "$f" ] || { echo "missing genome: $f" >&2; exit 2; }; done
[ -x "$BL/FastGA" ] || { echo "no FastGA in BL=$BL" >&2; exit 2; }

disk_guard() {  # refuse if free space on $1 < MIN_FREE_GB
  mkdir -p "$1"
  local avail; avail=$(df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc '0-9')
  [ -n "$avail" ] || { echo "disk guard: cannot stat $1" >&2; exit 3; }
  if [ "$avail" -lt "$MIN_FREE_GB" ]; then
    echo "disk guard: $1 has ${avail}GB free < ${MIN_FREE_GB}GB required — refusing" >&2
    exit 3
  fi
  echo "disk guard: ${avail}GB free on $1 (need ${MIN_FREE_GB}) — OK"
}

echo "=== pair '$LABEL': $(basename "$G1") x $(basename "$G2"), T=$T, $REPS reps ==="
disk_guard "$SCRATCH"
if [ "$DRY" = 1 ]; then
  echo "[dry-run] would run $REPS reps of: $BL/FastGA -v -T$T -P<tmp> -L:<log> g1 g2"
  echo "[dry-run] outputs -> $OUT/{logs,time}/"; exit 0
fi

mkdir -p "$OUT/logs" "$OUT/time"
# record measured genome sizes (bp) for the secondary divergence proxy
size_bp() { grep -v '^>' "$1" | tr -d '\n' | wc -c; }
{ echo -e "genome\tsize_bp"; echo -e "G1\t$(size_bp "$G1")"; echo -e "G2\t$(size_bp "$G2")"; } > "$OUT/meta.tsv"

for rep in $(seq 1 "$REPS"); do
  W="$SCRATCH/work"; TMP="$SCRATCH/tmp"; rm -rf "$W" "$TMP"; mkdir -p "$W" "$TMP"
  ln -sf "$G1" "$W/g1.fna"; ln -sf "$G2" "$W/g2.fna"
  ( cd "$W" && /usr/bin/time -v -o "$OUT/time/rep${rep}.time" \
      "$BL/FastGA" -v -T$T -P"$TMP" -L:"$OUT/logs/rep${rep}.Llog" g1.fna g2.fna \
      > /dev/null 2> "$OUT/logs/rep${rep}.stderr" )
  echo "[$(date +%H:%M:%S)] $LABEL rep$rep done"
  rm -rf "$W" "$TMP"
done
echo "PAIR_DONE $LABEL"
