#!/bin/bash
# Thread-scaling sweep for STOCK UPSTREAM FastGA on the EXAMPLE dataset
# (HAP1 vs HAP2, ~86 Mbp each). Writes results.tsv next to this script.
#
# The baseline binary must be a stock-upstream build (NOT optimize-memory /
# agent-optimization). Point FASTGA at it — see ../README.md "Reproduce":
#   FASTGA=/path/to/upstream/FastGA bash run_thread_scaling.sh
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"            # -> repo root (…/FASTGA)
FGA="${FASTGA:-$REPO/FastGA}"                    # override with FASTGA=…
export PATH="$(dirname "$FGA"):$PATH"           # FastGA shells out to GIXmake
RUN="$(mktemp -d)"; LOGS="$RUN/logs"; mkdir -p "$LOGS"
RESULTS="$HERE/results.tsv"

echo "Binary : $FGA"
echo "Workdir: $RUN"
cp -f "$REPO/EXAMPLE/HAP1.fasta.gz" "$REPO/EXAMPLE/HAP2.fasta.gz" "$RUN/"
cd "$RUN" || exit 1

clean_generated() {
  rm -f HAP1.1gdb HAP2.1gdb .HAP1.bps .HAP2.bps HAP1.gix HAP2.gix \
        .HAP1.ktab.* .HAP2.ktab.* .HAP1.post.* .HAP2.post.* \
        HAP1.1ano HAP2.1ano out.paf 2>/dev/null
}

echo -e "threads\trep\twall_s\tuser_s\tsys_s\tcpu_pct\tmaxrss_kb\tnonredundant_aln" > "$RESULTS"
THREADS="${THREADS:-1 2 4 8 16 32}"   # T=64 fails: FastGA hard-caps at 32
REPS="${REPS:-1 2 3}"

echo "START $(date)"
for T in $THREADS; do
  for R in $REPS; do
    clean_generated
    tag=$(printf "T%02d_rep%d" "$T" "$R"); tf="$LOGS/${tag}.time"; fl="$LOGS/${tag}.fastga_log"
    /usr/bin/time -v -o "$tf" \
      "$FGA" -v -T"$T" -P. -L:"$fl" HAP1.fasta.gz HAP2.fasta.gz > out.paf 2> "$LOGS/${tag}.stderr"
    wall=$(grep 'Elapsed' "$tf" | sed 's/.*: //')
    wall_s=$(echo "$wall" | awk -F: '{if(NF==3){print $1*3600+$2*60+$3}else if(NF==2){print $1*60+$2}else{print $1}}')
    user_s=$(grep 'User time' "$tf" | sed 's/.*: //')
    sys_s=$(grep 'System time' "$tf" | sed 's/.*: //')
    cpu=$(grep 'Percent of CPU' "$tf" | sed 's/.*: //; s/%//')
    rss=$(grep 'Maximum resident' "$tf" | sed 's/.*: //')
    nr=$(grep -oE '[0-9]+ non-redundant' "$fl" 2>/dev/null | grep -oE '[0-9]+' | head -1)
    echo -e "${T}\t${R}\t${wall_s}\t${user_s}\t${sys_s}\t${cpu}\t${rss}\t${nr}" >> "$RESULTS"
    echo "  done $tag  wall=${wall_s}s cpu=${cpu}% nr_aln=${nr}"
  done
done
clean_generated; rm -rf "$RUN"
echo "DONE $(date) -> $RESULTS"
