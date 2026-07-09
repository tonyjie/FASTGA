#!/bin/bash
# Thread scaling study for upstream FastGA (commit 10ebff7 = upstream ddeea32 + .gitignore)
# on the repo EXAMPLE dataset (HAP1 vs HAP2, ~86 Mbp each).
set -u
REPO=/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA
FGA=$REPO/FastGA
RUN=/tmp/claude-1564080/-work-shared-users-phd-jl4257-Project-genomics-agent-FASTGA/66f02e18-e698-4955-9f60-687e7a3925eb/scratchpad/tscaling
LOGS=$RUN/logs
mkdir -p $LOGS
cd $RUN || exit 1

# fresh copies of inputs
cp -f $REPO/EXAMPLE/HAP1.fasta.gz $REPO/EXAMPLE/HAP2.fasta.gz $RUN/

clean_generated() {
  rm -f $RUN/HAP1.1gdb $RUN/HAP2.1gdb $RUN/.HAP1.bps $RUN/.HAP2.bps \
        $RUN/HAP1.gix $RUN/HAP2.gix $RUN/.HAP1.ktab.* $RUN/.HAP2.ktab.* \
        $RUN/.HAP1.post.* $RUN/.HAP2.post.* $RUN/HAP1.1ano $RUN/HAP2.1ano \
        $RUN/out.paf 2>/dev/null
}

RESULTS=$RUN/results.tsv
echo -e "threads\trep\twall_s\tuser_s\tsys_s\tcpu_pct\tmaxrss_kb\tnonredundant_aln" > $RESULTS

THREADS="1 2 4 8 16 32 64"
REPS="1 2 3"

echo "START $(date)"
for T in $THREADS; do
  for R in $REPS; do
    clean_generated
    tag=$(printf "T%02d_rep%d" $T $R)
    tf=$LOGS/${tag}.time
    fl=$LOGS/${tag}.fastga_log
    /usr/bin/time -v -o $tf \
      $FGA -v -T$T -P. -L:$fl HAP1.fasta.gz HAP2.fasta.gz > $RUN/out.paf 2> $LOGS/${tag}.stderr
    # parse
    wall=$(grep 'Elapsed' $tf | sed 's/.*: //')
    # convert m:ss or h:mm:ss to seconds
    wall_s=$(echo "$wall" | awk -F: '{ if(NF==3){print $1*3600+$2*60+$3} else if(NF==2){print $1*60+$2} else {print $1} }')
    user_s=$(grep 'User time' $tf | sed 's/.*: //')
    sys_s=$(grep 'System time' $tf | sed 's/.*: //')
    cpu=$(grep 'Percent of CPU' $tf | sed 's/.*: //; s/%//')
    rss=$(grep 'Maximum resident' $tf | sed 's/.*: //')
    nr=$(grep -oE 'non-redundant aln[^,]*' $LOGS/${tag}.stderr 2>/dev/null | grep -oE '[0-9,]+ non' | grep -oE '[0-9,]+' | head -1)
    [ -z "$nr" ] && nr=$(grep -oE '[0-9]+ non-redundant' $fl 2>/dev/null | grep -oE '[0-9]+' | head -1)
    echo -e "${T}\t${R}\t${wall_s}\t${user_s}\t${sys_s}\t${cpu}\t${rss}\t${nr}" >> $RESULTS
    echo "  done $tag  wall=${wall_s}s cpu=${cpu}% nr_aln=${nr}"
  done
done
clean_generated
echo "DONE $(date)"
echo "=== results.tsv ==="
cat $RESULTS
