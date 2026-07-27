#!/bin/bash
WT=/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/.claude/worktrees/agent-opt
ALN=/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/output.1aln
SEEDS=$WT/gpu/wave.seeds
BIN=$WT/gpu/wave_bench_gpu
OUT=/home/jl4257/.claude/jobs/2237663e/tmp
export TMPDIR=$OUT
M=sm__warps_active.avg.pct_of_peak_sustained_active,\
launch__occupancy_limit_blocks,launch__occupancy_limit_registers,\
launch__occupancy_limit_shared_mem,launch__occupancy_limit_warps,\
launch__waves_per_multiprocessor,\
l1tex__t_sector_hit_rate.pct,lts__t_sector_hit_rate.pct,\
dram__bytes.sum.per_second,gpu__dram_throughput.avg.pct_of_peak_sustained_elapsed,\
sm__throughput.avg.pct_of_peak_sustained_elapsed,\
smsp__thread_inst_executed_per_inst_executed.ratio
cd $WT

echo "############ PROFILE 1: discovery forward_sweep_warp ############"
ncu --kill yes --launch-count 1 --launch-skip 4 --kernel-name forward_sweep_warp \
    --metrics $M $BIN $ALN $SEEDS --stage1-fit 2>&1 | grep -vE '^==WARNING==|^$'
echo; echo "############ PROFILE 2: trace wave_trace_warp @ TR_POOL=512 ############"
TR_POOL=512 ncu --kill yes --launch-count 1 --launch-skip 4 --kernel-name wave_trace_warp \
    --metrics $M $BIN $ALN $SEEDS --stage2-fit 2>&1 | grep -vE '^==WARNING==|^$'
echo; echo "############ PROFILE 3: trace wave_trace_warp @ TR_POOL=2160 ############"
TR_POOL=2160 ncu --kill yes --launch-count 1 --launch-skip 4 --kernel-name wave_trace_warp \
    --metrics $M $BIN $ALN $SEEDS --stage2-fit 2>&1 | grep -vE '^==WARNING==|^$'
echo "############ DONE ############"
