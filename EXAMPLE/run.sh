#!/bin/bash
# Usage: ./run.sh <num_thread>
num_thread=$1
echo "Running FastGA with $num_thread threads"

ROOT_DIR="/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/EXAMPLE"
THREAD_DIR="${ROOT_DIR}/thread_${num_thread}"

# remove all the intermediate files
# Note: GIX/GDB “hidden parts” are dotfiles (e.g. `.HAP1.bps`, `.HAP1.ktab.1`, `.HAP1.post.1`),
# and `*` does not match dotfiles by default. Use `.??*` to avoid matching `.` / `..`.
rm -f ${ROOT_DIR}/*.1aln
rm -f ${ROOT_DIR}/*.1gdb
rm -f ${ROOT_DIR}/*.gix
rm -f ${ROOT_DIR}/*.bps ${ROOT_DIR}/.??*.bps
rm -f ${ROOT_DIR}/*.ktab.* ${ROOT_DIR}/.??*.ktab.*
rm -f ${ROOT_DIR}/*.post.* ${ROOT_DIR}/.??*.post.*

mkdir -p ${THREAD_DIR}
# run FastGA
/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/FastGA \
  -vk \
  -v \
  -T${num_thread} \
  -L:${THREAD_DIR}/fastga.log \
  -P${THREAD_DIR} \
  -1:${THREAD_DIR}/H1vH2.1aln \
  ${ROOT_DIR}/HAP1.fasta.gz \
  ${ROOT_DIR}/HAP2.fasta.gz \
  2>&1 | tee ${THREAD_DIR}/fastga_stdout.log