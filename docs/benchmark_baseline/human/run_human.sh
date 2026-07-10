#!/bin/bash
# Human-genome profiling for STOCK UPSTREAM FastGA, T=32 (GRCh38 x CHM13, ~3.1 Gbp each).
#   - performance: 3 reps on a fast local disk; captures FastGA -L per-stage logs
#     (wall + CPU% per phase) and /usr/bin/time.  -> perf_data/logs/, perf_data/time/
#   - storage: 1 rep on an NFS mount (du must see open()-then-unlink()ed temp files;
#     NFS keeps them via silly-rename, local ext4/xfs do not); monitors work+tmp over
#     time.  -> storage_data/timeline.tsv, storage_data/storage.Llog
#
# Baseline binary: ddeea32 FastGA/GIXmake + a 5671357 FAtoGDB (ddeea32's FAtoGDB
# segfaults building CHM13's GDB — ANO regression). Pass BL=<dir with those>.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BL="${BL:?set BL=<baseline bin dir>}"
export PATH="$BL:$PATH"
G1="${G1:-/scratch/jl4257/seq_align/fastga_datasets/GRCh38/GCF_000001405.40_GRCh38.p14_genomic.fna}"
G2="${G2:-/scratch/jl4257/seq_align/fastga_datasets/CHM13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna}"
FAST="${FAST_SCRATCH:-/scratch/jl4257/human_perf}"      # local fast disk for perf
STORE="${STORE_SCRATCH:-/scratch/jl4257/human_store}"  # local disk for storage (see temp note)
T=32

PD="$HERE/perf_data"; SD="$HERE/storage_data"
mkdir -p "$PD/logs" "$PD/time" "$SD"

setup() { d=$1; rm -rf "$d"; mkdir -p "$d/work" "$d/tmp"; ln -sf "$G1" "$d/work/g1.fna"; ln -sf "$G2" "$d/work/g2.fna"; }

echo "=== PERFORMANCE: 3x T=$T on $FAST ==="
for rep in $( [ "${STORAGE_ONLY:-0}" = 1 ] && echo || echo 1 2 3 ); do
  setup "$FAST"; cd "$FAST/work"
  echo "[$(date +%H:%M:%S)] perf rep$rep start"
  /usr/bin/time -v -o "$PD/time/rep${rep}.time" \
    "$BL/FastGA" -v -T$T -P"$FAST/tmp" -L:"$PD/logs/rep${rep}.Llog" g1.fna g2.fna > /dev/null 2> "$PD/logs/rep${rep}.stderr"
  echo "[$(date +%H:%M:%S)] perf rep$rep done"
  cd "$HERE"; rm -rf "$FAST"
done

echo "=== STORAGE: 1x T=$T on $STORE ==="
# persistent = du(work) [real GDB/GIX files, correct on any FS];
# temp is open()-then-unlink()ed so `du` on local disk misses it — we sum the
# process's open fds pointing into tmp, deduped by inode, using st_SIZE (bytes
# written). st_size matches du-on-NFS; st_blocks would over-count by ~2x because
# XFS speculatively preallocates blocks for the actively-written seed files.
setup "$STORE"; cd "$STORE/work"
MON="$SD/timeline.tsv"; echo -e "elapsed_s\tpersistent_mb\ttemp_du_mb\ttemp_fd_mb" > "$MON"
echo "[$(date +%H:%M:%S)] storage run start"
"$BL/FastGA" -v -T$T -P"$STORE/tmp" -L:"$SD/storage.Llog" g1.fna g2.fna > /dev/null 2> "$SD/storage.stderr" &
FPID=$!; START=$(date +%s.%N)
while kill -0 "$FPID" 2>/dev/null; do
  W=$(du -sb "$STORE/work" 2>/dev/null|cut -f1)
  Tdu=$(du -sb "$STORE/tmp" 2>/dev/null|cut -f1)
  Tfd=$(for p in "$FPID" $(pgrep -P "$FPID" 2>/dev/null||true); do
          [ -d /proc/$p/fd ] || continue
          for fd in /proc/$p/fd/*; do [ -e "$fd" ] || continue
            case "$(readlink "$fd" 2>/dev/null||true)" in "$STORE/tmp"*) stat -L -c '%i %s' "$fd" 2>/dev/null;; esac
          done
        done | sort -u -k1,1 | awk '{s+=$2}END{print s+0}')
  echo -e "$(echo "$(date +%s.%N)-$START"|bc)\t$(echo "scale=1;${W:-0}/1048576"|bc)\t$(echo "scale=1;${Tdu:-0}/1048576"|bc)\t$(echo "scale=1;${Tfd:-0}/1048576"|bc)" >> "$MON"
  sleep 2
done
wait "$FPID"
echo "[$(date +%H:%M:%S)] storage run done"
cd "$HERE"; rm -rf "$STORE"
echo "ALL_DONE"
