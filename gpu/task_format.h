/* Shared binary task-file format for the FastGA GPU-alignment experiment.
 *
 * A "task" is one FastGA alignment's two aligned segments (as extracted the same
 * way ALNtoPAF does for its CIGAR: A = a-contig[abpos..aepos]; B = b-contig piece
 * [bbpos..bepos], reverse-complemented if the overlap is COMP), plus FastGA's
 * reference edit distance (path->diffs).  The GPU kernel recomputes the banded
 * unit-cost edit distance of (A,B) and is validated against `diffs`.
 *
 * Sequences are NUMERIC: one byte per base, values 0..3 (a,c,g,t).  Little-endian.
 *
 *   Header (16 bytes):
 *     uint32 magic   = 0x46474154   ("FGAT")
 *     uint32 ntasks
 *     uint32 wmax                    (segments with max(m,n) > wmax were skipped)
 *     uint32 reserved = 0
 *   Then ntasks records, each:
 *     uint32 m        (len of A segment = aepos-abpos)
 *     uint32 n        (len of B segment = bepos-bbpos)
 *     uint32 diffs    (FastGA path->diffs — the reference edit distance)
 *     uint8  A[m]     (0..3)
 *     uint8  B[n]     (0..3)
 */
#ifndef FGA_TASK_FORMAT_H
#define FGA_TASK_FORMAT_H
#include <stdint.h>
#define FGA_TASK_MAGIC 0x46474154u
typedef struct { uint32_t magic, ntasks, wmax, reserved; } FGATaskHeader;
#endif
