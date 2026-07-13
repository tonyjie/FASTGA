/* Trace-point task format: for validating a GPU kernel that emits FastGA trace-points
 * (not just endpoints/diffs).  For each FastGA alignment we dump the exact aligned
 * rectangle A[abpos..aepos] x B[bbpos..bepos] (NUMERIC bytes 0..3), the ORIGINAL contig
 * A-coordinate abpos (needed to phase trace-point panels on global multiples of tspace),
 * and FastGA's reference trace-point vector recomputed by Compute_Alignment(DIFF_TRACE)
 * on that same rectangle.
 *
 * The GPU aligns A[0..aw] vs B[0..bw] corner-to-corner (fixed endpoints) and must emit a
 * uint16 pair vector identical to ref_trace: per A-panel k = floor((abpos + a_local)/tspace),
 * even = diffs in panel (subs+ins+del), odd = Delta-b (B-bases consumed).  Contract:
 *   tlen = 2*(ceil(aepos/tspace) - floor(abpos/tspace)),
 *   sum(even) = ref_diffs,  sum(odd) = bw = bepos-bbpos.
 *
 * Little-endian.
 *   Header (16 bytes): uint32 magic=0x54474146 ("FAGT"), ntasks, tspace, reserved
 *   Per task:
 *     int32  abpos, aepos, bbpos, bepos     (original contig coords)
 *     uint32 aw, bw                         (aw=aepos-abpos, bw=bepos-bbpos)
 *     uint32 ref_diffs
 *     uint32 ref_tlen                       (# uint16 entries in ref_trace = 2*npanels)
 *     uint8  A[aw], uint8 B[bw]
 *     uint16 ref_trace[ref_tlen]
 */
#ifndef FGA_TRACE_FORMAT_H
#define FGA_TRACE_FORMAT_H
#include <stdint.h>
#define FGA_TRACE_MAGIC 0x54474146u
typedef struct { uint32_t magic, ntasks, tspace, reserved; } FGATraceHeader;
#endif
