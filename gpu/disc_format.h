/* Discovery-task format: for validating a GPU LOCAL-alignment kernel that must
 * rediscover the alignment endpoints from a seed (like FastGA's Local_Alignment),
 * not just compute edit distance on a pre-cut segment.
 *
 * For each FastGA alignment we dump a WINDOW = the aligned extent grown by `margin`
 * on all sides (clipped to contig), a seed anchor = the alignment midpoint, and the
 * reference endpoints (window-relative) + diffs.  A GPU local x-drop wave started at
 * the seed should rediscover ref_ab..ref_ae (A) and ref_bb..ref_be (B).
 *
 * NUMERIC sequences (bytes 0..3).  Little-endian.
 *
 *   Header (16 bytes): uint32 magic=0x46474144 ("DAGF"), ntasks, margin, reserved
 *   Per task:
 *     uint32 aw, bw                 (A/B window lengths)
 *     uint32 sa, sb                 (seed anchor, window-relative: 0<=sa<aw, 0<=sb<bw)
 *     int32  ref_ab, ref_ae         (reference A start/end, window-relative)
 *     int32  ref_bb, ref_be         (reference B start/end, window-relative)
 *     uint32 ref_diffs
 *     uint8  A[aw], uint8 B[bw]
 */
#ifndef FGA_DISC_FORMAT_H
#define FGA_DISC_FORMAT_H
#include <stdint.h>
#define FGA_DISC_MAGIC 0x46474144u
typedef struct { uint32_t magic, ntasks, margin, reserved; } FGADiscHeader;
#endif
