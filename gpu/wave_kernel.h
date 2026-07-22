/* wave_kernel.h -- C-callable GPU genome residency for the wave-parallelism characterization
 * study (Task 4). See docs/superpowers/plans/2026-07-22-gpu-wave-parallelism-characterization-plan.md
 * Task 4 and gpu/WAVE_PORT_NOTES.md for the sentinel/orientation contract this mirrors.
 *
 * Residency layout: BOTH genomes are uploaded as ONE flat resident NUMERIC (0-3) byte array
 * per genome side, with a per-contig base-offset table (byte offset of contig i's base 0
 * within the flat array). Contigs are laid out back-to-back separated by a single sentinel
 * byte (value 4) so that contig i's data spans flat[base[i]-1 .. base[i]+len_i] inclusive:
 *   flat[base[i]-1]        == 4                          (leading sentinel)
 *   flat[base[i] .. +len_i-1]  == contig i's NUMERIC bases
 *   flat[base[i]+len_i]    == 4                          (trailing sentinel, == next
 *                                                          contig's own leading sentinel)
 * This is exactly the align.c contract (New_Contig_Buffer/Get_Contig: seq[-1]==4,
 * seq[len]==4) applied across a whole concatenated genome, so a future wave kernel that
 * walks off a contig's end still hits a real `4` byte instead of another contig's bases.
 *
 * B needs BOTH orientations, mirroring wave_bench_cpu.c's CPU baseline: Bfwd (forward) and
 * Brev (reverse-complement, built once per contig via Complement_Seq) are uploaded as two
 * separate resident arrays sharing the SAME base-offset table (complementing in place does
 * not change contig boundaries/lengths). Callers pick Brev for COMP seeds, Bfwd otherwise --
 * exactly as wave_bench_cpu.c's `COMP(s->flags) ? Brev[...] : Bfwd[...]` does.
 *
 * Offsets/lengths are `long` (not `int`): a concatenated human-scale genome (~3.1 Gbp total
 * across all contigs) exceeds INT32_MAX (~2.147 Gbp), so genome-absolute addressing must be
 * 64-bit even though any *single* alignment's local extent fits comfortably in int.
 */
#ifndef WAVE_KERNEL_H
#define WAVE_KERNEL_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct wave_ctx wave_ctx;

wave_ctx *wave_open(void);
void      wave_close(wave_ctx *g);

/* Upload both genomes resident. A/Bfwd/Brev are flat NUMERIC (0-3) byte arrays laid out as
 * described above; Abase[0..nA) / Bbase[0..nB) are the per-contig base offsets; Alen/Blen are
 * the flat arrays' total byte lengths (>= 1 + sum(contig len) + ncontig, the sentinel-padded
 * total). Bfwd and Brev share Bbase/nB/Blen (same contig geometry, opposite orientation).
 * Safe to call more than once on the same ctx (reallocates device storage only if it grows).
 */
void wave_load_genomes(wave_ctx *g,
                        const unsigned char *A,    const long *Abase, int nA, long Alen,
                        const unsigned char *Bfwd, const unsigned char *Brev,
                                                    const long *Bbase, int nB, long Blen);

/* Read back [offset,offset+n) from the resident device array into a caller-allocated host
 * buffer `out`. Returns 0 on success, -1 if the range is out of bounds or nothing is loaded
 * yet. For validation (--selftest) and, later, for any host-side spot-check. */
int wave_readback_A(wave_ctx *g, long offset, long n, unsigned char *out);
int wave_readback_Bfwd(wave_ctx *g, long offset, long n, unsigned char *out);
int wave_readback_Brev(wave_ctx *g, long offset, long n, unsigned char *out);

#ifdef __cplusplus
}
#endif
#endif
