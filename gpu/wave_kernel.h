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

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct wave_ctx wave_ctx;

/* One seed to sweep, contig-relative (mirrors wave_bench_cpu.c's per-seed setup):
 *   aread/bread : contig indices into the resident A / B genomes
 *   alen/blen   : contig lengths (bound reads to [base,base+len))
 *   anti/diag   : seed anti-diagonal (a+b) and diagonal (a-b); the forward wave starts from
 *                 band [diag,diag] at anti-diagonal `anti` (exactly Local_Alignment's low=hgh=diag)
 *   comp        : 1 => use Brev (B reverse-complemented), 0 => use Bfwd -- the COMP(flags) pick */
typedef struct
  { int aread, bread;
    int alen, blen;
    int anti, diag;
    int comp;
  } wave_seed;

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

/* Stage-1 FORWARD sweep (Task 5): a faithful, warp-cooperative port of align.c's forward_wave.
 * One warp per seed; 32 lanes tile the diagonal band [low,hgh]; the `dif` loop is sequential
 * (strict antidiagonal recurrence) while the per-diagonal furthest-reach slide runs in parallel
 * across lanes, and a serial lane-0 pass reproduces forward_wave's exact k-descending
 * besta/trima trim tie-breaking (TABLE/SCORE suffix-positivity, TRIM_MLAG=250 continuation).
 *
 * Produces, per seed i, the FORWARD endpoint + forward diff count, contig-relative (identical
 * frame to align.c's apath->aepos/bepos/diffs, since forward_wave is an independent sweep from
 * the seed):
 *     ae[i]    = apath->aepos (= trimx)
 *     be[i]    = apath->bepos (= trimy = trima - trimx)
 *     fdiff[i] = apath->diffs (= trimd), the forward-half diff count
 * A seed whose band exceeds the per-warp scratch width is reported as ae[i]=be[i]=-2 (overflow;
 * count reported by the caller). `table`/`score` are the 32768-entry int16 suffix-positivity
 * tables from set_table (New_Align_Spec); path_ave is spec->ave_path. Genome residency must be
 * loaded (wave_load_genomes) first. Returns 0 on success, -1 on a fatal launch/alloc error. */
int wave_forward_batch(wave_ctx *g, int n, const wave_seed *seeds,
                        const int16_t *table, const int16_t *score, int path_ave,
                        int *ae, int *be, int *fdiff);

/* Stage-1 FULL discovery (Task 6): forward sweep + reverse sweep combined into the complete
 * Local_Alignment-equivalent endpoint set.  The reverse sweep is a faithful mirror of align.c's
 * reverse_wave (align.c:919-1459) -- opposite extension direction (x-=1), retired-diagonal
 * sentinel INT32_MAX (reverse minimizes reach, so a retired diagonal must read artificially
 * HIGH), reverse trim continuation `lasta <= besta + TRIM_MLAG`, minimizing best test `c<besta`,
 * and the swapped low/high clip roles (align.c:1316-1341).
 *
 * The reverse sweep starts from band [seed.diag, seed.diag] at the seed anti-diagonal, exactly
 * as Local_Alignment feeds reverse_wave(low,low,anti,...): because wave_bench always starts the
 * forward wave from a SINGLE diagonal [dg,dg], every forward trim-lineage roots at dg, so
 * forward_wave's *mind (align.c:871-875, cells[root].diag) is identically dg -- forward and
 * reverse are both anchored at the seed point.
 *
 * Produces, per seed i, the full endpoint quintet, contig-relative (same frame as .cpuref /
 * Local_Alignment's path):
 *     ab[i]    = apath->abpos  (reverse trimx)
 *     ae[i]    = apath->aepos  (forward trimx)
 *     bb[i]    = apath->bbpos  (reverse trima - trimx)
 *     be[i]    = apath->bepos  (forward trima - trimx)
 *     diffs[i] = apath->diffs  (forward trimd + reverse trimd)
 * A seed whose band overflows the per-warp scratch in either sweep is reported as all -2.
 * `table`/`score`/`path_ave` and residency requirements are as for wave_forward_batch.
 * Returns 0 on success, -1 on a fatal launch/alloc error. */
int wave_discover_batch(wave_ctx *g, int n, const wave_seed *seeds,
                        const int16_t *table, const int16_t *score, int path_ave,
                        int *ab, int *ae, int *bb, int *be, int *diffs);

#ifdef __cplusplus
}
#endif
#endif
