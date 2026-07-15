/* C-callable GPU local-alignment offload for FastGA (G4 integration).
 *
 * FastGA's align_contigs offloads its per-tube Local_Alignment DISCOVERY (find
 * endpoints+diffs from a seed) to the A100 in batches, then regenerates trace
 * points on the CPU (Compute_Trace_PTS) from the returned endpoints.
 *
 * Usage per contig-pair:
 *   gpu_ctx *g = gpu_open();
 *   gpu_load_seqs(g, A, alen, B, blen);      // NUMERIC (0-3); resident, reused per tube batch
 *   ... collect N tubes as seed anchors (sa[i], sb[i]) in A/B contig coords ...
 *   gpu_discover_batch(g, N, sa, sb, ab, ae, bb, be, diffs);   // fills the 5 output arrays
 *   ... CPU: filter + Compute_Trace_PTS + write .las ...
 *   gpu_close(g);
 */
#ifndef FASTGA_GPU_H
#define FASTGA_GPU_H
#ifdef __cplusplus
extern "C" {
#endif

typedef struct gpu_ctx gpu_ctx;

gpu_ctx *gpu_open(void);
void     gpu_close(gpu_ctx *g);

/* Copy the two contig sequences (NUMERIC, 1 byte/base 0-3) to the device; they
 * stay resident and are reused across every gpu_discover_batch for this pair. */
void gpu_load_seqs(gpu_ctx *g, const unsigned char *A, int alen,
                              const unsigned char *B, int blen);

/* Same, but the inputs are the 2-bit PACKED contigs (from .bps, (len+3)/4 bytes each):
 * the device unpacks them on-GPU (byte-identical to Uncompress_Read), so the CPU never
 * decompresses. Eliminates the 29%-of-runtime Uncompress_Read hotspot and shrinks H2D 4x. */
void gpu_load_seqs_2bit(gpu_ctx *g, const unsigned char *packedA, int alen,
                                    const unsigned char *packedB, int blen);

/* Discover N local alignments: from each seed anchor (sa[i],sb[i]) (A/B contig
 * coords) run forward+reverse x-drop wave -> endpoints (ab,ae,bb,be) + diffs.
 * Output arrays are caller-allocated, length >= N. Returns 0 on success. */
int gpu_discover_batch(gpu_ctx *g, int n,
                       const int *sa, const int *sb,
                       int *ab, int *ae, int *bb, int *be, int *diffs);

/* Emit FastGA trace-points for N alignments whose endpoints (ab,ae,bb,be) index
 * into the resident A/B contigs.  out_trace is caller-allocated as N fixed slots
 * of FGA_TRACE_MAX_PAIRS uint16 each (out_trace + i*FGA_TRACE_MAX_PAIRS); out_tlen[i]
 * is the used length (2*npanels), or negative on overflow: -1 band, -2 depth, -3 panels.
 * Trace pairs are (diffs, delta-b) per global-tspace panel (same as Compute_Alignment).
 * Returns 0 on success. */
#define FGA_TRACE_MAX_PAIRS 512
int gpu_trace_batch(gpu_ctx *g, int n,
                    const int *ab, const int *ae, const int *bb, const int *be,
                    int tspace, unsigned short *out_trace, int *out_tlen);

/* Genome-resident variants (A3b): A/B are whole genomes loaded via gpu_load_seqs_2bit; all
 * coords are genome-relative. discover takes per-tube contig bounds [aLo,aHi)x[bLo,bHi) so the
 * x-drop cannot cross a contig; trace takes per-tube base[i] = A-contig base offset for
 * contig-relative panel phasing. */
int gpu_discover_batch_g(gpu_ctx *g, int n, const int *sa, const int *sb,
                         const int *aLo, const int *aHi, const int *bLo, const int *bHi,
                         int *ab, int *ae, int *bb, int *be, int *diffs);
int gpu_trace_batch_g(gpu_ctx *g, int n,
                      const int *ab, const int *ae, const int *bb, const int *be,
                      const int *base, int tspace, unsigned short *out_trace, int *out_tlen);

#ifdef __cplusplus
}
#endif
#endif
