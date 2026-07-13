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

/* Discover N local alignments: from each seed anchor (sa[i],sb[i]) (A/B contig
 * coords) run forward+reverse x-drop wave -> endpoints (ab,ae,bb,be) + diffs.
 * Output arrays are caller-allocated, length >= N. Returns 0 on success. */
int gpu_discover_batch(gpu_ctx *g, int n,
                       const int *sa, const int *sb,
                       int *ab, int *ae, int *bb, int *be, int *diffs);

#ifdef __cplusplus
}
#endif
#endif
