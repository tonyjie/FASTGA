/* Fidelity-vs-divergence proxy for the FastGA GPU aligner (validation B).
 *
 * Reads a .disc file (real FastGA discovery tasks: seed anchor + window + reference
 * endpoints/diffs), runs the GPU x-drop discovery per task via the fastga_gpu library,
 * and reports the fraction of tasks whose 4 endpoints all land within a tolerance of
 * FastGA's -- BUCKETED BY per-task local divergence (ref_diffs / ref-A-length).
 *
 * This turns a single dataset into a fidelity-vs-divergence curve: the tail of high-diff
 * (tandem-repeat / low-complexity) tubes is exactly where the GPU x-drop is expected to
 * miss and the CPU fallback engages.
 *
 * Build: see Makefile target `disc_fidelity` (links gpu/fastga_gpu.o + -lcudart).
 * Run:   ./disc_fidelity <tasks.disc>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "disc_format.h"
#include "fastga_gpu.h"

// divergence buckets (percent), upper-exclusive; last is open-ended
static const double BUP[] = { 1, 2, 5, 10, 20, 1e9 };
static const char  *BNM[] = { "0-1%", "1-2%", "2-5%", "5-10%", "10-20%", ">20%" };
#define NB (int)(sizeof(BUP)/sizeof(BUP[0]))

static int absi(int x){ return x<0?-x:x; }

int main(int argc, char **argv)
{ if (argc < 2) { fprintf(stderr,"usage: %s <tasks.disc>\n",argv[0]); return 1; }
  FILE *f = fopen(argv[1],"rb");
  if (!f) { fprintf(stderr,"cannot open %s\n",argv[1]); return 1; }
  FGADiscHeader h;
  if (fread(&h,sizeof(h),1,f)!=1 || h.magic!=FGA_DISC_MAGIC)
    { fprintf(stderr,"bad .disc header\n"); return 1; }

  long n50[NB]={0}, n10[NB]={0}, ntot[NB]={0};
  long g50=0, g10=0, gtot=0;
  gpu_ctx *g = gpu_open();

  unsigned char *A=NULL, *B=NULL; size_t acap=0, bcap=0;
  for (uint32_t i=0; i<h.ntasks; i++)
    { uint32_t aw,bw,sa,sb,rd; int rab,rae,rbb,rbe;
      if (fread(&aw,4,1,f)!=1) break;
      fread(&bw,4,1,f); fread(&sa,4,1,f); fread(&sb,4,1,f);
      fread(&rab,4,1,f); fread(&rae,4,1,f); fread(&rbb,4,1,f); fread(&rbe,4,1,f); fread(&rd,4,1,f);
      if (aw>acap){ A=(unsigned char*)realloc(A,aw); acap=aw; }
      if (bw>bcap){ B=(unsigned char*)realloc(B,bw); bcap=bw; }
      if (aw) fread(A,1,aw,f);
      if (bw) fread(B,1,bw,f);

      int reflen = rae-rab; if (reflen<=0) continue;
      double div = 100.0*(double)rd/(double)reflen;
      int b=0; while (b<NB-1 && div>=BUP[b]) b++;

      int isa=(int)sa,isb=(int)sb,ab,ae,bb,be,df;
      gpu_load_seqs(g,A,aw,B,bw);
      gpu_discover_batch(g,1,&isa,&isb,&ab,&ae,&bb,&be,&df);
      int mx = absi(ab-rab); if(absi(ae-rae)>mx)mx=absi(ae-rae);
      if(absi(bb-rbb)>mx)mx=absi(bb-rbb); if(absi(be-rbe)>mx)mx=absi(be-rbe);

      ntot[b]++; gtot++;
      if (mx<=50){ n50[b]++; g50++; }
      if (mx<=10){ n10[b]++; g10++; }
    }
  gpu_close(g);
  fclose(f);

  printf("GPU discovery fidelity vs local divergence  (%s, %ld tasks)\n", argv[1], gtot);
  printf("  %-8s %10s %12s %12s\n","div","tasks","<=50bp","<=10bp");
  for (int b=0;b<NB;b++)
    if (ntot[b])
      printf("  %-8s %10ld %11.1f%% %11.1f%%\n", BNM[b], ntot[b],
             100.0*n50[b]/ntot[b], 100.0*n10[b]/ntot[b]);
  printf("  %-8s %10ld %11.1f%% %11.1f%%\n","ALL",gtot,
         100.0*g50/gtot, 100.0*g10/gtot);
  return 0;
}
