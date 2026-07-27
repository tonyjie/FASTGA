/* Deep-wave CPU baseline: time FastGA's REAL Local_Alignment (the actual hot-path aligner,
 * NOT Compute_Alignment) on the deep waves in a .disc file. Each task's margin-grown window
 * is treated as a mini-contig; the seed (sa,sb) is fed as diag=sa-sb, anti=sa+sb, exactly as
 * align_contigs does (FastGA.c:3169/3493 cross path, lbord=hbord=-1).
 *
 * Reports single-thread aln/s (per-core rate) and N-thread aln/s (real aggregate, dynamic
 * schedule since deep waves are imbalanced), plus endpoint/diffs agreement vs the .1aln
 * reference (so the paired GPU bench can be judged apples-to-apples).
 *
 * Build: see Makefile target deep_cpu_bench
 * Run:   ./gpu/deep_cpu_bench <tasks.disc> [nthreads=28]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include "align.h"
#include "disc_format.h"

typedef unsigned char u8;

int main(int argc, char **argv)
{ if (argc < 2) { fprintf(stderr,"usage: %s <tasks.disc> [nthreads=28]\n",argv[0]); return 1; }
  int NT = argc>=3 ? atoi(argv[2]) : 28;

  FILE *f = fopen(argv[1],"rb");
  if (!f) { fprintf(stderr,"open fail\n"); return 1; }
  FGADiscHeader h;
  if (fread(&h,sizeof(h),1,f)!=1 || h.magic!=FGA_DISC_MAGIC) { fprintf(stderr,"bad header\n"); return 1; }

  int   N   = (int) h.ntasks;
  int  *AW  = malloc(N*sizeof(int)), *BW = malloc(N*sizeof(int));
  int  *SA  = malloc(N*sizeof(int)), *SB = malloc(N*sizeof(int));
  long *AO  = malloc(N*sizeof(long)), *BO = malloc(N*sizeof(long));
  int  *RAB = malloc(N*sizeof(int)), *RAE = malloc(N*sizeof(int));
  int  *RBB = malloc(N*sizeof(int)), *RBE = malloc(N*sizeof(int));
  int  *RDF = malloc(N*sizeof(int));

  //  Contiguous sequence pools with a shared sentinel byte (value 4) between/around every task,
  //  so aseq[-1]==4 and aseq[len]==4 exactly like New_Contig_Buffer + Get_Contig produce.
  size_t acap = 1<<28, bcap = 1<<28, an = 0, bn = 0;
  u8 *Ab = malloc(acap), *Bb = malloc(bcap);
  Ab[an++] = 4;   // leading sentinel
  Bb[bn++] = 4;

  int Nok = 0;
  for (int i = 0; i < N; i++)
    { uint32_t aw,bw,sa,sb,rd; int rab,rae,rbb,rbe;
      if (fread(&aw,4,1,f)!=1) break;
      fread(&bw,4,1,f); fread(&sa,4,1,f); fread(&sb,4,1,f);
      fread(&rab,4,1,f); fread(&rae,4,1,f); fread(&rbb,4,1,f); fread(&rbe,4,1,f); fread(&rd,4,1,f);
      while (an+aw+1 > acap) { acap*=2; Ab=realloc(Ab,acap); }
      while (bn+bw+1 > bcap) { bcap*=2; Bb=realloc(Bb,bcap); }
      long ao = an, bo = bn;
      if (aw) { fread(Ab+an,1,aw,f); an+=aw; } Ab[an++]=4;   // task bytes + trailing sentinel
      if (bw) { fread(Bb+bn,1,bw,f); bn+=bw; } Bb[bn++]=4;
      if (aw==0||bw==0) continue;
      AW[Nok]=aw; BW[Nok]=bw; SA[Nok]=sa; SB[Nok]=sb; AO[Nok]=ao; BO[Nok]=bo;
      RAB[Nok]=rab; RAE[Nok]=rae; RBB[Nok]=rbb; RBE[Nok]=rbe; RDF[Nok]=rd;
      Nok++;
    }
  fclose(f);
  N = Nok;

  //  Window/depth distribution (reference diffs = a lower bound on wave depth).
  long tw=0; int mnw=1<<30,mxw=0, mnd=1<<30,mxd=0; double td=0;
  for (int i=0;i<N;i++){ int w=AW[i]>BW[i]?AW[i]:BW[i]; tw+=w; if(w<mnw)mnw=w; if(w>mxw)mxw=w;
                         if(RDF[i]<mnd)mnd=RDF[i]; if(RDF[i]>mxd)mxd=RDF[i]; td+=RDF[i]; }
  printf("loaded %d deep-wave tasks\n",N);
  printf("  window (max A/B):  min %d  mean %ld  max %d bp\n", mnw, tw/N, mxw);
  printf("  ref diffs (depth): min %d  mean %.0f  max %d\n", mnd, td/N, mxd);

  static float freq[4] = {0.25f,0.25f,0.25f,0.25f};

  //  ---- single-thread: per-core rate (skip with SKIP_SINGLE=1 for fast multi-thread reruns) ----
  if (getenv("SKIP_SINGLE") == NULL)
  { Align_Spec *spec = New_Align_Spec(0.7, 100, freq, 0);   // 1.-ALIGN_RATE(.3)=.7, TS=100
    Work_Data  *work = New_Work_Data();
    Alignment aln; Path path; aln.path=&path;
    long e50=0; double sumdiff_ratio=0;
    double t0 = omp_get_wtime();
    for (int i=0;i<N;i++)
      { aln.aseq=(char*)(Ab+AO[i]); aln.bseq=(char*)(Bb+BO[i]);
        aln.alen=AW[i]; aln.blen=BW[i]; aln.flags=0;
        int dg = SA[i]-SB[i], anti = SA[i]+SB[i];
        if (Local_Alignment(&aln,work,spec,dg,dg,anti,-1,-1)) { fprintf(stderr,"LA err %d\n",i); return 1; }
        int mx = abs(path.abpos-RAB[i]); int t;
        t=abs(path.aepos-RAE[i]); if(t>mx)mx=t; t=abs(path.bbpos-RBB[i]); if(t>mx)mx=t;
        t=abs(path.bepos-RBE[i]); if(t>mx)mx=t;
        if (mx<=50) e50++;
        if (RDF[i]>0) sumdiff_ratio += (double)path.diffs/RDF[i];
      }
    double s = omp_get_wtime()-t0;
    printf("\nCPU  1-thread Local_Alignment: %.3f s  (%.1f aln/s per core)\n", s, N/s);
    printf("  endpoint agreement vs .1aln: %ld/%d within 50bp (%.1f%%);  mean diffs(LA)/diffs(ref)=%.3f\n",
           e50, N, 100.0*e50/N, sumdiff_ratio/N);
    Free_Work_Data(work); Free_Align_Spec(spec);
  }

  //  ---- N-thread: real aggregate on the deep tail (dynamic schedule) ----
  { omp_set_num_threads(NT);
    double t0 = omp_get_wtime();
    #pragma omp parallel
    { Align_Spec *spec = New_Align_Spec(0.7, 100, freq, 0);
      Work_Data  *work = New_Work_Data();
      Alignment aln; Path path; aln.path=&path;
      #pragma omp for schedule(dynamic,16)
      for (int i=0;i<N;i++)
        { aln.aseq=(char*)(Ab+AO[i]); aln.bseq=(char*)(Bb+BO[i]);
          aln.alen=AW[i]; aln.blen=BW[i]; aln.flags=0;
          int dg = SA[i]-SB[i], anti = SA[i]+SB[i];
          Local_Alignment(&aln,work,spec,dg,dg,anti,-1,-1);
        }
      Free_Work_Data(work); Free_Align_Spec(spec);
    }
    double s = omp_get_wtime()-t0;
    printf("CPU %2d-thread Local_Alignment: %.3f s  (%.1f aln/s aggregate)\n", NT, s, N/s);
  }
  return 0;
}
