/* Clean wave-vs-wave baseline: time CPU Compute_Alignment(DIFF_TRACE) at 32 threads on the
 * SAME .trace alignments the GPU trace_bench uses. Compare aln/s directly (no sort/output
 * confound). Build: make cpu_trace_bench   Run: ./gpu/cpu_trace_bench <tasks.trace> [nthreads]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#include "align.h"
#include "trace_format.h"

int main(int argc,char**argv){
  if(argc<2){fprintf(stderr,"usage: %s <tasks.trace> [nthreads=32]\n",argv[0]);return 1;}
  int NT = argc>=3 ? atoi(argv[2]) : 32;
  FILE*f=fopen(argv[1],"rb"); if(!f){fprintf(stderr,"open fail\n");return 1;}
  FGATraceHeader h;
  if(fread(&h,sizeof(h),1,f)!=1||h.magic!=FGA_TRACE_MAGIC){fprintf(stderr,"bad header\n");return 1;}
  int ts=(int)h.tspace;

  // load: keep each task's A/B rectangle contiguous
  int cap=h.ntasks;
  int *AW=malloc(cap*sizeof(int)),*BW=malloc(cap*sizeof(int)),*AO=malloc(cap*sizeof(int)),*BO=malloc(cap*sizeof(int));
  size_t abytes=0,bbytes=0, acp=1<<28,bcp=1<<28;
  unsigned char *Ab=malloc(acp),*Bb=malloc(bcp);
  int N=0;
  for(uint32_t i=0;i<h.ntasks;i++){
    int32_t ab,ae,bb,be; uint32_t aw,bw,rd,rl;
    if(fread(&ab,4,1,f)!=1)break;
    fread(&ae,4,1,f);fread(&bb,4,1,f);fread(&be,4,1,f);
    fread(&aw,4,1,f);fread(&bw,4,1,f);fread(&rd,4,1,f);fread(&rl,4,1,f);
    while(abytes+aw>acp){acp*=2;Ab=realloc(Ab,acp);}
    while(bbytes+bw>bcp){bcp*=2;Bb=realloc(Bb,bcp);}
    if(aw)fread(Ab+abytes,1,aw,f); if(bw)fread(Bb+bbytes,1,bw,f);
    if(rl)fseek(f,(long)rl*2,SEEK_CUR);
    if(aw==0||bw==0)continue;
    AO[N]=abytes;BO[N]=bbytes;AW[N]=aw;BW[N]=bw;
    abytes+=aw;bbytes+=bw;N++;
  }
  fclose(f);
  printf("loaded %d alignments\n",N);

  omp_set_num_threads(NT);
  double t0=omp_get_wtime();
  long overflow=0;
  #pragma omp parallel reduction(+:overflow)
  { Work_Data *wk=New_Work_Data();
    Alignment aln; Path path;
    aln.path=&path;
    #pragma omp for schedule(dynamic,256)
    for(int i=0;i<N;i++){
      aln.aseq=(char*)(Ab+AO[i]); aln.bseq=(char*)(Bb+BO[i]);
      aln.alen=AW[i]; aln.blen=BW[i]; aln.flags=0;
      path.abpos=0; path.aepos=AW[i]; path.bbpos=0; path.bepos=BW[i];
      if(Compute_Alignment(&aln,wk,DIFF_TRACE,ts)) overflow++;
    }
    Free_Work_Data(wk);
  }
  double s=omp_get_wtime()-t0;
  printf("CPU %d-thread Compute_Alignment(DIFF_TRACE): %.3f s  (%.0f aln/s)  [errors=%ld]\n",
         NT,s,N/s,overflow);
  return 0;
}
