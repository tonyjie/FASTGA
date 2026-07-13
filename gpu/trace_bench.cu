/* A3b feasibility bench: how fast can the GPU trace ALL alignments of a part, batched?
 *
 * Packs every .trace task into one resident A/B buffer (phase-aligned, as trace_lib_test),
 * then times gpu_trace_batch over the whole set (H2D endpoints + kernel + D2H), repeated.
 * Reports GPU alignments/sec.  Compare against the CPU align-phase throughput
 * (EXAMPLE: 323,569 alignments in 16.7 s wall at T=8 -> ~19,400 aln/s) to decide whether the
 * search_seeds batching redesign (A3b) can beat the CPU aligner before investing in it.
 *
 * Build: make trace_bench    Run: ./gpu/trace_bench <tasks.trace>
 */
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include "trace_format.h"
#include "fastga_gpu.h"
typedef unsigned char u8;
using clk = std::chrono::high_resolution_clock;

int main(int argc,char**argv){
  if(argc<2){fprintf(stderr,"usage: %s <tasks.trace>\n",argv[0]);return 1;}
  FILE*f=fopen(argv[1],"rb"); if(!f){fprintf(stderr,"open fail\n");return 1;}
  FGATraceHeader h;
  if(fread(&h,sizeof(h),1,f)!=1||h.magic!=FGA_TRACE_MAGIC){fprintf(stderr,"bad header\n");return 1;}
  int ts=(int)h.tspace;

  std::vector<u8> Abuf,Bbuf; std::vector<int> AB,AE,BB,BE;
  long tot_aw=0;
  for(uint32_t i=0;i<h.ntasks;i++){
    int32_t ab,ae,bb,be; uint32_t aw,bw,rd,rl;
    if(fread(&ab,4,1,f)!=1)break;
    fread(&ae,4,1,f);fread(&bb,4,1,f);fread(&be,4,1,f);
    fread(&aw,4,1,f);fread(&bw,4,1,f);fread(&rd,4,1,f);fread(&rl,4,1,f);
    std::vector<u8> A(aw),B(bw); if(aw)fread(A.data(),1,aw,f); if(bw)fread(B.data(),1,bw,f);
    if(rl){ if(fseek(f,(long)rl*2,SEEK_CUR)){} }
    if(aw==0||bw==0)continue;
    int r=((ab%ts)+ts)%ts;
    while((int)(Abuf.size()%ts)!=r)Abuf.push_back(0);
    int ao=(int)Abuf.size();Abuf.insert(Abuf.end(),A.begin(),A.end());
    int bo=(int)Bbuf.size();Bbuf.insert(Bbuf.end(),B.begin(),B.end());
    AB.push_back(ao);AE.push_back(ao+(int)aw);BB.push_back(bo);BE.push_back(bo+(int)bw);
    tot_aw+=aw;
  }
  fclose(f);
  int N=(int)AB.size();
  printf("packed %d alignments; resident A=%.1f MB B=%.1f MB; mean A-len=%.0f\n",
         N,Abuf.size()/1e6,Bbuf.size()/1e6,(double)tot_aw/N);

  gpu_ctx*g=gpu_open();
  gpu_load_seqs(g,Abuf.data(),(int)Abuf.size(),Bbuf.data(),(int)Bbuf.size());
  std::vector<unsigned short> out((size_t)N*FGA_TRACE_MAX_PAIRS);
  std::vector<int> tlen(N);

  double best=1e30;
  for(int rep=0;rep<3;rep++){
    auto t0=clk::now();
    int rc=gpu_trace_batch(g,N,AB.data(),AE.data(),BB.data(),BE.data(),ts,out.data(),tlen.data());
    auto t1=clk::now();
    if(rc){fprintf(stderr,"gpu_trace_batch rc=%d\n",rc);return 1;}
    double s=std::chrono::duration<double>(t1-t0).count();
    if(s<best)best=s;
    printf("  rep %d: %.3f s  (%.0f aln/s)\n",rep,s,N/s);
  }
  long ov=0; for(int i=0;i<N;i++) if(tlen[i]<0)ov++;
  gpu_close(g);
  printf("BEST: %.3f s for %d alignments  =>  %.0f aln/s   [overflow=%ld]\n",best,N,N/best,ov);
  printf("CPU align phase (EXAMPLE T=8): 16.7 s for 323,569 alns => ~19,400 aln/s\n");
  printf("GPU trace-only speedup vs CPU align phase (trace is only part of A3b): %.1fx\n",
         16.7/best);
  return 0;
}
