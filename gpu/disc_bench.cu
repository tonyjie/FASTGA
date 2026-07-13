/* A3b feasibility: time batched gpu_discover_batch over all discovery tasks of a part.
 * Packs every .disc task's window into one resident buffer, seeds -> absolute coords,
 * times gpu_discover_batch. Reports GPU discoveries/sec (the OTHER GPU pass in A3b,
 * alongside the trace bench). Build: make disc_bench   Run: ./gpu/disc_bench <tasks.disc>
 */
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include "disc_format.h"
#include "fastga_gpu.h"
typedef unsigned char u8;
using clk=std::chrono::high_resolution_clock;

int main(int argc,char**argv){
  if(argc<2){fprintf(stderr,"usage: %s <tasks.disc>\n",argv[0]);return 1;}
  FILE*f=fopen(argv[1],"rb"); if(!f){fprintf(stderr,"open fail\n");return 1;}
  FGADiscHeader h;
  if(fread(&h,sizeof(h),1,f)!=1||h.magic!=FGA_DISC_MAGIC){fprintf(stderr,"bad header\n");return 1;}
  std::vector<u8> Abuf,Bbuf; std::vector<int> SA,SB;
  for(uint32_t i=0;i<h.ntasks;i++){
    uint32_t aw,bw,sa,sb,rd; int rab,rae,rbb,rbe;
    if(fread(&aw,4,1,f)!=1)break;
    fread(&bw,4,1,f);fread(&sa,4,1,f);fread(&sb,4,1,f);
    fread(&rab,4,1,f);fread(&rae,4,1,f);fread(&rbb,4,1,f);fread(&rbe,4,1,f);fread(&rd,4,1,f);
    std::vector<u8> A(aw),B(bw); if(aw)fread(A.data(),1,aw,f); if(bw)fread(B.data(),1,bw,f);
    if(aw==0||bw==0)continue;
    int ao=(int)Abuf.size();Abuf.insert(Abuf.end(),A.begin(),A.end());
    int bo=(int)Bbuf.size();Bbuf.insert(Bbuf.end(),B.begin(),B.end());
    SA.push_back(ao+(int)sa); SB.push_back(bo+(int)sb);
  }
  fclose(f);
  int N=(int)SA.size();
  printf("packed %d discovery tasks; resident A=%.1f MB B=%.1f MB\n",N,Abuf.size()/1e6,Bbuf.size()/1e6);
  gpu_ctx*g=gpu_open();
  gpu_load_seqs(g,Abuf.data(),(int)Abuf.size(),Bbuf.data(),(int)Bbuf.size());
  std::vector<int> ab(N),ae(N),bb(N),be(N),df(N);
  double best=1e30;
  for(int r=0;r<3;r++){
    auto t0=clk::now();
    int rc=gpu_discover_batch(g,N,SA.data(),SB.data(),ab.data(),ae.data(),bb.data(),be.data(),df.data());
    auto t1=clk::now();
    if(rc){fprintf(stderr,"rc=%d\n",rc);return 1;}
    double s=std::chrono::duration<double>(t1-t0).count();
    if(s<best)best=s;
    printf("  rep %d: %.3f s (%.0f disc/s)\n",r,s,N/s);
  }
  gpu_close(g);
  printf("BEST: %.3f s for %d discoveries => %.0f disc/s\n",best,N,N/best);
  return 0;
}
