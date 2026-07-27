/* Deep-wave GPU throughput: one warp per wave, band across 32 lanes, x-drop discovery.
 * This is a MEASUREMENT tool (borrows the wave writing-style from fastga_gpu.cu's extend_warp)
 * but lifts the two caps that would exclude deep waves:
 *   - positions are int (not short): deep waves reach x up to ~3M >> 32767.
 *   - MAXD is large: depth up to ~560K; the ping-pong keeps only 2 layers so this is free.
 * KB (band half-width) stays 320 -- the profiling shows active band ~128, so 320 is generous.
 *
 * Compares GPU deep-waves/sec against the paired CPU Local_Alignment bench on the SAME .disc.
 * Reports endpoint agreement vs the .1aln reference so the timing is judged apples-to-apples,
 * and reports both naive order and depth-sorted order (to separate raw throughput from the
 * warp-imbalance penalty of the heavy tail).
 *
 * Build: see Makefile target deep_gpu_bench    Run: ./gpu/deep_gpu_bench <tasks.disc>
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cuda_runtime.h>
#include "disc_format.h"

typedef unsigned char u8;
#define KB   320
#define DW   (2*KB+1)
#define MAXD 600000
#define XDROP 40

#define CK(c) do{cudaError_t e=(c); if(e!=cudaSuccess){fprintf(stderr,"CUDA %s:%d %s\n",__FILE__,__LINE__,cudaGetErrorString(e));exit(1);}}while(0)

// int-position furthest-reaching x-drop wave, band swept by the warp (argmax over lanes).
// Same algorithm as fastga_gpu.cu:extend_warp, positions widened short->int.
__device__ static void extend_warp_i(const u8*A,long a0,int da,long alen,
                                      const u8*B,long b0,int db,long blen,
                                      int xdrop,int*f,int*M,int lane,int*ox,int*oy,int*od){
  int *fp=f, *fc=f+DW, *Mp=M, *Mc=M+DW;
  for(int i=lane;i<DW;i+=32){fp[i]=-1;Mp[i]=0;fc[i]=-1;Mc[i]=0;}
  __syncwarp();
  int x0=0,mm0=0;
  if(lane==0){int x=0,mm=0;while(x<alen&&x<blen&&A[a0+da*(long)x]==B[b0+db*(long)x]){x++;mm++;}
              fp[KB]=x;Mp[KB]=mm;x0=x;mm0=mm;}
  __syncwarp();
  int bestS=__shfl_sync(0xffffffff,mm0,0);
  int bx=__shfl_sync(0xffffffff,x0,0), by=bx, bd=0, stop=0;
  for(int d=1;d<=MAXD && !stop;d++){
    for(int i=lane;i<DW;i+=32) fc[i]=-1;
    __syncwarp();
    int kmin=-KB>-d?-KB:-d, kmax=KB<d?KB:d;
    int locS=-1000000,lcx=0,lcy=0;
    for(int k=kmin+lane;k<=kmax;k+=32){
      int idx=k+KB,best=-1,bm=0,v;
      if(k-1>=-KB){v=fp[idx-1];if(v>=0&&v+1>best){best=v+1;bm=Mp[idx-1];}}
      if(k+1<=KB){v=fp[idx+1];if(v>=0&&v  >best){best=v;  bm=Mp[idx+1];}}
      v=fp[idx];if(v>=0&&v+1>best){best=v+1;bm=Mp[idx];}
      if(best<0)continue;
      int xx=best;if(xx>alen)xx=alen;if((long)xx-k>blen)xx=(int)(blen+k);if(xx<0)xx=0;int yy=xx-k;
      while(xx<alen&&yy<blen&&A[a0+da*(long)xx]==B[b0+db*(long)yy]){xx++;yy++;bm++;}
      fc[idx]=xx;Mc[idx]=bm;
      int S=bm-2*d; if(S>locS){locS=S;lcx=xx;lcy=yy;}
    }
    __syncwarp();
    int mS=locS,mx=lcx,my=lcy;
    for(int o=16;o>=1;o>>=1){
      int oS=__shfl_down_sync(0xffffffff,mS,o),oX=__shfl_down_sync(0xffffffff,mx,o),oY=__shfl_down_sync(0xffffffff,my,o);
      if(oS>mS){mS=oS;mx=oX;my=oY;}
    }
    int curbestS=__shfl_sync(0xffffffff,mS,0);
    int cx=__shfl_sync(0xffffffff,mx,0), cy=__shfl_sync(0xffffffff,my,0);
    if(lane==0 && curbestS>bestS){bestS=curbestS;bx=cx;by=cy;bd=d;}
    bestS=__shfl_sync(0xffffffff,bestS,0);
    if(curbestS<bestS-xdrop) stop=1;
    int*t; t=fp;fp=fc;fc=t; t=Mp;Mp=Mc;Mc=t;
  }
  if(lane==0){*ox=bx;*oy=by;*od=bd;}
}

// one warp per wave; SA/SB are absolute offsets into the packed A/B pools; aLo/aHi/bLo/bHi
// bound the wave to this wave's window so x-drop can't run into the neighbouring task.
__global__ void disc_deep(const u8*A,const u8*B,int n,int xdrop,
                          const long*SA,const long*SB,
                          const long*aLo,const long*aHi,const long*bLo,const long*bHi,
                          int*ab,int*ae,int*bb,int*be,int*df){
  extern __shared__ int sh[];
  int warp=(blockIdx.x*blockDim.x+threadIdx.x)>>5, lane=threadIdx.x&31, wib=threadIdx.x>>5;
  if(warp>=n)return;
  int*f=sh+(size_t)wib*(4*DW); int*M=f+2*DW;
  long sa=SA[warp],sb=SB[warp],aL=aLo[warp],aH=aHi[warp],bL=bLo[warp],bH=bHi[warp];
  int fx,fy,fd,rx,ry,rd;
  extend_warp_i(A,sa,  +1,aH-sa,B,sb,  +1,bH-sb,xdrop,f,M,lane,&fx,&fy,&fd);
  extend_warp_i(A,sa-1,-1,sa-aL,B,sb-1,-1,sb-bL,xdrop,f,M,lane,&rx,&ry,&rd);
  if(lane==0){ae[warp]=(int)(sa+fx-aL);be[warp]=(int)(sb+fy-bL);
              ab[warp]=(int)(sa-rx-aL);bb[warp]=(int)(sb-ry-bL);df[warp]=fd+rd;}
}

int main(int argc,char**argv){
  if(argc<2){fprintf(stderr,"usage: %s <tasks.disc>\n",argv[0]);return 1;}
  FILE*f=fopen(argv[1],"rb"); if(!f){fprintf(stderr,"open fail\n");return 1;}
  FGADiscHeader h; if(fread(&h,sizeof(h),1,f)!=1||h.magic!=FGA_DISC_MAGIC){fprintf(stderr,"bad header\n");return 1;}
  int N=(int)h.ntasks;
  std::vector<u8> Ah,Bh; std::vector<long> hSA,hSB,haLo,haHi,hbLo,hbHi;
  std::vector<int> RAB,RAE,RBB,RBE,RDF,AWv;
  for(int i=0;i<N;i++){
    uint32_t aw,bw,sa,sb,rd; int rab,rae,rbb,rbe;
    if(fread(&aw,4,1,f)!=1)break;
    fread(&bw,4,1,f);fread(&sa,4,1,f);fread(&sb,4,1,f);
    fread(&rab,4,1,f);fread(&rae,4,1,f);fread(&rbb,4,1,f);fread(&rbe,4,1,f);fread(&rd,4,1,f);
    long ao=(long)Ah.size(), bo=(long)Bh.size();
    Ah.resize(ao+aw); if(aw)fread(Ah.data()+ao,1,aw,f);
    Bh.resize(bo+bw); if(bw)fread(Bh.data()+bo,1,bw,f);
    if(aw==0||bw==0)continue;
    haLo.push_back(ao); haHi.push_back(ao+aw); hbLo.push_back(bo); hbHi.push_back(bo+bw);
    hSA.push_back(ao+sa); hSB.push_back(bo+sb);
    RAB.push_back(rab);RAE.push_back(rae);RBB.push_back(rbb);RBE.push_back(rbe);RDF.push_back(rd);AWv.push_back(aw);
    (void)rae;(void)rbb;
  }
  fclose(f);
  int n=(int)hSA.size();
  printf("packed %d deep waves; A=%.2f GB B=%.2f GB\n",n,Ah.size()/1e9,Bh.size()/1e9);

  // device pools
  u8*dA,*dB; CK(cudaMalloc(&dA,Ah.size())); CK(cudaMalloc(&dB,Bh.size()));
  CK(cudaMemcpy(dA,Ah.data(),Ah.size(),cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dB,Bh.data(),Bh.size(),cudaMemcpyHostToDevice));

  auto run = [&](std::vector<int>&order,int xdrop,const char*tag)->double{
    // reorder task metadata per `order`
    std::vector<long> SA(n),SB(n),aLo(n),aHi(n),bLo(n),bHi(n);
    for(int i=0;i<n;i++){int j=order[i];SA[i]=hSA[j];SB[i]=hSB[j];aLo[i]=haLo[j];aHi[i]=haHi[j];bLo[i]=hbLo[j];bHi[i]=hbHi[j];}
    long *dSA,*dSB,*daLo,*daHi,*dbLo,*dbHi; int *dab,*dae,*dbb,*dbe,*ddf;
    CK(cudaMalloc(&dSA,n*8));CK(cudaMalloc(&dSB,n*8));CK(cudaMalloc(&daLo,n*8));CK(cudaMalloc(&daHi,n*8));
    CK(cudaMalloc(&dbLo,n*8));CK(cudaMalloc(&dbHi,n*8));
    CK(cudaMalloc(&dab,n*4));CK(cudaMalloc(&dae,n*4));CK(cudaMalloc(&dbb,n*4));CK(cudaMalloc(&dbe,n*4));CK(cudaMalloc(&ddf,n*4));
    CK(cudaMemcpy(dSA,SA.data(),n*8,cudaMemcpyHostToDevice));CK(cudaMemcpy(dSB,SB.data(),n*8,cudaMemcpyHostToDevice));
    CK(cudaMemcpy(daLo,aLo.data(),n*8,cudaMemcpyHostToDevice));CK(cudaMemcpy(daHi,aHi.data(),n*8,cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dbLo,bLo.data(),n*8,cudaMemcpyHostToDevice));CK(cudaMemcpy(dbHi,bHi.data(),n*8,cudaMemcpyHostToDevice));
    int TPB=64, wpb=TPB/32, bl=(n+wpb-1)/wpb; size_t shmem=(size_t)wpb*4*DW*sizeof(int);
    cudaEvent_t e0,e1; cudaEventCreate(&e0);cudaEventCreate(&e1);
    double best=1e30; std::vector<int> ab(n),ae(n),bb(n),be(n),df(n);
    for(int r=0;r<3;r++){
      cudaEventRecord(e0);
      disc_deep<<<bl,TPB,shmem>>>(dA,dB,n,xdrop,dSA,dSB,daLo,daHi,dbLo,dbHi,dab,dae,dbb,dbe,ddf);
      cudaEventRecord(e1); CK(cudaEventSynchronize(e1));
      CK(cudaGetLastError());
      float ms; cudaEventElapsedTime(&ms,e0,e1); if(ms/1e3<best)best=ms/1e3;
    }
    CK(cudaMemcpy(ab.data(),dab,n*4,cudaMemcpyDeviceToHost));CK(cudaMemcpy(ae.data(),dae,n*4,cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(bb.data(),dbb,n*4,cudaMemcpyDeviceToHost));CK(cudaMemcpy(be.data(),dbe,n*4,cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(df.data(),ddf,n*4,cudaMemcpyDeviceToHost));
    // endpoint agreement vs reference (order[i] maps back to original task j)
    long e50=0; double dr=0; int drn=0;
    for(int i=0;i<n;i++){int j=order[i];
      int mx=abs(ab[i]-RAB[j]);int t;t=abs(ae[i]-RAE[j]);if(t>mx)mx=t;t=abs(bb[i]-RBB[j]);if(t>mx)mx=t;t=abs(be[i]-RBE[j]);if(t>mx)mx=t;
      if(mx<=50)e50++; if(RDF[j]>0){dr+=(double)df[i]/RDF[j];drn++;}}
    printf("GPU [%s]: %.3f s  (%.1f aln/s)  | endpoints %ld/%d within 50bp (%.1f%%)  mean df/ref=%.3f\n",
           tag,best,n/best,e50,n,100.0*e50/n,dr/drn);
    cudaFree(dSA);cudaFree(dSB);cudaFree(daLo);cudaFree(daHi);cudaFree(dbLo);cudaFree(dbHi);
    cudaFree(dab);cudaFree(dae);cudaFree(dbb);cudaFree(dbe);cudaFree(ddf);
    return best;
  };

  std::vector<int> naive(n); for(int i=0;i<n;i++)naive[i]=i;
  std::vector<int> sorted=naive;
  std::sort(sorted.begin(),sorted.end(),[&](int a,int b){return RDF[a]>RDF[b];}); // deepest first
  run(sorted,250,"sorted xdrop=250");
  cudaFree(dA);cudaFree(dB);
  return 0;
}
