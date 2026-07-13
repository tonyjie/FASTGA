/* C-callable GPU local-alignment offload for FastGA (see fastga_gpu.h).
 * Persistent device context: contig sequences stay resident; tube batches stream.
 * Build (library, into FastGA):  nvcc -O3 -arch=sm_80 -c fastga_gpu.cu
 * Build (self-test):  nvcc -O3 -arch=sm_80 -DGPU_LIB_TEST -Xcompiler -fopenmp \
 *                          -o fastga_gpu_test fastga_gpu.cu
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cuda_runtime.h>
#include "fastga_gpu.h"

typedef unsigned char u8;
#define KB 320
#define DW (2*KB+1)
#define MAXD 8192

#define CK(c) do{cudaError_t e=(c); if(e!=cudaSuccess){fprintf(stderr,"CUDA %s:%d %s\n",__FILE__,__LINE__,cudaGetErrorString(e));}}while(0)

// --- device local-alignment (same wave as gpu_discover.cu) ---
__device__ static void extend(const u8*A,int a0,int da,int alen,const u8*B,int b0,int db,int blen,
                              int*ox,int*oy,int*od){
  short f[2][DW], M[2][DW];
  int cur=0,prv=1;
  for(int i=0;i<DW;i++){f[prv][i]=-1;M[prv][i]=0;}
  int x=0,mm=0; while(x<alen&&x<blen&&A[a0+da*x]==B[b0+db*x]){x++;mm++;}
  for(int i=0;i<DW;i++){f[cur][i]=-1;M[cur][i]=0;}
  f[cur][KB]=(short)x; M[cur][KB]=mm;
  int bestS=mm,bx=x,by=x,bd=0; {int t=cur;cur=prv;prv=t;}
  for(int d=1;d<=MAXD;d++){
    for(int i=0;i<DW;i++) f[cur][i]=-1;
    int kmin=-KB>-d?-KB:-d,kmax=KB<d?KB:d,curbestS=-1000000,cx=0,cy=0;
    for(int k=kmin;k<=kmax;k++){
      int best=-1,bm=0,v;
      if(k-1>=-KB){v=f[prv][k-1+KB];if(v>=0&&v+1>best){best=v+1;bm=M[prv][k-1+KB];}}
      if(k+1<=KB){v=f[prv][k+1+KB];if(v>=0&&v>best){best=v;bm=M[prv][k+1+KB];}}
      v=f[prv][k+KB];if(v>=0&&v+1>best){best=v+1;bm=M[prv][k+KB];}
      if(best<0)continue;
      int xx=best;if(xx>alen)xx=alen;if(xx-k>blen)xx=blen+k;if(xx<0)xx=0;int yy=xx-k;
      while(xx<alen&&yy<blen&&A[a0+da*xx]==B[b0+db*yy]){xx++;yy++;bm++;}
      f[cur][k+KB]=(short)xx;M[cur][k+KB]=bm;
      int S=bm-2*d; if(S>curbestS){curbestS=S;cx=xx;cy=yy;}
    }
    if(curbestS>bestS){bestS=curbestS;bx=cx;by=cy;bd=d;}
    if(curbestS<bestS-40) break;
    int t=cur;cur=prv;prv=t;
  }
  *ox=bx;*oy=by;*od=bd;
}
// discover from a seed anchor (sa,sb) in the resident A/B contigs
__global__ void disc_batch(const u8*A,int alen,const u8*B,int blen,int n,
                           const int*sa,const int*sb,int*ab,int*ae,int*bb,int*be,int*df){
  int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=n)return;
  int SA=sa[i],SB=sb[i],fx,fy,fd,rx,ry,rd;
  extend(A,SA,+1,alen-SA,B,SB,+1,blen-SB,&fx,&fy,&fd);
  extend(A,SA-1,-1,SA,   B,SB-1,-1,SB,   &rx,&ry,&rd);
  ae[i]=SA+fx; be[i]=SB+fy; ab[i]=SA-rx; bb[i]=SB-ry; df[i]=fd+rd;
}

struct gpu_ctx { u8 *dA,*dB; int alen,blen; size_t acap,bcap;
                 int *dsa,*dsb,*dab,*dae,*dbb,*dbe,*ddf; int ncap; };

extern "C" gpu_ctx *gpu_open(void){
  gpu_ctx *g=(gpu_ctx*)calloc(1,sizeof(gpu_ctx)); return g;
}
extern "C" void gpu_close(gpu_ctx *g){
  if(!g)return;
  if(g->dA)cudaFree(g->dA); if(g->dB)cudaFree(g->dB);
  if(g->dsa){cudaFree(g->dsa);cudaFree(g->dsb);cudaFree(g->dab);cudaFree(g->dae);cudaFree(g->dbb);cudaFree(g->dbe);cudaFree(g->ddf);}
  free(g);
}
extern "C" void gpu_load_seqs(gpu_ctx *g,const u8*A,int alen,const u8*B,int blen){
  if((size_t)alen>g->acap){ if(g->dA)cudaFree(g->dA); CK(cudaMalloc(&g->dA,alen)); g->acap=alen; }
  if((size_t)blen>g->bcap){ if(g->dB)cudaFree(g->dB); CK(cudaMalloc(&g->dB,blen)); g->bcap=blen; }
  CK(cudaMemcpy(g->dA,A,alen,cudaMemcpyHostToDevice));
  CK(cudaMemcpy(g->dB,B,blen,cudaMemcpyHostToDevice));
  g->alen=alen; g->blen=blen;
}
static void ensure_n(gpu_ctx*g,int n){
  if(n<=g->ncap)return;
  if(g->dsa){cudaFree(g->dsa);cudaFree(g->dsb);cudaFree(g->dab);cudaFree(g->dae);cudaFree(g->dbb);cudaFree(g->dbe);cudaFree(g->ddf);}
  CK(cudaMalloc(&g->dsa,n*4));CK(cudaMalloc(&g->dsb,n*4));CK(cudaMalloc(&g->dab,n*4));
  CK(cudaMalloc(&g->dae,n*4));CK(cudaMalloc(&g->dbb,n*4));CK(cudaMalloc(&g->dbe,n*4));CK(cudaMalloc(&g->ddf,n*4));
  g->ncap=n;
}
extern "C" int gpu_discover_batch(gpu_ctx*g,int n,const int*sa,const int*sb,
                                  int*ab,int*ae,int*bb,int*be,int*diffs){
  ensure_n(g,n);
  CK(cudaMemcpy(g->dsa,sa,n*4,cudaMemcpyHostToDevice));
  CK(cudaMemcpy(g->dsb,sb,n*4,cudaMemcpyHostToDevice));
  int TPB=64,bl=(n+TPB-1)/TPB;
  disc_batch<<<bl,TPB>>>(g->dA,g->alen,g->dB,g->blen,n,g->dsa,g->dsb,g->dab,g->dae,g->dbb,g->dbe,g->ddf);
  if(cudaGetLastError()!=cudaSuccess||cudaDeviceSynchronize()!=cudaSuccess) return 1;
  CK(cudaMemcpy(ab,g->dab,n*4,cudaMemcpyDeviceToHost));CK(cudaMemcpy(ae,g->dae,n*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(bb,g->dbb,n*4,cudaMemcpyDeviceToHost));CK(cudaMemcpy(be,g->dbe,n*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(diffs,g->ddf,n*4,cudaMemcpyDeviceToHost));
  return 0;
}

#ifdef GPU_LIB_TEST
// validate the library reproduces discovery on the disc tasks (each task's window
// treated as a contig pair; seed = task seed) — same metric as gpu_discover.
#include <vector>
#include <algorithm>
#include "disc_format.h"
int main(int argc,char**argv){
  if(argc<2){fprintf(stderr,"usage: %s <disc-file>\n",argv[0]);return 1;}
  FILE*f=fopen(argv[1],"rb"); FGADiscHeader h; if(fread(&h,sizeof(h),1,f)!=1||h.magic!=FGA_DISC_MAGIC){fprintf(stderr,"bad\n");return 1;}
  gpu_ctx*g=gpu_open();
  long e50=0,tot=0; int N=h.ntasks;
  for(int i=0;i<N;i++){
    uint32_t aw,bw,sa,sb,rd; int rab,rae,rbb,rbe;
    fread(&aw,4,1,f);fread(&bw,4,1,f);fread(&sa,4,1,f);fread(&sb,4,1,f);
    fread(&rab,4,1,f);fread(&rae,4,1,f);fread(&rbb,4,1,f);fread(&rbe,4,1,f);fread(&rd,4,1,f);
    std::vector<u8> A(aw),B(bw); if(aw)fread(A.data(),1,aw,f); if(bw)fread(B.data(),1,bw,f);
    int isa=sa,isb=sb,ab,ae,bb,be,df;
    gpu_load_seqs(g,A.data(),aw,B.data(),bw);
    gpu_discover_batch(g,1,&isa,&isb,&ab,&ae,&bb,&be,&df);
    int mx=std::max(std::max(abs(ab-rab),abs(ae-rae)),std::max(abs(bb-rbb),abs(be-rbe)));
    if(mx<=50)e50++; tot++;
  }
  gpu_close(g); fclose(f);
  printf("fastga_gpu library self-test: %ld/%ld endpoints within 50bp (%.1f%%)\n",e50,tot,100.0*e50/tot);
  return 0;
}
#endif
