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

// ---- trace-point emission (A2): one warp per alignment, banded wave (store all d) +
//      lane-0 traceback binning (diff, delta-b) into global-tspace-phased panels.
//      Validated offline by gpu/trace_validate.cu: 100% Check_Trace_Points-valid, exact
//      edit distance vs FastGA. ----
#define TR_KBAND 384
#define TR_WCAP  1024        // max band width; wider -> overflow (-1)
#define TR_DCAP  2048        // max diffs stored; deeper -> overflow (-2)
#define TR_PCAP  256         // max panels/alignment in the output slot; more -> overflow (-3)

__global__ void trace_batch(const u8*A,const u8*B,int n,
                            const int*abA,const int*aeA,const int*bbA,const int*beA,
                            int tspace, short*Wave, unsigned short*Out, int*Otlen, int max_pairs){
  int warp = (blockIdx.x*blockDim.x + threadIdx.x) >> 5;
  int lane = threadIdx.x & 31;
  if (warp >= n) return;
  int abpos=abA[warp], aepos=aeA[warp], bbpos=bbA[warp], bepos=beA[warp];
  const u8*a = A + abpos; int aw = aepos-abpos;
  const u8*b = B + bbpos; int bw = bepos-bbpos;
  int drift = aw-bw;
  int klo=(drift<0?drift:0)-TR_KBAND, khi=(drift>0?drift:0)+TR_KBAND;
  int width = khi-klo+1;
  short *F = Wave + (size_t)warp*TR_DCAP*TR_WCAP;
  if (width>TR_WCAP){ if(lane==0)Otlen[warp]=-1; return; }

  for (int i=lane;i<width;i+=32) F[i]=-1;
  __syncwarp();
  if (lane==0){ int x=0; while(x<aw&&x<bw&&a[x]==b[x])x++; F[0-klo]=(short)x; }
  __syncwarp();
  int result=-1;
  if (drift==0 && F[0-klo]>=aw) result=0;
  for (int d=1; d<=TR_DCAP-1 && result<0; d++){
    short *fp=F+(size_t)(d-1)*TR_WCAP, *fc=F+(size_t)d*TR_WCAP;
    for (int i=lane;i<width;i+=32) fc[i]=-1;
    __syncwarp();
    int kmin=klo>-d?klo:-d, kmax=khi<d?khi:d;
    for (int k=kmin+lane;k<=kmax;k+=32){
      int idx=k-klo,best=-1,v;
      if(k-1>=klo){v=fp[idx-1];if(v>=0&&v+1>best)best=v+1;}
      if(k+1<=khi){v=fp[idx+1];if(v>=0&&v  >best)best=v;}
      v=fp[idx];if(v>=0&&v+1>best)best=v+1;
      if(best>=0){int x=best;if(x>aw)x=aw;if(x-k>bw)x=bw+k;if(x<0)x=0;int y=x-k;
        while(x<aw&&y<bw&&a[x]==b[y]){x++;y++;} fc[idx]=(short)x;}
    }
    __syncwarp();
    int xe=fc[drift-klo];
    if(xe>=aw&&xe-drift>=bw) result=d;
  }
  if (lane!=0) return;
  if (result<0){ Otlen[warp]=-2; return; }

  int p0=abpos/tspace;
  int npanel=(aepos-1)/tspace - p0 + 1;
  if (npanel>TR_PCAP || 2*npanel>max_pairs){ Otlen[warp]=-3; return; }
  unsigned short *o = Out + (size_t)warp*max_pairs;
  for (int p=0;p<2*npanel;p++) o[p]=0;
  #define TR_DBIN(al) o[2*(((abpos+(al))/tspace)-p0)]
  #define TR_BBIN(al) o[2*(((abpos+(al))/tspace)-p0)+1]
  int d=result, k=drift;
  int x=(int)F[(size_t)d*TR_WCAP+(k-klo)];
  while (d>0){
    short *fp=F+(size_t)(d-1)*TR_WCAP;
    int best=-1,op=0,xpred=-1,v;
    v=fp[k-klo];  if(v>=0&&v+1>best){best=v+1;op=2;xpred=v;}                  // sub
    if(k-1>=klo){v=fp[(k-1)-klo];if(v>=0&&v+1>best){best=v+1;op=0;xpred=v;}}  // del
    if(k+1<=khi){v=fp[(k+1)-klo];if(v>=0&&v  >best){best=v;  op=1;xpred=v;}}  // ins
    int x0=best;
    for (int al=x0; al<x; al++) TR_BBIN(al) += 1;      // match slide: 1 B each
    if (op==0){ TR_DBIN(xpred)+=1; }                   // del: +1 diff, 0 B
    else      { TR_DBIN(xpred)+=1; TR_BBIN(xpred)+=1; }// sub/ins: +1 diff, +1 B
    if (op==0) k=k-1; else if (op==1) k=k+1;
    x=xpred; d=d-1;
  }
  for (int al=0; al<x; al++) TR_BBIN(al) += 1;         // d=0 slide from (0,0)
  #undef TR_DBIN
  #undef TR_BBIN
  Otlen[warp]=2*npanel;
}

// ---- warp-parallel discovery (A3b): one warp per task, x-drop band across 32 lanes.
//      Same furthest-reaching wave + score=matches-2*diffs x-drop as extend(), but the band
//      is swept by the warp and the per-d best score is a warp argmax-reduction. ----
__device__ static void extend_warp(const u8*A,int a0,int da,int alen,const u8*B,int b0,int db,int blen,
                                   short*f,short*M,int lane,int*ox,int*oy,int*od){
  short *fp=f, *fc=f+DW, *Mp=M, *Mc=M+DW;
  for(int i=lane;i<DW;i+=32){fp[i]=-1;Mp[i]=0;fc[i]=-1;Mc[i]=0;}
  __syncwarp();
  int x0=0,mm0=0;
  if(lane==0){int x=0,mm=0;while(x<alen&&x<blen&&A[a0+da*x]==B[b0+db*x]){x++;mm++;}
              fp[KB]=(short)x;Mp[KB]=(short)mm;x0=x;mm0=mm;}
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
      int xx=best;if(xx>alen)xx=alen;if(xx-k>blen)xx=blen+k;if(xx<0)xx=0;int yy=xx-k;
      while(xx<alen&&yy<blen&&A[a0+da*xx]==B[b0+db*yy]){xx++;yy++;bm++;}
      fc[idx]=(short)xx;Mc[idx]=(short)bm;
      int S=bm-2*d; if(S>locS){locS=S;lcx=xx;lcy=yy;}
    }
    __syncwarp();
    int mS=locS,mx=lcx,my=lcy;               // warp argmax of the per-d best score
    for(int o=16;o>=1;o>>=1){
      int oS=__shfl_down_sync(0xffffffff,mS,o),oX=__shfl_down_sync(0xffffffff,mx,o),oY=__shfl_down_sync(0xffffffff,my,o);
      if(oS>mS){mS=oS;mx=oX;my=oY;}
    }
    int curbestS=__shfl_sync(0xffffffff,mS,0);
    int cx=__shfl_sync(0xffffffff,mx,0), cy=__shfl_sync(0xffffffff,my,0);
    if(lane==0 && curbestS>bestS){bestS=curbestS;bx=cx;by=cy;bd=d;}
    bestS=__shfl_sync(0xffffffff,bestS,0);
    if(curbestS<bestS-40) stop=1;
    short*t; t=fp;fp=fc;fc=t; t=Mp;Mp=Mc;Mc=t;
  }
  if(lane==0){*ox=bx;*oy=by;*od=bd;}
}
__global__ void disc_batch_warp(const u8*A,int alen,const u8*B,int blen,int n,
                                const int*sa,const int*sb,int*ab,int*ae,int*bb,int*be,int*df){
  extern __shared__ short sh[];
  int warp=(blockIdx.x*blockDim.x+threadIdx.x)>>5, lane=threadIdx.x&31, wib=threadIdx.x>>5;
  if(warp>=n)return;
  short*f=sh+(size_t)wib*(4*DW); short*M=f+2*DW;
  int SA=sa[warp],SB=sb[warp],fx,fy,fd,rx,ry,rd;
  extend_warp(A,SA,  +1,alen-SA,B,SB,  +1,blen-SB,f,M,lane,&fx,&fy,&fd);
  extend_warp(A,SA-1,-1,SA,     B,SB-1,-1,SB,     f,M,lane,&rx,&ry,&rd);
  if(lane==0){ae[warp]=SA+fx;be[warp]=SB+fy;ab[warp]=SA-rx;bb[warp]=SB-ry;df[warp]=fd+rd;}
}

// on-device 2-bit -> NUMERIC unpack: base[p] = (packed[p>>2] >> ((p&3)*2)) & 3
// (byte-identical to gene_core.c Uncompress_Read). Eliminates the CPU decompression that
// perf shows is 29% of the human align run, and shrinks the host->device copy 4x.
__global__ void unpack2bit(const u8* packed, u8* out, int len){
  int p = blockIdx.x*blockDim.x + threadIdx.x;
  if(p>=len) return;
  out[p] = (packed[p>>2] >> ((p&3)*2)) & 3;
}

struct gpu_ctx { u8 *dA,*dB; int alen,blen; size_t acap,bcap;
                 u8 *dPackA,*dPackB; size_t pcapA,pcapB;   // device 2-bit scratch
                 int *dsa,*dsb,*dab,*dae,*dbb,*dbe,*ddf; int ncap;
                 // trace-emission scratch/buffers
                 int *tab,*tae,*tbb,*tbe,*tOtlen; short *tWave; unsigned short *tOut;
                 int tncap, tchunkcap; };

extern "C" gpu_ctx *gpu_open(void){
  gpu_ctx *g=(gpu_ctx*)calloc(1,sizeof(gpu_ctx)); return g;
}
extern "C" void gpu_close(gpu_ctx *g){
  if(!g)return;
  if(g->dA)cudaFree(g->dA); if(g->dB)cudaFree(g->dB);
  if(g->dPackA)cudaFree(g->dPackA); if(g->dPackB)cudaFree(g->dPackB);
  if(g->dsa){cudaFree(g->dsa);cudaFree(g->dsb);cudaFree(g->dab);cudaFree(g->dae);cudaFree(g->dbb);cudaFree(g->dbe);cudaFree(g->ddf);}
  if(g->tab){cudaFree(g->tab);cudaFree(g->tae);cudaFree(g->tbb);cudaFree(g->tbe);cudaFree(g->tOtlen);cudaFree(g->tOut);}
  if(g->tWave)cudaFree(g->tWave);
  free(g);
}
extern "C" void gpu_load_seqs(gpu_ctx *g,const u8*A,int alen,const u8*B,int blen){
  if((size_t)alen>g->acap){ if(g->dA)cudaFree(g->dA); CK(cudaMalloc(&g->dA,alen)); g->acap=alen; }
  if((size_t)blen>g->bcap){ if(g->dB)cudaFree(g->dB); CK(cudaMalloc(&g->dB,blen)); g->bcap=blen; }
  CK(cudaMemcpy(g->dA,A,alen,cudaMemcpyHostToDevice));
  CK(cudaMemcpy(g->dB,B,blen,cudaMemcpyHostToDevice));
  g->alen=alen; g->blen=blen;
}
// Load 2-bit packed contigs (from .bps, NO CPU decompression) and unpack on-device.
// packedA/B hold (alen+3)/4 and (blen+3)/4 bytes. Resident dA/dB end up NUMERIC (0-3),
// identical to gpu_load_seqs' input -- but the CPU never runs Uncompress_Read.
extern "C" void gpu_load_seqs_2bit(gpu_ctx *g,const u8*packedA,int alen,const u8*packedB,int blen){
  size_t pa=(alen+3)/4, pb=(blen+3)/4;
  if((size_t)alen>g->acap){ if(g->dA)cudaFree(g->dA); CK(cudaMalloc(&g->dA,alen)); g->acap=alen; }
  if((size_t)blen>g->bcap){ if(g->dB)cudaFree(g->dB); CK(cudaMalloc(&g->dB,blen)); g->bcap=blen; }
  if(pa>g->pcapA){ if(g->dPackA)cudaFree(g->dPackA); CK(cudaMalloc(&g->dPackA,pa)); g->pcapA=pa; }
  if(pb>g->pcapB){ if(g->dPackB)cudaFree(g->dPackB); CK(cudaMalloc(&g->dPackB,pb)); g->pcapB=pb; }
  CK(cudaMemcpy(g->dPackA,packedA,pa,cudaMemcpyHostToDevice));
  CK(cudaMemcpy(g->dPackB,packedB,pb,cudaMemcpyHostToDevice));
  int TPB=256;
  unpack2bit<<<(int)((alen+TPB-1)/TPB),TPB>>>(g->dPackA,g->dA,alen);
  unpack2bit<<<(int)((blen+TPB-1)/TPB),TPB>>>(g->dPackB,g->dB,blen);
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
  int TPB=128, wpb=TPB/32, bl=(n+wpb-1)/wpb;
  size_t shmem=(size_t)wpb*4*DW*sizeof(short);
  disc_batch_warp<<<bl,TPB,shmem>>>(g->dA,g->alen,g->dB,g->blen,n,g->dsa,g->dsb,g->dab,g->dae,g->dbb,g->dbe,g->ddf);
  if(cudaGetLastError()!=cudaSuccess||cudaDeviceSynchronize()!=cudaSuccess) return 1;
  CK(cudaMemcpy(ab,g->dab,n*4,cudaMemcpyDeviceToHost));CK(cudaMemcpy(ae,g->dae,n*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(bb,g->dbb,n*4,cudaMemcpyDeviceToHost));CK(cudaMemcpy(be,g->dbe,n*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(diffs,g->ddf,n*4,cudaMemcpyDeviceToHost));
  return 0;
}

#define TR_CHUNK 2048     // max warps per trace launch; wave scratch = min(n,TR_CHUNK)*4 MB.
                          // Throughput is occupancy-bound: 256->24k, 2048->114k, 8192->306k aln/s
                          // on the A100 (bench). Raise for big batches if GPU memory allows.
static void ensure_trace(gpu_ctx*g,int n){
  if(n>g->tncap){
    if(g->tab){cudaFree(g->tab);cudaFree(g->tae);cudaFree(g->tbb);cudaFree(g->tbe);cudaFree(g->tOtlen);cudaFree(g->tOut);}
    CK(cudaMalloc(&g->tab,n*4));CK(cudaMalloc(&g->tae,n*4));CK(cudaMalloc(&g->tbb,n*4));CK(cudaMalloc(&g->tbe,n*4));
    CK(cudaMalloc(&g->tOtlen,n*4));
    CK(cudaMalloc(&g->tOut,(size_t)n*FGA_TRACE_MAX_PAIRS*sizeof(unsigned short)));
    g->tncap=n;
  }
  int chunk = n<TR_CHUNK ? n : TR_CHUNK;     // don't over-allocate wave scratch for small n
  if(chunk>g->tchunkcap){
    if(g->tWave)cudaFree(g->tWave);
    CK(cudaMalloc(&g->tWave,(size_t)chunk*TR_DCAP*TR_WCAP*sizeof(short)));
    g->tchunkcap=chunk;
  }
}
extern "C" int gpu_trace_batch(gpu_ctx*g,int n,
        const int*ab,const int*ae,const int*bb,const int*be,int tspace,
        unsigned short*out_trace,int*out_tlen){
  if(n<=0) return 0;
  ensure_trace(g,n);
  CK(cudaMemcpy(g->tab,ab,n*4,cudaMemcpyHostToDevice));CK(cudaMemcpy(g->tae,ae,n*4,cudaMemcpyHostToDevice));
  CK(cudaMemcpy(g->tbb,bb,n*4,cudaMemcpyHostToDevice));CK(cudaMemcpy(g->tbe,be,n*4,cudaMemcpyHostToDevice));
  for(int c=0;c<n;c+=TR_CHUNK){
    int m=(n-c<TR_CHUNK)?n-c:TR_CHUNK;
    int tpb=128, blk=(m*32+tpb-1)/tpb;
    trace_batch<<<blk,tpb>>>(g->dA,g->dB,m,g->tab+c,g->tae+c,g->tbb+c,g->tbe+c,
                             tspace,g->tWave,g->tOut+(size_t)c*FGA_TRACE_MAX_PAIRS,g->tOtlen+c,
                             FGA_TRACE_MAX_PAIRS);
  }
  if(cudaGetLastError()!=cudaSuccess||cudaDeviceSynchronize()!=cudaSuccess) return 1;
  CK(cudaMemcpy(out_tlen,g->tOtlen,n*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(out_trace,g->tOut,(size_t)n*FGA_TRACE_MAX_PAIRS*sizeof(unsigned short),cudaMemcpyDeviceToHost));
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
