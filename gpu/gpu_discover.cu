/* FastGA GPU local-alignment DISCOVERY kernel (toward G4 integration).
 *
 * Unlike gpu_align (edit distance of a pre-cut segment), this reproduces the core
 * of FastGA's Local_Alignment: from a seed anchor inside a window, extend forward
 * and backward with an x-drop furthest-reaching wave to DISCOVER the alignment
 * endpoints (abpos,aepos,bbpos,bepos) + diffs.  Validated against FastGA's endpoints
 * from the .1aln (via extract_disc).  Correctness-first: one thread per task.
 *
 * Build: nvcc -O3 -arch=sm_80 -Xcompiler -fopenmp -o gpu_discover gpu_discover.cu
 */
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <chrono>
#include <algorithm>
#include <omp.h>
#include <cuda_runtime.h>
#include "disc_format.h"

typedef unsigned char u8;
#define KB    320       // half band; M2: band<256
#define DW    (2*KB+1)
#define MAXD  8192
#define XDROP 40        // stop extending when frontier falls XDROP anti-diagonals behind best

#define CK(c) do{cudaError_t e=(c); if(e!=cudaSuccess){fprintf(stderr,"CUDA %s:%d %s\n",__FILE__,__LINE__,cudaGetErrorString(e));exit(1);}}while(0)

// Forward x-drop furthest-reaching (Levenshtein) extension from the origin along
// A[a0 + da*i] and B[b0 + db*i] (i>=0).  Returns the best (furthest anti-diagonal)
// reach as offsets (*ox,*oy) from the origin and its diff count *od.
__host__ __device__ static void extend(const u8* A,int a0,int da,int alen,
                                       const u8* B,int b0,int db,int blen,
                                       int* ox,int* oy,int* od)
{
  short f[2][DW], M[2][DW];              // f = furthest x; M = matches along that path
  int cur=0, prv=1;
  for (int i=0;i<DW;i++){ f[prv][i]=-1; M[prv][i]=0; }
  int x=0, mm=0; while (x<alen && x<blen && A[a0+da*x]==B[b0+db*x]){x++;mm++;}
  for (int i=0;i<DW;i++){ f[cur][i]=-1; M[cur][i]=0; }
  f[cur][0+KB]=(short)x; M[cur][0+KB]=mm;
  int bestS=mm, bx=x, by=x, bd=0;      // score = matches - diffs
  { int t=cur; cur=prv; prv=t; }
  for (int d=1; d<=MAXD; d++)
    { for (int i=0;i<DW;i++) f[cur][i]=-1;
      int kmin=-KB>-d?-KB:-d, kmax=KB<d?KB:d, curbestS=-1000000, cx=0, cy=0;
      for (int k=kmin;k<=kmax;k++)
        { int best=-1, bm=0, v;
          if (k-1>=-KB){ v=f[prv][k-1+KB]; if(v>=0&&v+1>best){best=v+1; bm=M[prv][k-1+KB];} } // horizontal
          if (k+1<= KB){ v=f[prv][k+1+KB]; if(v>=0&&v  >best){best=v;   bm=M[prv][k+1+KB];} } // vertical
          v=f[prv][k+KB]; if(v>=0&&v+1>best){best=v+1; bm=M[prv][k+KB];}                      // substitution
          if (best<0) continue;
          int xx=best; if(xx>alen)xx=alen; if(xx-k>blen)xx=blen+k; if(xx<0)xx=0;
          int yy=xx-k;
          while (xx<alen && yy<blen && A[a0+da*xx]==B[b0+db*yy]){xx++;yy++;bm++;}
          f[cur][k+KB]=(short)xx; M[cur][k+KB]=bm;
          int S = bm - 2*d;                    // local score at this point
          if (S>curbestS){ curbestS=S; cx=xx; cy=yy; }
        }
      if (curbestS > bestS){ bestS=curbestS; bx=cx; by=cy; bd=d; }
      if (curbestS < bestS - XDROP) break;   // x-drop on matches-minus-diffs score
      int t=cur; cur=prv; prv=t;
    }
  *ox=bx; *oy=by; *od=bd;
}

// discover local alignment: extend forward + backward from seed (sa,sb) in window
__host__ __device__ static void discover(const u8* A,int aw,const u8* B,int bw,
                                         int sa,int sb,
                                         int* ab,int* ae,int* bb,int* be,int* diffs)
{
  int fx,fy,fd, rx,ry,rd;
  extend(A, sa, +1, aw-sa, B, sb, +1, bw-sb, &fx,&fy,&fd);   // forward from seed
  extend(A, sa-1,-1, sa,    B, sb-1,-1, sb,    &rx,&ry,&rd);   // backward from seed
  *ae = sa+fx; *be = sb+fy; *ab = sa-rx; *bb = sb-ry; *diffs = fd+rd;
}

struct Disc { int n; std::vector<u8> A,B; std::vector<int> Ao,Bo,Aw,Bw,Sa,Sb,Rab,Rae,Rbb,Rbe,Rd; };

static void load(const char* p, Disc& D){
  FILE* f=fopen(p,"rb"); if(!f){fprintf(stderr,"open %s\n",p);exit(1);}
  FGADiscHeader h; if(fread(&h,sizeof(h),1,f)!=1||h.magic!=FGA_DISC_MAGIC){fprintf(stderr,"bad magic\n");exit(1);}
  D.n=h.ntasks; printf("disc tasks: %u  margin: %u\n",h.ntasks,h.margin);
  for(uint32_t i=0;i<h.ntasks;i++){
    uint32_t aw,bw,sa,sb,rd; int rab,rae,rbb,rbe;
    fread(&aw,4,1,f);fread(&bw,4,1,f);fread(&sa,4,1,f);fread(&sb,4,1,f);
    fread(&rab,4,1,f);fread(&rae,4,1,f);fread(&rbb,4,1,f);fread(&rbe,4,1,f);fread(&rd,4,1,f);
    D.Ao.push_back(D.A.size()); D.Aw.push_back(aw);
    D.Bo.push_back(D.B.size()); D.Bw.push_back(bw);
    D.Sa.push_back(sa);D.Sb.push_back(sb);D.Rab.push_back(rab);D.Rae.push_back(rae);
    D.Rbb.push_back(rbb);D.Rbe.push_back(rbe);D.Rd.push_back(rd);
    size_t a=D.A.size();D.A.resize(a+aw); if(aw&&fread(&D.A[a],1,aw,f)!=aw){fprintf(stderr,"truncA\n");exit(1);}
    size_t b=D.B.size();D.B.resize(b+bw); if(bw&&fread(&D.B[b],1,bw,f)!=bw){fprintf(stderr,"truncB\n");exit(1);}
  }
  fclose(f);
}

__global__ void disc_kernel(const u8*A,const u8*B,const int*Ao,const int*Bo,const int*Aw,const int*Bw,
                            const int*Sa,const int*Sb,int*oab,int*oae,int*obb,int*obe,int*od,int n){
  int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=n)return;
  discover(A+Ao[i],Aw[i],B+Bo[i],Bw[i],Sa[i],Sb[i],&oab[i],&oae[i],&obb[i],&obe[i],&od[i]);
}

int main(int argc,char**argv){
  // self-test: embed a known ~500bp alignment (5% diverged) in random flanks; discover it.
  { srand(7); int bad=0,tot=200;
    for(int t=0;t<tot;t++){
      int core=300+rand()%400, fl=100+rand()%150, aw=fl+core+fl, bw=fl+core+fl;
      std::vector<u8> A(aw),B(bw); for(auto&x:A)x=rand()&3; for(auto&x:B)x=rand()&3;
      // shared core with 5% substitutions at A[fl..fl+core), B[fl..fl+core)
      for(int i=0;i<core;i++){ u8 c=rand()&3; A[fl+i]=c; B[fl+i]=(rand()%100<5)?(u8)((c+1)&3):c; }
      int ab,ae,bb,be,df; discover(A.data(),aw,B.data(),bw, fl+core/2, fl+core/2, &ab,&ae,&bb,&be,&df);
      // endpoints should be near [fl, fl+core]
      if(abs(ab-fl)>30||abs(ae-(fl+core))>30){ if(bad<4)fprintf(stderr,"selftest t=%d core=%d got A[%d,%d] want~[%d,%d] df=%d\n",t,core,ab,ae,fl,fl+core,df); bad++; } }
    if(bad>tot/20){ fprintf(stderr,"DISC SELF-TEST: %d/%d off\n",bad,tot); }
    else printf("disc self-test: embedded-alignment endpoints recovered (%d/%d within 30bp)  OK\n",tot-bad,tot); }

  if(argc<2){ fprintf(stderr,"usage: %s <disc-file>\n",argv[0]); return 1; }
  Disc D; load(argv[1],D); int N=D.n;
  std::vector<int> gab(N),gae(N),gbb(N),gbe(N),gd(N);

  u8*dA,*dB; int *dAo,*dBo,*dAw,*dBw,*dSa,*dSb,*doab,*doae,*dobb,*dobe,*dod;
  CK(cudaMalloc(&dA,D.A.size()));CK(cudaMalloc(&dB,D.B.size()));
  CK(cudaMalloc(&dAo,N*4));CK(cudaMalloc(&dBo,N*4));CK(cudaMalloc(&dAw,N*4));CK(cudaMalloc(&dBw,N*4));
  CK(cudaMalloc(&dSa,N*4));CK(cudaMalloc(&dSb,N*4));
  CK(cudaMalloc(&doab,N*4));CK(cudaMalloc(&doae,N*4));CK(cudaMalloc(&dobb,N*4));CK(cudaMalloc(&dobe,N*4));CK(cudaMalloc(&dod,N*4));
  CK(cudaMemcpy(dA,D.A.data(),D.A.size(),cudaMemcpyHostToDevice));CK(cudaMemcpy(dB,D.B.data(),D.B.size(),cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dAo,D.Ao.data(),N*4,cudaMemcpyHostToDevice));CK(cudaMemcpy(dBo,D.Bo.data(),N*4,cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dAw,D.Aw.data(),N*4,cudaMemcpyHostToDevice));CK(cudaMemcpy(dBw,D.Bw.data(),N*4,cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dSa,D.Sa.data(),N*4,cudaMemcpyHostToDevice));CK(cudaMemcpy(dSb,D.Sb.data(),N*4,cudaMemcpyHostToDevice));
  int TPB=64,bl=(N+TPB-1)/TPB;
  auto t0=std::chrono::high_resolution_clock::now();
  disc_kernel<<<bl,TPB>>>(dA,dB,dAo,dBo,dAw,dBw,dSa,dSb,doab,doae,dobb,dobe,dod,N);
  CK(cudaGetLastError());CK(cudaDeviceSynchronize());
  auto t1=std::chrono::high_resolution_clock::now();
  CK(cudaMemcpy(gab.data(),doab,N*4,cudaMemcpyDeviceToHost));CK(cudaMemcpy(gae.data(),doae,N*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(gbb.data(),dobb,N*4,cudaMemcpyDeviceToHost));CK(cudaMemcpy(gbe.data(),dobe,N*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(gd.data(),dod,N*4,cudaMemcpyDeviceToHost));
  double ker=std::chrono::duration<double>(t1-t0).count();

  // validate endpoints vs FastGA reference (within tolerance) + coverage overlap
  long e10=0,e50=0; double covsum=0; long lediff=0;
  for(int i=0;i<N;i++){
    int dab=abs(gab[i]-D.Rab[i]),dae=abs(gae[i]-D.Rae[i]),dbb=abs(gbb[i]-D.Rbb[i]),dbe=abs(gbe[i]-D.Rbe[i]);
    int mx=std::max(std::max(dab,dae),std::max(dbb,dbe));
    if(mx<=10)e10++; if(mx<=50)e50++;
    // A-interval overlap fraction with reference
    int ib=std::max(gab[i],D.Rab[i]), ie=std::min(gae[i],D.Rae[i]);
    int refl=D.Rae[i]-D.Rab[i]; double ov = refl>0 ? (double)std::max(0,ie-ib)/refl : 1.0; covsum+=ov;
    if(gd[i]<=D.Rd[i]) lediff++;
  }
  printf("\n=== discovery correctness (vs FastGA endpoints) ===\n");
  printf("endpoints within 10bp (all 4)  : %.3f%%\n",100.0*e10/N);
  printf("endpoints within 50bp (all 4)  : %.3f%%\n",100.0*e50/N);
  printf("mean A-interval overlap w/ ref  : %.3f%%\n",100.0*covsum/N);
  printf("GPU diffs <= FastGA diffs       : %.3f%%\n",100.0*lediff/N);
  printf("\n=== perf ===\nGPU discovery kernel: %.3f s  (%.3f M aln/s)\n",ker,N/ker/1e6);
  return 0;
}
