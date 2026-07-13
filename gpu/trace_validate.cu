/* A1 validation: GPU emits FastGA trace-points, compared to Compute_Alignment's reference.
 *
 * For each task (gpu/trace_format.h): align the fixed rectangle A[0..aw] x B[0..bw]
 * corner-to-corner with a banded furthest-reaching O(ND) wave, STORING every d-layer, then
 * traceback (lane 0) recovering match/sub/ins/del ops and binning (diff, delta-b) into
 * trace-point panels phased on global multiples of tspace (panel = floor((abpos+a)/ts)).
 * Emit the uint16 pair vector and compare to ref_trace pair-by-pair.
 *
 * Build:  make trace_validate      Run: ./gpu/trace_validate <tasks.trace>
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <cuda_runtime.h>
#include "trace_format.h"

typedef unsigned char u8;
#define KBAND 384
#define WCAP  1024          // max band width; task wider -> overflow
#define DCAP  2048          // max diffs stored; task deeper -> overflow
#define PCAP  256           // max panels per task in the output buffer

#define CK(c) do{cudaError_t e=(c); if(e!=cudaSuccess){fprintf(stderr,"CUDA %s:%d %s\n",__FILE__,__LINE__,cudaGetErrorString(e));exit(1);}}while(0)

// ---- one warp per task: banded wave (store all d) + lane-0 traceback + panel binning ----
__global__ void trace_kernel(const u8 *Abuf, const u8 *Bbuf,
                             const int *Aoff, const int *Boff,
                             const int *AW, const int *BW, const int *ABP,
                             int tspace, short *Wave,   // Wave: warp-local slice DCAP*WCAP
                             unsigned short *Out, int *Otlen, int *Odiffs, int ntasks)
{
  int warp = (blockIdx.x*blockDim.x + threadIdx.x) >> 5;
  int lane = threadIdx.x & 31;
  if (warp >= ntasks) return;

  const u8 *A = Abuf + Aoff[warp]; int aw = AW[warp];
  const u8 *B = Bbuf + Boff[warp]; int bw = BW[warp];
  int abpos = ABP[warp];
  int drift = aw - bw;
  int klo = (drift<0?drift:0) - KBAND;
  int khi = (drift>0?drift:0) + KBAND;
  int width = khi - klo + 1;
  short *F = Wave + (size_t)warp * DCAP * WCAP;   // F[d*WCAP + (k-klo)]

  if (width > WCAP) { if(lane==0){ Otlen[warp]=-1; } return; }   // overflow: band

  // d = 0 : only k=0 reachable, slide matches from (0,0)
  for (int i = lane; i < width; i += 32) F[i] = -1;
  __syncwarp();
  if (lane == 0) { int x=0; while (x<aw && x<bw && A[x]==B[x]) x++; F[0-klo] = (short)x; }
  __syncwarp();

  int result = -1;
  {
    int reached0 = (drift==0 && F[0-klo] >= aw);
    if (reached0) result = 0;
  }
  for (int d = 1; d <= DCAP-1 && result < 0; d++)
    { short *fp = F + (size_t)(d-1)*WCAP;
      short *fc = F + (size_t)d*WCAP;
      for (int i = lane; i < width; i += 32) fc[i] = -1;
      __syncwarp();
      int kmin = klo>-d?klo:-d, kmax = khi<d?khi:d;
      for (int k = kmin+lane; k <= kmax; k += 32)
        { int idx=k-klo, best=-1, v;
          if (k-1>=klo){ v=fp[idx-1]; if(v>=0&&v+1>best)best=v+1; }   // horizontal (del)
          if (k+1<=khi){ v=fp[idx+1]; if(v>=0&&v  >best)best=v;   }   // vertical  (ins)
          v=fp[idx];     if(v>=0&&v+1>best)best=v+1;                  // substitution
          if (best>=0){ int x=best; if(x>aw)x=aw; if(x-k>bw)x=bw+k; if(x<0)x=0; int y=x-k;
            while (x<aw && y<bw && A[x]==B[y]){x++;y++;} fc[idx]=(short)x; }
        }
      __syncwarp();
      int xe = fc[drift-klo];
      if (xe>=aw && xe-drift>=bw) result = d;
    }

  if (lane != 0) return;
  if (result < 0) { Otlen[warp] = -2; return; }   // overflow: depth

  // ---- panel bins (global multiples of tspace) ----
  int aepos = abpos + aw;
  int p0    = abpos / tspace;
  int npanel = (aepos-1)/tspace - p0 + 1;         // == Check_Trace_Points count
  if (npanel > PCAP) { Otlen[warp] = -3; return; }
  unsigned short *o = Out + (size_t)warp * 2 * PCAP;
  for (int p = 0; p < 2*npanel; p++) o[p] = 0;
  #define DBIN(aloc) o[2*(((abpos+(aloc))/tspace) - p0)]
  #define BBIN(aloc) o[2*(((abpos+(aloc))/tspace) - p0) + 1]

  // ---- traceback from (aw,bw) i.e. (d=result, k=drift) back to (0,0) ----
  int d = result, k = drift;
  int x = (int) F[(size_t)d*WCAP + (k-klo)];       // == aw
  while (d > 0)
    { short *fp = F + (size_t)(d-1)*WCAP;
      // pick predecessor achieving the (canonical) furthest-reaching max; tie order tuned
      // to match FastGA's leftmost canonicalization.  SUB first (diagonal), then DEL, then INS.
      int best=-1, op=0, xpred=-1;   // op: 0=DEL(horiz), 1=INS(vert), 2=SUB
      int v;
      v=fp[k-klo];  if(v>=0&&v+1>best){best=v+1;op=2;xpred=v;}                          // sub
      if (k-1>=klo){ v=fp[(k-1)-klo]; if(v>=0&&v+1>best){best=v+1;op=0;xpred=v;} }       // del
      if (k+1<=khi){ v=fp[(k+1)-klo]; if(v>=0&&v  >best){best=v;  op=1;xpred=v;} }       // ins
      int x0 = best;                 // A-coord after the op, before the match slide
      // match slide [x0, x): each match consumes 1 A and 1 B at A-coord a
      for (int a = x0; a < x; a++) { BBIN(a) += 1; }
      // the single diff op at A-coord xpred
      if (op == 0)        { DBIN(xpred) += 1; }                 // DEL: +1 diff, no B
      else if (op == 1)   { DBIN(xpred) += 1; BBIN(xpred) += 1; } // INS: +1 diff, +1 B
      else                { DBIN(xpred) += 1; BBIN(xpred) += 1; } // SUB: +1 diff, +1 B
      // step to predecessor
      if (op == 0) k = k-1;
      else if (op == 1) k = k+1;
      x = xpred;
      d = d-1;
    }
  // d==0: slide matches from (0,0) to (x,x) -- consume B along the way
  for (int a = 0; a < x; a++) { BBIN(a) += 1; }

  int diffs = 0, bsum = 0;
  for (int p = 0; p < npanel; p++) { diffs += o[2*p]; bsum += o[2*p+1]; }
  Otlen[warp]  = 2*npanel;
  Odiffs[warp] = diffs;
  (void) bsum;
}

int main(int argc, char **argv)
{ if (argc < 2) { fprintf(stderr,"usage: %s <tasks.trace>\n",argv[0]); return 1; }
  FILE *f = fopen(argv[1],"rb");
  if (!f) { fprintf(stderr,"cannot open %s\n",argv[1]); return 1; }
  FGATraceHeader h;
  if (fread(&h,sizeof(h),1,f)!=1 || h.magic!=FGA_TRACE_MAGIC) { fprintf(stderr,"bad header\n"); return 1; }
  int tspace = (int) h.tspace;
  printf("tasks: %u   tspace: %d\n", h.ntasks, tspace);

  std::vector<u8>  Abuf, Bbuf;
  std::vector<int> Aoff, Boff, AW, BW, ABP;
  std::vector<int> refDiffs, refTlen, refOff;
  std::vector<unsigned short> refTr;
  int skipped_big = 0;
  for (uint32_t i=0;i<h.ntasks;i++)
    { int32_t ab,ae,bb,be; uint32_t aw,bw,rd,rl;
      if (fread(&ab,4,1,f)!=1) break;
      fread(&ae,4,1,f);fread(&bb,4,1,f);fread(&be,4,1,f);
      fread(&aw,4,1,f);fread(&bw,4,1,f);fread(&rd,4,1,f);fread(&rl,4,1,f);
      std::vector<u8> A(aw),B(bw); if(aw)fread(A.data(),1,aw,f); if(bw)fread(B.data(),1,bw,f);
      std::vector<unsigned short> T(rl); if(rl)fread(T.data(),2,rl,f);
      if (aw==0||bw==0){ continue; }
      Aoff.push_back((int)Abuf.size()); AW.push_back((int)aw); Abuf.insert(Abuf.end(),A.begin(),A.end());
      Boff.push_back((int)Bbuf.size()); BW.push_back((int)bw); Bbuf.insert(Bbuf.end(),B.begin(),B.end());
      ABP.push_back(ab);
      refDiffs.push_back((int)rd); refTlen.push_back((int)rl);
      refOff.push_back((int)refTr.size()); refTr.insert(refTr.end(),T.begin(),T.end());
    }
  fclose(f);
  int N = (int)AW.size();
  printf("loaded %d tasks (%d skipped empty)\n", N, skipped_big);

  // device buffers
  u8 *dA,*dB; int *dAo,*dBo,*dAW,*dBW,*dABP,*dOtlen,*dOdiffs; short *dWave; unsigned short *dOut;
  CK(cudaMalloc(&dA,Abuf.size())); CK(cudaMalloc(&dB,Bbuf.size()));
  CK(cudaMemcpy(dA,Abuf.data(),Abuf.size(),cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dB,Bbuf.data(),Bbuf.size(),cudaMemcpyHostToDevice));
  #define MK(p,n) CK(cudaMalloc(&p,(n)*sizeof(int)));
  auto up=[&](int*d,std::vector<int>&v){CK(cudaMemcpy(d,v.data(),v.size()*sizeof(int),cudaMemcpyHostToDevice));};
  MK(dAo,N);MK(dBo,N);MK(dAW,N);MK(dBW,N);MK(dABP,N);MK(dOtlen,N);MK(dOdiffs,N);
  up(dAo,Aoff);up(dBo,Boff);up(dAW,AW);up(dBW,BW);up(dABP,ABP);
  CK(cudaMalloc(&dOut,(size_t)N*2*PCAP*sizeof(unsigned short)));

  // process in chunks to bound the wave scratch (DCAP*WCAP shorts per warp)
  const int CHUNK = 1024;
  CK(cudaMalloc(&dWave,(size_t)CHUNK*DCAP*WCAP*sizeof(short)));
  std::vector<int> Otlen(N), Odiffs(N);
  std::vector<unsigned short> Out((size_t)N*2*PCAP);

  for (int c=0;c<N;c+=CHUNK)
    { int n = (N-c<CHUNK)?N-c:CHUNK;
      int tpb=128, blk=(n*32 + tpb-1)/tpb;
      trace_kernel<<<blk,tpb>>>(dA,dB,dAo+c,dBo+c,dAW+c,dBW+c,dABP+c,tspace,
                                dWave,dOut+(size_t)c*2*PCAP,dOtlen+c,dOdiffs+c,n);
      CK(cudaDeviceSynchronize());
    }
  CK(cudaMemcpy(Otlen.data(),dOtlen,N*sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(Odiffs.data(),dOdiffs,N*sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(Out.data(),dOut,(size_t)N*2*PCAP*sizeof(unsigned short),cudaMemcpyDeviceToHost));

  long exact=0, tlen_ok=0, diff_ok=0, bsum_ok=0, ctp_ok=0, valid=0, ov_band=0, ov_depth=0, ov_pan=0;
  for (int i=0;i<N;i++)
    { int tl=Otlen[i];
      if (tl==-1){ov_band++;continue;} if (tl==-2){ov_depth++;continue;} if (tl==-3){ov_pan++;continue;}
      valid++;
      int rtl=refTlen[i];
      int tok = (tl==rtl), dok=(Odiffs[i]==refDiffs[i]);
      // Δb sum over our panels vs bw (the Check_Trace_Points B invariant)
      const unsigned short *g = &Out[(size_t)i*2*PCAP];
      int gb=0, gd=0; for (int p=0;p<tl/2;p++){ gd+=g[2*p]; gb+=g[2*p+1]; }
      int bok = (gb==BW[i]);
      if (tok) tlen_ok++; if (dok) diff_ok++; if (bok) bsum_ok++;
      // Check_Trace_Points-valid: correct panel count, diffs sum to edit distance, Δb sum to bw
      if (tok && dok && bok && gd==Odiffs[i]) ctp_ok++;
      if (tok && dok)
        { int allpair=1;
          const unsigned short *r = &refTr[refOff[i]];
          for (int p=0;p<rtl;p++) if (g[p]!=r[p]){allpair=0;break;}
          if (allpair) exact++;
        }
    }
  printf("valid(computed): %ld / %d   [overflow band=%ld depth=%ld panels=%ld]\n",
         valid,N,ov_band,ov_depth,ov_pan);
  printf("  tlen match:        %ld (%.2f%% of valid)\n", tlen_ok, 100.0*tlen_ok/valid);
  printf("  total-diffs match: %ld (%.2f%% of valid)\n", diff_ok, 100.0*diff_ok/valid);
  printf("  Delta-b sum == bw: %ld (%.2f%% of valid)\n", bsum_ok, 100.0*bsum_ok/valid);
  printf("  Check_Trace_Points-VALID (tlen+diffsum+bsum): %ld (%.2f%% of valid)\n", ctp_ok, 100.0*ctp_ok/valid);
  printf("  EXACT bit match vs FastGA: %ld (%.2f%% of valid)\n", exact, 100.0*exact/valid);
  return 0;
}
