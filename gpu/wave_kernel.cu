/*******************************************************************************************
 *
 *  wave_kernel.cu -- GPU genome residency + (future) wave kernel for the wave-parallelism
 *                    characterization study (Task 4: context + loaders ONLY; no wave kernel
 *                    yet -- that is Tasks 5/6). See gpu/wave_kernel.h for the residency
 *                    layout contract and gpu/WAVE_PORT_NOTES.md for the align.c wave
 *                    contract this residency layout must satisfy (sentinel bytes at contig
 *                    boundaries).
 *
 *  Mirrors gpu/fastga_gpu.cu's residency pattern (gpu_open/gpu_close/gpu_load_seqs): a
 *  persistent device context, cudaMalloc-once-grow-as-needed device arrays, plain H2D
 *  cudaMemcpy. The difference from fastga_gpu.cu is scale (whole genomes, not one contig
 *  pair at a time) and that B is resident in BOTH orientations (Bfwd/Brev) to mirror
 *  wave_bench_cpu.c's CPU baseline -- see wave_kernel.h's header comment for why.
 *
 *  Build (library, into a driver):  nvcc -O3 -arch=sm_80 -Igpu -c gpu/wave_kernel.cu
 *
 *******************************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cuda_runtime.h>
#include "wave_kernel.h"

typedef unsigned char  u8;
typedef unsigned short uint16;

#define CK(c) do{ cudaError_t e=(c); \
  if(e!=cudaSuccess){ fprintf(stderr,"CUDA %s:%d %s\n",__FILE__,__LINE__,cudaGetErrorString(e)); } \
}while(0)

struct wave_ctx
  { u8   *dA, *dBfwd, *dBrev;      // resident flat NUMERIC genome arrays (device)
    long *dAbase, *dBbase;         // resident per-contig base-offset tables (device)
    size_t Acap, Bcap;             // device byte capacity of dA / dBfwd&dBrev (grow-only)
    size_t AbaseCap, BbaseCap;     // device element capacity (longs) of dAbase / dBbase
    int   nA, nB;                 // # contigs currently resident
    long  Alen, Blen;              // flat array lengths (bytes) currently resident
  };

extern "C" wave_ctx *wave_open(void)
{ return (wave_ctx *) calloc(1,sizeof(wave_ctx));
}

extern "C" void wave_close(wave_ctx *g)
{ if (g == NULL)
    return;
  if (g->dA)     cudaFree(g->dA);
  if (g->dBfwd)  cudaFree(g->dBfwd);
  if (g->dBrev)  cudaFree(g->dBrev);
  if (g->dAbase) cudaFree(g->dAbase);
  if (g->dBbase) cudaFree(g->dBbase);
  free(g);
}

extern "C" void wave_load_genomes(wave_ctx *g,
        const unsigned char *A,    const long *Abase, int nA, long Alen,
        const unsigned char *Bfwd, const unsigned char *Brev,
                                    const long *Bbase, int nB, long Blen)
{ if ((size_t) Alen > g->Acap)
    { if (g->dA) cudaFree(g->dA);
      CK(cudaMalloc(&g->dA,(size_t) Alen));
      g->Acap = (size_t) Alen;
    }
  if ((size_t) Blen > g->Bcap)
    { if (g->dBfwd) cudaFree(g->dBfwd);
      if (g->dBrev) cudaFree(g->dBrev);
      CK(cudaMalloc(&g->dBfwd,(size_t) Blen));
      CK(cudaMalloc(&g->dBrev,(size_t) Blen));
      g->Bcap = (size_t) Blen;
    }
  if ((size_t) nA > g->AbaseCap)
    { if (g->dAbase) cudaFree(g->dAbase);
      CK(cudaMalloc(&g->dAbase,(size_t) nA*sizeof(long)));
      g->AbaseCap = (size_t) nA;
    }
  if ((size_t) nB > g->BbaseCap)
    { if (g->dBbase) cudaFree(g->dBbase);
      CK(cudaMalloc(&g->dBbase,(size_t) nB*sizeof(long)));
      g->BbaseCap = (size_t) nB;
    }

  CK(cudaMemcpy(g->dA,    A,    (size_t) Alen,          cudaMemcpyHostToDevice));
  CK(cudaMemcpy(g->dBfwd, Bfwd, (size_t) Blen,          cudaMemcpyHostToDevice));
  CK(cudaMemcpy(g->dBrev, Brev, (size_t) Blen,          cudaMemcpyHostToDevice));
  CK(cudaMemcpy(g->dAbase,Abase,(size_t) nA*sizeof(long),cudaMemcpyHostToDevice));
  CK(cudaMemcpy(g->dBbase,Bbase,(size_t) nB*sizeof(long),cudaMemcpyHostToDevice));

  g->nA = nA;  g->nB = nB;
  g->Alen = Alen;  g->Blen = Blen;
}

static int readback(const u8 *dArr, long cap, long offset, long n, unsigned char *out)
{ if (dArr == NULL || offset < 0 || n < 0 || offset+n > cap)
    return -1;
  CK(cudaMemcpy(out,dArr+offset,(size_t) n,cudaMemcpyDeviceToHost));
  return 0;
}

extern "C" int wave_readback_A(wave_ctx *g, long offset, long n, unsigned char *out)
{ return readback(g->dA,g->Alen,offset,n,out);
}
extern "C" int wave_readback_Bfwd(wave_ctx *g, long offset, long n, unsigned char *out)
{ return readback(g->dBfwd,g->Blen,offset,n,out);
}
extern "C" int wave_readback_Brev(wave_ctx *g, long offset, long n, unsigned char *out)
{ return readback(g->dBrev,g->Blen,offset,n,out);
}

/*******************************************************************************************
 *
 *  Task 5: Stage-1 FORWARD sweep kernel -- a faithful warp-cooperative port of align.c's
 *  forward_wave (align.c:385-916).  See gpu/WAVE_PORT_NOTES.md for the contract.
 *
 *  Design (correctness first):
 *    - one WARP per seed; 32 lanes tile the diagonal band [low,hgh].
 *    - the `dif` (wave-number) loop is SEQUENTIAL: diagonal k at wave dif reads V/T/M[k-1,k,k+1]
 *      from wave dif-1 (align.c:687-716) -- a strict antidiagonal recurrence.  Within a wave the
 *      per-diagonal recurrence+furthest-reach slide (align.c:722-744) runs in PARALLEL across
 *      lanes; a serial lane-0 pass then reproduces forward_wave's exact k-DESCENDING
 *      besta/trima trim tie-breaking (align.c:680-831), so ties resolve to the same diagonal the
 *      CPU picks.
 *    - only the FORWARD ENDPOINT (trimx,trimy=trima-trimx,trimd) is produced, so the cells[]/
 *      HA/NA checkpoint + trace-point machinery (align.c:746-768,846-911) is dropped: it feeds
 *      only apath->trace, never the endpoint (REACH=0 => trimy=trima-trimx, align.c:852-859).
 *    - each wave is stored REBASED to buffer index (k - curbase); ping-pong prev/cur buffers.
 *      Out-of-[prev_lo,prev_hi] reads return the align.c retired-diagonal sentinel V=-1 (forward
 *      maximizes reach, so a retired diagonal reads artificially LOW -- WAVE_PORT_NOTES 3.(3)),
 *      with don't-care T=PATH_INT/M=PATH_LEN (a -1 diagonal never wins the recurrence max).
 *    - int positions throughout (deep waves reach x>2^15); no depth cap.
 *
 *******************************************************************************************/

#define BVEC       unsigned long long
#define TRIM_LEN   15
#define PATH_LEN   60
#define PATH_TOP   0x1000000000000000ULL     // 1 << PATH_LEN
#define PATH_INT   0x0fffffffffffffffULL      // PATH_TOP - 1
#define TRIM_MASK  0x7fff                      // (1<<TRIM_LEN)-1
#define TRIM_MLAG  250                          // align.c:179
#define WAVE_LAG   70                           // align.c:180
#define WV_INTMAX  2147483647

#define WV_SPAN    2048    // per-warp diagonal-band scratch width (band bounded by WAVE_LAG)
#define WV_POOL    8192    // max resident warps; seeds are grid-strided over the pool

//  Per-warp scalar state shared across the 32 lanes (one warp per block).
struct WvState
  { int  low, hgh;          // current band (final of the last completed wave)
    int  prevbase, prev_lo, prev_hi;   // prev buffer's rebase origin + valid k-range
    int  besta, bestx, lasta;
    int  trima, trimx, trimd;
    int  more;
    int  dif;
    int  overflow;
  };

//  Read prev-wave furthest-reach V for diagonal k (align.c retired-diagonal sentinel = -1).
__device__ __forceinline__ int prevV(const int *Vb, const WvState &s, int k)
{ if (k < s.prev_lo || k > s.prev_hi) return -1;
  return Vb[k - s.prevbase];
}

/*  One warp sweeps one seed's forward wave.  Lane = threadIdx.x (blockDim.x == 32). */
__global__ void forward_sweep_warp(
        const u8 *A,    const long *Abase,
        const u8 *Bfwd, const u8 *Brev, const long *Bbase,
        int n, const wave_seed *seeds,
        int *ae, int *be, int *fdiff,
        int *Vpool, BVEC *Tpool, int *Mpool, unsigned char *Cpool,
        const short *TABLE, const short *SCORE, int PATH_AVE)
{
  int lane = threadIdx.x;
  int warp = blockIdx.x;                 // one warp per block
  __shared__ WvState s;

  //  per-warp scratch (two V/T/M buffers for ping-pong, one clip buffer)
  int  *Vbuf0 = Vpool + ((size_t)warp*2 + 0)*WV_SPAN;
  int  *Vbuf1 = Vpool + ((size_t)warp*2 + 1)*WV_SPAN;
  BVEC *Tbuf0 = Tpool + ((size_t)warp*2 + 0)*WV_SPAN;
  BVEC *Tbuf1 = Tpool + ((size_t)warp*2 + 1)*WV_SPAN;
  int  *Mbuf0 = Mpool + ((size_t)warp*2 + 0)*WV_SPAN;
  int  *Mbuf1 = Mpool + ((size_t)warp*2 + 1)*WV_SPAN;
  unsigned char *Cbuf = Cpool + (size_t)warp*WV_SPAN;

  for (int si = warp; si < n; si += gridDim.x)
    { const wave_seed sd = seeds[si];

      const u8 *aseq = A + Abase[sd.aread];                       // aseq[-1]==4, aseq[alen]==4
      const u8 *bseq = (sd.comp ? Brev : Bfwd) + Bbase[sd.bread];
      int mida = sd.anti;

      //  ping-pong: `cur` is the wave we are computing, `prev` the one before it
      int *curV = Vbuf0, *prevVb = Vbuf1;
      BVEC *curT = Tbuf0, *prevTb = Tbuf1;
      int *curM = Mbuf0, *prevMb = Mbuf1;

      if (lane == 0)
        { int low = sd.diag, hgh = sd.diag;
          while (((mida - hgh) >> 1) < 0) hgh -= 1;      // align.c:1504-1505 anti clamp
          s.low = low; s.hgh = hgh;
          s.besta = s.trima = s.lasta = mida;
          s.bestx = s.trimx = (mida + hgh) >> 1;
          s.trimd = 0;
          s.more  = 1;
          s.dif   = 0;
          s.overflow = 0;
        }
      __syncwarp();

      int low = s.low, hgh = s.hgh;
      int curbase = low;

      //  overflow guard: 0-wave span (typically 1) must fit
      if (hgh - low + 1 > WV_SPAN)
        { if (lane == 0) s.overflow = 1; }
      else
        {
          //  ---- 0-wave (align.c:462-548): furthest reach per diagonal from the seed anti ----
          for (int k = low + lane; k <= hgh; k += 32)
            { int x = (mida + k) >> 1;
              int code = 0;                               // 0 none, 1 A-clip(d==4), 2 B-clip(c==4)
              while (1)
                { int cB = bseq[x - k];
                  if (cB == 4) { code = 2; break; }
                  int dA = aseq[x];
                  if (cB != dA) { if (dA == 4) code = 1; break; }
                  x += 1;
                }
              int c = (x << 1) - k;
              int idx = k - curbase;
              curV[idx] = c; curT[idx] = PATH_INT; curM[idx] = PATH_LEN; Cbuf[idx] = (unsigned char) code;
            }
          __syncwarp();

          //  ---- 0-wave serial trim (align.c:535-539, unconditional new-best; align.c:551-574 clip) ----
          if (lane == 0)
            { int aclip = WV_INTMAX, bclip = -WV_INTMAX;
              for (int k = hgh; k >= low; k--)
                { int idx = k - curbase;
                  int c = curV[idx]; int x = (c + k) >> 1; int code = Cbuf[idx];
                  if (code == 2) { s.more = 0; if (bclip < k) bclip = k; }
                  else if (code == 1) { s.more = 0; aclip = k; }
                  if (c > s.besta)
                    { s.besta = s.trima = s.lasta = c; s.bestx = s.trimx = x; }
                }
              if (s.more == 0)
                { if (bseq[s.besta - s.bestx] != 4 && aseq[s.bestx] != 4) s.more = 1;
                  if (hgh >= aclip) hgh = aclip - 1;
                  if (low <= bclip) low = bclip + 1;
                }
              s.low = low; s.hgh = hgh;
              //  prev range for wave 1 = this (clipped) 0-wave band
              s.prevbase = curbase; s.prev_lo = low; s.prev_hi = hgh;
            }
          __syncwarp();
          low = s.low; hgh = s.hgh;

          //  swap: the buffer we just filled becomes `prev`
          { int *tV=curV; curV=prevVb; prevVb=tV;
            BVEC *tT=curT; curT=prevTb; prevTb=tT;
            int *tM=curM; curM=prevMb; prevMb=tM; }

          //  ---- steady waves (align.c:583-844) ----
          while (s.more && s.lasta >= s.besta - TRIM_MLAG && !s.overflow)
            { low = s.low - 1; hgh = s.hgh + 1;          // align.c:590-591 (minp/maxp never clip here)
              curbase = low;
              if (hgh - low + 1 > WV_SPAN) { if (lane==0) s.overflow = 1; __syncwarp(); break; }
              int dif = s.dif + 1;

              //  phase 1: parallel per-diagonal recurrence + furthest-reach slide
              for (int k = low + lane; k <= hgh; k += 32)
                { int ap = prevV(prevVb,s,k+1);
                  int ac = prevV(prevVb,s,k);
                  int am = prevV(prevVb,s,k-1);
                  int c; int m; BVEC b;
                  if (ac < am)
                    { if (am < ap) { c=ap+1; int j=k+1; b=(j<s.prev_lo||j>s.prev_hi)?PATH_INT:prevTb[j-s.prevbase]; m=(j<s.prev_lo||j>s.prev_hi)?PATH_LEN:prevMb[j-s.prevbase]; }
                      else          { c=am+1; int j=k-1; b=(j<s.prev_lo||j>s.prev_hi)?PATH_INT:prevTb[j-s.prevbase]; m=(j<s.prev_lo||j>s.prev_hi)?PATH_LEN:prevMb[j-s.prevbase]; } }
                  else
                    { if (ac < ap) { c=ap+1; int j=k+1; b=(j<s.prev_lo||j>s.prev_hi)?PATH_INT:prevTb[j-s.prevbase]; m=(j<s.prev_lo||j>s.prev_hi)?PATH_LEN:prevMb[j-s.prevbase]; }
                      else          { c=ac+2; int j=k;   b=(j<s.prev_lo||j>s.prev_hi)?PATH_INT:prevTb[j-s.prevbase]; m=(j<s.prev_lo||j>s.prev_hi)?PATH_LEN:prevMb[j-s.prevbase]; } }

                  if ((b & PATH_TOP) != 0) m -= 1;
                  b <<= 1;

                  int x = (c + k) >> 1;
                  int code = 0;
                  while (1)
                    { int cB = bseq[x - k];
                      if (cB == 4) { code = 2; break; }
                      int dA = aseq[x];
                      if (cB != dA) { if (dA == 4) code = 1; break; }
                      x += 1;
                      if ((b & PATH_TOP) == 0) m += 1;
                      b = (b << 1) | 1;
                    }
                  c = (x << 1) - k;
                  int idx = k - curbase;
                  curV[idx]=c; curM[idx]=m; curT[idx]=b; Cbuf[idx]=(unsigned char)code;
                }
              __syncwarp();

              //  phase 2: serial trim/clip/prune (align.c:680-831), exact k-descending order
              if (lane == 0)
                { int aclip = WV_INTMAX, bclip = -WV_INTMAX;
                  for (int k = hgh; k >= low; k--)
                    { int idx = k - curbase;
                      int c = curV[idx]; int m = curM[idx]; BVEC b = curT[idx];
                      int x = (c + k) >> 1; int code = Cbuf[idx];
                      if (code == 2) { s.more = 0; if (bclip < k) bclip = k; }
                      else if (code == 1) { s.more = 0; aclip = k; }
                      if (c > s.besta)
                        { s.besta = c; s.bestx = x;
                          if (m >= PATH_AVE)
                            { s.lasta = c;
                              if (TABLE[b & TRIM_MASK] >= 0)
                                if (TABLE[(b >> TRIM_LEN) & TRIM_MASK] + SCORE[b & TRIM_MASK] >= 0)
                                  { s.trima = c; s.trimx = x; s.trimd = dif; }
                            }
                        }
                    }
                  s.dif = dif;
                  //  more==0 clip handling (align.c:796-821); morea/REACH branch dropped (REACH=0)
                  if (s.more == 0)
                    { if (bseq[s.besta - s.bestx] != 4 && aseq[s.bestx] != 4) s.more = 1;
                      if (hgh >= aclip) hgh = aclip - 1;
                      if (low <= bclip) low = bclip + 1;
                    }
                  //  WAVE_LAG band prune (align.c:823-831)
                  { int nlag = s.besta - WAVE_LAG;
                    while (hgh >= low)
                      { if (curV[hgh - curbase] < nlag) hgh -= 1;
                        else { while (curV[low - curbase] < nlag) low += 1; break; }
                      }
                  }
                  s.low = low; s.hgh = hgh;
                  s.prevbase = curbase; s.prev_lo = low; s.prev_hi = hgh;
                }
              __syncwarp();

              //  swap cur<->prev for the next wave
              { int *tV=curV; curV=prevVb; prevVb=tV;
                BVEC *tT=curT; curT=prevTb; prevTb=tT;
                int *tM=curM; curM=prevMb; prevMb=tM; }
            }
        }

      if (lane == 0)
        { if (s.overflow)
            { ae[si] = -2; be[si] = -2; fdiff[si] = -2; }
          else
            { ae[si] = s.trimx; be[si] = s.trima - s.trimx; fdiff[si] = s.trimd; }
        }
      __syncwarp();
    }
}

extern "C" int wave_forward_batch(wave_ctx *g, int n, const wave_seed *seeds,
                        const int16_t *table, const int16_t *score, int path_ave,
                        int *ae, int *be, int *fdiff)
{ if (g == NULL || n <= 0) return 0;

  wave_seed *dSeeds = NULL;
  int   *dAe=NULL, *dBe=NULL, *dFd=NULL;
  short *dTab=NULL, *dScr=NULL;
  int   *dV=NULL, *dM=NULL; BVEC *dT=NULL; unsigned char *dC=NULL;

  int pool = n < WV_POOL ? n : WV_POOL;
  size_t vtm = (size_t) pool * 2 * WV_SPAN;

  CK(cudaMalloc(&dSeeds,(size_t) n * sizeof(wave_seed)));
  CK(cudaMalloc(&dAe,(size_t) n * sizeof(int)));
  CK(cudaMalloc(&dBe,(size_t) n * sizeof(int)));
  CK(cudaMalloc(&dFd,(size_t) n * sizeof(int)));
  CK(cudaMalloc(&dTab,(size_t)(TRIM_MASK+1) * sizeof(short)));
  CK(cudaMalloc(&dScr,(size_t)(TRIM_MASK+1) * sizeof(short)));
  CK(cudaMalloc(&dV, vtm * sizeof(int)));
  CK(cudaMalloc(&dT, vtm * sizeof(BVEC)));
  CK(cudaMalloc(&dM, vtm * sizeof(int)));
  CK(cudaMalloc(&dC,(size_t) pool * WV_SPAN * sizeof(unsigned char)));
  if (dSeeds==NULL||dAe==NULL||dBe==NULL||dFd==NULL||dTab==NULL||dScr==NULL||
      dV==NULL||dT==NULL||dM==NULL||dC==NULL)
    { fprintf(stderr,"wave_forward_batch: device alloc failed\n"); return -1; }

  CK(cudaMemcpy(dSeeds,seeds,(size_t) n * sizeof(wave_seed),cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dTab,table,(size_t)(TRIM_MASK+1) * sizeof(short),cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dScr,score,(size_t)(TRIM_MASK+1) * sizeof(short),cudaMemcpyHostToDevice));

  forward_sweep_warp<<<pool,32>>>(g->dA,g->dAbase,g->dBfwd,g->dBrev,g->dBbase,
        n,dSeeds,dAe,dBe,dFd,dV,dT,dM,dC,dTab,dScr,path_ave);
  cudaError_t err = cudaGetLastError();
  if (err == cudaSuccess) err = cudaDeviceSynchronize();
  if (err != cudaSuccess)
    { fprintf(stderr,"wave_forward_batch: kernel error: %s\n",cudaGetErrorString(err));
      return -1;
    }

  CK(cudaMemcpy(ae,dAe,(size_t) n * sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(be,dBe,(size_t) n * sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(fdiff,dFd,(size_t) n * sizeof(int),cudaMemcpyDeviceToHost));

  cudaFree(dSeeds); cudaFree(dAe); cudaFree(dBe); cudaFree(dFd);
  cudaFree(dTab); cudaFree(dScr);
  cudaFree(dV); cudaFree(dT); cudaFree(dM); cudaFree(dC);
  return 0;
}

/*******************************************************************************************
 *
 *  Task 6: Stage-1 REVERSE sweep kernel -- a faithful warp-cooperative mirror of align.c's
 *  reverse_wave (align.c:919-1459).  The mirror of forward_sweep_warp above, differing only in
 *  the direction-dependent details called out in gpu/WAVE_PORT_NOTES.md section 3.(3):
 *    - extends LEFTWARD: the furthest-reach slide walks x-=1 (align.c:1037,1259), and A/B are
 *      indexed off pointers pre-decremented by 1 (align.c:921-922 `aseq = align->aseq - 1`).
 *    - reverse MINIMIZES reach, so the retired-diagonal sentinel is INT32_MAX (align.c:1181/
 *      1188/1195), the best test is `c < besta` (align.c:1290), the recurrence picks the MIN
 *      source and subtracts (am/ap-1, ac-2, align.c:1211-1236), and the WAVE_LAG prune drops
 *      diagonals whose V has risen above besta+WAVE_LAG (align.c:1343-1351).
 *    - trim continuation is `lasta <= besta + TRIM_MLAG` (align.c:1108).
 *    - the more==0 clip roles are swapped vs forward: A-clip bounds LOW, B-clip bounds HIGH
 *      (align.c:1319-1330), and bclip minimizes (align.c:1247 `bclip > k`).
 *    - the k-loop is ASCENDING (align.c:1200), so the serial trim pass visits k low->hgh and
 *      first (lowest-k) wins ties, matching the CPU.
 *  Endpoint-only, same as forward: cells[]/HA/NA trace machinery is dropped (REACH=0), and the
 *  reverse start band is [seed.diag,seed.diag] (see wave_kernel.h -- forward's *mind == dg).
 *
 *******************************************************************************************/

//  Read prev-wave furthest-reach V for reverse (retired-diagonal sentinel = INT32_MAX).
__device__ __forceinline__ int prevVr(const int *Vb, const WvState &s, int k)
{ if (k < s.prev_lo || k > s.prev_hi) return WV_INTMAX;
  return Vb[k - s.prevbase];
}

/*  One warp sweeps one seed's reverse wave.  Lane = threadIdx.x (blockDim.x == 32). */
__global__ void reverse_sweep_warp(
        const u8 *A,    const long *Abase,
        const u8 *Bfwd, const u8 *Brev, const long *Bbase,
        int n, const wave_seed *seeds,
        int *ab, int *bb, int *rdiff,
        int *Vpool, BVEC *Tpool, int *Mpool, unsigned char *Cpool,
        const short *TABLE, const short *SCORE, int PATH_AVE)
{
  int lane = threadIdx.x;
  int warp = blockIdx.x;
  __shared__ WvState s;

  int  *Vbuf0 = Vpool + ((size_t)warp*2 + 0)*WV_SPAN;
  int  *Vbuf1 = Vpool + ((size_t)warp*2 + 1)*WV_SPAN;
  BVEC *Tbuf0 = Tpool + ((size_t)warp*2 + 0)*WV_SPAN;
  BVEC *Tbuf1 = Tpool + ((size_t)warp*2 + 1)*WV_SPAN;
  int  *Mbuf0 = Mpool + ((size_t)warp*2 + 0)*WV_SPAN;
  int  *Mbuf1 = Mpool + ((size_t)warp*2 + 1)*WV_SPAN;
  unsigned char *Cbuf = Cpool + (size_t)warp*WV_SPAN;

  for (int si = warp; si < n; si += gridDim.x)
    { const wave_seed sd = seeds[si];

      //  reverse offsets both pointers by -1 (align.c:921-922); aseq[-1]/bseq[-1] are the
      //  contig's leading `4` sentinel, so walking left off base 0 clips exactly as the CPU.
      const u8 *aseq = A + Abase[sd.aread] - 1;
      const u8 *bseq = (sd.comp ? Brev : Bfwd) + Bbase[sd.bread] - 1;
      int mida = sd.anti;

      int *curV = Vbuf0, *prevVb = Vbuf1;
      BVEC *curT = Tbuf0, *prevTb = Tbuf1;
      int *curM = Mbuf0, *prevMb = Mbuf1;

      if (lane == 0)
        { int low = sd.diag, hgh = sd.diag;      // reverse band = [dg,dg] (== forward's *mind)
          s.low = low; s.hgh = hgh;
          s.besta = s.trima = s.lasta = mida;
          s.bestx = s.trimx = (mida + hgh) >> 1;
          s.trimd = 0;
          s.more  = 1;
          s.dif   = 0;
          s.overflow = 0;
        }
      __syncwarp();

      int low = s.low, hgh = s.hgh;
      int curbase = low;

      if (hgh - low + 1 > WV_SPAN)
        { if (lane == 0) s.overflow = 1; }
      else
        {
          //  ---- 0-wave (align.c:990-1076): furthest reach per diagonal, sliding LEFT ----
          for (int k = low + lane; k <= hgh; k += 32)
            { int x = (mida + k) >> 1;
              int code = 0;                               // 0 none, 1 A-clip(d==4), 2 B-clip(c==4)
              while (1)
                { int cB = bseq[x - k];
                  if (cB == 4) { code = 2; break; }
                  int dA = aseq[x];
                  if (cB != dA) { if (dA == 4) code = 1; break; }
                  x -= 1;
                }
              int c = (x << 1) - k;
              int idx = k - curbase;
              curV[idx] = c; curT[idx] = PATH_INT; curM[idx] = PATH_LEN; Cbuf[idx] = (unsigned char) code;
            }
          __syncwarp();

          //  ---- 0-wave serial trim (align.c:1062-1066 unconditional; 1078-1101 clip) ----
          if (lane == 0)
            { int aclip = -WV_INTMAX, bclip = WV_INTMAX;
              for (int k = low; k <= hgh; k++)               // reverse: ASCENDING (align.c:994)
                { int idx = k - curbase;
                  int c = curV[idx]; int x = (c + k) >> 1; int code = Cbuf[idx];
                  if (code == 2) { s.more = 0; if (bclip > k) bclip = k; }
                  else if (code == 1) { s.more = 0; aclip = k; }
                  if (c < s.besta)
                    { s.besta = s.trima = s.lasta = c; s.bestx = s.trimx = x; }
                }
              if (s.more == 0)
                { if (bseq[s.besta - s.bestx] != 4 && aseq[s.bestx] != 4) s.more = 1;
                  if (low <= aclip) low = aclip + 1;         // A-clip bounds LOW (align.c:1319)
                  if (hgh >= bclip) hgh = bclip - 1;         // B-clip bounds HIGH (align.c:1329)
                }
              s.low = low; s.hgh = hgh;
              s.prevbase = curbase; s.prev_lo = low; s.prev_hi = hgh;
            }
          __syncwarp();
          low = s.low; hgh = s.hgh;

          { int *tV=curV; curV=prevVb; prevVb=tV;
            BVEC *tT=curT; curT=prevTb; prevTb=tT;
            int *tM=curM; curM=prevMb; prevMb=tM; }

          //  ---- steady waves (align.c:1108-1364) ----
          while (s.more && s.lasta <= s.besta + TRIM_MLAG && !s.overflow)
            { low = s.low - 1; hgh = s.hgh + 1;          // align.c:1115-1116 (minp/maxp never clip)
              curbase = low;
              if (hgh - low + 1 > WV_SPAN) { if (lane==0) s.overflow = 1; __syncwarp(); break; }
              int dif = s.dif + 1;

              //  phase 1: parallel per-diagonal recurrence + furthest-reach slide (MIN, x-=1)
              for (int k = low + lane; k <= hgh; k += 32)
                { int ap = prevVr(prevVb,s,k+1);
                  int ac = prevVr(prevVb,s,k);
                  int am = prevVr(prevVb,s,k-1);
                  int c; int m; BVEC b;
                  if (ac > ap)
                    { if (ap > am) { c=am-1; int j=k-1; b=(j<s.prev_lo||j>s.prev_hi)?PATH_INT:prevTb[j-s.prevbase]; m=(j<s.prev_lo||j>s.prev_hi)?PATH_LEN:prevMb[j-s.prevbase]; }
                      else          { c=ap-1; int j=k+1; b=(j<s.prev_lo||j>s.prev_hi)?PATH_INT:prevTb[j-s.prevbase]; m=(j<s.prev_lo||j>s.prev_hi)?PATH_LEN:prevMb[j-s.prevbase]; } }
                  else
                    { if (ac > am) { c=am-1; int j=k-1; b=(j<s.prev_lo||j>s.prev_hi)?PATH_INT:prevTb[j-s.prevbase]; m=(j<s.prev_lo||j>s.prev_hi)?PATH_LEN:prevMb[j-s.prevbase]; }
                      else          { c=ac-2; int j=k;   b=(j<s.prev_lo||j>s.prev_hi)?PATH_INT:prevTb[j-s.prevbase]; m=(j<s.prev_lo||j>s.prev_hi)?PATH_LEN:prevMb[j-s.prevbase]; } }

                  if ((b & PATH_TOP) != 0) m -= 1;
                  b <<= 1;

                  int x = (c + k) >> 1;
                  int code = 0;
                  while (1)
                    { int cB = bseq[x - k];
                      if (cB == 4) { code = 2; break; }
                      int dA = aseq[x];
                      if (cB != dA) { if (dA == 4) code = 1; break; }
                      x -= 1;
                      if ((b & PATH_TOP) == 0) m += 1;
                      b = (b << 1) | 1;
                    }
                  c = (x << 1) - k;
                  int idx = k - curbase;
                  curV[idx]=c; curM[idx]=m; curT[idx]=b; Cbuf[idx]=(unsigned char)code;
                }
              __syncwarp();

              //  phase 2: serial trim/clip/prune (align.c:1200-1351), exact k-ASCENDING order
              if (lane == 0)
                { int aclip = -WV_INTMAX, bclip = WV_INTMAX;
                  for (int k = low; k <= hgh; k++)
                    { int idx = k - curbase;
                      int c = curV[idx]; int m = curM[idx]; BVEC b = curT[idx];
                      int x = (c + k) >> 1; int code = Cbuf[idx];
                      if (code == 2) { s.more = 0; if (bclip > k) bclip = k; }
                      else if (code == 1) { s.more = 0; aclip = k; }
                      if (c < s.besta)
                        { s.besta = c; s.bestx = x;
                          if (m >= PATH_AVE)
                            { s.lasta = c;
                              if (TABLE[b & TRIM_MASK] >= 0)
                                if (TABLE[(b >> TRIM_LEN) & TRIM_MASK] + SCORE[b & TRIM_MASK] >= 0)
                                  { s.trima = c; s.trimx = x; s.trimd = dif; }
                            }
                        }
                    }
                  s.dif = dif;
                  if (s.more == 0)
                    { if (bseq[s.besta - s.bestx] != 4 && aseq[s.bestx] != 4) s.more = 1;
                      if (low <= aclip) low = aclip + 1;
                      if (hgh >= bclip) hgh = bclip - 1;
                    }
                  //  WAVE_LAG band prune (align.c:1343-1351): reverse drops V > besta+WAVE_LAG
                  { int nlag = s.besta + WAVE_LAG;
                    while (hgh >= low)
                      { if (curV[hgh - curbase] > nlag) hgh -= 1;
                        else { while (curV[low - curbase] > nlag) low += 1; break; }
                      }
                  }
                  s.low = low; s.hgh = hgh;
                  s.prevbase = curbase; s.prev_lo = low; s.prev_hi = hgh;
                }
              __syncwarp();

              { int *tV=curV; curV=prevVb; prevVb=tV;
                BVEC *tT=curT; curT=prevTb; prevTb=tT;
                int *tM=curM; curM=prevMb; prevMb=tM; }
            }
        }

      if (lane == 0)
        { if (s.overflow)
            { ab[si] = -2; bb[si] = -2; rdiff[si] = -2; }
          else
            { ab[si] = s.trimx; bb[si] = s.trima - s.trimx; rdiff[si] = s.trimd; }
        }
      __syncwarp();
    }
}

extern "C" int wave_discover_batch(wave_ctx *g, int n, const wave_seed *seeds,
                        const int16_t *table, const int16_t *score, int path_ave,
                        int *ab, int *ae, int *bb, int *be, int *diffs)
{ if (g == NULL || n <= 0) return 0;

  wave_seed *dSeeds = NULL;
  int   *dAe=NULL, *dBe=NULL, *dFd=NULL;    // forward endpoints + forward diff
  int   *dAb=NULL, *dBb=NULL, *dRd=NULL;    // reverse endpoints + reverse diff
  short *dTab=NULL, *dScr=NULL;
  int   *dV=NULL, *dM=NULL; BVEC *dT=NULL; unsigned char *dC=NULL;

  int pool = n < WV_POOL ? n : WV_POOL;
  size_t vtm = (size_t) pool * 2 * WV_SPAN;

  CK(cudaMalloc(&dSeeds,(size_t) n * sizeof(wave_seed)));
  CK(cudaMalloc(&dAe,(size_t) n * sizeof(int)));
  CK(cudaMalloc(&dBe,(size_t) n * sizeof(int)));
  CK(cudaMalloc(&dFd,(size_t) n * sizeof(int)));
  CK(cudaMalloc(&dAb,(size_t) n * sizeof(int)));
  CK(cudaMalloc(&dBb,(size_t) n * sizeof(int)));
  CK(cudaMalloc(&dRd,(size_t) n * sizeof(int)));
  CK(cudaMalloc(&dTab,(size_t)(TRIM_MASK+1) * sizeof(short)));
  CK(cudaMalloc(&dScr,(size_t)(TRIM_MASK+1) * sizeof(short)));
  CK(cudaMalloc(&dV, vtm * sizeof(int)));
  CK(cudaMalloc(&dT, vtm * sizeof(BVEC)));
  CK(cudaMalloc(&dM, vtm * sizeof(int)));
  CK(cudaMalloc(&dC,(size_t) pool * WV_SPAN * sizeof(unsigned char)));
  if (dSeeds==NULL||dAe==NULL||dBe==NULL||dFd==NULL||dAb==NULL||dBb==NULL||dRd==NULL||
      dTab==NULL||dScr==NULL||dV==NULL||dT==NULL||dM==NULL||dC==NULL)
    { fprintf(stderr,"wave_discover_batch: device alloc failed\n"); return -1; }

  CK(cudaMemcpy(dSeeds,seeds,(size_t) n * sizeof(wave_seed),cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dTab,table,(size_t)(TRIM_MASK+1) * sizeof(short),cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dScr,score,(size_t)(TRIM_MASK+1) * sizeof(short),cudaMemcpyHostToDevice));

  //  forward sweep -> (ae,be,fdiff); the V/T/M/C pools are reused by the reverse launch below
  forward_sweep_warp<<<pool,32>>>(g->dA,g->dAbase,g->dBfwd,g->dBrev,g->dBbase,
        n,dSeeds,dAe,dBe,dFd,dV,dT,dM,dC,dTab,dScr,path_ave);
  //  reverse sweep -> (ab,bb,rdiff), starting from [seed.diag,seed.diag] at the seed anti
  reverse_sweep_warp<<<pool,32>>>(g->dA,g->dAbase,g->dBfwd,g->dBrev,g->dBbase,
        n,dSeeds,dAb,dBb,dRd,dV,dT,dM,dC,dTab,dScr,path_ave);
  cudaError_t err = cudaGetLastError();
  if (err == cudaSuccess) err = cudaDeviceSynchronize();
  if (err != cudaSuccess)
    { fprintf(stderr,"wave_discover_batch: kernel error: %s\n",cudaGetErrorString(err));
      return -1;
    }

  int *hFd = (int *) malloc((size_t) n * sizeof(int));
  int *hRd = (int *) malloc((size_t) n * sizeof(int));
  if (hFd==NULL || hRd==NULL) { fprintf(stderr,"wave_discover_batch: host alloc failed\n"); return -1; }

  CK(cudaMemcpy(ae, dAe,(size_t) n * sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(be, dBe,(size_t) n * sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hFd,dFd,(size_t) n * sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(ab, dAb,(size_t) n * sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(bb, dBb,(size_t) n * sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hRd,dRd,(size_t) n * sizeof(int),cudaMemcpyDeviceToHost));

  //  combine: diffs = fdiff + rdiff (align.c:1453 apath->diffs = apath->diffs + trimd); a band
  //  overflow in EITHER sweep marks the whole seed -2 (endpoints already carry -2 from the kernel)
  for (int i = 0; i < n; i++)
    { if (hFd[i] == -2 || hRd[i] == -2)
        { ab[i]=ae[i]=bb[i]=be[i]=diffs[i] = -2; }
      else
        diffs[i] = hFd[i] + hRd[i];
    }

  free(hFd); free(hRd);
  cudaFree(dSeeds); cudaFree(dAe); cudaFree(dBe); cudaFree(dFd);
  cudaFree(dAb); cudaFree(dBb); cudaFree(dRd);
  cudaFree(dTab); cudaFree(dScr);
  cudaFree(dV); cudaFree(dT); cudaFree(dM); cudaFree(dC);
  return 0;
}

/*******************************************************************************************
 *
 *  Task 7: CUDA-event-timed discovery + occupancy readout.  Body is wave_discover_batch's
 *  above, unchanged, with three cudaEvent windows dropped in around exactly the phases the
 *  Task 7 brief asks to separate (per-batch seed H2D / kernel-only / result D2H). The tiny
 *  one-time TABLE/SCORE upload is deliberately OUTSIDE the timed windows (a real deployment
 *  loads it once, like the genome) so it cannot inflate either reported basis.
 *
 *******************************************************************************************/

extern "C" int wave_discover_batch_timed(wave_ctx *g, int n, const wave_seed *seeds,
                        const int16_t *table, const int16_t *score, int path_ave,
                        int *ab, int *ae, int *bb, int *be, int *diffs,
                        float *ms_h2d, float *ms_kernel, float *ms_d2h)
{ if (g == NULL || n <= 0) return 0;

  wave_seed *dSeeds = NULL;
  int   *dAe=NULL, *dBe=NULL, *dFd=NULL;
  int   *dAb=NULL, *dBb=NULL, *dRd=NULL;
  short *dTab=NULL, *dScr=NULL;
  int   *dV=NULL, *dM=NULL; BVEC *dT=NULL; unsigned char *dC=NULL;

  int pool = n < WV_POOL ? n : WV_POOL;
  size_t vtm = (size_t) pool * 2 * WV_SPAN;

  CK(cudaMalloc(&dSeeds,(size_t) n * sizeof(wave_seed)));
  CK(cudaMalloc(&dAe,(size_t) n * sizeof(int)));
  CK(cudaMalloc(&dBe,(size_t) n * sizeof(int)));
  CK(cudaMalloc(&dFd,(size_t) n * sizeof(int)));
  CK(cudaMalloc(&dAb,(size_t) n * sizeof(int)));
  CK(cudaMalloc(&dBb,(size_t) n * sizeof(int)));
  CK(cudaMalloc(&dRd,(size_t) n * sizeof(int)));
  CK(cudaMalloc(&dTab,(size_t)(TRIM_MASK+1) * sizeof(short)));
  CK(cudaMalloc(&dScr,(size_t)(TRIM_MASK+1) * sizeof(short)));
  CK(cudaMalloc(&dV, vtm * sizeof(int)));
  CK(cudaMalloc(&dT, vtm * sizeof(BVEC)));
  CK(cudaMalloc(&dM, vtm * sizeof(int)));
  CK(cudaMalloc(&dC,(size_t) pool * WV_SPAN * sizeof(unsigned char)));
  if (dSeeds==NULL||dAe==NULL||dBe==NULL||dFd==NULL||dAb==NULL||dBb==NULL||dRd==NULL||
      dTab==NULL||dScr==NULL||dV==NULL||dT==NULL||dM==NULL||dC==NULL)
    { fprintf(stderr,"wave_discover_batch_timed: device alloc failed\n"); return -1; }

  //  one-time-in-a-real-deployment tables -- outside the timed windows, see header comment
  CK(cudaMemcpy(dTab,table,(size_t)(TRIM_MASK+1) * sizeof(short),cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dScr,score,(size_t)(TRIM_MASK+1) * sizeof(short),cudaMemcpyHostToDevice));

  cudaEvent_t e0,e1,e2,e3;
  cudaEventCreate(&e0); cudaEventCreate(&e1); cudaEventCreate(&e2); cudaEventCreate(&e3);

  //  ---- window 1: per-batch seed H2D ----
  cudaEventRecord(e0);
  CK(cudaMemcpy(dSeeds,seeds,(size_t) n * sizeof(wave_seed),cudaMemcpyHostToDevice));
  cudaEventRecord(e1);

  //  ---- window 2: kernel-only (both sweeps) ----
  forward_sweep_warp<<<pool,32>>>(g->dA,g->dAbase,g->dBfwd,g->dBrev,g->dBbase,
        n,dSeeds,dAe,dBe,dFd,dV,dT,dM,dC,dTab,dScr,path_ave);
  reverse_sweep_warp<<<pool,32>>>(g->dA,g->dAbase,g->dBfwd,g->dBrev,g->dBbase,
        n,dSeeds,dAb,dBb,dRd,dV,dT,dM,dC,dTab,dScr,path_ave);
  cudaEventRecord(e2);

  cudaError_t err = cudaGetLastError();
  if (err == cudaSuccess) err = cudaEventSynchronize(e2);
  if (err != cudaSuccess)
    { fprintf(stderr,"wave_discover_batch_timed: kernel error: %s\n",cudaGetErrorString(err));
      cudaEventDestroy(e0); cudaEventDestroy(e1); cudaEventDestroy(e2); cudaEventDestroy(e3);
      cudaFree(dSeeds); cudaFree(dAe); cudaFree(dBe); cudaFree(dFd);
      cudaFree(dAb); cudaFree(dBb); cudaFree(dRd);
      cudaFree(dTab); cudaFree(dScr);
      cudaFree(dV); cudaFree(dT); cudaFree(dM); cudaFree(dC);
      return -1;
    }

  int *hFd = (int *) malloc((size_t) n * sizeof(int));
  int *hRd = (int *) malloc((size_t) n * sizeof(int));
  if (hFd==NULL || hRd==NULL)
    { fprintf(stderr,"wave_discover_batch_timed: host alloc failed\n"); return -1; }

  //  ---- window 3: result D2H (all 6 result arrays) ----
  cudaEventRecord(e2);   // re-mark start of D2H precisely at kernel-end (already synced above)
  CK(cudaMemcpy(ae, dAe,(size_t) n * sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(be, dBe,(size_t) n * sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hFd,dFd,(size_t) n * sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(ab, dAb,(size_t) n * sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(bb, dBb,(size_t) n * sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hRd,dRd,(size_t) n * sizeof(int),cudaMemcpyDeviceToHost));
  cudaEventRecord(e3);
  cudaEventSynchronize(e3);

  float t_h2d=0, t_kernel=0, t_d2h=0;
  cudaEventElapsedTime(&t_h2d,   e0, e1);
  cudaEventElapsedTime(&t_kernel,e1, e2);
  cudaEventElapsedTime(&t_d2h,   e2, e3);
  if (ms_h2d)    *ms_h2d    = t_h2d;
  if (ms_kernel) *ms_kernel = t_kernel;
  if (ms_d2h)    *ms_d2h    = t_d2h;

  cudaEventDestroy(e0); cudaEventDestroy(e1); cudaEventDestroy(e2); cudaEventDestroy(e3);

  //  combine: diffs = fdiff + rdiff (align.c:1453); overflow in EITHER sweep marks -2
  for (int i = 0; i < n; i++)
    { if (hFd[i] == -2 || hRd[i] == -2)
        { ab[i]=ae[i]=bb[i]=be[i]=diffs[i] = -2; }
      else
        diffs[i] = hFd[i] + hRd[i];
    }

  free(hFd); free(hRd);
  cudaFree(dSeeds); cudaFree(dAe); cudaFree(dBe); cudaFree(dFd);
  cudaFree(dAb); cudaFree(dBb); cudaFree(dRd);
  cudaFree(dTab); cudaFree(dScr);
  cudaFree(dV); cudaFree(dT); cudaFree(dM); cudaFree(dC);
  return 0;
}

/*  Task 7: occupancy readout for the mechanism decomposition -- see wave_kernel.h contract. */
extern "C" void wave_query_occupancy(int *maxBlocksFwd, int *maxBlocksRev, int *smCount)
{ int dev = 0;
  cudaDeviceProp prop;

  cudaGetDevice(&dev);
  cudaGetDeviceProperties(&prop,dev);
  *smCount = prop.multiProcessorCount;

  int mbf = 0, mbr = 0;
  cudaOccupancyMaxActiveBlocksPerMultiprocessor(&mbf,forward_sweep_warp,32,0);
  cudaOccupancyMaxActiveBlocksPerMultiprocessor(&mbr,reverse_sweep_warp,32,0);
  *maxBlocksFwd = mbf;
  *maxBlocksRev = mbr;
}

/*******************************************************************************************
 *
 *  Task 8: Stage-2 SPARSE-CHECKPOINT trace-points.
 *
 *  A faithful warp-cooperative port of align.c's Local_Alignment (forward_wave + reverse_wave)
 *  that ALSO reproduces the Pebble cells[] checkpoints pushed DURING the sweep at each tspace
 *  boundary and the backward pointer-walk that emits (diff,Δb) per panel:
 *    - forward:  root/crossing pushes align.c:482-533 (0-wave), 746-768 (steady); backward walk
 *                & emit align.c:860-911.
 *    - reverse:  root(mark=x)/crossing pushes align.c:1010-1060 (0-wave), 1266-1288 (steady);
 *                backward walk, emit & seed-junction MERGE into the forward trace align.c:1366-1456.
 *
 *  Memory: O(checkpoints), NOT O(band*depth).  cells[] is a BOUNDED per-warp buffer (TR_CELLS_MAX);
 *  the per-seed trace is a bounded centred scratch (TR_TRACE_CAP).  Either overflowing marks the
 *  seed -2 (excluded by the validator, exactly like Stage-1's WV_SPAN band overflow).  Total scratch
 *  is capped by processing seeds over a fixed TR_POOL of concurrent warps (grid-stride).
 *
 *  DESIGN DECISION (disclosed): this is a SEPARATE kernel from the validated Stage-1 sweep kernels
 *  (forward_sweep_warp / reverse_sweep_warp), which are left BYTE-IDENTICAL -- so the Stage-1 gate
 *  cannot regress.  It is NOT a 2x re-sweep of Stage-1: it is self-contained (computes endpoints AND
 *  trace in its own single forward + single reverse pass, exactly as align.c does).  It reuses the
 *  SAME phase-1-parallel-slide / phase-2-serial-trim structure as Stage-1, so its endpoints match
 *  Stage-1's and the checkpoint pushes / trima tracking happen in the naturally-serial phase-2 pass
 *  (preserving align.c's exact avail ordering).  The only added per-diagonal state vs Stage-1 is the
 *  ping-ponged HA[]/NA[] (checkpoint-lineage index & next-boundary per diagonal) and a phase-1->2
 *  hand-off of the inherited lineage index (ha_src).
 *
 *******************************************************************************************/

#define TR_WV_SPAN    2048       // per-warp diagonal-band scratch (== Stage-1 WV_SPAN)
#define TR_CELLS_MAX  262144     // per-warp checkpoint capacity; overflow -> seed marked -2
#define TR_TRACE_CAP  32768      // per-seed centred trace scratch (uint16); TR_TRACE_CAP/2 each side
#define TR_POOL       512        // concurrent warps (bounds total device scratch)
#define TR_CENTER     (TR_TRACE_CAP/2)

//  Launch-grid (concurrency) cap.  Committed default is TR_POOL=512 (a self-imposed ~2 GB
//  checkpoint-scratch budget, NOT a hardware limit -- the kernel's own occupancy ceiling is
//  2160 warps).  Overridable at run time via env TR_POOL for the Stage-2 concurrency sweep
//  (measurement only -- the kernel algorithm is unchanged).  Returns the committed default
//  unless a positive TR_POOL is set in the environment.
static int tr_pool_cap(void)
{ static int cap = -1;
  if (cap < 0)
    { const char *e = getenv("TR_POOL");
      cap = (e != NULL && atoi(e) > 0) ? atoi(e) : TR_POOL;
    }
  return cap;
}

struct WvTrState
  { int low, hgh;
    int prevbase, prev_lo, prev_hi;
    int besta, bestx, lasta;
    int trima, trimx, trimd, trimha;
    int more, dif, overflow;
    int avail;                 // # checkpoints pushed so far (this direction)
    // emit outputs handed forward-pass -> reverse-pass:
    int fwd_tlen;
    int aepos, bepos, fdiff;   // forward endpoint + forward diff
  };

//  prev-wave furthest-reach reads for the trace state (mirror prevV/prevVr but for WvTrState).
__device__ __forceinline__ int prevV_tr(const int *Vb, const WvTrState &s, int k)
{ if (k < s.prev_lo || k > s.prev_hi) return -1; return Vb[k - s.prevbase]; }
__device__ __forceinline__ int prevVr_tr(const int *Vb, const WvTrState &s, int k)
{ if (k < s.prev_lo || k > s.prev_hi) return WV_INTMAX; return Vb[k - s.prevbase]; }

//  lane-0 checkpoint push; returns the new lineage index (== avail before increment), or -1 on
//  overflow (also sets s.overflow).  Mirrors "pb=cells+avail; ...; ha=avail++".
__device__ __forceinline__ int tr_push(WvTrState &s,
        int *cptr,int *cdiag,int *cdiff,int *cmark, int ptr,int diag,int diff,int mark)
{ int a = s.avail;
  if (a >= TR_CELLS_MAX) { s.overflow = 3; return -1; }   // 3 = cells-buffer overflow
  cptr[a]=ptr; cdiag[a]=diag; cdiff[a]=diff; cmark[a]=mark;
  s.avail = a+1;
  return a;
}

/*  One warp does the FULL forward+reverse trace for one seed (mirrors Local_Alignment). */
__global__ void wave_trace_warp(
        const u8 *A,    const long *Abase,
        const u8 *Bfwd, const u8 *Brev, const long *Bbase,
        int n, const wave_seed *seeds, int aoff, int tspace,
        int *o_ab,int *o_ae,int *o_bb,int *o_be,int *o_diffs,int *o_tlen,
        uint16 *o_trace, int trace_stride,
        int *Vpool,BVEC *Tpool,int *Mpool,int *HApool,int *NApool,
        unsigned char *Cpool,int *HSpool,
        int *CPptr,int *CPdiag,int *CPdiff,int *CPmark, uint16 *TRpool,
        const short *TABLE,const short *SCORE,int PATH_AVE)
{
  int lane = threadIdx.x;
  int warp = blockIdx.x;
  __shared__ WvTrState s;

  int  *V0 = Vpool + ((size_t)warp*2+0)*TR_WV_SPAN, *V1 = Vpool + ((size_t)warp*2+1)*TR_WV_SPAN;
  BVEC *T0 = Tpool + ((size_t)warp*2+0)*TR_WV_SPAN, *T1 = Tpool + ((size_t)warp*2+1)*TR_WV_SPAN;
  int  *M0 = Mpool + ((size_t)warp*2+0)*TR_WV_SPAN, *M1 = Mpool + ((size_t)warp*2+1)*TR_WV_SPAN;
  int  *H0 = HApool+ ((size_t)warp*2+0)*TR_WV_SPAN, *H1 = HApool+ ((size_t)warp*2+1)*TR_WV_SPAN;
  int  *N0 = NApool+ ((size_t)warp*2+0)*TR_WV_SPAN, *N1 = NApool+ ((size_t)warp*2+1)*TR_WV_SPAN;
  unsigned char *Cbuf = Cpool + (size_t)warp*TR_WV_SPAN;
  int  *HSbuf = HSpool + (size_t)warp*TR_WV_SPAN;
  int  *cptr = CPptr + (size_t)warp*TR_CELLS_MAX, *cdiag = CPdiag + (size_t)warp*TR_CELLS_MAX;
  int  *cdiff= CPdiff+ (size_t)warp*TR_CELLS_MAX, *cmark = CPmark + (size_t)warp*TR_CELLS_MAX;
  uint16 *trbuf = TRpool + (size_t)warp*TR_TRACE_CAP;

  for (int si = warp; si < n; si += gridDim.x)
    { const wave_seed sd = seeds[si];
      int mida = sd.anti;

      // ====================================================================================
      //  FORWARD  (align.c forward_wave: extends rightward, x+=1, retired sentinel V=-1)
      // ====================================================================================
      { const u8 *aseq = A + Abase[sd.aread];
        const u8 *bseq = (sd.comp ? Brev : Bfwd) + Bbase[sd.bread];

        int *curV=V0,*prevVb=V1; BVEC *curT=T0,*prevTb=T1; int *curM=M0,*prevMb=M1;
        int *curH=H0,*prevHb=H1; int *curN=N0,*prevNb=N1;

        if (lane==0)
          { int low=sd.diag, hgh=sd.diag;
            while (((mida-hgh)>>1) < 0) hgh -= 1;
            s.low=low; s.hgh=hgh;
            s.besta=s.trima=s.lasta=mida; s.bestx=s.trimx=(mida+hgh)>>1;
            s.trimd=0; s.trimha=0; s.more=1; s.dif=0; s.overflow=0; s.avail=0;
          }
        __syncwarp();
        int low=s.low, hgh=s.hgh, curbase=low;

        if (hgh-low+1 > TR_WV_SPAN) { if(lane==0) s.overflow=2; }
        else
        {
          //  ---- 0-wave slide (parallel) ----
          for (int k=low+lane; k<=hgh; k+=32)
            { int x=(mida+k)>>1, code=0;
              while (1)
                { int cB=bseq[x-k]; if(cB==4){code=2;break;}
                  int dA=aseq[x]; if(cB!=dA){ if(dA==4)code=1; break;} x+=1; }
              int c=(x<<1)-k, idx=k-curbase;
              curV[idx]=c; curT[idx]=PATH_INT; curM[idx]=PATH_LEN; Cbuf[idx]=(unsigned char)code;
            }
          __syncwarp();

          //  ---- 0-wave serial: root+crossing pushes, trim, set HA/NA (k descending) ----
          if (lane==0)
            { int aclip=WV_INTMAX, bclip=-WV_INTMAX;
              for (int k=hgh; k>=low; k--)
                { int idx=k-curbase;
                  int c=curV[idx]; int xf=(c+k)>>1; int code=Cbuf[idx];
                  int x0=(mida+k)>>1;                                   // seed x (pre-slide)
                  int na=((x0+(tspace-aoff))/tspace-1)*tspace+aoff;     // align.c:482
                  int ha=tr_push(s,cptr,cdiag,cdiff,cmark,-1,k,0,na);   // root (ptr=-1)
                  na+=tspace;
                  while (xf>=na && !s.overflow)                         // align.c:514-533
                    { ha=tr_push(s,cptr,cdiag,cdiff,cmark,ha,k,0,na); na+=tspace; }
                  if (code==2){ s.more=0; if(bclip<k)bclip=k; }
                  else if (code==1){ s.more=0; aclip=k; }
                  if (c>s.besta){ s.besta=s.trima=s.lasta=c; s.bestx=s.trimx=xf; s.trimha=ha; }
                  curH[idx]=ha; curN[idx]=na;
                }
              if (s.more==0)
                { if (bseq[s.besta-s.bestx]!=4 && aseq[s.bestx]!=4) s.more=1;
                  if (hgh>=aclip) hgh=aclip-1;
                  if (low<=bclip) low=bclip+1;
                }
              s.low=low; s.hgh=hgh;
              s.prevbase=curbase; s.prev_lo=low; s.prev_hi=hgh;
            }
          __syncwarp();
          low=s.low; hgh=s.hgh;
          { int *t;BVEC *tb; t=curV;curV=prevVb;prevVb=t; tb=curT;curT=prevTb;prevTb=tb;
            t=curM;curM=prevMb;prevMb=t; t=curH;curH=prevHb;prevHb=t; t=curN;curN=prevNb;prevNb=t; }

          //  ---- steady waves ----
          while (s.more && s.lasta >= s.besta-TRIM_MLAG && !s.overflow)
            { low=s.low-1; hgh=s.hgh+1; curbase=low;
              if (hgh-low+1 > TR_WV_SPAN){ if(lane==0)s.overflow=2; __syncwarp(); break; }
              int dif=s.dif+1;

              //  phase 1: parallel recurrence + slide; record ha_src of winning source
              for (int k=low+lane; k<=hgh; k+=32)
                { int ap=prevV_tr(prevVb,s,k+1), ac=prevV_tr(prevVb,s,k), am=prevV_tr(prevVb,s,k-1);
                  int c,m,jsrc; BVEC b;
                  if (ac<am){ if(am<ap){c=ap+1;jsrc=k+1;} else {c=am+1;jsrc=k-1;} }
                  else      { if(ac<ap){c=ap+1;jsrc=k+1;} else {c=ac+2;jsrc=k;} }
                  int inb=(jsrc>=s.prev_lo && jsrc<=s.prev_hi);
                  b=inb?prevTb[jsrc-s.prevbase]:PATH_INT;
                  m=inb?prevMb[jsrc-s.prevbase]:PATH_LEN;
                  int hs=inb?prevHb[jsrc-s.prevbase]:-1;
                  if((b&PATH_TOP)!=0)m-=1; b<<=1;
                  int x=(c+k)>>1, code=0;
                  while(1){ int cB=bseq[x-k]; if(cB==4){code=2;break;}
                    int dA=aseq[x]; if(cB!=dA){ if(dA==4)code=1; break;}
                    x+=1; if((b&PATH_TOP)==0)m+=1; b=(b<<1)|1; }
                  c=(x<<1)-k; int idx=k-curbase;
                  curV[idx]=c; curM[idx]=m; curT[idx]=b; Cbuf[idx]=(unsigned char)code; HSbuf[idx]=hs;
                }
              __syncwarp();

              //  phase 2: serial trim/pushes/clip/prune (k descending, align.c:680-831)
              if (lane==0)
                { int aclip=WV_INTMAX, bclip=-WV_INTMAX;
                  for (int k=hgh; k>=low; k--)
                    { int idx=k-curbase;
                      int c=curV[idx], m=curM[idx]; BVEC b=curT[idx];
                      int xf=(c+k)>>1, code=Cbuf[idx], ha=HSbuf[idx];
                      //  starting NA[k] carried from prev wave (edges inherit neighbour, align.c:656/663)
                      int na;
                      if (k < s.prev_lo)      na=prevNb[s.prev_lo - s.prevbase];
                      else if (k > s.prev_hi) na=prevNb[s.prev_hi - s.prevbase];
                      else                    na=prevNb[k - s.prevbase];
                      while (xf>=na && !s.overflow)                     // align.c:746-768
                        { if (cmark[ha] < na)
                            ha=tr_push(s,cptr,cdiag,cdiff,cmark,ha,k,dif,na);
                          na+=tspace;
                        }
                      if (code==2){ s.more=0; if(bclip<k)bclip=k; }
                      else if (code==1){ s.more=0; aclip=k; }
                      if (c>s.besta)
                        { s.besta=c; s.bestx=xf;
                          if (m>=PATH_AVE)
                            { s.lasta=c;
                              if (TABLE[b&TRIM_MASK]>=0)
                                if (TABLE[(b>>TRIM_LEN)&TRIM_MASK]+SCORE[b&TRIM_MASK]>=0)
                                  { s.trima=c; s.trimx=xf; s.trimd=dif; s.trimha=ha; }
                            }
                        }
                      curH[idx]=ha; curN[idx]=na;
                    }
                  s.dif=dif;
                  if (s.more==0)
                    { if (bseq[s.besta-s.bestx]!=4 && aseq[s.bestx]!=4) s.more=1;
                      if (hgh>=aclip) hgh=aclip-1;
                      if (low<=bclip) low=bclip+1;
                    }
                  { int nlag=s.besta-WAVE_LAG;
                    while (hgh>=low)
                      { if (curV[hgh-curbase]<nlag) hgh-=1;
                        else { while(curV[low-curbase]<nlag) low+=1; break; } }
                  }
                  s.low=low; s.hgh=hgh; s.prevbase=curbase; s.prev_lo=low; s.prev_hi=hgh;
                }
              __syncwarp();
              { int *t;BVEC *tb; t=curV;curV=prevVb;prevVb=t; tb=curT;curT=prevTb;prevTb=tb;
                t=curM;curM=prevMb;prevMb=t; t=curH;curH=prevHb;prevHb=t; t=curN;curN=prevNb;prevNb=t; }
            }

          //  ---- forward backward-walk + emit (align.c:860-911); REACH=0 so trimy=trima-trimx ----
          if (lane==0 && !s.overflow)
            { int a,b,k,h,d,e,atlen=0;
              int trimx=s.trimx, trimd=s.trimd, trimy=s.trima-s.trimx;
              a=-1;
              for (h=s.trimha; h>=0; h=b){ b=cptr[h]; cptr[h]=a; a=h; }
              h=a;
              k=cdiag[h]; b=(mida-k)>>1; e=0;
              for (h=cptr[h]; h>=0; h=cptr[h])
                { k=cdiag[h]; a=cmark[h]-k; d=cdiff[h];
                  if (atlen+2 > TR_CENTER){ s.overflow=4; break; }
                  trbuf[TR_CENTER+atlen++]=(uint16)(d-e);
                  trbuf[TR_CENTER+atlen++]=(uint16)(a-b);
                  b=a; e=d;
                }
              if (!s.overflow)
                { if (b+k != trimx)
                    { if (atlen+2>TR_CENTER) s.overflow=4;
                      else { trbuf[TR_CENTER+atlen++]=(uint16)(trimd-e);
                             trbuf[TR_CENTER+atlen++]=(uint16)(trimy-b); } }
                  else if (b != trimy)
                    { trbuf[TR_CENTER+atlen-1]=(uint16)(trbuf[TR_CENTER+atlen-1]+(trimy-b));
                      trbuf[TR_CENTER+atlen-2]=(uint16)(trbuf[TR_CENTER+atlen-2]+(trimd-e)); }
                }
              s.fwd_tlen=atlen; s.aepos=trimx; s.bepos=trimy; s.fdiff=trimd;
            }
          __syncwarp();
        }
      }

      // ====================================================================================
      //  REVERSE (align.c reverse_wave: extends leftward, x-=1, retired sentinel V=INT32_MAX)
      // ====================================================================================
      if (!s.overflow)
      { const u8 *aseq = A + Abase[sd.aread] - 1;
        const u8 *bseq = (sd.comp ? Brev : Bfwd) + Bbase[sd.bread] - 1;

        int *curV=V0,*prevVb=V1; BVEC *curT=T0,*prevTb=T1; int *curM=M0,*prevMb=M1;
        int *curH=H0,*prevHb=H1; int *curN=N0,*prevNb=N1;

        if (lane==0)
          { int low=sd.diag, hgh=sd.diag;
            s.low=low; s.hgh=hgh;
            s.besta=s.trima=s.lasta=mida; s.bestx=s.trimx=(mida+hgh)>>1;
            s.trimd=0; s.trimha=0; s.more=1; s.dif=0; s.avail=0;   // reuse cells (avail reset)
          }
        __syncwarp();
        int low=s.low, hgh=s.hgh, curbase=low;

        if (hgh-low+1 > TR_WV_SPAN) { if(lane==0) s.overflow=2; }
        else
        {
          //  ---- 0-wave slide (parallel, x-=1) ----
          for (int k=low+lane; k<=hgh; k+=32)
            { int x=(mida+k)>>1, code=0;
              while (1)
                { int cB=bseq[x-k]; if(cB==4){code=2;break;}
                  int dA=aseq[x]; if(cB!=dA){ if(dA==4)code=1; break;} x-=1; }
              int c=(x<<1)-k, idx=k-curbase;
              curV[idx]=c; curT[idx]=PATH_INT; curM[idx]=PATH_LEN; Cbuf[idx]=(unsigned char)code;
            }
          __syncwarp();

          //  ---- 0-wave serial (k ascending): root(mark=x)+crossing pushes, trim ----
          if (lane==0)
            { int aclip=-WV_INTMAX, bclip=WV_INTMAX;
              for (int k=low; k<=hgh; k++)
                { int idx=k-curbase;
                  int c=curV[idx]; int xf=(c+k)>>1; int code=Cbuf[idx];
                  int x0=(mida+k)>>1;
                  int na=((x0+(tspace-aoff)-1)/tspace-1)*tspace+aoff;   // align.c:1010
                  int ha=tr_push(s,cptr,cdiag,cdiff,cmark,-1,k,0,x0);   // root mark=x (align.c:1018)
                  while (xf<=na && !s.overflow)                         // align.c:1041-1060
                    { ha=tr_push(s,cptr,cdiag,cdiff,cmark,ha,k,0,na); na-=tspace; }
                  if (code==2){ s.more=0; if(bclip>k)bclip=k; }
                  else if (code==1){ s.more=0; aclip=k; }
                  if (c<s.besta){ s.besta=s.trima=s.lasta=c; s.bestx=s.trimx=xf; s.trimha=ha; }
                  curH[idx]=ha; curN[idx]=na;
                }
              if (s.more==0)
                { if (bseq[s.besta-s.bestx]!=4 && aseq[s.bestx]!=4) s.more=1;
                  if (low<=aclip) low=aclip+1;
                  if (hgh>=bclip) hgh=bclip-1;
                }
              s.low=low; s.hgh=hgh; s.prevbase=curbase; s.prev_lo=low; s.prev_hi=hgh;
            }
          __syncwarp();
          low=s.low; hgh=s.hgh;
          { int *t;BVEC *tb; t=curV;curV=prevVb;prevVb=t; tb=curT;curT=prevTb;prevTb=tb;
            t=curM;curM=prevMb;prevMb=t; t=curH;curH=prevHb;prevHb=t; t=curN;curN=prevNb;prevNb=t; }

          //  ---- steady waves (k ascending, MIN recurrence) ----
          while (s.more && s.lasta <= s.besta+TRIM_MLAG && !s.overflow)
            { low=s.low-1; hgh=s.hgh+1; curbase=low;
              if (hgh-low+1 > TR_WV_SPAN){ if(lane==0)s.overflow=2; __syncwarp(); break; }
              int dif=s.dif+1;

              for (int k=low+lane; k<=hgh; k+=32)
                { int ap=prevVr_tr(prevVb,s,k+1), ac=prevVr_tr(prevVb,s,k), am=prevVr_tr(prevVb,s,k-1);
                  int c,m,jsrc; BVEC b;
                  if (ac>ap){ if(ap>am){c=am-1;jsrc=k-1;} else {c=ap-1;jsrc=k+1;} }
                  else      { if(ac>am){c=am-1;jsrc=k-1;} else {c=ac-2;jsrc=k;} }
                  int inb=(jsrc>=s.prev_lo && jsrc<=s.prev_hi);
                  b=inb?prevTb[jsrc-s.prevbase]:PATH_INT;
                  m=inb?prevMb[jsrc-s.prevbase]:PATH_LEN;
                  int hs=inb?prevHb[jsrc-s.prevbase]:-1;
                  if((b&PATH_TOP)!=0)m-=1; b<<=1;
                  int x=(c+k)>>1, code=0;
                  while(1){ int cB=bseq[x-k]; if(cB==4){code=2;break;}
                    int dA=aseq[x]; if(cB!=dA){ if(dA==4)code=1; break;}
                    x-=1; if((b&PATH_TOP)==0)m+=1; b=(b<<1)|1; }
                  c=(x<<1)-k; int idx=k-curbase;
                  curV[idx]=c; curM[idx]=m; curT[idx]=b; Cbuf[idx]=(unsigned char)code; HSbuf[idx]=hs;
                }
              __syncwarp();

              if (lane==0)
                { int aclip=-WV_INTMAX, bclip=WV_INTMAX;
                  for (int k=low; k<=hgh; k++)
                    { int idx=k-curbase;
                      int c=curV[idx], m=curM[idx]; BVEC b=curT[idx];
                      int xf=(c+k)>>1, code=Cbuf[idx], ha=HSbuf[idx];
                      int na;
                      if (k < s.prev_lo)      na=prevNb[s.prev_lo - s.prevbase];
                      else if (k > s.prev_hi) na=prevNb[s.prev_hi - s.prevbase];
                      else                    na=prevNb[k - s.prevbase];
                      while (xf<=na && !s.overflow)                     // align.c:1266-1288
                        { if (cmark[ha] > na)
                            ha=tr_push(s,cptr,cdiag,cdiff,cmark,ha,k,dif,na);
                          na-=tspace;
                        }
                      if (code==2){ s.more=0; if(bclip>k)bclip=k; }
                      else if (code==1){ s.more=0; aclip=k; }
                      if (c<s.besta)
                        { s.besta=c; s.bestx=xf;
                          if (m>=PATH_AVE)
                            { s.lasta=c;
                              if (TABLE[b&TRIM_MASK]>=0)
                                if (TABLE[(b>>TRIM_LEN)&TRIM_MASK]+SCORE[b&TRIM_MASK]>=0)
                                  { s.trima=c; s.trimx=xf; s.trimd=dif; s.trimha=ha; }
                            }
                        }
                      curH[idx]=ha; curN[idx]=na;
                    }
                  s.dif=dif;
                  if (s.more==0)
                    { if (bseq[s.besta-s.bestx]!=4 && aseq[s.bestx]!=4) s.more=1;
                      if (low<=aclip) low=aclip+1;
                      if (hgh>=bclip) hgh=bclip-1;
                    }
                  { int nlag=s.besta+WAVE_LAG;
                    while (hgh>=low)
                      { if (curV[hgh-curbase]>nlag) hgh-=1;
                        else { while(curV[low-curbase]>nlag) low+=1; break; } }
                  }
                  s.low=low; s.hgh=hgh; s.prevbase=curbase; s.prev_lo=low; s.prev_hi=hgh;
                }
              __syncwarp();
              { int *t;BVEC *tb; t=curV;curV=prevVb;prevVb=t; tb=curT;curT=prevTb;prevTb=tb;
                t=curM;curM=prevMb;prevMb=t; t=curH;curH=prevHb;prevHb=t; t=curN;curN=prevNb;prevNb=t; }
            }

          //  ---- reverse backward-walk + emit + junction merge (align.c:1366-1456) ----
          if (lane==0 && !s.overflow)
            { int a,b,k,h,d,e,atlen=0;
              int trimx=s.trimx, trimd=s.trimd, trimy=s.trima-s.trimx;
              int fwd_tlen=s.fwd_tlen;
              uint16 *atrace = trbuf + TR_CENTER;
              a=-1;
              for (h=s.trimha; h>=0; h=b){ b=cptr[h]; cptr[h]=a; a=h; }
              h=a;
              k=cdiag[h]; b=cmark[h]-k; e=0;
              if ((b+k)%tspace != aoff)
                { h=cptr[h];
                  if (h<0){ a=trimy; d=trimd; }
                  else { k=cdiag[h]; a=cmark[h]-k; d=cdiff[h]; }
                  if (fwd_tlen==0)
                    { if (TR_CENTER+atlen-2 < 0){ s.overflow=4; }
                      else { atrace[--atlen]=(uint16)(b-a); atrace[--atlen]=(uint16)(d-e); } }
                  else
                    { atrace[1]=(uint16)(atrace[1]+(b-a)); atrace[0]=(uint16)(atrace[0]+(d-e)); }
                  b=a; e=d;
                }
              if (h>=0 && !s.overflow)
                { for (h=cptr[h]; h>=0; h=cptr[h])
                    { k=cdiag[h]; a=cmark[h]-k;
                      if (TR_CENTER+atlen-2 < 0){ s.overflow=4; break; }
                      atrace[--atlen]=(uint16)(b-a);
                      d=cdiff[h];
                      atrace[--atlen]=(uint16)(d-e);
                      b=a; e=d;
                    }
                  if (!s.overflow)
                    { if (b+k != trimx)
                        { if (TR_CENTER+atlen-2 < 0) s.overflow=4;
                          else { atrace[--atlen]=(uint16)(b-trimy);
                                 atrace[--atlen]=(uint16)(trimd-e); } }
                      else if (b != trimy)
                        { atrace[atlen+1]=(uint16)(atrace[atlen+1]+(b-trimy));
                          atrace[atlen]  =(uint16)(atrace[atlen]  +(trimd-e)); }
                    }
                }
              if (!s.overflow)
                { int tlen = s.fwd_tlen - atlen;      // atlen<=0 (grew backward)
                  int start = TR_CENTER + atlen;
                  s.trimx=trimx; s.trimd=trimd;       // reverse endpoint stash
                  //  write outputs
                  o_ab[si]=trimx; o_ae[si]=s.aepos; o_bb[si]=trimy; o_be[si]=s.bepos;
                  o_diffs[si]=s.fdiff+trimd; o_tlen[si]=tlen;
                  if (tlen > trace_stride) { s.overflow=4; }   // 4 = trace-buffer overflow
                  else
                    { uint16 *dst = o_trace + (size_t)si*trace_stride;
                      for (int q=0;q<tlen;q++) dst[q]=trbuf[start+q];
                    }
                }
            }
          __syncwarp();
        }
      }

      if (lane==0 && s.overflow)
        { //  endpoints/diffs carry -2; o_tlen carries the overflow CAUSE (-2 band, -3 cells,
          //  -4 trace) so the validator can attribute it.
          o_ab[si]=o_ae[si]=o_bb[si]=o_be[si]=o_diffs[si]=-2;
          o_tlen[si] = -s.overflow;
        }
      __syncwarp();
    }
}

extern "C" int wave_trace_batch(wave_ctx *g, int n, const wave_seed *seeds,
                     const int16_t *table, const int16_t *score, int path_ave,
                     int aoff, int tspace,
                     int *ab, int *ae, int *bb, int *be, int *diffs, int *tlen,
                     uint16_t *trace, int trace_stride)
{ if (g == NULL || n <= 0) return 0;

  wave_seed *dSeeds=NULL;
  int *dAb=NULL,*dAe=NULL,*dBb=NULL,*dBe=NULL,*dDf=NULL,*dTl=NULL;
  uint16 *dTrace=NULL;
  short *dTab=NULL,*dScr=NULL;
  int *dV=NULL,*dM=NULL,*dH=NULL,*dN=NULL,*dHS=NULL; BVEC *dT=NULL; unsigned char *dC=NULL;
  int *dCp=NULL,*dCd=NULL,*dCf=NULL,*dCm=NULL; uint16 *dTR=NULL;

  int tr_pool = tr_pool_cap();
  int pool = n < tr_pool ? n : tr_pool;
  size_t vspan = (size_t) pool * 2 * TR_WV_SPAN;
  size_t cells = (size_t) pool * TR_CELLS_MAX;

  CK(cudaMalloc(&dSeeds,(size_t)n*sizeof(wave_seed)));
  CK(cudaMalloc(&dAb,(size_t)n*sizeof(int))); CK(cudaMalloc(&dAe,(size_t)n*sizeof(int)));
  CK(cudaMalloc(&dBb,(size_t)n*sizeof(int))); CK(cudaMalloc(&dBe,(size_t)n*sizeof(int)));
  CK(cudaMalloc(&dDf,(size_t)n*sizeof(int))); CK(cudaMalloc(&dTl,(size_t)n*sizeof(int)));
  CK(cudaMalloc(&dTrace,(size_t)n*trace_stride*sizeof(uint16)));
  CK(cudaMalloc(&dTab,(size_t)(TRIM_MASK+1)*sizeof(short)));
  CK(cudaMalloc(&dScr,(size_t)(TRIM_MASK+1)*sizeof(short)));
  CK(cudaMalloc(&dV,vspan*sizeof(int))); CK(cudaMalloc(&dT,vspan*sizeof(BVEC)));
  CK(cudaMalloc(&dM,vspan*sizeof(int))); CK(cudaMalloc(&dH,vspan*sizeof(int)));
  CK(cudaMalloc(&dN,vspan*sizeof(int)));
  CK(cudaMalloc(&dC,(size_t)pool*TR_WV_SPAN*sizeof(unsigned char)));
  CK(cudaMalloc(&dHS,(size_t)pool*TR_WV_SPAN*sizeof(int)));
  CK(cudaMalloc(&dCp,cells*sizeof(int))); CK(cudaMalloc(&dCd,cells*sizeof(int)));
  CK(cudaMalloc(&dCf,cells*sizeof(int))); CK(cudaMalloc(&dCm,cells*sizeof(int)));
  CK(cudaMalloc(&dTR,(size_t)pool*TR_TRACE_CAP*sizeof(uint16)));
  if (dSeeds==NULL||dAb==NULL||dTrace==NULL||dV==NULL||dCp==NULL||dTR==NULL)
    { fprintf(stderr,"wave_trace_batch: device alloc failed\n"); return -1; }

  CK(cudaMemcpy(dSeeds,seeds,(size_t)n*sizeof(wave_seed),cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dTab,table,(size_t)(TRIM_MASK+1)*sizeof(short),cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dScr,score,(size_t)(TRIM_MASK+1)*sizeof(short),cudaMemcpyHostToDevice));

  wave_trace_warp<<<pool,32>>>(g->dA,g->dAbase,g->dBfwd,g->dBrev,g->dBbase,
        n,dSeeds,aoff,tspace, dAb,dAe,dBb,dBe,dDf,dTl,dTrace,trace_stride,
        dV,dT,dM,dH,dN,dC,dHS,dCp,dCd,dCf,dCm,dTR,dTab,dScr,path_ave);
  cudaError_t err = cudaGetLastError();
  if (err==cudaSuccess) err = cudaDeviceSynchronize();
  if (err!=cudaSuccess)
    { fprintf(stderr,"wave_trace_batch: kernel error: %s\n",cudaGetErrorString(err)); return -1; }

  CK(cudaMemcpy(ab,dAb,(size_t)n*sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(ae,dAe,(size_t)n*sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(bb,dBb,(size_t)n*sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(be,dBe,(size_t)n*sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(diffs,dDf,(size_t)n*sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(tlen,dTl,(size_t)n*sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(trace,dTrace,(size_t)n*trace_stride*sizeof(uint16),cudaMemcpyDeviceToHost));

  cudaFree(dSeeds); cudaFree(dAb); cudaFree(dAe); cudaFree(dBb); cudaFree(dBe);
  cudaFree(dDf); cudaFree(dTl); cudaFree(dTrace); cudaFree(dTab); cudaFree(dScr);
  cudaFree(dV); cudaFree(dT); cudaFree(dM); cudaFree(dH); cudaFree(dN);
  cudaFree(dC); cudaFree(dHS); cudaFree(dCp); cudaFree(dCd); cudaFree(dCf); cudaFree(dCm); cudaFree(dTR);
  return 0;
}

/*  Task 9: CUDA-event-timed trace-inclusive discovery.  Body == wave_trace_batch above, with the
 *  four timing windows (seed H2D / kernel / meta D2H / strided-trace D2H) dropped in exactly as
 *  wave_discover_batch_timed does for Stage-1.  See wave_kernel.h for the honesty contract on the
 *  strided-trace-D2H artifact. */
extern "C" int wave_trace_batch_timed(wave_ctx *g, int n, const wave_seed *seeds,
                     const int16_t *table, const int16_t *score, int path_ave,
                     int aoff, int tspace,
                     int *ab, int *ae, int *bb, int *be, int *diffs, int *tlen,
                     uint16_t *trace, int trace_stride,
                     float *ms_h2d, float *ms_kernel, float *ms_d2h_meta, float *ms_d2h_trace)
{ if (g == NULL || n <= 0) return 0;

  wave_seed *dSeeds=NULL;
  int *dAb=NULL,*dAe=NULL,*dBb=NULL,*dBe=NULL,*dDf=NULL,*dTl=NULL;
  uint16 *dTrace=NULL;
  short *dTab=NULL,*dScr=NULL;
  int *dV=NULL,*dM=NULL,*dH=NULL,*dN=NULL,*dHS=NULL; BVEC *dT=NULL; unsigned char *dC=NULL;
  int *dCp=NULL,*dCd=NULL,*dCf=NULL,*dCm=NULL; uint16 *dTR=NULL;

  int tr_pool = tr_pool_cap();
  int pool = n < tr_pool ? n : tr_pool;
  size_t vspan = (size_t) pool * 2 * TR_WV_SPAN;
  size_t cells = (size_t) pool * TR_CELLS_MAX;

  CK(cudaMalloc(&dSeeds,(size_t)n*sizeof(wave_seed)));
  CK(cudaMalloc(&dAb,(size_t)n*sizeof(int))); CK(cudaMalloc(&dAe,(size_t)n*sizeof(int)));
  CK(cudaMalloc(&dBb,(size_t)n*sizeof(int))); CK(cudaMalloc(&dBe,(size_t)n*sizeof(int)));
  CK(cudaMalloc(&dDf,(size_t)n*sizeof(int))); CK(cudaMalloc(&dTl,(size_t)n*sizeof(int)));
  CK(cudaMalloc(&dTrace,(size_t)n*trace_stride*sizeof(uint16)));
  CK(cudaMalloc(&dTab,(size_t)(TRIM_MASK+1)*sizeof(short)));
  CK(cudaMalloc(&dScr,(size_t)(TRIM_MASK+1)*sizeof(short)));
  CK(cudaMalloc(&dV,vspan*sizeof(int))); CK(cudaMalloc(&dT,vspan*sizeof(BVEC)));
  CK(cudaMalloc(&dM,vspan*sizeof(int))); CK(cudaMalloc(&dH,vspan*sizeof(int)));
  CK(cudaMalloc(&dN,vspan*sizeof(int)));
  CK(cudaMalloc(&dC,(size_t)pool*TR_WV_SPAN*sizeof(unsigned char)));
  CK(cudaMalloc(&dHS,(size_t)pool*TR_WV_SPAN*sizeof(int)));
  CK(cudaMalloc(&dCp,cells*sizeof(int))); CK(cudaMalloc(&dCd,cells*sizeof(int)));
  CK(cudaMalloc(&dCf,cells*sizeof(int))); CK(cudaMalloc(&dCm,cells*sizeof(int)));
  CK(cudaMalloc(&dTR,(size_t)pool*TR_TRACE_CAP*sizeof(uint16)));
  if (dSeeds==NULL||dAb==NULL||dTrace==NULL||dV==NULL||dCp==NULL||dTR==NULL)
    { fprintf(stderr,"wave_trace_batch_timed: device alloc failed\n"); return -1; }

  //  one-time-in-a-real-deployment tables -- outside the timed windows, as in Stage-1
  CK(cudaMemcpy(dTab,table,(size_t)(TRIM_MASK+1)*sizeof(short),cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dScr,score,(size_t)(TRIM_MASK+1)*sizeof(short),cudaMemcpyHostToDevice));

  cudaEvent_t e0,e1,e2,e3,e4;
  cudaEventCreate(&e0); cudaEventCreate(&e1); cudaEventCreate(&e2);
  cudaEventCreate(&e3); cudaEventCreate(&e4);

  //  ---- window 1: per-batch seed H2D ----
  cudaEventRecord(e0);
  CK(cudaMemcpy(dSeeds,seeds,(size_t)n*sizeof(wave_seed),cudaMemcpyHostToDevice));
  cudaEventRecord(e1);

  //  ---- window 2: kernel-only (single-pass discover+trace) ----
  wave_trace_warp<<<pool,32>>>(g->dA,g->dAbase,g->dBfwd,g->dBrev,g->dBbase,
        n,dSeeds,aoff,tspace, dAb,dAe,dBb,dBe,dDf,dTl,dTrace,trace_stride,
        dV,dT,dM,dH,dN,dC,dHS,dCp,dCd,dCf,dCm,dTR,dTab,dScr,path_ave);
  cudaEventRecord(e2);

  cudaError_t err = cudaGetLastError();
  if (err==cudaSuccess) err = cudaEventSynchronize(e2);
  if (err!=cudaSuccess)
    { fprintf(stderr,"wave_trace_batch_timed: kernel error: %s\n",cudaGetErrorString(err));
      cudaEventDestroy(e0); cudaEventDestroy(e1); cudaEventDestroy(e2);
      cudaEventDestroy(e3); cudaEventDestroy(e4);
      cudaFree(dSeeds); cudaFree(dAb); cudaFree(dAe); cudaFree(dBb); cudaFree(dBe);
      cudaFree(dDf); cudaFree(dTl); cudaFree(dTrace); cudaFree(dTab); cudaFree(dScr);
      cudaFree(dV); cudaFree(dT); cudaFree(dM); cudaFree(dH); cudaFree(dN);
      cudaFree(dC); cudaFree(dHS); cudaFree(dCp); cudaFree(dCd); cudaFree(dCf); cudaFree(dCm); cudaFree(dTR);
      return -1;
    }

  //  ---- window 3: meta D2H (6 small result arrays) ----
  cudaEventRecord(e2);   // re-mark start of D2H precisely at kernel-end (already synced)
  CK(cudaMemcpy(ab,dAb,(size_t)n*sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(ae,dAe,(size_t)n*sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(bb,dBb,(size_t)n*sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(be,dBe,(size_t)n*sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(diffs,dDf,(size_t)n*sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(tlen,dTl,(size_t)n*sizeof(int),cudaMemcpyDeviceToHost));
  cudaEventRecord(e3);

  //  ---- window 4: strided trace D2H (padding-dominated artifact -- see header) ----
  CK(cudaMemcpy(trace,dTrace,(size_t)n*trace_stride*sizeof(uint16),cudaMemcpyDeviceToHost));
  cudaEventRecord(e4);
  cudaEventSynchronize(e4);

  float t_h2d=0,t_kernel=0,t_meta=0,t_trace=0;
  cudaEventElapsedTime(&t_h2d,   e0,e1);
  cudaEventElapsedTime(&t_kernel,e1,e2);
  cudaEventElapsedTime(&t_meta,  e2,e3);
  cudaEventElapsedTime(&t_trace, e3,e4);
  if (ms_h2d)       *ms_h2d       = t_h2d;
  if (ms_kernel)    *ms_kernel    = t_kernel;
  if (ms_d2h_meta)  *ms_d2h_meta  = t_meta;
  if (ms_d2h_trace) *ms_d2h_trace = t_trace;

  cudaEventDestroy(e0); cudaEventDestroy(e1); cudaEventDestroy(e2);
  cudaEventDestroy(e3); cudaEventDestroy(e4);
  cudaFree(dSeeds); cudaFree(dAb); cudaFree(dAe); cudaFree(dBb); cudaFree(dBe);
  cudaFree(dDf); cudaFree(dTl); cudaFree(dTrace); cudaFree(dTab); cudaFree(dScr);
  cudaFree(dV); cudaFree(dT); cudaFree(dM); cudaFree(dH); cudaFree(dN);
  cudaFree(dC); cudaFree(dHS); cudaFree(dCp); cudaFree(dCd); cudaFree(dCf); cudaFree(dCm); cudaFree(dTR);
  return 0;
}

/*  Task 9: occupancy readout for the trace kernel (see wave_kernel.h contract). */
extern "C" void wave_query_trace_occupancy(int *maxBlocksTrace, int *smCount, int *gridCap)
{ int dev = 0;
  cudaDeviceProp prop;
  cudaGetDevice(&dev);
  cudaGetDeviceProperties(&prop,dev);
  *smCount = prop.multiProcessorCount;
  int mbt = 0;
  cudaOccupancyMaxActiveBlocksPerMultiprocessor(&mbt,wave_trace_warp,32,0);
  *maxBlocksTrace = mbt;
  *gridCap = tr_pool_cap();
}
