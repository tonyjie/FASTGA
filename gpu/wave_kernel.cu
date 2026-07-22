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

typedef unsigned char u8;

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
