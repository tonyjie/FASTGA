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
