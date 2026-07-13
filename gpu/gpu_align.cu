/* FastGA GPU alignment experiment — batched banded unit-cost (Levenshtein) edit
#include <utility>
#include <algorithm>
 * distance on A100, validated against FastGA's reference `diffs` and benchmarked
 * against a 32-thread CPU run of the SAME algorithm.
 *
 * Kernel target: the per-alignment work of FastGA's Phase-3 Myers O(nd) local wave.
 * We extract each alignment's two aligned segments (via extract_tasks) and recompute
 * the unit-cost edit distance (substitution=1, indel=1 — matching FastGA's wave,
 * whose V[k]+2 diagonal move is a cost-1 substitution). One CUDA thread per alignment.
 *
 * Build:  nvcc -O3 -arch=sm_80 -Xcompiler -fopenmp -o gpu_align gpu_align.cu
 * Run:    ./gpu_align <tasks-file>
 */
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <chrono>
#include <omp.h>
#include <cuda_runtime.h>
#include "task_format.h"

// band window + bounds.  M2 measured band < 256 across all divergence; KBAND gives margin.
#define MAXW  1024     // max k-window width (2*KBAND + |m-n|); overflow -> -2
#define KBAND 384
#define MAXD  8192     // max differences before giving up -> -1

typedef unsigned char u8;

#define CK(call) do{ cudaError_t e=(call); if(e!=cudaSuccess){ \
  fprintf(stderr,"CUDA error %s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(e)); exit(1);} }while(0)

// ---- shared edit-distance kernel (host + device) --------------------------
// Levenshtein edit distance between A[0..m) and B[0..n) via furthest-reaching
// O(ND) wave, restricted to a k-window around the m-n diagonal.
// Returns: edit distance d (>=0), or -1 (d>MAXD), or -2 (window>MAXW).
__host__ __device__ static int editdist(const u8* A, int m, const u8* B, int n)
{
  int drift = m - n;
  int klo = (drift < 0 ? drift : 0) - KBAND;
  int khi = (drift > 0 ? drift : 0) + KBAND;
  int width = khi - klo + 1;
  if (width > MAXW) return -2;

  short f[2][MAXW];
  int cur = 0, prv = 1;
  for (int i = 0; i < width; i++) f[prv][i] = -1;

  // d = 0: only diagonal k=0 reachable; slide matches from (0,0)
  { int x = 0;
    while (x < m && x < n && A[x] == B[x]) x++;
    for (int i = 0; i < width; i++) f[cur][i] = -1;
    f[cur][0 - klo] = (short) x;
    if (drift == 0 && x >= m) return 0;
    int t = cur; cur = prv; prv = t;              // swap: prv now holds d=0
  }

  for (int d = 1; d <= MAXD; d++)
    { for (int i = 0; i < width; i++) f[cur][i] = -1;
      int kmin = klo > -d ? klo : -d;
      int kmax = khi <  d ? khi :  d;
      for (int k = kmin; k <= kmax; k++)
        { int best = -1, v;
          if (k-1 >= klo) { v = f[prv][k-1-klo]; if (v >= 0 && v+1 > best) best = v+1; } // horizontal
          if (k+1 <= khi) { v = f[prv][k+1-klo]; if (v >= 0 && v   > best) best = v;   } // vertical
          v = f[prv][k-klo];      if (v >= 0 && v+1 > best) best = v+1;                  // substitution
          if (best < 0) continue;
          int x = best;
          if (x > m) x = m;
          if (x - k > n) x = n + k;     // clip so y = x-k <= n
          if (x < 0) x = 0;
          int y = x - k;
          while (x < m && y < n && A[x] == B[y]) { x++; y++; }
          f[cur][k-klo] = (short) x;
        }
      int xe = f[cur][drift-klo];
      if (xe >= m && xe - drift >= n) return d;
      int t = cur; cur = prv; prv = t;
    }
  return -1;
}

// ---- GPU kernel: one WARP per task (band across lanes, DP in shared mem) ---
#define WWIDTH 1024    // shared k-window per warp; task overflows if width exceeds
__global__ void align_kernel_warp(const u8* Abuf, const u8* Bbuf,
                                  const int* Aoff, const int* Boff,
                                  const int* Mlen, const int* Nlen,
                                  int* out, int ntasks)
{
  extern __shared__ short sh[];                 // (blockDim.x/32) * 2 * WWIDTH
  int warp = (blockIdx.x*blockDim.x + threadIdx.x) >> 5;   // one warp == one task
  int lane = threadIdx.x & 31;
  int wib  = threadIdx.x >> 5;
  if (warp >= ntasks) return;
  const u8* A = Abuf + Aoff[warp]; int m = Mlen[warp];
  const u8* B = Bbuf + Boff[warp]; int n = Nlen[warp];
  int drift = m - n;
  int klo = (drift<0?drift:0) - KBAND;
  int khi = (drift>0?drift:0) + KBAND;
  int width = khi - klo + 1;
  if (width > WWIDTH) { if (lane==0) out[warp] = -2; return; }
  short* fp = sh + wib*2*WWIDTH;
  short* fc = fp + WWIDTH;
  for (int i = lane; i < width; i += 32) fp[i] = -1;
  __syncwarp();
  if (lane == 0) { int x=0; while (x<m && x<n && A[x]==B[x]) x++; fp[0-klo] = (short)x; }
  __syncwarp();
  if (drift == 0 && fp[0-klo] >= m) { if (lane==0) out[warp]=0; return; }
  int result = -1;
  for (int d = 1; d <= MAXD; d++)
    { for (int i = lane; i < width; i += 32) fc[i] = -1;
      __syncwarp();
      int kmin = klo>-d?klo:-d, kmax = khi<d?khi:d;
      for (int k = kmin+lane; k <= kmax; k += 32)
        { int idx=k-klo, best=-1, v;
          if (k-1>=klo){ v=fp[idx-1]; if(v>=0&&v+1>best)best=v+1; }
          if (k+1<=khi){ v=fp[idx+1]; if(v>=0&&v  >best)best=v;   }
          v=fp[idx];     if(v>=0&&v+1>best)best=v+1;
          if (best>=0){ int x=best; if(x>m)x=m; if(x-k>n)x=n+k; if(x<0)x=0; int y=x-k;
            while (x<m && y<n && A[x]==B[y]){x++;y++;} fc[idx]=(short)x; }
        }
      __syncwarp();
      int xe = fc[drift-klo];
      if (xe>=m && xe-drift>=n) { result = d; break; }
      short* t=fp; fp=fc; fc=t;
    }
  if (lane==0) out[warp] = result;
}

// ---- GPU kernel: one thread per task --------------------------------------
__global__ void align_kernel(const u8* Abuf, const u8* Bbuf,
                             const int* Aoff, const int* Boff,
                             const int* Mlen, const int* Nlen,
                             int* out, int ntasks)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= ntasks) return;
  out[i] = editdist(Abuf + Aoff[i], Mlen[i], Bbuf + Boff[i], Nlen[i]);
}

// ---- task loading ---------------------------------------------------------
struct Tasks {
  int ntasks;
  std::vector<u8>  Abuf, Bbuf;
  std::vector<int> Aoff, Boff, Mlen, Nlen, Ref;
};

static void load(const char* path, Tasks& T)
{
  FILE* f = fopen(path, "rb");
  if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
  FGATaskHeader h;
  if (fread(&h, sizeof(h), 1, f) != 1 || h.magic != FGA_TASK_MAGIC)
    { fprintf(stderr, "bad task file / magic\n"); exit(1); }
  T.ntasks = h.ntasks;
  printf("tasks: %u   wmax: %u\n", h.ntasks, h.wmax);
  T.Aoff.reserve(h.ntasks); T.Boff.reserve(h.ntasks);
  T.Mlen.reserve(h.ntasks); T.Nlen.reserve(h.ntasks); T.Ref.reserve(h.ntasks);
  for (uint32_t i = 0; i < h.ntasks; i++)
    { uint32_t m, n, diffs;
      if (fread(&m,4,1,f)!=1 || fread(&n,4,1,f)!=1 || fread(&diffs,4,1,f)!=1) { fprintf(stderr,"trunc\n"); exit(1); }
      T.Aoff.push_back((int)T.Abuf.size()); T.Mlen.push_back((int)m);
      T.Boff.push_back((int)T.Bbuf.size()); T.Nlen.push_back((int)n);
      T.Ref.push_back((int)diffs);
      size_t as=T.Abuf.size(); T.Abuf.resize(as+m); if(m&&fread(&T.Abuf[as],1,m,f)!=m){fprintf(stderr,"truncA\n");exit(1);}
      size_t bs=T.Bbuf.size(); T.Bbuf.resize(bs+n); if(n&&fread(&T.Bbuf[bs],1,n,f)!=n){fprintf(stderr,"truncB\n");exit(1);}
    }
  fclose(f);
}

// brute-force full Levenshtein DP (correctness oracle for the self-test)
static int brute(const u8*A,int m,const u8*B,int n){
  std::vector<int> prev(n+1), cur(n+1);
  for(int j=0;j<=n;j++) prev[j]=j;
  for(int i=1;i<=m;i++){ cur[0]=i;
    for(int j=1;j<=n;j++){ int s=prev[j-1]+(A[i-1]!=B[j-1]);
      int d=prev[j]+1, ins=cur[j-1]+1; cur[j]=s<d?(s<ins?s:ins):(d<ins?d:ins);}
    std::swap(prev,cur);}
  return prev[n];
}

int main(int argc, char** argv)
{
  if (argc < 2) { fprintf(stderr,"usage: %s <tasks-file>\n",argv[0]); return 1; }

  // --- self-test editdist vs brute force on random small pairs ---
  { srand(12345); int bad=0;
    for(int t=0;t<2000;t++){ int m=rand()%80+1, n=rand()%80+1; u8 A[80],B[80];
      for(int i=0;i<m;i++)A[i]=rand()&3; for(int i=0;i<n;i++)B[i]=rand()&3;
      int e=editdist(A,m,B,n), b=brute(A,m,B,n);
      if(e!=b){ if(bad<5) fprintf(stderr,"SELFTEST mismatch m=%d n=%d ed=%d brute=%d\n",m,n,e,b); bad++; } }
    if(bad){ fprintf(stderr,"SELF-TEST FAILED (%d/2000)\n",bad); return 1; }
    printf("self-test: editdist == brute-force Levenshtein on 2000 random pairs  OK\n"); }

  Tasks T; load(argv[1], T);
  int N = T.ntasks;
  std::vector<int> cpuRes(N), gpuRes(N);

  // --- CPU: 32 threads, same editdist ---
  int nth = 32; omp_set_num_threads(nth);
  auto c0 = std::chrono::high_resolution_clock::now();
  #pragma omp parallel for schedule(dynamic,256)
  for (int i = 0; i < N; i++)
    cpuRes[i] = editdist(&T.Abuf[T.Aoff[i]], T.Mlen[i], &T.Bbuf[T.Boff[i]], T.Nlen[i]);
  auto c1 = std::chrono::high_resolution_clock::now();
  double cpu_s = std::chrono::duration<double>(c1-c0).count();

  // --- GPU ---
  u8 *dA,*dB; int *dAo,*dBo,*dM,*dN,*dO;
  CK(cudaMalloc(&dA,T.Abuf.size())); CK(cudaMalloc(&dB,T.Bbuf.size()));
  CK(cudaMalloc(&dAo,N*4)); CK(cudaMalloc(&dBo,N*4)); CK(cudaMalloc(&dM,N*4)); CK(cudaMalloc(&dN,N*4)); CK(cudaMalloc(&dO,N*4));
  auto t0 = std::chrono::high_resolution_clock::now();
  CK(cudaMemcpy(dA,T.Abuf.data(),T.Abuf.size(),cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dB,T.Bbuf.data(),T.Bbuf.size(),cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dAo,T.Aoff.data(),N*4,cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dBo,T.Boff.data(),N*4,cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dM,T.Mlen.data(),N*4,cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dN,T.Nlen.data(),N*4,cudaMemcpyHostToDevice));
  auto t1 = std::chrono::high_resolution_clock::now();
  int TPB=64, blocks=(N+TPB-1)/TPB;
  align_kernel<<<blocks,TPB>>>(dA,dB,dAo,dBo,dM,dN,dO,N);
  CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
  auto t2 = std::chrono::high_resolution_clock::now();
  CK(cudaMemcpy(gpuRes.data(),dO,N*4,cudaMemcpyDeviceToHost));
  auto t3 = std::chrono::high_resolution_clock::now();
  double gpu_h2d = std::chrono::duration<double>(t1-t0).count();
  double gpu_ker = std::chrono::duration<double>(t2-t1).count();
  double gpu_tot = std::chrono::duration<double>(t3-t0).count();

  // --- GPU warp-cooperative (band across lanes, DP in shared) ---
  std::vector<int> gpuResW(N);
  int wpb = 8, TPBw = wpb*32, blocksW = (N + wpb - 1)/wpb;
  size_t shmem = (size_t)wpb * 2 * WWIDTH * sizeof(short);
  CK(cudaFuncSetAttribute(align_kernel_warp, cudaFuncAttributeMaxDynamicSharedMemorySize, shmem));
  auto w0 = std::chrono::high_resolution_clock::now();
  align_kernel_warp<<<blocksW,TPBw,shmem>>>(dA,dB,dAo,dBo,dM,dN,dO,N);
  CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
  auto w1 = std::chrono::high_resolution_clock::now();
  CK(cudaMemcpy(gpuResW.data(),dO,N*4,cudaMemcpyDeviceToHost));
  double gpu_warp = std::chrono::duration<double>(w1-w0).count();

  // --- validate ---
  long same_cg=0, eq_ref=0, le_ref=0, ovf=0, wsame=0, weq=0, wle=0, wovf=0;
  for (int i=0;i<N;i++){
    if (gpuRes[i]==cpuRes[i]) same_cg++;
    if (gpuRes[i]<0) ovf++;
    else { if (gpuRes[i]==T.Ref[i]) eq_ref++; if (gpuRes[i]<=T.Ref[i]) le_ref++; }
    if (gpuResW[i]==cpuRes[i]) wsame++;
    if (gpuResW[i]<0) wovf++;
    else { if (gpuResW[i]==T.Ref[i]) weq++; if (gpuResW[i]<=T.Ref[i]) wle++; }
  }
  printf("\n=== correctness ===\n");
  printf("naive GPU == CPU (sanity)      : %.4f%%\n",100.0*same_cg/N);
  printf("warp  GPU == CPU (sanity)      : %.4f%%\n",100.0*wsame/N);
  printf("warp  GPU == FastGA diffs      : %.4f%%\n",100.0*weq/N);
  printf("warp  GPU <= FastGA diffs      : %.4f%%   (equal-or-better quality)\n",100.0*wle/N);
  printf("warp  overflow (band/d)        : %.4f%%  (%ld)\n",100.0*wovf/N,wovf);
  printf("\n=== performance (N=%d alignments) ===\n",N);
  printf("CPU 32-thread (same editdist)  : %.3f s   (%.3f M aln/s)   1.0x\n",cpu_s,N/cpu_s/1e6);
  printf("GPU naive (1 thread/aln) kernel: %.3f s   (%.3f M aln/s)   %.2fx\n",gpu_ker,N/gpu_ker/1e6,cpu_s/gpu_ker);
  printf("GPU warp  (1 warp/aln)   kernel: %.3f s   (%.3f M aln/s)   %.2fx\n",gpu_warp,N/gpu_warp/1e6,cpu_s/gpu_warp);
  printf("GPU warp  incl H2D+D2H         : %.3f s   (%.3f M aln/s)   %.2fx   (H2D %.3fs)\n",
         gpu_h2d+gpu_warp, N/(gpu_h2d+gpu_warp)/1e6, cpu_s/(gpu_h2d+gpu_warp), gpu_h2d);
  return 0;
}
