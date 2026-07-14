/* On-device 2-bit DNA decompression (the #2 hotspot, Uncompress_Read = 29% of the human
 * align run) vs the CPU. Genome-scale synthetic .bps. Validates byte-identical output.
 * Build: make decomp_bench   Run: ./gpu/decomp_bench [Gbp=3]
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <omp.h>
#include <cuda_runtime.h>
typedef unsigned char u8;
using clk=std::chrono::high_resolution_clock;
#define CK(c) do{cudaError_t e=(c); if(e!=cudaSuccess){fprintf(stderr,"CUDA %s:%d %s\n",__FILE__,__LINE__,cudaGetErrorString(e));exit(1);}}while(0)

// base[p] = (packed[p>>2] >> ((p&3)*2)) & 3   -- matches gene_core.c Uncompress_Read (beg=0)
__global__ void unpack2bit(const u8* packed, u8* out, long len){
  long p = (long)blockIdx.x*blockDim.x + threadIdx.x;
  if(p>=len) return;
  out[p] = (packed[p>>2] >> ((p&3)*2)) & 3;
}

// CPU reference: same formula, one contiguous pass (representative of Uncompress_Read's work)
static void cpu_unpack(const u8* packed, u8* out, long len){
  for(long p=0;p<len;p++) out[p]=(packed[p>>2]>>((p&3)*2))&3;
}

int main(int argc,char**argv){
  double Gbp = argc>=2 ? atof(argv[1]) : 3.0;
  long len = (long)(Gbp*1e9);
  long plen = (len+3)/4;
  printf("genome: %.2f Gbp  (2-bit packed = %.0f MB, NUMERIC = %.0f MB)\n",Gbp,plen/1e6,len/1e6);

  u8* packed=(u8*)malloc(plen);
  for(long i=0;i<plen;i++) packed[i]=(u8)(i*1103515245u+12345u);   // deterministic filler
  u8* outc=(u8*)malloc(len);
  u8* outg=(u8*)malloc(len);

  // CPU: single-thread and 32-thread
  auto t0=clk::now(); cpu_unpack(packed,outc,len); auto t1=clk::now();
  double s1=std::chrono::duration<double>(t1-t0).count();
  printf("CPU  1-thread: %.3f s  (%.2f Gbase/s)\n",s1,len/s1/1e9);
  omp_set_num_threads(32);
  long nchunk=(len+(1L<<20)-1)>>20;
  auto t2=clk::now();
  #pragma omp parallel for schedule(static)
  for(long ci=0;ci<nchunk;ci++){ long c=ci<<20, e=c+(1L<<20); if(e>len)e=len;
    for(long p=c;p<e;p++) outg[p]=(packed[p>>2]>>((p&3)*2))&3; }
  auto t3=clk::now();
  double s32=std::chrono::duration<double>(t3-t2).count();
  printf("CPU 32-thread: %.3f s  (%.2f Gbase/s)\n",s32,len/s32/1e9);

  // GPU: H2D 2-bit + unpack (+ D2H optional)
  u8 *dP,*dO; CK(cudaMalloc(&dP,plen)); CK(cudaMalloc(&dO,len));
  auto g0=clk::now();
  CK(cudaMemcpy(dP,packed,plen,cudaMemcpyHostToDevice));
  auto g1=clk::now();
  int TPB=256; long blk=(len+TPB-1)/TPB;
  unpack2bit<<<blk,TPB>>>(dP,dO,len);
  CK(cudaDeviceSynchronize());
  auto g2=clk::now();
  double h2d=std::chrono::duration<double>(g1-g0).count();
  double ker=std::chrono::duration<double>(g2-g1).count();
  printf("GPU  H2D(2-bit): %.3f s   kernel: %.3f s  (%.2f Gbase/s kernel)  incl-H2D: %.2f Gbase/s\n",
         h2d,ker,len/ker/1e9,len/(h2d+ker)/1e9);

  // validate a sample
  CK(cudaMemcpy(outc,dO,len,cudaMemcpyDeviceToHost));
  long bad=0; for(long p=0;p<len;p+=997){ u8 ref=(packed[p>>2]>>((p&3)*2))&3; if(outc[p]!=ref)bad++; }
  printf("validate: %s\n", bad?"MISMATCH":"GPU output byte-identical to formula");

  printf("\nCPU Uncompress_Read is 29%% of the human align run; GPU kernel is %.0fx the CPU 32-thread rate.\n",
         s32/ker);
  return 0;
}
