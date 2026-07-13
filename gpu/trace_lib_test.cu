/* A2 self-test: validate the BATCHED gpu_trace_batch library path against FastGA's
 * reference trace-points (gpu/trace_format.h from extract_trace).
 *
 * Each .trace task is an independent rectangle; we pack them all into ONE resident A/B
 * buffer and call gpu_trace_batch once.  Panel binning depends only on abpos mod tspace
 * (panel(al)=floor((abpos%ts+al)/ts)), so we pad each task's A-offset to be congruent to
 * its original abpos mod tspace -> identical panels to the reference.  B offset is a pure
 * buffer index (no phasing role).
 *
 * Build: make trace_lib_test    Run: ./gpu/trace_lib_test <tasks.trace>
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "trace_format.h"
#include "fastga_gpu.h"

typedef unsigned char u8;

int main(int argc,char**argv){
  if(argc<2){fprintf(stderr,"usage: %s <tasks.trace>\n",argv[0]);return 1;}
  FILE*f=fopen(argv[1],"rb");
  if(!f){fprintf(stderr,"cannot open %s\n",argv[1]);return 1;}
  FGATraceHeader h;
  if(fread(&h,sizeof(h),1,f)!=1||h.magic!=FGA_TRACE_MAGIC){fprintf(stderr,"bad header\n");return 1;}
  int ts=(int)h.tspace;
  printf("tasks: %u  tspace: %d\n",h.ntasks,ts);

  std::vector<u8> Abuf,Bbuf;
  std::vector<int> AB,AE,BB,BE,refTl,refOff;
  std::vector<unsigned short> refTr;
  for(uint32_t i=0;i<h.ntasks;i++){
    int32_t ab,ae,bb,be; uint32_t aw,bw,rd,rl;
    if(fread(&ab,4,1,f)!=1)break;
    fread(&ae,4,1,f);fread(&bb,4,1,f);fread(&be,4,1,f);
    fread(&aw,4,1,f);fread(&bw,4,1,f);fread(&rd,4,1,f);fread(&rl,4,1,f);
    std::vector<u8> A(aw),B(bw); if(aw)fread(A.data(),1,aw,f); if(bw)fread(B.data(),1,bw,f);
    std::vector<unsigned short> T(rl); if(rl)fread(T.data(),2,rl,f);
    if(aw==0||bw==0) continue;
    int r = ((ab % ts)+ts)%ts;
    while((int)(Abuf.size()%ts)!=r) Abuf.push_back(0);   // phase-align A offset
    int aoff=(int)Abuf.size(); Abuf.insert(Abuf.end(),A.begin(),A.end());
    int boff=(int)Bbuf.size(); Bbuf.insert(Bbuf.end(),B.begin(),B.end());
    AB.push_back(aoff); AE.push_back(aoff+(int)aw);
    BB.push_back(boff); BE.push_back(boff+(int)bw);
    refTl.push_back((int)rl); refOff.push_back((int)refTr.size());
    refTr.insert(refTr.end(),T.begin(),T.end());
  }
  fclose(f);
  int N=(int)AB.size();
  printf("packed %d tasks; resident A=%zu B=%zu bytes\n",N,Abuf.size(),Bbuf.size());

  gpu_ctx*g=gpu_open();
  gpu_load_seqs(g,Abuf.data(),(int)Abuf.size(),Bbuf.data(),(int)Bbuf.size());
  std::vector<unsigned short> out((size_t)N*FGA_TRACE_MAX_PAIRS);
  std::vector<int> tlen(N);
  int rc=gpu_trace_batch(g,N,AB.data(),AE.data(),BB.data(),BE.data(),ts,out.data(),tlen.data());
  gpu_close(g);
  if(rc){fprintf(stderr,"gpu_trace_batch failed rc=%d\n",rc);return 1;}

  long valid=0,exact=0,tlok=0,ov=0;
  for(int i=0;i<N;i++){
    if(tlen[i]<0){ov++;continue;}
    valid++;
    int rl=refTl[i];
    if(tlen[i]==rl){ tlok++;
      const unsigned short*gp=&out[(size_t)i*FGA_TRACE_MAX_PAIRS];
      const unsigned short*rp=&refTr[refOff[i]];
      int eq=1; for(int p=0;p<rl;p++) if(gp[p]!=rp[p]){eq=0;break;}
      if(eq)exact++;
    }
  }
  printf("valid: %ld/%d  [overflow=%ld]\n",valid,N,ov);
  printf("  tlen match:   %ld (%.2f%% of valid)\n",tlok,100.0*tlok/valid);
  printf("  EXACT match:  %ld (%.2f%% of valid)\n",exact,100.0*exact/valid);
  return 0;
}
