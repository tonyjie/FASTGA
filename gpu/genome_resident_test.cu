/* Validate the A3b genome-resident kernels: pack all .disc windows into ONE resident "genome"
 * with per-task contig bounds, run gpu_discover_batch_g once, and check the (genome-relative)
 * endpoints match the per-task gpu_discover_batch (already validated) exactly. Then the same
 * for gpu_trace_batch_g vs gpu_trace_batch on a .trace file.
 * Build: make genome_resident_test   Run: ./gpu/genome_resident_test <disc> <trace>
 */
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "disc_format.h"
#include "trace_format.h"
#include "fastga_gpu.h"
typedef unsigned char u8;

int main(int argc,char**argv){
  if(argc<3){fprintf(stderr,"usage: %s <disc> <trace>\n",argv[0]);return 1;}
  gpu_ctx*g=gpu_open();

  // ---- discovery: per-task vs genome-resident ----
  { FILE*f=fopen(argv[1],"rb"); FGADiscHeader h;
    if(fread(&h,sizeof(h),1,f)!=1||h.magic!=FGA_DISC_MAGIC){fprintf(stderr,"bad disc\n");return 1;}
    std::vector<u8> GA,GB; std::vector<int> SA,SB,aLo,aHi,bLo,bHi,Aoff,Boff;
    std::vector<std::vector<u8>> As,Bs; std::vector<int> sa0,sb0;
    int N=(int)h.ntasks; if(N>20000)N=20000;
    for(int i=0;i<N;i++){
      uint32_t aw,bw,sa,sb,rd; int rab,rae,rbb,rbe;
      if(fread(&aw,4,1,f)!=1)break; fread(&bw,4,1,f);fread(&sa,4,1,f);fread(&sb,4,1,f);
      fread(&rab,4,1,f);fread(&rae,4,1,f);fread(&rbb,4,1,f);fread(&rbe,4,1,f);fread(&rd,4,1,f);
      std::vector<u8> A(aw),B(bw); if(aw)fread(A.data(),1,aw,f); if(bw)fread(B.data(),1,bw,f);
      if(aw==0||bw==0)continue;
      int ao=(int)GA.size(), bo=(int)GB.size();
      Aoff.push_back(ao);Boff.push_back(bo);
      SA.push_back(ao+(int)sa);SB.push_back(bo+(int)sb);
      aLo.push_back(ao);aHi.push_back(ao+(int)aw);bLo.push_back(bo);bHi.push_back(bo+(int)bw);
      GA.insert(GA.end(),A.begin(),A.end());GB.insert(GB.end(),B.begin(),B.end());
      As.push_back(std::move(A));Bs.push_back(std::move(B));sa0.push_back(sa);sb0.push_back(sb);
    }
    fclose(f);
    int M=(int)SA.size();
    // genome-resident, one call
    std::vector<int> ab(M),ae(M),bb(M),be(M),df(M);
    gpu_load_seqs(g,GA.data(),(int)GA.size(),GB.data(),(int)GB.size());
    gpu_discover_batch_g(g,M,SA.data(),SB.data(),aLo.data(),aHi.data(),bLo.data(),bHi.data(),
                         ab.data(),ae.data(),bb.data(),be.data(),df.data());
    // per-task reference
    long mism=0;
    for(int i=0;i<M;i++){
      int pab,pae,pbb,pbe,pdf,s=sa0[i],t=sb0[i];
      gpu_load_seqs(g,As[i].data(),(int)As[i].size(),Bs[i].data(),(int)Bs[i].size());
      gpu_discover_batch(g,1,&s,&t,&pab,&pae,&pbb,&pbe,&pdf);
      if(ab[i]-Aoff[i]!=pab||ae[i]-Aoff[i]!=pae||bb[i]-Boff[i]!=pbb||be[i]-Boff[i]!=pbe||df[i]!=pdf) mism++;
    }
    printf("discovery genome-resident vs per-task: %d tasks, %ld mismatches -> %s\n",
           M,mism,mism?"FAIL":"PASS");
  }

  // ---- trace: per-task vs genome-resident ----
  { FILE*f=fopen(argv[2],"rb"); FGATraceHeader h;
    if(fread(&h,sizeof(h),1,f)!=1||h.magic!=FGA_TRACE_MAGIC){fprintf(stderr,"bad trace\n");return 1;}
    int ts=(int)h.tspace;
    std::vector<u8> GA,GB; std::vector<int> AB,AE,BB,BE,BASE,Aoff,Boff;
    std::vector<std::vector<u8>> As,Bs; std::vector<int> ab0_mod;
    int N=(int)h.ntasks; if(N>20000)N=20000;
    for(int i=0;i<N;i++){
      int32_t ab,ae,bb,be; uint32_t aw,bw,rd,rl;
      if(fread(&ab,4,1,f)!=1)break; fread(&ae,4,1,f);fread(&bb,4,1,f);fread(&be,4,1,f);
      fread(&aw,4,1,f);fread(&bw,4,1,f);fread(&rd,4,1,f);fread(&rl,4,1,f);
      std::vector<u8> A(aw),B(bw); if(aw)fread(A.data(),1,aw,f); if(bw)fread(B.data(),1,bw,f);
      if(rl)fseek(f,(long)rl*2,SEEK_CUR);
      if(aw==0||bw==0)continue;
      // pack; genome A-offset carries an arbitrary base; phasing uses base = ao - (ab mod ts)
      int r=((ab%ts)+ts)%ts;
      while((int)(GA.size()%ts)!=r) GA.push_back(0);   // keep ao ≡ ab (mod ts) so base is clean
      int ao=(int)GA.size(), bo=(int)GB.size();
      AB.push_back(ao);AE.push_back(ao+(int)aw);BB.push_back(bo);BE.push_back(bo+(int)bw);
      BASE.push_back(ao-r);                            // contig base so (ao-base)=r=ab mod ts
      Aoff.push_back(ao);Boff.push_back(bo);
      GA.insert(GA.end(),A.begin(),A.end());GB.insert(GB.end(),B.begin(),B.end());
      As.push_back(std::move(A));Bs.push_back(std::move(B));ab0_mod.push_back(r);
    }
    fclose(f);
    int M=(int)AB.size();
    std::vector<unsigned short> outG((size_t)M*FGA_TRACE_MAX_PAIRS); std::vector<int> tlG(M);
    gpu_load_seqs(g,GA.data(),(int)GA.size(),GB.data(),(int)GB.size());
    gpu_trace_batch_g(g,M,AB.data(),AE.data(),BB.data(),BE.data(),BASE.data(),ts,outG.data(),tlG.data());
    long mism=0;
    std::vector<unsigned short> outP(FGA_TRACE_MAX_PAIRS); int tlP;
    for(int i=0;i<M;i++){
      int r=ab0_mod[i], aw=(int)As[i].size();
      // per-task must phase on r AND read the window: place the window at offset r (r-pad).
      std::vector<u8> pbuf(r+aw,0); for(int k=0;k<aw;k++) pbuf[r+k]=As[i][k];
      int pab=r, pae=r+aw, pbb=0, pbe=(int)Bs[i].size();
      gpu_load_seqs(g,pbuf.data(),(int)pbuf.size(),Bs[i].data(),(int)Bs[i].size());
      gpu_trace_batch(g,1,&pab,&pae,&pbb,&pbe,ts,outP.data(),&tlP);
      if(tlG[i]!=tlP){mism++;continue;}
      for(int p=0;p<tlP;p++) if(outG[(size_t)i*FGA_TRACE_MAX_PAIRS+p]!=outP[p]){mism++;break;}
    }
    printf("trace genome-resident vs per-task: %d tasks, %ld mismatches -> %s\n",
           M,mism,mism?"FAIL":"PASS");
  }
  gpu_close(g);
  return 0;
}
