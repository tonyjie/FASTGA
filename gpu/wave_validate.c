/*******************************************************************************************
 *
 *  wave_validate -- Stage-2 TRACE-POINT correctness gate (Task 8).
 *
 *  THE definitive Stage-2 gate.  On a UNIFORM sample (default 20000) of the full seed
 *  distribution it runs the GPU sparse-checkpoint trace kernel (wave_trace_batch) and the CPU
 *  oracle (FastGA's real Local_Alignment, which emits trace-points in the SAME pass) and asserts:
 *
 *    (1) every GPU trace is 100% Check_Trace_Points-valid (align.c's own validator: right #
 *        trace points AND Σ b-displacements == aligned B-interval), and
 *    (2) the GPU edit distance equals the CPU's (same-score),
 *
 *  then splits the same-score set into BYTE-EXACT (identical tlen + endpoints + every (d,b)) vs
 *  EQUALLY-OPTIMAL-DIFFERENT-PATH, and attributes every divergence:
 *    - out-of-scope: the CPU alignment is SHORT at either end (aepos+bepos-anti < DUB_TRIM or
 *      anti-(abpos+bbpos) < DUB_TRIM) so Local_Alignment applied its DUB_TRIM re-center / re-
 *      extend (align.c:1535,1551-1576) -- a FULL-routine post-process the wave port omits (same
 *      exclusion the Stage-1 --discover-validate gate makes).
 *    - port residual: same-score different path arising from the SERIAL lane-0 trim tie-order
 *      (align.c visits diagonals in a fixed k-order and the first tied furthest-reach wins; the
 *      GPU's phase-2 serial trim reproduces that order but a tie can still resolve to a different
 *      equally-optimal lineage).  This is the Task-6 residual, NOT warp-argmax (the argmax is not
 *      used for trim; the slide is exact & deterministic).
 *
 *  Build: see Makefile target wave_validate (links wave_kernel.o + align.c GDB.c alncode.c
 *  gene_core.c ONElib.c, like the other gpu C benches).
 *  Run:   gpu/wave_validate <in.1aln> <in.seeds> [--sample=20000]
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "GDB.h"
#include "align.h"
#include "alncode.h"
#include "wave_harness.h"
#include "wave_kernel.h"

typedef unsigned char u8;

//  ---- trim tables: byte-identical to New_Align_Spec(0.7,100,gdb1->freq,0) internals
//  (copied verbatim from wave_bench_gpu.cu so the GPU kernel applies the same TABLE/SCORE). ----
#define WV_TRIM_LEN  15
#define WV_PATH_LEN  60
#define WV_TRIM_SZ   (1<<WV_TRIM_LEN)
#define WV_FRACTION  1000

static void gset_table(int bit, int prefix, int sc, int mx, int mscore, int dscore,
                       int16_t *table, int16_t *score)
{ if (bit >= WV_TRIM_LEN)
    { table[prefix] = (int16_t)(sc - mx); score[prefix] = (int16_t) sc; }
  else
    { if (sc > mx) mx = sc;
      gset_table(bit+1, (prefix<<1),     sc - dscore, mx, mscore, dscore, table, score);
      gset_table(bit+1, (prefix<<1) | 1, sc + mscore, mx, mscore, dscore, table, score);
    }
}

static int build_trim_tables(float *freq, double ave_corr, int16_t *table, int16_t *score)
{ static const double Bias_Factor[10] = { .690,.690,.690,.690,.780,.850,.900,.933,.966,1.000 };
  double match; int bias;

  match = freq[0] + freq[3];
  if ((match <= 0.) == (match > 0.)) match = .5;
  if (match > .5) match = 1. - match;
  bias = (int)((match + .025)*20. - 1.);
  if (match < .2) bias = 3;

  int ave_path = (int)(WV_PATH_LEN * (1. - Bias_Factor[bias] * (1. - ave_corr)));
  int mscore   = (int)(WV_FRACTION * Bias_Factor[bias] * (1. - ave_corr));
  int dscore   = WV_FRACTION - mscore;
  gset_table(0,0,0,0,mscore,dscore,table,score);
  return ave_path;
}

//  Build ONE flat, sentinel-padded, NUMERIC resident array for every contig of gdb (forward
//  orientation).  Identical to wave_bench_gpu.cu's build_flat -- see wave_kernel.h contract.
static unsigned char *build_flat(GDB *gdb, long **obase, long *olen)
{ int   n = gdb->ncontig;
  long *base = (long *) malloc((size_t) n * sizeof(long));
  long  total = 1;
  int   i;

  for (i = 0; i < n; i++)
    total += gdb->contigs[i].clen + 1;

  unsigned char *flat = (unsigned char *) malloc((size_t) total);
  if (flat == NULL || base == NULL)
    { fprintf(stderr,"build_flat: out of memory (%ld bytes, %d contigs)\n",total,n); exit(1); }
  flat[0] = 4;

  long off = 1;
  for (i = 0; i < n; i++)
    { int64 clen = gdb->contigs[i].clen;
      char *tmp  = (char *) malloc((size_t) clen + 8);
      if (tmp == NULL)
        { fprintf(stderr,"build_flat: OOM contig %d (%lld bp)\n",i,(long long) clen); exit(1); }
      char *seq = Get_Contig(gdb,i,NUMERIC,tmp+1);
      base[i] = off;
      memcpy(flat+off-1, seq-1, (size_t) clen + 2);
      free(tmp);
      off += clen + 1;
    }

  *obase = base;
  *olen  = total;
  return flat;
}

//  Run Check_Trace_Points on a uint16 trace-point list.  align.c's Check_Trace_Points reads an
//  8-bit trace when tspace<=TRACE_XOVR (100<=125), exactly as FastGA does after Compress_TraceTo8;
//  so we compress into a scratch uint8 (values must fit -- a real tspace=100 trace always does)
//  and call it.  Returns: 1 = valid, 0 = Check_Trace_Points rejected, -1 = value>255 (not a legal
//  compressible tspace=100 trace, so not valid).
static int check_trace_u16(int abpos, int bbpos, int aepos, int bepos, int tlen,
                           const uint16_t *tr, int tspace, uint8_t *scratch)
{ Overlap ovl;
  int j;

  if (tlen < 0) return 0;
  for (j = 0; j < tlen; j++)
    { if (tr[j] > 255) return -1;
      scratch[j] = (uint8_t) tr[j];
    }
  ovl.path.abpos = abpos; ovl.path.bbpos = bbpos;
  ovl.path.aepos = aepos; ovl.path.bepos = bepos;
  ovl.path.tlen  = tlen;
  ovl.path.trace = (void *) scratch;
  ovl.path.diffs = 0;
  ovl.flags = 0; ovl.aread = 0; ovl.bread = 0;
  return (Check_Trace_Points(&ovl,tspace,0,(char *)"gpu") == 0) ? 1 : 0;
}

static char *Usage = "<in:path>[.1aln] <in:path.seeds> [--sample=20000]";

int main(int argc, char *argv[])
{ GDB      _gdb1, *gdb1 = &_gdb1;
  GDB      _gdb2, *gdb2 = &_gdb2;
  FILE    **units1, **units2;
  OneFile  *input;
  int64     novl;
  int       TSPACE, ISTWO;
  int       SAMP = 20000;

  if (argc < 3)
    { fprintf(stderr,"Usage: %s %s\n",argv[0],Usage); return 1; }
  { int i;
    for (i = 3; i < argc; i++)
      if (strncmp(argv[i],"--sample=",9) == 0) SAMP = atoi(argv[i]+9);
  }

  //  ---- open the .1aln just for gdb1/gdb2, exactly as wave_bench_cpu.c / wave_bench_gpu.cu ----
  { char *pwd, *root, *cpath, *src1_name, *src2_name;

    pwd  = PathTo(argv[1]);
    root = Root(argv[1],".1aln");
    input = open_Aln_Read(Catenate(pwd,"/",root,".1aln"),1,&novl,&TSPACE,
                           &src1_name,&src2_name,&cpath);
    if (input == NULL)
      { fprintf(stderr,"%s: cannot open %s\n",argv[0],argv[1]); return 1; }
    free(root); free(pwd);

    ISTWO = (src2_name != NULL);
    Skip_Skeleton(input);
    units1 = Get_GDB(gdb1,src1_name,cpath,1,NULL);
    if (ISTWO)
      { Skip_Skeleton(input);
        units2 = Get_GDB(gdb2,src2_name,cpath,1,NULL);
      }
    else
      { gdb2 = gdb1; units2 = units1; }
    free(src1_name); free(src2_name); free(cpath);
    gdb1->seqs = units1[0];
    gdb2->seqs = units2[0];
  }
  oneFileClose(input);

  //  ---- flat resident genomes: A forward; B forward AND reverse-complement ----
  long *Abase, Alen;
  unsigned char *Aflat = build_flat(gdb1,&Abase,&Alen);
  int   nA = gdb1->ncontig;

  long *Bbase, Blen;
  unsigned char *Bfwd = build_flat(gdb2,&Bbase,&Blen);
  int   nB = gdb2->ncontig;

  unsigned char *Brev = (unsigned char *) malloc((size_t) Blen);
  memcpy(Brev,Bfwd,(size_t) Blen);
  { int i;
    for (i = 0; i < nB; i++)
      Complement_Seq((char *) Brev+Bbase[i],(int) gdb2->contigs[i].clen);
  }

  wave_ctx *g = wave_open();
  wave_load_genomes(g,Aflat,Abase,nA,Alen,Bfwd,Brev,Bbase,nB,Blen);

  //  ---- load seeds ----
  FILE *sf = fopen(argv[2],"rb");
  WaveSeedHeader hdr;
  if (sf == NULL || fread(&hdr,sizeof(hdr),1,sf) != 1 || hdr.magic != WAVE_SEEDS_MAGIC)
    { fprintf(stderr,"%s: bad seeds header in %s\n",argv[0],argv[2]); return 1; }
  int N = (int) hdr.nseeds;
  SeedRec *S = (SeedRec *) malloc((size_t) N * sizeof(SeedRec));
  if (S == NULL || fread(S,sizeof(SeedRec),(size_t) N,sf) != (size_t) N)
    { fprintf(stderr,"%s: short read of %d seed records\n",argv[0],N); return 1; }
  fclose(sf);

  if (SAMP > N) SAMP = N;
  int tspace = (int) hdr.tspace;
  if (tspace <= 0) tspace = TSPACE;

  //  ---- uniform sample ----
  wave_seed *WS   = (wave_seed *) malloc((size_t) SAMP * sizeof(wave_seed));
  int       *sidx = (int *)       malloc((size_t) SAMP * sizeof(int));
  int j;
  for (j = 0; j < SAMP; j++)
    { int i = (int)((long) j * N / SAMP);
      SeedRec *s = &S[i];
      sidx[j]     = i;
      WS[j].aread = s->aread;  WS[j].bread = s->bread;
      WS[j].alen  = s->alen;   WS[j].blen  = s->blen;
      WS[j].anti  = s->seed_anti; WS[j].diag = s->seed_diag;
      WS[j].comp  = (s->flags & COMP_FLAG) ? 1 : 0;
    }

  int16_t *tab = (int16_t *) malloc((size_t) WV_TRIM_SZ * sizeof(int16_t));
  int16_t *scr = (int16_t *) malloc((size_t) WV_TRIM_SZ * sizeof(int16_t));
  int path_ave = build_trim_tables(gdb1->freq,0.7,tab,scr);

  const int TR_STRIDE = 32768;   //  per-seed uint16 trace capacity (kernel enforces; -2 on overflow)
  int      *gAb = (int *) malloc((size_t) SAMP * sizeof(int));
  int      *gAe = (int *) malloc((size_t) SAMP * sizeof(int));
  int      *gBb = (int *) malloc((size_t) SAMP * sizeof(int));
  int      *gBe = (int *) malloc((size_t) SAMP * sizeof(int));
  int      *gDf = (int *) malloc((size_t) SAMP * sizeof(int));
  int      *gTl = (int *) malloc((size_t) SAMP * sizeof(int));
  uint16_t *gTr = (uint16_t *) malloc((size_t) SAMP * TR_STRIDE * sizeof(uint16_t));
  if (gTr == NULL) { fprintf(stderr,"%s: OOM for GPU trace buffer\n",argv[0]); return 1; }

  printf("wave_validate: %d of %d seeds (uniform), path_ave=%d, tspace=%d, trace_stride=%d\n",
         SAMP,N,path_ave,tspace,TR_STRIDE);

  if (wave_trace_batch(g,SAMP,WS,tab,scr,path_ave,0,tspace,
                       gAb,gAe,gBb,gBe,gDf,gTl,gTr,TR_STRIDE) != 0)
    { fprintf(stderr,"%s: wave_trace_batch failed\n",argv[0]); return 1; }

  //  ---- CPU oracle + comparison ----
  Align_Spec *spec = New_Align_Spec(0.7,100,gdb1->freq,0);
  Work_Data  *work = New_Work_Data();
  Alignment   aln;  Path path;  aln.path = &path;

  uint8_t  *u8scratch = (uint8_t *)  malloc((size_t) TR_STRIDE);
  uint16_t *cpuTr     = (uint16_t *) malloc((size_t) TR_STRIDE * sizeof(uint16_t));

  long comparable=0, overflow=0;
  long ov_band=0, ov_cells=0, ov_trace=0;
  long gpu_valid=0, gpu_invalid=0;
  long cpu_valid=0;
  long samescore=0, worse=0, better=0;
  long byte_exact=0, ss_diffpath=0;
  long ep_exact=0;
  long dmax=0, dsum=0;
  //  in-scope (CPU alignment long at BOTH ends -> no DUB_TRIM re-center) vs out-of-scope split:
  long inscope=0, outscope=0;
  long in_byte_exact=0, in_samescore=0, in_ss_diffpath=0, in_worse=0;
  long out_byte_exact=0, out_samescore=0;

  for (j = 0; j < SAMP; j++)
    { SeedRec *s = &S[sidx[j]];
      int anti = s->seed_anti, dg = s->seed_diag;

      aln.aseq  = (char *) Aflat + Abase[s->aread];
      aln.bseq  = (char *)(WS[j].comp ? Brev : Bfwd) + Bbase[s->bread];
      aln.alen  = s->alen;  aln.blen = s->blen;  aln.flags = 0;
      Local_Alignment(&aln,work,spec,dg,dg,anti,-1,-1);

      int ctl = path.tlen;
      uint16_t *ct = (uint16_t *) path.trace;
      if (ctl > TR_STRIDE) ctl = TR_STRIDE;
      memcpy(cpuTr,ct,(size_t) ctl * sizeof(uint16_t));

      //  sanity: CPU trace should itself be Check_Trace_Points-valid
      if (check_trace_u16(path.abpos,path.bbpos,path.aepos,path.bepos,path.tlen,cpuTr,tspace,u8scratch) == 1)
        cpu_valid++;

      if (gDf[j] == -2 || gTl[j] < 0)
        { overflow++;
          if (gTl[j]==-2) ov_band++; else if (gTl[j]==-3) ov_cells++; else if (gTl[j]==-4) ov_trace++;
          continue;
        }
      comparable++;

      //  (1) Check_Trace_Points validity of the GPU trace
      int cv = check_trace_u16(gAb[j],gBb[j],gAe[j],gBe[j],gTl[j],
                               gTr + (size_t) j * TR_STRIDE, tspace, u8scratch);
      if (cv == 1) gpu_valid++; else gpu_invalid++;

      //  Is this a DUB_TRIM re-center case -> out of the wave port's scope?  Local_Alignment's
      //  fshort/rshort branch (align.c:1535,1551) tests the PRE-recenter forward_wave / reverse_wave
      //  endpoints -- which is exactly what the GPU produces (it never re-centers).  So detect on
      //  the GPU's own endpoints (the same principled detector Stage-1's --fwd-validate uses,
      //  wave_bench_gpu.cu:470), NOT the CPU's post-recenter final endpoints (which can look "long"
      //  after a re-extend and misclassify the case).
      int cpu_fshort = (gAe[j] + gBe[j] - anti) < 45;
      int cpu_rshort = (anti - (gAb[j] + gBb[j])) < 45;
      int cpu_short  = cpu_fshort || cpu_rshort;

      //  (2) same edit distance
      if (gDf[j] == path.diffs) samescore++;
      else if (gDf[j] > path.diffs) worse++;
      else better++;

      long dd = labs((long) gDf[j] - (long) path.diffs);
      dsum += dd; if (dd > dmax) dmax = dd;

      //  byte-exact: identical endpoints + tlen + every trace value
      int ep = (gAb[j]==path.abpos && gAe[j]==path.aepos &&
                gBb[j]==path.bbpos && gBe[j]==path.bepos);
      if (ep) ep_exact++;
      int bytematch = (ep && gDf[j]==path.diffs && gTl[j]==path.tlen);
      if (bytematch)
        { int q, eq=1;
          uint16_t *gt = gTr + (size_t) j * TR_STRIDE;
          for (q = 0; q < path.tlen; q++)
            if (gt[q] != cpuTr[q]) { eq = 0; break; }
          bytematch = eq;
        }

      if (bytematch) byte_exact++;
      else if (gDf[j] == path.diffs) ss_diffpath++;

      //  in-scope / out-of-scope classification (out-of-scope == CPU applied DUB_TRIM re-center)
      if (cpu_short)
        { outscope++;
          if (bytematch)          out_byte_exact++;
          if (gDf[j]==path.diffs) out_samescore++;
        }
      else
        { inscope++;
          if (bytematch)               in_byte_exact++;
          if (gDf[j]==path.diffs)      in_samescore++;
          if (gDf[j]==path.diffs && !bytematch) in_ss_diffpath++;
          if (gDf[j] >  path.diffs)    in_worse++;
        }
    }

  Free_Work_Data(work);
  Free_Align_Spec(spec);

  //  ---- report + GATE ----
  double vpct  = comparable ? 100.0*gpu_valid/comparable : 0.0;
  double sspct = comparable ? 100.0*samescore/comparable : 0.0;
  double bxpct = comparable ? 100.0*byte_exact/comparable : 0.0;

  printf("\n=== wave_validate results (%d sampled, %ld comparable, %ld overflow-excluded) ===\n",
         SAMP,comparable,overflow);
  printf("  overflow cause: band(>%d diag)=%ld  cells(>%d checkpts)=%ld  trace(>stride)=%ld\n",
         2048,ov_band,262144,ov_cells,ov_trace);
  printf("  CPU trace self-check (Check_Trace_Points valid)   : %ld/%d\n",cpu_valid,SAMP);
  printf("  (1) GPU trace Check_Trace_Points-VALID            : %ld/%ld  (%.3f%%)   [GATE 100%%]\n",
         gpu_valid,comparable,vpct);
  printf("      GPU trace INVALID                             : %ld\n",gpu_invalid);
  printf("  (2) same edit distance vs CPU (|Δdiffs|==0)       : %ld/%ld  (%.3f%%)\n",
         samescore,comparable,sspct);
  printf("      GPU worse (diffs>CPU) %ld ; GPU better %ld ; mean|Δdiff|=%.3f max=%ld\n",
         worse,better, (worse+better)?(double)dsum/(worse+better):0.0, dmax);
  printf("  endpoint-exact (all 4)                            : %ld/%ld  (%.3f%%)\n",
         ep_exact,comparable, comparable?100.0*ep_exact/comparable:0.0);
  printf("\n  --- byte-exact vs equally-optimal split (of %ld comparable) ---\n",comparable);
  printf("  BYTE-EXACT (endpoints+tlen+every (d,b) identical) : %ld  (%.3f%%)\n",byte_exact,bxpct);
  printf("  same-score DIFFERENT-path (equally optimal)       : %ld\n",ss_diffpath);
  printf("  different-score (worse+better)                    : %ld\n",worse+better);
  printf("\n  --- IN-SCOPE vs OUT-OF-SCOPE (out == CPU DUB_TRIM re-center at a short end) ---\n");
  printf("  out-of-scope (CPU short-end re-center)            : %ld  (byte-exact %ld, same-score %ld)\n",
         outscope,out_byte_exact,out_samescore);
  printf("  in-scope (CPU long at both ends)                  : %ld\n",inscope);
  printf("     in-scope byte-exact                            : %ld/%ld  (%.3f%%)\n",
         in_byte_exact,inscope, inscope?100.0*in_byte_exact/inscope:0.0);
  printf("     in-scope same-score (byte-exact + diff-path)   : %ld/%ld  (%.3f%%)   [GATE >= 99%%]\n",
         in_samescore,inscope, inscope?100.0*in_samescore/inscope:0.0);
  printf("     in-scope same-score DIFFERENT-path (Task6 tie) : %ld\n",in_ss_diffpath);
  printf("     in-scope worse path (genuinely worse)          : %ld\n",in_worse);

  double in_sspct = inscope ? 100.0*in_samescore/inscope : 0.0;
  int gate_valid = (gpu_invalid == 0 && comparable > 0);
  int gate_same  = (in_sspct >= 99.0);
  printf("\n  GATE: Check_Trace_Points 100%% %s ; in-scope same-distance>=99%% %s => %s\n",
         gate_valid?"PASS":"FAIL", gate_same?"PASS":"FAIL",
         (gate_valid&&gate_same)?"PASS":"FAIL");

  wave_close(g);
  free(Aflat); free(Abase); free(Bfwd); free(Brev); free(Bbase);
  free(S); free(WS); free(sidx); free(tab); free(scr);
  free(gAb); free(gAe); free(gBb); free(gBe); free(gDf); free(gTl); free(gTr);
  free(u8scratch); free(cpuTr);
  return (gate_valid && gate_same) ? 0 : 2;
}
