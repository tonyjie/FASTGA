/*******************************************************************************************
 *
 *  wave_bench_gpu -- GPU driver skeleton for the wave-parallelism characterization study
 *                    (Task 4: genome residency + seed upload; NO wave kernel yet -- Tasks
 *                    5/6 add the actual forward/reverse sweep). Loads both GDBs, builds the
 *                    flat sentinel-padded NUMERIC genome arrays (A forward; B BOTH forward
 *                    (Bfwd) and reverse-complement (Brev), exactly as gpu/wave_bench_cpu.c's
 *                    CPU baseline does), uploads them resident via wave_kernel.cu's
 *                    wave_open/wave_load_genomes, and (if requested) runs a --selftest that
 *                    reads a handful of contig byte-ranges back off the device and checks
 *                    them against the host NUMERIC buffers -- including that the resident
 *                    Brev really is the complement-reverse of the resident Bfwd.
 *
 *  Also reads the .seeds file (Task 2) and builds genome-absolute seed-anchor arrays
 *  (contig-relative seed_anti/seed_diag -> absolute index via the contig base-offset
 *  tables), uploads them to the device, and sanity-checks every seed lands inside its own
 *  contig's span -- the "device seed arrays" Tasks 5/6 will consume.
 *
 *  Usage: wave_bench_gpu <in.1aln> <in.seeds> [--selftest]
 *
 *  Build: see Makefile target wave_bench_gpu (nvcc -O3 -arch=sm_80 -Igpu, links wave_kernel.o
 *  + the FastGA C objects GDB.c/align.c/alncode.c/gene_core.c/ONElib.c).
 *
 *******************************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <ctime>
#include <vector>
#include <cuda_runtime.h>

// GDB.h/align.h/alncode.h/gene_core.h have no extern "C" guards of their own (they are
// plain C headers used elsewhere only by gcc-compiled .c files); wrap them here so nvcc's
// C++ host pass emits C-linkage symbol names that match the gcc-compiled .o's we link
// against (GDB.o/align.o/alncode.o/gene_core.o/ONElib.o), per the Makefile's link line.
extern "C" {
#include "GDB.h"
#include "align.h"
#include "alncode.h"
}

#include "wave_harness.h"
#include "wave_kernel.h"

typedef unsigned char u8;

static double now(void)
{ struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC,&ts);
  return ts.tv_sec + ts.tv_nsec*1e-9;
}

//  ---- Replicate align.c's suffix-positivity trim tables + ave_path (New_Align_Spec/set_table,
//  align.c:197-269) so the GPU kernel can apply the identical TABLE/SCORE test.  We can't read
//  them off the opaque _Align_Spec, so we rebuild them from the same inputs (freq, ave_corr).
#define WV_TRIM_LEN  15
#define WV_PATH_LEN  60
#define WV_TRIM_SZ   (1<<WV_TRIM_LEN)     // 32768
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

//  Build the int16 table[32768] + score[32768] and return ave_path -- byte-identical to
//  New_Align_Spec(ave_corr,100,freq,0)'s internal spec->table/score/ave_path.
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

/*  Build ONE flat, sentinel-padded, NUMERIC resident array for every contig of `gdb`
 *  (forward orientation only -- callers wanting the reverse-complement copy it and
 *  Complement_Seq() each contig's span in place, see main()).
 *
 *  Layout (see wave_kernel.h for the full contract): contig i occupies
 *  flat[base[i]-1 .. base[i]+clen_i] inclusive (clen_i+2 bytes: leading sentinel, clen_i
 *  NUMERIC bases, trailing sentinel), contigs packed back to back so each contig's
 *  trailing sentinel doubles as the next contig's leading sentinel. Built via a private
 *  clen+8-padded scratch buffer per contig (Get_Contig's Uncompress_Read can write a few
 *  bytes past [len], per GDB.c's own New_Contig_Buffer margin -- see wave_bench_cpu.c's
 *  identical comment) so no contig's decompression scratch can corrupt its neighbor; only
 *  the confirmed-valid clen+2 byte range [-1..len] is copied into the shared flat array.
 */
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
  flat[0] = 4;   // global leading sentinel (overwritten redundantly by contig 0's own copy)

  long off = 1;
  for (i = 0; i < n; i++)
    { int64 clen = gdb->contigs[i].clen;
      char *tmp  = (char *) malloc((size_t) clen + 8);   // +8 margin, see header comment
      if (tmp == NULL)
        { fprintf(stderr,"build_flat: OOM loading contig %d (%lld bp)\n",i,(long long) clen); exit(1); }
      char *seq = Get_Contig(gdb,i,NUMERIC,tmp+1);        // seq[-1]==4, seq[clen]==4
      base[i] = off;
      memcpy(flat+off-1, seq-1, (size_t) clen + 2);
      free(tmp);
      off += clen + 1;
    }

  *obase = base;
  *olen  = total;
  return flat;
}

/*  Read back [base-1, base+clen+1) (the sentinel-inclusive span of one contig) from a
 *  device-resident array via `rb` and compare byte-for-byte against the host flat array.
 *  Returns 0 on match, 1 on mismatch/out-of-bounds. */
typedef int (*readback_fn)(wave_ctx *, long, long, unsigned char *);

static int check_contig(const char *label, readback_fn rb, wave_ctx *g,
                        const unsigned char *hostFlat, long base, long clen)
{ long lo = base - 1, n = clen + 2;
  unsigned char *buf = (unsigned char *) malloc((size_t) n);
  int bad = 0;

  if (rb(g,lo,n,buf) != 0)
    { fprintf(stderr,"%s: readback [%ld,%ld) out of bounds\n",label,lo,lo+n); free(buf); return 1; }
  for (long k = 0; k < n; k++)
    if (buf[k] != hostFlat[lo+k])
      bad++;
  free(buf);
  if (bad)
    fprintf(stderr,"%s: %d/%ld bytes mismatch in [%ld,%ld)\n",label,bad,n,lo,lo+n);
  return bad != 0;
}

/*  Verify the DEVICE-resident Brev really is the complement-reverse of the DEVICE-resident
 *  Bfwd for one contig (both sides read back from the device, not from host memory) --
 *  the literal "Brev = complement of Bfwd" acceptance check. */
static int check_brev_is_complement(wave_ctx *g, long base, long clen)
{ unsigned char *fwd = (unsigned char *) malloc((size_t) clen);
  unsigned char *rev = (unsigned char *) malloc((size_t) clen);
  int bad = 0;

  if (wave_readback_Bfwd(g,base,clen,fwd) != 0 || wave_readback_Brev(g,base,clen,rev) != 0)
    { fprintf(stderr,"check_brev_is_complement: readback out of bounds\n"); free(fwd); free(rev); return 1; }
  for (long k = 0; k < clen; k++)
    if (rev[k] != (unsigned char) (3 - fwd[clen-1-k]))
      bad++;
  free(fwd);
  free(rev);
  if (bad)
    fprintf(stderr,"check_brev_is_complement: %d/%ld positions violate Brev[k]==3-Bfwd[clen-1-k]\n",bad,clen);
  return bad != 0;
}

//  Pick up to 4 sample contig indices (first, ~1/3, ~2/3, last), deduplicated.
static int pick_samples(int n, int *out)
{ int cand[4] = { 0, n/3, (2*n)/3, n-1 };
  int m = 0, i, j;
  for (i = 0; i < 4; i++)
    { if (cand[i] < 0 || cand[i] >= n) continue;
      for (j = 0; j < m; j++) if (out[j] == cand[i]) break;
      if (j == m) out[m++] = cand[i];
    }
  return m;
}

//  ---- Task 7: (band,depth) stratification buckets -- IDENTICAL to wave_bench_cpu.c's, so
//  GPU and CPU per-bucket aln/s are directly comparable (the fit map). ----
#define G_NBAND 5
#define G_NDEPTH 5
static long G_BAND_HI[G_NBAND]   = {   8,   32,  128,   512, INT32_MAX };
static long G_DEPTH_HI[G_NDEPTH] = {  10,   50,  200,  1000, INT32_MAX };

static int gbucket_of(long v, long *hi, int n)
{ int i;
  for (i = 0; i < n; i++)
    if (v <= hi[i])
      return i;
  return n-1;
}

//  derive the .cpuref path from a .seeds path: strip ".seeds", append ".cpuref"
static void cpuref_path_from_seeds(const char *seeds_path, char *out, size_t outsz)
{ size_t L = strlen(seeds_path);
  const char *suf = ".seeds";  size_t sl = strlen(suf);
  if (L > sl && strcmp(seeds_path+L-sl,suf) == 0)
    { memcpy(out,seeds_path,L-sl); out[L-sl] = 0; }
  else
    { memcpy(out,seeds_path,L); out[L] = 0; }
  strncat(out,".cpuref",outsz-strlen(out)-1);
}

static char *Usage =
  "<in:path>[.1aln] <in:path.seeds> [--selftest] [--fwd-validate] [--discover-validate]\n"
  "       [--stage1-fit] [--cpu-1core=<aln/s>] [--cpu-32core=<aln/s>]";

int main(int argc, char *argv[])
{ GDB      _gdb1, *gdb1 = &_gdb1;
  GDB      _gdb2, *gdb2 = &_gdb2;
  FILE    **units1, **units2;
  OneFile  *input;
  int64     novl;
  int       TSPACE, ISTWO;
  int       selftest = 0;
  int       fwdvalidate = 0;
  int       discvalidate = 0;
  int       stage1fit = 0;
  //  defaults = Task 3's measured, unpinned numbers; override with --cpu-1core=/--cpu-32core=
  //  after a fresh pinned re-run (Task 7 brief: "re-confirm").
  double    cpu_1core_alnps  = 1212.0;
  double    cpu_32core_alnps = 37782.9;

  if (argc < 3)
    { fprintf(stderr,"Usage: %s %s\n",argv[0],Usage); return 1; }
  { int i;
    for (i = 3; i < argc; i++)
      { if (strcmp(argv[i],"--selftest") == 0)          selftest = 1;
        if (strcmp(argv[i],"--fwd-validate") == 0)      fwdvalidate = 1;
        if (strcmp(argv[i],"--discover-validate") == 0) discvalidate = 1;
        if (strcmp(argv[i],"--stage1-fit") == 0)        stage1fit = 1;
        if (strncmp(argv[i],"--cpu-1core=",12) == 0)    cpu_1core_alnps  = atof(argv[i]+12);
        if (strncmp(argv[i],"--cpu-32core=",13) == 0)   cpu_32core_alnps = atof(argv[i]+13);
      }
  }

  //  ---- open the .1aln just for gdb1/gdb2 (paths+skeletons), exactly as wave_bench_cpu.c ----

  { char *pwd, *root, *cpath, *src1_name, *src2_name;

    pwd  = PathTo(argv[1]);
    root = Root(argv[1],".1aln");
    input = open_Aln_Read(Catenate(pwd,"/",root,".1aln"),1,&novl,&TSPACE,
                           &src1_name,&src2_name,&cpath);
    if (input == NULL)
      { fprintf(stderr,"%s: cannot open %s\n",argv[0],argv[1]); return 1; }
    free(root);
    free(pwd);

    ISTWO = (src2_name != NULL);

    Skip_Skeleton(input);
    units1 = Get_GDB(gdb1,src1_name,cpath,1,NULL);

    if (ISTWO)
      { Skip_Skeleton(input);
        units2 = Get_GDB(gdb2,src2_name,cpath,1,NULL);
      }
    else
      { gdb2   = gdb1;
        units2 = units1;
      }

    free(src1_name);
    free(src2_name);
    free(cpath);

    gdb1->seqs = units1[0];
    gdb2->seqs = units2[0];
  }
  oneFileClose(input);

  //  ---- build the flat resident genome arrays: A forward; B forward AND reverse-complement ----

  double tl0 = now();

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

  double tl1 = now();
  printf("host genome build (Get_Contig, A + Bfwd + Brev): %.3f s "
         "(A: %d contigs, %lld bp; B: %d contigs, %lld bp, fwd+rev)\n",
         tl1-tl0, nA,(long long) gdb1->seqtot, nB,(long long) gdb2->seqtot);

  //  ---- upload resident, timing the one-time H2D residency cost separately ----

  wave_ctx *g = wave_open();

  double th0 = now();
  wave_load_genomes(g,Aflat,Abase,nA,Alen,Bfwd,Brev,Bbase,nB,Blen);
  cudaDeviceSynchronize();
  double th1 = now();

  double gb = (Alen + 2.0*Blen) / 1e9;
  printf("genome H2D upload: %.3f s  (%.3f GB resident: A=%.3fGB Bfwd=%.3fGB Brev=%.3fGB) "
         "-> %.2f GB/s\n", th1-th0, gb, Alen/1e9, Blen/1e9, Blen/1e9, gb/(th1-th0));

  //  ---- --selftest: read back a few contigs' byte-ranges and compare vs host NUMERIC ----

  int fail = 0;
  if (selftest)
    { int idxA[4], idxB[4], nsA, nsB, k;

      nsA = pick_samples(nA,idxA);
      nsB = pick_samples(nB,idxB);

      for (k = 0; k < nsA; k++)
        fail |= check_contig("A",wave_readback_A,g,Aflat,Abase[idxA[k]],gdb1->contigs[idxA[k]].clen);
      for (k = 0; k < nsB; k++)
        { fail |= check_contig("Bfwd",wave_readback_Bfwd,g,Bfwd,Bbase[idxB[k]],gdb2->contigs[idxB[k]].clen);
          fail |= check_contig("Brev",wave_readback_Brev,g,Brev,Bbase[idxB[k]],gdb2->contigs[idxB[k]].clen);
          fail |= check_brev_is_complement(g,Bbase[idxB[k]],gdb2->contigs[idxB[k]].clen);
        }

      printf("--selftest: %d A contig(s), %d B contig(s) checked (byte-match + Brev=complement(Bfwd)) -> %s\n",
             nsA,nsB, fail ? "FAIL" : "PASS");
    }

  //  ---- Step 3: load .seeds, build genome-absolute seed-anchor arrays, upload to device ----

  { FILE *sf = fopen(argv[2],"rb");
    WaveSeedHeader hdr;

    if (sf == NULL)
      { fprintf(stderr,"%s: cannot open %s\n",argv[0],argv[2]); wave_close(g); return 1; }
    if (fread(&hdr,sizeof(hdr),1,sf) != 1 || hdr.magic != WAVE_SEEDS_MAGIC)
      { fprintf(stderr,"%s: bad seeds header in %s\n",argv[0],argv[2]); fclose(sf); wave_close(g); return 1; }

    int       N = (int) hdr.nseeds;
    SeedRec  *S = (SeedRec *) malloc((size_t) N * sizeof(SeedRec));
    if (S == NULL || fread(S,sizeof(SeedRec),(size_t) N,sf) != (size_t) N)
      { fprintf(stderr,"%s: short read of %d seed records\n",argv[0],N); fclose(sf); wave_close(g); return 1; }
    fclose(sf);

    long *SA = (long *) malloc((size_t) N * sizeof(long));
    long *SB = (long *) malloc((size_t) N * sizeof(long));
    long  oob = 0;
    int   i;

    for (i = 0; i < N; i++)
      { SeedRec *s  = &S[i];
        long     sa = ((long) s->seed_anti + (long) s->seed_diag) / 2;   // contig-relative
        long     sb = ((long) s->seed_anti - (long) s->seed_diag) / 2;
        long     absA = Abase[s->aread] + sa;
        long     absB = Bbase[s->bread] + sb;

        SA[i] = absA;
        SB[i] = absB;
        if (sa < 0 || sa > s->alen || sb < 0 || sb > s->blen)
          oob++;
      }

    double ts0 = now();
    long *dSA, *dSB;
    cudaMalloc(&dSA,(size_t) N * sizeof(long));
    cudaMalloc(&dSB,(size_t) N * sizeof(long));
    cudaMemcpy(dSA,SA,(size_t) N * sizeof(long),cudaMemcpyHostToDevice);
    cudaMemcpy(dSB,SB,(size_t) N * sizeof(long),cudaMemcpyHostToDevice);
    double ts1 = now();

    printf("seeds: %d loaded from %s, %ld with contig-relative anchor out of [0,len] bounds, "
           "device seed-array upload %.4f s (%.2f MB)\n",
           N,argv[2],oob,ts1-ts0,(2.0*N*sizeof(long))/1e6);

    cudaFree(dSA);
    cudaFree(dSB);
    free(SA);
    free(SB);

    //  ---- Task 5 --fwd-validate: GPU forward endpoint vs align.c's forward endpoint ----
    //  Oracle = apath->aepos/bepos from the SAME Local_Alignment the CPU baseline runs
    //  (forward_wave is an independent sweep, so its endpoint is identical alone or as part of
    //  the full routine).  Fixed uniform 5000-seed sample; gate = forward-ENDPOINT agreement.

    if (fwdvalidate)
      { int SAMP = 5000; if (SAMP > N) SAMP = N;
        long stride = (long) N / SAMP;
        if (stride < 1) stride = 1;

        wave_seed *WS = (wave_seed *) malloc((size_t) SAMP * sizeof(wave_seed));
        int *sidx = (int *) malloc((size_t) SAMP * sizeof(int));
        int *gAe = (int *) malloc((size_t) SAMP * sizeof(int));
        int *gBe = (int *) malloc((size_t) SAMP * sizeof(int));
        int *gFd = (int *) malloc((size_t) SAMP * sizeof(int));
        int j;
        for (j = 0; j < SAMP; j++)
          { int i = (int)((long) j * N / SAMP);
            SeedRec *s = &S[i];
            sidx[j]      = i;
            WS[j].aread  = s->aread;  WS[j].bread = s->bread;
            WS[j].alen   = s->alen;   WS[j].blen  = s->blen;
            WS[j].anti   = s->seed_anti; WS[j].diag = s->seed_diag;
            WS[j].comp   = (s->flags & 0x1) ? 1 : 0;   // COMP_FLAG
          }

        //  trim tables (identical to New_Align_Spec(0.7,100,gdb1->freq,0) internals)
        int16_t *tab = (int16_t *) malloc((size_t) WV_TRIM_SZ * sizeof(int16_t));
        int16_t *scr = (int16_t *) malloc((size_t) WV_TRIM_SZ * sizeof(int16_t));
        int path_ave = build_trim_tables(gdb1->freq,0.7,tab,scr);
        printf("\n--fwd-validate: sampling %d of %d seeds (stride %ld), path_ave=%d\n",
               SAMP,N,stride,path_ave);

        double g0 = now();
        if (wave_forward_batch(g,SAMP,WS,tab,scr,path_ave,gAe,gBe,gFd) != 0)
          { fprintf(stderr,"%s: wave_forward_batch failed\n",argv[0]); }
        double g1 = now();

        //  host oracle: real Local_Alignment on the sampled seeds (same setup as wave_bench_cpu)
        Align_Spec *spec = New_Align_Spec(0.7,100,gdb1->freq,0);
        Work_Data  *work = New_Work_Data();
        Alignment   aln;  Path path;  aln.path = &path;

        long exact=0, within1=0, overflow=0, fshort_mismatch=0;
        long mism_le5=0, mism_gt5=0, maxdelta=0;
        long longset=0, longset_agree=0, mism_cpu_short=0;   // non-circular cross-checks
        double c0 = now();
        for (j = 0; j < SAMP; j++)
          { SeedRec *s = &S[sidx[j]];
            aln.aseq  = (char *) Aflat + Abase[s->aread];               // A always forward

            aln.bseq  = (char *)(WS[j].comp ? Brev : Bfwd) + Bbase[s->bread];
            aln.alen  = s->alen;  aln.blen = s->blen;  aln.flags = 0;
            int dg = s->seed_diag, anti = s->seed_anti;
            Local_Alignment(&aln,work,spec,dg,dg,anti,-1,-1);

            if (gAe[j] == -2) { overflow++; continue; }
            int dae = gAe[j]-path.aepos, dbe = gBe[j]-path.bepos;
            //  Non-circular subset: seeds whose CPU FINAL endpoint is clearly long
            //  (>=DUB_TRIM past the seed).  A long CPU alignment cannot be a fshort&&rshort
            //  point-collapse, so on it the forward endpoint is either preserved or re-extended;
            //  a GPU that truncated forward would show up here as a disagreement, so 100% here
            //  independently rules out a "truncate long alignments" bug.
            int cpu_long = (path.aepos + path.bepos - anti) >= 45;
            if (cpu_long) { longset++; if (dae==0 && dbe==0) longset_agree++; }
            if (dae==0 && dbe==0) exact++;
            if (abs(dae)<=1 && abs(dbe)<=1) within1++;
            else
              { //  Attribute the mismatch.  Local_Alignment re-centers SHORT alignments (fshort,
                //  align.c:1535,1553-1576, DUB_TRIM=45) as a FULL-routine post-process -- out of
                //  the forward-wave port's scope (WAVE_PORT_NOTES "out of scope").  fshort is
                //  tested on the FORWARD endpoint *before* reverse, which is exactly the GPU's
                //  output, so use it as the principled detector.
                if (!cpu_long) mism_cpu_short++;
                if ((gAe[j] + gBe[j]) - anti < 45) fshort_mismatch++;
                else
                  { int md = abs(dae) > abs(dbe) ? abs(dae) : abs(dbe);   // genuine port residual
                    if (md > maxdelta) maxdelta = md;
                    if (md <= 5) mism_le5++; else mism_gt5++;
                  }
              }
          }
        double c1 = now();

        Free_Work_Data(work);
        Free_Align_Spec(spec);

        long comparable = SAMP - overflow;
        printf("--fwd-validate results (%d seeds, %ld comparable, %ld band-overflow):\n",
               SAMP,comparable,overflow);
        printf("  forward-endpoint EXACT (ae&be)        : %ld/%ld  (%.2f%%)\n",
               exact,comparable, comparable?100.0*exact/comparable:0.0);
        printf("  forward-endpoint within 1bp (ae&be)   : %ld/%ld  (%.2f%%)\n",
               within1,comparable, comparable?100.0*within1/comparable:0.0);
        printf("  of the %ld non-agreeing: fshort(full-LA re-center, out-of-scope)=%ld;"
               "  genuine port residual: |delta|<=5=%ld, |delta|>5=%ld, max|delta|=%ld\n",
               comparable-within1, fshort_mismatch, mism_le5, mism_gt5, maxdelta);
        printf("  => port-attributable endpoint agreement (excl. fshort re-center): "
               "%ld/%ld (%.2f%%)\n",
               comparable-(mism_le5+mism_gt5), comparable,
               comparable?100.0*(comparable-(mism_le5+mism_gt5))/comparable:0.0);
        printf("  [non-circular] CPU-long subset (final aepos+bepos-anti>=45) EXACT: "
               "%ld/%ld (%.2f%%);  mismatches with CPU also short: %ld/%ld\n",
               longset_agree,longset, longset?100.0*longset_agree/longset:0.0,
               mism_cpu_short, comparable-within1);
        printf("  GPU batch %.4fs (%.0f seeds/s), host oracle %.4fs\n",
               g1-g0, SAMP/(g1-g0), c1-c0);

        free(WS); free(sidx); free(gAe); free(gBe); free(gFd); free(tab); free(scr);
      }

    //  ---- Task 6 --discover-validate: FULL fwd+rev discovery vs the CPU .cpuref (Task 3) ----
    //  THE definitive Stage-1 correctness gate: over the ENTIRE seed distribution, compare the
    //  GPU's ab/ae/bb/be/diffs (wave_discover_batch) to Local_Alignment's own output recorded in
    //  the .cpuref.  same-score (|Δdiffs|==0) is the faithfulness proof the forward-only endpoint
    //  check could not do; the depth-work ratio proves the GPU does ~the same wave work.

    if (discvalidate)
      { //  derive the .cpuref path from the .seeds path: strip ".seeds", append ".cpuref"
        char cpuref_path[4096];
        { size_t L = strlen(argv[2]);
          const char *suf = ".seeds";  size_t sl = strlen(suf);
          if (L > sl && strcmp(argv[2]+L-sl,suf) == 0)
            { memcpy(cpuref_path,argv[2],L-sl); cpuref_path[L-sl] = 0; }
          else
            { memcpy(cpuref_path,argv[2],L); cpuref_path[L] = 0; }
          strcat(cpuref_path,".cpuref");
        }

        FILE *cf = fopen(cpuref_path,"rb");
        CpuRefHeader ch;
        if (cf == NULL || fread(&ch,sizeof(ch),1,cf) != 1 || ch.magic != WAVE_CPUREF_MAGIC)
          { fprintf(stderr,"%s: cannot open/parse cpuref %s (run wave_bench_cpu first)\n",
                    argv[0],cpuref_path);
            if (cf) fclose(cf);
          }
        else if ((int) ch.nrecs != N)
          { fprintf(stderr,"%s: cpuref has %u recs, .seeds has %d -- mismatch\n",
                    argv[0],ch.nrecs,N);
            fclose(cf);
          }
        else
          { CpuRefRec *ref = (CpuRefRec *) malloc((size_t) N * sizeof(CpuRefRec));
            if (ref == NULL || fread(ref,sizeof(CpuRefRec),(size_t) N,cf) != (size_t) N)
              { fprintf(stderr,"%s: short read of %d cpuref recs\n",argv[0],N); fclose(cf); }
            else
              { fclose(cf);

                wave_seed *WS = (wave_seed *) malloc((size_t) N * sizeof(wave_seed));
                int *gAb=(int*)malloc((size_t)N*sizeof(int)), *gAe=(int*)malloc((size_t)N*sizeof(int));
                int *gBb=(int*)malloc((size_t)N*sizeof(int)), *gBe=(int*)malloc((size_t)N*sizeof(int));
                int *gDf=(int*)malloc((size_t)N*sizeof(int));
                int i;
                for (i = 0; i < N; i++)
                  { SeedRec *s = &S[i];
                    WS[i].aread = s->aread;  WS[i].bread = s->bread;
                    WS[i].alen  = s->alen;   WS[i].blen  = s->blen;
                    WS[i].anti  = s->seed_anti; WS[i].diag = s->seed_diag;
                    WS[i].comp  = (s->flags & 0x1) ? 1 : 0;   // COMP_FLAG
                  }

                int16_t *tab = (int16_t *) malloc((size_t) WV_TRIM_SZ * sizeof(int16_t));
                int16_t *scr = (int16_t *) malloc((size_t) WV_TRIM_SZ * sizeof(int16_t));
                int path_ave = build_trim_tables(gdb1->freq,0.7,tab,scr);

                printf("\n--discover-validate: FULL distribution, %d seeds, path_ave=%d, "
                       "cpuref=%s\n",N,path_ave,cpuref_path);

                double g0 = now();
                if (wave_discover_batch(g,N,WS,tab,scr,path_ave,gAb,gAe,gBb,gBe,gDf) != 0)
                  fprintf(stderr,"%s: wave_discover_batch failed\n",argv[0]);
                double g1 = now();

                //  ---- compare ----
                long comparable=0, overflow=0;
                long samescore=0, exact4=0, within2=0;
                double ratio_sum=0; long ratio_n=0;
                //  breakdown of same-score MISMATCHES (attribution):
                long ss_mm=0, ss_mm_cpushort=0;   // cpu short => fshort/rshort re-center (out of scope)
                long dmax=0; long dsum_abs=0;

                for (i = 0; i < N; i++)
                  { if (gDf[i] == -2) { overflow++; continue; }
                    comparable++;
                    CpuRefRec *r = &ref[i];

                    if (gDf[i] == r->diffs) samescore++;
                    else
                      { ss_mm++;
                        long dd = labs((long) gDf[i] - (long) r->diffs);
                        dsum_abs += dd; if (dd > dmax) dmax = dd;
                        //  CPU alignment short at either end => Local_Alignment re-centered/re-
                        //  extended (align.c:1535,1551-1576), which the endpoint-only port omits.
                        int cpu_fshort = (r->aepos + r->bepos - S[i].seed_anti) < 45;
                        int cpu_rshort = (S[i].seed_anti - (r->abpos + r->bbpos)) < 45;
                        if (cpu_fshort || cpu_rshort) ss_mm_cpushort++;
                      }

                    int e_ab = abs(gAb[i]-r->abpos), e_ae = abs(gAe[i]-r->aepos);
                    int e_bb = abs(gBb[i]-r->bbpos), e_be = abs(gBe[i]-r->bepos);
                    int emax = e_ab; if(e_ae>emax)emax=e_ae; if(e_bb>emax)emax=e_bb; if(e_be>emax)emax=e_be;
                    if (emax == 0) exact4++;
                    if (emax <= 2) within2++;

                    if (r->diffs > 0) { ratio_sum += (double) gDf[i] / r->diffs; ratio_n++; }
                  }

                double ss_pct   = comparable ? 100.0*samescore/comparable : 0.0;
                double ex_pct   = comparable ? 100.0*exact4/comparable     : 0.0;
                double w2_pct   = comparable ? 100.0*within2/comparable     : 0.0;
                double ratio_mn = ratio_n ? ratio_sum/ratio_n : 0.0;

                printf("--discover-validate results (%d seeds, %ld comparable, %ld band-overflow):\n",
                       N,comparable,overflow);
                printf("  same-score  (|Δdiffs|==0)            : %ld/%ld  (%.3f%%)   [GATE >= 95%%]\n",
                       samescore,comparable,ss_pct);
                printf("  depth-work ratio mean(GPUdiff/CPUdiff): %.4f            [GATE in 0.97..1.05]\n",
                       ratio_mn);
                printf("  exact endpoint (all 4, tol 0)         : %ld/%ld  (%.3f%%)\n",
                       exact4,comparable,ex_pct);
                printf("  endpoint within 2bp (all 4)           : %ld/%ld  (%.3f%%)\n",
                       within2,comparable,w2_pct);
                printf("  of %ld same-score mismatches: %ld have a SHORT CPU end (fshort/rshort "
                       "full-LA re-center, out-of-scope); mean|Δdiff|=%.2f max|Δdiff|=%ld\n",
                       ss_mm, ss_mm_cpushort, ss_mm?(double)dsum_abs/ss_mm:0.0, dmax);
                printf("  => port-attributable same-score (excl. CPU-short re-center): %ld/%ld (%.3f%%)\n",
                       comparable-(ss_mm-ss_mm_cpushort), comparable,
                       comparable?100.0*(comparable-(ss_mm-ss_mm_cpushort))/comparable:0.0);
                printf("  GPU discover batch %.4fs (%.0f seeds/s)\n", g1-g0, N/(g1-g0));

                int gate_ss = ss_pct >= 95.0;
                int gate_dr = (ratio_mn >= 0.97 && ratio_mn <= 1.05);
                printf("  GATE: same-score>=95%% %s ; depth-ratio in[0.97,1.05] %s => %s\n",
                       gate_ss?"PASS":"FAIL", gate_dr?"PASS":"FAIL",
                       (gate_ss && gate_dr)?"PASS":"FAIL");

                free(WS); free(gAb); free(gAe); free(gBb); free(gBe); free(gDf);
                free(tab); free(scr); free(ref);
              }
          }
      }

    //  ---- Task 7 --stage1-fit: THE headline benchmark -- full-distribution GPU-vs-CPU (both
    //  timing bases), the per-(band,depth) fit map, and the occupancy-based mechanism
    //  decomposition (concurrent warps vs the CPU's 32 cores). Self-contained: loads its own
    //  .cpuref (does not require --discover-validate to have also been passed).

    if (stage1fit)
      { char cpuref_path[4096];
        cpuref_path_from_seeds(argv[2],cpuref_path,sizeof(cpuref_path));

        FILE *cf = fopen(cpuref_path,"rb");
        CpuRefHeader ch;
        if (cf == NULL || fread(&ch,sizeof(ch),1,cf) != 1 || ch.magic != WAVE_CPUREF_MAGIC)
          { fprintf(stderr,"%s: --stage1-fit: cannot open/parse cpuref %s (run wave_bench_cpu first)\n",
                    argv[0],cpuref_path);
            if (cf) fclose(cf);
          }
        else if ((int) ch.nrecs != N)
          { fprintf(stderr,"%s: --stage1-fit: cpuref has %u recs, .seeds has %d -- mismatch\n",
                    argv[0],ch.nrecs,N);
            fclose(cf);
          }
        else
          { CpuRefRec *ref = (CpuRefRec *) malloc((size_t) N * sizeof(CpuRefRec));
            if (ref == NULL || fread(ref,sizeof(CpuRefRec),(size_t) N,cf) != (size_t) N)
              { fprintf(stderr,"%s: --stage1-fit: short read of %d cpuref recs\n",argv[0],N); if(cf) fclose(cf); }
            else
              { fclose(cf);

                wave_seed *WS = (wave_seed *) malloc((size_t) N * sizeof(wave_seed));
                int *gAb=(int*)malloc((size_t)N*sizeof(int)), *gAe=(int*)malloc((size_t)N*sizeof(int));
                int *gBb=(int*)malloc((size_t)N*sizeof(int)), *gBe=(int*)malloc((size_t)N*sizeof(int));
                int *gDf=(int*)malloc((size_t)N*sizeof(int));
                int i;
                for (i = 0; i < N; i++)
                  { SeedRec *s = &S[i];
                    WS[i].aread = s->aread;  WS[i].bread = s->bread;
                    WS[i].alen  = s->alen;   WS[i].blen  = s->blen;
                    WS[i].anti  = s->seed_anti; WS[i].diag = s->seed_diag;
                    WS[i].comp  = (s->flags & 0x1) ? 1 : 0;
                  }

                int16_t *tab = (int16_t *) malloc((size_t) WV_TRIM_SZ * sizeof(int16_t));
                int16_t *scr = (int16_t *) malloc((size_t) WV_TRIM_SZ * sizeof(int16_t));
                int path_ave = build_trim_tables(gdb1->freq,0.7,tab,scr);

                printf("\n--stage1-fit: FULL distribution, %d seeds, path_ave=%d, cpuref=%s\n",
                       N,path_ave,cpuref_path);

                //  ---- mechanism decomposition: occupancy -> concurrent warps ----
                int maxBlocksFwd=0, maxBlocksRev=0, smCount=0;
                wave_query_occupancy(&maxBlocksFwd,&maxBlocksRev,&smCount);
                int warpsPerSM = maxBlocksFwd < maxBlocksRev ? maxBlocksFwd : maxBlocksRev;
                long concurrentWarps = (long) warpsPerSM * smCount;
                printf("\noccupancy: forward_sweep_warp %d blocks/SM, reverse_sweep_warp %d blocks/SM"
                       " (block=32 threads=1 warp), %d SMs -> concurrent warps = min(%d,%d)*%d = %ld\n",
                       maxBlocksFwd,maxBlocksRev,smCount,maxBlocksFwd,maxBlocksRev,smCount,concurrentWarps);

                //  ---- headline: full-distribution, CUDA-event timed, best-of-N, both bases ----
                #define STAGE1_NREPS 5
                float h2d_ms[STAGE1_NREPS], ker_ms[STAGE1_NREPS], d2h_ms[STAGE1_NREPS];
                int rep;
                for (rep = 0; rep < STAGE1_NREPS; rep++)
                  { if (wave_discover_batch_timed(g,N,WS,tab,scr,path_ave,gAb,gAe,gBb,gBe,gDf,
                                                   &h2d_ms[rep],&ker_ms[rep],&d2h_ms[rep]) != 0)
                      fprintf(stderr,"%s: wave_discover_batch_timed failed on rep %d\n",argv[0],rep);
                  }
                float best_ker = ker_ms[0], best_tot = h2d_ms[0]+ker_ms[0]+d2h_ms[0];
                int   best_tot_rep = 0;
                for (rep = 1; rep < STAGE1_NREPS; rep++)
                  { if (ker_ms[rep] < best_ker) best_ker = ker_ms[rep];
                    float tot = h2d_ms[rep]+ker_ms[rep]+d2h_ms[rep];
                    if (tot < best_tot) { best_tot = tot; best_tot_rep = rep; }
                  }
                double gpu_alnps_i  = N / (best_ker/1000.0);    // basis (i): kernel-only (best single component across reps)
                double gpu_alnps_ii = N / (best_tot/1000.0);    // basis (ii): +per-batch H2D/D2H (best single REP's total)

                printf("\nGPU wave_discover_batch, full %d seeds, best-of-%d CUDA-event timing:\n",
                       N,STAGE1_NREPS);
                printf("  basis (i)  wave-engine/kernel-only : %7.1f ms  (%9.1f seeds/s)  "
                       "[genome H2D was one-time, reported above, NOT included]\n",
                       best_ker, gpu_alnps_i);
                printf("  basis (ii) realistic (+seed H2D+result D2H): %7.1f ms  (%9.1f seeds/s)  "
                       "(rep %d: h2d=%.2fms kernel=%.1fms d2h=%.2fms)\n",
                       best_tot, gpu_alnps_ii, best_tot_rep,
                       h2d_ms[best_tot_rep], ker_ms[best_tot_rep], d2h_ms[best_tot_rep]);

                //  ---- correctness re-check vs .cpuref (mirrors --discover-validate's gate) ----
                long comparable=0, overflow=0, samescore=0;
                double ratio_sum=0; long ratio_n=0;
                for (i = 0; i < N; i++)
                  { if (gDf[i] == -2) { overflow++; continue; }
                    comparable++;
                    if (gDf[i] == ref[i].diffs) samescore++;
                    if (ref[i].diffs > 0) { ratio_sum += (double) gDf[i] / ref[i].diffs; ratio_n++; }
                  }
                double ss_pct = comparable ? 100.0*samescore/comparable : 0.0;
                double ratio_mn = ratio_n ? ratio_sum/ratio_n : 0.0;
                printf("\ncorrectness re-check (matches --discover-validate's gate): same-score %ld/%ld"
                       " (%.3f%%), depth-ratio %.4f, band-overflow %ld\n",
                       samescore,comparable,ss_pct,ratio_mn,overflow);

                //  ---- fit map: bucket every comparable seed by the GPU's OWN (band,depth),
                //  same buckets as wave_bench_cpu.c, then re-run wave_discover_batch_timed on
                //  each bucket's subset (best-of-3) to get a directly comparable per-bucket
                //  GPU aln/s table (kernel-only basis, matching the CPU table's pure-compute
                //  methodology). ----
                std::vector<int> bktIdx[G_NBAND][G_NDEPTH];
                for (i = 0; i < N; i++)
                  { if (gDf[i] == -2) continue;
                    long band  = labs((long)(gAe[i]-gAb[i]) - (long)(gBe[i]-gBb[i]));
                    long depth = gDf[i];
                    int bb = gbucket_of(band,G_BAND_HI,G_NBAND);
                    int db = gbucket_of(depth,G_DEPTH_HI,G_NDEPTH);
                    bktIdx[bb][db].push_back(i);
                  }

                #define STAGE1_BREPS 3
                double bkt_alnps[G_NBAND][G_NDEPTH];
                long   bkt_n[G_NBAND][G_NDEPTH];
                { int bb, db;
                  for (bb = 0; bb < G_NBAND; bb++)
                    for (db = 0; db < G_NDEPTH; db++)
                      { std::vector<int> &idx = bktIdx[bb][db];
                        long cnt = (long) idx.size();
                        bkt_n[bb][db] = cnt;
                        if (cnt == 0) { bkt_alnps[bb][db] = 0; continue; }

                        wave_seed *sub = (wave_seed *) malloc((size_t) cnt * sizeof(wave_seed));
                        int *sab=(int*)malloc((size_t)cnt*sizeof(int)), *sae=(int*)malloc((size_t)cnt*sizeof(int));
                        int *sbb=(int*)malloc((size_t)cnt*sizeof(int)), *sbe=(int*)malloc((size_t)cnt*sizeof(int));
                        int *sdf=(int*)malloc((size_t)cnt*sizeof(int));
                        long k;
                        for (k = 0; k < cnt; k++) sub[k] = WS[idx[k]];

                        float bh[STAGE1_BREPS], bk[STAGE1_BREPS], bd[STAGE1_BREPS];
                        int r;
                        for (r = 0; r < STAGE1_BREPS; r++)
                          wave_discover_batch_timed(g,(int) cnt,sub,tab,scr,path_ave,
                                                     sab,sae,sbb,sbe,sdf,&bh[r],&bk[r],&bd[r]);
                        float bestk = bk[0];
                        for (r = 1; r < STAGE1_BREPS; r++) if (bk[r] < bestk) bestk = bk[r];
                        bkt_alnps[bb][db] = cnt / (bestk/1000.0);

                        free(sub); free(sab); free(sae); free(sbb); free(sbe); free(sdf);
                      }
                }

                printf("\nfit map: GPU aln/s (kernel-only basis), same (band,depth) buckets as "
                       "wave_bench_cpu.c\n");
                printf("%14s","depth\\band");
                { int bb;
                  for (bb=0; bb<G_NBAND; bb++)
                    printf("  <=%-8ld", G_BAND_HI[bb]==INT32_MAX?999999999L:G_BAND_HI[bb]);
                  printf("\n");
                }
                { int db;
                  for (db=0; db<G_NDEPTH; db++)
                    { printf("%12s<=%-1ld", "", G_DEPTH_HI[db]==INT32_MAX?999999999L:G_DEPTH_HI[db]);
                      int bb;
                      for (bb=0; bb<G_NBAND; bb++)
                        { if (bkt_n[bb][db] > 0)
                            printf("  %10.1f", bkt_alnps[bb][db]);
                          else
                            printf("  %10s", "-");
                        }
                      printf("   (n=");
                      long tot=0; for (bb=0; bb<G_NBAND; bb++) tot += bkt_n[bb][db];
                      printf("%ld)%s\n", tot, db==G_NDEPTH-1 ? "   <-- DIAGNOSTIC (deep tail), NOT the headline" : "");
                    }
                }

                //  ---- mechanism decomposition: per-warp vs per-core ----
                double gpu_per_warp = concurrentWarps > 0 ? gpu_alnps_i / concurrentWarps : 0.0;
                printf("\nmechanism decomposition (basis (i), kernel-only):\n");
                printf("  GPU: %.1f aln/s aggregate / %ld concurrent warps = %.4f aln/s per warp\n",
                       gpu_alnps_i, concurrentWarps, gpu_per_warp);
                printf("  CPU: %.1f aln/s per core (1-thread measured)\n", cpu_1core_alnps);
                printf("  per-core / per-warp = %.1fx  (%ld warps at %.4f aln/s/warp vs 32 cores at "
                       "%.1f aln/s/core is how the aggregate ratio arises: many-but-slow vs few-but-fast)\n",
                       gpu_per_warp>0 ? cpu_1core_alnps/gpu_per_warp : 0.0, concurrentWarps, gpu_per_warp,
                       cpu_1core_alnps);

                //  ---- final headline summary ----
                printf("\n=== HEADLINE: full distribution, %d seeds (GPU vs CPU) ===\n",N);
                printf("  GPU basis (i)  kernel-only : %10.1f seeds/s\n", gpu_alnps_i);
                printf("  GPU basis (ii) realistic   : %10.1f seeds/s\n", gpu_alnps_ii);
                printf("  CPU 1-core                 : %10.1f aln/s\n", cpu_1core_alnps);
                printf("  CPU 32-core                : %10.1f aln/s\n", cpu_32core_alnps);
                printf("  speedup vs 1-core : basis(i) %6.2fx   basis(ii) %6.2fx\n",
                       gpu_alnps_i/cpu_1core_alnps, gpu_alnps_ii/cpu_1core_alnps);
                printf("  speedup vs 32-core: basis(i) %6.2fx   basis(ii) %6.2fx\n",
                       gpu_alnps_i/cpu_32core_alnps, gpu_alnps_ii/cpu_32core_alnps);
                printf("  CAVEAT: this Stage-1 comparison is discovery (wave sweep: endpoints+diffs)"
                       " ONLY -- no trace-points -- whereas CPU Local_Alignment ALSO emits trace in the"
                       " same pass, so this slightly favors the GPU. Trace-inclusive: Stage 2 (Tasks 8-9).\n");

                free(WS); free(gAb); free(gAe); free(gBb); free(gBe); free(gDf);
                free(tab); free(scr); free(ref);
              }
          }
      }

    free(S);
  }

  wave_close(g);
  free(Aflat); free(Abase);
  free(Bfwd);  free(Brev);  free(Bbase);

  return fail;
}
