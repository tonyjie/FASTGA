/*******************************************************************************************
 *
 *  wave_bench_cpu -- CPU baseline for the GPU wave-parallelism characterization study
 *                    (Task 3).  Times FastGA's REAL Local_Alignment (the actual hot-path
 *                    aligner, NOT Compute_Alignment) run per-seed on GENOME-RESIDENT contigs,
 *                    over the FULL seed distribution extracted by extract_seeds (Task 2).
 *
 *  Unlike deep_cpu_bench (which reads pre-cut, pre-oriented margin-grown windows out of a
 *  .disc file), this harness loads every contig of BOTH genomes into memory once -- A forward
 *  only, B in BOTH orientations (forward + reverse-complement, built once via Complement_Seq)
 *  -- and then, for every seed, aligns full-contig-to-full-contig starting from the seed's
 *  (diag,anti), exactly as FastGA's align_contigs does for the real pipeline (cross path,
 *  lbord=hbord=-1).
 *
 *  Why both B orientations resident instead of complementing per seed: a chromosome-scale
 *  contig is hit by many thousands of seeds, and Complement_Seq is O(contig length) -- paying
 *  that cost once per contig (not once per seed) is the only way "genome-resident" is also
 *  fast.  It costs another ~1x B-genome's worth of RAM, which is negligible here (single-digit
 *  GB for human-scale genomes on this node).
 *
 *  Produces:
 *    - <out>.cpuref: one CpuRefRec per seed (same order as the .seeds file), the reference
 *      dump Tasks 6/7 validate the GPU wave kernel against (see wave_harness.h).
 *    - stdout: seed count, 1-core and N-core aln/s, endpoint-agreement vs the .1aln reference,
 *      and a (band-proxy x depth) stratified aln/s table.
 *
 *  Build: see Makefile target wave_bench_cpu (mirrors deep_cpu_bench's link line, + -fopenmp).
 *  Run:   ./gpu/wave_bench_cpu <in.1aln> <in.seeds> <out.cpuref> [nthreads=32]
 *         SKIP_SINGLE=1 skips the single-thread reference pass (and hence .cpuref + agreement
 *         reporting) for fast N-thread-only reruns, exactly as in deep_cpu_bench.
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>

#include "GDB.h"
#include "align.h"
#include "alncode.h"
#include "wave_harness.h"

typedef unsigned char u8;

static char *Usage = "<in:path>[.1aln] <in:path.seeds> <out:path.cpuref> [nthreads=32]";

//  Stratification buckets: band proxy = |(ae-ab)-(be-bb)| (net drift of the FOUND alignment,
//  a cheap proxy for the true active diagonal band the wave had to sweep -- see brief), and
//  depth = path.diffs (edit-distance work actually done).  Both log-scaled.

#define NBAND 5
#define NDEPTH 5
static long BAND_HI[NBAND]   = {   8,   32,  128,   512, INT32_MAX };
static long DEPTH_HI[NDEPTH] = {  10,   50,  200,  1000, INT32_MAX };

static int bucket_of(long v, long *hi, int n)
{ int i;
  for (i = 0; i < n; i++)
    if (v <= hi[i])
      return i;
  return n-1;
}

int main(int argc, char *argv[])
{ GDB      _gdb1, *gdb1 = &_gdb1;
  GDB      _gdb2, *gdb2 = &_gdb2;
  FILE    **units1, **units2;
  OneFile  *input;
  int64     novl;
  int       TSPACE, ISTWO;
  int       NT;

  if (argc < 4 || argc > 5)
    { fprintf(stderr,"Usage: %s %s\n",argv[0],Usage);
      return 1;
    }
  NT = argc >= 5 ? atoi(argv[4]) : 32;

  //  ---- Open the .1aln just to get gdb1/gdb2 (paths + skeletons), exactly as
  //  extract_seeds.c/extract_deep.c do.  We never read the .1aln's overlap records --
  //  the seed list comes from the separate .seeds file. ----

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
  oneFileClose(input);   //  done with the .1aln itself; only needed it for the GDB paths

  //  ---- Load the seed list ----

  int       N;
  SeedRec  *S;

  { FILE *sf = fopen(argv[2],"rb");
    WaveSeedHeader hdr;

    if (sf == NULL) { fprintf(stderr,"%s: cannot open %s\n",argv[0],argv[2]); return 1; }
    if (fread(&hdr,sizeof(hdr),1,sf) != 1 || hdr.magic != WAVE_SEEDS_MAGIC)
      { fprintf(stderr,"%s: bad seeds header in %s\n",argv[0],argv[2]); return 1; }

    N = (int) hdr.nseeds;
    S = malloc((size_t) N * sizeof(SeedRec));
    if (S == NULL) { fprintf(stderr,"%s: out of memory for %d seeds\n",argv[0],N); return 1; }

    if (fread(S,sizeof(SeedRec),(size_t) N,sf) != (size_t) N)
      { fprintf(stderr,"%s: short read of seed records\n",argv[0]); return 1; }
    fclose(sf);
  }
  printf("loaded %d seeds\n",N);

  //  ---- Genome residency: A forward per contig; B forward AND reverse-complement per
  //  contig (built once each, so per-seed orientation is a pointer pick, not a copy). ----

  double tload0 = omp_get_wtime();

  char **Aseq = malloc((size_t) gdb1->ncontig * sizeof(char*));
  { int i;
    for (i = 0; i < gdb1->ncontig; i++)
      { int64 clen = gdb1->contigs[i].clen;
        //  +8, not +2: mirrors New_Contig_Buffer's own padding (maxctg+8) -- Uncompress_Read's
        //  4-base unpacking can write a few bytes past s[len] for short/misaligned contigs,
        //  so a tight clen+2 allocation heap-corrupts. See GDB.c New_Contig_Buffer.
        char *buf = malloc((size_t) clen + 8);
        if (buf == NULL) { fprintf(stderr,"%s: OOM loading A contig %d\n",argv[0],i); return 1; }
        Aseq[i] = Get_Contig(gdb1,i,NUMERIC,buf+1);
      }
  }

  char **Bfwd = malloc((size_t) gdb2->ncontig * sizeof(char*));
  char **Brev = malloc((size_t) gdb2->ncontig * sizeof(char*));
  { int i;
    for (i = 0; i < gdb2->ncontig; i++)
      { int64 clen = gdb2->contigs[i].clen;
        char *fbuf = malloc((size_t) clen + 8);   //  +8: see Aseq loop comment above
        char *rbuf = malloc((size_t) clen + 8);
        if (fbuf == NULL || rbuf == NULL)
          { fprintf(stderr,"%s: OOM loading B contig %d\n",argv[0],i); return 1; }
        Bfwd[i] = Get_Contig(gdb2,i,NUMERIC,fbuf+1);
        //  Copy the meaningful sentinel-padded region (clen+2 bytes: [-1]..[clen]) then
        //  complement in place -- Complement_Seq only touches [0,clen), leaving the
        //  sentinels (already 4 at both ends, copied verbatim) untouched.  The extra
        //  padding bytes beyond [clen] are never read by anything downstream.
        memcpy(rbuf,Bfwd[i]-1,(size_t) clen + 2);
        Brev[i] = rbuf+1;
        Complement_Seq(Brev[i],(int) clen);
      }
  }
  double tload = omp_get_wtime()-tload0;
  printf("genome-resident: A %d contigs (%lld bp), B %d contigs (%lld bp, fwd+rev) -- loaded in %.1fs\n",
         gdb1->ncontig,(long long) gdb1->seqtot,gdb2->ncontig,(long long) gdb2->seqtot,tload);

  //  ---- .cpuref output (one record per seed, same order as .seeds) ----

  CpuRefRec *ref = malloc((size_t) N * sizeof(CpuRefRec));
  if (ref == NULL) { fprintf(stderr,"%s: OOM for %d cpuref records\n",argv[0],N); return 1; }

  //  Stratification accumulators (single-thread pass only -- see header comment).
  long   bkt_cnt[NBAND][NDEPTH];
  double bkt_time[NBAND][NDEPTH];
  memset(bkt_cnt,0,sizeof(bkt_cnt));
  memset(bkt_time,0,sizeof(bkt_time));

  //  ---- single-thread: per-core rate + the authoritative .cpuref/agreement pass
  //  (skip with SKIP_SINGLE=1 for fast multi-thread-only reruns, as in deep_cpu_bench) ----

  if (getenv("SKIP_SINGLE") == NULL)
    { Align_Spec *spec = New_Align_Spec(0.7,100,gdb1->freq,0);
      Work_Data  *work = New_Work_Data();
      Alignment   aln;
      Path        path;
      long        e50 = 0;
      double      sumdiff_ratio = 0;
      int         i;

      aln.path = &path;
      double t0 = omp_get_wtime();
      for (i = 0; i < N; i++)
        { SeedRec *s = &S[i];

          aln.aseq  = Aseq[s->aread];
          aln.bseq  = COMP(s->flags) ? Brev[s->bread] : Bfwd[s->bread];
          aln.alen  = s->alen;
          aln.blen  = s->blen;
          aln.flags = 0;   //  orientation already resolved by the Bfwd/Brev pick above

          double ts = omp_get_wtime();
          int dg = s->seed_diag, anti = s->seed_anti;
          if (Local_Alignment(&aln,work,spec,dg,dg,anti,-1,-1))
            fprintf(stderr,"%s: LA err on seed %d\n",argv[0],i);
          double te = omp_get_wtime();

          ref[i].abpos = path.abpos;  ref[i].aepos = path.aepos;
          ref[i].bbpos = path.bbpos;  ref[i].bepos = path.bepos;
          ref[i].diffs = path.diffs;

          int mx = abs(path.abpos - s->ref_ab);
          int t;
          t = abs(path.aepos - s->ref_ae); if (t>mx) mx=t;
          t = abs(path.bbpos - s->ref_bb); if (t>mx) mx=t;
          t = abs(path.bepos - s->ref_be); if (t>mx) mx=t;
          if (mx <= 50) e50++;
          if (s->ref_diffs > 0) sumdiff_ratio += (double) path.diffs / s->ref_diffs;

          long band  = labs((long)(path.aepos-path.abpos) - (long)(path.bepos-path.bbpos));
          long depth = path.diffs;
          int  bb = bucket_of(band,BAND_HI,NBAND);
          int  db = bucket_of(depth,DEPTH_HI,NDEPTH);
          bkt_cnt[bb][db]++;
          bkt_time[bb][db] += (te-ts);
        }
      double s_elapsed = omp_get_wtime()-t0;

      printf("\nCPU  1-thread Local_Alignment: %.3f s  (%.1f aln/s per core)\n",
             s_elapsed, N/s_elapsed);
      printf("endpoint agreement vs .1aln (contig-relative) within 50bp: %ld/%d (%.1f%%);"
             "  mean diffs(LA)/diffs(ref)=%.3f\n",
             e50, N, 100.0*e50/N, sumdiff_ratio/N);

      //  ---- stratified aln/s table ----
      printf("\nstratified aln/s  (rows=depth(path.diffs) bucket, cols=band-proxy"
             " |Δa-Δb| bucket; band proxy for the true active band, see brief)\n");
      printf("%14s","depth\\band");
      { int bb;
        for (bb=0; bb<NBAND; bb++)
          printf("  <=%-8ld", BAND_HI[bb]==INT32_MAX?999999999L:BAND_HI[bb]);
        printf("\n");
      }
      { int db;
        for (db=0; db<NDEPTH; db++)
          { printf("%12s<=%-1ld", "", DEPTH_HI[db]==INT32_MAX?999999999L:DEPTH_HI[db]);
            int bb;
            for (bb=0; bb<NBAND; bb++)
              { if (bkt_cnt[bb][db] > 0)
                  printf("  %10.1f", bkt_cnt[bb][db]/bkt_time[bb][db]);
                else
                  printf("  %10s", "-");
              }
            printf("   (n=");
            long tot=0; for (bb=0; bb<NBAND; bb++) tot += bkt_cnt[bb][db];
            printf("%ld)\n",tot);
          }
      }

      //  ---- write .cpuref ----
      { FILE *of = fopen(argv[3],"wb");
        CpuRefHeader h;
        if (of == NULL) { fprintf(stderr,"%s: cannot open %s for writing\n",argv[0],argv[3]); return 1; }
        h.magic = WAVE_CPUREF_MAGIC;
        h.nrecs = (uint32_t) N;
        h.reserved1 = 0; h.reserved2 = 0;
        if (fwrite(&h,sizeof(h),1,of) != 1 ||
            fwrite(ref,sizeof(CpuRefRec),(size_t) N,of) != (size_t) N)
          { fprintf(stderr,"%s: write error to %s\n",argv[0],argv[3]); return 1; }
        fclose(of);
        printf("\nwrote %d cpuref records to %s\n",N,argv[3]);
      }

      Free_Work_Data(work);
      Free_Align_Spec(spec);
    }
  else
    printf("SKIP_SINGLE set: skipping single-thread pass (no .cpuref/agreement this run)\n");

  //  ---- N-thread: real aggregate throughput (dynamic schedule; seeds are wildly
  //  imbalanced in work -- deep waves vs shallow ones) ----

  { omp_set_num_threads(NT);
    double t0 = omp_get_wtime();
    #pragma omp parallel
    { Align_Spec *spec = New_Align_Spec(0.7,100,gdb1->freq,0);
      Work_Data  *work = New_Work_Data();
      Alignment   aln;
      Path        path;
      aln.path = &path;
      #pragma omp for schedule(dynamic,16)
      for (int i = 0; i < N; i++)
        { SeedRec *s = &S[i];
          aln.aseq  = Aseq[s->aread];
          aln.bseq  = COMP(s->flags) ? Brev[s->bread] : Bfwd[s->bread];
          aln.alen  = s->alen;
          aln.blen  = s->blen;
          aln.flags = 0;
          int dg = s->seed_diag, anti = s->seed_anti;
          Local_Alignment(&aln,work,spec,dg,dg,anti,-1,-1);
        }
      Free_Work_Data(work);
      Free_Align_Spec(spec);
    }
    double s_elapsed = omp_get_wtime()-t0;
    printf("CPU %2d-thread Local_Alignment: %.3f s  (%.1f aln/s aggregate)\n",
           NT, s_elapsed, N/s_elapsed);
  }

  return 0;
}
