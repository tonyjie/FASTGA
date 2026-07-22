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

static char *Usage = "<in:path>[.1aln] <in:path.seeds> [--selftest]";

int main(int argc, char *argv[])
{ GDB      _gdb1, *gdb1 = &_gdb1;
  GDB      _gdb2, *gdb2 = &_gdb2;
  FILE    **units1, **units2;
  OneFile  *input;
  int64     novl;
  int       TSPACE, ISTWO;
  int       selftest = 0;

  if (argc < 3)
    { fprintf(stderr,"Usage: %s %s\n",argv[0],Usage); return 1; }
  { int i;
    for (i = 3; i < argc; i++)
      if (strcmp(argv[i],"--selftest") == 0)
        selftest = 1;
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
    free(S);
  }

  wave_close(g);
  free(Aflat); free(Abase);
  free(Bfwd);  free(Brev);  free(Bbase);

  return fail;
}
