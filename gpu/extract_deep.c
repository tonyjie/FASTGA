/*******************************************************************************************
 *
 *  extract_disc -- Extract LOCAL-ALIGNMENT DISCOVERY tasks from a FastGA .1aln file into a
 *                   binary "disc task file" for validating a GPU kernel that must rediscover
 *                   alignment endpoints from a seed (see disc_format.h).
 *
 *  Single-threaded.  Reading/setup/orientation logic is copied verbatim from extract_tasks.c
 *  (itself based on the sequence-extraction path of ALNtoPAF.c).  The only difference from
 *  extract_tasks is *what* is dumped per alignment: instead of the exact aligned segments, we
 *  dump a WINDOW grown by `margin` bases on all sides of the aligned extent (clipped to the
 *  contig), a seed anchor at the alignment midpoint (window-relative), and the reference
 *  endpoints (window-relative) + diffs.  A GPU local x-drop wave started at the seed should
 *  rediscover [ref_ab,ref_ae) in A and [ref_bb,ref_be) in B.
 *
 *  Usage: extract_disc <in.1aln> <out.disc> [max_tasks=100000] [margin=64] [wmax=8192]
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>

#include "GDB.h"
#include "align.h"
#include "alncode.h"
#include "disc_format.h"

static char *Usage = "<alignment:path>[.1aln] <out:path> [max_tasks=100000] [margin=64] [wmax=4000000] [dmin=0]";

//  ---------------------------------------------------------------------------------------
//  Self-check edit distance: banded Myers O(ND) furthest-reaching wave, ported from
//  gpu/gpu_align.cu's editdist() (the same algorithm the GPU kernel uses), for verifying
//  that the reference sub-segment really has edit distance == ref_diffs.
//  Returns d>=0, or -1 (d>MAXD), or -2 (window>MAXW).
//  ---------------------------------------------------------------------------------------

#define CHK_MAXW 4096
#define CHK_KBAND 1536
#define CHK_MAXD 8192

static int editdist_check(const unsigned char *A, int m, const unsigned char *B, int n)
{ int drift = m - n;
  int klo = (drift < 0 ? drift : 0) - CHK_KBAND;
  int khi = (drift > 0 ? drift : 0) + CHK_KBAND;
  int width = khi - klo + 1;
  int cur, prv, d, t;
  static int f[2][CHK_MAXW];

  if (width > CHK_MAXW)
    return -2;

  cur = 0;
  prv = 1;
  for (t = 0; t < width; t++)
    f[prv][t] = -1;

  { int x = 0;
    while (x < m && x < n && A[x] == B[x])
      x++;
    for (t = 0; t < width; t++)
      f[cur][t] = -1;
    f[cur][0-klo] = x;
    if (drift == 0 && x >= m)
      return 0;
    t = cur; cur = prv; prv = t;
  }

  for (d = 1; d <= CHK_MAXD; d++)
    { int kmin = klo > -d ? klo : -d;
      int kmax = khi <  d ? khi :  d;
      int k, xe;

      for (t = 0; t < width; t++)
        f[cur][t] = -1;

      for (k = kmin; k <= kmax; k++)
        { int best = -1, v, x, y;

          if (k-1 >= klo) { v = f[prv][k-1-klo]; if (v >= 0 && v+1 > best) best = v+1; }
          if (k+1 <= khi) { v = f[prv][k+1-klo]; if (v >= 0 && v   > best) best = v;   }
          v = f[prv][k-klo];      if (v >= 0 && v+1 > best) best = v+1;
          if (best < 0)
            continue;

          x = best;
          if (x > m) x = m;
          if (x - k > n) x = n + k;
          if (x < 0) x = 0;
          y = x - k;
          while (x < m && y < n && A[x] == B[y])
            { x++; y++; }
          f[cur][k-klo] = x;
        }

      xe = f[cur][drift-klo];
      if (xe >= m && xe - drift >= n)
        return d;
      t = cur; cur = prv; prv = t;
    }

  return -1;
}

int main(int argc, char *argv[])
{ GDB       _gdb1, *gdb1 = &_gdb1;
  GDB       _gdb2, *gdb2 = &_gdb2;
  FILE     **units1;
  FILE     **units2;
  OneFile   *input;
  int64      novl;
  int        TSPACE;
  int        ISTWO;

  int64      max_tasks;
  uint32_t   margin;
  uint32_t   wmax;
  uint32_t   dmin;
  FILE      *out;

  Prog_Name = Strdup("extract_disc","Allocating program name");

  if (argc < 3 || argc > 7)
    { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
      exit (1);
    }

  max_tasks = 100000;
  margin    = 64;
  wmax      = 4000000;
  dmin      = 0;
  if (argc >= 4)
    max_tasks = strtoll(argv[3],NULL,10);
  if (argc >= 5)
    margin = (uint32_t) strtoul(argv[4],NULL,10);
  if (argc >= 6)
    wmax = (uint32_t) strtoul(argv[5],NULL,10);
  if (argc >= 7)
    dmin = (uint32_t) strtoul(argv[6],NULL,10);

  //  Open the .1aln and the two source GDBs, exactly as extract_tasks.c / ALNtoPAF.c's
  //  main does when it needs sequence data (CIGAR || DIFFS path).

  { char *pwd, *root, *cpath;
    char *src1_name, *src2_name;
    char *head, *sptr, *eptr;
    int   s;

    pwd  = PathTo(argv[1]);
    root = Root(argv[1],".1aln");
    input = open_Aln_Read(Catenate(pwd,"/",root,".1aln"),1,&novl,&TSPACE,
                           &src1_name,&src2_name,&cpath);
    if (input == NULL)
      exit (1);
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

    //  Truncate GDB headers to first white-space (mirrors ALNtoPAF; harmless here,
    //  kept for fidelity in case headers are ever needed for debugging).

    head = gdb1->headers;
    for (s = 0; s < gdb1->nscaff; s++)
      { sptr = head + gdb1->scaffolds[s].hoff;
        for (eptr = sptr; *eptr != '\0'; eptr++)
          if (isspace((unsigned char) *eptr))
            break;
        *eptr = '\0';
      }

    if (ISTWO)
      { head = gdb2->headers;
        for (s = 0; s < gdb2->nscaff; s++)
          { sptr = head + gdb2->scaffolds[s].hoff;
            for (eptr = sptr; *eptr != '\0'; eptr++)
              if (isspace((unsigned char) *eptr))
                break;
            *eptr = '\0';
          }
      }

    gdb1->seqs = units1[0];
    gdb2->seqs = units2[0];
  }

  out = fopen(argv[2],"wb");
  if (out == NULL)
    { fprintf(stderr,"%s: Cannot open %s for writing\n",Prog_Name,argv[2]);
      exit (1);
    }

  //  Reserve space for the header, to be filled in at the end.

  { FGADiscHeader hdr;
    memset(&hdr,0,sizeof(hdr));
    if (fwrite(&hdr,sizeof(hdr),1,out) != 1)
      { fprintf(stderr,"%s: Cannot write header to %s\n",Prog_Name,argv[2]);
        exit (1);
      }
  }

  //  Sweep every alignment: build the margin-grown window (A from the full contig buffer,
  //  B via a re-fetched, correctly oriented contig piece covering the window), then dump
  //  the binary disc-task record.  Self-check the first ~200 tasks against a from-scratch
  //  edit-distance recomputation on the reference sub-segment.

  { Overlap   _ovl, *ovl = &_ovl;
    Alignment _aln, *aln = &_aln;
    Path     *path;
    char     *aseq, *bseq;
    int       alast;
    uint32_t  ntasks, nskip;
    int64     i;
    int       nchecked, nmatch;

    aseq = New_Contig_Buffer(gdb1);
    bseq = New_Contig_Buffer(gdb2);
    if (aseq == NULL || bseq == NULL)
      exit (1);

    aln->path = path = &(ovl->path);
    aln->aseq = aseq;

    ntasks = 0;
    nskip  = 0;
    alast  = -1;
    nchecked = 0;
    nmatch   = 0;

    for (i = 0; i < novl; i++)
      { int   acontig, bcontig;
        int   alen, blen;
        int   ab, ae, bb, be;
        int   a0, a1, b0, b1;
        int   aw, bw;
        int   bmin, bmax;
        char *bact;

        Read_Aln_Overlap(input,ovl);
        Skip_Aln_Trace(input);

        if (ntasks >= (uint32_t) max_tasks)
          continue;             //  keep sweeping only to advance the OneFile reader.

        acontig = ovl->aread;
        bcontig = ovl->bread;
        alen = gdb1->contigs[acontig].clen;
        blen = gdb2->contigs[bcontig].clen;
        aln->alen  = alen;
        aln->blen  = blen;
        aln->flags = ovl->flags;

        if ((uint32_t) path->diffs < dmin)   //  keep only deep waves (many differences)
          { nskip++;
            continue;
          }

        ab = path->abpos;  ae = path->aepos;
        bb = path->bbpos;  be = path->bepos;

        a0 = ab - (int) margin;  if (a0 < 0) a0 = 0;
        a1 = ae + (int) margin;  if (a1 > alen) a1 = alen;
        aw = a1 - a0;

        b0 = bb - (int) margin;  if (b0 < 0) b0 = 0;
        b1 = be + (int) margin;  if (b1 > blen) b1 = blen;
        bw = b1 - b0;

        if ((uint32_t) aw > wmax || (uint32_t) bw > wmax)
          { nskip++;
            continue;
          }

        if (acontig != alast)
          Get_Contig(gdb1,acontig,NUMERIC,aseq);
        alast = acontig;

        //  Fetch the B window, oriented exactly as extract_tasks does for the exact
        //  aligned extent, but for the margin-grown window [b0,b1) instead of [bb,be).

        if (COMP(aln->flags))
          { bmin = blen - b1;
            bmax = blen - b0;
          }
        else
          { bmin = b0;
            bmax = b1;
          }

        bact = Get_Contig_Piece(gdb2,bcontig,bmin,bmax,NUMERIC,bseq);
        if (COMP(aln->flags))
          { Complement_Seq(bact,bmax-bmin);
            aln->bseq = bact - (aln->blen-bmax);
          }
        else
          aln->bseq = bact - bmin;

        //  A window = aln->aseq[a0..a1), B window = aln->bseq[b0..b1) -- NUMERIC (0..3),
        //  already oriented/reverse-complemented per FastGA's COMP convention, same
        //  coordinate system extract_tasks uses for abpos/aepos/bbpos/bepos.

        { unsigned char *Awin = (unsigned char *) (aln->aseq + a0);
          unsigned char *Bwin = (unsigned char *) (aln->bseq + b0);
          uint32_t sa, sb;
          int32_t  ref_ab, ref_ae, ref_bb, ref_be;
          uint32_t ref_diffs;
          uint32_t rec[9];

          sa = (uint32_t) ((ab+ae)/2 - a0);
          sb = (uint32_t) ((bb+be)/2 - b0);
          ref_ab = ab - a0;  ref_ae = ae - a0;
          ref_bb = bb - b0;  ref_be = be - b0;
          ref_diffs = (uint32_t) path->diffs;

          rec[0] = (uint32_t) aw;
          rec[1] = (uint32_t) bw;
          rec[2] = sa;
          rec[3] = sb;
          rec[4] = (uint32_t) ref_ab;
          rec[5] = (uint32_t) ref_ae;
          rec[6] = (uint32_t) ref_bb;
          rec[7] = (uint32_t) ref_be;
          rec[8] = ref_diffs;

          if (fwrite(rec,sizeof(rec),1,out) != 1 ||
              fwrite(Awin,1,(size_t) aw,out) != (size_t) aw ||
              fwrite(Bwin,1,(size_t) bw,out) != (size_t) bw)
            { fprintf(stderr,"%s: Write error to %s\n",Prog_Name,argv[2]);
              exit (1);
            }

          //  Self-check: for the first ~200 tasks, verify the invariants and recompute
          //  the edit distance of the reference sub-segment from scratch.

          if (nchecked < 200)
            { int ok = 1;
              int d;

              if (!(sa < (uint32_t) aw)) ok = 0;
              if (!(sb < (uint32_t) bw)) ok = 0;
              if (!(ref_ab >= 0 && ref_ab <= ref_ae && ref_ae <= aw)) ok = 0;
              if (!(ref_bb >= 0 && ref_bb <= ref_be && ref_be <= bw)) ok = 0;
              if (!(sa >= (uint32_t) ref_ab && sa <= (uint32_t) ref_ae)) ok = 0;
              if (!(sb >= (uint32_t) ref_bb && sb <= (uint32_t) ref_be)) ok = 0;

              if (ok)
                { d = editdist_check(Awin+ref_ab,ref_ae-ref_ab,Bwin+ref_bb,ref_be-ref_bb);
                  if (d == (int) ref_diffs)
                    nmatch++;
                }

              nchecked++;
            }
        }

        ntasks++;
      }

    free(bseq-1);
    free(aseq-1);

    //  Patch in the final header.

    { FGADiscHeader hdr;

      hdr.magic    = FGA_DISC_MAGIC;
      hdr.ntasks   = ntasks;
      hdr.margin   = margin;
      hdr.reserved = 0;

      rewind(out);
      if (fwrite(&hdr,sizeof(hdr),1,out) != 1)
        { fprintf(stderr,"%s: Cannot patch header in %s\n",Prog_Name,argv[2]);
          exit (1);
        }
    }

    fclose(out);

    fprintf(stderr,"%s: %u alignments in file, %u tasks written, %u skipped (window > wmax=%u),"
                    " stopped-early=%s\n",
                    Prog_Name,(unsigned) novl,ntasks,nskip,wmax,
                    (ntasks >= (uint32_t) max_tasks && novl > ntasks + nskip) ? "yes" : "no");
    fprintf(stderr,"%s: self-check (edit distance of reference sub-segment == ref_diffs, "
                    "plus bound invariants): %d/%d match (%.2f%%)\n",
                    Prog_Name,nmatch,nchecked,
                    nchecked > 0 ? 100.0*nmatch/nchecked : 0.0);
  }

  if (ISTWO)
    Close_GDB(gdb2);
  Close_GDB(gdb1);

  oneFileClose(input);

  exit (0);
}
