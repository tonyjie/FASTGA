/*******************************************************************************************
 *
 *  extract_trace -- Extract TRACE-POINT validation tasks from a FastGA .1aln file into a
 *                    binary "trace task file" for validating a GPU kernel that must emit
 *                    FastGA trace-points (not just endpoints/diffs) (see trace_format.h).
 *
 *  Single-threaded.  Reading/setup/orientation logic is copied verbatim from extract_disc.c
 *  (itself based on extract_tasks.c / ALNtoPAF.c).  The differences from extract_disc are:
 *
 *    1. NO margin, NO seed, NO window growing -- the exact aligned rectangle
 *       A[abpos..aepos) x B[bbpos..bepos) is dumped.
 *    2. For each alignment, FastGA's reference trace-point vector is (re)computed by
 *       calling Compute_Alignment(aln,work,DIFF_TRACE,TSPACE) with path->abpos/aepos/
 *       bbpos/bepos set to the ORIGINAL contig coordinates (so panels are phased on
 *       global multiples of TSPACE exactly as FastGA does internally -- see FastGA.c's
 *       gpu_align_tube for the same idiom).
 *    3. The reference trace vector (uint16 pairs: diffs,delta-b per TSPACE-panel) is
 *       dumped alongside the rectangle.
 *
 *  Usage: extract_trace <in.1aln> <out.trace> [max_tasks=100000] [wmax=8192]
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
#include "trace_format.h"

static char *Usage = "<alignment:path>[.1aln] <out:path> [max_tasks=100000] [wmax=8192]";

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
  uint32_t   wmax;
  FILE      *out;
  Work_Data *work;

  Prog_Name = Strdup("extract_trace","Allocating program name");

  if (argc < 3 || argc > 5)
    { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
      exit (1);
    }

  max_tasks = 100000;
  wmax      = 8192;
  if (argc >= 4)
    max_tasks = strtoll(argv[3],NULL,10);
  if (argc >= 5)
    wmax = (uint32_t) strtoul(argv[4],NULL,10);

  //  Open the .1aln and the two source GDBs, exactly as extract_disc.c's main does.

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

    //  Truncate GDB headers to first white-space (mirrors extract_disc; harmless here).

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

  { FGATraceHeader hdr;
    memset(&hdr,0,sizeof(hdr));
    if (fwrite(&hdr,sizeof(hdr),1,out) != 1)
      { fprintf(stderr,"%s: Cannot write header to %s\n",Prog_Name,argv[2]);
        exit (1);
      }
  }

  work = New_Work_Data();
  if (work == NULL)
    exit (1);

  //  Sweep every alignment: fetch the exact aligned rectangle (A from the full contig
  //  buffer, B via a re-fetched, correctly oriented contig piece covering [bbpos,bepos)),
  //  recompute the reference trace-point vector with Compute_Alignment, and dump the
  //  binary trace-task record.  Self-check the first ~200 tasks against the trace-point
  //  invariants Check_Trace_Points enforces.

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

        ab = path->abpos;  ae = path->aepos;
        bb = path->bbpos;  be = path->bepos;

        aw = ae - ab;
        bw = be - bb;

        if ((uint32_t) aw > wmax || (uint32_t) bw > wmax)
          { nskip++;
            continue;
          }

        if (acontig != alast)
          Get_Contig(gdb1,acontig,NUMERIC,aseq);
        alast = acontig;

        //  Fetch the B piece for the EXACT aligned extent [bb,be), oriented exactly as
        //  extract_disc does for its margin-grown window, here with margin=0.

        if (COMP(aln->flags))
          { bmin = blen - be;
            bmax = blen - bb;
          }
        else
          { bmin = bb;
            bmax = be;
          }

        bact = Get_Contig_Piece(gdb2,bcontig,bmin,bmax,NUMERIC,bseq);
        if (COMP(aln->flags))
          { Complement_Seq(bact,bmax-bmin);
            aln->bseq = bact - (aln->blen-bmax);
          }
        else
          aln->bseq = bact - bmin;

        //  Set the path's absolute contig coordinates (unchanged from what
        //  Read_Aln_Overlap gave us) so Compute_Alignment indexes aln->aseq/aln->bseq
        //  correctly and phases trace-point panels on global multiples of TSPACE.

        path->abpos = ab;  path->aepos = ae;
        path->bbpos = bb;  path->bepos = be;

        if (Compute_Alignment(aln,work,DIFF_TRACE,TSPACE))
          { fprintf(stderr,"%s: Compute_Alignment failed on overlap %lld\n",
                            Prog_Name,(long long) i);
            exit (1);
          }

        { unsigned char *Aseg = (unsigned char *) (aln->aseq + ab);
          unsigned char *Bseg = (unsigned char *) (aln->bseq + bb);
          uint16        *trace = (uint16 *) path->trace;
          uint32_t       ref_diffs = (uint32_t) path->diffs;
          uint32_t       ref_tlen  = (uint32_t) path->tlen;
          int32_t        rec[4];
          uint32_t       urec[4];

          rec[0] = (int32_t) ab;
          rec[1] = (int32_t) ae;
          rec[2] = (int32_t) bb;
          rec[3] = (int32_t) be;
          urec[0] = (uint32_t) aw;
          urec[1] = (uint32_t) bw;
          urec[2] = ref_diffs;
          urec[3] = ref_tlen;

          if (fwrite(rec,sizeof(rec),1,out) != 1 ||
              fwrite(urec,sizeof(urec),1,out) != 1 ||
              fwrite(Aseg,1,(size_t) aw,out) != (size_t) aw ||
              fwrite(Bseg,1,(size_t) bw,out) != (size_t) bw ||
              fwrite(trace,sizeof(uint16),(size_t) ref_tlen,out) != (size_t) ref_tlen)
            { fprintf(stderr,"%s: Write error to %s\n",Prog_Name,argv[2]);
              exit (1);
            }

          //  Self-check: for the first ~200 tasks, verify the trace-point invariants
          //  that FastGA's Check_Trace_Points enforces.

          if (nchecked < 200)
            { int      ok = 1;
              uint32_t expect_tlen;
              uint32_t sum_diffs, sum_db;
              uint32_t k;

              expect_tlen = 2*(uint32_t) (((ae + TSPACE - 1)/TSPACE) - (ab/TSPACE));
              if (ref_tlen != expect_tlen)
                ok = 0;

              sum_diffs = 0;
              sum_db    = 0;
              for (k = 0; k < ref_tlen; k += 2)
                { sum_diffs += trace[k];
                  sum_db    += trace[k+1];
                }
              if (sum_diffs != ref_diffs)
                ok = 0;
              if (sum_db != (uint32_t) bw)
                ok = 0;

              if (ok)
                nmatch++;
              nchecked++;
            }
        }

        ntasks++;
      }

    free(bseq-1);
    free(aseq-1);

    //  Patch in the final header.

    { FGATraceHeader hdr;

      hdr.magic    = FGA_TRACE_MAGIC;
      hdr.ntasks   = ntasks;
      hdr.tspace   = (uint32_t) TSPACE;
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
    fprintf(stderr,"%s: self-check: %d/%d valid trace vectors\n",
                    Prog_Name,nmatch,nchecked);
  }

  Free_Work_Data(work);

  if (ISTWO)
    Close_GDB(gdb2);
  Close_GDB(gdb1);

  oneFileClose(input);

  exit (0);
}
