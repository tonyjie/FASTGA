/*******************************************************************************************
 *
 *  extract_tasks -- Extract aligned segment pairs from a FastGA .1aln file into a binary
 *                    "task file" for a GPU edit-distance experiment (see task_format.h).
 *
 *  Single-threaded.  Based directly on the sequence-extraction path of ALNtoPAF.c: for
 *  each alignment, the A segment is aln[abpos..aepos) of the (NUMERIC) A-contig, and the
 *  B segment is the (possibly reverse-complemented, for COMP overlaps) B-contig piece
 *  [bbpos..bepos), extracted exactly the way ALNtoPAF prepares aln->aseq / aln->bseq
 *  before it computes CIGAR/diffs.  FastGA's own path->diffs (read straight from the
 *  .1aln file) is stored as the reference edit distance.
 *
 *  Usage: extract_tasks <in.1aln> <out.tasks> [max_tasks=200000] [wmax=8192]
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
#include "task_format.h"

static char *Usage = "<alignment:path>[.1aln] <out:path> [max_tasks=200000] [wmax=8192]";

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

  Prog_Name = Strdup("extract_tasks","Allocating program name");

  if (argc < 3 || argc > 5)
    { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
      exit (1);
    }

  max_tasks = 200000;
  wmax      = 8192;
  if (argc >= 4)
    max_tasks = strtoll(argv[3],NULL,10);
  if (argc >= 5)
    wmax = (uint32_t) strtoul(argv[4],NULL,10);

  //  Open the .1aln and the two source GDBs, exactly as ALNtoPAF.c's main does when it
  //  needs sequence data (CIGAR || DIFFS path).

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

  { FGATaskHeader hdr;
    memset(&hdr,0,sizeof(hdr));
    if (fwrite(&hdr,sizeof(hdr),1,out) != 1)
      { fprintf(stderr,"%s: Cannot write header to %s\n",Prog_Name,argv[2]);
        exit (1);
      }
  }

  //  Sweep every alignment, extracting the A/B segments exactly as ALNtoPAF.c does for
  //  its CIGAR/diffs path (see ALNtoPAF.c lines ~256-276), then dump the binary task.

  { Overlap   _ovl, *ovl = &_ovl;
    Alignment _aln, *aln = &_aln;
    Path     *path;
    char     *aseq, *bseq;
    int       alast;
    uint32_t  ntasks, nskip;
    int64     i;

    aseq = New_Contig_Buffer(gdb1);
    bseq = New_Contig_Buffer(gdb2);
    if (aseq == NULL || bseq == NULL)
      exit (1);

    aln->path = path = &(ovl->path);
    aln->aseq = aseq;

    ntasks = 0;
    nskip  = 0;
    alast  = -1;

    for (i = 0; i < novl; i++)
      { int   acontig, bcontig;
        int   bmin, bmax;
        char *bact;
        int   m, n;

        Read_Aln_Overlap(input,ovl);
        Skip_Aln_Trace(input);

        if (ntasks >= (uint32_t) max_tasks)
          continue;             //  keep sweeping only to count remaining skips accurately? no:
                                 //  we stop extraction but must still advance the OneFile reader,
                                 //  which Read_Aln_Overlap/Skip_Aln_Trace above already did.
        acontig = ovl->aread;
        bcontig = ovl->bread;
        aln->alen  = gdb1->contigs[acontig].clen;
        aln->blen  = gdb2->contigs[bcontig].clen;
        aln->flags = ovl->flags;

        m = path->aepos - path->abpos;
        n = path->bepos - path->bbpos;

        if ((uint32_t) m > wmax || (uint32_t) n > wmax)
          { nskip++;
            continue;
          }

        if (acontig != alast)
          Get_Contig(gdb1,acontig,NUMERIC,aseq);
        alast = acontig;

        if (COMP(aln->flags))
          { bmin = (aln->blen - path->bepos);
            bmax = (aln->blen - path->bbpos);
          }
        else
          { bmin = path->bbpos;
            bmax = path->bepos;
          }

        bact = Get_Contig_Piece(gdb2,bcontig,bmin,bmax,NUMERIC,bseq);
        if (COMP(aln->flags))
          { Complement_Seq(bact,bmax-bmin);
            aln->bseq = bact - (aln->blen-bmax);
          }
        else
          aln->bseq = bact - bmin;

        //  A = aln->aseq[abpos..aepos), B = aln->bseq[bbpos..bepos) -- NUMERIC (0..3),
        //  already oriented/reverse-complemented per FastGA's COMP convention.

        { uint32_t rec[3];
          unsigned char *A = (unsigned char *) (aln->aseq + path->abpos);
          unsigned char *B = (unsigned char *) (aln->bseq + path->bbpos);

          rec[0] = (uint32_t) m;
          rec[1] = (uint32_t) n;
          rec[2] = (uint32_t) path->diffs;

          if (fwrite(rec,sizeof(rec),1,out) != 1 ||
              fwrite(A,1,(size_t) m,out) != (size_t) m ||
              fwrite(B,1,(size_t) n,out) != (size_t) n)
            { fprintf(stderr,"%s: Write error to %s\n",Prog_Name,argv[2]);
              exit (1);
            }
        }

        ntasks++;
      }

    free(bseq-1);
    free(aseq-1);

    //  Patch in the final header.

    { FGATaskHeader hdr;

      hdr.magic    = FGA_TASK_MAGIC;
      hdr.ntasks   = ntasks;
      hdr.wmax     = wmax;
      hdr.reserved = 0;

      rewind(out);
      if (fwrite(&hdr,sizeof(hdr),1,out) != 1)
        { fprintf(stderr,"%s: Cannot patch header in %s\n",Prog_Name,argv[2]);
          exit (1);
        }
    }

    fclose(out);

    fprintf(stderr,"%s: %u alignments in file, %u tasks written, %u skipped (> wmax=%u),"
                    " stopped-early=%s\n",
                    Prog_Name,(unsigned) novl,ntasks,nskip,wmax,
                    (ntasks >= (uint32_t) max_tasks && novl > ntasks + nskip) ? "yes" : "no");
  }

  if (ISTWO)
    Close_GDB(gdb2);
  Close_GDB(gdb1);

  oneFileClose(input);

  exit (0);
}
