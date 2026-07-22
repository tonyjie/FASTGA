/*******************************************************************************************
 *
 *  extract_seeds -- Extract the FULL real distribution of per-alignment seeds from a FastGA
 *                    .1aln file into a compact ".seeds" file (see wave_harness.h), for the
 *                    GPU wave-alignment characterization study (Task 2).
 *
 *  Unlike extract_disc/extract_tasks/extract_deep, this tool dumps NO sequence windows --
 *  genome residency (loading the sequence data the seeds refer to) is handled separately at
 *  consumption time (Tasks 3/4/6).  It also applies NO dmin/wmax filtering: every alignment
 *  record in the .1aln produces exactly one SeedRec, giving the unfiltered distribution.
 *
 *  Reading/setup logic (opening the .1aln + both source GDBs) is copied verbatim from
 *  gpu/extract_deep.c (itself based on extract_tasks.c / ALNtoPAF.c), minus everything that
 *  fetches sequence data, since none is needed here.
 *
 *  Seed convention (contig-relative, verified against deep_cpu_bench):
 *      sa = (path->abpos + path->aepos) / 2
 *      sb = (path->bbpos + path->bepos) / 2
 *      seed_diag = sa - sb
 *      seed_anti = sa + sb
 *
 *  Usage: extract_seeds <in.1aln> <out.seeds>
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "GDB.h"
#include "align.h"
#include "alncode.h"
#include "wave_harness.h"

static char *Usage = "<alignment:path>[.1aln] <out:path>";

int main(int argc, char *argv[])
{ GDB       _gdb1, *gdb1 = &_gdb1;
  GDB       _gdb2, *gdb2 = &_gdb2;
  FILE     **units1;
  FILE     **units2;
  OneFile   *input;
  int64      novl;
  int        TSPACE;
  int        ISTWO;
  FILE      *out;

  Prog_Name = Strdup("extract_seeds","Allocating program name");

  if (argc != 3)
    { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
      exit (1);
    }

  //  Open the .1aln and the two source GDBs, exactly as extract_deep.c's main does
  //  (minus the header truncation and sequence-buffer setup, which we don't need since
  //  we never fetch sequence data -- only contig lengths from gdb->contigs[].clen).

  { char *pwd, *root, *cpath;
    char *src1_name, *src2_name;

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

    gdb1->seqs = units1[0];
    gdb2->seqs = units2[0];
  }

  out = fopen(argv[2],"wb");
  if (out == NULL)
    { fprintf(stderr,"%s: Cannot open %s for writing\n",Prog_Name,argv[2]);
      exit (1);
    }

  //  Reserve space for the header, to be filled in at the end.

  { WaveSeedHeader hdr;
    memset(&hdr,0,sizeof(hdr));
    if (fwrite(&hdr,sizeof(hdr),1,out) != 1)
      { fprintf(stderr,"%s: Cannot write header to %s\n",Prog_Name,argv[2]);
        exit (1);
      }
  }

  //  Sweep every alignment record: one SeedRec per overlap, no filtering, no sequences.

  { Overlap  _ovl, *ovl = &_ovl;
    Path    *path;
    int64    i;
    uint32_t nseeds;

    path = &(ovl->path);
    nseeds = 0;

    for (i = 0; i < novl; i++)
      { int      acontig, bcontig;
        int      alen, blen;
        int      ab, ae, bb, be;
        int      sa, sb;
        SeedRec  rec;

        Read_Aln_Overlap(input,ovl);
        Skip_Aln_Trace(input);

        acontig = ovl->aread;
        bcontig = ovl->bread;
        alen = gdb1->contigs[acontig].clen;
        blen = gdb2->contigs[bcontig].clen;

        ab = path->abpos;  ae = path->aepos;
        bb = path->bbpos;  be = path->bepos;

        sa = (ab+ae)/2;
        sb = (bb+be)/2;

        rec.aread     = acontig;
        rec.bread     = bcontig;
        rec.flags     = (int32_t) ovl->flags;
        rec.alen      = alen;
        rec.blen      = blen;
        rec.seed_anti = sa+sb;
        rec.seed_diag = sa-sb;
        rec.ref_ab    = ab;
        rec.ref_ae    = ae;
        rec.ref_bb    = bb;
        rec.ref_be    = be;
        rec.ref_diffs = path->diffs;

        if (fwrite(&rec,sizeof(rec),1,out) != 1)
          { fprintf(stderr,"%s: Write error to %s\n",Prog_Name,argv[2]);
            exit (1);
          }

        nseeds++;
      }

    //  Patch in the final header.

    { WaveSeedHeader hdr;

      hdr.magic    = WAVE_SEEDS_MAGIC;
      hdr.nseeds   = nseeds;
      hdr.tspace   = (uint32_t) TSPACE;
      hdr.reserved = 0;

      rewind(out);
      if (fwrite(&hdr,sizeof(hdr),1,out) != 1)
        { fprintf(stderr,"%s: Cannot patch header in %s\n",Prog_Name,argv[2]);
          exit (1);
        }
    }

    fclose(out);

    fprintf(stderr,"%s: %u alignments in file, %u seeds written (tspace=%d)\n",
                    Prog_Name,(unsigned) novl,nseeds,TSPACE);
  }

  if (ISTWO)
    Close_GDB(gdb2);
  Close_GDB(gdb1);

  oneFileClose(input);

  exit (0);
}
