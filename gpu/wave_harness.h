/*******************************************************************************************
 *
 *  wave_harness.h -- On-disk format for the full-distribution seed extraction used by the
 *                     GPU wave-alignment characterization study (Task 2).  Produced by
 *                     gpu/extract_seeds.c, consumed by Tasks 3 and 6.
 *
 *  File layout:  WaveSeedHeader  followed by header.nseeds SeedRec's (no sequence data --
 *  genome residency is handled separately, at consumption time).
 *
 *******************************************************************************************/

#ifndef WAVE_HARNESS_H
#define WAVE_HARNESS_H

#include <stdint.h>

#define WAVE_SEEDS_MAGIC 0x53564157u   // "WAVS"

typedef struct
  { uint32_t magic, nseeds, tspace, reserved;
  } WaveSeedHeader;

typedef struct
  { int32_t aread, bread;            //  contig indices into gdb1 / gdb2
    int32_t flags;                   //  COMP bit = B reverse-complemented
    int32_t alen, blen;              //  contig lengths
    int32_t seed_anti, seed_diag;    //  midpoint seed: anti = sa+sb, diag = sa-sb
                                      //    where sa=(ref_ab+ref_ae)/2, sb=(ref_bb+ref_be)/2
    int32_t ref_ab, ref_ae, ref_bb, ref_be, ref_diffs;  //  .1aln reference, for validation
  } SeedRec;

/*******************************************************************************************
 *
 *  .cpuref format -- Task 3's CPU-baseline reference dump, produced by gpu/wave_bench_cpu.c,
 *  consumed by Tasks 6/7 to validate the GPU wave kernel.  One CpuRefRec per seed, in the
 *  SAME ORDER as the corresponding .seeds file (record i <-> SeedRec i), holding the
 *  contig-relative endpoints/diffs that FastGA's real Local_Alignment produced when started
 *  from that seed -- i.e. the actual algorithm output the GPU port must reproduce, as opposed
 *  to SeedRec's ref_* fields (the .1aln's own recorded truth).
 *
 *******************************************************************************************/

#define WAVE_CPUREF_MAGIC 0x55504357u   // "WCPU"

typedef struct
  { uint32_t magic, nrecs, reserved1, reserved2;
  } CpuRefHeader;

typedef struct
  { int32_t abpos, aepos, bbpos, bepos, diffs;   //  contig-relative, from Local_Alignment's path
  } CpuRefRec;

#endif
