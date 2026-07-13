# G4: wiring the GPU aligner into FastGA (design + status)

Goal: `FastGA -G` offloads the Phase-3 per-tube alignment to the A100 and produces a
`.1aln` of equivalent quality, faster end-to-end than the 32-thread CPU.

## Offload boundary (identified)
`align_contigs()` (FastGA.c ~3063), called per contig-pair from `search_seeds`. Its
`CALL_ALIGNER` block (~3323-3390) walks anti-diagonal tubes and per tube calls:
```
Local_Alignment(align, work, spec, dgmin, dgmax, amid, lbord, hbord)
  -> path->{abpos,aepos,bbpos,bepos,diffs} + trace-points (path.trace,path.tlen)
```
then filters (`rlen>=alnMin && alnRate*rlen>=diffs`) and writes `ovl`+trace to the
per-thread `.las`. `Local_Alignment` is the 80-98% bottleneck (M1).

## The GPU library (built, `fastga_gpu.{h,cu}` — compiles to `fastga_gpu.o`)
- `gpu_open/close` — persistent device context.
- `gpu_load_seqs(g, A, alen, B, blen)` — copy the contig-pair's NUMERIC (0-3) sequences
  resident on the GPU (reused across all its tubes).
- `gpu_discover_batch(g, n, sa, sb, ab,ae,bb,be,diffs)` — N tubes' seed anchors
  (contig coords) -> endpoints+diffs via the warp/thread x-drop wave (validated:
  84% of endpoints within 50bp of FastGA, self-test exact).

## Wiring plan (the remaining work)
1. **Seed anchor per tube**: from `(dgmin,dgmax,amid)`, seed diagonal `dg=(dgmin+dgmax)/2`,
   `sa=(amid+dg)/2`, `sb=(amid-dg)/2` (contig coords). Feed these as the batch.
2. **Batch, don't call per-tube**: restructure the `CALL_ALIGNER` do-loop into two passes —
   pass A collects all tube seed anchors for the contig-pair (or, better, buffer across the
   whole `search_seeds` part for GPU saturation); pass B calls `gpu_discover_batch` once,
   then iterates results.
3. **Sequences to GPU**: `align->aseq/bseq` are the (complemented if COMP) contig arrays;
   `gpu_load_seqs` once per contig-pair (they can be large — chr1 ~248 MB fits in 80 GB;
   reused across the pair's tubes).
4. **Trace-points**: GPU returns endpoints+diffs only. Regenerate trace-points on the CPU
   with `Compute_Trace_PTS(align, work, TSPACE, GREEDIEST, 1, -1)` over the returned
   rectangle (cheaper than the discovery wave), then `Compress_TraceTo8` + write `.las`.
   (Future: emit trace-points on the GPU to remove this CPU pass.)
5. **Flag**: `-G` selects the GPU path; without it, the CPU path is unchanged.
6. **Output** then flows through the existing redundancy removal + `la_sort`/`la_merge`
   unchanged -> `.1aln`.

## Validation plan
- Correctness/quality: compare the `-G` `.1aln` to the CPU `.1aln` — #alignments, bases
  covered, and per-alignment endpoint/diff agreement (accept non-bit-exact; require no
  coverage loss). The GPU wave finds equal-or-better edit distance (99.99% <= FastGA).
- Speed: end-to-end wall vs CPU 32-thread; and the sort+align phase specifically.

## Open risks (honest)
- **Endpoint fidelity 84% -> ~99%**: the x-drop uses `score = matches - 2*diffs`; FastGA
  uses its exact `TRIM_MLAG`/`M[]` criterion. Tune to match before the `.1aln` is trusted;
  the ~16% off are tandem-repeat/low-complexity tubes.
- **CPU trace regeneration cost**: `Compute_Trace_PTS` per alignment adds CPU work that
  eats into the GPU savings; measure, and move trace to GPU if it dominates.
- **Batch granularity**: per-contig-pair batches may be small (poor GPU use); buffering
  across a whole `search_seeds` part gives better saturation but more restructure.
- **Giant contigs**: loading chr-scale contigs per pair is fine on 80 GB but the H2D cost
  must be amortized across the pair's tubes.

## Status
Foundation done: offload boundary identified, GPU library built + linkable, kernels
validated (16.7x editdist; 84% discovery endpoints). Remaining: the `align_contigs`
two-pass batching + `-G` flag + CPU trace regen + `.1aln` quality/speed validation.
