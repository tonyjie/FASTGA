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

## Status (2026-07-13)

**Wired and building.** `FastGA.gpu` (make target; `-DGPU`, links `fastga_gpu.o` +
`-lcudart`, `--default-stream per-thread`) compiles; base `FastGA` unaffected (all
`#ifdef GPU`). The `-G` path is fully plumbed: flag parse, per-thread `gpu_ctx`
open/close, `gpu_load_seqs` once per contig-pair, `gpu_align_tube` replaces the
non-comp `Local_Alignment`, CPU `Compute_Trace_PTS` regenerates trace. It runs
end-to-end (EXAMPLE, 42 s, GPU engaged).

**Blocking bug (diagnosed): endpoint validity + no fallback.** The run aborts in
`Compute_Trace_PTS` ("Bad alignment between trace points") on the first bad GPU
alignment. Two compounding causes:
  (a) the seed formula `sa=(amid+dg)/2, sb=(amid-dg)/2` matches the aligner's
      `anti`/diagonal->index mapping, BUT it seeds from the tube-BOX centre, which can
      lie outside the true alignment (the box is >= the alignment) -> the x-drop wave
      finds a wrong/partial alignment;
  (b) `gpu_align_tube` feeds those endpoints straight into `Compute_Trace_PTS`, which
      hard-exits (not a return code) the moment a 100 bp panel doesn't align -> one bad
      tube kills the whole run. (The kernel itself is fine: 84% within 50 bp when seeded
      at the true midpoint via `extract_disc`.)

**To finish (bounded):**
1. **Seed from a real chain seed** (a sorted-seed entry inside the chain that triggered
   the tube) rather than the box centre -> guaranteed inside the alignment.
2. **Validate before trace / fall back**: check the GPU endpoints (monotone,
   `aepos<=alen`, `bepos<=blen`, `diffs` consistent with length) and for any that look
   off, fall back to CPU `Local_Alignment` for that tube instead of aborting. Or make a
   non-fatal `Compute_Trace` path.
3. Improve x-drop fidelity toward FastGA's `TRIM_MLAG` (84% -> ~99%).
4. Then batch (two-pass per contig-pair) for speed; validate `.1aln` coverage vs CPU +
   end-to-end timing.
