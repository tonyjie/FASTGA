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

## Status (2026-07-13) — END-TO-END WORKING, correctness-first

**`-G` produces a valid `.1aln` recovering 99% of alignments.** On EXAMPLE (hap1 vs
hap2, T=8):

| | records | aligned A-bases |
|---|---:|---:|
| CPU baseline | 323,569 | 632,119,471 |
| GPU (`-G`)   | 320,433 (**99.03%**) | 626,012,667 (**99.03%**) |

The full Phase-3 offload runs to completion, the output flows through the normal
redundancy-removal + `la_sort`/`la_merge`, and the `.1aln` is near-equivalent to the CPU
one. The ~1% gap is the tandem-repeat / low-complexity tubes where the GPU x-drop
endpoints differ slightly.

### What the two earlier blocking bugs actually were

1. **The seed mapping was never wrong.** A per-tube diagnostic (`SEED_DEBUG`) proved
   `sa=(amid+dg)/2, sb=(amid-dg)/2` lands *inside* FastGA's true alignment on 25/25
   tubes (`sa_in=sb_in=1`). The abort was not a coordinate-frame error.
2. **`Compute_Trace_PTS` was the wrong function.** It *refines an existing* trace-point
   vector (`path->trace`/`path->tlen` are INPUTS); `gpu_align_tube` never set them, so it
   read stale garbage → first a hard `exit(1)` ("Bad alignment between trace points"),
   then, once made non-fatal, a SIGSEGV. The function that builds a trace from *just
   endpoints* is **`Compute_Alignment(align, work, DIFF_TRACE, TSPACE)`** (it aligns the
   substrings via `split_nd` and sets `path->trace/tlen/diffs`).

### The fix (current code)

- `main`: under `-G`, set gene_core's `Error_Buffer` non-NULL so any internal
  `EXIT` in the aligner *returns* instead of `exit()`-ing the process (caught below).
- `gpu_align_tube`: GPU discovers endpoints → self-consistency guard (monotone, in-range,
  keep-filter, drift ≤ 512) → `Compute_Alignment(DIFF_TRACE)` builds the trace →
  **on any failure or a rejected guard, fall back to the exact CPU `Local_Alignment`**
  for that tube (bit-identical to no `-G`). No tube can abort the run or be silently lost.

### The catch: correctness-first is SLOW (no speedup yet)

EXAMPLE T=8, sort+align phase wall time:

| | seed-merge | **sort+align** | total |
|---|---:|---:|---:|
| CPU baseline | 41 s | **16.7 s** | 57.5 s |
| GPU (`-G`)   | 41 s | **44.5 min** (2671 s) | 45.2 min |

The align phase is **~160× slower**, not 2×: `Compute_Alignment` (`split_nd`, full
divide-and-conquer trace of each substring) is itself far more expensive than
`Local_Alignment`'s banded x-drop wave, *and* every aligning tube is now aligned twice
(GPU x-drop for endpoints + CPU `split_nd` for the trace), *and* the GPU runs one tube
per kernel launch (no batching). The GPU endpoint discovery is not on the critical path;
the CPU trace regeneration is. This version is a **correctness proof of the pipeline**,
not a speedup.

## Validation B: GPU discovery fidelity vs divergence (2026-07-13)

`gpu/disc_fidelity.c` (`make disc_fidelity`) replays real FastGA discovery tasks
(`extract_disc`) through the GPU x-drop and buckets the endpoint agreement by **per-task
local divergence** (`ref_diffs / ref-A-length`) — turning one dataset into a
fidelity-vs-divergence curve. EXAMPLE, 20,000 tasks:

| local divergence | tasks | ≤50 bp | ≤10 bp |
|---|---:|---:|---:|
| 0–1%   |    123 | **100.0%** | 93.5% |
| 1–2%   |    264 | **100.0%** | 94.7% |
| 2–5%   |  1,528 | **97.4%**  | 89.4% |
| 5–10%  |  5,380 | **95.1%**  | 80.4% |
| 10–20% | 10,761 | 81.4%      | 64.5% |
| >20%   |  1,944 | 52.5%      | 36.5% |
| **ALL**| 20,000 | **83.8%**  | 68.5% |

**Reading:** the GPU x-drop is near-exact (≥95% within 50 bp) up to ~10% local
divergence — i.e. on genuine homology, where FastGA alignments actually live — and only
falls off in the high-divergence tail (tandem-repeat / low-complexity / spurious tubes).
That tail is precisely where `gpu_align_tube`'s guard rejects the GPU endpoints and the
exact CPU `Local_Alignment` takes over, which is why the end-to-end `.1aln` keeps 99% of
alignments/bases despite the 83.8% overall endpoint match. The bucketing is a proxy:
divergence here is driven by within-dataset repeats rather than cross-species orthology,
but the x-drop responds to edit rate regardless of provenance, so the *relationship*
(fidelity vs divergence) transfers. A true high-divergence genome (e.g. mouse ~15%) would
populate the tail with orthology instead of repeats; measuring it needs a GIX build
(deferred — not part of the fast proxy).

### To get real speedup (next phase)

1. **Emit trace-points on the GPU.** The only way to remove the double-alignment: have the
   kernel produce the TSPACE=100 trace-point vector (diffs per panel) directly, so the CPU
   just packs it (`Compress_TraceTo8`) — no `Compute_Alignment`. This is the real win.
2. **Batch across tubes / contig-pairs** for GPU saturation (currently one tube per
   `gpu_discover_batch` call).
3. Improve x-drop fidelity toward FastGA's `TRIM_MLAG` (84% → ~99%) to shrink the CPU
   fallback rate.
4. Validate on divergent genomes (human/chimp/mouse) and measure end-to-end vs CPU T=32.
