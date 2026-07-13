# Lessons

## Verify the actual failure mode before assuming a cause (GPU -G integration)

**Pattern:** I spent a long time assuming the `-G` abort was a *seed coordinate-frame*
bug ("the tube-box centre maps to the wrong aseq/bseq position"). A 20-line `SEED_DEBUG`
print settled it in one run: `sa=(amid+dg)/2, sb=(amid-dg)/2` lands INSIDE FastGA's true
alignment on 25/25 tubes. The mapping was correct all along.

**Rule:** When a call aborts, first instrument the *inputs vs. the reference output*
empirically. Don't theorize about a coordinate frame for multiple edits when one print
resolves it.

## Read what a library function CONSUMES, not just what it produces

**Pattern:** I used `Compute_Trace_PTS` to "regenerate the trace from endpoints." Its
`align.h` doc says "Given ... trace point vector in `align->path.trace`" — it *refines an
existing* trace; `path->trace`/`path->tlen` are INPUTS. Feeding it only endpoints read
stale garbage → `exit(1)`, then SIGSEGV once I made it non-fatal. The correct
from-endpoints function is `Compute_Alignment(align, work, DIFF_TRACE, TSPACE)`
(aligns the substrings via `split_nd`, sets `path->trace/tlen/diffs`).

**Rule:** Before wiring a trace/alignment helper, read its align.h contract for required
INPUT fields, not only its outputs. `Compute_Trace_*` = refine existing trace;
`Compute_Alignment` = build from endpoints.

## gene_core EXIT is process-fatal unless Error_Buffer is non-NULL

`EXIT(x)` in gene_core.h: `if (Error_Buffer==NULL) exit(1); else return(x);`. To make an
aligner call non-fatal (so you can catch its error and fall back), set the global
`Error_Buffer` to a persistent non-NULL buffer first (see `ONEaln.c` for the
save/restore pattern). The message contents race across threads but are never read — only
the non-NULL-ness matters for control flow.

## The GPU speedup is gated by FastGA's sequential tube-walk, not the kernel

The GPU trace kernel is validated (100% valid trace-points) and end-to-end -G produces a
99%-coverage .1aln. But every -G variant is far SLOWER than CPU because the GPU is launched
ONE TUBE AT A TIME. The reason it can't be batched trivially: FastGA's tube walk has a
sequential dependency -- within a chain, each alignment's end anti-diagonal `eant` sets the
next tube's start (`alow = eant`, FastGA.c ~3494). Tube i+1 needs tube i's result. Only
first-tubes-of-chains and different contig-pairs are independent. So batching needs a
`search_seeds` REDESIGN (collect all chains' first tubes across a part -> one big
discover+trace batch -> iterate continuations), not a local edit. Lesson: before promising a
GPU speedup for a drop-in, check whether the host loop's granularity is even batchable --
a sequential data dependency in the caller can make a fast kernel useless.

## Correctness-first GPU offload ≠ speedup

`-G` with GPU endpoint discovery + CPU `Compute_Alignment` trace recovers 99% of
alignments but is ~10× SLOWER than CPU, because every aligning tube is aligned twice
(GPU x-drop for endpoints, then CPU `split_nd` for the trace). Real speedup requires the
GPU to emit the TSPACE trace-points directly. Prove the pipeline first, then move the
trace onto the GPU.
