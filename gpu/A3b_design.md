# A3b design: batching FastGA's align sweep for the GPU

## Why this is a rewrite, not a wrapper
`align_contigs` is a sequential sweep with **three nested levels of data-dependent state**,
each feeding the next tube's coordinates:

1. **Within a chain** (`CALL_ALIGNER` do-while, FastGA.c ~3388-3498): after each alignment,
   `eant = aepos+bepos` (non-comp) and `alow = (eant<=alow ? amid : eant)`; the next tube is
   `amid = alow + BUCK_ANTI`. Tube *i+1* cannot be placed without tube *i*'s endpoints.
2. **Within a diagonal-bucket triple** (`alast`, reset to -1 at ~3243): each completed chain
   sets `alast = alow`; the next chain runs only `if (ahgh > alast)` and clamps `alow` up to
   `alast`. Chain *i+1* depends on chain *i*'s result.
3. **Across triples and across contig-pairs:** independent — this is the only parallelism.

So a fast GPU batch can only be filled from *different triples / contig-pairs*, never from
the sequential interior of one. Batching therefore means **interleaving many triple-sweeps as
coroutines** and gathering their currently-pending tube from each into one batch.

## Architecture (coroutine batching)
Represent each triple-sweep as a resumable state machine `TripleJob { seeds cursor, alast,
current chain {dgmin,dgmax,alow,ahgh}, ovl/ctg context }`.

```
pool = all TripleJobs for a seed-sort part (across all contig-pairs in the part)
while pool not empty:
  batch = []
  for job in pool:                       # each job contributes its ONE pending tube
    tube = job.next_pending_tube()       # advances seeds/chains until a tube is ready
    if job.done: remove from pool; continue
    batch.append((job, tube seed sa,sb, endpoints-to-fill))
  gpu_discover_batch(batch)              # thousands of tubes -> endpoints  (warp-parallel*)
  gpu_trace_batch(passing)               # thousands -> trace-points        (warp-parallel)
  for (job, tube, result) in batch:
    apply guard + filter; write .las on pass
    job.advance(eant -> alow; alast; next chain)   # resumes the sequential interior
```
Most tubes resolve a chain in one step, so the pool drains in a few rounds, each a large batch
that saturates the A100 (the 15.8×/306k-aln/s regime from the feasibility bench).

\*Requires **warp-parallelizing the discovery kernel** first (currently 1 thread/task, the
feasibility bottleneck) using the trace kernel's band-across-lanes pattern.

## Correctness plan
- The per-job state machine reproduces the exact CPU walk (same `amid`, `eant`, `alow`,
  `alast`, `dgmin` mutation, filter `rlen>=alnMin && alnRate*rlen>=diffs`), so the SET of
  written alignments matches the CPU to within the trace kernel's equally-optimal-path
  differences (endpoints+diffs identical -> 99% coverage as in A3a).
- Validate against the CPU `.1aln` (records, aligned bases, identity) exactly as A3a was.
- `comp` and `self` tubes stay on the CPU (as today) or get their own jobs.

## Risks / open questions
- **Orchestration overhead:** the coroutine bookkeeping is CPU work that eats into the GPU
  win; must stay well under the ~14 s it replaces.
- **Load balance across jobs:** triples vary in chain count; a few heavy jobs can leave the
  GPU under-filled in late rounds (tail). Mitigate by pooling across the *whole part* (many
  contig-pairs), not one pair.
- **Memory:** trace wave scratch (`TR_CHUNK`×4 MB) + resident contigs; 32 GB gives 306k
  aln/s and fits the 80 GB A100.
- **Giant contigs:** H2D per contig-pair amortized across its tubes (already handled).

## Effort
Multi-day: the resumable triple-sweep state machine is the bulk (faithfully reproducing the
three-level sequential logic), plus warp-parallel discovery, plus validation. This is a
core-loop rewrite and must be validated bit-for-coverage against the CPU `.1aln` before trust.
