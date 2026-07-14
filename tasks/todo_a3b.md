# A3b: batched GPU align pipeline — implementation plan

**Goal:** a *measured* end-to-end GPU wall speedup on human, by batching FastGA's per-tube
wave (discovery+trace) onto the A100 while keeping the `.1aln` correct (≥99% coverage vs CPU).

**Architecture (chosen): dynamic batching via a shared queue + whole-genome-resident sequences.**
Rationale over the coroutine rewrite: each CPU thread keeps its EXACT sequential tube-walk
(`alow=eant`, `alast`, triple logic unchanged → the alignment *selection* is provably identical);
only the per-tube GPU *call* is batched. The whole genome stays 2-bit-resident on the GPU
(kills the 29% decompression + all per-pair transfer). Correctness surface is small and each
task validates against the CPU `.1aln` by md5 of the ALNshow summary.

**Global constraints:**
- Every task ends with: build clean, run on EXAMPLE, `ALNshow <out>.1aln | tail -n+3 | md5sum`
  equals the CPU baseline's (or, once GPU endpoints diverge, record count ≥99% + coverage check).
- Keep it all under `#ifdef GPU` / a `-G` flag; base `FastGA` untouched.
- No busy-wait; blocked threads must not burn CPU (condition variables).

---

## Task 1 — Whole-genome-resident sequences on the GPU (foundation)
**Files:** `gpu/fastga_gpu.{h,cu}`, a new offline test.
- Add `gpu_load_genome_2bit(g, packedA, baseA[], nA, packedB, ...)`: upload each genome's
  `.bps` once, unpack on-device to a resident NUMERIC genome (concatenated contigs), and keep a
  device table of per-contig base offsets. (Reuse `unpack2bit`.)
- Make `gpu_discover_batch` / `gpu_trace_batch` accept GENOME-relative coords (caller adds the
  contig base offset). No kernel change — they already index `A[coord]`.
- **Validate:** extract EXAMPLE tasks, offset every task's coords by its contig base, load the
  whole genome resident, run discover+trace → identical endpoints/traces to the per-contig runs.
- **Gate:** offline fidelity unchanged (83.8% disc, 100% trace-valid).

## Task 2 — The batching queue (single-thread plumbing)
**Files:** `FastGA.c` (`#ifdef GPU`), maybe `gpu/batch_queue.{h,c}`.
- A shared queue: `submit(tube request) -> blocks -> returns {endpoints, trace}`. A flush runs
  `gpu_discover_batch`+`gpu_trace_batch` on all pending, wakes submitters.
- Flush policy: when pending ≥ FLUSH_N **or** all active threads are blocked (deadlock guard).
- Wire `gpu_align_tube` to `submit()` instead of the per-tube n=1 calls.
- **Validate at T=1, FLUSH_N=1:** `-G` on EXAMPLE → ALNshow md5 == the current `-G` output
  (this proves the queue plumbing changes nothing).
- **Gate:** EXAMPLE `.1aln` record count within the A3a 99% band.

## Task 3 — Scale threads + adaptive flush for saturation
- Raise the worker count (oversubscribe; blocked threads are cheap) so batches reach ~512-2048.
- Adaptive flush: size + a short timeout; the "all blocked" guard prevents stalls.
- **Validate:** EXAMPLE + human `.1aln` coverage ≥99% vs CPU; no deadlock/hang across 3 runs.
- **Gate:** GPU utilization > 50% during the align phase (nvidia-smi), correctness held.

## Task 4 — End-to-end measurement
- Human GRCh38×CHM13, T=32-equivalent: wall vs CPU stock; `.1aln` coverage/identity vs CPU.
- Long-alignment tail: alignments exceeding the kernel caps (short-int x >32 kbp, band, depth)
  fall back to CPU `Local_Alignment` — count them; raise caps (int `x`) if they dominate.
- **Deliver:** the measured GPU wall speedup + a coverage table. Honest about the fallback tail
  and the ~2-3× realistic ceiling (wave is ~50% of wall; non-wave stays CPU).

---

## Risks / notes
- Memory: resident NUMERIC genome = ~3 GB each (6 GB); or keep 2-bit (750 MB each) and unpack
  per-access in-kernel to save 4×. Start with unpacked NUMERIC (simpler), optimize later.
- Deadlock: the "all-active-threads-blocked → force flush" guard is mandatory and must be tested.
- Batch size vs threads: batch ≤ live thread count; need ~512+ for good util (see feasibility
  table: 512→40k, 2048→114k aln/s).
- Comp strand: A is complemented per contig-pair; the resident genome is forward, so comp tubes
  need the complemented A — either a second resident (complemented A) or per-request handling.
  Decide in Task 2; keep comp on the CPU fallback initially if it complicates correctness.
