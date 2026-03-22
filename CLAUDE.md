## Workflow Orchestration

### 1. Plan Node Default
- Enter plan mode for ANY non-trivial task (3+ steps or architectural decisions)
- If something goes sideways, STOP and re-plan immediately – don't keep pushing
- Use plan mode for verification steps, not just building
- Write detailed specs upfront to reduce ambiguity

### 2. Subagent Strategy
- Use subagents liberally to keep main context window clean
- Offload research, exploration, and parallel analysis to subagents
- For complex problems, throw more compute at it via subagents
- One tack per subagent for focused execution

### 3. Self-Improvement Loop
- After ANY correction from the user: update `tasks/lessons.md` with the pattern
- Write rules for yourself that prevent the same mistake
- Ruthlessly iterate on these lessons until mistake rate drops
- Review lessons at session start for relevant project

### 4. Verification Before Done
- Never mark a task complete without proving it works
- Diff behavior between main and your changes when relevant
- Ask yourself: "Would a staff engineer approve this?"
- Run tests, check logs, demonstrate correctness

### 5. Demand Elegance (Balanced)
- For non-trivial changes: pause and ask "is there a more elegant way?"
- If a fix feels hacky: "Knowing everything I know now, implement the elegant solution"
- Skip this for simple, obvious fixes – don't over-engineer
- Challenge your own work before presenting it

### 6. Autonomous Bug Fixing
- When given a bug report: just fix it. Don't ask for hand-holding
- Point at logs, errors, failing tests – then resolve them
- Zero context switching required from the user
- Go fix failing CI tests without being told how

## Task Management

1. **Plan First**: Write plan to `tasks/todo.md` with checkable items
2. **Verify Plan**: Check in before starting implementation
3. **Track Progress**: Mark items complete as you go
4. **Explain Changes**: High-level summary at each step
5. **Document Results**: Add review section to `tasks/todo.md`
6. **Capture Lessons**: Update `tasks/lessons.md` after corrections

## Core Principles

- **Simplicity First**: Make every change as simple as possible. Impact minimal code.
- **No Laziness**: Find root causes. No temporary fixes. Senior developer standards.
- **Minimal Impact**: Changes should only touch what's necessary. Avoid introducing bugs.

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## Project Overview

**FastGA** is a fast whole-genome aligner by Gene Myers, Richard Durbin, and Chenxi Zhou. It finds all local DNA alignments between two high-quality genome assemblies. The core design principle is **cache coherence**: instead of hash-table lookups (which cause random memory accesses / cache misses), FastGA uses MSD radix sorts and linear merges throughout.

Reference paper: `paper/Myers and Durbin - FastGA Fast Genome Alignment.pdf`

### Repository Structure

This is a pure C codebase. All source files are in the root directory — no subdirectories for source code.

| Category | Key Files | Purpose |
|---|---|---|
| Core aligner | `FastGA.c` | Main pipeline: seed merge, sort, chain, align, output |
| Index builder | `GIXmake.c`, `MSDsort.c` | Build sorted k-mer index (GIX) from genome database |
| Genome database | `GDB.c`, `GDB.h`, `FAtoGDB.c`, `GDBtoFA.c` | Convert FASTA to 2-bit compressed genome database |
| K-mer table I/O | `libfastk.c`, `libfastk.h` | Stream/random access to partitioned k-mer tables |
| Alignment engine | `align.c`, `align.h` | Wave-based local alignment (LA-finder) |
| Alignment encoding | `alncode.c`, `alncode.h` | Trace-point encoding for compact alignment storage |
| Radix sort | `RSDsort.c` (reverse MSD), `MSDsort.c` (forward MSD) | In-place radix sort for seeds and k-mers |
| Alignment utilities | `ALNshow.c`, `ALNtoPAF.c`, `ALNtoPSL.c`, `ALNchain.c`, `ALNplot.c`, `ALNreset.c` | View, convert, filter alignments |
| Format converters | `PAFtoALN.c`, `PAFtoPSL.c` | Convert between PAF/PSL/ALN formats |
| Annotation | `ANO.c`, `ANO.h`, `ANOshow.c`, `ANOstat.c`, `ANOtoBED.c`, `BEDtoANO.c` | Soft masking via .1ano interval files |
| Index utilities | `GIXshow.c`, `GIXrm.c`, `GIXxfer.c` | Inspect, remove, copy/move GIX ensembles |
| ONEcode library | `ONElib.c`, `ONElib.h` | Binary file format used for .1gdb, .1aln, .1ano |
| Shared utilities | `gene_core.c`, `gene_core.h`, `hash.c`, `hash.h`, `select.c`, `select.h` | Memory allocation wrappers, hashing, selection |

## Build & Run

### Build
```bash
make          # builds all binaries
make clean    # remove binaries
make install  # copies binaries to ~/bin (set DEST_DIR in Makefile to change)
```
Requires: `gcc`, `zlib` (`-lz`), `pthreads` (`-lpthread`), `math` (`-lm`). No external dependencies.

### Basic Usage
```bash
# Simplest: compare two FASTA genomes, outputs PAF to stdout
FastGA genome1.fa genome2.fa

# With options
FastGA -v -T16 -k genome1 genome2          # verbose, 16 threads, keep intermediate files
FastGA -v -T8 -1:output.1aln genome1 genome2  # output as compact .1aln file
FastGA -v -T8 -pafx genome1 genome2         # PAF with CIGAR strings

# Self-comparison (finds repeats, haplotype homology)
FastGA genome1

# Step-by-step (manual pipeline)
FAtoGDB genome1.fa                          # FASTA -> .1gdb
FAtoGDB genome2.fa
GIXmake -v -T8 genome1                     # .1gdb -> .gix + .ktab + .post
GIXmake -v -T8 genome2
FastGA -v -T8 -k genome1.gix genome2.gix   # align using pre-built indices
ALNshow genome1 genome2 output.1aln         # view alignments
ALNtoPAF genome1 genome2 output.1aln        # convert to PAF
```

### Key Options
| Flag | Default | Description |
|---|---|---|
| `-T<int>` | 8 | Number of threads |
| `-P<dir>` | `$TMPDIR` or `.` | Directory for temporary files |
| `-k` | off | Keep intermediate GDB/GIX files after completion |
| `-S` | off | Symmetric mode (use adaptamers from both genomes) |
| `-M` | off | Use soft masks encoded in GIX |
| `-f<int>` | 10 | Max k-mer frequency threshold (adaptamer repetitiveness filter) |
| `-c<int>` | 85 | Minimum chain coverage (bases) |
| `-s<int>` | 1000 | Maximum seed spacing in a chain |
| `-l<int>` | 100 | Minimum alignment length |
| `-i<float>` | 0.7 | Minimum alignment identity |
| `-L:<path>` | none | Append performance log to file |

## Algorithm Pipeline

The pipeline has 4 major phases inside `FastGA.c`:

### Phase 1: Adaptive Seed Merge
- Merges two sorted GIX indices in a single linear pass
- For each k-mer (K=40, syncmer-filtered to ~50% of positions) in genome A, finds the longest prefix match in genome B's index — this is the "adaptamer"
- If the k-mer occurs <= `-f` (default 10) times in genome B, emits seed position pairs
- Seeds are written to temporary files partitioned by A-contig and strand (normal/complement)
- This is cache-coherent: both indices are swept linearly, no random lookups

### Phase 2: Reimport and Radix Sort
- Reads seed pairs back from temp files
- Transforms to (anti-diagonal, diagonal-bucket, j-contig) representation
- Sorts in-place using reverse MSD radix sort (`rmsd_sort()` in `RSDsort.c`)
- Groups seeds by contig pair and diagonal band for chaining

### Phase 3: Chaining and Alignment
- Linear scan over sorted seeds to build chains: seeds within diagonal band width 128 and spacing <= `-s` (1000)
- Chains covering >= `-c` (85) bases trigger alignment
- Calls `Local_Alignment()` (wave-based LA-finder from `align.c`) on each chain's tube
- Removes redundant/entwined alignments
- Writes results to per-thread `.las` temp files

### Phase 4: Sort and Merge Output
- Sorts each thread's overlap file by (aread, abpos, bread, comp)
- K-way merge of all threads into final `.1aln` output (ONEcode format)
- If PAF/PSL output requested, calls `ALNtoPAF`/`ALNtoPSL` via `system()`

## Data Structures and File Formats

### GDB (.1gdb)
- ONEcode binary + hidden `.bps` file with 2-bit compressed DNA
- Compact: ~genome_size/4 bytes for sequence data
- Supports random contig access without text parsing

### GIX (.gix + hidden files)
- Truncated suffix array of syncmer-filtered K-mers (K=40)
- Each entry: 10 bytes (2-bit packed 40-mer) + 2 bytes (contig #) + 4 bytes (signed position)
- Plus LCP (longest common prefix) array computed during MSD radix sort
- Hidden partition files: `.<root>.ktab.<1..NPARTS>` (k-mer table) + `.<root>.post.<1..NPARTS>` (position lists)
- **Size: ~14 GB per 1 Gbp of genome** (README figure; paper says ~11 GB/Gbp after syncmer filtering)
- The `.gix` stub file itself contains a 128MB prefix index (16M int64 entries)

### ALN (.1aln)
- ONEcode binary with trace-point encoding (delta=100bp panels)
- Very compact: e.g., 66 MB for 1.63 Gbp of aligned sequence (vs 1.03 GB in PAF with CIGAR)
- ~15-16x smaller than equivalent PAF output

## Storage Problem: Large Intermediate Files

This is the primary motivation for the `optimize-memory` branch.

### The Problem
When aligning large genomes, FastGA creates massive intermediate files that can consume hundreds of GBs of disk space. Running multiple FastGA processes on the same node compounds this, easily hitting multi-TB scratch limits and causing crashes.

### Where the Storage Goes

**1. GIX Index Files (~14 GB per Gbp of genome)**
- The dominant persistent storage cost
- A 3 Gbp human genome requires ~42 GB of index files
- Two genomes = ~84 GB just for indices
- With `-k` flag these persist for reuse; without it they are rebuilt each time

**2. Seed Pair Temp Files (in `$TMPDIR` or `-P` dir)**
- Created during Phase 1 (adaptive seed merge)
- `_pair.<pid>.<k>.N` and `.C` files: NPARTS x NTHREADS files for normal and complement strands
- Store (lcp, post1, post2) tuples for every seed hit
- Size scales with number of k-mer matches between genomes — can be 10-100x genome size for repetitive genomes
- Files are `unlink()`-ed immediately after creation (temp files), but still consume disk while open

**3. Alignment Temp Files (in `$TMPDIR` or `-P` dir)**
- `_uniq.<pid>.<p>.las` and `_algn.<pid>.<p>.las` per thread
- Scale with number of alignments found
- Also immediately unlinked

**4. GIXmake Distribution Files**
- `._post.<pid>.<k>.idx`: NTHREADS x NPARTS temp files during index construction
- Immediately unlinked after use

### Key Scaling Factors
- GIX size scales linearly with genome size
- Seed pair files scale with genome_size x average_k-mer_frequency (can be superlinear for repetitive genomes)
- NPARTS is auto-calculated to keep sort partitions near 4GB; larger genomes = more partitions = more files
- Memory during index construction: peak is the MSD radix sort array (~4GB per partition)
- Memory during alignment: sort array sized to max seeds in any partition panel, plus ~2GB I/O buffers

### Optimization Goals
1. **Reduce intermediate file sizes** — especially GIX indices and seed pair temp files
2. **Maintain alignment precision** — no significant loss in sensitivity or specificity
3. **Maintain runtime performance** — avoid major slowdowns
4. **Enable concurrent runs** — multiple FastGA processes on the same node without disk exhaustion
