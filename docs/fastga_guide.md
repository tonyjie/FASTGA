# FastGA: A Practical Guide

This guide explains FastGA from the ground up — what it does, why it's designed the way it is, and what happens to your data at every step.

---

## What Problem Does FastGA Solve?

Imagine you have two genome assemblies — say, a human reference genome and a newly assembled genome from a different individual. Each is about 3 billion DNA bases long. You want to find **all the regions where these two genomes match**, even approximately.

This is **whole-genome alignment**: find every local region in genome A that aligns to some region in genome B. The output is a list of alignments, each saying "bases 1000–5000 in contig X of genome A match bases 2000–6000 in contig Y of genome B, with 95% identity."

The challenge: doing this fast on genomes with billions of bases.

---

## The Key Idea: Cache Coherence

Most genome aligners (like minimap2) use **hash tables** to find matching sequences. You take a short sequence (a "k-mer") from genome A, hash it, and look up where it appears in genome B. This is conceptually simple, but it has a fundamental hardware problem:

**Hash table lookups cause random memory accesses.**

When your data structure is 30+ GB (as it is for human genomes), almost every lookup goes to a different region of RAM. The CPU cache — a small, fast memory buffer — can't help because consecutive lookups hit completely unrelated memory addresses. Each lookup stalls the CPU for ~200 clock cycles waiting for data from main memory. When you're doing billions of lookups, this adds up.

FastGA takes a radically different approach:

> **Instead of hashing and looking up, sort everything and sweep linearly.**

Both genomes' k-mers are pre-sorted. To find matches, FastGA simply merges two sorted lists in one pass — like merging two sorted decks of cards. Every memory access is sequential and predictable. The CPU cache stays warm. The prefetcher can anticipate what's needed next.

This is the "cache coherence" principle, and it drives every design decision in FastGA.

---

## The Big Picture: What Happens When You Run FastGA

```
FastGA genome1.fa genome2.fa
```

When you type this, FastGA runs a pipeline of four major stages:

```
                    genome1.fa          genome2.fa
                        │                   │
                   ┌────▼────┐         ┌────▼────┐
          Stage 1  │ Compress │         │ Compress │   FASTA → compact binary
                   │ (GDB)   │         │ (GDB)   │
                   └────┬────┘         └────┬────┘
                        │                   │
                   ┌────▼────┐         ┌────▼────┐
          Stage 2  │  Index  │         │  Index  │    Build sorted k-mer index
                   │  (GIX)  │         │  (GIX)  │
                   └────┬────┘         └────┬────┘
                        │                   │
                        └───────┬───────────┘
                           ┌────▼────┐
          Stage 3          │  Seed   │                 Merge indices → seed pairs
                           │  Merge  │
                           └────┬────┘
                                │
                           ┌────▼─────────┐
          Stage 4          │ Sort, Chain,  │           Seeds → chains → alignments
                           │   & Align    │
                           └────┬─────────┘
                                │
                                ▼
                           output.1aln
                          (or PAF/PSL)
```

Think of it as: **compress → index → find seeds → align**.

Let's walk through each stage.

---

## Stage 1: Compress the Genomes (FASTA → GDB)

### What it does

A FASTA file is a text file. Each DNA base is stored as one ASCII character (A, C, G, or T) — that's 8 bits per base. But DNA only has 4 possible values, so you only need 2 bits per base. FastGA converts the text FASTA into a compact binary format called a **GDB** (Genome DataBase).

### Intuition

Think of it like compressing a book that only uses 4 letters. Instead of storing each letter in a full byte, you pack 4 letters into a single byte. A 3-billion-base genome goes from ~3 GB of text to ~750 MB of binary.

### Files produced

| File | What it contains | Size (human genome, ~3 Gbp) |
|---|---|---:|
| `genome.1gdb` | Metadata: contig names, lengths, scaffold structure | ~92 KB |
| `.genome.bps` | 2-bit compressed DNA sequence data | ~750 MB |

The `.bps` file is a "hidden" file (note the leading dot). It contains all the actual sequence data, packed 4 bases per byte. The `.1gdb` file is a small index that says "contig chr1 starts at byte offset X in the .bps file and is Y bases long."

**Why this design?** Random access. When the aligner needs to fetch bases 1000–5000 of chromosome 3, it looks up the byte offset in the `.1gdb` metadata, seeks to that position in `.bps`, and reads the compressed bases directly. No text parsing needed.

**Storage cost:** ~0.25 GB per Gbp of genome. Negligible.

---

## Stage 2: Build the Sorted K-mer Index (GDB → GIX)

This is the most storage-intensive stage and the heart of FastGA's design.

### What it does

For each genome, FastGA extracts all 40-base subsequences (40-mers) and sorts them lexicographically. The result is a **GIX** (Genome IndeX): a sorted table of every k-mer and its position in the genome.

### Intuition

Imagine a book's index at the back: "apple appears on pages 3, 17, 42." Now imagine building that index for every possible 40-character phrase in the book, sorted alphabetically. That's what GIX is — except the "phrases" are DNA sequences and the "page numbers" are genome positions.

### But not every position — Syncmer filtering

Indexing every single position would be wasteful. Many adjacent positions share almost the same k-mer (shifted by one base). FastGA uses a technique called **syncmer filtering** to select only ~50% of positions. A syncmer is a k-mer where a specific small sub-pattern appears at a designated position within the k-mer. This deterministically selects a representative subset of positions while guaranteeing that any sufficiently long match between two genomes will share at least one selected position.

### Both strands

DNA is double-stranded. A match on the "reverse complement" strand is just as biologically meaningful as a match on the forward strand. So FastGA indexes both the forward and reverse complement of every selected k-mer. This roughly doubles the index size.

### Files produced

| File | What it contains | Size (human genome, ~3 Gbp) |
|---|---|---:|
| `genome.gix` | Prefix lookup table (first 24 bits → offset) | 128 MB (always) |
| `.genome.ktab.1` ... `.ktab.N` | Sorted k-mer table partitions | ~30 GB total |
| `.genome.post.1` ... `.post.N` | Position lists for each k-mer | ~2 GB total |

**Total GIX size: ~10.1–10.5 GB per Gbp of genome.**

For two human genomes, that's **~63 GB of index files**.

### Why is the GIX so large?

Let's do the math for a 3 Gbp genome:

1. **Number of k-mers**: 3 billion positions × 50% syncmer rate × 2 strands = ~3 billion k-mers
2. **Per k-mer storage**:
   - 10 bytes: the 40-mer itself (2 bits × 40 = 80 bits, packed into 10 bytes)
   - 2 bytes: which contig this k-mer is in
   - 4 bytes: position within the contig (signed, to indicate strand)
   - ~1 byte: LCP (longest common prefix with the previous k-mer in sorted order)
3. **Total**: ~3 billion × ~11 bytes ≈ **33 GB**

The `.gix` stub file is always exactly 128 MB. It contains a prefix index: for each possible 24-bit prefix (16 million entries), it stores the offset into the sorted table where k-mers with that prefix begin. This enables fast lookup during the merge.

### Why are the k-mer tables partitioned?

Instead of one giant file, the sorted k-mers are split across N partition files (`.ktab.1` through `.ktab.N`). This is a natural consequence of the MSD radix sort: the sort distributes k-mers into buckets by their high-order bits, and each bucket becomes a partition file. It also enables parallel I/O and keeps individual files under OS limits.

### Why this design?

The k-mers must be **sorted** for the linear merge in Stage 3. A hash table would allow point lookups but would cause cache misses. A sorted table enables a single sweep. The cost is that building the sorted index takes time and disk space — but the payoff comes in Stage 3, where the merge is extremely fast.

---

## Stage 3: Seed Merge (Finding Matching Positions)

### What it does

This is where FastGA finds all positions where the two genomes share a common k-mer. It takes the two sorted GIX indices and merges them in a single linear pass, like merging two sorted lists.

### Intuition

Imagine you and a friend each have an alphabetically sorted word list. To find words that appear on both lists, you don't need to look up each word in the other's list. You simply walk through both lists together: if your current word is "apple" and theirs is "banana", advance your list (since "apple" < "banana"). If both lists say "cherry", you have a match — record it and advance both. One pass, no random lookups.

FastGA does exactly this, but with a twist: instead of requiring exact matches, it finds the **longest prefix match** for each k-mer. If genome A has the 40-mer `ACGTACGT...` and genome B has `ACGTACGG...` (matching for 38 out of 40 bases), that's still a seed. This "adaptive" matching is why they're called **adaptamers** — the match length adapts to local sequence similarity.

### Frequency filtering

If a k-mer appears thousands of times in genome B (because it's in a repetitive region), using it as a seed would generate too many spurious matches. FastGA skips any k-mer that appears more than `-f` times (default: 10) in genome B.

### Files produced

Seeds are written to temporary files partitioned by which contig in genome A they came from:

| File | What it contains | Size (human genomes) |
|---|---|---:|
| `_pair.<pid>.<k>.N` | Seed pairs (normal strand) | ~7 GB total |
| `_pair.<pid>.<k>.C` | Seed pairs (complement strand) | (included above) |

Each seed pair record contains:
- 1 byte: how many bases of the k-mer matched (the LCP)
- ~4 bytes: position in genome A
- ~4 bytes: position in genome B

**Total: ~9 bytes per seed pair.**

For two human genomes, ~1.3 billion seeds are found → ~7 GB of temp files.

### Important detail: these files are invisible

FastGA `open()`s these temp files and then immediately `unlink()`s them. In Unix, this means the file is deleted from the directory listing but remains on disk as long as the process has it open. So `ls` won't show them, `du` won't count them, but they still consume disk space. This can be surprising when you're trying to understand why your disk is full.

### Runtime

The seed merge is extremely fast: **5.3 seconds for two human genomes** at 32 threads. This is the payoff of the sorted-index design — billions of comparisons happen in a single sequential sweep.

### Critical storage observation

**After the seed merge completes, the GIX files are never read again.** The remaining stages only need the `.bps` files (for fetching sequence during alignment) and the `_pair.*` temp files (the seeds). But currently, FastGA keeps the GIX files on disk until the program exits. This means ~63 GB of GIX data sits unused during the 8-minute alignment phase.

---

## Stage 4: Sort, Chain, and Align

This stage is the most computationally intensive. It turns raw seed positions into verified, scored alignments.

### Step 4a: Radix Sort the Seeds

**What:** The seed pairs from Stage 3 are loaded into memory and re-sorted using an in-place radix sort.

**Why re-sort?** The seeds from Stage 3 are grouped by contig-in-genome-A (because that's how the merge was parallelized). But for chaining, we need them sorted by their geometric relationship: which contigs they connect and where they fall relative to each other.

FastGA transforms each seed pair into a coordinate system that makes chaining natural:

- **Anti-diagonal** = position_in_A − position_in_B (seeds on the same diagonal have the same offset between genomes — they represent a co-linear match)
- **Diagonal bucket** = (position_in_A + position_in_B) >> 7 (groups seeds along each diagonal into 128-base buckets)
- **B-contig** = which contig in genome B

After sorting by these three keys, seeds that could form a chain end up adjacent in the array.

**Why radix sort?** For the same cache-coherence reason. Comparison-based sorts (quicksort, mergesort) jump around in memory unpredictably. Radix sort processes data in sequential sweeps, one byte at a time. For billions of seeds, this is significantly faster.

### Step 4b: Chaining

**What:** Scan through the sorted seeds and group them into **chains** — sequences of seeds that are roughly co-linear (on similar diagonals, in order along both genomes).

**Intuition:** Imagine scattering dots on a 2D plot where the x-axis is position in genome A and the y-axis is position in genome B. True homologous regions show up as diagonal lines of dots. Chaining finds these lines.

A chain is formed when consecutive seeds (in sorted order) satisfy:
- They connect the same pair of contigs
- They're within 128 bases of each other in diagonal space
- They're within 1000 bases of each other along the anti-diagonal (the `-s` parameter)

Chains covering at least 85 bases (the `-c` parameter) are passed to the alignment stage.

### Step 4c: Local Alignment

**What:** For each chain, FastGA runs a full local alignment using a wave-based algorithm.

**Intuition:** The chain gives us a rough "tube" — a narrow diagonal band where the alignment probably lives. The aligner extends from this tube in both directions, matching bases and allowing insertions, deletions, and mismatches, until the similarity drops below the identity threshold (default 70%, the `-i` parameter).

The alignment engine uses **trace-point encoding**: instead of recording every single match/mismatch/insertion/deletion (like a CIGAR string does), it records a checkpoint every 100 bases. Each checkpoint says "at this point in A, we're at position Y in B, and there were Z differences in this 100-base panel." This is far more compact than a full CIGAR and can be expanded to a full alignment on demand.

### Step 4d: Redundancy Elimination and Output

Overlapping or crossing alignments are removed, keeping the best non-redundant set. The surviving alignments are sorted and written to the output file.

### Files produced

| File | What it contains | Size (human genomes) |
|---|---|---:|
| `output.1aln` | All alignments in ONEcode trace-point format | **157 MB** |

If you requested PAF output instead, the alignments are converted to text PAF format (~1.6 GB with CIGAR strings — about 10× larger than the `.1aln`).

### Why is the .1aln output so compact?

An alignment between two 5000-base regions might have a CIGAR string like `"500M2I300M1D..."` — hundreds of characters. But trace-point encoding reduces this to just 50 checkpoints (one per 100 bases), each stored as 2–4 bytes. That's ~150 bytes instead of ~500+ bytes.

For human genomes, FastGA finds ~518K alignments totaling 1.63 Gbp of aligned sequence, compressed into just 157 MB.

---

## Summary: The Complete File Inventory

Here's everything that exists on disk at each moment during a human genome comparison:

### After Stage 1 (GDB creation)
| File | Size | Cumulative |
|---|---:|---:|
| genome1.1gdb + .bps | ~750 MB | |
| genome2.1gdb + .bps | ~750 MB | |
| **Total** | | **~1.5 GB** |

### After Stage 2 (GIX build) — **Peak storage**
| File | Size | Cumulative |
|---|---:|---:|
| GDB files (both) | 1.5 GB | |
| GIX files (genome1) | ~32.5 GB | |
| GIX files (genome2) | ~30.2 GB | |
| **Total** | | **~64 GB** |

### During Stage 3 (Seed merge) — **Absolute peak**
| File | Size | Cumulative |
|---|---:|---:|
| GDB files (both) | 1.5 GB | |
| GIX files (both) | ~62.7 GB | |
| Seed pair temps | ~7 GB | |
| **Total** | | **~71 GB** |

### During Stage 4 (Sort & align) — GIX no longer needed
| File | Size | Cumulative |
|---|---:|---:|
| GDB files (both) | 1.5 GB | |
| GIX files (both, **unused**) | ~62.7 GB | |
| Seed pair temps | ~7 GB | |
| Alignment temps | small | |
| **Total** | | **~71 GB** |

### After completion
| File | Size |
|---|---:|
| output.1aln | **157 MB** |
| GDB files (if `-k` flag used) | 1.5 GB |

The key observation: **88% of peak storage is GIX files that are only read for 5 seconds**, during a process that takes nearly 10 minutes total.

---

## How FastGA Compares to Other Aligners

| Aspect | minimap2 | MUMmer/nucmer | FastGA |
|---|---|---|---|
| **Index type** | Hash table (minimizers) | Suffix array/tree | Sorted k-mer table |
| **Memory access** | Random (cache-unfriendly) | Mixed | Sequential (cache-friendly) |
| **Seed finding** | Hash lookup per k-mer | Suffix tree traversal | Linear merge of sorted lists |
| **Sort algorithm** | Comparison-based | N/A | Radix sort (no comparisons) |
| **Output format** | PAF (text) | delta/coords (text) | .1aln trace-points (binary) |
| **Output size** | Large (with CIGAR) | Medium | Very compact (15-16× smaller) |

FastGA's win comes from a fundamental architectural choice: by paying upfront cost to sort everything, every subsequent operation becomes a linear sweep. On modern CPUs where cache misses dominate performance, this matters enormously.

---

## Quick Reference: Key Parameters

| Parameter | Default | What it controls |
|---|---|---|
| `-T<int>` | 8 | Number of threads |
| `-f<int>` | 10 | Max k-mer frequency (filters repetitive seeds) |
| `-c<int>` | 85 | Min bases in a chain to trigger alignment |
| `-s<int>` | 1000 | Max gap between seeds in a chain |
| `-l<int>` | 100 | Min alignment length to report |
| `-i<float>` | 0.7 | Min alignment identity (0.7 = 70%) |
| `-k` | off | Keep intermediate files (GDB, GIX) after completion |
| `-P<dir>` | `$TMPDIR` | Where to put temporary files |
