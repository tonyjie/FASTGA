# FastGA: A Practical Guide

This guide explains FastGA from the ground up — what it does, why it's designed the way it is, and what happens to your data at every step.

---

## What Problem Does FastGA Solve?

Imagine you have two genome assemblies — say, a human reference genome and a newly assembled genome from a different individual. Each is about 3 billion DNA bases long. You want to find **all the regions where these two genomes match**, even approximately.

This is **whole-genome alignment**: find every local region in genome A that aligns to some region in genome B. The output is a list of alignments, each saying "bases 1000–5000 in contig X of genome A match bases 2000–6000 in contig Y of genome B, with 95% identity."

The challenge: doing this fast on genomes with billions of bases.

---

## Background: K-mers and How Aligners Find Matches

Before diving into FastGA's design, let's establish two foundational concepts.

### What is a k-mer?

A genome is a long string of bases: `ACGTTAGCACGT...`, billions of characters long, using only the 4-letter alphabet {A, C, G, T}. A **k-mer** is simply a substring of length k. For example, in the sequence `ACGTTAGC`:

```
Position 0:  ACGTT    ← a 5-mer
Position 1:   CGTTA   ← the next 5-mer (shifted by 1)
Position 2:    GTTAG  ← and so on
Position 3:     TTAGC
```

Every position in the genome defines a k-mer starting there. A genome with N bases has N−k+1 k-mers of length k.

**Why k-mers matter for alignment:** If two genomes share a region of high similarity, they will share many identical k-mers in that region. So the first step in any aligner is: find all k-mers that appear in both genomes. These shared k-mers are called **seeds** — they are the starting points from which full alignments are built.

FastGA uses k=40 (40-mers). Longer k-mers are more specific — a random 40-base match is astronomically unlikely, so nearly every shared 40-mer reflects genuine biological homology.

### How traditional aligners find seeds (hash tables)

The standard approach, used by tools like minimap2, works like this:

1. **Build a hash table** from genome B: for each k-mer in genome B, compute a hash (a number derived from the k-mer's sequence), and store the k-mer's position in a table slot determined by that hash.
2. **Query with genome A**: for each k-mer in genome A, compute the same hash and look up the table to see if genome B has a matching k-mer.

This is fast in theory — O(1) per lookup. But in practice, it has a serious hardware problem.

---

## The Key Idea: Cache Coherence

### Why hash table lookups cause random memory access

To understand the problem, you need to know how modern CPUs access memory. There are multiple layers:

```
┌──────────────┐
│   CPU Core   │   Registers: 1 cycle, bytes
├──────────────┤
│   L1 Cache   │   ~4 cycles,   64 KB     ← very fast, very small
├──────────────┤
│   L2 Cache   │   ~12 cycles,  256 KB
├──────────────┤
│   L3 Cache   │   ~40 cycles,  32 MB     ← shared across cores
├──────────────┤
│   Main RAM   │   ~200 cycles, 64+ GB    ← huge but slow
└──────────────┘
```

When the CPU reads a memory address, it loads a **cache line** (64 bytes) around that address into L1. If the next read is nearby (within the same or adjacent cache lines), it's served from cache — fast. If it's far away, the CPU stalls for ~200 cycles waiting for main memory.

Now consider a hash table for a human genome. The table holds billions of entries and occupies 30+ GB — far larger than any cache. When you look up k-mer #1, the hash sends you to, say, byte offset 2,400,000,000. The CPU loads that cache line from RAM (~200 cycles). Then k-mer #2 hashes to offset 18,700,000,000 — a completely unrelated location. Another cache miss, another 200 cycles. K-mer #3 goes to offset 5,100,000,000. And so on.

**Every lookup jumps to a random location in a 30 GB structure.** The CPU cache is useless because consecutive lookups are uncorrelated — by the time you revisit a region, its cache line was evicted long ago. You're paying the full ~200-cycle RAM penalty on virtually every access.

When you're doing this billions of times, those 200-cycle stalls dominate the runtime. The CPU spends most of its time waiting for data, not computing.

### FastGA's alternative: sort and sweep

FastGA avoids this entirely by never using hash tables. Instead:

1. **Pre-sort** all k-mers from each genome, lexicographically (alphabetically by their DNA sequence).
2. **Merge** the two sorted lists in a single linear pass.

With sorted data, memory access is perfectly sequential: you read element 1, then element 2, then element 3, always moving forward. The CPU's **hardware prefetcher** detects this pattern and loads upcoming cache lines *before* you need them. Every access hits cache. No stalls.

> **Instead of hashing and looking up, sort everything and sweep linearly.**

The cost is building the sorted index (which takes time and disk space). The payoff is that the actual seed-finding step becomes blazingly fast — **5.3 seconds for two entire human genomes** — because the CPU is never waiting for memory.

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

For each genome, FastGA extracts 40-base subsequences (40-mers) from selected positions and sorts them lexicographically. The result is a **GIX** (Genome IndeX): a sorted table mapping k-mer sequences to their genome positions.

### Intuition

Imagine a book's index at the back: "apple appears on pages 3, 17, 42." Now imagine building that index for every possible 40-character phrase in the book, sorted alphabetically. That's what GIX is — except the "phrases" are DNA sequences and the "page numbers" are genome positions.

### Why not index every position?

A naive approach would extract a 40-mer starting at every position in the genome. But consider what that looks like:

```
Position 0:  ACGTTAGCACGTNNNN...NNNN    (40 bases)
Position 1:   CGTTAGCACGTNNNN...NNNX    (shifted by 1 base)
Position 2:    GTTAGCACGTNNNN...NNXY    (shifted by 2)
```

Adjacent 40-mers overlap by 39 out of 40 bases. They are almost identical — they carry almost the same information. If position 0's k-mer matches something in genome B, position 1's k-mer will almost certainly match too. Indexing both is redundant work and wasted storage.

For a 3 Gbp genome, indexing every position on both strands would produce ~6 billion entries. At ~14 bytes each, that's ~84 GB per genome — impractical.

### Syncmer filtering: smart subsampling

FastGA solves this by only indexing a carefully chosen subset of positions. The technique is called **syncmer filtering**. It's important to be precise about what's happening at each level, so let's build it up step by step.

**The question being answered:** For each genome position, should we index the 40-mer starting here, or skip it?

**The answer uses a small "fingerprint" test:** Look at just 12 bases (a "TMER") at that position. Within those 12 bases, find the smallest 8-base subsequence (an "SMER"). If that smallest 8-mer happens to be at the left or right edge of the 12-mer, this position passes — index the full 40-mer here.

Let's be concrete. FastGA defines:
- **SMER = 8**: the small fingerprint
- **TMER = 12**: the test window (SMER + 4 = 12, so there are 5 places an 8-mer can sit inside a 12-mer)
- **K = 40**: the actual k-mer that gets indexed

**Step by step at one genome position:**

```
Genome:  ...NNNNNACGTTAGCACGTNNNN...NNNN...
                 ↑
                 Position i

Step 1: Look at the 12-mer starting at position i:

         A C G T T A G C A C G T        (12 bases = TMER)

Step 2: Find ALL 8-mers (SMER) within this 12-mer.
        There are 5 possible positions (12 - 8 + 1 = 5):

         Position 0:  [A C G T T A G C]  A C G T
         Position 1:   A[C G T T A G C A] C G T
         Position 2:   A C[G T T A G C A C] G T
         Position 3:   A C G[T T A G C A C G] T
         Position 4:   A C G T[T A G C A C G T]
                       ↑                       ↑
                    left edge              right edge

Step 3: Rank these five 8-mers by a canonical ordering
        (considering both forward and reverse complement,
         so it's strand-independent).

        Suppose position 0's 8-mer ranks smallest.
        Position 0 IS the left edge → PASS ✓

        → Index the full 40-mer starting at position i.

Step 4: If instead position 2's 8-mer ranked smallest:
        Position 2 is the middle — neither edge → FAIL ✗

        → Skip this position. Don't index the 40-mer here.
```

**To be clear about what is what:**
- The **syncmer** is the 40-mer that gets indexed. It's the thing that goes into the GIX.
- The **12-mer (TMER)** is the test window used to decide whether to select this position.
- The **8-mer (SMER)** is the small fingerprint used within the test.
- The 12-mer and 8-mer are NOT stored in the index. They are only used for the pass/fail decision.

**Why the edge rule works — guaranteed coverage:**

Imagine two genomes share a long identical region. As we slide position-by-position through this region, the 12-mer window changes by one base at each step. The minimum 8-mer within the window can only change when either:
- A new 8-mer enters from the right and happens to be smaller → new minimum at the **right edge** → syncmer selected
- The old minimum falls off the left end → the departing minimum was at the **left edge** → syncmer was already selected at that step

So between any two "minimum changes," at least one syncmer position is selected. This guarantees that within any shared stretch of ~12+ identical bases, both genomes will select at least one syncmer in the same position. No true match can be missed.

**Density:** Roughly **35–40% of positions are selected** as syncmers. (Theory for SMER=8 / TMER=12: the 12-mer window holds 5 overlapping 8-mers, and the position is kept only when the minimal 8-mer sits at one of the two ends — about 2/5 of the time.) Since every selected position emits **two** entries (forward + reverse complement, see below), the index ends up with ~0.7 entries per genome base — i.e. the total number of index entries is roughly **70% of the genome size**, which for GRCh38 is the measured ~2.2 billion entries.

> Note: FastGA's verbose output prints a `Sampled: N (X%)` line reporting ~75–81%. That percentage is **entries per position counting both strands**, *not* the per-position selection rate — with ~2 entries per selected position it works out to ~0.7–0.8 per position, consistent with the ~35–40% single-strand selection above. Don't read the `Sampled %` as "fraction of positions kept."

### Both strands

DNA is double-stranded. The two strands are reverse complements of each other: if one strand reads `ACGT→`, the other reads `←TGCA` (A↔T, C↔G, reversed). A biological match can occur on either strand.

For every position that passes the syncmer filter, FastGA creates **two** index entries: one for the forward k-mer and one for its reverse complement. The contig ID field uses a sign bit to distinguish them (positive = forward, negative = reverse complement).

### Files produced

Now, what exactly goes into the files? The GIX consists of three types of files:

#### 1. The `.gix` stub file (always 128 MB)

This is a **prefix jump table** — NOT a hash table. It's a simple array indexed by the first 12 bases of a k-mer, giving the offset into the sorted ktab where k-mers with that prefix begin.

Think of it like the thumb tabs on a dictionary. A dictionary is sorted alphabetically (like the ktab), and the thumb tabs let you jump directly to "M" without flipping from the beginning. The `.gix` file is the thumb tab — except instead of 26 letter tabs, it has 16 million tabs (one for every possible 12-base DNA prefix).

```
Structure of genome.gix:

┌─────────────────────────────────┐
│ Header: KMER (40), NPARTS, etc. │  16 bytes
├─────────────────────────────────┤
│ Prefix jump table: 16M entries  │  16,777,216 × 8 bytes = 128 MB
│                                 │
│  Index encodes first 12 bases of the k-mer as a 24-bit number:
│  A=00, C=01, G=10, T=11, so AAAAAAAAAAAA = 0x000000, etc.
│                                 │
│  table[0x000000] = 0            │  "k-mers starting with AAAAAAAAAAAA begin at ktab offset 0"
│  table[0x000001] = 1547         │  "k-mers starting with AAAAAAAAAAAC begin at offset 1547"
│  ...                            │
│  table[0xFFFFFF] = ...          │  "k-mers starting with TTTTTTTTTTTT begin at offset ..."
├─────────────────────────────────┤
│ Metadata: PostBytes, ContBytes, │  small
│ contig permutation, etc.        │
└─────────────────────────────────┘
```

**How it's used during the merge (Stage 3):** During the linear sweep, FastGA doesn't actually need to look up random k-mers — it sweeps through both sorted tables sequentially. The prefix table is mainly used during the GIX build to efficiently partition and reassemble the sorted k-mers.

**Why exactly 128 MB?** 12 bases × 2 bits/base = 24 bits → 2^24 = 16,777,216 possible prefixes × 8 bytes per offset = 128 MB. This is a fixed cost regardless of genome size — the same 128 MB whether the genome is 1 Mbp or 3 Gbp.

#### 2. The `.ktab.*` partition files (~30 GB total for human genome)

These contain the actual sorted k-mer data — the core of the index. Before diving into the byte layout, let's clarify three genomics terms that appear in every entry.

**What is a contig?** A genome assembly is not one continuous string. Current sequencing technology produces the genome in pieces. A **contig** (short for "contiguous sequence") is one such piece — an unbroken stretch of DNA whose sequence is known. The human genome assembly has ~700 contigs. The largest (parts of chromosome 1) are ~249 million bases long. The smallest might be a few thousand bases. Think of contigs as chapters in a book — the genome is the whole book, but you read it one chapter at a time.

**What is a strand?** DNA is a double helix — two strands running in opposite directions, with complementary bases (A pairs with T, C pairs with G). If the forward strand reads `ACGT`, the reverse complement strand reads `ACGT` backwards with each base flipped: `ACGT → TGCA`. A match between two genomes can occur on either strand. When we say a k-mer is from the "forward strand" at position 500, we mean the 40 bases starting at position 500 read left-to-right. The "reverse complement" at position 500 means the 40 bases ending at position 500, read right-to-left with bases flipped.

**What is position?** The offset (in bases) from the start of a contig. Position 0 is the first base of the contig, position 1000 means 1000 bases in. Combined with the contig ID and strand, this uniquely identifies where in the genome a k-mer came from.

Now, each ktab entry on disk (original FastGA code, before our optimizations):

```
Per ktab entry (human genome, original code):

┌──────────────────────┬────────────┬───────────┬──────────────┬────────────────┐
│  K-mer suffix        │  Mask byte │  LCP      │  Position    │  Contig+Strand │
│  (remaining bases    │  (always 0 │           │  (offset     │  (which contig,│
│   after prefix)      │   when no  │           │   within     │   which strand)│
│                      │   masking) │           │   contig)    │                │
├──────────────────────┼────────────┼───────────┼──────────────┼────────────────┤
│  7 bytes             │  1 byte    │  1 byte   │  4 bytes     │  2 bytes       │
│  (28 bases × 2 bits  │            │           │  (supports   │  (supports up  │
│   = 56 bits, packed  │            │           │   contigs    │   to ~32K      │
│   into 7 bytes)      │            │           │   up to 4 GB)│   contigs,     │
│                      │            │           │              │   sign=strand) │
└──────────────────────┴────────────┴───────────┴──────────────┴────────────────┘
  In-memory during sort: 17 bytes (adds LCP byte at front + 2 prefix bytes)
  On disk in .ktab file: 15 bytes per entry
```

*(Note: this 15-byte layout is the **original** upstream code. Our optimizations shrink it:
Optimization 3 removes the always-zero mask byte when masking is off (→ 14 bytes), and
Optimization 4 additionally drops the LCP byte and recomputes it on the fly (→ 13 bytes).
The original code always stores all 15.)*

Let's go through each field:

**K-mer suffix (7 bytes):** The full 40-mer is 80 bits = 10 bytes. But not all 10 bytes need to be stored in each entry:
- The first 4 bases (8 bits) are encoded by which **partition** the entry was sorted into (the bucket index).
- The next 8 bases (16 bits = 2 bytes) are absorbed into the **prefix jump table** within each partition.
- That leaves 28 bases (56 bits) — which fits exactly in 7 bytes.

So 12 bases are implicit in the table structure, and only the remaining 28 bases are stored per entry.

**Mask byte (1 byte):** This byte is meant to store soft-masking information (indicating repetitive regions in the genome). In the original code, this byte is **always allocated** — even when masking is disabled (the `-M` flag is not set), every entry carries a mask byte set to 0. This is a design choice that simplifies the code (every entry has the same width regardless of masking mode), but wastes 1 byte per entry. Across billions of entries, that's gigabytes of wasted disk.

**Position (4 bytes):** Which base in the contig this k-mer starts at. The code computes `PostBytes` dynamically — for human chromosomes up to ~249 Mbp, 4 bytes (max ~4.3 billion) is sufficient. For a tiny genome with short contigs, this might be only 2 or 3 bytes.

**Contig+Strand (2 bytes):** Identifies which contig and which strand. The sign bit encodes strand: positive = forward, negative = reverse complement. For the human genome's ~700 contigs, 2 bytes (max ~32K values) is plenty. Also computed dynamically (`ContBytes`).

**LCP (1 byte):**

LCP stands for **Longest Common Prefix**. For each entry in the sorted table, the LCP records how many bases this k-mer shares with the **previous** entry in sorted order.

Example in a sorted k-mer table:
```
Entry 100:  ACGTTAGCACGT...  (40 bases)
Entry 101:  ACGTTAGCCCGT...  (40 bases)    LCP = 32  (first 32 bases identical)
Entry 102:  ACGTTAGCCCTT...  (40 bases)    LCP = 36  (first 36 bases identical)
Entry 103:  ACTTGGGCACGT...  (40 bases)    LCP = 3   (only ACT in common)
```

The LCP is computed as a byproduct of the MSD radix sort — it falls out naturally during sorting. During the compress step before writing to disk, the LCP is repositioned: in memory during sorting it sits at byte 0 (before the k-mer), but on disk it's placed between the k-mer suffix and the position field. Entries with identical k-mers (duplicates from different positions) get an LCP of 40 (= full match with predecessor).

**Why the LCP is crucial for the algorithm:**

The LCP is what makes the Stage 3 merge fast. During the merge, FastGA sweeps through both sorted tables (genome A and genome B) simultaneously. When it finds a k-mer in genome A, it needs to find all k-mers in genome B that share a long prefix with it. The LCP values let it do this without comparing full 40-mers character by character.

Here's the key insight: in a sorted list, all k-mers sharing a prefix form a contiguous block. The LCP values tell you exactly where each block begins and ends. When the LCP drops below some value N, you know you've left the block of k-mers sharing an N-base prefix.

```
Scanning genome B's sorted table during merge:

  Entry:   ...ACGTTAGCACGT...    LCP=40  ─┐
  Entry:   ...ACGTTAGCACGT...    LCP=40   │ All share 38-base prefix
  Entry:   ...ACGTTAGCACTT...    LCP=38   │ with genome A's current k-mer
  Entry:   ...ACGTTAGCCCGT...    LCP=32  ─┘
  Entry:   ...ACTTGGGCACGT...    LCP=3   ← LCP drops below 32 → block ends
```

Without LCP values, you'd have to compare every k-mer byte-by-byte to know when you've exited the matching block. With LCP values, you just check one number. This is what enables the "adaptamer" matching — FastGA finds the longest prefix match length and all positions sharing that prefix, in a single forward scan.

### How ktab entries are used in the algorithm

The sorted ktab is consumed during **Stage 3 (Seed Merge)**. Here's how:

```
Genome A's ktab (sorted):          Genome B's ktab (sorted):
  AAACCCGGGTTT... pos=100 chr1      AAACCCGGGTTT... pos=200 chr3
  AAACCCGGGTTX... pos=500 chr1      AAACCCGGGUUU... pos=800 chr5
  AAACCCTTTTTT... pos=300 chr2      AAACCCTTTTTT... pos=150 chr2
  ...                                ...

Merge sweep (both advance left to right):
  ─────────────────────────────────────────────────────
  A's current k-mer: AAACCCGGGTTT...
  B's matching block: entries with prefix AAACCCGGG... (found via LCP)
     → B has 2 entries matching with prefix length ≥ 36
     → Both have freq ≤ 10, so emit seed pairs:
        (A: chr1, pos=100)  ↔  (B: chr3, pos=200)
        (A: chr1, pos=100)  ↔  (B: chr5, pos=800)

  These seed pairs mean: "position 100 in A's chr1 and position 200
  in B's chr3 share a 36+ base prefix — probably a real homologous region.
  Investigate further in Stage 4."
  ─────────────────────────────────────────────────────
```

The position and contig+strand fields are what make the seeds useful: they tell Stage 4 exactly *where* in each genome to look for a full alignment. The k-mer suffix tells you *how long* the match is (via LCP comparison). Together, they define a candidate alignment anchor.

**Why partitioned into multiple files?** The MSD radix sort processes k-mers in batches that fit in ~4 GB of RAM. Each batch, once sorted, is written to its own `.ktab.N` file. The number of partitions (N) depends on genome size — typically 8 to 64. This keeps memory usage bounded and enables parallel I/O.

#### 3. The `.post.*` position list files (~2 GB total for human genome)

These are used during the GIX build phase as intermediate data — they store the raw (unsorted) positions of syncmer k-mers, which are then read back and sorted into the ktab files. In the final GIX, the position information is embedded in the ktab entries, so the `.post.*` files are primarily a build artifact. They use variable-length encoding (1–4 bytes per position, delta-compressed) to save space.

### Storage math: why ~10 GB per Gbp

Let's trace the numbers for a 3 Gbp human genome:

```
Genome size:           3,100,000,000 bases (GRCh38)
Syncmer selection:     ~35% of positions kept (per strand)
Entries per position:  × 2 (forward + reverse complement)
Total index entries:   3.1B × 0.35 × 2 ≈ 2.2 billion entries
                       (measured: ~2.2 billion for GRCh38)

  (Each entry is stored individually — even if two positions share
   the same k-mer sequence, both get their own entry with their
   own position/contig fields. Duplicates get LCP = 40.)

Per entry on disk:     15 bytes (7 suffix + 1 mask + 1 LCP + 4 position + 2 contig)
K-mer table total:     2.2B × 15 ≈ 33 GB
Plus .gix stub:        0.128 GB

Total GIX:             ~32.5 GB per genome (measured: 32.5 GB for GRCh38)
Per Gbp:               ~10.5 GB/Gbp
```

For **two** human genomes: ~63 GB of index files. This is 88% of FastGA's peak disk usage.

### Why this design?

The k-mers must be **sorted** for the linear merge in Stage 3. A hash table would allow O(1) point lookups but would cause cache misses (as explained earlier). A sorted table enables a single sequential sweep. The cost is building the sorted index (time and disk space) — but the payoff comes in Stage 3, where the merge is extremely fast: 5.3 seconds for two entire human genomes.

---

## Stage 3: Seed Merge (Finding Matching Positions)

### What it does

This is where FastGA finds all positions where the two genomes share a common k-mer. It takes the two sorted GIX indices and merges them in a single linear pass, like merging two sorted lists.

### Intuition

Imagine you and a friend each have an alphabetically sorted word list. To find words that appear on both lists, you don't need to look up each word in the other's list. You simply walk through both lists together: if your current word is "apple" and theirs is "banana", advance your list (since "apple" < "banana"). If both lists say "cherry", you have a match — record it and advance both. One pass, no random lookups.

FastGA does exactly this, but with a twist: instead of requiring exact matches, it finds the **longest prefix match** for each k-mer. If genome A has the 40-mer `ACGTACGT...` and genome B has `ACGTACGG...` (matching for 38 out of 40 bases), that's still a seed. This "adaptive" matching is why they're called **adaptamers** — the match length adapts to local sequence similarity.

### How matching works: longest prefix + frequency filtering

For each k-mer in genome A, the merge finds the **longest prefix match** against genome B's sorted table. The algorithm extends the match base by base — starting from the 12-base prefix group — until T2 has no entries that agree on the next base. The resulting match length `plen` can be anything from 12 (only the prefix group matches) to 40 (exact full k-mer match).

**There is no minimum match length requirement.** Instead, FastGA uses a single filter: **frequency**. It counts how many entries in genome B share the longest prefix match. If that count exceeds `-f` (default: 10), the k-mer is skipped. The match length `plen` is recorded in the seed pair output, so Stage 4 knows how long the match was.

**Why frequency filtering implicitly handles match quality:**

At first this seems strange — wouldn't short matches (e.g., 13 bases) produce meaningless seeds? In practice, the frequency filter takes care of this automatically:

```
Match length vs. expected frequency (human genome, ~2.2B indexed k-mers):

  40-base match:  4^40 ≈ 10^24 possibilities
                  → Almost every 40-mer is globally unique
                  → Frequency typically 1 → KEPT (strong seed)

  25-base match:  4^25 ≈ 10^15 possibilities
                  → Still astronomically specific
                  → Frequency usually 1-3 → KEPT (good seed)

  13-base match:  4^13 ≈ 67 million possibilities
                  → 2.2 billion entries / 67 million prefixes ≈ 33 per prefix
                  → Frequency typically ~30 → SKIPPED (exceeds -f 10)

  12-base match:  4^12 ≈ 16.7 million possibilities
                  → ~130 entries per prefix on average → SKIPPED
```

Longer prefix matches are inherently more specific, so they naturally have low frequency and pass the filter. Short prefix matches are inherently common, so they naturally have high frequency and get filtered out. The frequency threshold (`-f 10`) acts as both a repetitiveness filter AND an implicit match-quality filter — one parameter handles both.

The rare exception: a short match (say, 15 bases) in a truly unique region of the genome *could* pass the frequency filter. But this is uncommon, and even if such a seed slips through, Stage 4's chaining will discard it unless it lines up with other nearby seeds.

### Files produced

Seeds are written to temporary files, partitioned by which contig in genome A they came from and by **strand orientation**:

| File | What it contains | Size (human genomes) |
|---|---|---:|
| `_pair.<pid>.<k>.N` | Seed pairs, same-direction matches | ~7 GB total |
| `_pair.<pid>.<k>.C` | Seed pairs, inverted matches | (included above) |

**What do `.N` and `.C` mean?**

Recall that each genome position is indexed on both strands (forward and reverse complement). When a seed pair is found, the two k-mers might come from the same strand orientation or opposite orientations:

```
Normal (.N) — same direction:
  Genome A:  ───────ACGTTAGC──────→    (forward)
  Genome B:  ───────ACGTTAGC──────→    (forward)
  The matching regions run in the same direction in both genomes.

Complement (.C) — inverted:
  Genome A:  ───────ACGTTAGC──────→    (forward)
  Genome B:  ←──────GCTAACGT───────    (reverse complement)
  One region is flipped relative to the other.
```

This separation matters for Stage 4: same-direction matches and inverted matches require different chaining geometry. By writing them to separate files, Stage 4 can process each orientation independently without needing to sort or filter by strand.

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

**Why re-sort?** The seeds from Stage 3 are grouped by contig-in-genome-A (because that's how the merge was parallelized). But for chaining, we need seeds from the same alignment to be *adjacent* in memory. This requires a different sort order based on geometric relationships.

**First, some intuition.** Think of a seed pair as a dot on a 2D plot:

```
Position in B
     ▲
     │          • seed 3
     │        • seed 2          ← These 4 seeds form a diagonal line.
     │      • seed 1               They probably represent a real
     │    • seed 0                 alignment between the genomes.
     │
     │                  • seed X  ← This seed is far off — different
     │                               alignment or noise.
     └──────────────────────────→ Position in A
```

Seeds from a real alignment sit along a **diagonal line** (slope ≈ 1), because both genomes advance at roughly the same rate through homologous regions. Our goal is to sort seeds so that dots on the same diagonal line end up next to each other in the array.

**The coordinate transformation.** FastGA converts each (posA, posB) seed pair into three sort keys:

**Key 1 — Diagonal** = posA − posB (for normal-strand matches):

```
Position in B
     ▲
     │  /  /  /  /  /          Each "/" line is a diagonal where
     │ /  /  /  /  /           posA - posB = constant.
     │/  /  /  /  /
     │  /  /  /  /             Seeds on the same diagonal have the
     │ /  /  /  /              same offset between the two genomes.
     │/  /  /  /
     └──────────────→ Position in A

  seed 0: posA=1000, posB=2000 → diagonal = -1000
  seed 1: posA=1040, posB=2038 → diagonal = -998   ← almost same diagonal
  seed 2: posA=1080, posB=2079 → diagonal = -999   ← almost same diagonal
  seed X: posA=5000, posB=100  → diagonal = 4900   ← completely different diagonal
```

Seeds from the same alignment have *nearly* the same diagonal value (not exactly, because of insertions/deletions). To handle this, FastGA groups diagonals into **buckets** of width 64 (= 2^BUCK_SHIFT). Seeds within 64 bases of each other in diagonal space land in the same bucket or an adjacent one.

The reimport stores:
- The **diagonal bucket** = diagonal >> 6 (which 64-base band)
- The **diagonal remainder** = diagonal mod 64 (position within the band)

**Key 2 — Anti-diagonal** = posA + posB:

```
Position in B
     ▲
     │\  \  \  \  \           Each "\" line is an anti-diagonal where
     │ \  \  \  \  \          posA + posB = constant.
     │  \  \  \  \  \
     │   \  \  \  \  \        The anti-diagonal measures progress
     │    \  \  \  \  \       ALONG a diagonal.
     └──────────────────→ Position in A

  seed 0: posA=1000, posB=2000 → anti-diagonal = 3000
  seed 1: posA=1040, posB=2038 → anti-diagonal = 3078   ← 78 bases further along
  seed 2: posA=1080, posB=2079 → anti-diagonal = 3159   ← 159 bases further along
```

While the diagonal tells you *which* line the seed is on, the anti-diagonal tells you *where along* that line the seed sits. Sorting by anti-diagonal within a diagonal bucket orders seeds by their position along the alignment.

**Key 3 — B-contig:** Which contig in genome B. Seeds connecting different pairs of contigs can never be part of the same alignment, so they must be separated.

**The sort order:** Seeds are sorted by (A-contig, B-contig, diagonal bucket, anti-diagonal). After sorting, seeds from the same alignment — same contig pair, same diagonal band, ordered along the alignment — end up contiguous in memory. This makes chaining a simple linear scan.

**Why radix sort?** For the same cache-coherence reason. Comparison-based sorts (quicksort, mergesort) jump around in memory unpredictably. Radix sort processes data in sequential sweeps, one byte at a time. For billions of seeds, this is significantly faster.

### Step 4b: Chaining

**What:** Scan through the sorted seeds and group them into **chains** — sequences of seeds that are co-linear (on similar diagonals, in order along both genomes).

Thanks to Step 4a's sort, seeds from the same alignment are now adjacent. Chaining is a simple linear scan: walk through the sorted array and group consecutive seeds that are "close enough."

A chain grows by accumulating seeds one at a time. The current chain extends as long as each next seed satisfies:
- Same pair of contigs
- Same or adjacent diagonal buckets (within 128 bases in diagonal space, i.e., 2 × BUCK_WIDTH)
- Anti-diagonal gap ≤ 1000 bases from the previous seed (the `-s` parameter) — no large gap along the alignment

When a seed violates any of these, the current chain ends and a new chain begins.

A chain can contain **any number of seeds** — 2, 5, 20, or more. The chain's **coverage** is the sum of the match lengths (`plen` values from Stage 3) of all its seeds. Only chains with coverage ≥ 85 bases (the `-c` parameter) trigger alignment:

```
Chain with 5 seeds (all in the same diagonal band):

  Seed 0: anti-diag=3000, plen=38  ─┐
  Seed 1: anti-diag=3078, plen=40   │  gap=78   ✓ (< 1000)
  Seed 2: anti-diag=3159, plen=35   │  gap=81   ✓
  Seed 3: anti-diag=3300, plen=40   │  gap=141  ✓
  Seed 4: anti-diag=3502, plen=36  ─┘  gap=202  ✓

  Total coverage: 38+40+35+40+36 = 189 bases  ✓ (≥ 85)
  → Trigger alignment over this region

Chain with 1 seed:

  Seed X: anti-diag=8000, plen=40  ─┐
  (next seed is gap > 1000)         ─┘

  Total coverage: 40 bases  ✗ (< 85)
  → Discard — not enough evidence
```

This is why the `-c 85` threshold matters: a single seed (even a perfect 40-base match) isn't enough evidence to justify the expensive alignment step. You need multiple seeds lining up consistently to be confident there's a real homologous region worth aligning.

### Step 4c: Local Alignment

For each qualifying chain, FastGA runs a full local alignment. The individual seeds are NOT aligned one by one. Instead, the chain defines a **tube** — a narrow diagonal band in the posA × posB space — and the aligner searches for the best alignment within that tube.

#### From chain to tube

The chain's seeds collectively define the tube:

```
Position in B
     ▲
     │        ╔══════════════════╗
     │        ║           • s4   ║
     │        ║         • seed 3 ║    Width  = dgmax - dgmin (diagonal range)
     │        ║      • seed 2    ║  ← The "tube": a narrow diagonal band
     │        ║    • seed 1      ║    defined by the chain's extent.
     │        ║  • seed 0        ║    Length = ahgh - alow  (anti-diagonal range)
     │        ╚══════════════════╝
     │
     └──────────────────────────────→ Position in A

Seeds progress from lower-left to upper-right: both posA and posB increase,
so anti-diagonal (posA + posB) increases along the chain, while diagonal
(posA - posB) stays roughly constant (same diagonal band).
```

The tube is defined by four numbers:
- **Diagonal range** [dgmin, dgmax]: the width of the band (from the chain's min/max diagonal values).
- **Anti-diagonal range** [alow, ahgh]: the length of the tube (how far the chain extends).

Once the tube is defined, the individual seeds are discarded — they've served their purpose. What happens next operates on the **raw contig sequences** directly.

#### What sequences does the aligner work on?

The entire contigs (from the `.bps` files) are loaded into memory:

```
Contig in genome A:  ──────────────────────────────────────  (e.g., 249 million bases)
Contig in genome B:  ────────────────────────────────────    (e.g., 200 million bases)

But the alignment only touches this small region:

Contig A:  ──────────[====tube====]─────────────────────────
Contig B:  ─────────────[====tube====]──────────────────────
                         ↑
                    ~5000 bases wide, ~128 diagonals
```

The full contigs must be in memory because the aligner needs random access to bases within the tube region. But the algorithm is constrained to the tube's diagonal band — it never visits cells outside it. For a 249M-base contig, the tube might cover only 5000 bases × 128 diagonals = 640K cells, an infinitesimal fraction of the full 249M × 200M matrix.

#### The DP matrix and why we need something smarter than Smith-Waterman

The alignment problem is: given two sequences and a diagonal band, find the best local alignment within that band. This is a dynamic programming problem on a 2D grid:

```
         Sequence A →
         A  C  G  T  T  A  G  C  A  C  G  T
     ┌──────────────────────────────────────
  A  │  ╲
  C  │     ╲           The optimal alignment is a path
  G  │        ╲        through this grid.
  T  │           ╲
  T  │              ╲     ╲ = match (diagonal move)
  A  │                 ╲   | = insertion (move down)
  G  │                    ╲  — = deletion (move right)
  C  │                       ╲
  A  │                          ╲
  C  │
```

**Smith-Waterman** would fill every cell in the band, row by row, regardless of how similar the sequences are. For a 5000-base alignment in a 128-diagonal band, that's 640K cells — always, whether the sequences are 99% identical or 70% identical.

FastGA uses a fundamentally different approach: **Gene Myers' O(nd) wave-front algorithm** (1986). Instead of sweeping through *positions*, it sweeps through *edit distances*.

#### How the wave-front algorithm works

The key idea: **start at a point and greedily extend exact matches. Only compute further when you hit a mismatch.**

```
Wave d=0 (no errors allowed):
  Starting from the tube's anti-diagonal midpoint, walk along the diagonal
  comparing bases: A[i] == B[j]? A[i+1] == B[j+1]? ...
  Extend as far as possible with zero mismatches.

  ────────────•••••••••••••|────────
              ↑ start      ↑ first mismatch after 13 exact matches

Wave d=1 (one error allowed):
  From every point reached by wave 0, try each possible error:
    - mismatch: skip one base in both A and B
    - insertion: skip one base in B only
    - deletion: skip one base in A only
  Then greedily extend exact matches again from each new position.

  ──────•••••••••••••••••••••••|────
        ↑ reached further (one error, then more exact matches)

Wave d=2 (two errors allowed):
  Same process: from wave 1's endpoints, try one more error, extend.

  ──••••••••••••••••••••••••••••••••
    ↑ even further

  ...continues until the error rate exceeds the identity threshold.
```

**Why this is faster than Smith-Waterman:** Each wave is cheap — it's mostly just comparing bases in a loop (the greedy exact-match extension). The number of waves equals the number of differences `d`, not the alignment length `n`. For high-identity alignments (which is the common case in whole-genome comparison):

```
  95% identity, 5000 bases:  d ≈ 250 errors
    Smith-Waterman (banded): 5000 × 128 = 640K cells      (always)
    Wave-front:              O(n × d) ≈ 5000 × 250 = ~1.25M ops
                             BUT most of each wave is greedy extension
                             (just base comparisons), so much faster in practice

  99% identity, 5000 bases:  d ≈ 50 errors
    Smith-Waterman (banded): still 640K cells               (same as before!)
    Wave-front:              O(5000 × 50) ≈ 250K ops       (5× less work)
```

The wave-front naturally does less work when sequences are more similar — exactly the situation in genome alignment. Smith-Waterman does the same work regardless of similarity.

#### Two passes: forward and reverse

`Local_Alignment` starts at a point within the tube (derived from the chain's anti-diagonal range) and runs two passes:

1. **`forward_wave`**: extends forward (increasing posA and posB) to find where the alignment ends.
2. **`reverse_wave`**: extends backward (decreasing posA and posB) to find where the alignment begins.

```
Along the anti-diagonal (posA + posB increasing →):

                      starting point
                          ↓
   ◄── reverse_wave ──── • ──── forward_wave ──►

   Result: abpos,bbpos          Result: aepos,bepos
   (alignment start)            (alignment end)
```

The result may extend beyond the chain's seeds (the aligner discovered more homology) or be shorter (part of the region had too many differences):

```
Chain seeds:              ••••••••
Actual alignment:  ════════════════════════
                   ↑                      ↑
               Extended beyond         Extended beyond
               first seed              last seed
```

#### Quality filtering

The resulting alignment is accepted only if:
- Length ≥ `-l` (default 100 bases)
- Identity ≥ `-i` (default 0.7 = 70%)

#### How the result is encoded — trace points

Rather than recording every match/mismatch/insertion/deletion (like a CIGAR string), the aligner records a compact **trace-point encoding** — a checkpoint every 100 bases in A (the trace spacing, TS=100):

```
Alignment from A[1000..5000] to B[2000..6020]:

  Trace points (one every 100 bases in A):
    at A[1100]: B is at 2105, 3 differences in this panel
    at A[1200]: B is at 2202, 2 differences in this panel
    at A[1300]: B is at 2298, 4 differences in this panel
    ...
    at A[4900]: B is at 5915, 1 difference in this panel

  Stored as pairs: (diffs, b_offset) = (3,105), (2,97), (4,96), ..., (1,95)
  Each pair is 2 unsigned shorts = 4 bytes per 100 bp of alignment.
```

This is extremely compact: a 4000-base alignment needs only 40 trace points × 4 bytes = 160 bytes, vs. a CIGAR string that might be 500+ characters.

**But how can you reconstruct the exact alignment from just these summaries?**

The trace points do NOT record where individual mismatches, insertions, or deletions are. They only record the *count* of differences per panel and the B-position at each checkpoint. So how do you recover the base-by-base alignment when you need it?

The answer: **re-align each panel independently.** This is what `Compute_Trace` does:

```
Trace point says:
  Panel: A[1100..1200] ↔ B[2105..2202], 2 differences

To find WHERE those 2 differences are:
  Run Smith-Waterman on just this panel:
    A[1100..1200]  (100 bases)
    B[2105..2202]  (97 bases)
    → 100 × 97 = ~10K DP cells  (trivial)

  Result: exact alignment, e.g.:
    A: ACGTTAGC-ACGT...    (100 bases)
    B: ACGTTAGCAACGT...    ( 97 bases)
                ↑
           1 insertion + 1 mismatch somewhere = 2 diffs ✓
```

This works because the trace points give you the **exact boundaries** of each panel — you know precisely which ~100 bases in A align to which ~100 bases in B. So each reconstruction is a tiny global alignment (fixed start and end), not an open-ended search. A 5000-base alignment becomes 50 independent ~100-base panels, each taking ~10K cells.

**The design tradeoff:**

```
Full CIGAR:     Store everything explicitly.
                Large files (500+ chars per alignment).
                Instant access to any operation.

Trace points:   Store only summaries (4 bytes per 100 bases).
                Tiny files (15-16× smaller).
                Reconstruct details on demand (~10K cells per panel).
```

FastGA chooses trace points because most downstream analysis doesn't need base-level detail for every alignment — and when it does, reconstruction is fast.

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
