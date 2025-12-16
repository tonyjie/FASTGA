# FastGA Execution Phases Explained

## Overview
FastGA performs whole genome alignment between two genomes (HAP1 and HAP2) through several distinct phases. This document explains each phase, its timing, and output meaning based on a complete run from scratch.

## Command Used
```bash
FastGA -vk -1:H1vH2 HAP1 HAP2
```

**Flags:**
- `-v`: Verbose mode (shows progress)
- `-k`: Keep intermediate files (GDB and GIX files)
- `-1:H1vH2`: Output format is ONEcode ALN format to file `H1vH2.1aln`

---

## Complete Phase Breakdown

### **Phase 1: GDB Creation for HAP1 (FAtoGDB)**
**Status:** ✅ Executed

**What it does:**
- Converts FASTA file (`HAP1.fasta.gz`) to Genome Database (`HAP1.1gdb`) format
- Creates a binary representation of the genome assembly
- Stores sequences in 2-bit compressed format in hidden `.HAP1.bps` file
- Extracts scaffold/contig information, names, and sizes
- Enables random access to sequences with 4x less I/O than FASTA

**Output files:**
- `HAP1.1gdb` - Main database file (~1.3KB metadata)
- `.HAP1.bps` - Hidden file with compressed sequences (~21MB)

**Timing:**
- **User time:** 0.424s
- **System time:** 0.027s
- **Wall time:** 0.454s
- **CPU utilization:** 99.4%
- **Peak memory:** 0MB (reported)

**Resource format:** `user_time system_time wall_time CPU_utilization%`

---

### **Phase 2: GIX Creation for HAP1 (GIXmake)**
**Status:** ✅ Executed (2 sub-phases)

#### **Sub-phase 2a: K-mer Partitioning**
**What it does:**
- Partitions k-mers via position lists into 8 parts for parallel processing
- Prepares k-mer data for sorting

**Timing:**
- **User time:** 1.273s
- **System time:** 0.116s
- **Wall time:** 0.295s
- **CPU utilization:** 471.2% (using ~4.7 cores)

#### **Sub-phase 2b: Sort & Index Output**
**What it does:**
- Sorts each of the 8 parts
- Compresses sorted k-mer data
- Outputs sorted k-mer tables to hidden `.ktab.N` files
- Creates position lists in hidden `.post.N` files
- Samples k-mers (75% kept for indexing)

**Process steps:**
1. Sorting part 1-8 (parallel)
2. Compressing part 1-8 (parallel)
3. Outputting part 1-8 (parallel)

**Output files:**
- `HAP1.gix` - Main index file (proxy, ~129MB)
- `.HAP1.ktab.1` through `.HAP1.ktab.8` - Hidden k-mer table files (~100MB each)
- `.HAP1.post.1` through `.HAP1.post.8` - Hidden position list files

**Timing:**
- **User time:** 6.148s
- **System time:** 1.041s
- **Wall time:** 2.083s
- **CPU utilization:** 345.1% (using ~3.5 cores)

**Statistics:**
- **Sampled:** 64,661,214 k-mers/positions (75.0% of total)

#### **Total Phase 2 (HAP1):**
- **User time:** 7.422s
- **System time:** 1.158s
- **Wall time:** 2.378s
- **CPU utilization:** 360.8%

---

### **Phase 3: GDB Creation for HAP2 (FAtoGDB)**
**Status:** ✅ Executed

**What it does:**
- Same as Phase 1, but for HAP2 genome
- Converts `HAP2.fasta.gz` to `HAP2.1gdb`

**Output files:**
- `HAP2.1gdb` - Main database file (~1.3KB)
- `.HAP2.bps` - Hidden file with compressed sequences (~21MB)

**Timing:**
- **User time:** 0.425s
- **System time:** 0.026s
- **Wall time:** 0.454s
- **CPU utilization:** 99.5%

---

### **Phase 4: GIX Creation for HAP2 (GIXmake)**
**Status:** ✅ Executed (2 sub-phases)

#### **Sub-phase 4a: K-mer Partitioning**
**Timing:**
- **User time:** 1.303s
- **System time:** 0.104s
- **Wall time:** 0.280s
- **CPU utilization:** 502.6% (using ~5 cores)

#### **Sub-phase 4b: Sort & Index Output**
**Timing:**
- **User time:** 6.198s
- **System time:** 1.136s
- **Wall time:** 2.128s
- **CPU utilization:** 344.6%

**Statistics:**
- **Sampled:** 65,045,862 k-mers/positions (75.0% of total)

#### **Total Phase 4 (HAP2):**
- **User time:** 7.501s
- **System time:** 1.241s
- **Wall time:** 2.409s
- **CPU utilization:** 362.9%

---

### **Phase 5: Adaptive Seed Merge**
**Status:** ✅ Executed

**What it does:**
- Merges k-mer tables from both genomes to find matching seeds
- Uses "adaptamer" algorithm (Martin Frith's method)
- Finds k-mer matches between G1 (HAP1) and G2 (HAP2) genomes
- Outputs seed position pairs for alignment candidates
- Performs cache-coherent merge of sorted k-mer tables
- Uses 8 threads in parallel

**Process:**
1. Divides k-mer tables into thread partitions
2. Performs parallel merge of sorted k-mer tables
3. Identifies matching k-mers between genomes
4. Outputs seed pairs to temporary files
5. Shows progress 0-100%

**Timing:**
- **User time:** 4.103s
- **System time:** 26.291s (high I/O - reading/writing temp files)
- **Wall time:** 8.296s
- **CPU utilization:** 366.4%

**Output statistics:**
- **Total seeds:** 51,082,720
- **Average seed length:** 28.5 bp
- **Seeds per genome position:** 0.6

**Note:** The high system time (26.3s) indicates significant I/O operations reading/writing temporary files during the merge process.

---

### **Phase 6: Seed Sort and Alignment Search**
**Status:** ✅ Executed

**What it does:**
- Sorts seeds by genomic coordinates
- Groups seeds into alignment candidates
- Performs local alignment using wave-based algorithm
- Filters alignments by minimum length (default 85bp)
- Chains overlapping alignments
- Removes redundant alignments
- Processes 16 parts in parallel

**Process for each part:**
1. **Loading:** Loads seeds for the part
2. **Sorting:** Sorts seeds by diagonal/anti-diagonal and position
3. **Searching:** Performs alignment search for each seed cluster

**Timing:**
- **User time:** 1:36.449s (96.4 seconds)
- **System time:** 1.263s
- **Wall time:** 19.574s
- **CPU utilization:** 499.2% (using ~5 cores efficiently)

**Output statistics:**
- **Total hits over 85bp:** 338,700
- **Total alignments:** 376,250
- **Non-redundant alignments:** 323,569
- **Average alignment length:** 1,953 bp

**Why it takes longest:** This phase performs the actual sequence alignment, which is computationally intensive. Uses 8 threads efficiently (~500% CPU utilization).

---

### **Phase 7: Sorting and Merging Alignments**
**Status:** ✅ Executed (included in Phase 6 timing)

**What it does:**
- Sorts all alignments by:
  1. Source1 contig number
  2. Source2 contig number  
  3. Start coordinate in source1
- Merges alignments from multiple threads
- Writes final output in ONEcode ALN format

**Timing:**
- Included in Phase 6 timing (no separate measurement)

**Output file:**
- `H1vH2.1aln` - Binary alignment file (~19MB)

---

## Complete Timing Summary

### Phase-by-Phase Breakdown

| Phase | Description | User Time | System Time | Wall Time | CPU % | % of Total Wall |
|-------|------------|-----------|-------------|-----------|-------|-----------------|
| 1 | GDB Creation (HAP1) | 0.424s | 0.027s | 0.454s | 99.4% | 1.6% |
| 2a | GIX Partitioning (HAP1) | 1.273s | 0.116s | 0.295s | 471.2% | 1.1% |
| 2b | GIX Sort & Index (HAP1) | 6.148s | 1.041s | 2.083s | 345.1% | 7.5% |
| **2 Total** | **GIX Creation (HAP1)** | **7.422s** | **1.158s** | **2.378s** | **360.8%** | **8.5%** |
| 3 | GDB Creation (HAP2) | 0.425s | 0.026s | 0.454s | 99.5% | 1.6% |
| 4a | GIX Partitioning (HAP2) | 1.303s | 0.104s | 0.280s | 502.6% | 1.0% |
| 4b | GIX Sort & Index (HAP2) | 6.198s | 1.136s | 2.128s | 344.6% | 7.6% |
| **4 Total** | **GIX Creation (HAP2)** | **7.501s** | **1.241s** | **2.409s** | **362.9%** | **8.6%** |
| 5 | Adaptive Seed Merge | 4.103s | 26.291s | 8.296s | 366.4% | 29.8% |
| 6 | Seed Sort & Alignment | 96.449s | 1.263s | 19.574s | 499.2% | 70.2% |
| 7 | Sort & Merge | (included) | (included) | (included) | - | - |
| **TOTAL** | **Complete Run** | **116.34s** | **30.05s** | **33.65s** | **435%** | **100%** |

### Key Observations

1. **GDB Creation** (Phases 1 & 3): Very fast (~0.45s each), single-threaded
2. **GIX Creation** (Phases 2 & 4): Moderate (~2.4s each), parallelized (~360% CPU)
3. **Adaptive Seed Merge** (Phase 5): Moderate (~8.3s), high I/O (26s system time)
4. **Alignment Search** (Phase 6): Longest (~19.6s), highly parallelized (~500% CPU)

### Time Distribution

- **Indexing (GDB+GIX):** ~5.7s (17% of total)
- **Seed Finding:** ~8.3s (25% of total)
- **Alignment:** ~19.6s (58% of total)

---

## System Resource Usage (from /usr/bin/time)

**Overall Statistics:**
- **User time:** 116.34 seconds
- **System time:** 30.05 seconds
- **Wall clock time:** 33.65 seconds
- **CPU utilization:** 435% (using ~4.4 cores effectively)
- **Peak memory:** 668,716 KB (~653 MB)
- **Page faults:** 651,555 minor, 0 major
- **File system I/O:** 
  - Inputs: 3,489,744 blocks
  - Outputs: 5,311,752 blocks

---

## Output File: H1vH2.1aln

**Format:** ONEcode binary alignment format

**Contents:**
- 323,569 non-redundant alignments
- Each alignment encoded using trace points (space-efficient)
- File size: ~19MB

**Conversion utilities:**
- `ALNtoPAF H1vH2` - Convert to PAF format
- `ALNtoPSL H1vH2` - Convert to PSL format
- `ALNshow H1vH2` - Display alignments

---

## Understanding Resource Metrics

**Format:** `user_time system_time wall_time CPU_utilization%`

- **user_time (u):** CPU time spent in user mode (actual computation)
- **system_time (s):** CPU time spent in system/kernel mode (I/O, system calls)
- **wall_time (w):** Real elapsed time (wall clock time)
- **CPU_utilization%:** (user_time + system_time) / wall_time × 100

**Example:** `1:36.449u  1.263s  19.574w  499.2%`
- User: 1 minute 36.449 seconds
- System: 1.263 seconds  
- Wall: 19.574 seconds
- Utilization: 499.2% (using ~5 CPU cores effectively)

---

## Key Insights

1. **GDB/GIX Reuse:** When `-k` flag is used, GDB and GIX files are kept and reused in subsequent runs, saving ~5.7 seconds (17% of total time).

2. **Parallelization:** FastGA efficiently uses multiple threads (default 8), achieving 360-500% CPU utilization across phases.

3. **Most Time-Consuming Phase:** Seed sort and alignment search takes ~58% of total time, as it performs the actual sequence alignment.

4. **I/O Intensive:** Adaptive seed merge phase has high system time (26s) due to reading/writing large temporary files.

5. **Space Efficiency:** The ONEcode ALN format is highly compressed (~19MB for 323K alignments vs potentially GBs in text formats).

6. **Memory Usage:** Peak memory usage is ~653MB, which is reasonable for whole genome alignment.

7. **Adaptive Seeds:** Uses adaptamer algorithm which finds seeds adaptively based on local sequence context, making it more sensitive than fixed k-mer approaches.

---

## Comparison: With vs Without Intermediate Files

| Scenario | Wall Time | Time Saved |
|----------|-----------|------------|
| **Full run (from scratch)** | 33.65s | - |
| **With existing GDB/GIX** | ~27s | ~6.7s (20%) |

When GDB and GIX files already exist (from previous run with `-k` flag), FastGA skips Phases 1-4, saving approximately 20% of total runtime.

---

## References

- FastGA codebase: `/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/`
- Main source: `FastGA.c` - Contains phase orchestration
- GDB creation: `FAtoGDB.c` - Phase 1 & 3 implementation
- GIX creation: `GIXmake.c` - Phase 2 & 4 implementation
- Alignment: `align.c`, `ONEaln.c` - Phases 5-7 implementation
- Full run log: `EXAMPLE/fastga_full_run.log`


