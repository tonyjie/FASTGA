# FastGA 算法高级概述

本文档从高层次解释 FastGA 的整体思路，结合论文、代码和日志输出进行讲解。

---

## 1. FastGA 核心思想

根据论文 **"FastGA: Fast Genome Alignment"** (Myers, Zhou, Durbin)：

> "FastGA finds alignments between two genome sequences more than an order of magnitude faster than previous methods that have comparable sensitivity."

### 1.1 三个关键创新

1. **Cache-coherent Architecture**: 只使用 cache-coherent 的 MSD radix sort 和 merge 操作
2. **Adaptive Seed (Adaptamer) Algorithm**: 在线性时间内通过 merge 两个排序的 k-mer 表找到 adaptive seed
3. **Wave-based Local Alignment**: Myers 的 adaptive wave 算法变体，可检测高达 25-30% 变异的比对

### 1.2 整体 Pipeline 概述

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                           FastGA Pipeline Overview                              │
├─────────────────────────────────────────────────────────────────────────────────┤
│                                                                                 │
│  ┌─────────────┐     ┌─────────────┐     ┌─────────────┐     ┌─────────────┐   │
│  │   FASTA     │     │    GDB      │     │    GIX      │     │    ALN      │   │
│  │  (Genome)   │ ──► │  (Database) │ ──► │   (Index)   │ ──► │ (Alignments)│   │
│  └─────────────┘     └─────────────┘     └─────────────┘     └─────────────┘   │
│        │                   │                   │                   │           │
│        │   FAtoGDB         │    GIXmake        │      FastGA       │           │
│        └───────────────────┴───────────────────┴───────────────────┘           │
│                                                                                 │
└─────────────────────────────────────────────────────────────────────────────────┘
```

---

## 2. Pipeline 各阶段详解

### 阶段 1: FASTA → GDB (FAtoGDB)

**日志输出对应**:
```
Creating genome data base (GDB) GCF_000001405.40_GRCh38.p14_genomic.1gdb in directory ...
Total Resources:  6.606u  2.187s  8.794w  100.0%  0MB
```

**代码位置**: `FAtoGDB.c`

**功能**:
- 将 FASTA 文件转换为 Genome Database (GDB) 格式
- 创建 2-bit 压缩的序列存储（`.bps` 隐藏文件）
- 提取 scaffold/contig 信息、名称和大小
- 支持随机访问，I/O 效率比 FASTA 提高 4 倍

**输出文件**:
- `*.1gdb` - ONEcode 二进制元数据文件
- `.*.bps` - 隐藏的 2-bit 压缩序列文件

**论文描述**:
> "The FASTA files are converted to **genome databases** with extension **.1gdb** that are a ONEcode binary file and associated hidden file containing the ASCII DNA sequences in 2-bit compressed form."

---

### 阶段 2: GDB → GIX (GIXmake)

**日志输出对应**:
```
Creating genome genome index (GIX) GCF_000001405.40_GRCh38.p14_genomic.gix in directory ...

Partitioning K-mers via pos-lists into 16 parts
Resources for phase:  49.565u  1.539s  7.880w  648.5%

Starting sort & index output of each part
  Sorting part 1  Compressing part 1  Outputing part 1  ...
  Done

Sampled:    2316875829 ( 73.9%) kmers/positions
Resources for phase:  3:45.285u  16.744s  47.762w  506.7%
Total Resources:  4:34.850u  18.284s  55.643w  526.8%  2MB
```

**代码位置**: `GIXmake.c`

```1:10:GIXmake.c
/*******************************************************************************************
 *
 *  Produce a k-mer index of genome's contigs suitable for finding adaptamer seed matches.
 *    As such only k-mers whose count is less than the adaptamer frequency cutoff are in
 *    the index.
 *
 *  Author  :  Gene Myers
 ...
```

**功能**:
- 创建 Genome Index (GIX)，本质上是一个截断的后缀数组
- 提取所有以 (12,8) syncmer 开头的 k-mers (默认 k=40)
- 按字典顺序排序 k-mers
- 每个 k-mer 关联其基因组位置

**子阶段**:

| 子阶段 | 描述 | 日志对应 |
|--------|------|----------|
| 2a | K-mer 分区 | `Partitioning K-mers via pos-lists into N parts` |
| 2b | 排序、压缩、输出 | `Sorting part X`, `Compressing part X`, `Outputing part X` |

**"Parts" 的含义**:
- **Parts = NTHREADS × 2** (对于 8 线程，parts = 16)
- 将 k-mer 数据分成多个部分以并行处理
- 每个 part 包含一定范围的 k-mer 前缀

**输出文件**:
- `*.gix` - 主索引文件（代理文件）
- `.*.ktab.N` - 隐藏的 k-mer 表文件（N 个）
- `.*.post.N` - 隐藏的位置列表文件（N 个）

**论文描述**:
> "A **genome index** with extension **.gix** is then built for each genome that is basically a truncated suffix array. One of the things that makes FastGA fast is that it **compares these two indices against each other directly** rather than looking up sequences of one genome in the index of the other."

---

### 阶段 3: Adaptive Seed Merge

**日志输出对应**:
```
Using 8 threads

Starting adaptive seed merge for G1
  Completed   0%  Completed   1% ... Completed 100%

Total seeds = 1298273853, ave. len = 39.0, seeds per genome position = 0.4
Resources for phase:  2:01.381u  13.447s  18.947w  711.6%
```

**代码位置**: `FastGA.c` 第 2280-2500 行

```2280:2290:FastGA.c
static void adaptamer_merge(Kmer_Stream *T1, Kmer_Stream *T2,
                            Post_List *P1,   Post_List *P2, int64 g1len)
{ SP         parm[NTHREADS];
#ifndef DEBUG_MERGE
  pthread_t  threads[NTHREADS];
#endif
  uint8     *cache;
  int64      nhits, tseed;
  int        i;
```

主程序调用（第 5121-5124 行）:
```5121:5124:FastGA.c
    if (SELF)
      self_adaptamer_merge(T1,P1,gdb1->seqtot);
    else
      adaptamer_merge(T1,T2,P1,P2,gdb1->seqtot);
```

**功能**:
- **核心算法**: 线性时间 merge 两个已排序的 k-mer 表
- 找到两个基因组之间匹配的 k-mers
- 输出 seed position pairs (种子位置对)

**Adaptamer (Adaptive Seed) 概念**:

根据论文：
> "An **adaptive seed** at a given position p of source1, is the **longest string beginning at that position that is also somewhere in source2**. If the number of occurrences of an adaptamer in source2 exceeds the frequency cutoff (e.g., 10), then the adaptamer is deemed *repetitive* and is not considered."

**算法思路**:
1. 两个基因组的 k-mer 表已按字典顺序排序
2. 使用 **线性 merge** 同时遍历两个表
3. 当找到匹配的 k-mer 时，检查其出现频率
4. 如果频率 ≤ cutoff，输出所有 (position_in_G1, position_in_G2) 对

**日志输出解释**:
- **Total seeds = 1,298,273,853**: 找到的总种子数（位置对）
- **ave. len = 39.0**: 平均种子长度（bp）
- **seeds per genome position = 0.4**: 每个基因组位置的平均种子数

**论文公式**:
种子匹配条件：对于位置 p 的 adaptamer，如果在 G2 中的出现次数 ≤ frequency cutoff，则输出所有 (p, q) 对，其中 q 是 G2 中的匹配位置。

---

### 阶段 4: Seed Sort and Alignment Search

**日志输出对应**:
```
Starting seed sort and alignment search, 18 parts
  Loading seeds for part 1  Sorting seeds for part 1  Searching seeds for part 1
  Loading seeds for part 2  Sorting seeds for part 2  Searching seeds for part 2
  ...
  Loading seeds for part 18  Sorting seeds for part 18  Searching seeds for part 18
  Done

Total hits over 85bp = 843675, 791660 aln's, 518037 non-redundant aln's of ave len 28415
```

**代码位置**: `FastGA.c` 第 4134-4410 行

```4134:4155:FastGA.c
static void pair_sort_search(GDB *gdb1, GDB *gdb2)
{ uint8 *sarray;
  int    swide;
  int64  nels;

  RP     rarm[NTHREADS];
  TP     tarm[NTHREADS];
  pthread_t threads[NTHREADS];
  int64    *panel;
  Range     range[NTHREADS];

  IOBuffer *unit[2], *nu;
  int       nused;
  int       i, p, j, u;

  if (VERBOSE)
    { fprintf(stderr,"\n  Starting seed sort and alignment search, %d parts\n",2*NPARTS);
      fflush(stderr);
    }
```

主程序调用（第 5171 行）:
```5171:5171:FastGA.c
    pair_sort_search(gdb1,gdb2);
```

**功能**:
这是 FastGA 中**最耗时的阶段**，包含三个子步骤：

| 子步骤 | 描述 | 日志对应 |
|--------|------|----------|
| Loading | 加载种子数据 | `Loading seeds for part X` |
| Sorting | 按对角线/反对角线排序 | `Sorting seeds for part X` |
| Searching | 执行实际比对 | `Searching seeds for part X` |

**"Parts" 的含义**:
- **Parts = 2 × NPARTS** (对于 8 线程，NPARTS=9，所以 parts=18)
- 每个 part 对应一组 contigs
- 分两组：正向（N_Units）和反向互补（C_Units）

**算法步骤**:

1. **Seed Chaining**: 
   - 将种子按对角线（diagonal）排序
   - 找到在 128bp 宽的对角线带内的种子链
   - 相邻种子间距 < 1000bp

   论文描述:
   > "FastGA then searches for runs or chains of adaptamer seed hits that (a) all lie within a diagonal band of width 128, (b) the spacing between every pair of consecutive seeds is less than -s(1000)"

2. **Alignment Verification**:
   - 使用 Myers 的 wave-based local alignment 算法
   - 可以检测高达 25-30% 的序列差异
   - 只保留长度 ≥ 85bp（-c 参数）的比对

   论文描述:
   > "a variant of the Myers adaptive wave algorithm to find alignments around a chain of seed hits that detects alignments with up to 25-30% variation"

**日志输出解释**:
- **Total hits over 85bp = 843,675**: 长度超过 85bp 的比对命中
- **791,660 aln's**: 总比对数
- **518,037 non-redundant aln's**: 去冗余后的比对数
- **ave len 28415**: 平均比对长度（bp）

---

### 阶段 5: Sorting and Merging Alignments

**日志输出对应**:
```
Sorting and merging alignments
Resources for phase:  42:26.326u  2:10.026s  14:30.410w  307.5%
```

**代码位置**: `FastGA.c` 第 4391-4410 行

```4391:4410:FastGA.c
  if (VERBOSE)
    { fprintf(stderr,"\n  Sorting and merging alignments\n");
      fflush(stderr);
    }
  if (LOG_FILE)
    fprintf(LOG_FILE,"\n  Sorting and merging alignments\n");

#ifdef DEBUG_LASORT
  for (p = 0; p < NTHREADS; p++)
    la_sort(tarm+p);
#else
  for (p = 1; p < NTHREADS; p++)
    pthread_create(threads+p,NULL,la_sort,tarm+p);
  la_sort(tarm);
  for (p = 1; p < NTHREADS; p++)
    pthread_join(threads[p],NULL);
#endif

  if (la_merge(tarm))
    Clean_Exit(1);
```

**功能**:
- 排序所有比对（按 contig1 号、contig2 号、起始位置）
- 合并多线程产生的比对结果
- 输出最终的 ONEcode ALN 格式文件

**输出文件**:
- `*.1aln` - 使用 trace point 编码的高效比对文件

论文描述:
> "FastGA records all the alignments it finds in a ONEcode binary file we refer to here as a ALN-formatted file with extension **.1aln** that uses a very space efficient **trace point encoding** of each alignment."

---

## 3. 代码架构总览

```
FastGA.c
├── main() [第 4433 行]
│   ├── 解析命令行参数
│   ├── 调用 FAtoGDB (如需要)
│   ├── 调用 GIXmake (如需要)
│   ├── Open_Post_List(), Open_Kmer_Stream() - 加载索引
│   ├── adaptamer_merge() 或 self_adaptamer_merge() [第 5121-5124 行]
│   │   └── 线性 merge 两个 k-mer 表，输出种子对
│   ├── pair_sort_search() [第 5171 行]
│   │   ├── 加载种子
│   │   ├── 按对角线排序
│   │   ├── Seed chaining
│   │   └── Wave-based alignment
│   └── la_merge() - 合并输出比对结果
│
├── FAtoGDB.c - FASTA → GDB 转换
├── GIXmake.c - GDB → GIX 索引构建
├── align.c   - Wave-based 比对算法
└── alncode.c - ONEcode ALN 编码/解码
```

---

## 4. 关键数据结构

### 4.1 K-mer Table (Kmer_Stream)
```c
// 存储排序的 k-mers 和对应位置
Kmer_Stream *T1, *T2;  // 两个基因组的 k-mer 流
```

### 4.2 Position List (Post_List)
```c
// 存储每个 k-mer 在基因组中的位置
Post_List *P1, *P2;
```

### 4.3 IOBuffer
```c
// 用于种子对的 I/O 缓冲
typedef struct {
  uint8 *bufr;   // 缓冲区
  uint8 *btop;   // 顶部指针
  uint8 *bend;   // 结束指针
  int64 *buck;   // 桶计数
  int    file;   // 文件描述符
  int    inum;   // 索引号
} IOBuffer;

static IOBuffer *N_Units;  // 正向种子对
static IOBuffer *C_Units;  // 反向互补种子对
```

---

## 5. 日志输出 vs 代码对应表

| 日志输出 | 代码位置 | 函数/阶段 |
|----------|----------|-----------|
| `Creating genome data base (GDB)` | FAtoGDB.c | FAtoGDB 主函数 |
| `Creating genome index (GIX)` | GIXmake.c | GIXmake 主函数 |
| `Partitioning K-mers via pos-lists` | GIXmake.c | 分区阶段 |
| `Sorting part N` | GIXmake.c | 排序阶段 |
| `Compressing part N` | GIXmake.c | 压缩阶段 |
| `Outputing part N` | GIXmake.c | 输出阶段 |
| `Sampled: N kmers/positions` | GIXmake.c | 采样统计 |
| `Using N threads` | FastGA.c:5043 | 线程初始化 |
| `Starting adaptive seed merge` | FastGA.c:2364/2505 | adaptamer_merge() |
| `Completed N%` | FastGA.c | merge 进度 |
| `Total seeds = N` | FastGA.c:2451 | 种子统计 |
| `Starting seed sort and alignment search` | FastGA.c:4150 | pair_sort_search() |
| `Loading seeds for part N` | FastGA.c | 加载阶段 |
| `Sorting seeds for part N` | FastGA.c | 排序阶段 |
| `Searching seeds for part N` | FastGA.c | 搜索阶段 |
| `Total hits over Nbp` | FastGA.c:4376 | 比对统计 |
| `Sorting and merging alignments` | FastGA.c:4391 | la_sort/la_merge |

---

## 6. 性能特点总结

根据论文和代码分析：

1. **Cache-coherent 设计**: 所有操作都是 MSD radix sort 和 linear merge，最大化缓存利用率

2. **线性时间 seed finding**: 通过同时遍历两个排序的 k-mer 表，在 O(n+m) 时间内找到所有匹配

3. **高效的比对编码**: 使用 trace point encoding，比 CIGAR 字符串节省大量空间

4. **并行化**: 
   - Seed merge: 高并行度（CPU 利用率可达 2000%+）
   - Alignment search: 受 I/O 限制，并行度较低

---

## 7. 下一步：Deep Dive

后续可以深入探讨：
1. **Adaptamer merge 算法细节** - 如何在线性时间内完成 merge
2. **Wave-based alignment** - Myers 的 adaptive wave 算法
3. **Trace point encoding** - ONEcode ALN 格式的编码方式
4. **Parallelization strategy** - 如何分配工作给多线程

---

## 参考

- 论文: "FastGA: Fast Genome Alignment" - Gene Myers, Chenxi Zhou, Richard Durbin
- 代码仓库: https://github.com/thegenemyers/FASTGA
- ONEcode: https://github.com/thegenemyers/ONEcode




