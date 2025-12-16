# FastGA 性能优化与正确性验证：方案1 vs 方案2

本文件记录我们关于 **Plan 1（exact match）** 与 **Plan 2（metrics match）** 的讨论，并给出建议的执行路线。

---

## 背景与目标

你的 overall goal：
- **优化 FastGA 性能**（例如提高并行利用率、减少 I/O、改进数据结构/内存布局）
- 同时通过可自动化的验证手段，确保 **functional correctness**

我们已经在 GRCh38 vs CHM13 的实验中验证了一个关键事实：
- 不同线程数得到的 `.1aln` **byte-level 不一致**（sha256 不同）
- 但 `ALNtoPAF` 转换后的输出（尤其是 `-x` 带 CIGAR）在 **排序后（order-independent）hash 完全一致**
  - 说明 alignment set 与 base-level alignment encoding（CIGAR X/=）在不同线程下保持一致

这为 Plan 1 提供了很强的回归判据。

---

## Plan 1：强制 exact match（推荐作为第一阶段）

### 定义
- 优化 FastGA 的实现/工程性能，但要求最终输出与 baseline **完全等价**。

### “完全等价”的判据（建议）
- **推荐强判据**：
  - `ALNtoPAF -x <.1aln> | sort | sha256sum` 必须一致
  - 这同时固定了：
    - alignment 集合
    - 精确 base-level 对齐（CIGAR X/=）
  - 由于 `sort` 去除了顺序差异影响，因此适合作为多线程/并行实现的等价判据。

- **补充判据（快速 sanity check）**：
  - `ALNshow <.1aln> | head -1` 的 records 数一致

### 适合做的优化类型（理论上不影响输出）
- 内存分配/内存池、减少 malloc/free
- I/O 批处理、减少 syscall、减少临时文件读写开销
- 更好的线程调度、负载均衡、减少线程间争用
- 更 cache-friendly 的数据结构布局
- 更确定性的 merge / tie-breaker（避免非确定性）

### 风险与注意
- `.1aln` 二进制文件本身很可能无法 byte-level 一致（内部元数据/写入顺序等）。
- 所以 **不要用 `.1aln` 的 sha256** 作为 Plan 1 的主判据，而应使用 **转换为标准文本表示（PAF+CIGAR）后再比较**。

---

## Plan 2：允许不 exact match，用 coverage/sensitivity 等指标判定

### 定义
- 允许 LLM/工程改动调整 FastGA 的 heuristics 或近似策略，只要最终质量指标达标。

### 可调的范围

#### 2.1 参数 knob（最稳健、可控）
- `-f`（seed frequency cutoff）
- `-S`（symmetric seeding）
- `-c/-s`（seed chain coverage / chain break）
- `-l/-i`（min alignment length / min identity）
- `-M`（soft-mask）

这些改变会 **直接影响输出**（alignment 数量、覆盖度、重复区域行为等）。

#### 2.2 算法级近似策略（风险更高）
- 更 aggressive 的 seed downsample（不是 hard cutoff）
- 更强的 early-pruning（chaining 阶段提前剪枝）
- 两阶段验证：cheap filter → full alignment
- 改 band/bucket 粒度或策略
- 改去冗余/merge 策略

### 质量指标（建议）
- Coverage A/B（Mb 或比例）
- 总 aligned bp
- alignment 数（total / non-redundant）
- divergence 分布（dv:f）
- （可选）与其他 aligner 的 overlap coverage 交叉验证

### 推荐实践
- Plan 2 最好作为 **Plan 1 之后的第二阶段/单独分支**，因为它需要更严格的评估体系和多数据集回归。

---

## 推荐路线

1. **先做 Plan 1**：
   - 以 `ALNtoPAF -x | sort | sha256sum` 作为 hard gate
   - 允许尽可能多的工程级优化而不改变结果

2. **再评估 Plan 2**：
   - 将 heuristics/算法变体封装为明确 profile
   - 建立 coverage/sensitivity 报告与数据集回归

---

## Plan 1 的验证脚本

我们将提供一个脚本（见 `docs_created_by_ai/verify_exact_match_plan1.sh`），用于比较 baseline 与 candidate `.1aln` 是否 exact-match（按 PAF+CIGAR 的 order-independent hash）。
