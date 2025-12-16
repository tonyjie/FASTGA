# FastGA 性能优化策略：方案1（Exact-match） vs 方案2（Metric-match）

## 背景
你的总体目标：**优化 FastGA 性能**，同时用 LLM agent 保证功能正确性（functional correctness）。

FastGA 本身包含两类可动点：
- **工程实现层**（并行、内存、I/O、数据结构、调度）——通常可以做到输出等价。
- **启发式/算法决策层**（seed/chaining/alignment 过滤阈值与策略）——会改变输出集合与统计。

---

## 方案1：强制 Exact Match（推荐先做）

### 定义
- 不要求 `.1aln` 二进制文件 byte-level 相同（多线程/写入顺序/元数据可能导致不同）。
- 以 **语义等价（alignment set 等价）** 为标准：
  - `ALNtoPAF` 输出排序后 hash 完全一致。
  - （更强）`ALNtoPAF -x`（带 CIGAR）输出排序后 hash 完全一致。

### 为什么适合 FastGA
- 我们已验证：对 GRCh38 vs CHM13，T=1/8/32 的 `.1aln` **file hash 不同**，但：
  - `ALNtoPAF | sort | sha256sum` **一致**
  - `ALNtoPAF -x | sort | sha256sum` **一致**
- 这给了你一个非常强、可自动化的回归标准，适合作为 LLM agent 的“硬门槛”。

### 可做的优化类型（不改变输出）
- 线程/任务划分与负载均衡（只改变调度，不改变决策）
- 内存分配优化（arena/pool，减少 malloc/free）
- I/O 合并与缓冲（减少 syscall、减少小文件写入，但保持内容等价）
- cache-friendly layout、减少拷贝
- 确定性增强（固定 tie-breaker / stable merge），使输出顺序更稳定（可选）

### 风险
- 少量工程优化也可能不小心改变对齐边界/过滤逻辑 → 通过 exact-match gate 及时发现。

---

## 方案2：不强制 Exact Match，用 Metrics 判定

### 定义
- 允许输出集合变化，只要总体指标达标，例如：
  - coverage A/B、总 aligned bp、alignment 数量
  - divergence/identity 分布
  - 与 baseline 或其它工具（minimap2/lastz）对照的覆盖率

论文中也用 coverage/对齐总量等指标对比不同工具（见 FastGA paper）。

### 适用场景
- 你愿意接受结果集合变化，以换取更大性能收益。
- 你能维护一套数据集与指标体系（回归集），并接受“tradeoff”。

### 可动的启发式/算法点
**用户参数 knob（最常见）**：
- `-f`（seed 频率 cutoff）、`-S`（symmetric seeding）
- `-s`（chain break）、`-c`（min chain coverage）
- `-l`（min alignment length）、`-i`（min identity）、`-M`（softmask）

**更深的算法级近似（风险更高）**：
- seed 的截断/抽样策略（不仅仅是 `-f` hard cutoff）
- 更 aggressive 的 chain pruning（early-exit）
- 两阶段 alignment（cheap filter + full wave alignment）
- 修改 band/bucket 规则（影响对大 indel/漂移的敏感度）
- 修改 non-redundant 的去冗余定义

### 风险
- 很容易“变快但变差”，需要更强的评估体系与数据覆盖。

---

## 推荐路线图

1. **先做方案1**：把 `sorted PAF (+ optional CIGAR) hash` 作为硬门槛，进行纯工程优化。
2. **再做方案2（可选）**：在一个独立分支/独立配置 profile 里做启发式/算法 tradeoff，建立 metrics report。

---

## 建议的方案1回归标准（最强）

- 基准：
  - `ALNtoPAF -T<N> baseline.1aln | sort | sha256sum`
  - `ALNtoPAF -x -T<N> baseline.1aln | sort | sha256sum`
- 优化后：同样两条 hash 必须一致。

