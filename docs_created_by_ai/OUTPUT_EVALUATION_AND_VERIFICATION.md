# FastGA 输出评估与验证指南

本文档帮助你理解如何评估 FastGA 输出的质量以及验证结果的正确性。这对于代码优化后的验证至关重要。

---

## 1. 输出格式概述

FastGA 可以输出三种格式：

| 格式 | 扩展名 | 说明 | 命令选项 |
|------|--------|------|----------|
| ONEcode ALN | `.1aln` | 二进制格式，最高效 | `-1:output.1aln` |
| PAF | 标准输出 | 文本格式，可加 CIGAR | `-paf`, `-pafx`, `-pafm` |
| PSL | 标准输出 | BLAT 兼容格式 | `-psl` |

---

## 2. 核心质量指标

### 2.1 ALNshow 输出解读

运行 `ALNshow <alignment_file>` 查看比对详情：

```bash
ALNshow H1vH2 | head -20
```

输出示例：
```
H1vH2: 380,294 records
 1.01n   1.01n   [         2..     5,982] x [       483..     6,497]  ~   4.35%   (86,251,610 x 86,746,603 bps,      261 diffs,     60 trace pts)
```

**字段解释**：

| 字段 | 含义 | 示例值 |
|------|------|--------|
| `1.01n` | Contig ID + 方向 (n=normal, c=complement) | `1.01n` = scaffold 1, contig 1, normal |
| `[2..5,982]` | Query 区间 | 从位置 2 到 5982 |
| `[483..6,497]` | Target 区间 | 从位置 483 到 6497 |
| `~4.35%` | **Divergence (差异率)** | 4.35% 的位置不同 |
| `261 diffs` | **差异数量** | 261 个 substitutions + indels |
| `60 trace pts` | Trace points 数量 | 用于编码比对的点数 |

### 2.2 差异率计算公式

从代码 `ALNshow.c:607-608`：

```c
printf("  ~  %5.2f%% ",(200.*ovl->path.diffs) /
       ((ovl->path.aepos - ovl->path.abpos) + (ovl->path.bepos - ovl->path.bbpos)) );
```

**公式**:
```
Divergence = (2 × diffs) / (query_length + target_length) × 100%
```

**论文解释**:
> "The divergence is defined as the edit distance divided by the average length of the two aligned segments."

### 2.3 PAF 输出质量指标

运行 `ALNtoPAF <alignment_file>` 查看 PAF 格式：

```bash
ALNtoPAF H1vH2 | head -5
```

输出示例：
```
>SUPER_2  86251610  2  5982  +  >SUPER_2  86746603  483  6497  5866  5997  255  dv:f:0.0191  df:i:261
```

**关键 SAM tags**：

| Tag | 含义 | 示例 |
|-----|------|------|
| `dv:f:0.0191` | **Divergence** (差异率) | 1.91% |
| `df:i:261` | **Diffs** (差异数量) | 261 个差异 |

**PAF 标准字段**：

| 列 | 含义 | 示例 |
|---|------|------|
| 10 | Number of matches | 5866 |
| 11 | Alignment block length | 5997 |
| 12 | Mapping quality | 255 |

---

## 3. 验证正确性的方法

### 3.1 方法 1: 比较比对数量和统计

**基准测试**：保存原始运行的统计信息

```bash
# 原始运行
FastGA -vk -1:original.1aln genome1.fna genome2.fna 2>&1 | tee original.log

# 优化后运行
FastGA -vk -1:optimized.1aln genome1.fna genome2.fna 2>&1 | tee optimized.log
```

**比较关键统计**：

```bash
# 提取关键统计
grep "Total seeds" original.log optimized.log
grep "non-redundant aln's" original.log optimized.log
grep "ave len" original.log optimized.log
```

**期望结果**：
- `Total seeds` 应该相同
- `non-redundant aln's` 应该相同（或非常接近）
- `ave len` 应该相似

### 3.2 方法 2: 比较 PAF 输出

```bash
# 生成 PAF 并排序
ALNtoPAF original.1aln | sort > original.paf
ALNtoPAF optimized.1aln | sort > optimized.paf

# 比较
diff original.paf optimized.paf | head -50

# 或者比较关键列
cut -f1-12 original.paf | sort > original_key.paf
cut -f1-12 optimized.paf | sort > optimized_key.paf
diff original_key.paf optimized_key.paf
```

### 3.3 方法 3: 验证 CIGAR 字符串

```bash
# 生成带 CIGAR 的 PAF
ALNtoPAF -x original.1aln | sort > original_cigar.paf
ALNtoPAF -x optimized.1aln | sort > optimized_cigar.paf

# 比较 CIGAR 字符串
diff original_cigar.paf optimized_cigar.paf
```

### 3.4 方法 4: 统计汇总比较

创建验证脚本 `verify_fastga.sh`:

```bash
#!/bin/bash
# verify_fastga.sh - 验证 FastGA 输出正确性

ORIGINAL=$1
OPTIMIZED=$2

echo "=== 比对数量 ==="
echo -n "Original:  "
ALNshow $ORIGINAL 2>/dev/null | head -1
echo -n "Optimized: "
ALNshow $OPTIMIZED 2>/dev/null | head -1

echo ""
echo "=== 差异率分布 ==="
echo "Original divergence stats:"
ALNtoPAF $ORIGINAL 2>/dev/null | awk -F'\t' '{for(i=1;i<=NF;i++)if($i~/dv:f:/){split($i,a,":"); print a[3]}}' | \
  awk '{sum+=$1; count++} END {printf "  Mean: %.4f, Count: %d\n", sum/count, count}'

echo "Optimized divergence stats:"
ALNtoPAF $OPTIMIZED 2>/dev/null | awk -F'\t' '{for(i=1;i<=NF;i++)if($i~/dv:f:/){split($i,a,":"); print a[3]}}' | \
  awk '{sum+=$1; count++} END {printf "  Mean: %.4f, Count: %d\n", sum/count, count}'

echo ""
echo "=== Coverage 分布 ==="
echo "Original:"
ALNtoPAF $ORIGINAL 2>/dev/null | awk '{sum+=$11} END {printf "  Total aligned bases: %d\n", sum}'

echo "Optimized:"
ALNtoPAF $OPTIMIZED 2>/dev/null | awk '{sum+=$11} END {printf "  Total aligned bases: %d\n", sum}'
```

运行：
```bash
chmod +x verify_fastga.sh
./verify_fastga.sh original.1aln optimized.1aln
```

---

## 4. 论文中的质量评估指标

### 4.1 Coverage（覆盖度）

论文 Table 2 中使用的指标：

> **Coverage A (Mb)**: 基因组 A 被比对覆盖的碱基数
> **Coverage B (Mb)**: 基因组 B 被比对覆盖的碱基数

计算方法：
```bash
# Coverage A (query)
ALNtoPAF alignment.1aln | awk '{sum+=($4-$3)} END {print sum/1000000 " Mb"}'

# Coverage B (target)
ALNtoPAF alignment.1aln | awk '{sum+=($9-$8)} END {print sum/1000000 " Mb"}'
```

### 4.2 与其他工具比较

论文比较了 FastGA 与以下工具：
- **minimap2**
- **LASTZ**
- **GSAlign**
- **AnchorWave**

**关键发现**（Table 2）：
- FastGA 的 coverage 与其他工具相当
- FastGA 速度最快（通常快 10-100 倍）
- FastGA 内存使用最低

---

## 5. 代码中的质量阈值

### 5.1 比对过滤参数

从 `FastGA.c` 中的默认值：

```c
// FastGA.c:4448-4452
FREQ = 10;           // -f: 最大种子频率（过滤重复 k-mers）
CHAIN_BREAK = 2000;  // -s: 链断裂阈值
CHAIN_MIN   =  170;  // -c: 最小链覆盖度 (85bp × 2)
ALIGN_MIN   =  100;  // -l: 最小比对长度
ALIGN_RATE  = .3;    // -i: 最大差异率 (1 - 0.7 = 30%)
```

### 5.2 比对验证代码

```c
// FastGA.c:3264-3265 和 3283-3284
if (rlen >= alnMin && alnRate*rlen >= path->diffs)
  { /* 接受比对 */ }

if (rlen >= ALIGN_MIN && ALIGN_RATE*rlen >= path->diffs)
  { /* 最终验证 */ }
```

**验证条件**：
1. `rlen >= ALIGN_MIN`：比对长度 ≥ 100bp（默认）
2. `ALIGN_RATE × rlen >= diffs`：差异率 ≤ 30%（默认）

---

## 6. 完整验证工作流

### 6.1 优化前的基准测试

```bash
# 1. 运行原始版本并保存所有输出
./FastGA_original -vk -1:baseline.1aln genome1.fna genome2.fna 2>&1 | tee baseline.log

# 2. 生成各种格式的输出
ALNshow baseline.1aln > baseline_alnshow.txt
ALNtoPAF baseline.1aln > baseline.paf
ALNtoPAF -x baseline.1aln > baseline_cigar.paf
ALNtoPSL baseline.1aln > baseline.psl

# 3. 提取统计信息
grep -E "Total seeds|non-redundant|ave len|Total Resources" baseline.log > baseline_stats.txt
```

### 6.2 优化后的验证

```bash
# 1. 运行优化版本
./FastGA_optimized -vk -1:test.1aln genome1.fna genome2.fna 2>&1 | tee test.log

# 2. 生成输出
ALNshow test.1aln > test_alnshow.txt
ALNtoPAF test.1aln > test.paf
ALNtoPAF -x test.1aln > test_cigar.paf

# 3. 比较
echo "=== ALNshow header ==="
head -5 baseline_alnshow.txt test_alnshow.txt

echo "=== PAF comparison ==="
diff <(sort baseline.paf) <(sort test.paf) | head -20

echo "=== Statistics ==="
diff baseline_stats.txt <(grep -E "Total seeds|non-redundant|ave len|Total Resources" test.log)
```

### 6.3 验证通过标准

| 指标 | 要求 | 说明 |
|------|------|------|
| 比对数量 | 完全相同 | `non-redundant aln's` 必须相同 |
| 种子数量 | 完全相同 | `Total seeds` 必须相同 |
| PAF 内容 | 完全相同 | 排序后的 PAF 文件应完全一致 |
| 差异率分布 | 完全相同 | 每个比对的 `dv:f:` 值应相同 |
| 覆盖度 | 完全相同 | Coverage A/B 应相同 |

---

## 7. 常见问题

### 7.1 比对数量不同

**可能原因**：
- 随机因素（浮点精度、线程调度）
- 种子频率阈值变化
- 链覆盖度计算差异

**检查方法**：
```bash
# 检查种子统计
grep "Total seeds" baseline.log test.log
grep "ave. len" baseline.log test.log
```

### 7.2 差异率略有不同

**可能原因**：
- 比对边界处理差异
- trace point 计算精度

**容忍度**：
- 差异率差异 < 0.01% 通常可接受
- 如果差异 > 0.1%，需要调查原因

### 7.3 性能 vs 准确性权衡

如果需要牺牲一些准确性来提升性能：
- 可以调整 `-f` 参数（种子频率阈值）
- 可以调整 `-c` 参数（最小链覆盖度）
- 但必须记录并验证影响范围

---

## 8. 推荐验证命令

### 快速验证

```bash
# 比对数量检查
ALNshow baseline.1aln 2>/dev/null | head -1
ALNshow test.1aln 2>/dev/null | head -1

# 应该显示相同的 records 数量
```

### 完整验证

```bash
# 完整 PAF 比较
ALNtoPAF baseline.1aln | sort | md5sum
ALNtoPAF test.1aln | sort | md5sum

# MD5 应该完全相同
```

### 深度验证

```bash
# 随机抽样比较实际比对内容
ALNshow -a baseline.1aln | head -1000 > baseline_sample.aln
ALNshow -a test.1aln | head -1000 > test_sample.aln
diff baseline_sample.aln test_sample.aln
```

---

## 总结

验证 FastGA 输出正确性的关键步骤：

1. **比较统计信息**：种子数、比对数、平均长度
2. **比较 PAF 输出**：排序后应完全一致
3. **比较 CIGAR 字符串**：确保比对细节正确
4. **检查覆盖度**：确保没有遗漏重要区域

优化代码后，必须确保以上所有指标与原始版本一致，才能确认优化没有影响结果正确性。

