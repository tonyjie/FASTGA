# FastGA 输出一致性对比：T=1 vs T=8 vs T=32

## 对比对象

- T=1: `/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/thread_1/GRCh38_vs_CHM13.1aln`
- T=8: `/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/thread_8/GRCh38_vs_CHM13.1aln`
- T=32: `/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/thread_32/GRCh38_vs_CHM13.1aln`

## 结论（先说结论）

- **(1) `.1aln` 二进制文件是否 byte-level 完全一致？**
  - **不一致**（sha256 不同）。

- **(2) alignment 结果集合是否完全一致（不受输出顺序影响）？**
  - **一致**：对三份 `.1aln` 用 `ALNtoPAF` 转换后，做 **排序（order-independent）** 的 PAF 流 sha256，三者 **完全相同**。

这意味着：
- **线程数改变不会改变最终 alignment set**（至少对这组输入/参数而言）。
- `.1aln` 的差异更可能来自 **文件内部元数据 / 写入细节 / trace 压缩或记录顺序** 等（但转换成标准 PAF 后，alignment 内容一致）。

---

## 证据 1：ALNshow 记录数与开头摘要一致

三份文件 `ALNshow | sed -n '1,5p'`：
- 都是 `GRCh38_vs_CHM13: 518,037 records`
- 前几条摘要行内容一致

---

## 证据 2：`.1aln` 文件 sha256（byte-level 不一致）

```bash
sha256sum \
  /scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/thread_1/GRCh38_vs_CHM13.1aln \
  /scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/thread_8/GRCh38_vs_CHM13.1aln \
  /scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/thread_32/GRCh38_vs_CHM13.1aln
```

输出（你的环境里实际观测）：
- T=1: `dfd2c77d...313e3d`
- T=8: `1ad6717c...fb1f83`
- T=32: `fbc73776...3adbe8`

---

## 证据 3：order-independent 的 PAF 内容 sha256（完全一致）

关键思路：不同线程可能导致输出顺序不同，所以我们把 PAF 结果按整行排序后再 hash。

```bash
cd /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA
mkdir -p /scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/_tmp_sort

for t in 1 8 32; do
  ALN=/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/thread_${t}/GRCh38_vs_CHM13.1aln
  echo "=== sorted-PAF sha256 t=${t} ==="
  ./ALNtoPAF -T32 "$ALN" \
    | LC_ALL=C sort -T /scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/_tmp_sort \
    | sha256sum
  echo
 done
```

输出（你的环境里实际观测）：
- T=1:  `54635008...8eb1ff5`
- T=8:  `54635008...8eb1ff5`
- T=32: `54635008...8eb1ff5`

**完全相同**。

---

## 对性能优化的意义

你现在可以把 **“sorted PAF hash 相同”** 作为一个非常强的回归标准：
- 只要优化后仍能保证 `ALNtoPAF | sort | sha256sum` 不变，就基本说明 alignment 输出集合没有变化。

同时你也应该预期：
- `.1aln` 二进制文件 **不一定**能做到 byte-level 完全一致（尤其在多线程/不同 run 之间）。

如果你希望进一步确认 `.1aln` 内部的“record 顺序是否一致”，可以再做：
- `./ALNtoPAF file.1aln | sha256sum`（不排序）

但这会把“输出顺序差异”也算作不一致；对验证正确性来说通常不如 **排序后 hash** 稳健。


---

## CIGAR string（`ALNtoPAF -x`）一致性对比

### 结论
- **排序后（order-independent）hash：一致** → 说明不仅 alignment 区间一致，**连精确的 base-level 对齐（CIGAR: X/=）也一致**。
- **不排序 hash：不一致** → 输出顺序仍然会随线程数变化。

### 证据：`-x` 且排序后的 PAF sha256
- T=1 / T=8 / T=32：`7fcd0e3e78559743322f686b3aa0e397f5573c4d38e1b7e4447a82762d1a904e`

### 证据：`-x` 且不排序的 PAF sha256（顺序敏感）
- T=1：`40038b24006358488b08de4508de39ac637d1f557973f83dec5cb7a728e4a0b5`
- T=8：`4939902f63cfb4d7da57deb35d6a686f63a5886665d02df84ab4051371cdf43d`
- T=32：`5be7ed8db2b2836bc0f36fe51ab6d128fae3fcce33953dfc0117a1f6bc863af4`

### 命令

```bash
cd /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA
mkdir -p /scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/_tmp_sort

# 排序后 hash（推荐回归标准）
for t in 1 8 32; do
  ALN=/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/thread_${t}/GRCh38_vs_CHM13.1aln
  ./ALNtoPAF -x -T32 "$ALN"     | LC_ALL=C sort -T /scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/_tmp_sort     | sha256sum
 done

# 不排序 hash（顺序敏感）
for t in 1 8 32; do
  ALN=/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/thread_${t}/GRCh38_vs_CHM13.1aln
  ./ALNtoPAF -x -T32 "$ALN" | sha256sum
 done
```
