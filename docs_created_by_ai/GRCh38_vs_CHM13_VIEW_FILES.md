# GRCh38 vs CHM13: 可视化/查看文件（由 ALNshow / ALNplot 生成）

## 输入
- **ALN 文件**: `/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/thread_32/GRCh38_vs_CHM13.1aln`

## 生成的可查看文件位置
都放在：
- `/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/thread_32/view/`

当前包含：
- `ALNshow_head500.txt`
  - **内容**: `ALNshow` 的前 500 条 alignment 记录（只看摘要行）
  - **用途**: 快速感受整体输出格式、diffs、trace pts 等统计

- `ALNshow_example_alignment_w120.txt`
  - **内容**: 从 `@1.2:10002-@1.2:32028` 区间抽取的 alignment（`ALNshow -a`，宽度 120，border 20，uppercase）
  - **用途**: 直接肉眼检查一个代表性 alignment 的字符级别正确性

- `ALNplot_overview.pdf`
  - **内容**: 全局 collinearity dotplot（默认最多 100,000 条用于画图，`-L` 关闭标签）
  - **用途**: 看大尺度 synteny/重排/反向互补条带

- `ALNplot_long_id90_len100k.pdf`
  - **内容**: 过滤后的 dotplot（`-l 100000 -i 0.90 -n 100000`）
  - **用途**: 只看更长、identity 更高的对齐，图更干净

## 如何复现这些文件（命令）

> 注：`ALNplot` 内部会先写 EPS 再用 `ps2pdf/epstopdf` 转成 PDF；在本环境里生成的 PDF 可能先出现在运行目录，再移动到目标目录。

```bash
OUTDIR=/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/thread_32/view
ALN=/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13/thread_32/GRCh38_vs_CHM13.1aln
FASTGA_BIN=/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA

mkdir -p "$OUTDIR"

# 1) ALNshow：前 500 条摘要
cd "$FASTGA_BIN"
./ALNshow "$ALN" | head -500 > "$OUTDIR/ALNshow_head500.txt" || true

# 2) ALNshow：抽取一个区间的完整 alignment 显示
./ALNshow -a -w120 -b20 -U "$ALN" "@1.2:10002-@1.2:32028" > "$OUTDIR/ALNshow_example_alignment_w120.txt"

# 3) ALNplot：全局图（PDF）
./ALNplot -L -T4 -p:ALNplot_overview.pdf "$ALN" > /dev/null
mv -f ALNplot_overview.pdf "$OUTDIR/"

# 4) ALNplot：过滤后的图（PDF）
./ALNplot -L -T4 -l100000 -i0.90 -n100000 -p:ALNplot_long_id90_len100k.pdf "$ALN" > /dev/null
mv -f ALNplot_long_id90_len100k.pdf "$OUTDIR/"
```
