# TODO: upstream storage comparison + optimize-memory rebase

Context:新旧 upstream 存储画像实测一致(旧 5671357 vs 新 ddeea32,EXAMPLE 默认跑法逐字节相同,仅 T=1 有 20B 采样噪声)。现在做三件事。

## 1. 文档 + 图:新旧 upstream 存储等价
- [ ] 写绘图脚本 `benchmarks/plot_upstream_storage_comparison.py`(读两份 audit 的 monitor tsv + summary）
- [ ] 生成图 `docs/benchmark_storage_upstream/storage_old_vs_new.png`
      - Panel A: temp 磁盘占用 over time(T=8,旧 vs 新 叠加 → 重合）
      - Panel B: peak_total vs threads(旧/新，flat）
- [ ] 写文档 `docs/benchmark_storage_upstream.md`(方法、结果表、图、代码 diff 论证、结论）
- [ ] 和 `docs/benchmark_thread_scaling_upstream.md` 配套(交叉引用）

## 2. 把 optimize-memory 的优化重建到新 upstream 并做真实对比
- [ ] 在独立 worktree 新建分支 `optimize-memory-ddeea32`（基于 optimize-memory）
- [ ] `git rebase --onto ddeea32 5671357`（把 18 个 commit 重放到新 upstream）
- [ ] 解决冲突(Opt1/Opt3 vs upstream ANO 改动，预计 FastGA.c/GIXmake.c 有冲突)
      - 若冲突非平凡 → STOP，报告用户，不擅自语义合并
- [ ] `make`，正确性验证(ONEview diff vs upstream bit-exact)
- [ ] 跑 audit → `benchmarks/storage_audit_optmem_ddeea32/`
- [ ] 对比 优化版 vs 新 upstream，量化 Opt1/Opt3 到底省多少(注意 -k 会抑制 Opt1）
- [ ] NOT push/commit to remote without approval

## 3. 清理临时 worktree
- [ ] 移除 `../FASTGA-old-upstream`、`../FASTGA-main-upstream`、rebase worktree
- [ ] 保留 benchmarks/ 下的结果数据

## Review (完成 2026-07-09)

### 1. 新旧 upstream 存储等价 ✅
- 图 `docs/benchmark_storage_upstream/storage_old_vs_new.png` + 文档 `docs/benchmark_storage_upstream.md`
- 结论:旧(5671357)/新(ddeea32)存储画像逐字节一致,仅 T=1 有 20B 采样噪声;代码 diff 证明写盘核心未变
- 已与 `benchmark_thread_scaling_upstream.md` 交叉引用

### 2. optimize-memory rebase 到 ddeea32 + 真实对比 ✅
- 分支 `optimize-memory-ddeea32`(local only,未 push);rebase **零冲突**
- 正确性:对齐输出 **bit-exact**(1,449,232 行全同)
- **Opt3**:GIX ktab **−7.69%**(840,595,878 → 775,934,664 B)——注:之前 stale binary 测不出
- **Opt1**:sort+align 长尾 **−93%**(2003 → 146 MB);峰值不变(峰在 seed merge)
- 图+文档 `docs/benchmark_optmem_vs_upstream.md`

### 3. 清理 worktree ✅
- 移除 FASTGA-old-upstream / FASTGA-main-upstream
- 保留 FASTGA-optmem-ddeea32(rebase 成果,可继续工作)
- optimize-memory 原分支 + 未提交改动全程未动
