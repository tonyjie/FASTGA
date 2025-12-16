# GRCh38 vs CHM13 FastGA 运行时间分析

## 输入文件
- **GRCh38**: `/scratch/jl4257/seq_align/fastga_datasets/GRCh38/GCF_000001405.40_GRCh38.p14_genomic.fna`
- **CHM13**: `/scratch/jl4257/seq_align/fastga_datasets/CHM13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna`

## 运行配置
- **线程数**: 1 thread (`-T1`), 8 threads (`-T8`), 和 32 threads (`-T32`)
- **输出目录**: `/scratch/jl4257/seq_align/fastga_datasets/GRCh38_vs_CHM13`
- **注意**: 1线程和32线程运行时，GDB和GIX文件已缓存，只运行了FastGA主程序

---

## 阶段时间分解

### Phase 1: GDB Creation (GRCh38)
**命令**: `FAtoGDB` for GRCh38

| 指标 | 时间 |
|------|------|
| User time | 6.606s |
| System time | 2.187s |
| **Wall time** | **8.794s** |
| CPU utilization | 100.0% |
| Peak memory | 0MB |

**说明**: 将GRCh38 FASTA文件转换为Genome Database格式

---

### Phase 2: GIX Creation (GRCh38)
**命令**: `GIXmake` for GRCh38

#### Sub-phase 2a: K-mer Partitioning
| 指标 | 时间 |
|------|------|
| User time | 49.565s |
| System time | 1.539s |
| **Wall time** | **7.880s** |
| CPU utilization | 648.5% (使用~6.5核心) |

**说明**: 将k-mers分区为16个部分进行并行处理

#### Sub-phase 2b: Sorting & Index Output
| 指标 | 时间 |
|------|------|
| User time | 3:45.285s (225.285s) |
| System time | 16.744s |
| **Wall time** | **47.762s** |
| CPU utilization | 506.7% (使用~5.1核心) |

**统计**:
- **Sampled k-mers**: 2,316,875,829 (73.9% of total)
- **处理了16个部分** (Sorting, Compressing, Outputting)

#### Phase 2 Total (GRCh38)
| 指标 | 时间 |
|------|------|
| User time | 4:34.850s (274.850s) |
| System time | 18.284s |
| **Wall time** | **55.643s** |
| CPU utilization | 526.8% |
| Peak memory | 2MB |

---

### Phase 3: GDB Creation (CHM13)
**命令**: `FAtoGDB` for CHM13

| 指标 | 时间 |
|------|------|
| User time | 6.498s |
| System time | 2.203s |
| **Wall time** | **8.703s** |
| CPU utilization | 100.0% |
| Peak memory | 0MB |

**说明**: 将CHM13 FASTA文件转换为Genome Database格式

---

### Phase 4: GIX Creation (CHM13)
**命令**: `GIXmake` for CHM13

#### Sub-phase 4a: K-mer Partitioning
| 指标 | 时间 |
|------|------|
| User time | 48.044s |
| System time | 1.277s |
| **Wall time** | **7.982s** |
| CPU utilization | 617.9% (使用~6.2核心) |

#### Sub-phase 4b: Sorting & Index Output
| 指标 | 时间 |
|------|------|
| User time | 3:46.288s (226.288s) |
| System time | 15.467s |
| **Wall time** | **47.579s** |
| CPU utilization | 508.1% (使用~5.1核心) |

**统计**:
- **Sampled k-mers**: 2,305,788,869 (74.0% of total)
- **处理了16个部分** (Sorting, Compressing, Outputting)

#### Phase 4 Total (CHM13)
| 指标 | 时间 |
|------|------|
| User time | 4:34.333s (274.333s) |
| System time | 16.745s |
| **Wall time** | **55.562s** |
| CPU utilization | 523.9% |
| Peak memory | 2MB |

---

### Phase 5: Adaptive Seed Merge
**命令**: `FastGA` - Adaptive seed merge for G1

| 指标 | 时间 |
|------|------|
| User time | 2:01.381s (121.381s) |
| System time | 13.447s |
| **Wall time** | **18.947s** |
| CPU utilization | 711.6% (使用~7.1核心) |

**统计**:
- **Total seeds**: 1,298,273,853
- **Average seed length**: 39.0 bp
- **Seeds per genome position**: 0.4

**说明**: 合并两个基因组的k-mer表，找到匹配的种子

---

### Phase 6: Seed Sort and Alignment Search
**命令**: `FastGA` - Seed sort and alignment search (18 parts)

| 指标 | 时间 |
|------|------|
| User time | 42:26.326s (2546.326s) |
| System time | 2:10.026s (130.026s) |
| **Wall time** | **14:30.410s (870.410s)** |
| CPU utilization | 307.5% (使用~3.1核心) |

**统计**:
- **Total hits over 85bp**: 843,675
- **Total alignments**: 791,660
- **Non-redundant alignments**: 518,037
- **Average alignment length**: 28,415 bp

**说明**: 
- 这是**最耗时的阶段**，占总时间的约85.5%
- 处理了18个部分（Loading, Sorting, Searching seeds）
- 执行实际的序列对齐操作

---

### Phase 7: Sorting and Merging Alignments
**命令**: `FastGA` - Final sorting and merging

**说明**: 此阶段的时间包含在Phase 6中，没有单独测量

---

## 总时间汇总

### FastGA主程序总时间
| 指标 | 时间 |
|------|------|
| User time | 44:27.707s (2667.707s) |
| System time | 2:23.486s (143.486s) |
| **Wall time** | **14:49.370s (889.370s)** |
| CPU utilization | 316.1% |
| Peak memory | 6MB |

### 完整流程总时间
| 指标 | 时间 |
|------|------|
| FAtoGDB (GRCh38) | 8.794s |
| GIXmake (GRCh38) | 55.643s |
| FAtoGDB (CHM13) | 8.703s |
| GIXmake (CHM13) | 55.562s |
| FastGA | 889.370s |
| **Total Wall Time** | **1018.072s (16:58)** |

---

## 完整时间线

### 8线程运行时间线

| Phase | 程序 | 描述 | Wall Time | % of Total | CPU Time | CPU % |
|-------|------|------|-----------|------------|----------|-------|
| 1 | FAtoGDB | GDB Creation (GRCh38) | 8.79s | 0.9% | 8.79s | 100% |
| 2a | GIXmake | GIX Partitioning (GRCh38) | 7.88s | 0.8% | 51.10s | 648% |
| 2b | GIXmake | GIX Sort & Index (GRCh38) | 47.76s | 4.7% | 242.03s | 507% |
| **2 Total** | **GIXmake** | **GIX Creation (GRCh38)** | **55.64s** | **5.5%** | **293.13s** | **527%** |
| 3 | FAtoGDB | GDB Creation (CHM13) | 8.70s | 0.9% | 8.70s | 100% |
| 4a | GIXmake | GIX Partitioning (CHM13) | 7.98s | 0.8% | 49.32s | 618% |
| 4b | GIXmake | GIX Sort & Index (CHM13) | 47.58s | 4.7% | 241.76s | 508% |
| **4 Total** | **GIXmake** | **GIX Creation (CHM13)** | **55.56s** | **5.5%** | **291.08s** | **524%** |
| 5 | FastGA | Adaptive Seed Merge | 18.95s | 1.9% | 134.83s | 712% |
| 6 | FastGA | Seed Sort & Alignment | 870.41s | 85.5% | 2676.35s | 308% |
| 7 | FastGA | Sort & Merge | (included) | - | - | - |
| **FastGA Total** | **FastGA** | **FastGA Complete** | **889.37s** | **87.4%** | **2811.18s** | **316%** |
| **TOTAL** | **All Programs** | **Complete Run** | **1018.07s** | **100%** | **3421.70s** | **336%** |

**8线程总运行时间**: 16分58秒 (1018.07秒)

### 1线程运行时间线 (FastGA only, GDB/GIX cached)

| Phase | 程序 | 描述 | Wall Time | % of Total | CPU Time | CPU % |
|-------|------|------|-----------|------------|----------|-------|
| 5 | FastGA | Adaptive Seed Merge | 136.37s | 4.9% | 136.37s | 100% |
| 6 | FastGA | Seed Sort & Alignment | 2645.03s | 95.1% | 2645.03s | 100% |
| 7 | FastGA | Sort & Merge | (included) | - | - | - |
| **FastGA Total** | **FastGA** | **FastGA Complete** | **2781.40s** | **100%** | **2781.40s** | **100%** |

**1线程FastGA运行时间**: 46分21秒 (2781.40秒)

### 32线程运行时间线 (FastGA only, GDB/GIX cached)

| Phase | 程序 | 描述 | Wall Time | % of Total | CPU Time | CPU % |
|-------|------|------|-----------|------------|----------|-------|
| 5 | FastGA | Adaptive Seed Merge | 7.68s | 1.3% | 177.92s | 2315% |
| 6 | FastGA | Seed Sort & Alignment | 584.18s | 98.7% | 3232.58s | 553% |
| 7 | FastGA | Sort & Merge | (included) | - | - | - |
| **FastGA Total** | **FastGA** | **FastGA Complete** | **592.05s** | **100%** | **3410.68s** | **576%** |

**32线程FastGA运行时间**: 9分52秒 (592.05秒)

### 性能提升总结

| 配置 | FastGA Wall Time | 相比1线程 | 相比8线程 | 1→N加速比 | 8→N加速比 |
|------|-----------------|-----------|-----------|-----------|-----------|
| 1 thread | 2781.40s (46:21) | - | - | 1.00x | - |
| 8 threads | 889.37s (14:49) | -1892.03s | - | **3.13x** | 1.00x |
| 32 threads | 592.05s (9:52) | -2189.35s | -297.32s | **4.70x** | **1.50x** |

**关键观察**:
- **1线程 → 8线程**: 加速3.13倍，效率很高
- **8线程 → 32线程**: 加速1.50倍，效率下降但仍有效
- **1线程 → 32线程**: 加速4.70倍，总体效果显著

---

## 关键观察

### 1. 时间分布
- **FAtoGDB** (Phases 1, 3): ~17秒 (1.7%)
  - GRCh38: 8.79秒
  - CHM13: 8.70秒
- **GIXmake** (Phases 2, 4): ~111秒 (10.9%)
  - GRCh38: 55.64秒
  - CHM13: 55.56秒
- **FastGA** (Phases 5-7): ~889秒 (87.4%) ⭐ **最耗时**
  - Adaptive Seed Merge: 18.95秒 (1.9%)
  - Seed Sort & Alignment: 870.41秒 (85.5%)

### 2. CPU利用率
- **最高**: Adaptive Seed Merge (711.6%) - 使用~7.1核心
- **最低**: Seed Sort & Alignment (307.5%) - 使用~3.1核心
- **平均**: ~384% (使用~3.8核心)

### 3. 基因组规模
- **GRCh38 k-mers**: 2,316,875,829 (73.9% sampled)
- **CHM13 k-mers**: 2,305,788,869 (74.0% sampled)
- **总种子数**: 1,298,273,853
- **最终对齐数**: 518,037 non-redundant alignments

### 4. 对齐质量
- **平均对齐长度**: 28,415 bp (非常长！)
- **总对齐数**: 791,660 (包含冗余)
- **非冗余对齐数**: 518,037
- **最小长度阈值**: 85bp

### 5. 性能瓶颈
- **FastGA程序**是最主要的瓶颈，占用约87.4%的总时间
  - 其中Phase 6 (对齐搜索) 占用约85.5%的总时间
- CPU利用率相对较低 (307.5%)，可能受I/O限制
- 可以考虑增加线程数或优化I/O来提升性能

---

## 与示例数据对比

| 指标 | GRCh38 vs CHM13 | HAP1 vs HAP2 (示例) | 比例 |
|------|-----------------|---------------------|------|
| Wall time | 889s (14:49) | 33.65s | **26.4x** |
| CPU time | 3414s | 146s | **23.4x** |
| K-mers (GRCh38) | 2.32B | 64.7M | **35.8x** |
| K-mers (CHM13) | 2.31B | 65.0M | **35.5x** |
| Seeds | 1.30B | 51.1M | **25.4x** |
| Alignments | 518K | 324K | **1.6x** |
| Avg alignment len | 28,415bp | 1,953bp | **14.6x** |

**说明**: GRCh38和CHM13是完整的人类基因组参考序列，比示例数据大得多，因此运行时间显著增加是预期的。

---

---

## 1线程运行结果对比

### 1线程配置运行结果

**运行配置**: 1 thread (`-T1`), GDB和GIX文件已缓存

#### Phase 5: Adaptive Seed Merge (1 thread)
| 指标 | 时间 |
|------|------|
| User time | 1:56.083s (116.083s) |
| System time | 20.255s |
| **Wall time** | **2:16.369s (136.369s)** |
| CPU utilization | 100.0% (单线程) |

**统计**:
- **Total seeds**: 1,298,273,853 (与其他配置相同)
- **Average seed length**: 39.0 bp
- **Seeds per genome position**: 0.4

#### Phase 6: Seed Sort and Alignment Search (1 thread)
| 指标 | 时间 |
|------|------|
| User time | 41:46.108s (2506.108s) |
| System time | 2:18.224s (138.224s) |
| **Wall time** | **44:05.032s (2645.032s)** |
| CPU utilization | 100.0% (单线程) |

**统计**:
- **Total hits over 85bp**: 843,675 (与其他配置相同)
- **Total alignments**: 791,660
- **Non-redundant alignments**: 518,037
- **Average alignment length**: 28,415 bp
- **处理了4个部分** (vs 8线程的18个部分，32线程的66个部分)

#### FastGA Total (1 thread)
| 指标 | 时间 |
|------|------|
| User time | 43:42.191s (2622.191s) |
| System time | 2:38.479s (158.479s) |
| **Wall time** | **46:21.402s (2781.402s)** |
| CPU utilization | 100.0% |
| Peak memory | 13MB |

---

## 32线程运行结果对比

### 32线程配置运行结果

**运行配置**: 32 threads (`-T32`), GDB和GIX文件已缓存

#### Phase 5: Adaptive Seed Merge (32 threads)
| 指标 | 时间 |
|------|------|
| User time | 2:30.015s (150.015s) |
| System time | 27.909s |
| **Wall time** | **7.684s** |
| CPU utilization | 2315.3% (使用~23.2核心) |

**统计**:
- **Total seeds**: 1,298,273,853 (与8线程相同)
- **Average seed length**: 39.0 bp
- **Seeds per genome position**: 0.4

#### Phase 6: Seed Sort and Alignment Search (32 threads)
| 指标 | 时间 |
|------|------|
| User time | 46:06.869s (2766.869s) |
| System time | 7:45.713s (465.713s) |
| **Wall time** | **9:44.182s (584.182s)** |
| CPU utilization | 553.4% (使用~5.5核心) |

**统计**:
- **Total hits over 85bp**: 843,675 (与8线程相同)
- **Total alignments**: 791,660
- **Non-redundant alignments**: 518,037
- **Average alignment length**: 28,415 bp
- **处理了66个部分** (vs 8线程的18个部分)

#### FastGA Total (32 threads)
| 指标 | 时间 |
|------|------|
| User time | 48:36.884s (2916.884s) |
| System time | 8:13.798s (493.798s) |
| **Wall time** | **9:52.046s (592.046s)** |
| CPU utilization | 576.1% |
| Peak memory | 20MB |

---

## 1线程 vs 8线程 vs 32线程性能对比

### FastGA主程序时间对比

| 阶段 | T=1 | | | T=8 | | | T=32 | | |
|------|-----|-|-|-----|-|-|------|-|-|
| | Wall Time | Speedup | CPU % | Wall Time | Speedup | CPU % | Wall Time | Speedup | CPU % |
| **Adaptive Seed Merge** | 136.37s | 1.00x | 100.0% | 18.95s | **7.20x** | 711.6% | 7.68s | **17.75x** | 2315.3% |
| **Seed Sort & Alignment** | 2645.03s | 1.00x | 100.0% | 870.41s | **3.04x** | 307.5% | 584.18s | **4.53x** | 553.4% |
| | (4 parts) | | | (18 parts) | | | (66 parts) | | |
| **FastGA Total** | 2781.40s | 1.00x | 100.0% | 889.37s | **3.13x** | 316.1% | 592.05s | **4.70x** | 576.1% |
| **Peak Memory** | 13MB | - | - | 6MB | - | - | 20MB | - | - |

### 关键发现

1. **总体性能提升**: 
   - **1线程 → 8线程**: Wall time从2781秒降至889秒，**加速3.13倍**
   - **8线程 → 32线程**: Wall time从889秒降至592秒，**加速1.50倍**
   - **1线程 → 32线程**: Wall time从2781秒降至592秒，**加速4.70倍**
   - 从1线程到32线程，节省了约36.5分钟（2189秒）

2. **Adaptive Seed Merge阶段**:
   - **1线程 → 8线程**: 从136.4秒降至18.9秒，**加速7.20倍** ⭐ 最佳加速
   - **8线程 → 32线程**: 从18.9秒降至7.7秒，**加速2.47倍**
   - **1线程 → 32线程**: 从136.4秒降至7.7秒，**加速17.75倍** ⭐ 最大加速
   - CPU利用率：100% → 712% → 2315%（使用~23核心）
   - 这是**最受益于多线程的阶段**

3. **Seed Sort & Alignment阶段**:
   - **1线程 → 8线程**: 从2645秒降至870秒，**加速3.04倍**
   - **8线程 → 32线程**: 从870秒降至584秒，**加速1.49倍**
   - **1线程 → 32线程**: 从2645秒降至584秒，**加速4.53倍**
   - CPU利用率：100% → 308% → 553%（但仍只使用~5.5核心）
   - 处理部分数：4 → 18 → 66（更好的并行度）
   - 系统时间：138s → 130s → 466s（32线程时I/O成为瓶颈）

4. **资源使用**:
   - CPU利用率：100% → 316% → 576%
   - 内存使用：13MB → 6MB → 20MB（32线程时最高）
   - 系统时间：158s → 143s → 494s（32线程时显著增加）

5. **瓶颈分析**:
   - **1线程 → 8线程**: 线性加速效果显著，CPU利用率提升至316%
   - **8线程 → 32线程**: 
     - Adaptive Seed Merge阶段充分利用了32线程（2315% CPU）
     - Seed Sort & Alignment阶段CPU利用率仍较低（553%），可能受：
       - I/O瓶颈（系统时间从130s增至466s）
       - 内存带宽限制
       - 算法本身的串行部分

6. **并行效率**:
   - **Adaptive Seed Merge**: 从1到8线程效率很高（7.20x加速），8到32线程仍有提升（2.47x）
   - **Seed Sort & Alignment**: 从1到8线程效率中等（3.04x加速），8到32线程提升有限（1.49x）
   - 总体：从1到8线程效率高（3.13x），8到32线程效率下降（1.50x）

---

## 建议优化

### 已验证的优化
1. ✅ **增加线程数**: 
   - **1线程 → 8线程**: 总体加速3.13倍，效率很高
   - **8线程 → 32线程**: 总体加速1.50倍，效率下降但仍有效
   - **1线程 → 32线程**: 总体加速4.70倍，效果显著
   - Adaptive Seed Merge阶段受益最大（1→32: 17.75x加速）
   - Seed Sort & Alignment阶段适度加速（1→32: 4.53x加速）

### 进一步优化建议
1. **I/O优化**: 
   - 对齐阶段系统时间显著增加（130s → 466s），I/O成为瓶颈
   - 考虑使用SSD或增加内存缓存
   - 使用更快的存储系统

2. **内存优化**:
   - 内存使用从6MB增至20MB，仍在合理范围
   - 可以考虑预加载数据到内存

3. **并行度调优**:
   - 32线程时处理66个部分（vs 8线程的18个部分）
   - 可以监控各部分负载均衡，进一步优化

4. **混合策略**:
   - Adaptive Seed Merge阶段：使用更多线程（32+）
   - Seed Sort & Alignment阶段：可能需要平衡线程数和I/O性能

