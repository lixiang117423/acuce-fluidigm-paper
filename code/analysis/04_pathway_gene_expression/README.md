# 04_pathway_gene_expression — 运行record

pathway水平的geneexpression重analysis：defence_pathways（Figure 3A/3B、4A）、69 分蘖期 BR+AUX（3C）、87 灌浆期 BR+SA+Drought（3D）、
独脚金内酯（4B）、高糖pathway（4C），以及pathway_enrichment_test（Table S3/S4）and高糖–防御相关性（Table S5）。

## 运行环境

- R：`/opt/homebrew/bin/R`（`R version 4.6.1 (2026-06-24)`），`/opt/homebrew/bin/Rscript "<脚本绝对路径>"`
- 关键包version：tidyverse 2.0.0、ggplot2 4.0.3、bioRtools 1.95.1、readxl 1.5.0.1、writexl 2.0.1
- 固定随机种子 `set.seed(12345)`

## 脚本and运行命令（2026-09-22，须按序运行）

```bash
cd "<repository_root>"
/opt/homebrew/bin/Rscript "code/04_pathway_gene_expression/01_pathway_gene_statistics.R"
/opt/homebrew/bin/Rscript "code/04_pathway_gene_expression/02_pathway_enrichment_test_TableS3S4.R"
/opt/homebrew/bin/Rscript "code/04_pathway_gene_expression/03_Figure3A_3B_defence_pathways.R"
/opt/homebrew/bin/Rscript "code/04_pathway_gene_expression/04_Figure3C_3D_focal_pathways.R"
/opt/homebrew/bin/Rscript "code/04_pathway_gene_expression/05_Figure4_cross_stage_pathways.R"
/opt/homebrew/bin/Rscript "code/04_pathway_gene_expression/06_TableS5_high_sugar_defence_correlation.R"
```

6 个脚本退出码均为 0；日志见 `result/04_pathway_gene_expression/00_运行日志_*.log` and `00_运行错误_*.log`。
（02 依赖 01 的 `stage_focal_level`，03–06 依赖 01 的 `clean_long` / `stage_level`。）

## 口径and实现notes

- **gene/pathway**：92 个gene / 20 条pathway（引物面板 96 个位点 = 92 gene + 4 内参；见主题 03 README）。
- **离群值**：按 stage × focal × neighbor × gene 用 `bioRtools::find_outliers()` 剔除，and主题 03 一致。
- **Fold change**：`log2FC = log2(mean_Mix / mean_Pure)`（线性尺度均值比）。
- **两种粒度**：① stage 层面（4 个 focal 合并，and旧版 Figure 4A/4C 口径一致）；② stage × focal 层面（用于 3C/3D）。
- **testandcorrection**：Wilcoxon 秩和为主、Welch t 为对照；BH-FDR 按test族分别correction（stage 层面按时期、stage × focal 层面按 时期 × focal）。
- **缩放**：图内expression为gene内 0–1 缩放（and旧版 `bioRtools::scale01()` 一致），使不同gene可同图展示。
- **热图实现**：and主题 03 相同，用 `ggplot2::geom_tile` 绘制（`pheatmapColor` 未安装、`ggplotify` 转换报 gpar 错误）。

### OsXTR1 and旧版剔除gene的处理

- **OsXTR1**：Mix 条件下expression近似 0（多数sample为 0），`log2FC` 为 −6.46（69 分蘖期）and −Inf（87 灌浆期）。
  旧版用 `filter(gene != "OsXTR1")` 直接剔除、未notes原因；本次**同样不纳入图内**，但在statistics表中保留该gene行，
  并在 README 写明原因（比值为 0 时 log2FC and可视化均无意义）。
- **OsGH9B9 / OsPR1b / OsPR5 / OsSPL7**：旧版在 Figure 4A 中无理由剔除；本次**全部纳入 Figure 4A**（12 个防御gene全部展示）。
  其中 OsGH9B9、OsPR1b 在 Mix 条件下大量为 0（log2FC 分别低至 −8.8 and −12.9），正文引用时须注意。

## 输出

| file | 内容 |
|---|---|
| `01_表_pathwaygene差异.xlsx` | 5 sheet：stage_level、stage_focal_level、combo_means、clean_long（清洗后长表）、gene_annotation |
| `02_表_pathway富集_TableS3S4.xlsx` | **Table S3/S4**：Fisher 富集（判据 A 上调 / 判据 B 显著）× 全部 focal × 时期 |
| `03_Figure3A_防御gene热图.pdf` / `.png` | Defense+JA+SA 17 gene × 32 组合（10 × 6 in） |
| `04_Figure3B_防御gene整体差异.pdf` / `.png` | 防御gene整体comparison（8 × 6 in） |
| `05_表_防御整体test.xlsx` | Figure 3B 的statistics表 |
| `06_Figure3C_69分蘖期_BR_AUX.pdf` / `.png`、`07_Figure3D_87灌浆期_BR_SA_Drought.pdf` / `.png` | 焦点纯系pathwaygene（8 × 8 in） |
| `08_Figure4A_防御genecross_stage.pdf` / `.png` | 12 个防御gene × 2 时期（15 × 6 in） |
| `09_Figure4B_独脚金内酯cross_stage.pdf` / `.png` | 3 个独脚金内酯gene × 2 时期（6 × 6 in，**旧版缺失**） |
| `10_Figure4C_高糖pathwaycross_stage.pdf` / `.png` | 5 个高糖gene × 2 时期（13.3 × 5 in） |
| `11_表_cross_stagepathway_gene_statistics.xlsx` | 3 sheet：defense、strigolactone、high_sugar |
| `12_表_TableS5_high_sugar_defence_correlation.xlsx` | **Table S5**：高糖 × 防御gene相关（170 个组合） |
| `13_FigureS5_高糖and防御相关性热图.pdf` / `.png` | 相关性热图（8 × 8 in） |

## results要点（2026-09-22）

### defence_pathways（3A/3B、4A）

- **整体**（gene内 0–1 缩放后comparison）：tillering 混种显著高于纯种（Wilcoxon P = 0.0042，n = 908 vs 351）；
  filling 无差异（P = 0.804）。and"分蘖期激活、灌浆期减弱"的叙事一致。
- **分蘖期gene水平**（stage 层面）：OsGNS1（log2FC = 1.66，P = 4.0e-4，**FDR = 0.018**）、OsPAL（1.34，P = 0.0099，FDR = 0.066）、
  OsNPR1（0.68，P = 0.034）、OsPOX223（0.93，P = 0.040）、OsSPL7（0.80，P = 0.048）方向为上调；
  OsPBZ1（P = 0.125）、OsRBBI2（P = 0.310）、OsPR5（P = 0.507）**未达显著**。
- **灌浆期gene水平**：12 个gene全部不显著（最小 P = 0.085，OsPAL，方向为下调）。

### focal_pathways（3C/3D）

- **69 分蘖期 BR + AUX**：OsBLE1 显著上调（log2FC = 1.10 即 ≈2.1 倍，P = 0.034）；OsIAA4 上调（log2FC = 3.20，P = 0.010）；
  OsGH3.2（log2FC = 0.99，P = 0.095）、OsIAA1（0.99，P = 0.078）为趋势上调；OsARF1、OsBLE2 不显著。
- **87 灌浆期 BR + SA + Drought**：OsWRKY7（log2FC = 1.26，P = 0.0077）、OsWRKY45（0.88，P = 0.049）显著上调；
  OsDHN1（1.75）、OsWRKY82（2.50）、OsWRKY114（2.16）方向上调但不显著；
  **OsBLE1 在 87 灌浆期为 +0.38（≈1.3 倍），P = 0.178，未达显著**（旧框架称 P = 0.05）。

### 独脚金内酯（4B，旧版缺失）

- 3 个gene在两个时期均不显著：tillering 的 OsD14-1 log2FC = +0.12、P = 0.077（FDR = 0.22）、OsD53 +0.02、OsD53like-1 −0.20；
  filling 各gene P = 0.22–0.55。**"分蘖期特异激活"的说法在本次重analysis中不成立**，只能表述为方向性趋势。

### 高糖pathway（4C）

- **分蘖期**：OsGBSSII 显著上调（log2FC = 1.22，P = 0.0013，**FDR = 0.039**）；
  OsHXK5（+0.47，P = 0.011）、OsHXK6（+0.40，P = 0.012）为趋势（FDR = 0.066）；OsCRCT、OsHXK2 不显著。
- **灌浆期**：5 个gene全部不显著（P = 0.46–0.99），但除 OsHXK5（−0.04）外方向仍为上调。and"高糖信号cross_stage持续、灌浆期减弱"的表述基本一致。

### pathway富集（Table S3/S4）—— 重要负面results

- 用 Fisher 精确test（因每条pathway仅 1–12 个gene，卡方不适用）在两个判据下（上调gene富集 / 显著gene富集）
  对 **全部 focal × 时期** 组合test，**没有任何pathway达到 FDR < 0.05**；
  69 分蘖期最低的原始 P 为 Auxin（4/4 gene上调，P = 0.134，FDR = 1.00），BR 为 P = 1.00；
  87 灌浆期同样无显著pathway。
- 结论：**paper框架里"BR P = 0.041、AUX P = 0.005、87 灌浆期 BR P = 0.023 / drought P = 0.010 / SA P = 0.022"无法复现**
  （旧代码里也没有任何富集test），正文不宜再使用这些数字。建议改为报告gene水平的 log2FC and FDR，
  并用"方向一致的gene比例"描述pathway趋势。

### 高糖–防御相关性（Table S5）

- 170 个组合中 169 个可计算（1 个因expression无变异为 NA）；**47 个组合 FDR < 0.05**（未correction P < 0.05 的 64 个），
  且集中在**灌浆期**。最强相关：OsHXK2–OsWRKY82（r = 0.79，FDR = 7.7e-15）、OsHXK5–OsWRKY82（r = 0.55，FDR = 1.2e-05）、
  OsGBSSII–OsWRKY82（r = 0.52，FDR = 6.4e-05）；OsHXK2/OsHXK5 and OsPBZ1、OsPR5、OsPOX223 亦显著正相关。
- 这为"高糖信号and防御维持同步"提供了sample级证据，但属**相关**results，不能直接推因果。

## and旧版 / paper框架的交叉核对（须作者确认）

| 项目 | 旧版或框架 | 本次重analysis | 结论 |
|---|---|---|---|
| 防御gene整体 | 分蘖期上调、灌浆期减弱 | tillering P = 0.0042；filling P = 0.804 | **一致** |
| 分蘖期显著防御gene | OsGNS1、OsPBZ1、OsPOX223、OsPR5、OsRBBI2 | OsGNS1（FDR 0.018）、OsPAL、OsNPR1、OsPOX223、OsSPL7；OsPBZ1/OsPR5/OsRBBI2 不显著 | **部分不一致**，正文gene名单需按新results核对 |
| OsBLE1（69 分蘖期） | FC = 2.22、P = 0.001 | log2FC = 1.10（FC ≈ 2.14）、P = 0.034 | **FC 一致、P 更保守** |
| OsBLE1（87 灌浆期） | 显著上调 P = 0.05（全文最强分子证据） | log2FC = +0.38、P = 0.178，**不显著** | **不一致，影响核心结论**：跨品系、cross_stage一致上调的表述需改写 |
| OsWRKY45（87 灌浆期） | FC = 1.86、P = 0.01 | log2FC = 0.88（FC ≈ 1.84）、P = 0.049 | FC 一致、P 更保守（仍显著） |
| 独脚金内酯（分蘖期特异激活） | OsD14-1 显著上调 | OsD14-1 P = 0.077（FDR 0.22），全部不显著 | **不一致**，只能表述为趋势 |
| 高糖pathway | 分蘖期 3 个gene显著、灌浆期持续上调 | 分蘖期 OsGBSSII FDR 0.039（另 2 个趋势）；灌浆期全部不显著但方向多上调 | **基本一致** |
| pathway富集卡方 | BR/AUX/SA/drought 若干显著 | **旧代码无此analysis**；Fisher test下无任何pathway FDR < 0.05 | **旧数字不可复现，须删除或重写** |
| 剔除的gene | OsXTR1（3C/3D）、OsGH9B9/OsPR1b/OsPR5/OsSPL7（4A）无理由剔除 | OsXTR1 仍不画（Mix 近似 0，已写明原因）；4A 全部gene均展示 | 图内容调整，已notes理由 |

## 备注

- 所有 P 值在图内统一用定点写法（`P < 0.001` / `P = 0.012`），避免 R and Excel 的科学计数法被误读为显著性字母。
- `01_表_pathwaygene差异.xlsx` 的 `clean_long` 为清洗后长表（供plotting复用），`combo_means` 为 32 个组合的均值表。

---

## 追加：统一口径为 limma + camera（2026-09-22，脚本 07–08）

```bash
/opt/homebrew/bin/Rscript "code/04_pathway_gene_expression/07_limma_unified_statistics.R"
/opt/homebrew/bin/Rscript "code/04_pathway_gene_expression/08_pathway_gene_plotting_limma_version.R"
```

### 为什么换口径

• 逐gene Wilcoxon（脚本 01）：每个gene独立估方差 + 92 个gene做 BH correction → 分蘖期仅 8 个gene FDR < 0.05、灌浆期 0 个。

• 逐pathway Fisher（脚本 02）：每条pathway只有 1–12 个gene（多数 2–4 个），2×2 表几乎无检出力 → 全部阴性。

• limma（`~ 0 + group + focal`，focal 作固定区组 + 方差收缩，**并显式构造 Mix − Pure 对比矩阵**）：
在本data上**两个时期均无gene达到 FDR < 0.05**（分蘖期 14 个gene P < 0.05，最强 OsGNS1 logFC = 1.15、P = 0.0024；
灌浆期仅 1 个 P < 0.05）。结论：limma and Wilcoxon 的口径**并不冲突**，两者都只能给出名义显著信号。

> **⚠️ 2026-09-22 更正（重要）**：本 README and CHANGELOG 早期version曾写"limma 下分蘖期 53/92、灌浆期 59/92 个gene
> FDR < 0.05、OsBLE1 FDR ≈ 1e-11"，**该results是错误编码所致**——当时用 `topTable(fit, coef = "groupMix")` 取系数，
> 而在 `~ 0 + group`（无截距 + 单因子）的写法下**系数是组均值本身**（Mix 组在参考 focal 的对数均值），
> 并非 Mix − Pure 的差异，实际test的是"该gene是否expression"这一平凡假设。改正为显式对比矩阵后，
> 正确results是：**无gene通过 FDR < 0.05**（详见下表）。已重跑脚本 07–09 and主题 06 脚本 06，并同步更正全部record。

• camera（limma 体系内的**竞争性gene集test**，会correctiongene间相关）：用于替代 Fisher 做pathwaytest，适合 2–15 个gene的小pathway。

### resultsand结论变化

| 项目 | Wilcoxon（脚本 01） | limma（脚本 07，**已更正**） | 对正文的影响 |
|---|---|---|---|
| 单gene显著性（FDR < 0.05） | 分蘖期 8 个、灌浆期 0 个 | **分蘖期 0 个、灌浆期 0 个** | 两套方法一致：**没有gene能通过多重correction**，单gene结论只能写"名义趋势（P < 0.05，未correction）" |
| 单gene名义显著性（P < 0.05，未correction） | —— | 分蘖期 14 个、灌浆期 1 个 | 可用于描述性表述，须注明未correction |
| `OsBLE1`（69 分蘖期，focal 内） | P = 0.034 | logFC = 0.72，P = 0.024、FDR = 0.33 | 两法一致：**名义显著、correction后不显著**；不能称"显著上调" |
| `OsBLE1`（87 灌浆期，focal 内） | P = 0.178 | logFC = 0.10，P ≈ 0.49 | **不显著**：manuscript"跨品系cross_stage一致上调"的表述**不能成立**，只能写方向一致的趋势 |
| `OsD14-1`（69 分蘖期，focal 内） | P = 0.077 | logFC = −0.28，P = 0.35 | 不显著；"分蘖期独脚金内酯上调"**不能写** |
| 防御gene整体 | P = 0.0042（分蘖期"显著上调"） | **P = 0.172**（分蘖期）、0.560（灌浆期） | **必须改口**：原 P 值把 1259 个"gene × sample"当独立观测（gene层面的伪重复），每sample一个分数才是正确口径 → defence_pathways整体不能称显著，只能写一致的上升趋势 |
| pathway水平（Table S3/S4） | Fisher 全部阴性 | camera **同样无 FDR < 0.05**（分蘖期最小 Oxidative stress P = 0.069；灌浆期最小 Drought P = 0.051） | 两套方法一致：pathway只能写方向andeffect量，不能写"显著富集" |
| 关键gene方向（已更正） | —— | 69 分蘖期（focal 内 limma）：OsGBSSII logFC = 1.20（P = 0.014）、OsGNS1 1.65（P = 0.028）、OsBLE1 0.72（P = 0.024）、OsGH3.2 0.85（P = 0.078）、OsIAA1 0.81（P = 0.136）、OsWRKY45 0.48（P = 0.129）、OsD14-1 −0.28（P = 0.35）| 分蘖期以高糖、生长素、防御相关gene为主的名义上调趋势；**均未通过 FDR correction** |

### 输出（new_definition）

| file | 内容 |
|---|---|
| `14_表_limma_unified_statistics.xlsx` | 6 sheet：stage 层面、stage × focal 层面、camera pathwaytest、防御分数test、关键gene、pathway定义 |
| `15_Figure3B_防御gene整体分数.pdf` / `.png` | 防御pathway_scores（每sample一个值），8 × 6 in |
| `16_Figure3C_69分蘖期_limma_version.pdf` / `.png`、`17_Figure3D_87灌浆期_limma_version.pdf` / `.png` | 焦点纯系 logFC 点图（focal 内 limma），7 × 4.5 / 7 × 6 in |
| `18_Figure4A_防御gene_limma_version.pdf` / `.png` | 12 个防御gene × 2 时期（13 × 7 in） |
| `19_Figure4B_独脚金内酯_limma_version.pdf` / `.png` | 独脚金内酯 3 gene × 2 时期（8 × 4 in） |
| `20_Figure4C_高糖pathway_limma_version.pdf` / `.png` | 高糖 5 gene × 2 时期（9 × 5 in） |

### 两点必须写进方法/讨论的局限

• 本表**没有芯片/批次变量**，无法直接correction上机批次；4 个 focal 方向 100% 一致可降低（但不能排除）批次混杂的可能。

• limma 的 logFC 是 **log2(expression + 1) 尺度、经 focal correction**后的effect，and旧版"均值比 logFC"不是同一个量
（例如 OsGNS1 分蘖期：limma 1.15 vs 均值比 1.66），正文引用数值时必须注明尺度。

### 建议的finalversion（投稿用）

• statisticsand图：以 `14_表_limma_unified_statistics.xlsx` and `15–20` 图为准；脚本 01–02 的 Wilcoxon/Fisher version作为方法学对照保留（可放补充材料notes"两种口径结论方向一致、limma 检出力更高"）。

• 仍需作者确认的三条：① defence_pathways整体不再称"显著上调"；② pathway富集不写 P 值，只写方向；③ 全部gene水平数值标注 limma and均值比两种尺度。

---

## 追加：three_stages扩展（2026-10-18，脚本 10–11）

按作者决定把花期（flowering）接入主题 05。判定性analysis见主题 `08_three_stagesexpression`（结论：不影响主线，纳入）。

```bash
/opt/homebrew/bin/Rscript "code/04_pathway_gene_expression/10_limma_unified_statistics_three_stages.R"
/opt/homebrew/bin/Rscript "code/04_pathway_gene_expression/11_pathway_gene_plotting_limma_version_three_stages.R"
```

- 输入：`result/03_overall_gene_expression/11_表_expression矩阵_three_stages.xlsx`（sample级 254 × 92）；口径and脚本 07/08 完全一致，仅时期数 2 → 3；输出为**新编号file**，未覆盖两时期version。

### 10 limma 统一口径statistics（three_stages）

输出 `result/04_pathway_gene_expression/23_表_limma_unified_statistics_three_stages.xlsx`（6 sheet：gene_stage_level / gene_stage_focal_level / pathway_camera / defense_score / key_genes / pathway_definition）。

| 项目 | tillering（n = 80） | flowering（n = 86） | filling（n = 88） |
|---|---|---|---|
| stage 层面名义 P < 0.05 | **14** | 1 | 1 |
| stage 层面 FDR < 0.05 | 0 | 0 | 0 |
| 防御pathway_scores（Mix − Pure） | **+0.217**（P = 0.172） | +0.075（P = 0.640） | −0.079（P = 0.560） |
| camera 最小 FDR | 0.419 | 0.702 | 0.461 |

- **and两时期口径完全一致**（分蘖期 14 个名义、灌浆期 1 个、三期均无 FDR 显著、防御分数三点单调下降）。
- **stage × focal 层面的新细节（须注意）**：three_stages下 **flowering × focal 49 有 4 个gene通过该族内的 FDR**（OsNADP-isocitrate logFC = −0.202、P = 1.7e-4、FDR = 0.016；OsNIA1 −0.003、4.8e-4、0.022；OsGlut_dehydr −0.079、1.5e-3、0.047；OsHXK5 −0.095、2.0e-3、0.047），而两时期version在**任何 时期 × focal 组合里都没有 FDR < 0.05**。这批命中的**effect量极小**（其中一个 |logFC| = 0.003），属"statistics上可检出、生物学上可忽略"，正文若引用必须同时给 logFC andsample量，不能只写"显著"；它也**不改变** stage 层面的主结论（三期均无gene过 FDR）。

### 11 pathway图（three_stages，面板级）

输出 `result/04_pathway_gene_expression/24–28_*_three_stages.pdf / .png`（PDF 页面尺寸已逐张 `pdfinfo` 核对）：

| 编号 | 内容 | 尺寸 |
|---|---|---|
| `24_Figure3A_防御gene热图_three_stages` | 17 个防御相关gene（Defense 12 + Jasmonic acid 3 + Salicylic acid 2）× 48 组合，`facet_grid(process + stage ~ group)` | 11 × 11 in（792×792 pt） |
| `25_Figure3B_防御gene整体分数_three_stages` | 每sample一个防御分数、按三个时期分面；图注给出各期 n and `lm(score ~ group + focal)` 口径 | 11 × 4.5 in（792×324 pt） |
| `26_Figure4A_防御gene_limma_version_three_stages` | 12 个防御gene跨three_stages（stage 层面 limma） | 13 × 7 in（936×504 pt） |
| `27_Figure4B_独脚金内酯_limma_version_three_stages` | 3 个独脚金内酯gene跨three_stages（按作者决定该图进补充材料） | 9 × 4 in（648×288 pt） |
| `28_Figure4C_高糖pathway_limma_version_three_stages` | 5 个高糖gene跨three_stages | 10 × 5 in（720×360 pt） |

- **Figure 3C / 3D 未重做**：它们分别是"69 分蘖期"and"87 灌浆期"的焦点纯系图，and是否加入花期无关，继续沿用两时期版（`16` / `17` 号file）。
- 本步只产出**面板级**file；final组图按 2026-10-17 规则由 VectorCraft 拼版。

### 12 stage_specific_comparison（three_stages）

输出 `result/04_pathway_gene_expression/29_表_stage_specific_comparison_three_stages.xlsx` and `30_/30b_FigureS15_*_three_stages.pdf / png`（各 8×6 in）。

- 模型：`grp = stage × group`（6 个水平）+ `~ 0 + grp + focal`（focal 固定区组）；对比 = 三期各自的 (Mix − Pure) and三个两两交互（花期−分蘖、灌浆−分蘖、灌浆−花期）。
- **results**：在该decomposition模型下 **tillering 有 17 个gene FDR < 0.05、flowering 0 个、filling 0 个、交互项 0 个**；重叠结构只有 "Tillering only"（17 个）and "Neither"（75 个）。
- effect量：中位 |logFC| = 0.039（分蘖）/ 0.017（花期）/ 0.024（灌浆）；分蘖 vs 花期 Wilcoxon P = 4.3e-10、分蘖 vs 灌浆 P = 5.2e-06、花期 vs 灌浆 P = 1.5e-04；三期 logFC 之间呈弱负相关（r ≈ −0.26 ~ −0.30）。
- **须注意的口径差异（重要）**：本脚本（联合 6 组模型）给出"分蘖期 17 个 FDR 显著"，而脚本 10 的**逐期模型**（每期单独 `~ 0 + group + focal`）给出 **0 个**。两者差别来自design/方差收缩方式（联合模型把所有时期的残差合并、并按 6 组重新做经验贝叶斯收缩）。两时期版同一现象已存在（S15 = 8 个 vs 主文 0 个）。**正文必须写明用的是哪一个模型**，不能同时用两种口径的数字。

### 13 / 14 Table S5（three_stages）and一次口径correction

```bash
/opt/homebrew/bin/Rscript "code/04_pathway_gene_expression/13_TableS5_high_sugar_defence_correlation_three_stages.R"
/opt/homebrew/bin/Rscript "code/04_pathway_gene_expression/14_TableS5_high_sugar_defence_correlation_three_stages_robust_definition.R"
```

- **脚本 13（沿用旧口径：未过滤 Pearson）**：`31_表_*_three_stages.xlsx` + `32_FigureS5_*_three_stages` —— three_stages FDR < 0.05 组合数 = tillering 26 / flowering 43 / filling 11（各期 85 个组合）。
- **发现并correction一个假相关问题（重要）**：脚本 13 里出现 `OsCRCT–OsPR1b` 在花期 **r = 0.9999915** 的组合。查证：`OsPR1b` 在分蘖期 80 个sample中 **72 个为 0**（最大 0.88）、`OsGH9B9` 检出率仅 5%——这类"绝大多数为零"的gene做 Pearson，相关系数完全由一两个非零点决定。**旧版两时期 Table S5 使用同一口径，存在同类风险**。
- **脚本 14（robust_definition）**：先按 **检出率（expression > 0）≥ 50%（每期各自计算）** 过滤gene，再用 **Spearman 秩相关**，BH-FDR 按时期correction。输出 `33_表_*_robust_definition.xlsx` + `34_FigureS5_*_robust_definition`（14×5.5 in）。
  - 被过滤的gene：三个时期都剔除 **OsGH9B9 / OsPR1b / OsJAZ5**，分蘖期另剔除 **OsWRKY82**（检出率 0.35）。
  - correction后的results：每期 65–70 个组合、最大 |rho| = **0.63**（不再是 0.99）；FDR < 0.05 = **tillering 17 / flowering 7 / filling 12**；最强组合仍是 **OsHXK2–OsWRKY82（灌浆期 rho = 0.63、FDR = 3.5e-9）**，and两时期版报告的最强组合一致。
  - **建议**：投稿版 Table S5 / Figure S5 使用脚本 14 的robust_definition；脚本 13 的outputs保留作方法学对照（并在正文/补充notes里写清过滤and秩相关的理由）。
