# 03_overall_gene_expression — 运行record

Fluidigm expressiondata的整体格局analysis：expression_matrix_build、PERMANOVA（Table S2）、PLS-DA（Figure 2B–2E）、
整体expression热图（Figure 2A）、VIP geneand遗传距离的相关（Figure S3）。

## 运行环境

- R：`/opt/homebrew/bin/R`（`R version 4.6.1 (2026-06-24)`），`/opt/homebrew/bin/Rscript "<脚本绝对路径>"`
- 关键包version：tidyverse 2.0.0、ggplot2 4.0.3、ropls 1.44.0、vegan 2.7.6、bioRtools 1.95.1、
  ggrepel 0.9.8、readxl 1.5.0.1、writexl 2.0.1
- 固定随机种子 `set.seed(12345)`（含 PERMANOVA 的 999 次置换）

## 脚本and运行命令（2026-09-22，须按序运行）

```bash
cd "<repository_root>"
/opt/homebrew/bin/Rscript "code/03_overall_gene_expression/01_expression_matrix_build.R"
/opt/homebrew/bin/Rscript "code/03_overall_gene_expression/02_PERMANOVA_TableS2.R"
/opt/homebrew/bin/Rscript "code/03_overall_gene_expression/03_PLSDA.R"
/opt/homebrew/bin/Rscript "code/03_overall_gene_expression/04_VIP_and_genetic_distance_correlation.R"
/opt/homebrew/bin/Rscript "code/03_overall_gene_expression/05_gene_expression_plotting_Figure2.R"
```

5 个脚本退出码均为 0；日志见 `result/03_overall_gene_expression/00_运行日志_*.log` and `00_运行错误_*.log`。

## 输入anddata口径

- 输入：`data/processed_inputs/04_fluidigm_normalized_expression.xlsx`（15456 行 × 10 列）
- **Fluidigm 覆盖 4 个焦点纯系的全部 16 个 focal × neighbor 组合（含 4 个纯种），每个时期 16 个组合，合计 32 个**；
  sample级共 168 个sample（tillering 80、filling 88）。
- genetics 表：`data/02_raw_data/05_遗传距离/月亮谷纯系的遗传距离.xlsx`（取 A15/A49/A69/A87 的 4 × 4 子矩阵）。
- 离群值：`bioRtools::find_outliers()` 按 stage × focal × neighbor × gene 判定，共标记 996 个（6.44%）；
  组合均值矩阵用剔除后的data，sample级矩阵保留completedata（每个 (stage, sample, gene) 只有 1 个值，剔除后会出现 NA）。

### 两个必须记住的技术点

1. **矩阵 sheet 的读取方式**：`bioRtools::read_data()` 的列名修复（`.name_repair = "universal"`）会把gene名中的连字符
   改成点号（`OsD14-1` → `OsD14.1`、`OsCPS4-1` → `OsCPS4.1`、`OsNADP-isocitrate` → `OsNADP.isocitrate`），
   因此读取 `01_表_expression矩阵.xlsx` 的矩阵 sheet **必须**用 `readxl::read_excel(..., .name_repair = "minimal")`。
   同理，矩阵/得分矩阵的**行名必须显式设置为组合标签**，否则 `opls()`/`prcomp()` 的得分行名退化为 1..n，and元信息 join 后全为 NA。
2. **ropls 的 `orthoI` 参数**：`ropls::opls(x, y, predI = 3)`（默认 `orthoI = 0`）and旧版一致，可建出模型；
   若写 `orthoI = NA`（自动决定正交成分），ropls 会先做交叉验证，本data"第一预测成分即不显著"→ **直接返回空模型**
   （无 p1/p2 得分）。这本身notes Pure/Mix 在组合均值层面的分离较弱，写作时须谨慎表述。

## 方法

| 脚本 | 方法 |
|---|---|
| `01` | build（a）组合均值矩阵 32 × 92（已剔离群值，供 PLS-DA）and（b）sample级矩阵 168 × 92（completedata，供 PERMANOVA） |
| `02` | `vegan::adonis2`：`d ~ stage + group + focal + stage:group`，`by = "terms"` and `by = "margin"` 各一次，999 置换；距离用"gene Z-score 后的欧氏距离"and"log2(x+1) 的 Bray–Curtis"两种；另按时期分层建模 `~ group + focal` |
| `03` | `ropls::opls(predI = 3)`；输出 scores / loadings / VIP / 每成分 R2X and交叉验证statistics；主轴test用 Welch t and Wilcoxon（分时期comparison Pure vs Mix，以及时期之间comparison），并给 Cohen's d |
| `04` | VIP 前 10 gene在 12 个混种组合中的expression vs 遗传距离，Pearson + FDR（BH，按时期） |
| `05` | Figure 2A 用 `ggplot2::geom_tile` 直接绘制热图（原因见下）；2B–2E 用 `theme_prism()` + `write_data(units="in")` |

**热图实现的调整**：旧版用 `pheatmapColor::pheatmap_bio()`（本机未安装）；改用 `pheatmap` + `ggplotify::as.ggplot()`
会报 `'gpar' 成分 'fill' 的长度不能为零`（version兼容问题），故 2A 直接用 `geom_tile` 绘制：
行内 0–1 缩放、gene按pathway + 载荷排序、列按 时期 × 种植模式 分面，内容and旧图一致。

## 输出

| file | 内容 |
|---|---|
| `01_表_expression矩阵.xlsx` | 5 sheet：gene_annotation、sample_info、combo_matrix（32×92+元信息）、sample_matrix（168×92）、combo_long |
| `02_表_PERMANOVA_TableS2.xlsx` | **Table S2**：permanova（4 套data × terms/margin × 分层）、sample_info |
| `03_表_PLSDA.xlsx` | 6 sheet：scores、loadings、vip、vip_top10、variance_explained、axis_tests、cv_summary |
| `04_表_VIP_and_genetic_distance_correlation.xlsx` | 2 sheet：vip_gene_expression、correlation |
| `05_Figure2A_整体expression热图.pdf` / `.png` | 10 × 12 in |
| `06_Figure2B_PLSDA得分图.pdf` / `.png` | 7 × 5.5 in |
| `07_Figure2C_Dim1分布.pdf`、`08_Figure2D_Dim2分布.pdf` / `.png` | 7 × 6 in |
| `09_Figure2E_VIP载荷图.pdf` / `.png` | 9 × 6 in |
| `10_FigureS3_VIPand遗传距离.pdf` / `.png` | 14 × 6 in |

## results要点（2026-09-22）

### PERMANOVA（Table S2，`by = "terms"`，主模型 `~ stage + group + focal + stage:group`）

| data / 距离 | stage | group | focal | stage:group |
|---|---|---|---|---|
| 组合均值 · 欧氏（Z-score） | R² = 28.7%，P = 0.001 | R² = 4.25%，P = 0.062 | R² = 9.39%，P = 0.078 | R² = 2.79%，P = 0.225 |
| 组合均值 · Bray–Curtis | R² = 76.5%，P = 0.001 | R² = 1.10%，P = 0.221 | R² = 4.87%，P = 0.056 | R² = 0.83%，P = 0.250 |
| sample级 · 欧氏（Z-score） | R² = 10.2%，P = 0.001 | R² = 1.09%，P = 0.051 | R² = 2.37%，P = 0.034 | R² = 0.51%，P = 0.470 |
| sample级 · Bray–Curtis | R² = 52.0%，P = 0.001 | R² = 0.60%，P = 0.100 | R² = 2.84%，P = 0.001 | R² = 0.39%，P = 0.172 |

- **发育时期是最主要的变异来源**（4 套analysis中 stage 始终 P = 0.001）；种植模式（group）and品系（focal）effect小且多数处在边缘（P ≈ 0.05–0.10）。
- **分时期分层建模支持"分蘖期是主要响应窗口"**：组合均值 · Bray–Curtis 下，tillering 的 group effect R² = 18.4%、P = 0.025，
  而 filling 的 group effect R² = 3.5%、P = 0.628；sample级 · 欧氏距离下 group 的 P 在 tillering 也小于 filling。这一条andpaper叙事一致。

### PLS-DA（Figure 2B–2E）

- 模型：R²X(cum) = 0.538、R²Y(cum) = 0.638、Q²(cum) = 0.224；每成分 R2X = 19.3% / 23.6% / 10.9%。
- 得分主轴test（n = 4 Pure vs 12 Mix，注意不平衡）：
  - Dim 1：tillering vs filling Wilcoxon P = 8.7e-05（Welch P = 0.052）；
    时期内 Pure vs Mix 的 Wilcoxon P = 0.0022（tillering）、0.0011（filling），而 Welch t 分别为 0.178 and 0.204。
  - Dim 2：tillering vs filling P = 6.4e-07；filling 内 Pure vs Mix Welch P = 0.033、Wilcoxon P = 0.0011。
- VIP 前 10 gene：`OsCPS4-1`、`OsJAZ5`、`OsXTR1`、`OsGH9B9`、`OsKS4`、`OsPR1b`、`OsIPA1`、`OsWRKY45`、`OsAOS2`、`OsCitr_synth`。

### VIP geneand遗传距离（Figure S3）

- 12 个混种组合中，10 个 VIP gene × 2 时期共 20 个相关analysis**全部不显著**（最大 |r| = 0.33，最小 P = 0.29，FDR < 0.05 的组合数为 0），
  and表型侧的阴性results一致：互作的分子effect同样不能由遗传距离预测。

## and旧版 / paper框架的交叉核对（须作者确认）

| 项目 | 旧版或框架 | 本次重analysis | 结论 |
|---|---|---|---|
| PERMANOVA | 框架称 stage R² = 0.53（P = 0.001）、focal R² = 0.02（P = 0.003）、互作 R² = 0.02（P = 0.57） | **旧代码里根本没有 PERMANOVA**；本次 4 套模型results见上表（stage 10%–77%、focal 2.4%–9.4%、group/interaction 均不显著） | 旧数字无出处；本次结论方向一致（stage 最强），但具体 R²/P 需按本次results改写 |
| PLS-DA 轴标签 | Figure 2B/2E 硬编码 "Dim 1 (26.00%)"、"Dim 2 (8.00%)" | 实算每成分 R2X = 19.3% / 23.6% | **旧标签and模型不符**，已在图中改为实算值 |
| Figure 2C/2D 的 P 值 | 硬编码 "P = 0.01 / 0.20 / 0.02 / 0.25 / 0.04 / P < 0.0001" | 实算见上（Welch and Wilcoxon results差异较大） | 旧 P 值不可复现；建议正文改用 Wilcoxon + Cohen's d，并notes Pure : Mix = 4 : 12 的不平衡 |
| VIP 前 10 gene | 框架列的名单含 OsPOX223、OsWRKY114、OsAMT2.1、OsGSTU24、OsHXK6 等 | 本次为 OsCPS4-1、OsJAZ5、OsXTR1、OsGH9B9、OsKS4、OsPR1b、OsIPA1、OsWRKY45、OsAOS2、OsCitr_synth（and旧名单仅 3 个重合） | **需作者确认**：旧图 VIP 名单无法用当前模型复现，正文涉及具体 VIP gene的表述要按新名单核对 |
| 弱分离提示 | —— | `orthoI = NA` 时 ropls 交叉验证直接拒绝建模（第一预测成分不显著）；Q²(cum) 仅 0.224 | 建议正文避免把 PLS-DA 分离表述过强 |

## 备注

- **表 04 的gene/pathway口径**：实测为 92 个gene / 20 条pathway；引物面板 96 个位点 = 92 个目标gene + 4 个内参
  （OsAct / OseEF-1α / OseIF-4a / OsUBQ5，其中 OsAct 在stage 0 已剔除）。paper里的"96 个gene / 21 条pathway"
  建议改写为"96 个检测位点（92 个gene + 4 个内参）/ 20 条pathway（另含内参 HKG）"。
- `bioRtools::read_data()` 会把 `R2Y(cum)`、`Q2(cum)` 这类含括号的列名改成 `R2Y.cum.`、`Q2.cum.`，
  下游引用时注意列名变化。

---

## 追加：three_stages扩展（2026-10-18，脚本 06–07）

按作者决定把 Fluidigm 面板的**第三个time点（花期 / flowering）**接入主题 03。判定性analysis见主题 `08_three_stagesexpression`（结论：不影响主线，纳入）。

```bash
/opt/homebrew/bin/Rscript "code/03_overall_gene_expression/06_expression_matrix_build_three_stages.R"
/opt/homebrew/bin/Rscript "code/03_overall_gene_expression/07_PERMANOVA_three_stages.R"
```

- 输入：`data/02_raw_data/04_Fluidigmexpression/00_引物and共享表/fluidigm.normaled.expression.3 times.xlsx`（three_stagesnormalized长表，23368 行）。
- 口径and脚本 01/02 完全一致，仅时期数由 2 变 3；输出为**新编号file**，未覆盖两时期version的 `01_表_expression矩阵.xlsx` / `02_表_PERMANOVA_TableS2.xlsx`。

### 06 expression_matrix_build（three_stages）

| 项目 | results |
|---|---|
| 输出 | `result/03_overall_gene_expression/11_表_expression矩阵_three_stages.xlsx`（6 sheet：gene_annotation / sample_info / combo_matrix / sample_matrix / combo_long / outlier_summary） |
| 规模 | 组合均值矩阵 **48 × 92**（3 期 × 16 组合）；sample级矩阵 **254 × 92**（缺失 0） |
| sample分组 | tillering Mix 58 / Pure 22；flowering Mix 66 / Pure 20；filling Mix 67 / Pure 21 |
| 离群值（stage × focal × neighbor × gene 内判定） | 合计 1411 个标记，占 6.04%（tillering 6.83%、flowering 5.25%、filling 6.09%） |

**一处datacleanup**：three_stages长表里 `OsWRKY82` 有两行gene注释（唯一差别是 `id` 字段带换行残留 `_x000D_`）。脚本 06 先去残留再按 gene 去重为**一gene一行（92 行）**，避免下游 `left_join` 产生 many-to-many 重复行。主题 08 的 `gene_limma` 因此曾多出 1 行（每期 93 行、92 个gene），但**该gene不在任何一期的名义显著名单里，报告的计数and结论不受影响**。

### 07 PERMANOVA（three_stages）

输出 `result/03_overall_gene_expression/12_表_PERMANOVA_three_stages.xlsx`。模型 `~ stage + group + focal + stage:group`，999 次置换，`by = terms` and `by = margin` 各一次，并按时段分层。

| data × 距离 | stage（R²%） | group（R²%） | focal（R²%） | stage:group（R²%） |
|---|---|---|---|---|
| 组合均值 · 列normalized欧氏 | **37.96**（P = 0.001） | 1.77（P = 0.185） | 5.87（P = 0.051） | 3.62（P = 0.151） |
| 组合均值 · Bray–Curtis | **83.37**（P = 0.001） | 0.45（P = 0.275） | 2.49（P = 0.026） | 0.76（P = 0.325） |
| sample级 · 列normalized欧氏 | **13.65**（P = 0.001） | 0.40（P = 0.267） | 1.68（P = 0.013） | 0.77（P = 0.309） |
| sample级 · Bray–Curtis | **60.76**（P = 0.001） | 0.24（P = 0.139） | 1.42（P = 0.004） | 0.41（P = 0.219） |

**分时期分层（by = terms，组合均值）**：tillering group R² = 18.43%（P = 0.019，Bray–Curtis）/ 11.26%（P = 0.112，欧氏）；flowering group R² = 3.62%（P = 0.722）/ 3.20%（P = 0.672）；filling group R² = 3.54%（P = 0.609）/ 7.53%（P = 0.210）。

**and两时期version对照**：结论完全一致——发育时期是主要变异来源（stage 全部 P = 0.001）；种植模式effect小、仅在分蘖期边缘显著（两时期：tillering group R² = 18.4%、P = 0.025；本版 18.43%、P = 0.019）；**`stage:group` 交互在两时期andthree_stages都不显著**。品系（focal）effect量小但对距离口径敏感（P 0.004–0.051）。

### 08 PLS-DA（three_stages）

输出 `result/03_overall_gene_expression/13_表_PLSDA_three_stages.xlsx`（scores / loadings / vip / vip_top10 / variance_explained / axis_tests / cv_summary）。

```bash
/opt/homebrew/bin/Rscript "code/03_overall_gene_expression/08_PLSDA_three_stages.R"
```

- 输入：`11_表_expression矩阵_three_stages.xlsx` 的 `combo_matrix`（**48 行** = 3 期 × 16 组合，Pure:Mix = 4:12 每期）。
- 口径：`ropls::opls(predI = 3, orthoI = 0)`，and两时期脚本一致；**时期整体comparison改用单因素test**（one-way ANOVA、Kruskal–Wallis）＋三对两两 Wilcoxon（BH），因为 t/Wilcoxon 只接受两组；期内 Pure vs Mix test不变。
- **模型statistics**（关键）：**R²X(cum) = 0.463、R²Y(cum) = 0.485、Q²(cum) = −0.023、pR2Y = 0.15、pQ2 = 0.15**；每成分 R2X = 19.1% / 19.0% / 8.2%。
  → **Q² 为负、置换 P = 0.15**，即该 PLS-DA 模型**没有预测能力、未通过置换test**（两时期version Q²(cum) = 0.224，本版更低）。
  **据此建议：PLS-DA 只能作为探索性展示，正文主证据改用 PCA + PERMANOVA（Table S2 / 12 表）。**
- **时期分离**：Dim 1 —— ANOVA P = 0.0036、Kruskal–Wallis P = 6.3e-06，两两 Wilcoxon（BH）全显著（flowering vs tillering P = 1.2e-05、filling vs tillering P = 1.1e-04、filling vs flowering P = 0.021）；Dim 2 —— ANOVA P = 1.0e-07、KW P = 1.7e-06，flowering vs tillering P = 1.0e-08、filling vs tillering P = 1.1e-04、filling vs flowering P = 0.138。
- **期内 Pure vs Mix（Wilcoxon，n = 4 Pure vs 12 Mix）**：Dim 1 —— tillering P = 0.0022（Cohen's d = −1.24）、flowering P = 0.52（d = 0.25）、filling P = 0.0011（d = −1.51）；Dim 2 —— tillering P = 0.078（d = −0.81）、flowering P = 0.107（d = −1.20）。**注意**：这些test是在"为区分组而优化的有监督轴"上做的，属循环论证（审稿意见 M5），不能作为独立证据。
- VIP 前 10（three_stages模型）：OsCPS4-1、OsXTR1、OsGH9B9、OsPR1b、OsWRKY45、OsWRKY102、OsJAZ5、OsAOS2、OsKS4、OsIPA1（and两时期名单 9/10 重合，差异为 OsCitr_synth ↔ OsWRKY102）。

### 09 Figure 2（three_stages，面板级出图）

输出 `result/03_overall_gene_expression/14–18_Figure2A–2E_*_three_stages.pdf / .png`（2A 11×13 in、2B 7.5×5.5 in、2C/2D 8×6 in、2E 9×6 in；PDF 页面尺寸已用 `pdfinfo` 逐张核对一致）。

- 2A 热图：48 个组合 × 92 gene，行内 0–1 缩放，`facet_grid(stage ~ group)` 6 个分面；
- 2B 得分图：6 组（3 期 × Pure/Mix）着色，椭圆用 `type = "norm"`（Pure 组每期仅 4 点，默认的稳健估计报收敛失败）；
- 2C/2D 分布图：标注改为"**stage（Kruskal–Wallis）整体 P + 各期 Mix vs Pure P**"（两时期version标的是 tillering vs filling）；
- 2E 载荷图：three_stages模型的 VIP 前 10。

**注意**：本步只产出**面板级**file；按 2026-10-17 规则，Figure 2 的final组图必须由 VectorCraft 拼版（见 `logs/03_旧版manuscriptoutputs_20261019/03_图表/03_组图方式核查_VectorCraft重做manifest_20261018.md`）。

### 10–11 盆号mappingand盆effectsensitivity（2026-10-18，A 桶：把盆号接进模型）

```bash
/opt/homebrew/bin/Rscript "code/03_overall_gene_expression/10_sample_pot_id_mapping_three_stages.R"
/opt/homebrew/bin/Rscript "code/03_overall_gene_expression/11_pot_id_effect_sensitivity_three_stages.R"
```

对应审稿意见 **M3**："同一 Mix 盆的两株分别作为不同 focal，但共享同一盆环境，应在模型中体现"。

**10 样品 → 盆号mapping**（输出 `19_表_sample_pot_id_mapping_three_stages.xlsx`，sheet：map_long / pot_summary / coverage / sample_level）

- data来源：三批样品的**标签表**——第一批 `第一批96个样品选择表.xlsx`（`label` 形如 `2(15_15)`，`tag` 即analysis用样品 ID）；第二/三批的孔位表（`XL1第二批RNA样品编号.xlsx`、`samples.xlsx`），每 3 行为一组（盆标签行 / `Focal` 样品 ID 行 / `Neighbor` 样品 ID 行）。取值列按"标签行里含 `(` 的列"动态判定（两个file列偏移不同）。
- 2026-10-18 已把两份上游标签表复制进项目只读副本目录：`data/02_raw_data/04_Fluidigmexpression/02_第二批/XL1第二批RNA样品编号.xlsx`（上游 `/Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/field试验/fluidigm/第二批/`）、`…/03_第三批/samples.xlsx`（上游 `…/fluidigm/第三批/`）。
- **results：254/254 个analysissample 100% 匹配到盆号**。盆结构要分"标签表"and"进入analysis的sample"两层看（2026-10-18 复核）：

  | 时期 | 标签表：Pure 盆 | 标签表：Mix 盆 | analysis集：Pure 盆（sample/盆） | analysis集：Mix 盆（sample/盆） |
  |---|---|---|---|---|
  | tillering | 22 盆（每盆 2 个样品 ID） | 51 盆（每盆 2 个样品 ID） | 20 盆：18×1、**2×2** | 48 盆：38×1、**10×2** |
  | flowering | 40 盆 × 2 | 60 盆 × 2 | 20 盆：**全部 1** | 35 盆：4×1、**31×2** |
  | filling | 40 盆 × 2 | 60 盆 × 2 | 21 盆：**全部 1** | 35 盆：3×1、**32×2** |

  - 即：**Fluidigm 板对 Pure 盆同样record了 Focal / Neighbor 两个样品 ID**（例如第三批 Pure 盆 `76(69_69)` 的两个 ID 是 401 and 402），但**进入本次three_stagesanalysis的sample只取了 Focal 那一株**（402、408、424、430 等 Neighbor ID 均不在analysis集内；analysis集里 Pure 组合的样品 ID 全为奇数，即 Focal 行）。
  - Mix 盆的情况相反：同一盆的两株分属互为 focal 的两个组合，因此两株**都会进入analysis**（flowering 31 盆、filling 32 盆、tillering 10 盆）。
  - **待作者确认**：design上 Pure 盆是否应当"两株都纳入analysis"？若是，则需把 Neighbor 行的样品也接进analysis矩阵（第三批 Neighbor ID 存在于标签表，但当前analysis输入表只含 Focal 行）；这会改变 Pure 组的sample量（每期 Pure 从 20–22 增至约 40），进而影响全部expressionstatistics，属**口径变更**，故本轮先如实record、未擅自纳入。

**11 盆随机截距sensitivityanalysis**（输出 `20_表_pot_id_effect_sensitivity_three_stages.xlsx`）

- 两种口径逐gene对照（每期 92 个gene）：① 现行 `lm(log2(expr+1) ~ group + focal)`（不含盆）；② `lmerTest::lmer(… + (1|pot))`（含盆随机截距）。BH-FDR 分别correction。
- **结论：加盆随机截距后，三期的 FDR 结论完全不变——两版都没有gene FDR < 0.05。**名义 P < 0.05 的gene数变化很小：tillering 14 → **12**、flowering 2 → 2、filling 1 → **0**。
- 盆方差占比（ICC，限"检出率 ≥ 50%"的gene，n = 68 / 70 / 71）：中位数 **0.434（tillering）/ 0.000（flowering）/ 0.020（filling）**，最大 0.998 / 0.701 / 0.996。
- **两点必须写进方法/局限**：① 盆随机截距的**奇异拟合比例很高**（17/92、40/92、37/92），原因是多数盆内只有 1 个sample（tillering 56/68 盆仅 1 株），盆方差识别不稳；② ICC 对**低检出gene**会假性接近 1（残差方差趋近 0），因此上表的 ICC 只取检出率 ≥ 50% 的gene。
- 因此：**现有结论（expression层无gene通过 FDR、effect以分蘖期可检出为主）在计入盆effect后依然成立**，可作为对 M3 的直接回应；ICC 的具体数值建议只作辅助陈述。

### 12 per_gene_variance_decomposition（three_stages）

输出 `result/03_overall_gene_expression/21_表_per_gene_variance_decomposition_three_stages.xlsx`（gene_variance / terms_summary / per_stage_group / summary_by_stage / top_by_stage）。

- 方法：响应 `log2(expression + 1)`；全data模型 `lm(y ~ stage * group + focal)` 用 `car::Anova(type = 2)` 偏平方和换算成**偏 R²**；另按时期拟合 `lm(y ~ group + focal)` 给出该期内 group 的偏 R²。
- **results（254 sample，92 gene的偏 R² 中位数）**：**stage 15.1%**、group **0.18%**、stage:group 0.57%、focal 1.49%、残差 81.9%。
- **分时期 group（邻株effect）偏 R²**：中位数 **tillering 2.86%（13 个gene > 5%，最大 OsGNS1 10.8%、OsARF1 10.2%）**、flowering 0.44%（0 个 > 5%）**、filling 0.49%（1 个 > 5%，OsCitr_synth 6.5%）**。
- 结论：**邻株effect对单geneexpression方差的解释力很小（中位 < 3%），且集中在分蘖期**——and"分蘖期可检出、随后减弱"的主线一致，正好回答审稿意见 M5 的"邻株effect各期解释多少变异"。偏 R² and PERMANOVA 的 R² 口径不同，需分别陈述。

### 13 outlier_sensitivity（three_stages）—— 一处需要作者拍板的重大口径问题

输出 `result/03_overall_gene_expression/22_表_outlier_sensitivity_three_stages.xlsx`（outlier_cells / outlier_summary / limma_compare / compare_summary）。

```bash
/opt/homebrew/bin/Rscript "code/03_overall_gene_expression/13_outlier_sensitivity_three_stages.R"
```

- **剔除规则**（and主题 03 脚本 01 的组合均值矩阵一致）：`bioRtools::find_outliers()`，按 **stage × focal × neighbor × gene** 分组判定；被标记比例 **tillering 6.83% / flowering 5.25% / filling 6.09%**。被标记最多的gene：OsHXK2（26）、OsDREB1B（25）、OsAMT2.1 / OsPR5 / OsRBBI2 / OsTP24（各 24）。
- **现行主statistics（sample-level limma）用的是未剔离群值的complete矩阵**；sensitivityanalysis把这些单元格置为 NA 后用**同一个 limma 模型**重跑：

  | 时期 | 名义 P<0.05（剔前 → 剔后） | **FDR<0.05（剔前 → 剔后）** | logFC 相关（Pearson） | P 值相关（Spearman） |
  |---|---|---|---|---|
  | tillering | 14 → **19** | 0 → **9** | 0.969 | 0.570 |
  | flowering | 1 → **9** | 0 → 0 | 0.568 | 0.033 |
  | filling | 1 → 2 | 0 → 0 | 0.751 | 0.462 |

- **剔除离群值后，分蘖期有 9 个gene通过 FDR < 0.05**：OsGNS1（FDR = 0.014）、OsIAA1（0.020）、OsGH3.2（0.031）、OsARF1（0.032）、OsHDR（0.036）、**OsBLE1（0.036，logFC = +0.449）**、OsPyr_dehydra（0.036）、OsGBSSII（0.047）、OsXTH17（0.047）。
- **这意味着**：manuscript现在写的"**经过 FDR correction没有任何gene存在差异**"是**口径依赖**的结论——如果采用"剔除离群值"的口径（这正是组合均值矩阵、Figure 2A/3A/4A 热图所用的口径），分蘖期就有一批通过 FDR 的gene，其中包含paper主线关心的 OsBLE1、OsGBSSII、OsGNS1 and生长素gene。**两种口径必须由作者选定其一为主口径**，因为正文、图注、摘要都要跟着改。
- 另需注意：`find_outliers` 在每个 (stage × focal × neighbor × gene) 只有 4–6 个sample时判定较激进（剔掉约 7% 的单元格），所以"剔后 9 个 FDR"也要谨慎解读；建议在方法里写清规则and标记比例，并把两版results都放进补充材料。

### 14 组合均值矩阵的"不剔除"口径and稳健聚合（2026-10-18，作者决定）

**决定**：作者 2026-10-18 明确"**先不剔除**"——**sample级statistics（limma/camera/PERMANOVA/PLS-DA）本来就是用未剔除的complete矩阵**，因此statistics口径不变；需要跟着改的是**图侧的组合均值矩阵**（`11_表` 的组合均值是剔除后算的，and_statistics口径不一致，正是审稿意见 M6"离群值"条目的要害）。

```bash
/opt/homebrew/bin/Rscript "code/03_overall_gene_expression/14_expression_matrix_build_three_stages_without_outlier_removal.R"
```

输出 `result/03_overall_gene_expression/23_表_expression矩阵_three_stages_without_outlier_removal.xlsx`（7 sheet：gene_annotation / sample_info / **combo_matrix** / **combo_matrix_median** / sample_matrix / combo_long / combo_diff、combo_diff_median）。

| 对照 | 逐格差异 | 最大绝对差 | 备注 |
|---|---|---|---|
| 不剔除 **均值** vs 剔除后均值 | **1305 / 4416 = 29.55%**（其中 1218 格相对差 > 10%） | **1366.9** | 极端例子：`flowering_69_15 × OsHDS` 剔除后 0.031 vs 不剔除 1366.9（4 个数量级） |
| 不剔除 **中位数** vs 剔除后均值 | 3789 / 4416 = 85.8% | 502.2（但**中位差仅 0.006**） | 差异遍在但都很小，不会被单个极端sample带飞 |

**建议（供作者采纳）**：图侧组合均值矩阵改用**中位数**——同样"不剔除任何data"，但热图不会被个别极端sample拉成二值图（旧版 Fig.2A"分蘖期 Mix 一列把数十个gene拉到 1、其余压到蓝色"正是这个毛病，审稿意见也点了）。若坚持用**不剔除的均值**，则 29.55% 的格子会变、个别格子差 4 个数量级，热图可读性差。

**下游需按new_definition重跑（图侧）**：Figure 2A 热图、Figure 3A 热图、PLS-DA（组合均值输入）、以及任何用组合均值的plotting；sample级statistics不受影响。
