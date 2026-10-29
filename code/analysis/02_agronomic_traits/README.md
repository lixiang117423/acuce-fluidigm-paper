# 02_agronomic_traits — 运行record

室外全生育期实验的agronomic_traits重analysis：height_change_index、tiller_number、单株yield、表型 PCA and欧氏距离、
遗传距离—表型相关性，另补充fieldplant_heightandthousand_grain_weight。对应 Figure 1A–1F and Figure S4。

## 运行环境

- R：`/opt/homebrew/bin/R`（`R version 4.6.1 (2026-06-24)`），运行方式 `/opt/homebrew/bin/Rscript "<脚本绝对路径>"`
- 关键包version：tidyverse 2.0.0、ggplot2 4.0.3、lme4 2.0.6、lmerTest 3.2.1、emmeans 2.0.4、
  multcomp 1.4.32、rstatix 1.1.0、bioRtools 1.95.1、readxl 1.5.0.1、writexl 2.0.1
- 固定随机种子 `set.seed(12345)`

## 脚本and运行命令（2026-09-22）

```bash
cd "<repository_root>"
/opt/homebrew/bin/Rscript "code/02_agronomic_traits/01_height_change_index_and_statistics.R"
/opt/homebrew/bin/Rscript "code/02_agronomic_traits/02_tiller_and_yield_statistics.R"
/opt/homebrew/bin/Rscript "code/02_agronomic_traits/03_phenotype_PCA_and_euclidean_distance.R"
/opt/homebrew/bin/Rscript "code/02_agronomic_traits/04_genetic_distance_and_phenotype_correlation.R"
/opt/homebrew/bin/Rscript "code/02_agronomic_traits/05_agronomic_traits_plotting_Figure1.R"
/opt/homebrew/bin/Rscript "code/02_agronomic_traits/06_field_height_and_thousand_grain_weight_supplement.R"
```

6 个脚本退出码均为 0；日志 `result/02_agronomic_traits/00_运行日志_*.log` / `00_运行错误_*.log`。
（04 依赖 03 的输出，05 依赖 01–04 的输出，须按上述顺序运行。）

## 输入anddata口径

| 输入 | notes |
|---|---|
| `data/processed_inputs/03_agronomic_traits_height_tiller_yield.xlsx` | 8095 行 × 6 列；plant_height 15 个time点（20220711–20221031）、tiller_number（20220927）andyield（20221227）各 1 个时点 |
| `data/02_raw_data/05_遗传距离/月亮谷纯系的遗传距离.xlsx` | 99 个纯系的遗传距离矩阵（取出 A15/A49/A69/A87 的 4 × 4 子矩阵） |
| `data/processed_inputs/02_field_disease_survey.xlsx` | 取 Category = PH（3564 行）and KGW（1140 行）用于补充analysis |

**必须记住的三点data事实**（实测核对 2026-09-22）：

1. plant_height表 2518 行 `time` 为空，analysis一律剔除；**focal 36 / 74 的plant_heightrecord time 全部为空**，故不进入本主题任何analysis。
2. Pure 组的 `neighbor` 已统一等于 `focal`（stage 0 prepare时设定），因此纯种均值按 `focal × time` 聚合即可。
3. tiller_numberandyield为单时点：纯种 20 株、每个混种组合 10 株（yield个别组合 7–9 株），**不存在批次/区组等随机结构**，
   故主analysis用单因素 ANOVA + Tukey，而不是混合模型。

## 方法

| 脚本 | 方法 |
|---|---|
| `01` | height_change_index `index = value − 同 focal × time 的 Pure 均值`；总体test = 混种 index 对 0 的单sample t test（另做 Wilcoxon）；每个 focal 同时用旧口径 `bioRtools::anova_posthoc()`（ANOVA + Tukey 字母）and `lmerTest::lmer(index ~ neighbor + (1 \| time))` + `emmeans`(Sidak) 两套并列；LMM 的随机截距为time点 |
| `02` | tiller_number / yield：每个 phenotype × focal 的 ANOVA + Tukey 字母；另用 `lm(value ~ focal * neighbor)` 给整体test |
| `03` | 三性状组合均值（16 个 focal × neighbor 组合）→ Z-score 后 `prcomp`；欧氏距离在normalized性状空间按"混种组合到该 focal 纯种点"计算 |
| `04` | 遗传距离and三性状的 Pearson 相关（主analysis = 12 个混种组合；另存含纯种的 16 组合version对照旧图），Bonferroni correction（α = 0.05/3） |
| `05` | Figure 1A–1F；`bioRtools::theme_prism()`、`write_data(units = "in", dpi = 500)`；P 值统一 `P < 0.001` / `P = 0.012` 定点写法 |
| `06` | field PH / KGW：`lmerTest::lmer(Value ~ condition + FocalPlant + (1 \| Year))` + `emmeans` 各纯系均值and 95% CI |

注意：`03` 中必须给矩阵设行名（组合标签），否则 `scale()`/`prcomp()` 的 scores 会丢失组合标识，
导致 Figure 1D 的点and连线全部落到 NA（本次已修复并复跑验证）。

## 输出

| file | 内容 |
|---|---|
| `01_表_height_change_index.xlsx` | 6 sheet：index_long、pure_mean、overall_test、anova_posthoc、lmm_emmeans、lmm_vs_ctrl |
| `02_表_tiller_and_yield_statistics.xlsx` | 3 sheet：anova_posthoc、descriptive、two_way_anova |
| `03_表_表型PCA.xlsx` | 4 sheet：combo_means、variance_explained、scores、loadings |
| `04_表_表型欧氏距离.xlsx` | 2 sheet：distance_long、focal_summary |
| `05_表_genetic_distance_and_phenotype_correlation.xlsx` | 3 sheet：genetic_distance_submatrix、data_joined、correlation |
| `06–11_Figure1A–1F_*.pdf` / `.png` | 1A height_change_index、1B tiller_number、1C 单株yield、1D 表型 PCA、1E 三性状整体comparison（10 × 5 in）、1F 遗传距离and表型（10 × 4.5 in）；1A–1D 为 8 × 6 in |
| `12_表_fieldplant_heightandthousand_grain_weight.xlsx` | 3 sheet：overall_model、emmeans、descriptive |
| `12_FigureS4_fieldplant_heightandthousand_grain_weight.pdf` / `.png` | 6 纯系 × Pure/Mix 的估计均值and 95% CI（10 × 5 in） |

## results要点（2026-09-22）

- **plant_height**：混种整体促进plant_height，变化指数均值 +0.789 cm（95% CI 0.326–1.252；单sample t test P = 0.00084，Wilcoxon P = 6.5e-05）；
  Pure vs Mix 两sample t test P = 0.0060（见 Figure 1E），andmanuscript P = 0.006 一致。
- **tiller_number**：整体不显著（P = 0.178），特定组合显著（见 `02_表_tiller_and_yield_statistics.xlsx`）。
- **yield**：整体不显著（P = 0.993）。
- **表型 PCA**：PC1 = 49.31%、PC2 = 38.85%（and旧版 Figure 1D 硬编码的轴标签完全一致，可作为 PCA 复现的直接证据）。
- **互作影响强度（normalized性状空间欧氏距离，混种到纯种）**：87（2.54）> 49（2.24）> 15（2.00）> 69（1.44），
  即 **87 受影响最大、69 最稳定**。
- **遗传距离—表型**：全部不显著（Mix 组合 n = 12：plant_height r = −0.12、P = 0.699；tiller_number r = −0.09、P = 0.772；
  yield r = −0.33、P = 0.298；Bonferroni correction后 P ≥ 0.89）。阴性results稳固，支持"互作effect不能由遗传距离预测"。
- **field补充（新）**：混种显著提高fieldplant_height（+0.98 cm，P = 0.0115）andthousand_grain_weight（+0.55 g，P = 4.6e-08），
  and室外实验"互作促进生长"的方向一致，且为独立year/独立性状的旁证。

## and旧版 / manuscript的交叉核对

| 项目 | 旧版或manuscript | 本次重analysis | 结论 |
|---|---|---|---|
| 表型 PCA 方差解释率 | Figure 1D 轴标签：PC1 49.31%、PC2 38.85% | PC1 0.4931、PC2 0.3885 | **完全一致** |
| plant_height整体促进 | manuscript P = 0.006 | 单sample t test P = 0.00084；Pure vs Mix t test P = 0.0060 | 方向一致；**manuscript的 0.006 对应 Pure vs Mix 两sampletest**，已在 README record口径 |
| 遗传距离—表型 | 无显著相关（阴性） | 三项均不显著 | **一致** |
| Figure 1F 口径 | 含 4 个纯种点（距离 = 0），共 16 点 | 主analysis改为 12 个混种组合；16 组合version同表另存 | 口径已改进，两者均已存档可比 |
| 互作影响强度排序 | paper框架称"15 受互作影响最大（以 69 为邻体最强）、87 次之、69 最稳定" | 87 > 49 > 15 > 69 | **不一致，需作者确认**：按normalized性状空间的欧氏距离，受影响最大的是 87 而非 15；"69 最稳定"一致 |
| 旧"PCAresults欧氏距离.xlsx" | 曾被当作表型距离 | 实际是**expression空间**的 3 时期 × 4 纯系距离矩阵（tillering/flowering/grain-filling） | 两者不是同一量，**不能互相替代或对照**；已在 `项目filerecord.md` 标注 |

## 备注

- `bioRtools::anova_posthoc()` 返回的字母列为 `tukey_signif`（不是 `signif`），且其输出的 ANOVA P 值为整组共用值；
  本主题按旧口径沿用该函数以保证and旧图可比，LMM results作为稳健性对照另存。
- emmeans 在观测数 > 3000 时会自动关闭自由度计算（列名变为 `asymp.LCL`/`asymp.UCL`），`06` 脚本已做列名归一化处理。

---

## 追加：planned_contrasts_and_raw_scale_effects（2026-10-18，脚本 07，A 桶第 2 项）

```bash
/opt/homebrew/bin/Rscript "code/02_agronomic_traits/07_planned_contrasts_and_raw_scale_effects.R"
```

对应审稿意见 **M4**（statistics方法前后不一致 + effect量尺度）and **M10**（原始尺度effect大小and 95% CI）。

- **改动要点**：① tiller_number / 单株yield由"各邻体两两 Tukey"改为**计划对比**——每个 Mix 邻体组合 vs **该 focal 自己的单作对照**（emmeans 对比向量 + Dunnett correction）；② height_change_index的test单元由"每行（15 调查 × 多株，伪重复）"改为**先按每次调查取均值、再以 15 个time点为单元**做 vs 0 的test（纵向模型的推断结论仍以主题 06 的 Figure S9/S12 为准）；③ inoculationdata用and主题 06 一致的盆级模型 `lesion ~ mode + (1|batch) + (1|pot)`，把effect换算成**原始尺度差值and百分比变化 + 95% CI**。
- **输出**：`result/02_agronomic_traits/13_表_planned_contrasts_and_raw_scale_effects.xlsx`（4 sheet：planned_contrasts_traits / height_index_contrasts / disease_effects / raw_scale_summary）。
- **关键results（原始尺度 + 百分比变化）**：
  - **胡麻叶斑病（盆级）**：−1.195 个病斑/cm²（95% CI −2.183 ~ −0.208），**−8.87%**（95% CI −16.19% ~ −1.54%），P = 0.018（n = 1043 叶 / 96 盆）——这正是审稿人要的"原始尺度百分比变化"。
  - **稻瘟病（盆级）**：−0.183 个病斑/cm²（−0.992 ~ +0.625），−5.43%（−29.36% ~ +18.50%），P = 0.65（n = 616 叶 / 64 盆）。
  - **单株yield**：focal 15 × 邻体 69 → +4.63 g/株（0.69 ~ 8.56），**+60.0%**（+8.9% ~ +111.1%），Dunnett P = 0.017（BH-FDR = 0.41）；其余组合均不显著。
  - **tiller_number**：focal 15 × 邻体 69 → +3.15 个/株（0.08 ~ 6.22），+23.6%，Dunnett P = 0.043（BH-FDR = 0.41）；其余不显著。
  - **height_change_index（按 15 次调查为单元）**：最强的是 **focal 87 × 邻体 15 → −8.23 cm**（95% CI −10.84 ~ −5.61，P = 4.6e-21），focal 49 × 87 → +6.07 cm、focal 49 × 69 → +7.06 cm；方向性and主题 06 的纵向模型（速率effect +0.027 cm/天）需在正文一并交代。
- **图注/正文要跟上的三处**：① Figure 1B/1C 的字母若仍用两两 Tukey，需and正文的计划对比分开表述；② tiller_numberandyield的**显著性只在 Dunnett 口径下有个别组合**（BH correction后均不显著），正文不宜写成"整体显著"；③ plant_height指数图的误差and显著性必须标"以 15 次调查为单元"，不能再用行数当 n。

---

## 追加：Figure 1 按new_definition重绘（2026-10-18，脚本 08）

```bash
/opt/homebrew/bin/Rscript "code/02_agronomic_traits/08_agronomic_traits_plotting_Figure1_new_definition.R"
```

- **输出（面板级，`result/02_agronomic_traits/`）**：
  | 编号 | 内容 | 尺寸（PDF 已核对） |
  |---|---|---|
  | `14_Figure1A_height_change_index_new_definition` | height_change_index（cm），**以 15 次调查为单元的 95% CI** + 各组合 P | 12×4.8 in（864×345 pt） |
  | `15_Figure1B_tiller_number_new_definition` | tiller_number原始尺度箱线 + 抖动散点，标注**and自身单作的计划对比 P**（`*` = BH-FDR < 0.05） | 12×4.8 in |
  | `16_Figure1C_单株yield_new_definition` | 单株yield（g/株），同上 | 12×4.8 in |
  | `17_Figure1D_表型PCA_new_definition` | 三性状 PCA（16 组合；Pure/Mix 用shape区分、focal 用颜色，虚线指向各 focal 重心） | 8×6 in（576×432 pt） |
  | `18_Figure1E_三性状整体comparison_new_definition` | 三性状 Pure vs Mix（**每个分面的条带带单位**：cm / per plant / g per plant） | 13×5.2 in（936×374 pt） |
  | `19_FigureS_遗传距离and表型_new_definition` | 遗传距离—表型（**按 M12 移出main_figures**，作探索性supplementary_figures，图注写明 n = 12、非独立、Bonferroni 后不显著） | 11×4.5 in（792×324 pt） |
- **改动对照（对应审稿意见）**：① 1A 的误差棒/显著性由"每行当 n"改为 15 次调查为单元（M3/M4）；② 1B/1C 由两两 Tukey 字母改为**计划对比 P**（M4）；③ 各分面complete单位（图规范条目）；④ 配色由图内mapping直接generate图例，不再用文字描述颜色（图例and图内颜色不符的条目）；⑤ 遗传距离—表型移出main_figures（M12）。
- **仍待同步的两处（属manuscript侧，D 桶）**：① 方法and_figure_legends要把"PCA 输入只含 3 个agronomic_traits、不含病斑密度"写清楚（M11）；② 正文对tiller_number/yield的表述要按 BH-FDR results收紧。
