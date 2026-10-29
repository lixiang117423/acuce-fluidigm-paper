# 06_design_decomposition_and_robustness — 运行record

本主题**不属于旧版流程**，是重analysis时新加入的独立analysis层：针对旧版在实验design、推断单位and_statistics方法上的
薄弱环节，用不同的analysis思路重新提问，再and旧results对照。

## 五个独立问题（旧版未涉及）

| 编号 | 问题 | 旧版做法 | 本次做法 |
|---|---|---|---|
| Q1 | inoculation实验里同一盆的多株植物能否当独立观测？ | 直接按株建模（df ≈ 600–1000） | 加盆（image）随机截距；并用盆聚合值复核 |
| Q2 | 邻体身份本身贡献多少？互作是否非加性？ | 按 focal 逐组切片comparison，无整体decomposition | `focal + neighbor + focal:neighbor` decomposition，看主effectand交互 |
| Q3 | 同一对品种的双向响应是否对称（reciprocity）？ | 未涉及 | 12 个无序品种对，comparison i→j and j→i 的effect（相关 + 对称性置换test） |
| Q4 | plant_height互作是速率effect还是水平effect？effect何时出现？ | 15 个time点压成 1 个指数做一次 ANOVA | 纵向模型（group × day）+ 逐time点test + FDR 定 onset |
| Q5 | pathway水平该怎么test？ | 逐genetest + 卡方/Fisher（每pathway 2–3 个gene） | pathway整体分数（pathway内gene Z 均值）+ 品系随机截距 + 置换test |

## 运行环境and命令（2026-09-22）

R version：`/opt/homebrew/bin/R`（4.6.1）；关键包：tidyverse 2.0.0、ggplot2 4.0.3、lme4 2.0.6、
lmerTest 3.2.1、emmeans 2.0.4、vegan 2.7.6、patchwork 1.3.2、bioRtools 1.95.1；固定 `set.seed(12345)`。

```bash
cd "<repository_root>"
/opt/homebrew/bin/Rscript "code/06_design_decomposition_and_robustness/01_inoculation_伪重复anddesigndecomposition.R"
/opt/homebrew/bin/Rscript "code/06_design_decomposition_and_robustness/02_plant_height_growth_curve_and_effect_onset.R"
/opt/homebrew/bin/Rscript "code/06_design_decomposition_and_robustness/03_transcriptome_pathway_scores_and_shared_accessions.R"
```

3 个脚本退出码均 0；日志见 `result/06_design_decomposition_and_robustness/00_运行日志_*.log`。

## results：新发现and对旧结论的correction

### 01 inoculation实验（Figure S8）

• Q1 盆内clusteringcorrection（含 2026-09-22 的一处更正）：采样结构经用户确认并用data复核为
**"每一株就是一个个体，一株只取一片叶"**，表中一行 = 一株植物，image 内 leaf_number 唯一、
一个 image 只对应一个 focal × neighbor 组合（每 image 6–14 株）。**因此不存在"同一株重复采叶"式的伪重复**
（本 README and脚本早期对此的"伪重复"表述已作废，见文末更正notes）。真正的问题是：同一盆内的多株植物
共享土壤and邻体环境、并非独立观测——实测胡麻叶斑病的盆间 ICC = 0.12，确实存在clustering。
加盆随机截距后，胡麻叶斑病的 P 值从 **0.00021**（按株，df ≈ 1039）变为 **0.018**（df ≈ 92），
按盆聚合复核 P = 0.021；稻瘟病从 0.491 变为 0.652（按盆聚合 0.646）。
结论方向不变（混种降低胡麻叶斑病病斑密度），但要理解成**有效独立单位是盆而不是株**：
1043 株其实只来自 96 盆，旧版按株计 df 会高估精度，effect证据强度比旧版弱约两个数量级。

• Q2 designdecomposition（新发现）：稻瘟病的 focal（F = 8.75，P = 1.2e-5）、neighbor（F = 5.40，P = 0.0012）and
**focal:neighbor 交互项（F = 4.19，P = 0.00098）**均显著；胡麻叶斑病则是 focal 主effect主导
（F = 16.1，P = 3.5e-10）、neighbor 较弱（F = 3.81，P = 0.0099）、交互项不显著（P = 0.35）。
即：稻瘟病表现为**组合特异的非加性互作**，胡麻叶斑病近似**加性**。这是"互作effect组合特异"的直接statistics证据，
同时notes两种病害的互作结构不同。

• Q3 reciprocity（新发现）：12 个无序品种对中，i 受 j 影响and j 受 i 影响既不相关也不对称——
稻瘟病 r = 0.52（P = 0.29）、平均绝对effect差 0.28（对称性置换 P = 0.90）；胡麻叶斑病 r = −0.36（P = 0.49）、
平均绝对effect差 0.33（置换 P = 0.88）。互作是**单向、不对等**的，一对品种中谁受影响更大and另一方无关。

### 02 plant_height纵向analysis（Figure S9）

• 速率 vs 水平：`group` 主effect P = 0.763（第 0 天无水平差异），而 **`group:day` = +0.0266 cm/天（t = 2.49，P = 0.0129）**，
notes互作改变的是**生长速率**（约比纯种快 4%），不是整体抬高水平。分 focal 只有 49 显著（+0.048 cm/天，P = 0.028），
15 / 69 / 87 方向一致但未达显著。

• effectonsettime：15 个time点中**只有第 0 天（2022-07-11）在 FDR 后显著**（混种高 1.83 cm，P = 0.0032，FDR = 0.048）；
第 14 天原始 P = 0.016（FDR = 0.081）；其余 13 个time点均不显著。旧版"互作整体显著促进plant_height（P = 0.006）"应改写为：
**effect在移栽后第一次调查即出现、幅度小（全季累计约 3 cm）、随后被个体间变异淹没**。

### 03 pathway整体分数（Figure S10）

• 方法：pathway内gene按时期做 Z-score，取pathway_scores（均值），再用 `lmer(score ~ group + (1|focal))` test
（品系随机截距以处理 12 个混种组合shared_accessions导致的非独立性），并用 999 次标签置换and BH-FDR 交叉验证。

• results：20 条pathway × 2 个时期**均未达到 FDR < 0.05**；最小原始 P = 0.044（灌浆期 Momilactone，仅 1 个gene；
置换 P = 0.050），分蘖期最小 P = 0.106。但effect方向清晰：**分蘖期以正向为主**（Auxin +0.26、High sugar +0.20、
Defense +0.14、Brassinosteroid +0.09 Z 单位，P = 0.17–0.56），**灌浆期多为负向**（Defense −0.19、BR −0.15、
Jasmonic acid −0.13、Strigolactone −0.06）。

• 结论：换成检出力更高的pathway_scores法后仍无显著pathway，notes在本data的sample量下**pathway水平的差异无法达到statistics显著**；
正文应改为报告effect量and方向（分蘖期上调趋势、灌浆期回落趋势），不能写"某pathway显著富集"。

## and旧版结论的对照

| 结论 | 旧版 | 本次独立analysis | 变化 |
|---|---|---|---|
| 混种降低胡麻叶斑病 | P = 0.0002，作为强证据 | pot 级 P = 0.018–0.021 | 方向不变，强度大幅下调 |
| 稻瘟病无差异 | P = 0.491 | pot 级 P = 0.65 | 一致 |
| 互作组合特异 | 由分 focal 字母间接推断 | 稻瘟病交互项 P = 0.00098（直接证据）；胡麻叶斑病无交互 | 新增机制层面结论：两种病害不同 |
| 双方响应关系 | 未考虑 | 无 reciprocity（P = 0.88–0.90） | 新增：互作不对等 |
| plant_height整体促进生长 | P = 0.006（水平effect） | group:day P = 0.013（速率effect）；仅第 0 天 FDR 显著 | 从水平effect改写为早期速率effect |
| pathway显著富集 | 框架列出 BR / AUX 等 P 值 | 卡方 / Fisher andpathway_scores法均无 FDR < 0.05 | 旧数字不可复现，须改为effect量描述 |

## 三个技术教训（已写入脚本注释）

• `aplot::insert_right()` 会强制对齐两张子图的 y 轴，导致 y 值超出另一图范围的点被**静默删除**
（实测 12 个点被删 7 个）；量纲不同的图拼接请改用 `patchwork::wrap_plots()`。

• 系数表里若存在名为 `df` 的列，在其上继续 `mutate(df$xxx)` 会被**data掩码**解析成该列
（报错 `$ operator is invalid for atomic vectors`）；需要先把sample量、组数取成普通标量。

• `group` 若以字符向量给出，`lmer` 按字母序把 `Mix` 设为参照水平，系数名会变成 `groupPure`；
必须显式 `factor(..., levels = c("Pure", "Mix"))`，否则按 `groupMix` 取值会下标越界。

---

## 追加analysis（Q6–Q9，2026-09-22 续做）

```bash
/opt/homebrew/bin/Rscript "code/06_design_decomposition_and_robustness/04_transcriptome_partner_effect_variance_decomposition.R"
/opt/homebrew/bin/Rscript "code/06_design_decomposition_and_robustness/05_plant_height_effect_shape_and_change_points.R"
/opt/homebrew/bin/Rscript "code/06_design_decomposition_and_robustness/06_expression_empirical_Bayes_differential_test.R"
/opt/homebrew/bin/Rscript "code/06_design_decomposition_and_robustness/07_field_year_and_field_stability.R"
```

4 个脚本退出码均 0；`04` 含 99 次置换、`05` 含 199 次自助法，均设 `set.seed(12345)`。

### Q6 transcriptome的partner_effectdecomposition（Figure S11）

• 方法：先用交叉随机effect模型 `(1|focal) + (1|neighbor) + (1|combo)` 做方差成分decomposition，
results **184 个模型中 174 个奇异拟合**、残差占比 96%–97% —— 在本sample量下方差成分不可识别。
改用两因素方差analysisdecomposition（响应 = log2(expression + 1)，Type III SS，逐gene）：
`y ~ focal * neighbor`，输出三项 partial η²（自身gene型 / 伙伴gene型 / 特定组合）。

• 逐gene FDR < 0.05 的gene数：分蘖期 focal 2 个、neighbor 0 个、交互 1 个（OsHDS，η² = 0.37，FDR = 0.029）；
灌浆期三项均为 0。

• **关键校准（重要）**：交互项的 η² 中位数看似最大（分蘖期 0.123、灌浆期 0.084，高于两个主effect），
但交互项自由度为 9（主effect各 3），直接comparison会偏向交互项。用 99 次置换打乱配对关系（neighbor 在 focal 内重排）做参照后：
观测值 0.123 vs 零分布 0.120（P = 0.42，分蘖期）、0.084 vs 0.104（P = 0.97，灌浆期）——
**交互项and"纯加性"下的期望水平没有差别**。

• 结论：transcriptome层面**没有可检出的特定组合（SCA）effect**，and主题 06 Q2 中"稻瘟病表型存在显著 focal × neighbor 交互"形成对照：
非加性effect在病害表型上可见，但在 92 个gene的expression水平上检不出。

### Q7 plant_heighteffect的shapeandchange_points（Figure S12）

• 方法：植物水平模型 `value ~ group * post + combo`（post = I(day ≥ tau)），对 14 个候选change_points做网格搜索，
似然比给 95% 区间；再comparison恒定差 / 线性差 / 阶跃差三种shape的 AIC；199 次分层自助法给出不确定性。

• 最优change_points tau = 36 天（2022-08-16），95% 似然比集合仅含 36；该点的阶跃effect +1.31 cm（P = 0.171，不显著）。

• **shapecomparison（决定性）**：线性差模型 AIC = 40273（权重 100%），阶跃差模型 42947，恒定差模型 48191；
199 次自助法中 **100% 选择线性差模型**，仅 13.6% 的自助sample给出显著阶跃。

• 结论：Mix and Pure 的plant_height差随time是**线性变化**（and Q4 的 `group:day` = +0.0266 cm/天一致，全季累计约 3 cm），
既不是"某时点之后恒定抬高"，也不是"全程恒定差"。此前"change_points在 36 天"只能理解为分段拟合的伪最优，
真正的信号是斜率差异。

### Q8 经验贝叶斯（limma）替代 Wilcoxon（Figure S13）

• 方法：响应 = log2(expression + 1)，design `~ 0 + group + focal`（focal 作固定区组），`limma::eBayes(robust = TRUE)`。

• results（**2026-09-22 更正后**）：改用**显式 Mix − Pure 对比矩阵**后，**分蘖期and灌浆期均无gene达到 FDR < 0.05**；
分蘖期有 14 个gene P < 0.05（未correction，最强 OsGNS1 logFC = 1.15、P = 0.0024），灌浆期仅 1 个。
方向一致性（and Wilcoxon 口径comparison）：分蘖期 79.3%、灌浆期 67.4%。

> **⚠️ 更正notes（重要）**：本 README and CHANGELOG 早期version写"分蘖期 53/92、灌浆期 59/92 个gene FDR < 0.05、
> OsBLE1 在 87 灌浆期 FDR = 6.7e-23"，**该results是错误编码所致**——当时用 `topTable(fit, coef = "groupMix")`，
> 而在 `~ 0 + group`（无截距 + 单因子）写法下系数是**组均值**而非 Mix − Pure 差异，
> test的其实是"该gene是否expression"。更正后：无gene通过 FDR correction；**OsBLE1 在 87 灌浆期确实不显著**
> （logFC = 0.10、P ≈ 0.49），manuscript"cross_stage一致上调"的表述**不成立**。

• 需要写明的两点局限：① 该表**没有芯片/批次变量**，虽然 4 个 focal 一致性能降低批次混杂的可能，
但不能完全排除；② limma 的 logFC 是 log2(x+1) 尺度、经 focal correction后的effect，and主题 04 中"均值比 log2FC"
不是同一个量（例如 OsGNS1 分蘖期 limma logFC = 5.99 vs 均值比 log2FC = 1.66），正文引用时必须注明尺度。

### Q9 field处理的跨年stability（Figure S14）

• 方法：疾病用 `glmer(cbind(diseased, 10 − diseased) ~ condition * Year + (1|Field), binomial)`，
性状用 `lmer(Value ~ condition * Year + (1|Field))`；另给逐年effectand留一年法（leave-one-year-out）。
只用 4 个焦点纯系。

• **condition × Year 交互在全部 5 个类别上均显著**（Blast P = 3.9e-28、BrownSpot P = 0.0195、Smut P = 1.6e-37、
PH P = 1.3e-180、KGW 亦显著，FDR 全部 < 0.05）——处理effect**随year变化**。

• 逐年effect（疾病为 OR，性状为均值差）：
 - Blast：2019 OR = 0.70（P = 0.058）、2020 OR = 1.41（P = 0.55）→ **方向反转**；
 - BrownSpot：2019 OR = 0.54（P = 5.4e-06）、2020 OR = 0.83（P = 0.16）→ 方向一致但 2020 减弱；
 - Smut：2019 OR = 1.31、2020 OR = 0.64 → 方向反转；
 - PH：2019 +3.20 cm（P = 0.0015）、**2020 −1.36 cm（P = 0.097）**、2021 +1.83 cm（P = 0.0065）→ 2020 反转；
 - KGW：2019 +0.57 g（P = 0.0065）、2020 +1.11 g（P = 2.3e-11）、2021 −0.25 g（P = 0.55）→ 2/3 年为正。

• 需注意的data事实：**病害（Blast / BrownSpot / Smut）只有 2019 and 2020 两年record**（每类 246 行），
PH / KGW 有 3 年；因此疾病侧的"stability"只能在两年间评估。

• 结论：**混种effect并不跨年稳定**，只有胡麻叶斑病在两年间方向一致（且 2020 明显减弱）。
正文若声称"全生育期持续、可泛化"的互作effect，必须补上这一year异质性notes，
否则结论会被year驱动（2020 年 PH 甚至反号）。

### 追加analysisresults汇总表

| analysis | 方法 | results | 对正文的影响 |
|---|---|---|---|
| Q6 transcriptomepartner_effect | 两因素 ANOVA + 99 次置换校准 | 交互项 η² and纯加性零分布无差别（P = 0.42 / 0.97） | 不可声称expression层面存在组合特异effect |
| Q7 plant_heightshape | change_points网格 + AIC + 199 次自助 | 线性差模型权重 100%，tau = 36 天但阶跃不显著 | 表述为"速率差异"而非"某时点后的跃升" |
| Q8 limma（已更正） | focal 区组 + 经验贝叶斯 + **显式 Mix − Pure 对比** | 两个时期均无gene FDR < 0.05；分蘖期 14 个名义显著 | 单gene结论只能写"名义趋势"；**OsBLE1 两个时期均未通过correction** |
| Q9 跨年stability | condition × Year + 逐年 + 留一年法 | 5 个类别交互项全部显著；Blast / Smut / PH 出现year反转 | 必须写明year异质性，避免"持续、可泛化"表述 |

### 输出

| file | 内容 |
|---|---|
| `04_表_partner_effect_variance_decomposition.xlsx`、`04_FigureS11_partner_effect_variance_decomposition.pdf` / `.png` | Q6：逐gene η²、置换参照、pathway层面（13 × 10 in） |
| `05_表_plant_height_effect_shape_and_change_points.xlsx`、`05_FigureS12_plant_heighteffect_shape_and_change_points.pdf` / `.png` | Q7：change_points网格、shapecomparison、自助results（15 × 5 in） |
| `06_表_expression_empirical_Bayes_differential_test.xlsx`、`06_FigureS13_经验贝叶斯火山图.pdf` / `.png` | Q8：limma results、and Wilcoxon 对比、focal 一致性（11 × 6 in） |
| `07_表_field_year_and_field_stability.xlsx`、`07_FigureS14_fieldyearstability.pdf` / `.png` | Q9：主模型、逐年effect、方向一致性、留一年法（12 × 6 in） |

---

## 更正notes（2026-09-22，Q1 的表述）

Q1 最初写作"叶片级伪重复（同一盆 6–14 片叶当独立重复）"。用户随后指出采样方式为
**"每一株就是一个个体，一个盆里有多株，采样时一株只取一片叶"**；据此复核data得到一致结论：
一个 image（盆）内有 6–14 株、`leaf_number` 在 image 内唯一、一个 image 只对应一个 focal × neighbor 组合。

因此 Q1 的**statistics做法and数字不变**（加盆随机截距后胡麻叶斑病 P 由 0.00021 变为 0.018、按盆聚合复核 0.021；
盆间 ICC = 0.12），但**性质需要更正**：

* 不是"同一株被重复采样的伪重复"，而是"同一盆内多株植物共享环境、彼此不独立"的**clustering**；
* 有效独立单位是盆（96 盆）而不是株（1043 株），这才是旧版 df ≈ 1039 偏高的原因；
* 正文描述实验时应写明"每个组合 4 盆、每盆多株、每株取一片叶"，并把重复数按盆计。

按"改写旧规则时追加新条目、不删除旧条目"的原则，原表述保留在 CHANGELOG 的历史条目中，此处仅作更正notes。
