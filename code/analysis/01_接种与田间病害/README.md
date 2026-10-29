# 01_接种与田间病害 — 运行记录

温室接种实验（稻瘟病 blast / 胡麻叶斑病 brown spot）与元阳梯田田间病害调查的重分析。
对应 Figure S2A–S2H 与田间病害调查小节（正文引用）。

## 运行环境

- R：`/opt/homebrew/bin/R`（`R version 4.6.1 (2026-06-24)`，Platform `aarch64-apple-darwin27.0.0`）
- 运行方式：`/opt/homebrew/bin/Rscript "<脚本绝对路径>"`；运行位置：本地 MacBook（本项目无超算端）
- 关键包版本：tidyverse 2.0.0、ggplot2 4.0.3、lme4 2.0.6、lmerTest 3.2.1、emmeans 2.0.4、
  multcomp 1.4.32、bioRtools 1.95.1、readxl 1.5.0.1、writexl 2.0.1
- 固定随机种子：各脚本 `set.seed(12345)`（本主题模型均为确定性拟合，种子用于将来加入重抽样时保持一致）

## 脚本与运行命令（2026-09-22）

```bash
cd "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
/opt/homebrew/bin/Rscript "code/01_接种与田间病害/01_接种数据统计.R"
/opt/homebrew/bin/Rscript "code/01_接种与田间病害/02_接种_各focal组合统计.R"
/opt/homebrew/bin/Rscript "code/01_接种与田间病害/03_接种绘图_FigureS2A_F.R"
/opt/homebrew/bin/Rscript "code/01_接种与田间病害/04_田间病害_GLMM与字母.R"
/opt/homebrew/bin/Rscript "code/01_接种与田间病害/05_田间病害绘图.R"
```

运行日志：`result/01_接种与田间病害/00_运行日志_<脚本名>.log`（stdout）、`00_运行错误_<脚本名>.log`（stderr）；
5 个脚本退出码均为 0。

## 输入

| 输入 | 说明 |
|---|---|
| `data/01_分析输入表/01_接种结果_合并表.xlsx` | 1659 行 × 13 列；4 focal × 4 neighbor × 2 病害；blast 2 批次、brown spot 3 批次；PHE 160 张图片、每组合 4 张图 |
| `data/01_分析输入表/02_元阳梯田田间病害调查表.xlsx` | 5658 行 × 10 列；2019–2021 × 3 田块 × 6 纯系（15/36/49/69/74/87）；Blast/BrownSpot/Smut 各 318 行（每小区 10 株口径），另有 KGW 1140 行、PH 3564 行（本主题不用） |

## 方法与关键参数

1. **温室接种总体比较**（`01`）：`lmerTest::lmer(lesion_number_per_cm2 ~ group + (1 | batch))`，分病害拟合；
   组间检验取固定效应 contrast（`groupMix`），df 用 Satterthwaite；Z-score 按 `disease × batch` 标准化（与旧版一致）。
2. **各 focal 邻体组合比较**（`02`）：`lmer(lesion_number_per_cm2 ~ combination + (1 | batch))`，逐 病害 × focal；
   `emmeans` + `multcomp::cld(adjust = "sidak")` 出字母，另出「与单作对照」的 `trt.vs.ctrl` 比较表。
3. **田间病害**（`04`）：`glmer(cbind(diseased, 10 − diseased) ~ combination + (1 | Year), family = binomial)`，
   逐 纯系 × 病害；`emmeans(type = "response")` 给发病率 + 95% 渐近 CI，字母在 link 尺度用 `cld(adjust = "sidak")` 计算后对应到 Prob 尺度
   （Prob 与 logit 单调对应，可按 combination 直接合并）。
4. **出图**（`03`、`05`）：`bioRtools::theme_prism()` + `bioRtools::scale_fill_research()` / `scale_color_research()`，
   `bioRtools::write_data(..., units = "in", dpi = 500)`；同时输出 PDF（矢量）与 PNG（500 dpi）。

### 统计口径须记住的三点

- **组均值（emmean）的 df 极小**：batch 只有 2–3 个水平且批次方差远大于残差（brown spot 的 ICC = 0.84），
  故 emmean 的 df ≈ 1–2、CI 极宽，**不能用于组间比较**；比较一律看 contrast（df = 613 / 1039）。
  因此 S2C–S2F 沿用旧版"箱线图 + 字母"的形式，不画误差棒；S2G/S2H 用的是 GLMM（渐近 z 推断），CI 正常可用。
- **batch 与 combination 不混杂**：每个 combination 出现在 2–3 个批次中（已实测核对），不存在"组合只在一个批次里"的设计缺陷。
- **多重比较校正**：emmeans 对单一 factor 的全部两两比较会把 `adjust = "tukey"` 自动改成 `"sidak"`
  （旧版脚本同样如此，故旧版字母实际是 Sidak 校正）。本主题显式写 `"sidak"`；已实测 **Tukey 与 Sidak 给出完全相同的字母**，
  故该选择不影响任何结论。

## 输出

| 文件 | 内容 |
|---|---|
| `01_表_接种LMM总体检验.xlsx` | 5 个 sheet：fixed_effect、anova_type3、model_info（含 ICC）、emmean_group_mean、descriptive |
| `02_表_接种Zscore中间表.xlsx` | 1659 行 + `lesion_norm`（按 disease × batch Z-score），S2A–S2F 的绘图输入 |
| `03_表_接种各focal组合emmeans字母.xlsx` | 3 个 sheet：emmean_letters、`pure_vs_mix_contrasts`（与单作对照比较，Sidak）、pairwise_contrasts |
| `04_表_田间病害主要GLMM.xlsx` | 3 个 sheet：emmean_letters_main（4 focal × Blast/BrownSpot）、overall_condition（OR + 95% CI）、descriptive |
| `05_表_田间病害补充_Smut与36_74.xlsx` | emmean_letters_all（6 纯系 × Blast/BrownSpot/Smut，18 个 GLMM） |
| `07–12_FigureS2A–S2F_*.pdf` / `.png` | 接种数据：标准化前、标准化后、各 focal 组合比较（8 × 6 in） |
| `13_FigureS2G_田间病害调查.pdf` / `.png` | 田间病害主图：4 焦点纯系 × Blast/BrownSpot（12 × 6 in） |
| `14_FigureS2H_田间病害调查_全部纯系与Smut.pdf` / `.png` | 补充图：6 纯系 × 3 病害（16 × 9 in） |

## 结果要点（2026-09-22）

- **温室接种总体**：brown spot 混种显著降低病斑密度（Estimate = −1.229，SE = 0.330，P = 0.000208，n = 1043）；
  blast 无显著差异（Estimate = −0.168，SE = 0.244，P = 0.491，n = 616）。
- **田间病害总体（4 焦点纯系）**：brown spot 混种发病率显著低于单作（OR = 0.678，95% CI 0.558–0.823，P = 8.8e-05）；
  blast 方向一致但未达显著（OR = 0.720，95% CI 0.509–1.018，P = 0.063）。与温室结果方向一致。
- 田间模型有 1 个边界奇异拟合（BrownSpot 焦点 15，Year 方差被估为 0），已在输出 `singular` 列标记，不影响该组组合比较结论。

## 与旧版 / 手稿的交叉核对

| 项目 | 旧版（手稿 v1 正文） | 本次重分析 | 结论 |
|---|---|---|---|
| brown spot 总体 LMM | Estimate = 1.229，SE = 0.33，P = 0.0002 | Estimate = −1.2294，SE = 0.3303，P = 0.000208 | **完全一致**（符号差异仅为 Mix 参照方向） |
| blast 总体 LMM | P = 0.491 | P = 0.4911 | **完全一致** |
| 模型设定、字母调整方法 | `lmer(~ group + (1\|batch))`、`cld(adjust = "tukey")` | 同一模型；`tukey` 与 `sidak` 实测字母相同 | 一致 |

**须提请作者注意的三处正文表述**（按 `pure_vs_mix_contrasts`（Sidak 校正）逐条核对）：

1. "混种 49 提高 15 对胡麻叶斑病的抗性" —— 邻体 49：−1.97，P = 0.048（**成立但为边缘显著**）；
   若按全部两两比较（6 个对比）的字母体系则与单作共享字母、不显著。建议正文改为"边缘显著"或补注校正口径。
2. "与 15 混种降低 49 对稻瘟病的抗性" —— 邻体 15：+1.26，P = 0.056，**未达显著**，旧版正文表述偏强。
3. "15 和 87 分别提高 69 对稻瘟病和胡麻叶斑病的抗性" —— 稻瘟病（邻体 15）−1.70，P = 0.0072 **成立**；
   胡麻叶斑病（邻体 87）−1.83，P = 0.073 **未达显著**。
   另：胡麻叶斑病下 49 的各邻体均不显著（P = 0.36 / 0.96 / 0.85），与"其余三个纯系均无显著影响"的表述基本相符。

## 备注

- 旧版"总体检验"用 `lme4::lmer() + summary()`，本身不输出 P 值；手稿中的 P 值与本次 Satterthwaite 结果一致，说明当时另有一次带 lmerTest 的运行。
- `lesion_numner` 为历史误拼字段（原 `lesion_nb`），按既定约定保留。
- 图内 P 值统一用定点写法（`P < 0.001` / `P = 0.491`）。R 默认会把小于约 0.001 的数值输出成科学计数法（如 `2.1e-04`），
  Excel 窄列也会显示成 `2.75E-05`，易被误读为显著性字母，故统一规避。
