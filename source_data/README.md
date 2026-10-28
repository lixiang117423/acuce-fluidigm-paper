# Source Data 清单（manuscript/source_data/）

> 本目录**只放数据与清单**：`Source_Data.xlsx`（总工作簿，一个子图一张 sheet）＋ 本 README。
> **绘图代码统一放在 `manuscript/code/`**（一图一脚本：`Fig1.R`…`Fig5.R`、`FigS1.R`、`FigS2.R`…`FigS19.R`），
> 代码只读本工作簿、不再读原始大表；`manuscript/code/build_source_data.R` 是**生成这本工作簿的组装脚本**
> （读 `result/` 与 `data/` 下既有产物，不重跑任何分析）。

## 一、sheet 命名与布局

- 命名：多面板图 = 图号 + 面板小写字母（`Fig1a`、`Fig1b`…、`FigS2a`…）；**单面板图 = 图号**（`FigS3`、`FigS4`…，不带字母）。
- 顺序：图序 → 面板序（Fig1a、Fig1b、Fig2a … Fig5c、FigS2a … FigS19b）。
- 每个 sheet 的布局：**第 1 行 = 元信息**（`figure_panel=` / `code_file=` / `generated=` / `software=` / `note=`）；
  第 2 行空行；**第 3 行 = 表头**；第 4 行起 = 数据。代码读取方式：`readxl::read_excel("source_data/Source_Data.xlsx", sheet = "Fig1a", skip = 2)`。
- 多面板图的 sheet 内均含 `Panel` 列（值与该 sheet 名后缀一致）；单面板图无后缀、也无 `Panel` 列。

## 二、子图 ↔ sheet ↔ 代码对照

| 图 | 子图 | sheet | 代码 | 说明 |
|---|---|---|---|---|
| Fig. 1 | a–e | `Fig1a`–`Fig1e` | `code/Fig1.R` | 农艺性状（株高变化指数、分蘖数、产量、表型 PCA、三性状整体比较） |
| Fig. 2 | a–e | `Fig2a`–`Fig2e` | `code/Fig2.R` | 通路水平热图、PLS-DA 得分图、Dim1/Dim2 分布、VIP 载荷 |
| Fig. 3 | a、b | `Fig3a`、`Fig3b` | `code/Fig3.R` | 防御基因热图、每样本防御分数 |
| Fig. 4 | a、b | `Fig4a`、`Fig4b` | `code/Fig4.R` | 跨阶段 limma（防御基因、高糖通路基因） |
| Fig. 5 | a–c | `Fig5a`–`Fig5c` | `code/Fig5.R` | 结构分解、reciprocity、逐年田间效应 |
| Fig. S1 | — | — | `code/FigS1.R` | **示意图，无 Source Data**（实验设计与流程；按规则登记例外） |
| Fig. S2 | a–h | `FigS2a`–`FigS2h` | `code/FigS2.R` | 温室接种（原始/Z 标准化、各 focal 邻体比较）与田间病害 |
| Fig. S3 | 单面板 | `FigS3` | `code/FigS3.R` | VIP 前 10 基因 × 遗传距离（三时期 + 相关统计） |
| Fig. S4 | 单面板 | `FigS4` | `code/FigS4.R` | 田间株高与千粒重（emmeans ± 95% CI） |
| Fig. S5 | 单面板 | `FigS5` | `code/FigS5.R` | 高糖–防御基因 Spearman 相关热图 |
| Fig. S6 | a、b | `FigS6a`、`FigS6b` | `code/FigS6.R` | 激素类别与物质检出率 |
| Fig. S7 | 单面板 | `FigS7` | `code/FigS7.R` | 全检出物质的每样本值（含描述统计与 KW 检验） |
| Fig. S8 | a、b | `FigS8a`、`FigS8b` | `code/FigS8.R` | 分析单元比较、reciprocity |
| Fig. S9 | a、b | `FigS9a`、`FigS9b` | `code/FigS9.R` | 株高生长曲线、逐日 P 值曲线 |
| Fig. S10 | 单面板 | `FigS10` | `code/FigS10.R` | 通路分数森林图（三时期） |
| Fig. S11 | a–c | `FigS11a`–`FigS11c` | `code/FigS11.R` | 伙伴效应方差分解（η² 分布、散点、基因排序） |
| Fig. S12 | a–c | `FigS12a`–`FigS12c` | `code/FigS12.R` | 株高效应形状与变化点（差值曲线、剖面似然、自助分布） |
| Fig. S13 | 单面板 | `FigS13` | `code/FigS13.R` | limma 火山图（三时期） |
| Fig. S14 | 单面板 | `FigS14` | `code/FigS14.R` | 田间效应逐年稳定性 |
| Fig. S15 | a、b | `FigS15a`、`FigS15b` | `code/FigS15.R` | 阶段特异性（平行坐标、各对比计数） |
| Fig. S16 | 单面板 | `FigS16` | `code/FigS16.R` | 遗传距离 × 农艺性状（12 个混种组合） |
| Fig. S17 | 单面板 | `FigS17` | `code/FigS17.R` | 独脚金内酯通路基因（三时期） |
| Fig. S18 | 单面板 | `FigS18` | `code/FigS18.R` | 基因级 92 基因热图（整页） |
| Fig. S19 | a、b | `FigS19a`、`FigS19b` | `code/FigS19.R` | 焦点纯系 69 分蘖期 / 87 灌浆期通路基因 |

## 三、怎么跑

```bash
cd "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript"
Rscript code/Fig1.R      # 逐图运行；每个脚本内部按子图读各自的 sheet
```

- 代码里的路径一律是**相对 `manuscript/` 的相对路径**（`source_data/Source_Data.xlsx`、`figures/...`），换机器可用；
- 每个 `Fig*.R` 会把面板写入 `figures/panels/`，再用 `code/_vectorcraft.R`（VectorCraft CLI 薄封装）拼成整图 → `figures/FigN.pdf`（矢量）＋ `FigN.png`（600 dpi，由同一 PDF 转）；
- 组图需要 VectorCraft CLI 0.12.0；路径可用环境变量 `VECTORCRAFT_CLI` 覆盖（默认 `/Users/lixiang/tools/vectorcraft/cli-0.12.0/vectorcraft-cli-0.12.0-macos-universal/vectorcraft-cli`）。

## 四、自检（子图清单 ↔ sheet ↔ 代码 三者对账）

```bash
cd "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript"
Rscript code/check_source_data.R
```

检查项：① sheet 名单与上表子图清单**一一对应**（扣掉登记的例外 Fig. S1）；② 多面板 sheet 不得出现裸 `FigN`、单面板 sheet 不得带字母，且 `Panel` 列与 sheet 名后缀一致；③ 每个 sheet 至少一行数据、无空 sheet；④ 每个子图都有对应的 `FigN.R`。
