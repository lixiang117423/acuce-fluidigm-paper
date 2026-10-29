# 00_数据整理 — 运行记录

把 `../../data/02_原始数据/` 的上游原始表整理成 `../../data/01_分析输入表/` 的分析入口表。

## 运行环境

- R：`/opt/homebrew/bin/R`（`R version 4.6.1 (2026-06-24)`，Platform `aarch64-apple-darwin27.0.0`）
- 运行方式：`/opt/homebrew/bin/Rscript "<脚本绝对路径>"`
- 运行位置：本地 MacBook（本项目无超算端）
- 关键包版本：tidyverse 2.0.0（dplyr 1.2.1、tidyr 1.3.2、stringr 1.6.0）、readxl 1.5.0.1、writexl 2.0.1、`fluidgr` 0.1（`/Users/lixiang/R/library`，2026-09-22 由用户安装）

## 2026-09-22 09:07 首次运行记录

> 2026-09-22 09:10 补充：用户完成 `fluidgr` 0.1 安装后，第 04 步已补齐并跑通（见下表最后一行），阶段 0 五个脚本全部完成。

运行命令（逐条执行，stdout 落 `../../result/00_数据整理/00_运行日志_<脚本名>.log`，stderr 落 `00_运行错误_<脚本名>.log`）：

```bash
/opt/homebrew/bin/Rscript "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/code/00_数据整理/01_整理接种结果表.R"
/opt/homebrew/bin/Rscript "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/code/00_数据整理/02_整理田间病害调查表.R"
/opt/homebrew/bin/Rscript "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/code/00_数据整理/03_整理农艺性状表.R"
/opt/homebrew/bin/Rscript "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/code/00_数据整理/05_整理激素_迈维表.R"
```

| 脚本 | 退出码 | 输入 | 输出 | 行 × 列 | 与旧版比对 |
|---|---|---|---|---|---|
| `01_整理接种结果表.R` | 0 | `02_原始数据/01_温室接种/blast.inoculation.result.xlsx`（616 行）、`brown.spot.inoculation.result.xlsx`（1043 行） | `01_分析输入表/01_接种结果_合并表.xlsx` | 1659 × 13 | 与旧表**逐格一致**（13 列 `all.equal` 全 TRUE） |
| `02_整理田间病害调查表.R` | 0 | `02_原始数据/02_田间病害调查/All_Data_in_YYT.xlsx` | `01_分析输入表/02_元阳梯田田间病害调查表.xlsx` | 5658 × 10 | 与旧表**逐格一致**（10 列 `all.equal` 全 TRUE） |
| `03_整理农艺性状表.R` | 0 | `02_原始数据/03_农艺性状/` 分蘖数（200 行）+ 单株产量（195 行）+ 株高调查结果（7700 行） | `01_分析输入表/03_分蘖数_株高_单株产量.xlsx` | 8095 × 6 | 与旧表**逐格一致**（6 列 `all.equal` 全 TRUE） |
| `04_整理fluidigm表达量表.R` | 0 | `02_原始数据/04_Fluidigm表达/`：引物表 + 第一批下机 csv（9120 行 / 96 样本）+ 第三批下机 csv（8835 行 / 93 样本） | `01_分析输入表/04_fluidigm标准化表达量.xlsx` | 15456 × 10 | 与旧表**逐格一致**（10 列不一致单元格数全为 0，`expression` 差值恒为 0） |
| `05_整理激素_迈维表.R` | 0 | `02_原始数据/06_激素_迈维/0.Data/data.levels.xlsx` + `sample_info.xlsx` | `01_分析输入表/05_激素_迈维含量表.xlsx` | 1308 × 12 | 新增脚本，无旧版可比 |

## 运行中发现的问题（须在后续分析中处理）

1. **1 / 2 / 3 号表内容与旧版完全一致，字节不同**：`.xlsx` 的 zip 时间戳与写入顺序导致 MD5 变化，内容经逐列 `all.equal` 核验一致。
2. **迈维激素数据不含油菜素甾醇**：109 个物质属生长素 / 细胞分裂素 / 赤霉素 / 茉莉酸 / 水杨酸 / 脱落酸 / 褪黑素 / 独角金内酯 / 乙烯类，**没有 BR**。
   论文中 BR 的证据来自 Fluidigm 基因表达（如 `OsBLE1`）。使用该数据时不得写成 BR 含量测定结果。
3. **迈维激素数据检出率低**：1308 个样本 × 物质组合中仅 254 个非缺失（约 19%），1054 个为 NA。
   全检出（12/12）的仅 12 个物质：`12-OH-JA`、`ACC`、`BAP`、`IAA`、`IAA-Phe`、`IAA-Val-Me`、`IAN`、`IP`、`TRP`、`GA4`（11/12）、`H2JA`（10/12）、`mT9G`（10/12）。
   分析口径（只用全检出物质 / 按类别汇总 / 调阈值）需与用户确认。

## 输出校验值（2026-09-22 09:07）

| 文件 | MD5 |
|---|---|
| `01_分析输入表/01_接种结果_合并表.xlsx` | `35cca7f9f22b6680a369ef6238799d4d` |
| `01_分析输入表/02_元阳梯田田间病害调查表.xlsx` | `5e7485fe75edc23c45ec2f17b0068c65` |
| `01_分析输入表/03_分蘖数_株高_单株产量.xlsx` | `94001bb8205e6c38bc56b138bb19ff97` |
| `01_分析输入表/04_fluidigm标准化表达量.xlsx` | `cb0abee6751d1b67ec2a79c41b494dd3`（2026-09-22 09:10 重跑后） |
| `01_分析输入表/05_激素_迈维含量表.xlsx` | `c1f1f609c87a1f69626208722ac04bd0` |
