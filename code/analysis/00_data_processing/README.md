# 00_data_processing — 运行record

把 `../../data/02_raw_data/` 的上游原始表prepare成 `../../data/processed_inputs/` 的analysis入口表。

## 运行环境

- R：`/opt/homebrew/bin/R`（`R version 4.6.1 (2026-06-24)`，Platform `aarch64-apple-darwin27.0.0`）
- 运行方式：`/opt/homebrew/bin/Rscript "<脚本绝对路径>"`
- 运行位置：本地 MacBook（本项目无超算端）
- 关键包version：tidyverse 2.0.0（dplyr 1.2.1、tidyr 1.3.2、stringr 1.6.0）、readxl 1.5.0.1、writexl 2.0.1、`fluidgr` 0.1（`/Users/lixiang/R/library`，2026-09-22 由用户安装）

## 2026-09-22 09:07 首次运行record

> 2026-09-22 09:10 补充：用户完成 `fluidgr` 0.1 安装后，第 04 步已complete并跑通（见下表最后一行），stage 0 五个脚本全部完成。

运行命令（item_by_item执行，stdout 落 `../../result/00_data_processing/00_运行日志_<脚本名>.log`，stderr 落 `00_运行错误_<脚本名>.log`）：

```bash
/opt/homebrew/bin/Rscript "<repository_root>/code/00_data_processing/01_prepare_inoculation_results.R"
/opt/homebrew/bin/Rscript "<repository_root>/code/00_data_processing/02_prepare_field_disease_survey.R"
/opt/homebrew/bin/Rscript "<repository_root>/code/00_data_processing/03_prepare_agronomic_traits.R"
/opt/homebrew/bin/Rscript "<repository_root>/code/00_data_processing/05_prepare_phytohormone_MetWare_table.R"
```

| 脚本 | 退出码 | 输入 | 输出 | 行 × 列 | and旧版比对 |
|---|---|---|---|---|---|
| `01_prepare_inoculation_results.R` | 0 | `02_raw_data/01_温室inoculation/blast.inoculation.result.xlsx`（616 行）、`brown.spot.inoculation.result.xlsx`（1043 行） | `processed_inputs/01_inoculation_merged.xlsx` | 1659 × 13 | and旧表**逐格一致**（13 列 `all.equal` 全 TRUE） |
| `02_prepare_field_disease_survey.R` | 0 | `02_raw_data/02_field_disease调查/All_Data_in_YYT.xlsx` | `processed_inputs/02_field_disease_survey.xlsx` | 5658 × 10 | and旧表**逐格一致**（10 列 `all.equal` 全 TRUE） |
| `03_prepare_agronomic_traits.R` | 0 | `02_raw_data/03_agronomic_traits/` tiller_number（200 行）+ 单株yield（195 行）+ plant_height调查results（7700 行） | `processed_inputs/03_agronomic_traits_height_tiller_yield.xlsx` | 8095 × 6 | and旧表**逐格一致**（6 列 `all.equal` 全 TRUE） |
| `04_prepare_fluidigm_expression.R` | 0 | `02_raw_data/04_Fluidigmexpression/`：引物表 + 第一批下机 csv（9120 行 / 96 sample）+ 第三批下机 csv（8835 行 / 93 sample） | `processed_inputs/04_fluidigm_normalized_expression.xlsx` | 15456 × 10 | and旧表**逐格一致**（10 列不一致单元格数全为 0，`expression` 差值恒为 0） |
| `05_prepare_phytohormone_MetWare_table.R` | 0 | `02_raw_data/06_激素_迈维/0.Data/data.levels.xlsx` + `sample_info.xlsx` | `processed_inputs/05_激素_迈维含量表.xlsx` | 1308 × 12 | 新增脚本，无旧版可比 |

## 运行中发现的问题（须在后续analysis中处理）

1. **1 / 2 / 3 号表内容and旧版完全一致，字节不同**：`.xlsx` 的 zip time戳and写入顺序导致 MD5 变化，内容经逐列 `all.equal` validation一致。
2. **迈维激素data不含油菜素甾醇**：109 个物质属生长素 / 细胞分裂素 / 赤霉素 / 茉莉酸 / 水杨酸 / 脱落酸 / 褪黑素 / 独角金内酯 / 乙烯类，**没有 BR**。
   paper中 BR 的证据来自 Fluidigm geneexpression（如 `OsBLE1`）。使用该data时不得写成 BR 含量测定results。
3. **迈维激素data检出率低**：1308 个sample × 物质组合中仅 254 个非缺失（约 19%），1054 个为 NA。
   全检出（12/12）的仅 12 个物质：`12-OH-JA`、`ACC`、`BAP`、`IAA`、`IAA-Phe`、`IAA-Val-Me`、`IAN`、`IP`、`TRP`、`GA4`（11/12）、`H2JA`（10/12）、`mT9G`（10/12）。
   analysis口径（只用全检出物质 / 按类别汇总 / 调阈值）需and用户确认。

## 输出校验值（2026-09-22 09:07）

| file | MD5 |
|---|---|
| `processed_inputs/01_inoculation_merged.xlsx` | `35cca7f9f22b6680a369ef6238799d4d` |
| `processed_inputs/02_field_disease_survey.xlsx` | `5e7485fe75edc23c45ec2f17b0068c65` |
| `processed_inputs/03_agronomic_traits_height_tiller_yield.xlsx` | `94001bb8205e6c38bc56b138bb19ff97` |
| `processed_inputs/04_fluidigm_normalized_expression.xlsx` | `cb0abee6751d1b67ec2a79c41b494dd3`（2026-09-22 09:10 重跑后） |
| `processed_inputs/05_激素_迈维含量表.xlsx` | `c1f1f609c87a1f69626208722ac04bd0` |
