# ============================================================
# 01_prepare_inoculation_results.R
# 作用：把温室inoculation（稻瘟病 blast + 胡麻叶斑病 brown spot）原始results合并成analysis输入表
# 运行：/opt/homebrew/bin/Rscript "<root>/code/00_data_processing/01_prepare_inoculation_results.R"
#
# 输入（本项目内只读副本）：
#   <root>/data/02_raw_data/01_温室inoculation/blast.inoculation.result.xlsx
#   <root>/data/02_raw_data/01_温室inoculation/brown.spot.inoculation.result.xlsx
# 输入上游原始路径（副本来源，未修改）：
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/02.fluidigm文章/data/blast.inoculation.result.xlsx
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/02.fluidigm文章/data/brown.spot.inoculation.result.xlsx
# 输出：<root>/data/processed_inputs/01_inoculation_merged.xlsx
# 其中 <root> = <repository_root>
#
# 字段notes：lesion_numner 沿用旧版拼写（原字段为 lesion_nb，旧脚本误拼为 numner，为and历史results一致未更正）；
#          lesion_number_per_cm2 = 单位面积病斑数（每 cm2）；lesion_percent = 病斑占比（%）。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：本步耗时 < 5 s 且输出确定，不做"输出已存在则跳过"判断，每次运行都重新generate，便于复现validation。
# versionrecord：脚本末尾打印 sessionInfo()；正式运行record写在 code/00_data_processing/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

project_root <- normalizePath(getwd())

path_in_blast <- file.path(project_root, "data", "02_raw_data", "01_温室inoculation", "blast.inoculation.result.xlsx")
path_in_spot <- file.path(project_root, "data", "02_raw_data", "01_温室inoculation", "brown.spot.inoculation.result.xlsx")
path_out <- file.path(project_root, "data", "processed_inputs", "01_inoculation_merged.xlsx")

stopifnot(file.exists(path_in_blast), file.exists(path_in_spot))

cat("[01] 输入：", path_in_blast, "\n")
cat("[01] 输入：", path_in_spot, "\n")

# 稻瘟病（blast）
readxl::read_excel(path_in_blast) %>%
  dplyr::mutate(disease = "blast") %>%
  dplyr::rename(batch = group) %>%
  dplyr::mutate(batch = as.character(batch)) -> df.blast

# 胡麻叶斑病（brown spot）
readxl::read_excel(path_in_spot) %>%
  dplyr::mutate(disease = "brown spot") %>%
  dplyr::rename(batch = group) %>%
  dplyr::mutate(batch = as.character(batch)) -> df.spot

cat("[01] blast 原始行数：", nrow(df.blast), "；brown spot 原始行数：", nrow(df.spot), "\n")

# 合并两种病害，并统一字段顺序and命名
dplyr::bind_rows(df.blast, df.spot) %>%
  dplyr::rename(group = group.2) %>%
  dplyr::select(1, 2, 15, 10, 11, 3, 13, 4:9) %>%
  dplyr::rename_with(~ gsub("\\.", "_", .x)) %>%
  dplyr::rename(lesion_numner = lesion_nb) %>%
  dplyr::rename(lesion_number_per_cm2 = lesion_per_area) %>%
  dplyr::rename(lesion_percent = pourcent_lesions) -> df.out

cat("[01] 合并后：", nrow(df.out), "行 x", ncol(df.out), "列\n")
cat("[01] 列名：", paste(colnames(df.out), collapse = ", "), "\n")
cat("[01] disease x group 计数：\n")
print(table(df.out$disease, df.out$group))

writexl::write_xlsx(df.out, path_out)
cat("[01] 输出：", path_out, "\n")

cat("\n[01] sessionInfo()\n")
sessionInfo()
