# ============================================================
# 03_prepare_agronomic_traits.R
# 作用：把tiller_number、单株yield、plant_height三张原始调查表合并成一张长表analysis输入表
# 运行：/opt/homebrew/bin/Rscript "<root>/code/00_data_processing/03_prepare_agronomic_traits.R"
#
# 输入（本项目内只读副本）：
#   <root>/data/02_raw_data/03_agronomic_traits/tiller_number.xlsx
#   <root>/data/02_raw_data/03_agronomic_traits/单株yield.xlsx
#   <root>/data/02_raw_data/03_agronomic_traits/plant_height调查results.xlsx
# 输入上游原始路径（副本来源，未修改）：
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/02.fluidigm文章/data/tiller_number.xlsx
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/02.fluidigm文章/data/单株yield.xlsx
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/02.fluidigm文章/data/plant_height调查results.xlsx
# 输出：<root>/data/processed_inputs/03_agronomic_traits_height_tiller_yield.xlsx
# 其中 <root> = <repository_root>
#
# 字段notes：focal / neighbor / group（Pure 或 Mix）/ time（YYYYMMDD）/ value / phenotype
#          （tiller number、yield、plant height）。
# 处理约定（沿用旧版）：
#   1) 单株yield表无 time 列，统一记调查日期 20221227，并丢弃辅助列 temp；
#   2) plant_height表中 Pure 组的 neighbor 统一改写为 focal 自身，代表同纯系对照；
#   3) 三张表统一为tiller_number表的列顺序。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算；time 为无分隔日期，勿当数值参and运算。
# 断点续传：本步耗时 < 5 s 且输出确定，每次运行都重新generate，便于复现validation。
# versionrecord：脚本末尾打印 sessionInfo()；正式运行record写在 code/00_data_processing/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

project_root <- normalizePath(getwd())

path_in_tiller <- file.path(project_root, "data", "02_raw_data", "03_agronomic_traits", "tiller_number.xlsx")
path_in_yield <- file.path(project_root, "data", "02_raw_data", "03_agronomic_traits", "单株yield.xlsx")
path_in_height <- file.path(project_root, "data", "02_raw_data", "03_agronomic_traits", "plant_height调查results.xlsx")
path_out <- file.path(project_root, "data", "processed_inputs", "03_agronomic_traits_height_tiller_yield.xlsx")

stopifnot(file.exists(path_in_tiller), file.exists(path_in_yield), file.exists(path_in_height))

cat("[03] 输入：", path_in_tiller, "\n")
cat("[03] 输入：", path_in_yield, "\n")
cat("[03] 输入：", path_in_height, "\n")

# tiller_number：作为列顺序的基准表
readxl::read_excel(path_in_tiller) %>%
  dplyr::mutate(phenotype = "tiller number") -> df.tiller

# 单株yield：无 time 列，补调查日期 20221227
readxl::read_excel(path_in_yield) %>%
  dplyr::mutate(phenotype = "yield") %>%
  dplyr::select(-temp) %>%
  dplyr::mutate(time = 20221227) %>%
  dplyr::select(dplyr::all_of(colnames(df.tiller))) -> df.yield

# plant_height：Pure 组用 focal 自身作为 neighbor
readxl::read_excel(path_in_height) %>%
  dplyr::mutate(phenotype = "plant height") %>%
  dplyr::mutate(neighbor = dplyr::case_when(group == "Pure" ~ focal, TRUE ~ neighbor)) %>%
  dplyr::select(dplyr::all_of(colnames(df.tiller))) -> df.height

cat("[03] tiller_number：", nrow(df.tiller), "行；单株yield：", nrow(df.yield), "行；plant_height：", nrow(df.height), "行\n")

dplyr::bind_rows(df.tiller, df.yield, df.height) -> df.out

cat("[03] 合并后：", nrow(df.out), "行 x", ncol(df.out), "列\n")
cat("[03] 列名：", paste(colnames(df.out), collapse = ", "), "\n")
cat("[03] phenotype x group 计数：\n")
print(table(df.out$phenotype, df.out$group))

writexl::write_xlsx(df.out, path_out)
cat("[03] 输出：", path_out, "\n")

cat("\n[03] sessionInfo()\n")
sessionInfo()
