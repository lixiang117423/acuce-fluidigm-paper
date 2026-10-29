# ============================================================
# 03_整理农艺性状表.R
# 作用：把分蘖数、单株产量、株高三张原始调查表合并成一张长表分析输入表
# 运行：/opt/homebrew/bin/Rscript "<root>/code/00_数据整理/03_整理农艺性状表.R"
#
# 输入（本项目内只读副本）：
#   <root>/data/02_原始数据/03_农艺性状/分蘖数.xlsx
#   <root>/data/02_原始数据/03_农艺性状/单株产量.xlsx
#   <root>/data/02_原始数据/03_农艺性状/株高调查结果.xlsx
# 输入上游原始路径（副本来源，未修改）：
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/02.fluidigm文章/data/分蘖数.xlsx
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/02.fluidigm文章/data/单株产量.xlsx
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/02.fluidigm文章/data/株高调查结果.xlsx
# 输出：<root>/data/01_分析输入表/03_分蘖数_株高_单株产量.xlsx
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 字段说明：focal / neighbor / group（Pure 或 Mix）/ time（YYYYMMDD）/ value / phenotype
#          （tiller number、yield、plant height）。
# 处理约定（沿用旧版）：
#   1) 单株产量表无 time 列，统一记调查日期 20221227，并丢弃辅助列 temp；
#   2) 株高表中 Pure 组的 neighbor 统一改写为 focal 自身，代表同纯系对照；
#   3) 三张表统一为分蘖数表的列顺序。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算；time 为无分隔日期，勿当数值参与运算。
# 断点续传：本步耗时 < 5 s 且输出确定，每次运行都重新生成，便于复现核验。
# 版本记录：脚本末尾打印 sessionInfo()；正式运行记录写在 code/00_数据整理/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"

path_in_tiller <- file.path(project_root, "data", "02_原始数据", "03_农艺性状", "分蘖数.xlsx")
path_in_yield <- file.path(project_root, "data", "02_原始数据", "03_农艺性状", "单株产量.xlsx")
path_in_height <- file.path(project_root, "data", "02_原始数据", "03_农艺性状", "株高调查结果.xlsx")
path_out <- file.path(project_root, "data", "01_分析输入表", "03_分蘖数_株高_单株产量.xlsx")

stopifnot(file.exists(path_in_tiller), file.exists(path_in_yield), file.exists(path_in_height))

cat("[03] 输入：", path_in_tiller, "\n")
cat("[03] 输入：", path_in_yield, "\n")
cat("[03] 输入：", path_in_height, "\n")

# 分蘖数：作为列顺序的基准表
readxl::read_excel(path_in_tiller) %>%
  dplyr::mutate(phenotype = "tiller number") -> df.tiller

# 单株产量：无 time 列，补调查日期 20221227
readxl::read_excel(path_in_yield) %>%
  dplyr::mutate(phenotype = "yield") %>%
  dplyr::select(-temp) %>%
  dplyr::mutate(time = 20221227) %>%
  dplyr::select(dplyr::all_of(colnames(df.tiller))) -> df.yield

# 株高：Pure 组用 focal 自身作为 neighbor
readxl::read_excel(path_in_height) %>%
  dplyr::mutate(phenotype = "plant height") %>%
  dplyr::mutate(neighbor = dplyr::case_when(group == "Pure" ~ focal, TRUE ~ neighbor)) %>%
  dplyr::select(dplyr::all_of(colnames(df.tiller))) -> df.height

cat("[03] 分蘖数：", nrow(df.tiller), "行；单株产量：", nrow(df.yield), "行；株高：", nrow(df.height), "行\n")

dplyr::bind_rows(df.tiller, df.yield, df.height) -> df.out

cat("[03] 合并后：", nrow(df.out), "行 x", ncol(df.out), "列\n")
cat("[03] 列名：", paste(colnames(df.out), collapse = ", "), "\n")
cat("[03] phenotype x group 计数：\n")
print(table(df.out$phenotype, df.out$group))

writexl::write_xlsx(df.out, path_out)
cat("[03] 输出：", path_out, "\n")

cat("\n[03] sessionInfo()\n")
sessionInfo()
