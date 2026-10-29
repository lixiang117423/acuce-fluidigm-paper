# ============================================================
# 02_prepare_field_disease_survey.R
# 作用：把元阳梯田field_disease调查原始表转成analysis输入表（内容and上游一致，仅统一存放位置and读取入口）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/00_data_processing/02_prepare_field_disease_survey.R"
#
# 输入（本项目内只读副本）：
#   <root>/data/02_raw_data/02_field_disease调查/All_Data_in_YYT.xlsx
# 输入上游原始路径（副本来源，未修改）：
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/field试验/All_Data_in_YYT.xlsx
# 输出：<root>/data/processed_inputs/02_field_disease_survey.xlsx
# 其中 <root> = <repository_root>
#
# 字段notes：row num / FocalPlant / Neighbor / condition / Field / Year / Month / Day / Category / Value；
#          Category 含 Blast、BrownSpot、Smut 等病害类型，Value = 调查值。
# 处理约定：本表and上游 All_Data_in_YYT.xlsx 内容一致（已核对 5658 行 x 10 列逐格相同），不做筛选或改写；
#          analysis脚本读入后自行 dplyr::select(-1) 去掉 row num。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：本步耗时 < 5 s 且输出确定，每次运行都重新generate，便于复现validation。
# versionrecord：脚本末尾打印 sessionInfo()；正式运行record写在 code/00_data_processing/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

project_root <- normalizePath(getwd())

path_in <- file.path(project_root, "data", "02_raw_data", "02_field_disease调查", "All_Data_in_YYT.xlsx")
path_out <- file.path(project_root, "data", "processed_inputs", "02_field_disease_survey.xlsx")

stopifnot(file.exists(path_in))

cat("[02] 输入：", path_in, "\n")

readxl::read_excel(path_in) -> df

# complete性校验：行数and列名必须and上游一致，防止读到截断file
stopifnot(
  nrow(df) == 5658,
  identical(
    colnames(df),
    c("row num", "FocalPlant", "Neighbor", "condition", "Field", "Year", "Month", "Day", "Category", "Value")
  )
)

cat("[02] 读入：", nrow(df), "行 x", ncol(df), "列\n")
cat("[02] Category 计数：\n")
print(table(df$Category, useNA = "ifany"))
cat("[02] FocalPlant 计数：\n")
print(table(df$FocalPlant, useNA = "ifany"))
cat("[02] condition 计数：\n")
print(table(df$condition, useNA = "ifany"))

writexl::write_xlsx(df, path_out)
cat("[02] 输出：", path_out, "\n")

cat("\n[02] sessionInfo()\n")
sessionInfo()
