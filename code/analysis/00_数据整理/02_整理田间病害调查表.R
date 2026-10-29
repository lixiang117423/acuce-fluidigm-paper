# ============================================================
# 02_整理田间病害调查表.R
# 作用：把元阳梯田田间病害调查原始表转成分析输入表（内容与上游一致，仅统一存放位置与读取入口）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/00_数据整理/02_整理田间病害调查表.R"
#
# 输入（本项目内只读副本）：
#   <root>/data/02_原始数据/02_田间病害调查/All_Data_in_YYT.xlsx
# 输入上游原始路径（副本来源，未修改）：
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/田间试验/All_Data_in_YYT.xlsx
# 输出：<root>/data/01_分析输入表/02_元阳梯田田间病害调查表.xlsx
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 字段说明：row num / FocalPlant / Neighbor / condition / Field / Year / Month / Day / Category / Value；
#          Category 含 Blast、BrownSpot、Smut 等病害类型，Value = 调查值。
# 处理约定：本表与上游 All_Data_in_YYT.xlsx 内容一致（已核对 5658 行 x 10 列逐格相同），不做筛选或改写；
#          分析脚本读入后自行 dplyr::select(-1) 去掉 row num。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：本步耗时 < 5 s 且输出确定，每次运行都重新生成，便于复现核验。
# 版本记录：脚本末尾打印 sessionInfo()；正式运行记录写在 code/00_数据整理/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"

path_in <- file.path(project_root, "data", "02_原始数据", "02_田间病害调查", "All_Data_in_YYT.xlsx")
path_out <- file.path(project_root, "data", "01_分析输入表", "02_元阳梯田田间病害调查表.xlsx")

stopifnot(file.exists(path_in))

cat("[02] 输入：", path_in, "\n")

readxl::read_excel(path_in) -> df

# 完整性校验：行数与列名必须与上游一致，防止读到截断文件
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
