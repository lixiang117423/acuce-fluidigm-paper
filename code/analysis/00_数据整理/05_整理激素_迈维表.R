# ============================================================
# 05_整理激素_迈维表.R
# 作用：把迈维（Metware）水稻根系分泌物激素靶向检测结果转成长表分析输入表
# 运行：/opt/homebrew/bin/Rscript "<root>/code/00_数据整理/05_整理激素_迈维表.R"
#
# 输入（本项目内只读副本）：
#   <root>/data/02_原始数据/06_激素_迈维/0.Data/data.levels.xlsx  含量表：109 物质 x (12 样本 + 3 QC)
#   <root>/data/02_原始数据/06_激素_迈维/0.Data/sample_info.xlsx  样本分组：A=15、B=69、C=1569
#   <root>/data/02_原始数据/06_激素_迈维/0.Data/unit.xlsx         单位：hormones = ng/mL
# 输入上游原始路径（副本来源，未修改）：
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/02.fluidigm文章/data/hormone-迈维/0.Data/
# 输出：<root>/data/01_分析输入表/05_激素_迈维含量表.xlsx
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 字段说明：sample（样本ID，如 15-1）/ group（A/B/C）/ focal（15 / 69 / 1569）/ replicate（1-4）/
#          index（物质简称）/ compound_en / compound_cn / class_en / class_cn / cpd_id / kegg_map /
#          value（含量，单位 ng/mL）
# 设计说明：12 个根系分泌物样本 = 纯系 15（A）x4、纯系 69（B）x4、混种 15+69（C）x4；
#          QC 样本（MWXS23186201a_R1/R2/R3）与标曲信息不进入分析输入表。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：本步耗时 < 5 s 且输出确定，每次运行都重新生成，便于复现核验。
# 版本记录：脚本末尾打印 sessionInfo()；正式运行记录写在 code/00_数据整理/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"

dir_in <- file.path(project_root, "data", "02_原始数据", "06_激素_迈维", "0.Data")
path_in_levels <- file.path(dir_in, "data.levels.xlsx")
path_in_sample <- file.path(dir_in, "sample_info.xlsx")
path_out <- file.path(project_root, "data", "01_分析输入表", "05_激素_迈维含量表.xlsx")

stopifnot(file.exists(path_in_levels), file.exists(path_in_sample))

cat("[05] 输入：", path_in_levels, "\n")
cat("[05] 输入：", path_in_sample, "\n")

# 样本信息：A/B/C 依次对应纯系 15 / 纯系 69 / 混种 15+69
readxl::read_excel(path_in_sample) %>%
  dplyr::mutate(
    replicate = as.integer(stringr::str_extract(Sample_ID, "[0-9]+$")),
    focal = dplyr::recode(Group, "A" = "15", "B" = "69", "C" = "1569")
  ) %>%
  dplyr::select(sample = Sample_ID, group = Group, focal, replicate) -> df.sample

cat("[05] 样本：", nrow(df.sample), "个\n")
print(as.data.frame(df.sample))

# 含量表：样本列固定为 sample_info 中的 12 个样本（其后为 mix/QC 与注释列）
readxl::read_excel(path_in_levels) -> df.level

sample_cols <- df.sample$sample
stopifnot(all(sample_cols %in% colnames(df.level)))

df.level %>%
  dplyr::select(
    index = Index,
    compound_en = Compounds,
    compound_cn = 物质,
    class_en = Class,
    class_cn = 物质类别,
    cpd_id = cpd_ID,
    kegg_map = kegg_map,
    dplyr::all_of(sample_cols)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_cols),
    names_to = "sample",
    values_to = "value"
  ) %>%
  dplyr::left_join(df.sample, by = "sample") %>%
  dplyr::mutate(value = as.numeric(value)) %>%
  dplyr::select(
    sample, group, focal, replicate,
    index, compound_en, compound_cn, class_en, class_cn, cpd_id, kegg_map, value
  ) -> df.out

cat("[05] 长表：", nrow(df.out), "行 x", ncol(df.out), "列（应为 109 物质 x 12 样本 = 1308 行）\n")
cat("[05] 物质类别计数：\n")
print(table(df.out$class_cn))
cat("[05] 缺失值个数：", sum(is.na(df.out$value)), "；0 值个数：", sum(df.out$value == 0, na.rm = TRUE), "\n")

writexl::write_xlsx(df.out, path_out)
cat("[05] 输出：", path_out, "\n")

cat("\n[05] sessionInfo()\n")
sessionInfo()
