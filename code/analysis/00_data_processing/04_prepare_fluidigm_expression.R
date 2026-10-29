# ============================================================
# 04_prepare_fluidigm_expression.R
# 作用：把 Fluidigm 高通量 qPCR 下机data（第一批 tillering + 第三批 filling）normalized并合并成analysis输入表
# 运行：/opt/homebrew/bin/Rscript "<root>/code/00_data_processing/04_prepare_fluidigm_expression.R"
#
# 输入（本项目内只读副本）：
#   <root>/data/02_raw_data/04_Fluidigmexpression/00_引物and共享表/FLUIDIGM-Primer.xlsx
#   <root>/data/02_raw_data/04_Fluidigmexpression/01_第一批/20220923第一批FLUIDIGM.csv
#   <root>/data/02_raw_data/04_Fluidigmexpression/01_第一批/第一批96个样品选择表.xlsx
#   <root>/data/02_raw_data/04_Fluidigmexpression/03_第三批/第三批Cq值.csv
#   <root>/data/02_raw_data/04_Fluidigmexpression/03_第三批/samples.information.xlsx
# 输入上游原始路径（副本来源，未修改）：
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/field试验/FLUIDIGM-Primer.xlsx
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/field试验/20220923第一批FLUIDIGM.csv
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/field试验/第一批96个样品选择表.xlsx
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/field试验/fluidigm/第三批/第三批Cq值.csv
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/field试验/fluidigm/第三批/samples.information.xlsx
# 输出：<root>/data/processed_inputs/04_fluidigm_normalized_expression.xlsx
# 其中 <root> = <repository_root>
#
# 算法（沿用旧版，用 fluidgr 0.1 实现）：
#   1) 内参gene OseEF-1 / OseIF-4a / OsUBQ5 在每个sample内的 Ct 几何均值 norm_geom_mean；
#   2) 非内参gene expression = 2^(-(ct_value - norm_geom_mean))，保留 5 位小数（fluidgr::normalize）；
#      旧脚本另调用了 fluidgr::scale_fluidigm()，但它只新增 scaled_expression / scaled_01_expression 两列，
#      旧脚本随后 dplyr::select(... expression ...) 取的是未缩放的 expression；本脚本保留该调用以对齐流程，
#      但不使用其新增列；
#   3) 只保留 3 <= norm_geom_mean <= 21 的sample-gene组合；
#   4) 拼上样品信息（focal / neighbor）andgene注释（id / response / process / reference）。
# data格式notes：下机 csv 为逗号分隔，前 12 行为仪器信息表头，第 13 行起为data（13 列），故 fread(skip = 12)；
#              样品 0 and 431 为空白/异常孔，第三批中剔除；OsAct 为额外肌动蛋白，两批均剔除。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：本步输出确定，不做"输出已存在则跳过"判断，每次运行都重新generate，便于复现validation。
# versionrecord：脚本末尾打印 sessionInfo() and fluidgr version；正式运行record写在 code/00_data_processing/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

project_root <- normalizePath(getwd())

dir_fluidigm <- file.path(project_root, "data", "02_raw_data", "04_Fluidigmexpression")
path_primer <- file.path(dir_fluidigm, "00_引物and共享表", "FLUIDIGM-Primer.xlsx")
path_cq_1st <- file.path(dir_fluidigm, "01_第一批", "20220923第一批FLUIDIGM.csv")
path_sample_1st <- file.path(dir_fluidigm, "01_第一批", "第一批96个样品选择表.xlsx")
path_cq_3rd <- file.path(dir_fluidigm, "03_第三批", "第三批Cq值.csv")
path_sample_3rd <- file.path(dir_fluidigm, "03_第三批", "samples.information.xlsx")
path_out <- file.path(project_root, "data", "processed_inputs", "04_fluidigm_normalized_expression.xlsx")

stopifnot(
  file.exists(path_primer), file.exists(path_cq_1st), file.exists(path_sample_1st),
  file.exists(path_cq_3rd), file.exists(path_sample_3rd)
)

cat("[04] 输入：", path_primer, "\n")
cat("[04] 输入：", path_cq_1st, "\n")
cat("[04] 输入：", path_sample_1st, "\n")
cat("[04] 输入：", path_cq_3rd, "\n")
cat("[04] 输入：", path_sample_3rd, "\n")

# Fluidigm 下机 csv 的 13 列字段名（第 13 行起为data）
col_names_fluidigm <- c(
  "chamber_id", "sample_name", "sample_type", "sample_rconc",
  "target_name", "target_type", "ct_value", "ct_quality", "ct_call",
  "ct_threshold", "tm_inrange", "tm_outrange", "tm_peakratio"
)

# gene名称归一（下机data的gene名and引物表不完全一致，统一为 Os 前缀写法）
fix_target_name <- function(x) {
  dplyr::case_when(
    x == "OsGlut dehydr" ~ "OsGlut_dehydr",
    x == "OsPR1B" ~ "OsPR1b",
    x == "OsMT2" ~ "OsMT2b",
    x == "OsActine" ~ "OsAct",
    x == "OsLIPoxy" ~ "OsLIPOxy",
    x == "OsFd-GOGAT =  GLU1" ~ "OsGLU1",
    x == "Os08g0126300" ~ "OsGAPC3",
    x == "Os09g0315700" ~ "OsPEPC",
    x == "Os02g0601300" ~ "OsGAPC1",
    x == "Os04g0479200" ~ "OsNAD-isocitrate2",
    x == "Os02g0194100" ~ "OsCitr_synth",
    x == "Os02g0739600" ~ "OsPyr_dehydra",
    x == "Os01g0654500" ~ "OsNADP-isocitrate",
    TRUE ~ x
  )
}

# gene注释：Gene ID / Gene name / Expected response / process2 / Reference
readxl::read_excel(path_primer, sheet = 1) %>%
  dplyr::select(1:3, 5, 8) %>%
  magrittr::set_names(c("id", "name", "response", "process", "reference")) %>%
  dplyr::mutate(
    name = stringr::str_replace(name, "Os", ""),
    name = dplyr::case_when(
      name == "AOS2" ~ "AOS2",
      TRUE ~ stringr::str_replace(name, "OS", "")
    ),
    name = paste0("Os", name)
  ) %>%
  dplyr::mutate(name = dplyr::case_when(
    name == "OseEF-1α" ~ "OseEF-1",
    TRUE ~ name
  )) %>%
  dplyr::rename(target_name = name) -> df.gene.info

cat("[04] gene注释：", nrow(df.gene.info), "个gene\n")

# 内参gene
normalizers <- c("OseEF-1", "OseIF-4a", "OsUBQ5")
stopifnot(all(normalizers %in% df.gene.info$target_name))

# ---------------- 第一批（tillering 期） ----------------

readxl::read_excel(path_sample_1st, sheet = 1) %>%
  dplyr::select(tag, focal, neighbor, BioRep) %>%
  dplyr::rename(sample_name = tag) -> df.sample.info.1st

data.table::fread(path_cq_1st, header = FALSE, skip = 12) %>%
  as.data.frame() %>%
  magrittr::set_names(col_names_fluidigm) %>%
  dplyr::mutate(
    target_name = stringr::str_replace(target_name, "Os", ""),
    target_name = paste0("Os", target_name)
  ) %>%
  dplyr::mutate(target_name = fix_target_name(target_name)) %>%
  dplyr::filter(target_name != "OsAct") -> df.cq.1st

cat("[04] 第一批下机data：", nrow(df.cq.1st), "行；sample数：", dplyr::n_distinct(df.cq.1st$sample_name), "\n")

df.cq.1st %>%
  fluidgr::normalize(normalizers = normalizers) %>%
  fluidgr::scale_fluidigm(.group = target_name) %>%
  dplyr::filter(norm_geom_mean <= 21, norm_geom_mean >= 3) %>%
  dplyr::select(sample_name, target_name, expression) %>%
  dplyr::left_join(df.sample.info.1st, by = "sample_name") %>%
  dplyr::left_join(df.gene.info, by = "target_name") %>%
  dplyr::select(1, 4, 5, 6, 2, 7:10, 3) %>%
  dplyr::rename(sample = sample_name, gene = target_name) %>%
  dplyr::mutate(focal = as.character(focal), neighbor = as.character(neighbor)) -> df.norm.1st

# ---------------- 第三批（filling 期） ----------------

readxl::read_excel(path_sample_3rd, sheet = 1) %>%
  dplyr::rename(sample_name = sample) %>%
  dplyr::mutate(sample_name = as.character(sample_name)) -> df.sample.info.3rd

data.table::fread(path_cq_3rd, header = FALSE, skip = 12) %>%
  as.data.frame() %>%
  magrittr::set_names(col_names_fluidigm) %>%
  dplyr::mutate(
    sample_name = as.character(sample_name),
    target_name = stringr::str_replace(target_name, "Os", ""),
    target_name = paste0("Os", target_name)
  ) %>%
  dplyr::mutate(target_name = fix_target_name(target_name)) %>%
  dplyr::filter(target_name != "OsAct", !sample_name %in% c("0", "431")) -> df.cq.3rd

cat("[04] 第三批下机data：", nrow(df.cq.3rd), "行；sample数：", dplyr::n_distinct(df.cq.3rd$sample_name), "\n")

df.cq.3rd %>%
  fluidgr::normalize(normalizers = normalizers) %>%
  fluidgr::scale_fluidigm(.group = target_name) %>%
  dplyr::filter(norm_geom_mean <= 21, norm_geom_mean >= 3) %>%
  dplyr::select(sample_name, target_name, expression) %>%
  dplyr::left_join(df.sample.info.3rd, by = "sample_name") %>%
  dplyr::left_join(df.gene.info, by = "target_name") %>%
  dplyr::select(1, 4, 5, 6, 2, 7:9, 3) %>%
  dplyr::rename(sample = sample_name, gene = target_name) %>%
  dplyr::mutate(focal = as.character(focal), neighbor = as.character(neighbor)) -> df.norm.3rd

# ---------------- 合并两个批次并导出 ----------------

df.norm.1st %>%
  dplyr::select(dplyr::all_of(colnames(df.norm.3rd))) %>%
  dplyr::mutate(stage = "tillering") %>%
  dplyr::select(stage, dplyr::all_of(colnames(df.norm.3rd))) %>%
  dplyr::mutate(sample = as.character(sample)) -> df.1

df.norm.3rd %>%
  dplyr::mutate(stage = "filling") %>%
  dplyr::select(dplyr::all_of(colnames(df.1))) %>%
  dplyr::mutate(sample = as.character(sample)) %>%
  dplyr::bind_rows(df.1) -> df.out

cat("[04] tillering：", nrow(df.1), "行；filling：", nrow(df.norm.3rd), "行\n")
cat("[04] 合并后：", nrow(df.out), "行 x", ncol(df.out), "列\n")
cat("[04] 列名：", paste(colnames(df.out), collapse = ", "), "\n")
cat("[04] stage x sample数 / gene数：\n")
print(df.out %>%
  dplyr::group_by(stage) %>%
  dplyr::summarise(n_sample = dplyr::n_distinct(sample), n_gene = dplyr::n_distinct(gene), .groups = "drop") %>%
  as.data.frame())

writexl::write_xlsx(df.out, path_out)
cat("[04] 输出：", path_out, "\n")

cat("\n[04] fluidgr version：", as.character(utils::packageVersion("fluidgr")), "\n")
cat("[04] sessionInfo()\n")
sessionInfo()
