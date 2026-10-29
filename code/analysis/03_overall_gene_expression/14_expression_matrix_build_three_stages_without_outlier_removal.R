# ============================================================
# 03_overall_gene_expression / 14_expression_matrix_build_three_stages_without_outlier_removal.R
# 作用：按作者 2026-10-18 的决定——**主口径不剔除离群值**——重建"组合均值矩阵"，
#       使**图（热图/PCA/PLS-DA）and_statistics（limma）口径一致**，
#       并量化"剔除 vs 不剔除"两种组合均值矩阵的差异（供方法/补充材料引用）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_overall_gene_expression/14_expression_matrix_build_three_stages_without_outlier_removal.R"
#
# 输入：<root>/data/02_raw_data/04_Fluidigmexpression/00_引物and共享表/fluidigm.normaled.expression.3 times.xlsx
#       <root>/result/03_overall_gene_expression/11_表_expression矩阵_three_stages.xlsx（用于and"剔除后"version对照）
# 输出：<root>/result/03_overall_gene_expression/23_表_expression矩阵_three_stages_without_outlier_removal.xlsx
#       sheet：gene_annotation / sample_info / combo_matrix / sample_matrix / combo_long /
#              combo_diff（and剔除后组合均值的逐格差异）
# 其中 <root> = <repository_root>
#
# 背景：主题 03/05 的**sample级statistics（limma/camera/PERMANOVA/PLS-DA 的sample级部分）本来就使用未剔除的complete矩阵**；
#   但 `11_表_expression矩阵_three_stages.xlsx` 的**组合均值矩阵**是按 `find_outliers()` 剔除后算的（沿用旧版口径），
#   导致热图（Figure 2A/3A/4A）and_statistics口径不一致——这正是审稿意见 M6"离群值"条目要求notes的地方。
#   本脚本给出"不剔除"version的组合均值矩阵，作为**图侧的主口径**。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_in <- file.path(project_root, "data", "02_raw_data", "04_Fluidigmexpression",
                     "00_引物and共享表", "fluidigm.normaled.expression.3 times.xlsx")
path_removed <- file.path(project_root, "result", "03_overall_gene_expression", "11_表_expression矩阵_three_stages.xlsx")
out_dir <- file.path(project_root, "result", "03_overall_gene_expression")
path_out <- file.path(out_dir, "23_表_expression矩阵_three_stages_without_outlier_removal.xlsx")

stopifnot(file.exists(path_in), file.exists(path_removed))

stage_map <- c("1st" = "tillering", "2nd" = "flowering", "3rd" = "filling")
stage_levels <- c("tillering", "flowering", "filling")

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_in, "Sheet1") %>%
  magrittr::set_names(c("time", "sample", "focal", "neighbor", "biorep", "gene", "id",
                        "response", "process", "reference", "expression")) %>%
  dplyr::mutate(stage = factor(unname(stage_map[time]), levels = stage_levels),
                sample = as.character(sample),
                gene = stringr::str_remove_all(gene, "_x000D_")) -> d3

cat("[14] three_stages长表：", nrow(d3), "行；sample", dplyr::n_distinct(d3$sample),
    "；gene", dplyr::n_distinct(d3$gene), "\n")

# ---- gene注释（一gene一行，remove换行残留）----
d3 %>%
  dplyr::mutate(id = stringr::str_remove_all(id, "_x000D_")) %>%
  dplyr::select(gene, id, response, process, reference) %>%
  dplyr::distinct() %>%
  dplyr::arrange(gene, id, process) %>%
  dplyr::group_by(gene) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() -> df.gene
stopifnot(nrow(df.gene) == 92)

# ---- sample级矩阵（不剔除；and 11 表一致）----
d3 %>%
  dplyr::select(stage, sample, focal, neighbor, gene, expression) %>%
  tidyr::pivot_wider(names_from = gene, values_from = expression) -> mat.sample

gene_cols <- setdiff(colnames(mat.sample), c("stage", "sample", "focal", "neighbor"))
stopifnot(nrow(mat.sample) == 254, length(gene_cols) == 92, sum(is.na(mat.sample)) == 0)

d3 %>%
  dplyr::distinct(stage, sample, focal, neighbor) %>%
  dplyr::mutate(group = dplyr::if_else(focal == neighbor, "Pure", "Mix")) -> df.sample.info

# ---- 组合均值矩阵（不剔除）----
# 同时给出 mean and median 两种聚合：**中位数同样"不剔除任何data"**，
# 但不会被单个极端sample带飞（实测 29.55% 的格子and"剔除后均值"不同，极端差 4 个数量级）。
d3 %>%
  dplyr::group_by(stage, focal, neighbor, gene) %>%
  dplyr::summarise(mean_expression = mean(expression),
                   median_expression = stats::median(expression),
                   n_obs = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(combo = paste0(stage, "_", focal, "_", neighbor),
                group = dplyr::if_else(focal == neighbor, "Pure", "Mix")) -> df.combo.long

df.combo.long %>%
  dplyr::select(combo, stage, focal, neighbor, group, gene, mean_expression) %>%
  tidyr::pivot_wider(names_from = gene, values_from = mean_expression) -> mat.combo

stopifnot(nrow(mat.combo) == 48, !anyNA(mat.combo))
cat("[14] 组合均值矩阵（不剔除）：", nrow(mat.combo), "行 ×", length(gene_cols), "gene\n")

df.combo.long %>%
  dplyr::select(combo, stage, focal, neighbor, group, gene, median_expression) %>%
  tidyr::pivot_wider(names_from = gene, values_from = median_expression) -> mat.combo.median
stopifnot(nrow(mat.combo.median) == 48, !anyNA(mat.combo.median))

# ---- and"剔除后"version逐格对照 ----
read_sheet_raw(path_removed, "combo_matrix") -> mat.combo.removed

long_new <- as.data.frame(mat.combo) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(gene_cols), names_to = "gene", values_to = "mean_no_removal")
long_old <- as.data.frame(mat.combo.removed) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(gene_cols), names_to = "gene", values_to = "mean_removed")

combo_diff <- long_new %>%
  dplyr::left_join(long_old, by = c("combo", "gene")) %>%
  dplyr::mutate(diff = mean_no_removal - mean_removed,
                rel_diff = dplyr::if_else(mean_removed == 0, NA_real_, diff / mean_removed),
                same = is.na(diff) | abs(diff) < 1e-12)

combo_diff_summary <- combo_diff %>%
  dplyr::summarise(n_cell = dplyr::n(),
                   n_different = sum(!same),
                   pct_different = round(100 * mean(!same), 2),
                   median_abs_diff = stats::median(abs(diff), na.rm = TRUE),
                   max_abs_diff = max(abs(diff), na.rm = TRUE),
                   n_rel_diff_gt10pct = sum(abs(rel_diff) > 0.1, na.rm = TRUE))

# 中位数version vs 剔除后均值：量化"稳健聚合"and"剔除后均值"的差异
long_med <- as.data.frame(mat.combo.median) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(gene_cols), names_to = "gene", values_to = "median_no_removal")

combo_diff_median <- long_med %>%
  dplyr::left_join(long_old, by = c("combo", "gene")) %>%
  dplyr::mutate(diff = median_no_removal - mean_removed)

median_diff_summary <- combo_diff_median %>%
  dplyr::summarise(n_cell = dplyr::n(),
                   n_different = sum(abs(diff) > 1e-12),
                   pct_different = round(100 * mean(abs(diff) > 1e-12), 2),
                   median_abs_diff = stats::median(abs(diff), na.rm = TRUE),
                   max_abs_diff = max(abs(diff), na.rm = TRUE))

cat("\n[14] 组合均值：不剔除 vs 剔除 的逐格差异：\n")
print(as.data.frame(combo_diff_summary))

cat("\n[14] 组合中位数（不剔除） vs 剔除后均值 的逐格差异：\n")
print(as.data.frame(median_diff_summary))

cat("\n[14] 差异最大的 10 个 组合×gene：\n")
print(as.data.frame(combo_diff %>% dplyr::arrange(dplyr::desc(abs(diff))) %>%
                      dplyr::select(combo, gene, mean_no_removal, mean_removed, diff, rel_diff) %>%
                      dplyr::slice_head(n = 10)))

bioRtools::write_data(
  list(gene_annotation = df.gene, sample_info = df.sample.info,
       combo_matrix = mat.combo, combo_matrix_median = mat.combo.median,
       sample_matrix = mat.sample, combo_long = df.combo.long,
       combo_diff = combo_diff, combo_diff_median = combo_diff_median),
  file = path_out
)

cat("\n[14] 输出：", path_out, "\n")
cat("[14] notes：sample级矩阵and 11 表一致（均未剔除）；差异只在组合均值矩阵（图侧口径）。\n")

cat("\n[14] sessionInfo()\n")
sessionInfo()
