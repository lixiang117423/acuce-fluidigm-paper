# ============================================================
# 03_overall_gene_expression / 13_outlier_sensitivity_three_stages.R
# 作用：回答审稿意见 M6 的"离群值"条目——notes剔除规则，并做**剔除前后的sensitivityanalysis**
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_overall_gene_expression/13_outlier_sensitivity_three_stages.R"
#
# 输入：<root>/result/03_overall_gene_expression/11_表_expression矩阵_three_stages.xlsx（sample_matrix）
# 输出：<root>/result/03_overall_gene_expression/22_表_outlier_sensitivity_three_stages.xlsx
#       sheet：outlier_cells（被标记的单元格）、outlier_summary（按 时期/gene/组合 汇总）、
#              limma_compare（逐gene 剔除前 vs 剔除后）、compare_summary（每期汇总）
# 其中 <root> = <repository_root>
#
# 口径：
#   - 现行statistics（limma）用的是**未剔离群值**的sample级矩阵（complete、无缺失）；
#   - 离群值规则：`bioRtools::find_outliers()`，按 **stage × focal × neighbor × gene** 分组判定
#     （and主题 03 脚本 01 的组合均值矩阵口径一致）；
#   - sensitivityanalysis：把这些单元格置为 NA 后**重跑同一 limma 模型**（`~ 0 + group + focal`，
#     显式 Mix − Pure 对比，`eBayes(robust = TRUE)`；limma 允许逐gene缺失），comparison两版results。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_mat <- file.path(project_root, "result", "03_overall_gene_expression", "11_表_expression矩阵_three_stages.xlsx")
out_dir <- file.path(project_root, "result", "03_overall_gene_expression")
path_out <- file.path(out_dir, "22_表_outlier_sensitivity_three_stages.xlsx")

stopifnot(file.exists(path_mat))

stage_levels <- c("tillering", "flowering", "filling")

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_mat, "sample_matrix") -> mat.sample
meta_cols <- c("stage", "sample", "focal", "neighbor")
gene_cols <- setdiff(colnames(mat.sample), meta_cols)

d <- mat.sample %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels),
                group = factor(dplyr::if_else(focal == neighbor, "Pure", "Mix"),
                               levels = c("Pure", "Mix")),
                focal = factor(focal))

# ---- 标记离群单元格（stage × focal × neighbor × gene）----
long <- d %>%
  dplyr::select(stage, sample, focal, neighbor, group, dplyr::all_of(gene_cols)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(gene_cols), names_to = "gene", values_to = "expression") %>%
  dplyr::group_by(stage, focal, neighbor, gene) %>%
  dplyr::mutate(is_outlier = suppressWarnings(bioRtools::find_outliers(expression)) != "no") %>%
  dplyr::ungroup()

outlier_cells <- long %>% dplyr::filter(is_outlier)

outlier_summary <- long %>%
  dplyr::group_by(stage) %>%
  dplyr::summarise(n_cell = dplyr::n(), n_outlier = sum(is_outlier),
                   pct = round(100 * mean(is_outlier), 2), .groups = "drop")

cat("[13] 离群单元格（规则：stage × focal × neighbor × gene 内 find_outliers）：\n")
print(as.data.frame(outlier_summary))

cat("\n[13] 被标记最多的 10 个gene：\n")
print(as.data.frame(outlier_cells %>% dplyr::count(gene, sort = TRUE) %>% dplyr::slice_head(n = 10)))

# ---- 剔除后矩阵（置 NA）----
wide_mark <- long %>%
  dplyr::mutate(expression_masked = dplyr::if_else(is_outlier, NA_real_, expression)) %>%
  dplyr::select(stage, sample, focal, neighbor, gene, expression_masked) %>%
  tidyr::pivot_wider(names_from = gene, values_from = expression_masked)

# ---- 两版 limma 对照（逐期、逐gene）----
run_limma <- function(dd, cols) {
  expr <- t(log2(as.matrix(dd[, cols]) + 1))
  design <- stats::model.matrix(~ 0 + group + focal, data = dd)
  colnames(design) <- make.names(colnames(design))
  con <- matrix(0, nrow = ncol(design), ncol = 1,
                dimnames = list(colnames(design), "mix_vs_pure"))
  con["groupMix", "mix_vs_pure"] <- 1
  con["groupPure", "mix_vs_pure"] <- -1
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr, design), con), robust = TRUE)
  limma::topTable(fit, coef = "mix_vs_pure", number = Inf, sort.by = "P")
}

compare <- purrr::map_dfr(stage_levels, function(st) {
  d1 <- d %>% dplyr::filter(stage == st)
  d2 <- wide_mark %>% dplyr::filter(stage == st) %>%
    dplyr::mutate(stage = factor(stage, levels = stage_levels),
                  group = factor(dplyr::if_else(focal == neighbor, "Pure", "Mix"),
                                 levels = c("Pure", "Mix")),
                  focal = factor(focal))

  t1 <- run_limma(d1, gene_cols)
  t2 <- run_limma(d2, gene_cols)

  tibble::tibble(
    stage = st, gene = rownames(t1),
    logFC_full = t1$logFC, p_full = t1$P.Value, fdr_full = t1$adj.P.Val,
    logFC_noOut = t2[rownames(t1), "logFC"], p_noOut = t2[rownames(t1), "P.Value"],
    fdr_noOut = t2[rownames(t1), "adj.P.Val"]
  )
})

compare_summary <- compare %>%
  dplyr::group_by(stage) %>%
  dplyr::summarise(
    n_gene = dplyr::n(),
    n_nominal_full = sum(p_full < 0.05),
    n_nominal_noOut = sum(p_noOut < 0.05, na.rm = TRUE),
    n_fdr_full = sum(fdr_full < 0.05),
    n_fdr_noOut = sum(fdr_noOut < 0.05, na.rm = TRUE),
    pearson_r_logFC = stats::cor(logFC_full, logFC_noOut, use = "complete.obs"),
    spearman_r_p = stats::cor(p_full, p_noOut, method = "spearman", use = "complete.obs"),
    median_abs_logFC_diff = stats::median(abs(logFC_full - logFC_noOut), na.rm = TRUE),
    .groups = "drop"
  )

cat("\n[13] 剔除前 vs 剔除后（逐期）：\n")
print(as.data.frame(compare_summary))

cat("\n[13] 被剔除离群值后显著性状态改变的gene：\n")
print(as.data.frame(compare %>%
                      dplyr::mutate(change = dplyr::case_when(
                        p_full < 0.05 & (is.na(p_noOut) | p_noOut >= 0.05) ~ "lost nominal",
                        p_full >= 0.05 & !is.na(p_noOut) & p_noOut < 0.05 ~ "gained nominal",
                        TRUE ~ "unchanged")) %>%
                      dplyr::filter(change != "unchanged") %>%
                      dplyr::select(stage, gene, logFC_full, p_full, p_noOut, change)))

bioRtools::write_data(
  list(outlier_cells = outlier_cells %>% dplyr::select(stage, sample, focal, neighbor, gene, expression),
       outlier_summary = outlier_summary,
       limma_compare = compare, compare_summary = compare_summary),
  file = path_out
)

cat("\n[13] 输出：", path_out, "\n")

cat("\n[13] sessionInfo()\n")
sessionInfo()
