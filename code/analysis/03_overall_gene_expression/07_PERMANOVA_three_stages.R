# ============================================================
# 03_overall_gene_expression / 07_PERMANOVA_three_stages.R
# 作用：把 PERMANOVA（Table S2）扩展到**三个发育时期**（分蘖 / 花期 / 灌浆）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_overall_gene_expression/07_PERMANOVA_three_stages.R"
#
# 输入：<root>/result/03_overall_gene_expression/11_表_expression矩阵_three_stages.xlsx（combo_matrix / sample_matrix / sample_info）
# 输出：<root>/result/03_overall_gene_expression/12_表_PERMANOVA_three_stages.xlsx
# 其中 <root> = <repository_root>
#
# 口径：and脚本 02（两时期 Table S2）完全一致，仅时期数由 2 变 3：
#   1) 两套data：组合均值矩阵（48 行 = 3 时期 × 16 组合，已剔离群值）andsample级矩阵（254 行，completedata）；
#   2) 两种距离：列（gene）Z-score 后的欧氏距离；log2(x+1) 后的 Bray–Curtis；
#   3) 模型 `~ stage + group + focal + stage:group`，by = "terms" and by = "margin" 各跑一次，permutations = 999；
#   4) 另按时期分层单独建模（`~ group + focal`）。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_in <- file.path(project_root, "result", "03_overall_gene_expression", "11_表_expression矩阵_three_stages.xlsx")
out_dir <- file.path(project_root, "result", "03_overall_gene_expression")
path_out <- file.path(out_dir, "12_表_PERMANOVA_three_stages.xlsx")

stopifnot(file.exists(path_in))

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_in, "combo_matrix") -> mat.combo
read_sheet_raw(path_in, "sample_matrix") -> mat.sample
bioRtools::read_data(path_in, sheet = "sample_info") -> df.info

cat("[07] 组合均值矩阵：", dim(mat.combo), "；sample级矩阵：", dim(mat.sample), "\n")

make_dist <- function(mat, method = c("euclidean_scaled", "bray_log")) {
  method <- match.arg(method)
  x <- as.matrix(mat)
  if (method == "euclidean_scaled") {
    x <- scale(x, center = TRUE, scale = TRUE)
    stats::dist(x, method = "euclidean")
  } else {
    vegan::vegdist(log2(x + 1), method = "bray")
  }
}

run_adonis <- function(dist_obj, meta, rhs, by = c("terms", "margin")) {
  by <- match.arg(by)
  d <- dist_obj
  df <- as.data.frame(meta)
  form <- stats::as.formula(paste("d ~", rhs))
  res <- vegan::adonis2(form, data = df, permutations = 999, by = by)
  as.data.frame(res) %>%
    tibble::rownames_to_column(var = "term") %>%
    dplyr::mutate(by = by, formula = paste("d ~", rhs), n_obs = nrow(df))
}

df.combo.meta <- mat.combo %>%
  dplyr::select(combo, stage, focal, neighbor, group) %>%
  as.data.frame()

df.sample.meta <- mat.sample %>%
  dplyr::select(stage, sample, focal, neighbor) %>%
  dplyr::mutate(group = dplyr::if_else(focal == neighbor, "Pure", "Mix")) %>%
  as.data.frame()

stopifnot(nrow(df.sample.meta) == nrow(mat.sample))

dists <- list(
  combo_euclidean = make_dist(dplyr::select(mat.combo, -combo, -stage, -focal, -neighbor, -group),
                              "euclidean_scaled"),
  combo_bray = make_dist(dplyr::select(mat.combo, -combo, -stage, -focal, -neighbor, -group),
                         "bray_log"),
  sample_euclidean = make_dist(dplyr::select(mat.sample, -stage, -sample, -focal, -neighbor),
                               "euclidean_scaled"),
  sample_bray = make_dist(dplyr::select(mat.sample, -stage, -sample, -focal, -neighbor),
                          "bray_log")
)

res_list <- list()

for (nm in names(dists)) {
  meta <- if (grepl("^combo", nm)) df.combo.meta else df.sample.meta
  d <- dists[[nm]]

  res_list[[paste0(nm, "_terms")]] <- run_adonis(d, meta, "stage + group + focal + stage:group", "terms") %>%
    dplyr::mutate(dataset = nm, distance = nm)
  res_list[[paste0(nm, "_margin")]] <- run_adonis(d, meta, "stage + group + focal + stage:group", "margin") %>%
    dplyr::mutate(dataset = nm, distance = nm)

  for (st in unique(meta$stage)) {
    sub_idx <- which(meta$stage == st)
    sub_dist <- as.matrix(d)[sub_idx, sub_idx] %>% stats::as.dist()
    sub_meta <- meta[sub_idx, , drop = FALSE]
    res_list[[paste0(nm, "_", st)]] <- run_adonis(sub_dist, sub_meta, "group + focal", "terms") %>%
      dplyr::mutate(dataset = paste0(nm, " | ", st), distance = nm, stage_subset = st)
  }
}

permanova_tab <- purrr::list_rbind(res_list) %>%
  dplyr::select(dataset, distance, by, formula, stage_subset, term, Df, SumOfSqs, R2, F, `Pr(>F)`, n_obs) %>%
  dplyr::rename(p_value = `Pr(>F)`) %>%
  dplyr::mutate(R2_percent = 100 * R2)

cat("\n[07] 主results（by = terms，两套data × 两种距离）：\n")
print(permanova_tab %>%
        dplyr::filter(by == "terms", is.na(stage_subset)) %>%
        dplyr::select(dataset, term, Df, R2_percent, F, p_value) %>%
        as.data.frame(), row.names = FALSE)

cat("\n[07] 边际test（by = margin）：\n")
print(permanova_tab %>%
        dplyr::filter(by == "margin", is.na(stage_subset)) %>%
        dplyr::select(dataset, term, Df, R2_percent, F, p_value) %>%
        as.data.frame(), row.names = FALSE)

cat("\n[07] 分时期分层模型（by = terms）：\n")
print(permanova_tab %>%
        dplyr::filter(!is.na(stage_subset)) %>%
        dplyr::select(dataset, term, Df, R2_percent, F, p_value) %>%
        as.data.frame(), row.names = FALSE)

bioRtools::write_data(list(permanova = permanova_tab, sample_info = df.sample.meta), file = path_out)
cat("\n[07] 输出：", path_out, "\n")

cat("\n[07] sessionInfo()\n")
sessionInfo()
