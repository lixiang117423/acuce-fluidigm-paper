# ============================================================
# 03_overall_gene_expression / 02_PERMANOVA_TableS2.R
# 作用：PERMANOVA（vegan::adonis2）test发育时期、种植模式and品系对整体expression格局的影响 —— Table S2
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_overall_gene_expression/02_PERMANOVA_TableS2.R"
#
# 输入：<root>/result/03_overall_gene_expression/01_表_expression矩阵.xlsx（combo_matrix / sample_matrix / sample_info）
# 输出：<root>/result/03_overall_gene_expression/02_表_PERMANOVA_TableS2.xlsx
# 其中 <root> = <repository_root>
#
# 重要notes：旧版脚本**从未运行过 PERMANOVA**（两份 qmd 全文都没有 adonis/vegan 调用），
#   paper框架里"stage R² = 0.53、focal R² = 0.02、互作 R² = 0.02 且不显著"等数字在代码中没有出处。
#   本脚本按 2026-09-22 确认的方案首次实现，结论以本脚本输出为准。
# 方法：
#   1) 两套data：组合均值矩阵（32 行 = 2 时期 × 16 组合，已剔离群值，and PLS-DA 输入一致）andsample级矩阵（168 行，completedata）；
#   2) 两种距离：gene先按行（sample）normalized后取欧氏距离；以及 log2(x+1) 后的 Bray–Curtis 距离；
#   3) 模型：`~ stage + group + focal + stage:group`，`by = "terms"`（顺序test）and `by = "margin"`（边际test）各跑一次，
#      permutations = 999 且固定随机种子；
#   4) 由于 Pure : Mix sample数不平衡（组合均值 8 : 24、sample级 43 : 125），另按时期分层单独建模（`~ group + focal`）。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# 读取notes：矩阵 sheet 必须用 readxl::read_excel(..., .name_repair = "minimal") 读取，
#   否则gene名里的连字符会被改写成点号（OsD14-1 → OsD14.1），下游按gene名匹配会失败。
# 断点续传：单步耗时 < 120 s，每次运行都重新generate。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/03_overall_gene_expression/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_in <- file.path(project_root, "result", "03_overall_gene_expression", "01_表_expression矩阵.xlsx")
out_dir <- file.path(project_root, "result", "03_overall_gene_expression")
path_out <- file.path(out_dir, "02_表_PERMANOVA_TableS2.xlsx")

stopifnot(file.exists(path_in))

# 矩阵 sheet 用 minimal 列名修复，保留gene名原样
read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_in, "combo_matrix") -> mat.combo
read_sheet_raw(path_in, "sample_matrix") -> mat.sample
bioRtools::read_data(path_in, sheet = "sample_info") -> df.info

cat("[02] 组合均值矩阵：", dim(mat.combo), "；sample级矩阵：", dim(mat.sample), "\n")

# ---- 距离矩阵工具 ----
make_dist <- function(mat, method = c("euclidean_scaled", "bray_log")) {
  method <- match.arg(method)
  x <- as.matrix(mat)
  if (method == "euclidean_scaled") {
    # gene（列）先做 Z-score，避免高expressiongene主导距离
    x <- scale(x, center = TRUE, scale = TRUE)
    stats::dist(x, method = "euclidean")
  } else {
    # Bray–Curtis 用 vegan::vegdist（stats::dist 不支持 bray）
    vegan::vegdist(log2(x + 1), method = "bray")
  }
}

# ---- PERMANOVA 工具 ----
# 注意：adonis2 若收到预先算好的距离矩阵，公式必须带左侧（d ~ ...），
#   且 d 要在 as.formula() 所处的环境里可见，否则会报 "找不到对象 'group'"。
run_adonis <- function(dist_obj, meta, rhs, by = c("terms", "margin")) {
  by <- match.arg(by)
  d <- dist_obj
  df <- as.data.frame(meta)
  form <- stats::as.formula(paste("d ~", rhs))
  res <- vegan::adonis2(form, data = df, permutations = 999, by = by)
  as.data.frame(res) %>%
    tibble::rownames_to_column(var = "term") %>%
    dplyr::mutate(
      by = by,
      formula = paste("d ~", rhs),
      n_obs = nrow(df)
    )
}

# 组合均值矩阵的元信息
df.combo.meta <- mat.combo %>%
  dplyr::select(combo, stage, focal, neighbor, group) %>%
  as.data.frame()

# sample级矩阵的元信息（and sample_info 对齐）
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
  if (grepl("^combo", nm)) {
    meta <- df.combo.meta
  } else {
    meta <- df.sample.meta
  }
  d <- dists[[nm]]

  res_list[[paste0(nm, "_terms")]] <- run_adonis(d, meta, "stage + group + focal + stage:group", "terms") %>%
    dplyr::mutate(dataset = nm, distance = nm)

  res_list[[paste0(nm, "_margin")]] <- run_adonis(d, meta, "stage + group + focal + stage:group", "margin") %>%
    dplyr::mutate(dataset = nm, distance = nm)

  # 分时期分层模型
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

cat("\n[02] PERMANOVA 主results（按 terms 顺序test，两套data × 两种距离）：\n")
permanova_tab %>%
  dplyr::filter(by == "terms", is.na(stage_subset)) %>%
  dplyr::select(dataset, term, Df, R2_percent, F, p_value) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n[02] 分时期分层模型：\n")
permanova_tab %>%
  dplyr::filter(!is.na(stage_subset)) %>%
  dplyr::select(dataset, term, Df, R2_percent, F, p_value) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

bioRtools::write_data(
  list(permanova = permanova_tab, sample_info = df.sample.meta),
  file = path_out
)
cat("[02] 输出：", path_out, "\n")

cat("\n[02] sessionInfo()\n")
sessionInfo()
