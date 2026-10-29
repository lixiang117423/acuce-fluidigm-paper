# ============================================================
# 03_overall_gene_expression / 12_per_gene_variance_decomposition_three_stages.R
# 作用：逐gene做variance_decomposition，量化"发育时期 / 种植模式（邻株effect）/ 两者交互 / 品系"各解释多少expression变异
#       —— 直接回答审稿意见 M5："邻株effect在分蘖期和灌浆期各解释了多少expression变异"
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_overall_gene_expression/12_per_gene_variance_decomposition_three_stages.R"
#
# 输入：<root>/result/03_overall_gene_expression/11_表_expression矩阵_three_stages.xlsx（sample_matrix）
# 输出：<root>/result/03_overall_gene_expression/21_表_per_gene_variance_decomposition_three_stages.xlsx
#       sheet：gene_variance（每个gene各effect项的类型 II 偏 R²）、terms_summary（各effect项的分布）、
#              per_stage_group（每期各gene的 group R²）、summary_by_stage
# 其中 <root> = <repository_root>
#
# 方法：响应 = log2(expression + 1)
#   ① 全data模型 `lm(y ~ stage * group + focal)`，用 car::Anova(type = 2) 的偏平方和换算成
#      "该effect项解释的变异比例"（偏 R²），方便跨genecomparison；
#   ② 分时期模型 `lm(y ~ group + focal)`，给出该期内 group 的 R²（邻株effect在期内的解释力）；
#   ③ 另给出 stage 主effectand stage:group 交互的偏 R²，用于判断"响应是否随发育时期改变"。
#   notes：偏 R² and PERMANOVA 的 R² 不是同一口径（前者按单gene、后者按多gene距离），
#   两者应分别陈述、不要混用。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_mat <- file.path(project_root, "result", "03_overall_gene_expression", "11_表_expression矩阵_three_stages.xlsx")
out_dir <- file.path(project_root, "result", "03_overall_gene_expression")
path_out <- file.path(out_dir, "21_表_per_gene_variance_decomposition_three_stages.xlsx")

stopifnot(file.exists(path_mat))

stage_levels <- c("tillering", "flowering", "filling")

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_mat, "sample_matrix") -> mat.sample
meta_cols <- c("stage", "sample", "focal", "neighbor")
gene_cols <- setdiff(colnames(mat.sample), meta_cols)
stopifnot(length(gene_cols) == 92, nrow(mat.sample) == 254)

d <- mat.sample %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels),
                group = factor(dplyr::if_else(focal == neighbor, "Pure", "Mix"),
                               levels = c("Pure", "Mix")),
                focal = factor(focal))

# ---- ① 全data：类型 II 偏平方和 → 偏 R² ----
gene_variance <- purrr::map_dfr(gene_cols, function(g) {
  dd <- data.frame(y = log2(d[[g]] + 1), stage = d$stage, group = d$group, focal = d$focal)
  m <- stats::lm(y ~ stage * group + focal, data = dd)
  a <- as.data.frame(car::Anova(m, type = 2))
  a$term <- rownames(a)
  a <- a[!is.na(a$`Sum Sq`), ]
  ss <- stats::setNames(a$`Sum Sq`, a$term)
  total <- sum(ss)                       # 含残差
  tibble::tibble(
    gene = g,
    r2_stage = 100 * ss[["stage"]] / total,
    r2_group = 100 * ss[["group"]] / total,
    r2_stage_group = 100 * (if ("stage:group" %in% names(ss)) ss[["stage:group"]] else 0) / total,
    r2_focal = 100 * ss[["focal"]] / total,
    r2_residual = 100 * ss[["Residuals"]] / total,
    ss_stage = ss[["stage"]], ss_group = ss[["group"]],
    ss_stage_group = if ("stage:group" %in% names(ss)) ss[["stage:group"]] else 0,
    ss_focal = ss[["focal"]], ss_residual = ss[["Residuals"]], ss_total = total
  )
}) %>%
  dplyr::left_join(read_sheet_raw(path_mat, "gene_annotation") %>% dplyr::select(gene, process),
                   by = "gene")

terms_summary <- gene_variance %>%
  dplyr::summarise(
    dplyr::across(dplyr::starts_with("r2_"),
                  list(median = ~ stats::median(.x),
                       mean = ~ mean(.x),
                       q25 = ~ stats::quantile(.x, 0.25),
                       q75 = ~ stats::quantile(.x, 0.75),
                       max = ~ max(.x)),
                  .names = "{.col}__{.fn}")
  ) %>%
  tidyr::pivot_longer(dplyr::everything(), names_to = "stat", values_to = "value") %>%
  tidyr::separate(stat, into = c("term", "stat"), sep = "__") %>%
  tidyr::pivot_wider(names_from = stat, values_from = value)

cat("[12] 全data（254 sample）逐gene偏 R²（%）分布：\n")
print(as.data.frame(terms_summary))

cat("\n[12] group（邻株effect）偏 R² 最大的 10 个gene：\n")
print(as.data.frame(gene_variance %>% dplyr::arrange(dplyr::desc(r2_group)) %>%
                      dplyr::select(gene, process, r2_stage, r2_group, r2_stage_group, r2_focal) %>%
                      dplyr::slice_head(n = 10)))

# ---- ② 分时期：group 的 R² ----
per_stage_group <- purrr::map_dfr(stage_levels, function(st) {
  dd0 <- d %>% dplyr::filter(stage == st)
  purrr::map_dfr(gene_cols, function(g) {
    dd <- data.frame(y = log2(dd0[[g]] + 1), group = dd0$group, focal = dd0$focal)
    m0 <- stats::lm(y ~ 1, data = dd)
    m1 <- stats::lm(y ~ group + focal, data = dd)
    r2 <- 100 * (sum(stats::residuals(m0)^2) - sum(stats::residuals(m1)^2)) /
      sum(stats::residuals(m0)^2)
    # group 单独的偏 R²（在已含 focal 的模型里加入 group）
    m2 <- stats::lm(y ~ focal, data = dd)
    r2_group <- 100 * (sum(stats::residuals(m2)^2) - sum(stats::residuals(m1)^2)) /
      sum(stats::residuals(m0)^2)
    tibble::tibble(stage = st, gene = g, r2_model = r2, r2_group = r2_group,
                   n_sample = nrow(dd))
  })
})

summary_by_stage <- per_stage_group %>%
  dplyr::group_by(stage) %>%
  dplyr::summarise(n_gene = dplyr::n(),
                   median_r2_group = stats::median(r2_group),
                   mean_r2_group = mean(r2_group),
                   q75_r2_group = stats::quantile(r2_group, 0.75),
                   max_r2_group = max(r2_group),
                   n_gene_r2_gt5 = sum(r2_group > 5),
                   .groups = "drop")

cat("\n[12] 分时期：group 的偏 R²（%）分布：\n")
print(as.data.frame(summary_by_stage))

cat("\n[12] 每期 group R² 最大的 5 个gene（含检出率提示）：\n")
top_by_stage <- per_stage_group %>%
  dplyr::group_by(stage) %>%
  dplyr::arrange(dplyr::desc(r2_group)) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(gene_variance %>% dplyr::select(gene, process), by = "gene")
print(as.data.frame(top_by_stage))

bioRtools::write_data(
  list(gene_variance = gene_variance, terms_summary = terms_summary,
       per_stage_group = per_stage_group, summary_by_stage = summary_by_stage,
       top_by_stage = top_by_stage),
  file = path_out
)

cat("\n[12] 输出：", path_out, "\n")

cat("\n[12] sessionInfo()\n")
sessionInfo()
