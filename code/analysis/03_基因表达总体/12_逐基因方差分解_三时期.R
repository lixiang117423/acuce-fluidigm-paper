# ============================================================
# 03_基因表达总体 / 12_逐基因方差分解_三时期.R
# 作用：逐基因做方差分解，量化"发育时期 / 种植模式（邻株效应）/ 两者交互 / 品系"各解释多少表达变异
#       —— 直接回答审稿意见 M5："邻株效应在分蘖期和灌浆期各解释了多少表达变异"
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_基因表达总体/12_逐基因方差分解_三时期.R"
#
# 输入：<root>/result/03_基因表达总体/11_表_表达矩阵_三时期.xlsx（sample_matrix）
# 输出：<root>/result/03_基因表达总体/21_表_逐基因方差分解_三时期.xlsx
#       sheet：gene_variance（每个基因各效应项的类型 II 偏 R²）、terms_summary（各效应项的分布）、
#              per_stage_group（每期各基因的 group R²）、summary_by_stage
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 方法：响应 = log2(expression + 1)
#   ① 全数据模型 `lm(y ~ stage * group + focal)`，用 car::Anova(type = 2) 的偏平方和换算成
#      "该效应项解释的变异比例"（偏 R²），方便跨基因比较；
#   ② 分时期模型 `lm(y ~ group + focal)`，给出该期内 group 的 R²（邻株效应在期内的解释力）；
#   ③ 另给出 stage 主效应与 stage:group 交互的偏 R²，用于判断"响应是否随发育时期改变"。
#   说明：偏 R² 与 PERMANOVA 的 R² 不是同一口径（前者按单基因、后者按多基因距离），
#   两者应分别陈述、不要混用。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_mat <- file.path(project_root, "result", "03_基因表达总体", "11_表_表达矩阵_三时期.xlsx")
out_dir <- file.path(project_root, "result", "03_基因表达总体")
path_out <- file.path(out_dir, "21_表_逐基因方差分解_三时期.xlsx")

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

# ---- ① 全数据：类型 II 偏平方和 → 偏 R² ----
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

cat("[12] 全数据（254 样本）逐基因偏 R²（%）分布：\n")
print(as.data.frame(terms_summary))

cat("\n[12] group（邻株效应）偏 R² 最大的 10 个基因：\n")
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

cat("\n[12] 每期 group R² 最大的 5 个基因（含检出率提示）：\n")
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
