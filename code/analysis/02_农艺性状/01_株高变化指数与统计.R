# ============================================================
# 02_农艺性状 / 01_株高变化指数与统计.R
# 作用：株高变化指数（Mix − Pure，按 focal × time 的纯种均值定义）的计算与统计检验
#       —— Figure 1A 的数据与统计输入
# 运行：/opt/homebrew/bin/Rscript "<root>/code/02_农艺性状/01_株高变化指数与统计.R"
#
# 输入：<root>/data/01_分析输入表/03_分蘖数_株高_单株产量.xlsx（8095 行 × 6 列）
# 输出：<root>/result/02_农艺性状/01_表_株高变化指数.xlsx（多 sheet）
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 数据口径（实测核对 2026-09-22）：
#   1) 株高调查共 15 个时间点（20220711–20221031）；表中 2518 行 time 为空，分析一律剔除（与旧版一致）；
#      focal 36 / 74 的全部株高记录 time 均为空，因此自动不进入本主题分析；
#   2) Pure 组的 neighbor 已统一等于 focal（整理时设定），故"纯种均值"按 focal × time 聚合；
#   3) 株高变化指数 index = value − 同 focal × time 的 Pure 均值；Pure 组自身 index 的均值恒为 0，
#      作为对照水平保留在模型与图中（与旧版一致）。
# 统计口径：① 总体检验 —— 混种（Mix）index 是否显著偏离 0（单样本 t 检验 + Wilcoxon 对照）；
#          ② 每个 focal 的邻体比较 —— 旧口径 `bioRtools::anova_posthoc()`（单因素 ANOVA + Tukey 字母，保证与旧图可比）
#          与 LMM `lmer(index ~ neighbor + (1 | time))` + emmeans(Sidak) 两套并列，作为稳健性对照；
#          ③ LMM 的随机截距为时间点（15 个水平），用于吸收时间维度残余变异。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：单步耗时 < 60 s，每次运行都重新生成。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录见 code/02_农艺性状/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_in <- file.path(project_root, "data", "01_分析输入表", "03_分蘖数_株高_单株产量.xlsx")
out_dir <- file.path(project_root, "result", "02_农艺性状")
path_out <- file.path(out_dir, "01_表_株高变化指数.xlsx")

stopifnot(file.exists(path_in))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("[01] 输入：", path_in, "\n")

bioRtools::read_data(path_in) -> d3
stopifnot(nrow(d3) == 8095, ncol(d3) == 6)

d3 %>%
  dplyr::filter(phenotype == "plant height", !is.na(time)) %>%
  dplyr::mutate(
    focal = as.factor(focal),
    neighbor = as.factor(neighbor),
    time = as.factor(time)
  ) -> df

cat("[01] 株高（有时间）行数：", nrow(df), "；focal：", paste(levels(df$focal), collapse = ","),
    "；时间点数：", dplyr::n_distinct(df$time), "\n")

# ---- 纯种均值（focal × time）----
df %>%
  dplyr::filter(group == "Pure") %>%
  dplyr::group_by(focal, time) %>%
  dplyr::summarise(mean_pure = mean(value), n_pure = dplyr::n(), .groups = "drop") -> df.pure

stopifnot(!anyNA(df.pure$mean_pure))

# ---- 变化指数 ----
df %>%
  dplyr::left_join(df.pure, by = c("focal", "time")) %>%
  dplyr::mutate(index = value - mean_pure) -> df.idx

stopifnot(!anyNA(df.idx$index))

cat("[01] index 取值范围：", round(min(df.idx$index), 2), " ~ ", round(max(df.idx$index), 2), "\n")

# ---- 总体检验：Mix 的 index 是否偏离 0 ----
df.idx %>%
  dplyr::filter(group == "Mix") -> df.mix

t.test(df.mix$index, mu = 0) -> tt
wilcox.test(df.mix$index, mu = 0) -> wt

tibble::tibble(
  test = c("one_sample_t", "wilcoxon"),
  n = c(nrow(df.mix), nrow(df.mix)),
  mean_index = c(mean(df.mix$index), mean(df.mix$index)),
  estimate = c(unname(tt$estimate), unname(wt$estimate)),
  ci_lower = c(tt$conf.int[1], NA_real_),
  ci_upper = c(tt$conf.int[2], NA_real_),
  statistic = c(unname(tt$statistic), unname(wt$statistic)),
  df = c(unname(tt$parameter), NA_real_),
  p_value = c(tt$p.value, wt$p.value),
  method = c(tt$method, wt$method)
) -> overall_tab

cat("\n[01] 总体检验（Mix 的株高变化指数是否偏离 0）：\n")
print(as.data.frame(overall_tab))

# ---- 每个 focal：旧口径 anova_posthoc（Tukey 字母）----
df.idx %>%
  dplyr::group_by(focal) %>%
  bioRtools::anova_posthoc(value = "index", group = "neighbor") %>%
  dplyr::ungroup() -> posthoc_tab

cat("\n[01] 旧口径 anova_posthoc（Tukey 字母）：\n")
print(as.data.frame(posthoc_tab))

# ---- 每个 focal：LMM + emmeans(Sidak) ----
fit_focal <- function(f) {
  d <- df.idx %>%
    dplyr::filter(focal == f) %>%
    droplevels()

  m <- lmerTest::lmer(index ~ neighbor + (1 | time), data = d)

  emmeans::emmeans(m, "neighbor") %>%
    multcomp::cld(alpha = 0.05, Letters = letters, adjust = "sidak", decreasing = TRUE) %>%
    as.data.frame() %>%
    dplyr::select(neighbor, emmean, SE, df, lower.CL, upper.CL, signif = .group) %>%
    dplyr::mutate(focal = f) -> emm

  # 与纯种对照比较（邻体是否显著改变株高）
  emmeans::contrast(
    emmeans::emmeans(m, "neighbor"), "trt.vs.ctrl",
    ref = which(levels(d$neighbor) == as.character(f)),
    adjust = "sidak"
  ) %>%
    as.data.frame() %>%
    dplyr::mutate(focal = f) %>%
    dplyr::select(focal, contrast, estimate, SE, df, t.ratio, p.value) -> ctrl

  list(emmeans = emm, vs_ctrl = ctrl)
}

focal_res <- purrr::map(levels(df.idx$focal), fit_focal)
lmm_emm <- purrr::map_dfr(focal_res, "emmeans")
lmm_ctrl <- purrr::map_dfr(focal_res, "vs_ctrl")

cat("\n[01] LMM + emmeans(Sidak) 各 focal 邻体均值：\n")
print(as.data.frame(lmm_emm))

bioRtools::write_data(
  list(
    index_long = df.idx,
    pure_mean = df.pure,
    overall_test = overall_tab,
    anova_posthoc = posthoc_tab,
    lmm_emmeans = lmm_emm,
    lmm_vs_ctrl = lmm_ctrl
  ),
  file = path_out
)
cat("[01] 输出：", path_out, "\n")

cat("\n[01] sessionInfo()\n")
sessionInfo()
