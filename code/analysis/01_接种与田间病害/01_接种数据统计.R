# ============================================================
# 01_接种与田间病害 / 01_接种数据统计.R
# 作用：温室接种数据（稻瘟病 blast / 胡麻叶斑病 brown spot）的线性混合模型总体检验
#       与 Z-score 标准化中间表（Figure S2A / S2B 的统计与绘图输入）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/01_接种与田间病害/01_接种数据统计.R"
#
# 输入：<root>/data/01_分析输入表/01_接种结果_合并表.xlsx（1659 行 × 13 列）
# 输出：<root>/result/01_接种与田间病害/01_表_接种LMM总体检验.xlsx（多 sheet）
#       <root>/result/01_接种与田间病害/02_表_接种Zscore中间表.xlsx
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 模型：lesion_number_per_cm2 ~ group + (1 | batch)
#       group = Pure / Mix（固定效应）；batch = 实验批次（随机截距）
# 统计口径说明（重要）：
#   1) 组间检验取固定效应 contrast（groupMix vs Pure），df 用 Satterthwaite（lmerTest）；
#   2) 各组的 emmean 因 batch 只有 2–3 个水平、批次方差远大于残差，其 df ≈ 1–2、CI 极宽，
#      不能用于组间比较（本表仅作记录），比较一律看 contrast（contrast 的 df 为 613 / 1039）；
#   3) 旧版脚本用 lme4::lmer() + summary()，不输出 P 值，旧版"总体检验"本无 P 值可比，
#      本次以 lmerTest Satterthwaite 结果为准；
#   4) batch 与 combination 未完全混杂（每个 combination 出现在 2–3 个批次），已核实。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：单步耗时 < 10 s，每次运行都重新生成，不做"输出已存在则跳过"判断。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录写在 code/01_接种与田间病害/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

# 本脚本不含随机过程；固定种子以便将来加入重抽样时行为一致
set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_in <- file.path(project_root, "data", "01_分析输入表", "01_接种结果_合并表.xlsx")
out_dir <- file.path(project_root, "result", "01_接种与田间病害")
path_out_lmm <- file.path(out_dir, "01_表_接种LMM总体检验.xlsx")
path_out_z <- file.path(out_dir, "02_表_接种Zscore中间表.xlsx")

stopifnot(file.exists(path_in))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("[01] 输入：", path_in, "\n")
cat("[01] 输出目录：", out_dir, "\n")

bioRtools::read_data(path_in) -> df.raw

# ---- 输入校验 ----
stopifnot(nrow(df.raw) == 1659, ncol(df.raw) == 13)
need_cols <- c("focal", "neighbor", "disease", "batch", "group", "lesion_number_per_cm2")
stopifnot(all(need_cols %in% colnames(df.raw)))
stopifnot(!anyNA(df.raw$lesion_number_per_cm2))

cat("[01] 输入表：", nrow(df.raw), "行 x", ncol(df.raw), "列\n")
cat("[01] disease x group 计数：\n")
print(table(df.raw$disease, df.raw$group))

df.raw %>%
  dplyr::mutate(
    batch = as.factor(batch),
    group = factor(group, levels = c("Pure", "Mix")),
    disease = factor(disease, levels = c("blast", "brown spot"))
  ) -> df

diseases <- levels(df$disease)

# ---- 逐病害拟合 LMM 并汇总 ----
fit_by_disease <- function(ds) {
  d <- df %>%
    dplyr::filter(disease == ds) %>%
    droplevels()

  m <- lmerTest::lmer(lesion_number_per_cm2 ~ group + (1 | batch), data = d)

  # 固定效应（Satterthwaite df，含组间对比的 P 值）
  summary(m)$coefficients %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "term") %>%
    dplyr::rename(std_error = `Std. Error`, t_value = `t value`, p_value = `Pr(>|t|)`) %>%
    dplyr::select(term, estimate = Estimate, std_error, df, t_value, p_value) %>%
    dplyr::mutate(disease = ds, n_obs = nrow(d), n_batch = dplyr::n_distinct(d$batch)) %>%
    dplyr::select(disease, term, estimate, std_error, df, t_value, p_value, n_obs, n_batch) -> fixed

  # Type III ANOVA
  as.data.frame(anova(m)) %>%
    tibble::rownames_to_column(var = "term") %>%
    dplyr::rename(sum_sq = `Sum Sq`, mean_sq = `Mean Sq`, num_df = NumDF,
                  den_df = DenDF, f_value = `F value`, p_value = `Pr(>F)`) %>%
    dplyr::select(term, num_df, den_df, f_value, p_value) %>%
    dplyr::mutate(disease = ds) %>%
    dplyr::select(disease, term, num_df, den_df, f_value, p_value) -> aov_tab

  # 随机效应方差与 ICC
  vc <- lme4::VarCorr(m)
  var_batch <- as.numeric(vc$batch)
  var_resid <- attr(vc, "sc")^2
  tibble::tibble(
    disease = ds,
    n_obs = nrow(d),
    n_batch = dplyr::n_distinct(d$batch),
    var_batch = var_batch,
    var_resid = var_resid,
    icc_batch = var_batch / (var_batch + var_resid)
  ) -> info_tab

  # 各组 emmean（仅作记录：df 极小、CI 极宽，不用于组间比较）
  emmeans::emmeans(m, "group", lmer.df = "satterthwaite") %>%
    as.data.frame() %>%
    dplyr::mutate(disease = ds) %>%
    dplyr::select(disease, group, emmean, SE, df, lower.CL, upper.CL) -> emmean_tab

  # 描述性统计
  d %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      n_leaf = dplyr::n(),
      mean = mean(lesion_number_per_cm2),
      sd = sd(lesion_number_per_cm2),
      median = median(lesion_number_per_cm2),
      q1 = quantile(lesion_number_per_cm2, 0.25),
      q3 = quantile(lesion_number_per_cm2, 0.75),
      min = min(lesion_number_per_cm2),
      max = max(lesion_number_per_cm2),
      .groups = "drop"
    ) %>%
    dplyr::mutate(disease = ds) %>%
    dplyr::select(disease, group, n_leaf, mean, sd, median, q1, q3, min, max) -> desc_tab

  list(fixed = fixed, anova = aov_tab, info = info_tab, emmean = emmean_tab, desc = desc_tab)
}

res <- purrr::map(diseases, fit_by_disease)
names(res) <- diseases

fixed_all <- purrr::map_dfr(res, "fixed")
aov_all <- purrr::map_dfr(res, "anova")
info_all <- purrr::map_dfr(res, "info")
emmean_all <- purrr::map_dfr(res, "emmean")
desc_all <- purrr::map_dfr(res, "desc")

cat("\n[01] 总体 Pure vs Mix 固定效应（Satterthwaite）：\n")
fixed_all %>%
  dplyr::filter(term == "groupMix") %>%
  dplyr::select(disease, estimate, std_error, df, t_value, p_value) %>%
  as.data.frame() %>%
  print()

cat("\n[01] Type III ANOVA：\n")
print(as.data.frame(aov_all))

cat("\n[01] 随机效应方差与 ICC：\n")
print(as.data.frame(info_all))

cat("\n[01] 描述性统计：\n")
print(as.data.frame(desc_all))

# ---- Z-score 标准化中间表（按 disease × batch 标准化，与旧版一致）----
df %>%
  dplyr::group_by(disease, batch) %>%
  dplyr::mutate(lesion_norm = as.vector(scale(lesion_number_per_cm2))) %>%
  dplyr::ungroup() -> df.z

stopifnot(!anyNA(df.z$lesion_norm))

cat("\n[01] Z-score 表：", nrow(df.z), "行 x", ncol(df.z), "列；lesion_norm 范围 ",
    min(df.z$lesion_norm), " ~ ", max(df.z$lesion_norm), "\n")

# ---- 输出 ----
bioRtools::write_data(
  list(
    fixed_effect = fixed_all,
    anova_type3 = aov_all,
    model_info = info_all,
    emmean_group_mean = emmean_all,
    descriptive = desc_all
  ),
  file = path_out_lmm
)
cat("[01] 输出：", path_out_lmm, "\n")

bioRtools::write_data(df.z, file = path_out_z)
cat("[01] 输出：", path_out_z, "\n")

cat("\n[01] sessionInfo()\n")
sessionInfo()
