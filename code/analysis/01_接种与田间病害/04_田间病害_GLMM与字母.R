# ============================================================
# 01_接种与田间病害 / 04_田间病害_GLMM与字母.R
# 作用：元阳梯田田间病害调查（3 年 × 3 田块 × 6 纯系组合）的二项 GLMM 分析
#       —— 田间病害调查小节（正文引用）+ Figure S2G / S2H 的统计输入
# 运行：/opt/homebrew/bin/Rscript "<root>/code/01_接种与田间病害/04_田间病害_GLMM与字母.R"
#
# 输入：<root>/data/01_分析输入表/02_元阳梯田田间病害调查表.xlsx（5658 行 × 10 列）
# 输出：<root>/result/01_接种与田间病害/04_表_田间病害主要GLMM.xlsx（多 sheet）
#       <root>/result/01_接种与田间病害/05_表_田间病害补充_Smut与36_74.xlsx（多 sheet）
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 模型：cbind(diseased_count, healthy_count) ~ combination + (1 | Year)，family = binomial
#       combination = FocalPlant-Neighbor；Year 为随机截距（2019–2021 调查年份）
#       口径：每个小区调查 10 株，diseased_count = Value，healthy_count = 10 − Value
#       （实测 Blast 最大 8、BrownSpot 最大 10、Smut 最大 7，与 10 株口径一致）
# 范围说明（按用户 2026-09-22 确认的原则）：
#       主分析 = 焦点纯系 15 / 49 / 69 / 87 × Blast / BrownSpot（与旧版一致，保证可比）；
#       补充分析 = Smut（4 个焦点纯系）+ 36 / 74 两个额外纯系的 Blast / BrownSpot / Smut。
# 多重比较：字母用 multcomp::cld(adjust = "sidak")；Prob 尺度与 link 尺度单调对应，
#       故字母（link 尺度）可直接对应到 Prob 估计值。GLMM 采用渐近（z）推断，无小 df 问题。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：单步耗时 < 60 s，每次运行都重新生成。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录见 code/01_接种与田间病害/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_in <- file.path(project_root, "data", "01_分析输入表", "02_元阳梯田田间病害调查表.xlsx")
out_dir <- file.path(project_root, "result", "01_接种与田间病害")
path_main <- file.path(out_dir, "04_表_田间病害主要GLMM.xlsx")
path_supp <- file.path(out_dir, "05_表_田间病害补充_Smut与36_74.xlsx")

stopifnot(file.exists(path_in))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("[04] 输入：", path_in, "\n")

bioRtools::read_data(path_in) -> df.raw

stopifnot(nrow(df.raw) == 5658)
need_cols <- c("FocalPlant", "Neighbor", "condition", "Field", "Year", "Category", "Value")
stopifnot(all(need_cols %in% colnames(df.raw)))

cat("[04] Category 计数：\n")
print(table(df.raw$Category))
cat("[04] FocalPlant 计数：\n")
print(table(df.raw$FocalPlant))

df.raw %>%
  dplyr::select(-1) %>%
  dplyr::mutate(
    FocalPlant = as.character(FocalPlant),
    Neighbor = as.character(Neighbor),
    Year = as.factor(Year),
    condition = factor(stringr::str_to_title(condition), levels = c("Pure", "Mix")),
    combination = paste0(FocalPlant, "-", Neighbor),
    disease = dplyr::recode(Category, "Blast" = "Blast", "BrownSpot" = "BrownSpot", "Smut" = "Smut")
  ) -> df

# 三种病害按 10 株口径展开成二项计数
df %>%
  dplyr::filter(Category %in% c("Blast", "BrownSpot", "Smut")) %>%
  dplyr::mutate(
    diseased_count = Value,
    healthy_count = 10 - Value
  ) -> df.dis

stopifnot(all(df.dis$healthy_count >= 0), !anyNA(df.dis$diseased_count))
stopifnot(nrow(df.dis) == 954)

cat("[04] 病害数据：", nrow(df.dis), "行；focal x disease x combination 组合数：",
    nrow(dplyr::count(df.dis, FocalPlant, Category, combination)), "\n")

main_focals <- c("15", "49", "69", "87")

# ---- 单个 focal × disease 的 GLMM ----
fit_one <- function(f, ds) {
  d <- df.dis %>%
    dplyr::filter(FocalPlant == f, Category == ds) %>%
    dplyr::mutate(combination = factor(combination)) %>%
    droplevels()

  m <- lme4::glmer(
    cbind(diseased_count, healthy_count) ~ combination + (1 | Year),
    data = d, family = binomial,
    control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )

  # 字母在 link 尺度计算（Prob 与 logit 单调对应，可直接对应到 Prob 估计值）
  emmeans::emmeans(m, "combination") %>%
    multcomp::cld(alpha = 0.05, Letters = letters, adjust = "sidak", decreasing = TRUE) %>%
    as.data.frame() %>%
    dplyr::select(combination, signif = .group) -> letters_tab

  emmeans::emmeans(m, "combination", type = "response") %>%
    as.data.frame() %>%
    dplyr::select(combination, prob, SE, asymp.LCL, asymp.UCL) %>%
    dplyr::left_join(letters_tab, by = "combination") %>%
    dplyr::mutate(
      focal = f,
      disease = ds,
      condition = dplyr::if_else(combination == paste0(f, "-", f), "Pure", "Mix"),
      neighbor = stringr::str_remove(combination, "^[0-9]+-"),
      n_obs = nrow(d),
      conv_code = m@optinfo$conv$opt,
      singular = lme4::isSingular(m)
    ) %>%
    dplyr::select(disease, focal, combination, neighbor, condition, prob, SE,
                  asymp.LCL, asymp.UCL, signif, n_obs, conv_code, singular)
}

grid_task <- tidyr::expand_grid(focal = unique(df.dis$FocalPlant), disease = unique(df.dis$Category)) %>%
  dplyr::mutate(in_main = focal %in% main_focals & disease %in% c("Blast", "BrownSpot"))

cat("[04] 待拟合 GLMM 数：", nrow(grid_task), "\n")

res <- purrr::pmap(list(grid_task$focal, grid_task$disease), fit_one)
all_tab <- purrr::list_rbind(res)

cat("[04] 收敛情况（conv_code 非 0 或 singular = TRUE 需注意）：\n")
print(as.data.frame(dplyr::count(all_tab, disease, focal, conv_code, singular)))

all_tab %>%
  dplyr::filter(focal %in% main_focals, disease %in% c("Blast", "BrownSpot")) -> main_tab

# ---- 总体 Pure vs Mix 的固定效应（主分析 4 个焦点纯系）----
fit_overall <- function(ds) {
  d <- df.dis %>%
    dplyr::filter(FocalPlant %in% main_focals, Category == ds,
                  condition %in% c("Pure", "Mix")) %>%
    dplyr::mutate(condition = factor(condition, levels = c("Pure", "Mix"))) %>%
    droplevels()

  m <- lme4::glmer(
    cbind(diseased_count, healthy_count) ~ condition + FocalPlant + (1 | Year),
    data = d, family = binomial,
    control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )

  summary(m)$coefficients %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "term") %>%
    dplyr::rename(std_error = `Std. Error`, z_value = `z value`, p_value = `Pr(>|z|)`) %>%
    dplyr::filter(term == "conditionMix") %>%
    dplyr::mutate(
      disease = ds,
      n_obs = nrow(d),
      odds_ratio = exp(Estimate),
      or_lower = exp(Estimate - 1.96 * std_error),
      or_upper = exp(Estimate + 1.96 * std_error)
    ) %>%
    dplyr::select(disease, term, estimate = Estimate, std_error, z_value, p_value,
                  odds_ratio, or_lower, or_upper, n_obs)
}

overall_tab <- c("Blast", "BrownSpot") %>% purrr::map(fit_overall) %>% purrr::list_rbind()

cat("\n[04] 总体 Pure vs Mix（田间，4 个焦点纯系）：\n")
print(as.data.frame(overall_tab))

# ---- 描述性统计 ----
desc_tab <- df.dis %>%
  dplyr::filter(FocalPlant %in% main_focals, Category %in% c("Blast", "BrownSpot")) %>%
  dplyr::mutate(incidence = Value / 10) %>%
  dplyr::group_by(Category, condition) %>%
  dplyr::summarise(
    n_plot = dplyr::n(),
    mean_incidence = mean(incidence),
    sd_incidence = sd(incidence),
    .groups = "drop"
  ) %>%
  dplyr::rename(disease = Category)

cat("\n[04] 描述性统计（发病率 = Value / 10）：\n")
print(as.data.frame(desc_tab))

cat("\n[04] 主分析 emmeans（Prob 尺度 + 字母）：\n")
print(as.data.frame(main_tab), row.names = FALSE)

# ---- 输出 ----
bioRtools::write_data(
  list(
    emmean_letters_main = main_tab,
    overall_condition = overall_tab,
    descriptive = desc_tab
  ),
  file = path_main
)
cat("[04] 输出：", path_main, "\n")

bioRtools::write_data(
  list(emmean_letters_all = all_tab),
  file = path_supp
)
cat("[04] 输出：", path_supp, "\n")

cat("\n[04] sessionInfo()\n")
sessionInfo()
