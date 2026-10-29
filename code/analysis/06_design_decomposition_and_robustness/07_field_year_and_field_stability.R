# ============================================================
# 06_design_decomposition_and_robustness / 07_field_year_and_field_stability.R
# 独立analysis（续 06 主题，Q9）：field处理的"跨年stability"and"田块随机effect"
#   旧版只把 Year 当随机截距、报告一个总体effect；主题 01 的analysis也只给总体 OR。
#   但field实验跨 3 年 × 3 田块，真正关键的问题是：混种effect是否**每年都成立**，还是由某一年驱动？
# 方法（对 Blast / BrownSpot / Smut 三个病害and PH / KGW 两个agronomic_traits各做一遍）：
#   1) 主模型（含交互）：疾病用 glmer(cbind(diseased, 10−diseased) ~ condition * Year + (1|Field), binomial)，
#      性状用 lmer(Value ~ condition * Year + (1|Field))；test condition × Year 交互（effect是否随year变化）；
#   2) 逐年effect：每年单独拟合 condition effect（疾病给出 OR and 95% CI，性状给出均值差and 95% CI）；
#   3) 留一年法（leave-one-year-out）：去掉某一年后的合并effect，判断结论是否依赖单一year；
#   4) 方向一致性：3 年中effect方向相同的year数。
#   主analysis只用 4 个焦点纯系（15/49/69/87），36/74 的固定配对组合在主题 01 中已作补充analysis。
#
# 运行：/opt/homebrew/bin/Rscript "<root>/code/06_design_decomposition_and_robustness/07_field_year_and_field_stability.R"
# 输入：<root>/data/processed_inputs/02_field_disease_survey.xlsx
# 输出：<root>/result/06_design_decomposition_and_robustness/07_表_field_year_and_field_stability.xlsx
#       <root>/result/06_design_decomposition_and_robustness/07_FigureS14_fieldyearstability.pdf / .png
# 其中 <root> = <repository_root>
#
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/06_design_decomposition_and_robustness/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_in <- file.path(project_root, "data", "processed_inputs", "02_field_disease_survey.xlsx")
out_dir <- file.path(project_root, "result", "06_design_decomposition_and_robustness")
path_out <- file.path(out_dir, "07_表_field_year_and_field_stability.xlsx")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stopifnot(file.exists(path_in))

bioRtools::read_data(path_in) -> d2

focus_lines <- c("15", "49", "69", "87")

d2 %>%
  dplyr::select(-1) %>%
  dplyr::mutate(
    FocalPlant = as.character(FocalPlant),
    Year = as.factor(Year),
    Field = as.factor(Field),
    condition = factor(stringr::str_to_title(condition), levels = c("Pure", "Mix"))
  ) %>%
  dplyr::filter(FocalPlant %in% focus_lines, !is.na(Value)) -> df

disease_cats <- c("Blast", "BrownSpot", "Smut")
trait_cats <- c("PH", "KGW")

cat("[07] 输入：", nrow(df), " 行；类别：", paste(sort(unique(df$Category)), collapse = ", "),
    "；year：", paste(sort(unique(as.character(df$Year))), collapse = ", "), "\n")
cat("[07] 每类别每年的sample量：\n")
print(as.data.frame(df %>% dplyr::count(Category, Year) %>%
                      tidyr::pivot_wider(names_from = Year, values_from = n, values_fill = 0)))

# ---------- 主模型：condition x Year 交互 ----------
fit_main <- function(cat) {
  d <- df %>% dplyr::filter(Category == cat) %>% droplevels()

  if (cat %in% disease_cats) {
    d <- d %>% dplyr::mutate(diseased = Value, healthy = 10 - Value)
    m0 <- lme4::glmer(cbind(diseased, healthy) ~ condition + (1 | Field),
                      data = d, family = binomial,
                      control = lme4::glmerControl(optimizer = "bobyqa"))
    m1 <- lme4::glmer(cbind(diseased, healthy) ~ condition * Year + (1 | Field),
                      data = d, family = binomial,
                      control = lme4::glmerControl(optimizer = "bobyqa"))
    cmp <- stats::anova(m0, m1)
    co <- summary(m1)$coefficients
    tibble::tibble(
      category = cat, type = "binomial GLMM",
      n_obs = nrow(d),
      effect_scale = "odds ratio (Mix vs Pure, Year 2019 baseline)",
      effect = exp(co["conditionMix", "Estimate"]),
      ci_lower = exp(co["conditionMix", "Estimate"] - 1.96 * co["conditionMix", "Std. Error"]),
      ci_upper = exp(co["conditionMix", "Estimate"] + 1.96 * co["conditionMix", "Std. Error"]),
      p_effect = co["conditionMix", "Pr(>|z|)"],
      lrt_interaction = cmp$Chisq[2],
      df_interaction = cmp$Df[2],
      p_interaction = cmp$`Pr(>Chisq)`[2],
      singular = lme4::isSingular(m1)
    )
  } else {
    m0 <- lmerTest::lmer(Value ~ condition + (1 | Field), data = d)
    m1 <- lmerTest::lmer(Value ~ condition * Year + (1 | Field), data = d)
    cmp <- stats::anova(m0, m1)
    co <- summary(m1)$coefficients
    tibble::tibble(
      category = cat, type = "LMM",
      n_obs = nrow(d),
      effect_scale = "mean difference (Mix - Pure, Year 2019 baseline)",
      effect = co["conditionMix", "Estimate"],
      ci_lower = co["conditionMix", "Estimate"] - 1.96 * co["conditionMix", "Std. Error"],
      ci_upper = co["conditionMix", "Estimate"] + 1.96 * co["conditionMix", "Std. Error"],
      p_effect = co["conditionMix", "Pr(>|t|)"],
      lrt_interaction = cmp$Chisq[2],
      df_interaction = cmp$Df[2],
      p_interaction = cmp$`Pr(>Chisq)`[2],
      singular = lme4::isSingular(m1)
    )
  }
}

main_tab <- purrr::map_dfr(c(disease_cats, trait_cats), fit_main) %>%
  dplyr::mutate(fdr_interaction = p.adjust(p_interaction, method = "BH"))

cat("\n[07] 主模型（含 condition x Year 交互）：\n")
print(as.data.frame(main_tab %>% dplyr::select(category, type, n_obs, effect, ci_lower, ci_upper,
                                              p_effect, p_interaction, fdr_interaction, singular)),
      row.names = FALSE)

# ---------- 逐年effect ----------
fit_year <- function(cat, yr) {
  d <- df %>% dplyr::filter(Category == cat, Year == yr) %>% droplevels()
  if (dplyr::n_distinct(d$condition) < 2) return(NULL)

  if (cat %in% disease_cats) {
    d <- d %>% dplyr::mutate(diseased = Value, healthy = 10 - Value)
    # 单一year内仍可能有 3 个田块：田块数 > 1 时用 glmer 保留田块随机截距，
    # 否则退回 glm（glmer 会报 "No random effects terms specified in formula"）
    if (dplyr::n_distinct(d$Field) > 1) {
      m <- suppressWarnings(lme4::glmer(cbind(diseased, healthy) ~ condition + (1 | Field),
                                        data = d, family = binomial))
    } else {
      m <- stats::glm(cbind(diseased, healthy) ~ condition, data = d, family = binomial)
    }
    co <- summary(m)$coefficients
    z_col <- if ("Pr(>|z|)" %in% colnames(co)) "Pr(>|z|)" else "Pr(>|t|)"
    tibble::tibble(
      category = cat, year = yr, scale = "OR",
      estimate = exp(co["conditionMix", "Estimate"]),
      ci_lower = exp(co["conditionMix", "Estimate"] - 1.96 * co["conditionMix", "Std. Error"]),
      ci_upper = exp(co["conditionMix", "Estimate"] + 1.96 * co["conditionMix", "Std. Error"]),
      p_value = co["conditionMix", z_col], n_obs = nrow(d)
    )
  } else {
    m <- stats::lm(Value ~ condition, data = d)
    co <- summary(m)$coefficients
    tibble::tibble(
      category = cat, year = yr, scale = "mean difference",
      estimate = co["conditionMix", "Estimate"],
      ci_lower = co["conditionMix", "Estimate"] - 1.96 * co["conditionMix", "Std. Error"],
      ci_upper = co["conditionMix", "Estimate"] + 1.96 * co["conditionMix", "Std. Error"],
      p_value = co["conditionMix", "Pr(>|t|)"], n_obs = nrow(d)
    )
  }
}

year_grid <- tidyr::expand_grid(category = c(disease_cats, trait_cats),
                                year = sort(unique(df$Year)))
year_tab <- purrr::pmap_dfr(list(year_grid$category, year_grid$year), fit_year)

# 方向一致性（以 OR < 1（疾病降低）或均值差 > 0 为该类别的"预期方向"由总体effect决定）
overall_dir <- main_tab %>% dplyr::select(category, effect) %>% dplyr::rename(overall_effect = effect)

year_dir <- year_tab %>%
  dplyr::left_join(overall_dir, by = "category") %>%
  dplyr::mutate(same_direction = sign(estimate - ifelse(scale == "OR", 1, 0)) ==
                  sign(overall_effect - ifelse(scale == "OR", 1, 0))) %>%
  dplyr::group_by(category) %>%
  dplyr::summarise(n_year = dplyr::n(), n_same_direction = sum(same_direction),
                   n_significant = sum(p_value < 0.05), .groups = "drop")

# ---------- 留一年法 ----------
fit_loo <- function(cat, drop_year) {
  d <- df %>% dplyr::filter(Category == cat, Year != drop_year) %>% droplevels()

  if (cat %in% disease_cats) {
    d <- d %>% dplyr::mutate(diseased = Value, healthy = 10 - Value)
    m <- suppressWarnings(lme4::glmer(cbind(diseased, healthy) ~ condition + (1 | Field),
                                      data = d, family = binomial))
    co <- summary(m)$coefficients
    tibble::tibble(category = cat, dropped_year = drop_year, scale = "OR",
                   estimate = exp(co["conditionMix", "Estimate"]),
                   ci_lower = exp(co["conditionMix", "Estimate"] - 1.96 * co["conditionMix", "Std. Error"]),
                   ci_upper = exp(co["conditionMix", "Estimate"] + 1.96 * co["conditionMix", "Std. Error"]),
                   p_value = co["conditionMix", "Pr(>|z|)"])
  } else {
    m <- lmerTest::lmer(Value ~ condition + (1 | Field), data = d)
    co <- summary(m)$coefficients
    tibble::tibble(category = cat, dropped_year = drop_year, scale = "mean difference",
                   estimate = co["conditionMix", "Estimate"],
                   ci_lower = co["conditionMix", "Estimate"] - 1.96 * co["conditionMix", "Std. Error"],
                   ci_upper = co["conditionMix", "Estimate"] + 1.96 * co["conditionMix", "Std. Error"],
                   p_value = co["conditionMix", "Pr(>|t|)"])
  }
}

loo_grid <- tidyr::expand_grid(category = c(disease_cats, trait_cats), drop_year = levels(df$Year))
loo_tab <- purrr::pmap_dfr(list(loo_grid$category, loo_grid$drop_year), fit_loo)

cat("\n[07] 逐年effect：\n")
print(as.data.frame(year_tab), row.names = FALSE)
cat("\n[07] 方向一致性：\n")
print(as.data.frame(year_dir), row.names = FALSE)
cat("\n[07] 留一年法：\n")
print(as.data.frame(loo_tab), row.names = FALSE)

bioRtools::write_data(
  list(main_model = main_tab, by_year = year_tab, direction_consistency = year_dir,
       leave_one_year_out = loo_tab),
  file = path_out
)
cat("[07] 输出：", path_out, "\n")

# ---------- 图：逐年effect森林图 ----------
ref_lines <- tibble::tibble(
  category = factor(c(disease_cats, trait_cats), levels = c(disease_cats, trait_cats)),
  ref = c(rep(1, length(disease_cats)), rep(0, length(trait_cats)))
)

p.forest <- year_tab %>%
  dplyr::mutate(
    category = factor(category, levels = c(disease_cats, trait_cats)),
    year = factor(year)
  ) %>%
  ggplot2::ggplot(aes(x = estimate, y = year, color = year)) +
  ggplot2::geom_vline(data = ref_lines, aes(xintercept = ref),
                      linetype = "dashed", color = "grey40") +
  ggplot2::geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::facet_wrap(~category, scales = "free_x") +
  ggplot2::scale_y_discrete(limits = rev(levels(factor(year_tab$year)))) +
  bioRtools::scale_color_research() +
  ggplot2::labs(
    x = "Effect estimate per year (diseases: odds ratio; traits: mean difference)",
    y = "Year", color = "Year",
    caption = paste0("Dash line at 1 for odds ratios (no effect). Interaction P values: ",
                     paste0(main_tab$category, " ", formatC(main_tab$p_interaction, format = "f", digits = 3),
                            collapse = "; "))
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "none",
    plot.caption = ggplot2::element_text(size = 7, hjust = 0)
  ) -> p.forest

bioRtools::write_data(p.forest, file = file.path(out_dir, "07_FigureS14_fieldyearstability.pdf"),
                      width = 12, height = 6, units = "in")
bioRtools::write_data(p.forest, file = file.path(out_dir, "07_FigureS14_fieldyearstability.png"),
                      width = 12, height = 6, units = "in", dpi = 500)
cat("[07] 输出：07_FigureS14_fieldyearstability.pdf / .png\n")

cat("\n[07] sessionInfo()\n")
sessionInfo()
