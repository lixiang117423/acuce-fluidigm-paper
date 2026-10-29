# ============================================================
# 02_agronomic_traits / 06_field_height_and_thousand_grain_weight_supplement.R
# 作用：field实验（2019–2021、3 田块、6 纯系）的plant_height（PH）andthousand_grain_weight（KGW）补充analysis
#       —— 作为主题 02 agronomic_traits的独立year/独立性状验证（Figure S4）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/02_agronomic_traits/06_field_height_and_thousand_grain_weight_supplement.R"
#
# 输入：<root>/data/processed_inputs/02_field_disease_survey.xlsx
#       （Category 中 PH 3564 行、KGW 1140 行，此前旧版未analysis）
# 输出：<root>/result/02_agronomic_traits/12_表_fieldplant_heightandthousand_grain_weight.xlsx
#       <root>/result/02_agronomic_traits/12_FigureS4_fieldplant_heightandthousand_grain_weight.pdf / .png
#
# 模型：value ~ condition + FocalPlant + (1 | Year)，分 PH / KGW 各拟合一次（condition = Pure / Mix）；
#       同时给出每个 focal 的 condition effect（emmeans 估计值 + 95% 渐近 CI），用于判断混种是否改变fieldplant_heightandthousand_grain_weight。
# 口径notes：36 / 74 两个纯系在field仅and固定伙伴配对（36-49、69-74），故其 condition 对比同样可比；
#       本图按 6 个纯系全部展示，主结论以 4 个焦点纯系（15/49/69/87）为准。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：单步耗时 < 60 s，每次运行都重新generate。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/02_agronomic_traits/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_in <- file.path(project_root, "data", "processed_inputs", "02_field_disease_survey.xlsx")
out_dir <- file.path(project_root, "result", "02_agronomic_traits")
path_out <- file.path(out_dir, "12_表_fieldplant_heightandthousand_grain_weight.xlsx")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stopifnot(file.exists(path_in))

format_p <- function(p, digits = 3) {
  dplyr::if_else(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, format = "f", digits = digits)))
}

# emmeans 在观测数 > 3000 时会关闭自由度计算，此时置信区间列名为 asymp.LCL / asymp.UCL，
# 统一改成 lower.CL / upper.CL，便于下游plottingand选列
norm_ci <- function(x) {
  if (!"lower.CL" %in% names(x) && "asymp.LCL" %in% names(x)) {
    x <- dplyr::rename(x, lower.CL = asymp.LCL, upper.CL = asymp.UCL)
  }
  x
}

bioRtools::read_data(path_in) -> df.raw

df.raw %>%
  dplyr::filter(Category %in% c("PH", "KGW")) %>%
  dplyr::mutate(
    FocalPlant = as.character(FocalPlant),
    Neighbor = as.character(Neighbor),
    Year = as.factor(Year),
    Field = as.factor(Field),
    condition = factor(stringr::str_to_title(condition), levels = c("Pure", "Mix"))
  ) -> df

cat("[06] 输入：PH", sum(df$Category == "PH"), "行；KGW", sum(df$Category == "KGW"), "行\n")
cat("[06] PH/KGW data覆盖（纯系 x year）：\n")
print(table(df$Category, df$FocalPlant, df$Year))

# ---- 总体 condition effect ----
fit_overall <- function(ct) {
  d <- df %>% dplyr::filter(Category == ct, !is.na(Value)) %>% droplevels()
  m <- lmerTest::lmer(Value ~ condition + FocalPlant + (1 | Year), data = d)
  s <- summary(m)$coefficients
  tibble::tibble(
    category = ct,
    n_obs = nrow(d),
    estimate = s["conditionMix", "Estimate"],
    std_error = s["conditionMix", "Std. Error"],
    df = s["conditionMix", "df"],
    t_value = s["conditionMix", "t value"],
    p_value = s["conditionMix", "Pr(>|t|)"],
    mean_pure = mean(d$Value[d$condition == "Pure"]),
    mean_mix = mean(d$Value[d$condition == "Mix"])
  )
}

overall_tab <- purrr::map_dfr(c("PH", "KGW"), fit_overall)

cat("\n[06] 总体 condition effect（LMM）：\n")
print(as.data.frame(overall_tab))

# ---- 每个纯系的 condition 均值and CI ----
emm_tab <- purrr::map_dfr(c("PH", "KGW"), function(ct) {
  d <- df %>% dplyr::filter(Category == ct, !is.na(Value)) %>% droplevels()
  m <- lmerTest::lmer(Value ~ condition + FocalPlant + (1 | Year), data = d)
  emmeans::emmeans(m, ~ condition | FocalPlant, lmer.df = "satterthwaite") %>%
    as.data.frame() %>%
    norm_ci() %>%
    dplyr::mutate(category = ct) %>%
    dplyr::select(category, focal = FocalPlant, condition, emmean, SE, df, lower.CL, upper.CL)
})

cat("\n[06] 各纯系 condition 均值：\n")
print(as.data.frame(emm_tab))

# ---- 描述性statistics ----
df %>%
  dplyr::filter(!is.na(Value)) %>%
  dplyr::group_by(Category, FocalPlant, condition, Year) %>%
  dplyr::summarise(n = dplyr::n(), mean = mean(Value), sd = sd(Value), .groups = "drop") -> desc_tab

bioRtools::write_data(
  list(overall_model = overall_tab, emmeans = emm_tab, descriptive = desc_tab),
  file = path_out
)
cat("[06] 输出：", path_out, "\n")

# ---- Figure S4 ----
emm_tab %>%
  dplyr::mutate(
    category = factor(category, levels = c("PH", "KGW"),
                      labels = c("Plant height (cm)", "1000-grain weight (g)")),
    focal = factor(focal, levels = c("15", "36", "49", "69", "74", "87"))
  ) %>%
  ggplot2::ggplot(aes(x = focal, y = emmean, color = condition)) +
  ggplot2::geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.18,
                         position = ggplot2::position_dodge(width = 0.5)) +
  ggplot2::geom_point(size = 2.6, position = ggplot2::position_dodge(width = 0.5)) +
  ggplot2::geom_text(
    aes(y = upper.CL, label = focal), alpha = 0, show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~category, scales = "free_y") +
  bioRtools::scale_color_research() +
  ggplot2::labs(
    x = "Focal accession",
    y = "Estimated mean (95% CI)",
    color = "Planting mode",
    caption = paste0(
      "Field trial 2019-2021 (3 years x 3 fields). Overall condition effect: ",
      paste0(overall_tab$category, " ", format_p(overall_tab$p_value), collapse = "; ")
    )
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    legend.position = "bottom",
    strip.text = ggplot2::element_text(face = "bold"),
    plot.caption = ggplot2::element_text(size = 8, hjust = 0)
  ) -> p.s4

bioRtools::write_data(p.s4, file = file.path(out_dir, "12_FigureS4_fieldplant_heightandthousand_grain_weight.pdf"),
                      width = 10, height = 5, units = "in")
bioRtools::write_data(p.s4, file = file.path(out_dir, "12_FigureS4_fieldplant_heightandthousand_grain_weight.png"),
                      width = 10, height = 5, units = "in", dpi = 500)
cat("[06] 输出：12_FigureS4_fieldplant_heightandthousand_grain_weight.pdf / .png\n")

cat("\n[06] sessionInfo()\n")
sessionInfo()
