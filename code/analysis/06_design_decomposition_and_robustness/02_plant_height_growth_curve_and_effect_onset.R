# ============================================================
# 06_design_decomposition_and_robustness / 02_plant_height_growth_curve_and_effect_onset.R
# 独立analysis（不属于旧版流程）：把plant_height当纵向性状处理，回答两个旧版没问的问题
#   Q1 互作effect改变的是生长速率还是final水平（group x day 交互）
#   Q2 effect从第几次调查开始出现（逐time点test + FDR 控制，给出 onset time）
# 旧版把 15 个time点压成一个变化指数再做一次 ANOVA，纵向信息全部丢弃，既无法区分速率effectand水平effect，
# 也给不出effect出现的time窗口。
#
# 运行：/opt/homebrew/bin/Rscript "<root>/code/06_design_decomposition_and_robustness/02_plant_height_growth_curve_and_effect_onset.R"
# 输入：<root>/data/processed_inputs/03_agronomic_traits_height_tiller_yield.xlsx
# 输出：<root>/result/06_design_decomposition_and_robustness/02_表_plant_height_growth_curveandonset点.xlsx
#       <root>/result/06_design_decomposition_and_robustness/02_FigureS9_plant_heightgrowth_curve_and_effect_onset.pdf / .png
# 其中 <root> = <repository_root>
#
# 口径：只用有 time 的plant_heightrecord（focal 36/74 无 time，自动排除）；time 为 YYYYMMDD 整数，
#   换算为距首次调查（20220711）的天数作为连续变量；Pure 组 neighbor = focal。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/06_design_decomposition_and_robustness/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_in <- file.path(project_root, "data", "processed_inputs", "03_agronomic_traits_height_tiller_yield.xlsx")
out_dir <- file.path(project_root, "result", "06_design_decomposition_and_robustness")
path_out <- file.path(out_dir, "02_表_plant_height_growth_curveandonset点.xlsx")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stopifnot(file.exists(path_in))

bioRtools::read_data(path_in) -> d3

d3 %>%
  dplyr::filter(phenotype == "plant height", !is.na(time)) %>%
  dplyr::mutate(
    focal = as.factor(focal),
    neighbor = as.factor(neighbor),
    group = factor(group, levels = c("Pure", "Mix")),
    date = as.Date(as.character(time), format = "%Y%m%d"),
    day = as.numeric(date - min(date)),
    combo = paste0(focal, "-", neighbor)
  ) -> df

cat("[02] plant_heightdata：", nrow(df), " 行；time点 ", dplyr::n_distinct(df$day),
    " 个（", min(df$day), "-", max(df$day), " 天）；组合 ", dplyr::n_distinct(df$combo), " 个\n")

# ---------- Q1：速率 vs 水平 ----------
rate_tab <- local({
  # 随机斜率（1 + day | combo）在 16 个组合下不可识别（large eigenvalue / max|grad 超限），
  # 故随机effect只保留组合截距；速率effect由 group:day 固定effect承担。
  # 注意：系数表里有一列名为 df，若在 mutate() 里直接用 df$combo 会被data掩码解析成该列，
  # 因此这里先把sample量and组合数取成普通标量。
  n_obs_val <- nrow(df)
  n_combo_val <- dplyr::n_distinct(df$combo)
  m <- lmerTest::lmer(value ~ group * day + (1 | combo), data = df,
                      control = lme4::lmerControl(optimizer = "bobyqa"))
  co <- summary(m)$coefficients
  as.data.frame(co) %>%
    tibble::rownames_to_column(var = "term") %>%
    dplyr::rename(std_error = `Std. Error`, t_value = `t value`, p_value = `Pr(>|t|)`) %>%
    dplyr::select(term, estimate = Estimate, std_error, df, t_value, p_value) %>%
    dplyr::mutate(scope = "all focals", n_obs = n_obs_val, n_combo = n_combo_val)
})

cat("\n[02] Q1 速率 vs 水平（group x day，随机截距+随机斜率 | combo）：\n")
print(as.data.frame(rate_tab), row.names = FALSE)

rate_by_focal <- purrr::map_dfr(levels(df$focal), function(f) {
  d <- df %>% dplyr::filter(focal == f) %>% droplevels()
  m <- lmerTest::lmer(value ~ group * day + (1 | combo), data = d,
                      control = lme4::lmerControl(optimizer = "bobyqa"))
  co <- summary(m)$coefficients
  tibble::tibble(
    focal = f,
    effect_groupMix = co["groupMix", "Estimate"],
    p_groupMix = co["groupMix", "Pr(>|t|)"],
    effect_groupMix_day = co["groupMix:day", "Estimate"],
    p_groupMix_day = co["groupMix:day", "Pr(>|t|)"],
    n_obs = nrow(d)
  )
})

cat("\n[02] 各 focal 的速率effect（groupMix:day = 每日plant_height差的变化率）：\n")
print(as.data.frame(rate_by_focal), row.names = FALSE)

# ---------- Q2：逐time点test + onset ----------
time_tab <- df %>%
  dplyr::group_by(day, date) %>%
  dplyr::group_modify(~ {
    tt <- stats::t.test(value ~ group, data = .x)
    tibble::tibble(
      n_pure = sum(.x$group == "Pure"),
      n_mix = sum(.x$group == "Mix"),
      mean_pure = mean(.x$value[.x$group == "Pure"]),
      mean_mix = mean(.x$value[.x$group == "Mix"]),
      diff_mix_pure = mean(.x$value[.x$group == "Mix"]) - mean(.x$value[.x$group == "Pure"]),
      t_p_value = tt$p.value
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(t_fdr = p.adjust(t_p_value, method = "BH")) %>%
  dplyr::arrange(day)

onset_raw <- time_tab %>% dplyr::filter(t_p_value < 0.05) %>% dplyr::slice_head(n = 1)
onset_fdr <- time_tab %>% dplyr::filter(t_fdr < 0.05) %>% dplyr::slice_head(n = 1)

cat("\n[02] Q2 逐time点test（前 6 行）：\n")
print(as.data.frame(time_tab %>% dplyr::slice_head(n = 6)), row.names = FALSE)
cat("首次原始 P < 0.05：", ifelse(nrow(onset_raw) > 0, paste0("第 ", onset_raw$day, " 天（", onset_raw$date, "）"), "无"),
    "；首次 FDR < 0.05：", ifelse(nrow(onset_fdr) > 0, paste0("第 ", onset_fdr$day, " 天（", onset_fdr$date, "）"), "无"), "\n")
cat("FDR < 0.05 的time点数：", sum(time_tab$t_fdr < 0.05), " / ", nrow(time_tab), "\n")
cat("原始 P < 0.05 的time点数：", sum(time_tab$t_p_value < 0.05), " / ", nrow(time_tab), "\n")

bioRtools::write_data(
  list(rate_model_all = rate_tab, rate_model_by_focal = rate_by_focal, time_point_tests = time_tab),
  file = path_out
)
cat("[02] 输出：", path_out, "\n")

# ---------- 图 ----------
trace_tab <- df %>%
  dplyr::group_by(focal, neighbor, group, day) %>%
  dplyr::summarise(mean_value = mean(value), se = sd(value) / sqrt(dplyr::n()), .groups = "drop")

p.trace <- trace_tab %>%
  ggplot2::ggplot(aes(x = day, y = mean_value, color = neighbor, linetype = group)) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::facet_wrap(~focal, ncol = 4) +
  bioRtools::scale_color_research() +
  ggplot2::labs(
    x = "Days after first survey",
    y = "Plant height (cm)",
    color = "Neighbor",
    linetype = "Planting mode"
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "right"
  ) -> p.trace

p.pval <- time_tab %>%
  ggplot2::ggplot(aes(x = day, y = t_p_value)) +
  ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
  ggplot2::geom_line(color = "#1e4e8c", linewidth = 0.7) +
  ggplot2::geom_point(aes(color = t_fdr < 0.05), size = 2.2) +
  ggplot2::scale_color_manual(
    values = c(`TRUE` = "#a63446", `FALSE` = "grey55"),
    labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "n.s."), name = NULL
  ) +
  ggplot2::labs(x = "Days after first survey", y = "P value (Mix vs Pure, all focals pooled)") +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "bottom") -> p.pval

p.s9 <- patchwork::wrap_plots(p.trace, p.pval, ncol = 1, heights = c(1.4, 1))

bioRtools::write_data(p.s9, file = file.path(out_dir, "02_FigureS9_plant_heightgrowth_curve_and_effect_onset.pdf"),
                      width = 14, height = 9, units = "in")
bioRtools::write_data(p.s9, file = file.path(out_dir, "02_FigureS9_plant_heightgrowth_curve_and_effect_onset.png"),
                      width = 14, height = 9, units = "in", dpi = 500)
cat("[02] 输出：02_FigureS9_plant_heightgrowth_curve_and_effect_onset.pdf / .png\n")

cat("\n[02] sessionInfo()\n")
sessionInfo()
