# ============================================================
# 06_设计分解与稳健性 / 02_株高_生长曲线与效应起始.R
# 独立分析（不属于旧版流程）：把株高当纵向性状处理，回答两个旧版没问的问题
#   Q1 互作效应改变的是生长速率还是最终水平（group x day 交互）
#   Q2 效应从第几次调查开始出现（逐时间点检验 + FDR 控制，给出 onset 时间）
# 旧版把 15 个时间点压成一个变化指数再做一次 ANOVA，纵向信息全部丢弃，既无法区分速率效应与水平效应，
# 也给不出效应出现的时间窗口。
#
# 运行：/opt/homebrew/bin/Rscript "<root>/code/06_设计分解与稳健性/02_株高_生长曲线与效应起始.R"
# 输入：<root>/data/01_分析输入表/03_分蘖数_株高_单株产量.xlsx
# 输出：<root>/result/06_设计分解与稳健性/02_表_株高_生长曲线与起始点.xlsx
#       <root>/result/06_设计分解与稳健性/02_FigureS9_株高生长曲线与效应起始.pdf / .png
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 口径：只用有 time 的株高记录（focal 36/74 无 time，自动排除）；time 为 YYYYMMDD 整数，
#   换算为距首次调查（20220711）的天数作为连续变量；Pure 组 neighbor = focal。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录见 code/06_设计分解与稳健性/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_in <- file.path(project_root, "data", "01_分析输入表", "03_分蘖数_株高_单株产量.xlsx")
out_dir <- file.path(project_root, "result", "06_设计分解与稳健性")
out_pub_dir <- file.path(project_root, "result", "08_投稿规格图")  # 投稿规格面板输出目录
dir.create(out_pub_dir, showWarnings = FALSE, recursive = TRUE)
path_out <- file.path(out_pub_dir, "98_重出副本_02_表_株高_生长曲线与起始点.xlsx")  # 投稿规格版：表写到这里，不动原目录

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

cat("[02] 株高数据：", nrow(df), " 行；时间点 ", dplyr::n_distinct(df$day),
    " 个（", min(df$day), "-", max(df$day), " 天）；组合 ", dplyr::n_distinct(df$combo), " 个\n")

# ---------- Q1：速率 vs 水平 ----------
rate_tab <- local({
  # 随机斜率（1 + day | combo）在 16 个组合下不可识别（large eigenvalue / max|grad 超限），
  # 故随机效应只保留组合截距；速率效应由 group:day 固定效应承担。
  # 注意：系数表里有一列名为 df，若在 mutate() 里直接用 df$combo 会被数据掩码解析成该列，
  # 因此这里先把样本量与组合数取成普通标量。
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

cat("\n[02] 各 focal 的速率效应（groupMix:day = 每日株高差的变化率）：\n")
print(as.data.frame(rate_by_focal), row.names = FALSE)

# ---------- Q2：逐时间点检验 + onset ----------
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

cat("\n[02] Q2 逐时间点检验（前 6 行）：\n")
print(as.data.frame(time_tab %>% dplyr::slice_head(n = 6)), row.names = FALSE)
cat("首次原始 P < 0.05：", ifelse(nrow(onset_raw) > 0, paste0("第 ", onset_raw$day, " 天（", onset_raw$date, "）"), "无"),
    "；首次 FDR < 0.05：", ifelse(nrow(onset_fdr) > 0, paste0("第 ", onset_fdr$day, " 天（", onset_fdr$date, "）"), "无"), "\n")
cat("FDR < 0.05 的时间点数：", sum(time_tab$t_fdr < 0.05), " / ", nrow(time_tab), "\n")
cat("原始 P < 0.05 的时间点数：", sum(time_tab$t_p_value < 0.05), " / ", nrow(time_tab), "\n")

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

# 投稿规格：两个面板分别导出（组图交 VectorCraft）
bioRtools::write_data(p.trace, file = file.path(out_pub_dir, "02a_FigureS9_生长曲线.pdf"),
                      width = 17.0, height = 9.5, units = "cm")
bioRtools::write_data(p.trace, file = file.path(out_pub_dir, "02a_FigureS9_生长曲线.png"),
                      width = 17.0, height = 9.5, units = "cm", dpi = 600)
bioRtools::write_data(p.pval, file = file.path(out_pub_dir, "02b_FigureS9_效应起始P曲线.pdf"),
                      width = 17.0, height = 6.8, units = "cm")
bioRtools::write_data(p.pval, file = file.path(out_pub_dir, "02b_FigureS9_效应起始P曲线.png"),
                      width = 17.0, height = 6.8, units = "cm", dpi = 600)
cat("[10] 输出：02a / 02b_FigureS9 面板（PDF + PNG）\n")

cat("\n[02] sessionInfo()\n")
sessionInfo()
