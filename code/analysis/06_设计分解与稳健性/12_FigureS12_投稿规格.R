# ============================================================
# 06_设计分解与稳健性 / 05_株高_效应形状与变化点.R
# 独立分析（续 06 主题，Q7）：株高互作效应的"时间形状"是什么，以及效应何时出现
#   旧版只给一个总变化指数；主题 06 的 Q4 给出"速率效应 + 首次调查即显著"，
#   本脚本进一步用**变化点（change point）模型**回答：效应是阶跃型（某时点后恒定抬高）、
#   线性衰减型（早期差异逐渐消失）还是恒定型（全程不变）？
# 方法：
#   1) 植物水平模型（combo 作为固定区组，吸收组合均值差异）：
#        value ~ group * post + combo，其中 post = I(day >= tau)
#      对 15 个候选 tau（观测时点）逐一拟合，按对数似然挑最优 tau；
#   2) 变化点的 95% 置信区间：似然比统计量 2*(logLik_max − logLik(tau)) < 3.84 的 tau 集合；
#   3) 形状比较：M0 恒定差（value ~ group + combo）、M1 线性差（value ~ group * day + combo）、
#      M2 阶跃差（group * post，取最优 tau），比较 AIC；
#   4) 自助法（199 次，按 时点 × 组合 × 处理 分层重抽样）给出 tau 的不确定性与模型选择稳定性。
#
# 运行：/opt/homebrew/bin/Rscript "<root>/code/06_设计分解与稳健性/05_株高_效应形状与变化点.R"
# 输入：<root>/data/01_分析输入表/03_分蘖数_株高_单株产量.xlsx
# 输出：<root>/result/06_设计分解与稳健性/05_表_株高_效应形状与变化点.xlsx
#       <root>/result/06_设计分解与稳健性/05_FigureS12_株高效应形状与变化点.pdf / .png
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
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
path_out <- file.path(out_pub_dir, "98_重出副本_05_表_株高_效应形状与变化点.xlsx")  # 投稿规格版：表写到这里，不动原目录

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stopifnot(file.exists(path_in))

bioRtools::read_data(path_in) -> d3

d3 %>%
  dplyr::filter(phenotype == "plant height", !is.na(time)) %>%
  dplyr::mutate(
    group = factor(group, levels = c("Pure", "Mix")),
    date = as.Date(as.character(time), format = "%Y%m%d"),
    day = as.numeric(date - min(date)),
    combo = factor(paste0(focal, "-", neighbor))
  ) -> df

# 候选变化点必须严格大于首次调查日：若 tau = min(day)，则 post = I(day >= tau) 恒为 1，
# 阶跃项与交互项被 lm 别名化（系数不存在，取值会下标出界）。
# 注意 I(day < tau) 与 I(day >= tau) 只是同一划分的两种写法，故只需在一个方向上枚举 tau，
# 其中 tau = 第 2 个时点即代表"仅首次调查差异"这一划分。
candidate_tau <- sort(unique(df$day))
candidate_tau <- candidate_tau[candidate_tau > min(candidate_tau)]
cat("[05] 数据：", nrow(df), " 行；候选变化点 ", length(candidate_tau), " 个：",
    paste(candidate_tau, collapse = ", "), "\n")

# ---- 变化点网格搜索 ----
fit_tau <- function(d, tau) {
  d$post <- as.integer(d$day >= tau)
  m <- stats::lm(value ~ group * post + combo, data = d)
  s <- summary(m)$coefficients
  ll <- stats::logLik(m)
  tibble::tibble(
    tau = tau,
    logLik = as.numeric(ll),
    df_model = attr(ll, "df"),
    aic = stats::AIC(m),
    b_group = s["groupMix", "Estimate"],
    p_group = s["groupMix", "Pr(>|t|)"],
    b_group_post = s["groupMix:post", "Estimate"],
    p_group_post = s["groupMix:post", "Pr(>|t|)"],
    n_obs = nrow(d)
  )
}

grid_tab <- purrr::map_dfr(candidate_tau, ~ fit_tau(df, .x)) %>%
  dplyr::arrange(dplyr::desc(logLik))

tau_hat <- grid_tab$tau[1]
ll_max <- grid_tab$logLik[1]

grid_tab <- grid_tab %>%
  dplyr::mutate(
    lr_stat = 2 * (ll_max - logLik),
    in_ci95 = lr_stat < 3.84
  ) %>%
  dplyr::arrange(tau)

ci_set <- grid_tab %>% dplyr::filter(in_ci95) %>% dplyr::pull(tau)

cat("\n[04->05] 最优变化点 tau =", tau_hat, " 天（", as.character(min(df$date) + tau_hat), "）\n")
cat("[05] 变化点 95% 似然比区间覆盖的候选时点：", paste(ci_set, collapse = ", "), "\n")
cat("[05] 最优 tau 下的阶跃效应：group:post = ",
    round(grid_tab$b_group_post[grid_tab$tau == tau_hat], 3),
    "，P = ", formatC(grid_tab$p_group_post[grid_tab$tau == tau_hat], format = "f", digits = 4), "\n")

# ---- 三种形状的模型比较 ----
fit_shape <- function(d) {
  d$post <- as.integer(d$day >= tau_hat)
  m0 <- stats::lm(value ~ group + combo, data = d)                       # 恒定差
  m1 <- stats::lm(value ~ group * day + combo, data = d)                 # 线性差
  m2 <- stats::lm(value ~ group * post + combo, data = d)                # 阶跃差

  tibble::tibble(
    shape = c("constant difference", "linear difference", "step change"),
    model = c("value ~ group + combo", "value ~ group * day + combo",
              paste0("value ~ group * post + combo (tau = ", tau_hat, ")")),
    aic = c(stats::AIC(m0), stats::AIC(m1), stats::AIC(m2)),
    logLik = c(as.numeric(stats::logLik(m0)), as.numeric(stats::logLik(m1)),
               as.numeric(stats::logLik(m2))),
    n_par = c(attr(stats::logLik(m0), "df"), attr(stats::logLik(m1), "df"),
              attr(stats::logLik(m2), "df")),
    # 关键效应量：恒定差模型给组间总差；线性模型给 group:day 斜率；阶跃模型给 group:post 阶跃
    key_effect = c(
      summary(m0)$coefficients["groupMix", "Estimate"],
      summary(m1)$coefficients["groupMix:day", "Estimate"],
      summary(m2)$coefficients["groupMix:post", "Estimate"]
    ),
    key_effect_p = c(
      summary(m0)$coefficients["groupMix", "Pr(>|t|)"],
      summary(m1)$coefficients["groupMix:day", "Pr(>|t|)"],
      summary(m2)$coefficients["groupMix:post", "Pr(>|t|)"]
    )
  ) %>%
    dplyr::mutate(
      delta_aic = aic - min(aic),
      aic_weight = exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
    )
}

shape_tab <- fit_shape(df)

cat("\n[05] 三种时间形状的模型比较：\n")
print(as.data.frame(shape_tab), row.names = FALSE)

# ---- 自助法：tau 与形状的不确定性 ----
boot_one <- function(b) {
  set.seed(1000 + b)
  d <- df %>%
    dplyr::group_by(day, combo, group) %>%
    dplyr::slice_sample(prop = 1, replace = TRUE) %>%
    dplyr::ungroup()

  g <- purrr::map_dfr(candidate_tau, ~ fit_tau(d, .x)) %>% dplyr::arrange(dplyr::desc(logLik))
  tau_b <- g$tau[1]

  d$post <- as.integer(d$day >= tau_b)
  a_m2 <- stats::AIC(stats::lm(value ~ group * post + combo, data = d))
  a_m1 <- stats::AIC(stats::lm(value ~ group * day + combo, data = d))
  a_m0 <- stats::AIC(stats::lm(value ~ group + combo, data = d))

  tibble::tibble(
    boot = b,
    tau = tau_b,
    best_shape = c("constant difference", "linear difference", "step change")[which.min(c(a_m0, a_m1, a_m2))],
    b_group_post = g$b_group_post[1],
    p_group_post = g$p_group_post[1]
  )
}

n_boot <- 199
boot_tab <- purrr::map_dfr(seq_len(n_boot), boot_one)

boot_summary <- tibble::tibble(
  n_boot = n_boot,
  tau_mode = names(sort(table(boot_tab$tau), decreasing = TRUE))[1],
  tau_q025 = quantile(boot_tab$tau, 0.025),
  tau_q500 = quantile(boot_tab$tau, 0.5),
  tau_q975 = quantile(boot_tab$tau, 0.975),
  pct_step = 100 * mean(boot_tab$best_shape == "step change"),
  pct_linear = 100 * mean(boot_tab$best_shape == "linear difference"),
  pct_constant = 100 * mean(boot_tab$best_shape == "constant difference"),
  pct_significant_step = 100 * mean(boot_tab$p_group_post < 0.05),
  median_b_group_post = median(boot_tab$b_group_post)
)

cat("\n[05] 自助法（", n_boot, " 次）结果：\n", sep = "")
print(as.data.frame(boot_summary), row.names = FALSE)

bioRtools::write_data(
  list(change_point_grid = grid_tab, shape_comparison = shape_tab,
       bootstrap_runs = boot_tab, bootstrap_summary = boot_summary),
  file = path_out
)
cat("[05] 输出：", path_out, "\n")

# ---- 图 ----
diff_tab <- df %>%
  dplyr::group_by(day, date, group) %>%
  dplyr::summarise(mean_value = mean(value), se = sd(value) / sqrt(dplyr::n()), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = group, values_from = c(mean_value, se)) %>%
  dplyr::mutate(
    diff = mean_value_Mix - mean_value_Pure,
    se_diff = sqrt(se_Mix^2 + se_Pure^2)
  )

p.diff <- diff_tab %>%
  ggplot2::ggplot(aes(x = day, y = diff)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  ggplot2::geom_vline(xintercept = tau_hat, linetype = "dotted", color = "#a63446") +
  ggplot2::geom_errorbar(aes(ymin = diff - 1.96 * se_diff, ymax = diff + 1.96 * se_diff),
                         width = 2, color = "grey35") +
  ggplot2::geom_point(size = 2.4, color = "#1e4e8c") +
  ggplot2::labs(
    x = "Days after first survey",
    y = "Plant height difference (Mix - Pure, cm)",
    subtitle = paste0("Change point: day ", tau_hat, "; 95% set: days ",
                      paste(range(ci_set), collapse = "-"))
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 7),
                 axis.title = ggplot2::element_text(size = 7)) -> p.diff

p.grid <- grid_tab %>%
  ggplot2::ggplot(aes(x = tau, y = lr_stat)) +
  ggplot2::geom_hline(yintercept = 3.84, linetype = "dashed", color = "grey40") +
  ggplot2::geom_line(color = "#1e4e8c", linewidth = 0.7) +
  ggplot2::geom_point(aes(color = in_ci95), size = 2.2) +
  ggplot2::scale_color_manual(values = c(`TRUE` = "#a63446", `FALSE` = "grey60"),
                              labels = c(`TRUE` = "in 95% set", `FALSE` = "outside"), name = NULL) +
  # 2026-10-19：旋转 y 轴标题在面板顶端被切（缺 2.2 pt）→ 缩短 + 缩小字号 + 留边距
  ggplot2::labs(x = "Candidate change point (day)", y = "Likelihood ratio") +
  bioRtools::theme_prism() +
  ggplot2::theme(axis.title = ggplot2::element_text(size = 6.5),
                 plot.margin = ggplot2::margin(3, 2, 3, 2)) +
  ggplot2::theme(legend.position = "bottom") -> p.grid

p.boot <- boot_tab %>%
  ggplot2::ggplot(aes(x = tau)) +
  ggplot2::geom_histogram(binwidth = 1, fill = "#1e4e8c", alpha = 0.85) +
  ggplot2::labs(x = "Bootstrap change point (day)", y = "Count",
                subtitle = paste0("Model preference: step ", round(boot_summary$pct_step, 1),
                                  "%, linear ", round(boot_summary$pct_linear, 1),
                                  "%, constant ", round(boot_summary$pct_constant, 1), "%")) +
  bioRtools::theme_prism() +
  ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8)) -> p.boot

# 投稿规格：三个面板分别导出（组图交 VectorCraft）
bioRtools::write_data(p.diff, file = file.path(out_pub_dir, "05a_FigureS12_效应形状.pdf"),
                      width = 8.3, height = 6.5, units = "cm")
bioRtools::write_data(p.diff, file = file.path(out_pub_dir, "05a_FigureS12_效应形状.png"),
                      width = 8.3, height = 6.5, units = "cm", dpi = 600)
bioRtools::write_data(p.grid, file = file.path(out_pub_dir, "05b_FigureS12_变化点网格.pdf"),
                      width = 8.3, height = 6.5, units = "cm")
bioRtools::write_data(p.grid, file = file.path(out_pub_dir, "05b_FigureS12_变化点网格.png"),
                      width = 8.3, height = 6.5, units = "cm", dpi = 600)
bioRtools::write_data(p.boot, file = file.path(out_pub_dir, "05c_FigureS12_自助分布.pdf"),
                      width = 17.0, height = 6.0, units = "cm")
bioRtools::write_data(p.boot, file = file.path(out_pub_dir, "05c_FigureS12_自助分布.png"),
                      width = 17.0, height = 6.0, units = "cm", dpi = 600)
cat("[12] 输出：05a / 05b / 05c_FigureS12 面板（PDF + PNG）\n")

cat("\n[05] sessionInfo()\n")
sessionInfo()
