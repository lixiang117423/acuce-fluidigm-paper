# ============================================================
# FigS12.R —— Figure S12（Shape of the plant height response and change point；3 个子图 a–c）
# Source Data：只读 data/source_data/Source_Data_v45.xlsx（sheet：FigS12a、FigS12b、FigS12c）
# 运行：Rscript code/figure_generation/FigS12.R
# 输出：figures/FigS12.pdf / .png；面板 figures/panels/FigS12a–c.pdf
# 本机绝对路径：<repository_root>/manuscript/code/FigS12.R
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/figure_generation/_vectorcraft.R")
sd_file <- "data/source_data/Source_Data_v45.xlsx"
rd <- function(sheet) readxl::read_excel(sd_file, sheet = sheet, skip = 2)
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)
save_panel <- function(p, name, w_cm, h_cm) {
  bioRtools::write_data(p, file = file.path("figures/panels", paste0(name, ".pdf")),
                        width = w_cm, height = h_cm, units = "cm")
}

d12a <- rd("FigS12a")
tau_hat <- d12a$tau_hat[1]
p12a <- d12a %>%
  ggplot2::ggplot(aes(x = day, y = diff_mix_minus_pure)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  ggplot2::geom_vline(xintercept = tau_hat, linetype = "dotted", color = "#a63446") +
  ggplot2::geom_errorbar(aes(ymin = diff_mix_minus_pure - 1.96 * se_diff,
                             ymax = diff_mix_minus_pure + 1.96 * se_diff),
                         width = 1.5, color = "grey35") +
  ggplot2::geom_point(size = 2.2, color = "#1e4e8c") +
  ggplot2::labs(x = "Days after first survey", y = "Plant height difference (Mix - Pure, cm)",
                subtitle = paste0("Change point: day ", tau_hat, "; 95% set: days ",
                                  d12a$ci95_low[1], "-", d12a$ci95_high[1])) +
  bioRtools::theme_prism() +
  ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 7),
                 axis.title = ggplot2::element_text(size = 7))
save_panel(p12a, "FigS12a", 8.5, 6.5)

d12b <- rd("FigS12b")
p12b <- d12b %>%
  ggplot2::ggplot(aes(x = tau, y = b_group)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  ggplot2::geom_line(color = "#1e4e8c", linewidth = 0.7) +
  ggplot2::geom_point(size = 1.6, color = "#1e4e8c") +
  # 2026-10-19：旋转 y 轴title顶端被切 → 缩短 + 缩小字号 + 留边距
  ggplot2::labs(x = "Candidate change point (day)", y = "Group effect after change") +
  bioRtools::theme_prism() +
  # 注意：theme() 必须写在 theme_prism() 之后，否则会被 theme_prism 的轴设置覆盖
  ggplot2::theme(axis.title = ggplot2::element_text(size = 6.5),
                 plot.margin = ggplot2::margin(3, 2, 3, 2))
save_panel(p12b, "FigS12b", 8.5, 6.5)

d12c <- rd("FigS12c")
pct_lab <- paste0("Model preference: step ", round(d12c$pct_step[1], 1), "%, linear ",
                  round(d12c$pct_linear[1], 1), "%, constant ", round(d12c$pct_constant[1], 1), "%")
p12c <- d12c %>%
  ggplot2::ggplot(aes(x = tau)) +
  ggplot2::geom_histogram(binwidth = 1, fill = "#1e4e8c", alpha = 0.85) +
  ggplot2::labs(x = "Bootstrap change point (day)", y = "Count", subtitle = pct_lab) +
  bioRtools::theme_prism() +
  ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8))
save_panel(p12c, "FigS12c", 17.0, 6.0)

compose_vectorcraft(
  list(
    list(path = "figures/panels/FigS12a.pdf", rect = c(2, 8, 85, 65), label = "A", label_pos = c(2, 3.5)),
    list(path = "figures/panels/FigS12b.pdf", rect = c(89, 8, 85, 65), label = "B", label_pos = c(89, 3.5)),
    list(path = "figures/panels/FigS12c.pdf", rect = c(2, 78, 170, 60), label = "C", label_pos = c(2, 73.5))
  ),
  width_mm = 174, height_mm = 140, out_pdf = "figures/FigS12.pdf")

cat("\n[FigS12] sessionInfo()\n")
sessionInfo()
