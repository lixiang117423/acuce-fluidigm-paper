# ============================================================
# FigS9.R —— Figure S9（Plant height trajectories and onset of the interaction effect；2 个子图 a、b）
# Source Data：只读 source_data/Source_Data.xlsx（sheet：FigS9a、FigS9b）
# 运行：cd manuscript && Rscript code/FigS9.R
# 输出：figures/FigS9.pdf / .png；面板 figures/panels/FigS9a–b.pdf
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/FigS9.R
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/_vectorcraft.R")
sd_file <- "source_data/Source_Data.xlsx"
rd <- function(sheet) readxl::read_excel(sd_file, sheet = sheet, skip = 2)
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)
save_panel <- function(p, name, w_cm, h_cm) {
  bioRtools::write_data(p, file = file.path("figures/panels", paste0(name, ".pdf")),
                        width = w_cm, height = h_cm, units = "cm")
}

d9a <- rd("FigS9a")
p9a <- d9a %>%
  dplyr::select(day, mean_pure, mean_mix) %>%
  tidyr::pivot_longer(cols = -day, names_to = "series", values_to = "mean_height") %>%
  dplyr::mutate(series = factor(series, levels = c("mean_pure", "mean_mix"),
                                labels = c("Monoculture", "Mixture"))) %>%
  ggplot2::ggplot(aes(x = day, y = mean_height, color = series)) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 1.4) +
  bioRtools::scale_color_research() +
  ggplot2::labs(x = "Days after first survey", y = "Plant height (cm)", color = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "top", legend.key.size = ggplot2::unit(0.35, "cm"))
save_panel(p9a, "FigS9a", 17.0, 9.5)

d9b <- rd("FigS9b")
p9b <- d9b %>%
  ggplot2::ggplot(aes(x = day, y = t_p_value)) +
  ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey40") +
  ggplot2::geom_line(color = "#1e4e8c", linewidth = 0.7) +
  ggplot2::geom_point(size = 1.4, color = "#1e4e8c") +
  # 2026-10-19：旋转 y 轴标题比绘图区还长、顶端被切（缺 132 px 墨迹）→ 缩短 + 缩小字号 + 留边距
  ggplot2::labs(x = "Days after first survey", y = "P value (Mix vs Pure)") +
  bioRtools::theme_prism() +
  # theme() 必须写在 theme_prism() 之后，否则会被 theme_prism 的轴设置覆盖
  ggplot2::theme(axis.title = ggplot2::element_text(size = 6.5),
                 plot.margin = ggplot2::margin(3, 2, 3, 2))
save_panel(p9b, "FigS9b", 17.0, 6.8)

compose_vectorcraft(
  list(
    list(path = "figures/panels/FigS9a.pdf", rect = c(2, 8, 170, 95), label = "A", label_pos = c(2, 3.5)),
    list(path = "figures/panels/FigS9b.pdf", rect = c(2, 108, 170, 68), label = "B", label_pos = c(2, 103.5))
  ),
  width_mm = 174, height_mm = 178, out_pdf = "figures/FigS9.pdf")

cat("\n[FigS9] sessionInfo()\n")
sessionInfo()
