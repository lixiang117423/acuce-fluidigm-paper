# ============================================================
# FigS5.R —— Figure S5（High-sugar vs defence correlation heatmap；单面板）
# Source Data：只读 data/source_data/Source_Data_v45.xlsx（sheet：FigS5）
# 运行：Rscript code/figure_generation/FigS5.R
# 输出：figures/FigS5.pdf / .png
# 本机绝对路径：<repository_root>/manuscript/code/FigS5.R
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/figure_generation/_vectorcraft.R")
sd_file <- "data/source_data/Source_Data_v45.xlsx"
d <- readxl::read_excel(sd_file, sheet = "FigS5", skip = 2) %>%
  dplyr::mutate(stage = factor(stage, levels = c("tillering", "flowering", "filling")))
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)

p <- d %>%
  ggplot2::ggplot(aes(x = defense_gene, y = sugar_gene, fill = rho)) +
  ggplot2::geom_tile() +
  ggplot2::facet_wrap(~stage, ncol = 2) +
  ggplot2::scale_fill_gradient2(low = "#3B1F4E", mid = "#FFFFFF", high = "#780018",
                                midpoint = 0, name = "Spearman rho") +
  ggplot2::labs(x = "Defence-related gene", y = "High-sugar gene") +
  bioRtools::theme_prism() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
                 axis.text.y = ggplot2::element_text(size = 6, face = "italic"),
                 strip.text = ggplot2::element_text(size = 8, face = "bold"),
                 panel.grid = ggplot2::element_blank(),
                 legend.key.width = ggplot2::unit(0.3, "cm"),
                 legend.title = ggplot2::element_text(size = 6.5, angle = 90),
                 legend.text = ggplot2::element_text(size = 6))
bioRtools::write_data(p, file = "figures/panels/FigS5.pdf", width = 17.4, height = 12.0, units = "cm")

compose_vectorcraft(
  list(list(path = "figures/panels/FigS5.pdf", rect = c(2, 8, 170, 120), label = "", label_pos = c(2, 3.5))),
  width_mm = 174, height_mm = 128, out_pdf = "figures/FigS5.pdf")

cat("\n[FigS5] sessionInfo()\n")
sessionInfo()
