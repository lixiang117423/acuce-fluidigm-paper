# ============================================================
# FigS18.R —— Figure S18（Gene-level heatmap of the 92-gene panel；单面板）
# Source Data：只读 data/source_data/Source_Data_v45.xlsx（sheet：FigS18）
# 运行：Rscript code/figure_generation/FigS18.R
# 输出：figures/FigS18.pdf / .png
# 本机绝对路径：<repository_root>/manuscript/code/FigS18.R
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/figure_generation/_vectorcraft.R")
sd_file <- "data/source_data/Source_Data_v45.xlsx"
d <- readxl::read_excel(sd_file, sheet = "FigS18", skip = 2) %>%
  dplyr::mutate(stage = factor(stage, levels = c("tillering", "flowering", "filling")),
                gene = factor(gene, levels = rev(levels(factor(gene)))),
                combo_label = factor(combo_label, levels = unique(combo_label)))
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)

p <- d %>%
  ggplot2::ggplot(aes(x = combo_label, y = gene, fill = scaled_expression)) +
  ggplot2::geom_tile() +
  ggplot2::facet_grid(. ~ stage, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_gradientn(colours = c("#f2f2f0", "#a8bedb", "#1e4e8c"),
                                name = "Scaled median expression (0-1)") +
  ggplot2::labs(x = "Focal-neighbour combination", y = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
                 axis.text.y = ggplot2::element_text(size = 5.2, face = "italic"),
                 axis.title.x = ggplot2::element_text(size = 7.5),
                 strip.text = ggplot2::element_text(size = 8, face = "bold"),
                 panel.grid = ggplot2::element_blank(),
                 legend.key.width = ggplot2::unit(0.3, "cm"),
                 legend.key.height = ggplot2::unit(1.2, "cm"),
                 legend.title = ggplot2::element_text(size = 6.5, angle = 90),
                 legend.text = ggplot2::element_text(size = 6))
bioRtools::write_data(p, file = "figures/panels/FigS18.pdf", width = 17.0, height = 22.0, units = "cm")

compose_vectorcraft(
  list(list(path = "figures/panels/FigS18.pdf", rect = c(2, 8, 170, 220), label = "", label_pos = c(2, 3.5))),
  width_mm = 174, height_mm = 228, out_pdf = "figures/FigS18.pdf")

cat("\n[FigS18] sessionInfo()\n")
sessionInfo()
