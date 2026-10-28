# ============================================================
# FigS7.R —— Figure S7（Exploratory comparison of fully detected hormone compounds；单面板）
# Source Data：只读 source_data/Source_Data.xlsx（sheet：FigS7；含每样本原值）
# 运行：cd manuscript && Rscript code/FigS7.R
# 输出：figures/FigS7.pdf / .png
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/FigS7.R
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/_vectorcraft.R")
sd_file <- "source_data/Source_Data.xlsx"
d <- readxl::read_excel(sd_file, sheet = "FigS7", skip = 2) %>%
  dplyr::mutate(group = factor(group), index = factor(index, levels = unique(index)))
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)

p <- d %>%
  ggplot2::ggplot(aes(x = index, y = value, fill = group)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.7, size = 1) +
  ggplot2::facet_wrap(~index, scales = "free_y", nrow = 2) +
  bioRtools::scale_fill_research() +
  ggplot2::labs(x = NULL, y = "Hormone concentration (a.u.)", fill = "Planting combination") +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 7, face = "bold"),
                 axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                 legend.position = "bottom", legend.key.size = ggplot2::unit(0.3, "cm"),
                 legend.text = ggplot2::element_text(size = 7))
bioRtools::write_data(p, file = "figures/panels/FigS7.pdf", width = 17.4, height = 17.4, units = "cm")

compose_vectorcraft(
  list(list(path = "figures/panels/FigS7.pdf", rect = c(2, 8, 170, 170), label = "", label_pos = c(2, 3.5))),
  width_mm = 174, height_mm = 178, out_pdf = "figures/FigS7.pdf")

cat("\n[FigS7] sessionInfo()\n")
sessionInfo()
