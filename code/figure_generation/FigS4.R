# ============================================================
# FigS4.R —— Figure S4（Field plant height and 1000-grain weight；单面板）
# Source Data：只读 data/source_data/Source_Data_v45.xlsx（sheet：FigS4）
# 运行：Rscript code/figure_generation/FigS4.R
# 输出：figures/FigS4.pdf / .png
# 本机绝对路径：<repository_root>/manuscript/code/FigS4.R
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/figure_generation/_vectorcraft.R")
sd_file <- "data/source_data/Source_Data_v45.xlsx"
d <- readxl::read_excel(sd_file, sheet = "FigS4", skip = 2) %>%
  dplyr::mutate(category = factor(category), focal = factor(as.character(focal)),
                condition = factor(condition, levels = c("Pure", "Mix")))
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)

p <- d %>%
  ggplot2::ggplot(aes(x = focal, y = emmean, fill = condition)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.62) +
  ggplot2::geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                         position = ggplot2::position_dodge(width = 0.7), width = 0.15) +
  ggplot2::facet_wrap(~category, scales = "free_y") +
  bioRtools::scale_fill_research() +
  ggplot2::labs(x = "Focal accession", y = NULL, fill = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 8, face = "bold"),
                 legend.position = "bottom", legend.key.size = ggplot2::unit(0.3, "cm"))
bioRtools::write_data(p, file = "figures/panels/FigS4.pdf", width = 17.4, height = 8.7, units = "cm")

compose_vectorcraft(
  list(list(path = "figures/panels/FigS4.pdf", rect = c(2, 8, 170, 87), label = "", label_pos = c(2, 3.5))),
  width_mm = 174, height_mm = 95, out_pdf = "figures/FigS4.pdf")

cat("\n[FigS4] sessionInfo()\n")
sessionInfo()
