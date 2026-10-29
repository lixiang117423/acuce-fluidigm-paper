# ============================================================
# FigS3.R —— Figure S3（VIP genes and genetic distance；单面板，无面板字母）
# Source Data：只读 source_data/Source_Data.xlsx（sheet：FigS3）
# 运行：cd manuscript && Rscript code/FigS3.R
# 输出：figures/FigS3.pdf / .png；面板 figures/panels/FigS3.pdf
# 统计：每基因 × 每时期对遗传距离做 Pearson 相关（12 个混种组合），FDR 按时期校正
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/FigS3.R
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/_vectorcraft.R")
sd_file <- "source_data/Source_Data.xlsx"
d <- readxl::read_excel(sd_file, sheet = "FigS3", skip = 2) %>%
  dplyr::mutate(stage = factor(stage, levels = c("tillering", "flowering", "filling")),
                gene = factor(gene, levels = unique(gene))) %>%
  dplyr::mutate(lab = paste0("r = ", formatC(r, format = "f", digits = 2), ", ",
                             dplyr::if_else(p_value < 0.001, "P<0.001",
                                            paste0("P=", formatC(p_value, format = "f", digits = 2)))))
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)

p <- d %>%
  ggplot2::ggplot(aes(x = genetic_distance, y = expression)) +
  ggplot2::geom_point(size = 1.5, color = "#1e4e8c", alpha = 0.85) +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "grey35", linewidth = 0.5) +
  ggplot2::geom_text(data = d %>% dplyr::distinct(stage, gene, lab),
                     aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.6,
                     size = 2.2, inherit.aes = FALSE) +
  ggplot2::facet_wrap(~ stage + gene, ncol = 5, nrow = 6, scales = "free_y") +
  ggplot2::labs(x = "Genetic distance", y = "Relative expression (2^-dCt)") +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 6.5, face = "bold"),
                 axis.text = ggplot2::element_text(size = 6.5),
                 axis.title = ggplot2::element_text(size = 7.5))
bioRtools::write_data(p, file = "figures/panels/FigS3.pdf", width = 17.4, height = 22.0, units = "cm")

compose_vectorcraft(
  list(list(path = "figures/panels/FigS3.pdf", rect = c(2, 8, 170, 220), label = "", label_pos = c(2, 3.5))),
  width_mm = 174, height_mm = 228, out_pdf = "figures/FigS3.pdf")

cat("\n[FigS3] sessionInfo()\n")
sessionInfo()
