# ============================================================
# FigS10.R —— Figure S10（Pathway score differences；单面板）
# Source Data：只读 data/source_data/Source_Data_v45.xlsx（sheet：FigS10）
# 运行：Rscript code/figure_generation/FigS10.R
# 输出：figures/FigS10.pdf / .png
# 本机绝对路径：<repository_root>/manuscript/code/FigS10.R
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/figure_generation/_vectorcraft.R")
sd_file <- "data/source_data/Source_Data_v45.xlsx"
d <- readxl::read_excel(sd_file, sheet = "FigS10", skip = 2) %>%
  dplyr::mutate(stage = factor(stage, levels = c("tillering", "flowering", "filling")),
                sig = dplyr::case_when(lmm_fdr < 0.05 ~ "FDR < 0.05",
                                       lmm_p < 0.05 ~ "P < 0.05",
                                       TRUE ~ "n.s."),
                sig = factor(sig, levels = c("FDR < 0.05", "P < 0.05", "n.s.")))
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)

p <- d %>%
  ggplot2::ggplot(aes(x = lmm_estimate, y = forcats::fct_reorder(process, lmm_estimate), color = sig)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  ggplot2::geom_errorbarh(aes(xmin = lmm_estimate - 1.96 * lmm_se, xmax = lmm_estimate + 1.96 * lmm_se),
                          height = 0.25, linewidth = 0.5) +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::facet_wrap(~stage, scales = "free_x") +
  ggplot2::scale_color_manual(values = c("FDR < 0.05" = "#a63446", "P < 0.05" = "#e0a458",
                                         "n.s." = "grey55"), name = NULL) +
  ggplot2::labs(x = "Pathway score difference (Mix - Pure, mean Z)", y = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8),
                 strip.text = ggplot2::element_text(face = "bold", size = 8),
                 legend.position = "bottom", legend.key.size = ggplot2::unit(0.35, "cm"),
                 legend.text = ggplot2::element_text(size = 7))
bioRtools::write_data(p, file = "figures/panels/FigS10.pdf", width = 17.4, height = 11.5, units = "cm")

compose_vectorcraft(
  list(list(path = "figures/panels/FigS10.pdf", rect = c(2, 8, 170, 115), label = "", label_pos = c(2, 3.5))),
  width_mm = 174, height_mm = 123, out_pdf = "figures/FigS10.pdf")

cat("\n[FigS10] sessionInfo()\n")
sessionInfo()
