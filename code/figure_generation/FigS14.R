# ============================================================
# FigS17.R —— Figure S17（Strigolactone pathway genes across the three stages；单面板）
# Source Data：只读 source_data/Source_Data.xlsx（sheet：FigS17）
# 运行：cd manuscript && Rscript code/FigS17.R
# 输出：figures/FigS17.pdf / .png
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/FigS17.R
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/_vectorcraft.R")
sd_file <- "source_data/Source_Data.xlsx"
d <- readxl::read_excel(sd_file, sheet = "FigS17", skip = 2) %>%
  dplyr::mutate(stage = factor(stage, levels = c("tillering", "flowering", "filling")),
                sig = dplyr::case_when(limma_fdr < 0.05 ~ "FDR < 0.05",
                                       limma_p < 0.05 ~ "P < 0.05",
                                       TRUE ~ "n.s."),
                sig = factor(sig, levels = c("FDR < 0.05", "P < 0.05", "n.s.")))
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)

p <- d %>%
  ggplot2::ggplot(aes(x = limma_logFC, y = forcats::fct_reorder(gene, limma_logFC))) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  ggplot2::geom_point(aes(color = sig), size = 3) +
  ggplot2::geom_text(aes(label = formatC(limma_logFC, format = "f", digits = 2)),
                     hjust = -0.25, size = 2.6, color = "grey25") +
  ggplot2::scale_color_manual(values = c("FDR < 0.05" = "#a63446", "P < 0.05" = "#e0a458",
                                         "n.s." = "grey62"), name = NULL) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.22)) +
  # 2026-10-19：关掉面板 clip——数值标签贴到分面右界时会被裁（实测 "0.11" 少画 3.36 pt）
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::facet_wrap(~stage, scales = "free_y") +
  ggplot2::labs(x = "Moderated log2 fold change (Mix vs Pure)", y = NULL,
                title = "Strigolactone genes (3)") +
  bioRtools::theme_prism() +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8, face = "italic"),
                 legend.position = "bottom", legend.key.size = ggplot2::unit(0.35, "cm"),
                 legend.text = ggplot2::element_text(size = 7),
                 plot.title = ggplot2::element_text(size = 8, face = "bold"))
bioRtools::write_data(p, file = "figures/panels/FigS17.pdf", width = 17.0, height = 6.5, units = "cm")

compose_vectorcraft(
  list(list(path = "figures/panels/FigS17.pdf", rect = c(2, 8, 170, 65), label = "", label_pos = c(2, 3.5))),
  width_mm = 174, height_mm = 73, out_pdf = "figures/FigS17.pdf")

cat("\n[FigS17] sessionInfo()\n")
sessionInfo()
