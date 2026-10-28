# ============================================================
# FigS13.R —— Figure S13（Volcano plot of the limma analysis；单面板）
# Source Data：只读 source_data/Source_Data.xlsx（sheet：FigS13）
# 运行：cd manuscript && Rscript code/FigS13.R
# 输出：figures/FigS13.pdf / .png
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/FigS13.R
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/_vectorcraft.R")
sd_file <- "source_data/Source_Data.xlsx"
d <- readxl::read_excel(sd_file, sheet = "FigS13", skip = 2) %>%
  dplyr::mutate(stage = factor(stage, levels = c("tillering", "flowering", "filling")),
                sig = dplyr::case_when(limma_fdr < 0.05 ~ "FDR < 0.05",
                                       limma_p < 0.05 ~ "P < 0.05",
                                       TRUE ~ "n.s."),
                sig = factor(sig, levels = c("FDR < 0.05", "P < 0.05", "n.s.")))
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)

# 诊断图数值自检：检查 P 值范围与显著性计数（形状结论必须由数字支撑）
cat("[FigS13] 每期 P<0.05 基因数：\n")
print(d %>% dplyr::group_by(stage) %>%
        dplyr::summarise(n = dplyr::n(), n_p05 = sum(limma_p < 0.05),
                         n_fdr05 = sum(limma_fdr < 0.05), min_p = min(limma_p), .groups = "drop") %>%
        as.data.frame())

p <- d %>%
  ggplot2::ggplot(aes(x = limma_logFC, y = -log10(limma_p))) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  ggplot2::geom_point(aes(color = sig), size = 1.8, alpha = 0.85) +
  ggrepel::geom_text_repel(data = d %>% dplyr::filter(limma_fdr < 0.05 | limma_p < 0.01),
                           aes(label = gene), size = 2.4, fontface = "italic",
                           max.overlaps = 25, seed = 12345) +
  ggplot2::facet_wrap(~stage) +
  ggplot2::scale_color_manual(values = c("FDR < 0.05" = "#a63446", "P < 0.05" = "#e0a458",
                                         "n.s." = "grey60"), name = NULL) +
  ggplot2::labs(x = "log2 fold change (Mix vs Pure, focal-adjusted)",
                y = "-log10 P (moderated t)") +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 8),
                 legend.position = "bottom", legend.key.size = ggplot2::unit(0.35, "cm"),
                 legend.text = ggplot2::element_text(size = 7),
                 axis.title = ggplot2::element_text(size = 7.5))
bioRtools::write_data(p, file = "figures/panels/FigS13.pdf", width = 17.0, height = 8.5, units = "cm")

compose_vectorcraft(
  list(list(path = "figures/panels/FigS13.pdf", rect = c(2, 8, 170, 85), label = "", label_pos = c(2, 3.5))),
  width_mm = 174, height_mm = 93, out_pdf = "figures/FigS13.pdf")

cat("\n[FigS13] sessionInfo()\n")
sessionInfo()
