# ============================================================
# Fig4.R —— Figure 4（Cross-stage expression: defence genes + high-sugar genes）
# Source Data：只读 source_data/Source_Data.xlsx（sheet：Fig4a、Fig4b）
# 运行：cd manuscript && Rscript code/Fig4.R
# 输出：figures/Fig4.pdf / .png；面板 figures/panels/Fig4a–b.pdf
# 软件：R 4.6.1；ggplot2；bioRtools::theme_prism()；set.seed(12345)
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/Fig4.R
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
stage_levels <- c("tillering", "flowering", "filling")
sig_colors <- c("FDR < 0.05" = "#a63446", "P < 0.05" = "#e0a458", "n.s." = "grey62")

panel_logfc <- function(sheet, title, name, w_cm, h_cm) {
  d <- rd(sheet) %>%
    dplyr::mutate(stage = factor(stage, levels = stage_levels),
                  sig = dplyr::case_when(limma_fdr < 0.05 ~ "FDR < 0.05",
                                         limma_p < 0.05 ~ "P < 0.05",
                                         TRUE ~ "n.s."),
                  sig = factor(sig, levels = names(sig_colors)))
  p <- d %>%
    ggplot2::ggplot(aes(x = limma_logFC, y = forcats::fct_reorder(gene, limma_logFC))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
    ggplot2::geom_point(aes(color = sig), size = 3) +
    ggplot2::geom_text(aes(label = formatC(limma_logFC, format = "f", digits = 2)),
                       hjust = -0.25, size = 2.2, color = "grey25") +
    ggplot2::scale_color_manual(values = sig_colors, name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.22)) +
    ggplot2::facet_wrap(~stage, scales = "free_y") +
    ggplot2::labs(x = "Moderated log2 fold change (Mix vs Pure)", y = NULL, title = title) +
    bioRtools::theme_prism() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6.5, face = "italic"),
                   panel.spacing = ggplot2::unit(0.5, "cm"),
                   legend.position = "bottom", legend.key.size = ggplot2::unit(0.35, "cm"),
                   legend.text = ggplot2::element_text(size = 7),
                   plot.title = ggplot2::element_text(size = 8, face = "bold"))
  save_panel(p, name, w_cm, h_cm)
}
panel_logfc("Fig4a", "Defence genes (12)", "Fig4a", 17.0, 9.0)
panel_logfc("Fig4b", "High-sugar genes (5)", "Fig4b", 17.0, 6.0)

compose_vectorcraft(
  list(
    list(path = "figures/panels/Fig4a.pdf", rect = c(2, 8, 170, 90), label = "A", label_pos = c(2, 3.5)),
    list(path = "figures/panels/Fig4b.pdf", rect = c(2, 103, 170, 60), label = "B", label_pos = c(2, 98.5))
  ),
  width_mm = 174, height_mm = 165, out_pdf = "figures/Fig4.pdf")

cat("\n[Fig4] sessionInfo()\n")
sessionInfo()
