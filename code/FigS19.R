# ============================================================
# FigS19.R —— Figure S19（Accession-level expression of candidate pathway genes；2 个子图 a、b）
# Source Data：只读 source_data/Source_Data.xlsx（sheet：FigS19a、FigS19b）
# 运行：cd manuscript && Rscript code/FigS19.R
# 输出：figures/FigS19.pdf / .png；面板 figures/panels/FigS19a–b.pdf
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/FigS19.R
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/_vectorcraft.R")
sd_file <- "source_data/Source_Data.xlsx"
rd <- function(sheet) readxl::read_excel(sd_file, sheet = sheet, skip = 2)
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)
sig_colors <- c("FDR < 0.05" = "#a63446", "P < 0.05" = "#e0a458", "n.s." = "grey62")

dot_panel <- function(sheet, title, name) {
  d <- rd(sheet) %>%
    dplyr::mutate(sig = dplyr::case_when(focal_fdr < 0.05 ~ "FDR < 0.05",
                                         focal_p < 0.05 ~ "P < 0.05",
                                         TRUE ~ "n.s."),
                  sig = factor(sig, levels = names(sig_colors)))
  p <- d %>%
    ggplot2::ggplot(aes(x = focal_logFC, y = forcats::fct_reorder(gene, focal_logFC))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
    ggplot2::geom_point(aes(color = sig), size = 3) +
    ggplot2::geom_text(aes(label = formatC(focal_logFC, format = "f", digits = 2)),
                       hjust = -0.25, size = 2.2, color = "grey25") +
    ggplot2::scale_color_manual(values = sig_colors, name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.22)) +
    ggplot2::labs(x = "Moderated log2 fold change (Mix vs Pure)", y = NULL, title = title) +
    bioRtools::theme_prism() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6.3, face = "italic"),
                   legend.position = "bottom", legend.key.size = ggplot2::unit(0.35, "cm"),
                   legend.text = ggplot2::element_text(size = 7),
                   plot.title = ggplot2::element_text(size = 8, face = "bold"))
  bioRtools::write_data(p, file = file.path("figures/panels", paste0(name, ".pdf")),
                        width = 17.0, height = 6.5, units = "cm")
}
dot_panel("FigS19a", "Focal 69, tillering: BR + auxin", "FigS19a")
dot_panel("FigS19b", "Focal 87, grain filling: BR + SA + drought", "FigS19b")

compose_vectorcraft(
  list(
    list(path = "figures/panels/FigS19a.pdf", rect = c(2, 8, 170, 65), label = "A", label_pos = c(2, 3.5)),
    list(path = "figures/panels/FigS19b.pdf", rect = c(2, 78, 170, 65), label = "B", label_pos = c(2, 73.5))
  ),
  width_mm = 174, height_mm = 148, out_pdf = "figures/FigS19.pdf")

cat("\n[FigS19] sessionInfo()\n")
sessionInfo()
