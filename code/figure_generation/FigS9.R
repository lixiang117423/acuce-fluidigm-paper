# ============================================================
# FigS11.R —— Figure S11（Partner-effect variance decomposition；3 个子图 a–c）
# Source Data：只读 data/source_data/Source_Data_v45.xlsx（sheet：FigS11a、FigS11b、FigS11c）
# 运行：Rscript code/figure_generation/FigS11.R
# 输出：figures/FigS11.pdf / .png；面板 figures/panels/FigS11a–c.pdf
# 本机绝对路径：<repository_root>/manuscript/code/FigS11.R
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/figure_generation/_vectorcraft.R")
sd_file <- "data/source_data/Source_Data_v45.xlsx"
rd <- function(sheet) readxl::read_excel(sd_file, sheet = sheet, skip = 2)
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)
save_panel <- function(p, name, w_cm, h_cm) {
  bioRtools::write_data(p, file = file.path("figures/panels", paste0(name, ".pdf")),
                        width = w_cm, height = h_cm, units = "cm")
}
stage_levels <- c("tillering", "flowering", "filling")

d11a <- rd("FigS11a") %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels)) %>%
  dplyr::select(stage, gene, eta2_focal, eta2_neighbor, eta2_interaction) %>%
  tidyr::pivot_longer(cols = dplyr::starts_with("eta2_"), names_to = "term", values_to = "eta2") %>%
  dplyr::mutate(term = factor(term, levels = c("eta2_focal", "eta2_neighbor", "eta2_interaction"),
                              labels = c("Own genotype", "Partner genotype", "Specific pairing")))
p11a <- d11a %>%
  ggplot2::ggplot(aes(x = term, y = eta2, fill = term)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  ggplot2::facet_wrap(~stage) +
  bioRtools::scale_fill_bioinfo() +
  ggplot2::labs(x = NULL, y = expression("Partial " * eta^2 * " per gene"), fill = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "none",
                 axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1, size = 7),
                 strip.text = ggplot2::element_text(face = "bold", size = 8))
save_panel(p11a, "FigS11a", 8.5, 7.0)

d11b <- rd("FigS11b") %>% dplyr::mutate(stage = factor(stage, levels = stage_levels))
p11b <- d11b %>%
  ggplot2::ggplot(aes(x = eta2_focal, y = eta2_neighbor)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::geom_point(aes(color = eta2_interaction), size = 2, alpha = 0.85) +
  ggplot2::facet_wrap(~stage) +
  ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
  ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
  ggplot2::scale_color_gradientn(colours = c("#d9d9d6", "#a63446"),
                                 name = expression("Interaction " * eta^2)) +
  ggplot2::labs(x = expression("Own genotype " * eta^2), y = expression("Partner genotype " * eta^2)) +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 8),
                 axis.text = ggplot2::element_text(size = 5.5),
                 axis.title = ggplot2::element_text(size = 7),
                 legend.key.size = ggplot2::unit(0.35, "cm"),
                 legend.text = ggplot2::element_text(size = 6.5),
                 legend.title = ggplot2::element_text(size = 6.5))
save_panel(p11b, "FigS11b", 8.5, 7.0)

d11c <- rd("FigS11c") %>% dplyr::mutate(stage = factor(stage, levels = stage_levels))
p11c <- d11c %>%
  ggplot2::ggplot(aes(x = eta2_interaction, y = forcats::fct_reorder(gene, eta2_interaction))) +
  ggplot2::geom_col(fill = "#1e4e8c", width = 0.7) +
  ggplot2::facet_wrap(~stage, scales = "free_y") +
  ggplot2::labs(x = expression("Interaction " * eta^2 * " (focal x neighbor)"), y = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7, face = "italic"),
                 strip.text = ggplot2::element_text(face = "bold", size = 8))
save_panel(p11c, "FigS11c", 17.0, 7.5)

compose_vectorcraft(
  list(
    list(path = "figures/panels/FigS11a.pdf", rect = c(2, 8, 85, 70), label = "A", label_pos = c(2, 3.5)),
    list(path = "figures/panels/FigS11b.pdf", rect = c(89, 8, 85, 70), label = "B", label_pos = c(89, 3.5)),
    list(path = "figures/panels/FigS11c.pdf", rect = c(2, 83, 170, 75), label = "C", label_pos = c(2, 78.5))
  ),
  width_mm = 174, height_mm = 160, out_pdf = "figures/FigS11.pdf")

cat("\n[FigS11] sessionInfo()\n")
sessionInfo()
