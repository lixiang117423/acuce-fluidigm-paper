# ============================================================
# FigS8.R —— Figure S8（Pot-level clustering and design decomposition；2 个子图 a、b）
# Source Data：只读 source_data/Source_Data.xlsx（sheet：FigS8a、FigS8b）
# 运行：cd manuscript && Rscript code/FigS8.R
# 输出：figures/FigS8.pdf / .png；面板 figures/panels/FigS8a–b.pdf
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/FigS8.R
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
fmt_p <- function(p, digits = 2) dplyr::if_else(p < 0.001, "P<0.001", paste0("P=", formatC(p, format = "f", digits = digits)))

d8a <- rd("FigS8a") %>%
  dplyr::select(disease, p_leaf, p_pot_random, p_pot_agg) %>%
  tidyr::pivot_longer(cols = -disease, names_to = "model", values_to = "p_value") %>%
  dplyr::mutate(model = factor(model, levels = c("p_leaf", "p_pot_random", "p_pot_agg"),
                               labels = c("Per plant", "+ pot RE", "Pot-level")),
                disease = factor(disease, levels = c("blast", "brown spot"),
                                 labels = c("Rice blast", "Brown spot")),
                lab = fmt_p(p_value))
p8a <- d8a %>%
  ggplot2::ggplot(aes(x = model, y = p_value, fill = disease)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.65) +
  ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
  ggplot2::geom_text(aes(label = lab), position = ggplot2::position_dodge(width = 0.7),
                     vjust = -0.4, size = 2.5) +
  bioRtools::scale_fill_research() +
  # 2026-10-19：旋转 y 轴标题顶端被切 → 缩短 + 缩小字号 + 留边距
  ggplot2::labs(x = NULL, y = "P value (Pure vs Mix)", fill = "Disease") +
  bioRtools::theme_prism() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7, angle = 30, hjust = 1, vjust = 1),
                 axis.title.y = ggplot2::element_text(size = 6.5),
                 plot.margin = ggplot2::margin(3, 3, 3, 2),
                 legend.position = "bottom", legend.key.size = ggplot2::unit(0.3, "cm"),
                 legend.text = ggplot2::element_text(size = 7))
save_panel(p8a, "FigS8a", 8.3, 7.0)

d8b <- rd("FigS8b") %>% dplyr::mutate(disease = factor(disease, levels = c("blast", "brown spot"),
                                                       labels = c("Rice blast", "Brown spot")))
p8b <- d8b %>%
  ggplot2::ggplot(aes(x = effect_focal, y = effect_neighbor)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey45") +
  ggplot2::geom_point(size = 2.4, color = "#1e4e8c") +
  ggrepel::geom_text_repel(aes(label = pair), size = 2.4, max.overlaps = 20, seed = 12345) +
  ggplot2::facet_wrap(~disease, scales = "free") +
  ggplot2::labs(x = "Effect on focal (Z)", y = "Effect on the partner (Z)") +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 8, face = "bold"),
                 axis.title = ggplot2::element_text(size = 6.5),
                 plot.margin = ggplot2::margin(3, 4, 3, 2))
save_panel(p8b, "FigS8b", 8.7, 7.0)

compose_vectorcraft(
  list(
    list(path = "figures/panels/FigS8a.pdf", rect = c(2, 8, 83, 70), label = "A", label_pos = c(2, 3.5)),
    list(path = "figures/panels/FigS8b.pdf", rect = c(87, 8, 87, 70), label = "B", label_pos = c(87, 3.5))
  ),
  width_mm = 174, height_mm = 80, out_pdf = "figures/FigS8.pdf")

cat("\n[FigS8] sessionInfo()\n")
sessionInfo()
