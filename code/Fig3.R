# ============================================================
# Fig3.R —— Figure 3（Pathway-level responses: defence genes + defence score）
# Source Data：只读 source_data/Source_Data.xlsx（sheet：Fig3a、Fig3b）
# 运行：cd manuscript && Rscript code/Fig3.R
# 输出：figures/Fig3.pdf / .png；面板 figures/panels/Fig3a–b.pdf
# 软件：R 4.6.1；ggplot2；bioRtools::theme_prism()；set.seed(12345)
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/Fig3.R
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

# 3A：防御基因热图（combo × gene，值是每个组合内的 0-1 缩放表达）
d3a <- rd("Fig3a") %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels),
                process = factor(process, levels = c("Defense", "Jasmonic acid", "Salicylic acid")),
                gene = factor(gene, levels = rev(sort(unique(gene)))),
                combo_label = factor(combo_label, levels = unique(combo_label[order(group, focal, neighbor)])))
p3a <- d3a %>%
  ggplot2::ggplot(aes(x = combo_label, y = gene, fill = scaled_expression)) +
  ggplot2::geom_tile() +
  ggplot2::facet_grid(process ~ stage, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_gradientn(colours = c("#f2f2f0", "#a8bedb", "#1e4e8c"),
                                name = "Scaled median expression (0-1)") +
  ggplot2::labs(x = "Focal-neighbour combination", y = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
                 axis.text.y = ggplot2::element_text(size = 5.2, face = "italic"),
                 # 2026-10-19：轴标题左端曾越出纸边 2.6 pt（"Fo" 被切）→ 字号 7.5 → 7 + 左/下边距
                 axis.title.x = ggplot2::element_text(size = 7),
                 plot.margin = ggplot2::margin(2, 2, 3, 8),
                 strip.text = ggplot2::element_text(size = 6, face = "bold"),
                 panel.grid = ggplot2::element_blank(),
                 legend.key.width = ggplot2::unit(0.3, "cm"),
                 legend.title = ggplot2::element_text(size = 6, angle = 90),
                 legend.text = ggplot2::element_text(size = 5.5))
save_panel(p3a, "Fig3a", 17.0, 14.5)

# 3B：每样本防御分数（Mix vs Pure，三时期）
d3b <- rd("Fig3b") %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels),
                group = factor(group, levels = c("Pure", "Mix")))
pval_lab <- d3b %>% dplyr::distinct(stage, p_value) %>%
  dplyr::mutate(lab = dplyr::if_else(p_value < 0.001, "P<0.001",
                                     paste0("P=", formatC(p_value, format = "f", digits = 3))))
p3b <- d3b %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(y_pos = max(defence_score, na.rm = TRUE) + 0.06 * diff(range(defence_score, na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot(aes(x = group, y = defence_score, fill = group)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  ggplot2::geom_text(data = pval_lab, aes(x = 1.5, y = Inf, label = lab),
                     inherit.aes = FALSE, vjust = 1.6, size = 2.6) +
  ggplot2::facet_wrap(~stage) +
  ggplot2::scale_x_discrete(limits = c("Pure", "Mix")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.14))) +
  bioRtools::scale_fill_research() +
  # 2026-10-19：旋转 y 轴标题两头被切 → 缩短 + 缩小字号 + 上下留边距
  ggplot2::labs(x = NULL, y = "Defence pathway score (Z)") +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(face = "bold", size = 8),
                 axis.title.y = ggplot2::element_text(size = 6.5),
                 plot.margin = ggplot2::margin(3, 2, 3, 2))
save_panel(p3b, "Fig3b", 17.0, 5.2)

compose_vectorcraft(
  list(
    list(path = "figures/panels/Fig3a.pdf", rect = c(2, 8, 170, 145), label = "A", label_pos = c(2, 3.5)),
    list(path = "figures/panels/Fig3b.pdf", rect = c(2, 158, 170, 52), label = "B", label_pos = c(2, 153.5))
  ),
  width_mm = 174, height_mm = 212, out_pdf = "figures/Fig3.pdf")

cat("\n[Fig3] sessionInfo()\n")
sessionInfo()
