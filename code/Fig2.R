# ============================================================
# Fig2.R —— Figure 2（Transcriptional landscape of the 92-gene panel）
# Source Data：只读 source_data/Source_Data.xlsx（sheet：Fig2a–Fig2e）
# 运行：cd manuscript && Rscript code/Fig2.R
# 输出：figures/Fig2.pdf / .png；面板 figures/panels/Fig2a–e.pdf
# 软件：R 4.6.1；ggplot2；bioRtools::theme_prism()；set.seed(12345)
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/Fig2.R
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

d2a <- rd("Fig2a") %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels),
                pathway = factor(pathway, levels = rev(sort(unique(pathway)))),
                combo_label = factor(combo_label, levels = unique(combo_label[order(group, focal, neighbor)])))
p2a <- d2a %>%
  ggplot2::ggplot(aes(x = combo_label, y = pathway, fill = mean_z)) +
  ggplot2::geom_tile() +
  ggplot2::facet_grid(. ~ stage, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_gradient2(low = "#3B1F4E", mid = "#FFFFFF", high = "#780018",
                                midpoint = 0, name = "Scaled expression (Z)") +
  ggplot2::labs(x = "Focal-neighbour combination", y = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
                 axis.text.y = ggplot2::element_text(size = 5.5),
                 axis.title.x = ggplot2::element_text(size = 7.5),
                 strip.text = ggplot2::element_text(size = 8, face = "bold"),
                 panel.grid = ggplot2::element_blank(),
                 legend.key.width = ggplot2::unit(0.3, "cm"),
                 legend.key.height = ggplot2::unit(1.0, "cm"),
                 legend.title = ggplot2::element_text(size = 5.5, angle = 0),
                 legend.text = ggplot2::element_text(size = 6))
save_panel(p2a, "Fig2a", 17.0, 7.5)

d2b <- rd("Fig2b") %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels),
                plotgroup = factor(paste(stage, group), levels = as.vector(outer(stage_levels, c("Pure", "Mix"), paste))))
p2b <- d2b %>%
  ggplot2::ggplot(aes(x = p1, y = p2, color = plotgroup)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::stat_ellipse(level = 0.6, type = "norm", linewidth = 0.5) +
  ggplot2::geom_point(size = 3.2) +
  bioRtools::scale_color_research() +
  # 轴标签取自 PLS-DA 的 R2X（单成分解释率）：p1 = 29.4%、p2 = 10.8%、p3 = 19.1%
  # （2026-10-19 修：原写 19.1% / 19.0% 是两时期旧模型的残留数值）
  ggplot2::labs(x = "Dim 1 (29.4%)", y = "Dim 2 (10.8%)", color = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "right", legend.key.size = ggplot2::unit(0.3, "cm"),
                 legend.text = ggplot2::element_text(size = 6), legend.margin = ggplot2::margin(0, 0, 0, 0))
save_panel(p2b, "Fig2b", 8.3, 5.0)

panel_axis <- function(sheet, ylab, name) {
  d <- rd(sheet) %>%
    dplyr::mutate(stage = factor(stage, levels = stage_levels),
                  plotgroup = factor(paste(stage, group), levels = as.vector(outer(stage_levels, c("Pure", "Mix"), paste))))
  grp_p <- d %>% dplyr::distinct(stage, plotgroup, wilcoxon_p) %>%
    dplyr::mutate(y = -Inf, lab = dplyr::if_else(wilcoxon_p < 0.001, "P<0.001",
                                                 paste0("P=", formatC(wilcoxon_p, format = "f", digits = 2))))
  kruskal_p <- d$kruskal_p[1]
  p <- d %>%
    ggplot2::ggplot(aes(x = plotgroup, y = score, fill = plotgroup)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.55, size = 1.2) +
    ggplot2::geom_text(data = grp_p, aes(x = plotgroup, y = y, label = lab),
                       vjust = -0.8, size = 2.6, inherit.aes = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.14, 0.16))) +
    bioRtools::scale_fill_research() +
    ggplot2::labs(x = NULL, y = ylab,
                  subtitle = paste0("Stage effect (Kruskal-Wallis): P=", formatC(kruskal_p, format = "f", digits = 2))) +
    bioRtools::theme_prism() +
    ggplot2::theme(legend.position = "none", plot.subtitle = ggplot2::element_text(size = 7, hjust = 0),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1, size = 7))
  save_panel(p, name, 8.3, 5.0)
}
panel_axis("Fig2c", "Dim 1 score", "Fig2c")
panel_axis("Fig2d", "Dim 2 score", "Fig2d")

d2e <- rd("Fig2e")
p2e <- d2e %>%
  ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_segment(aes(x = 0, y = 0, xend = p1, yend = p2, color = process),
                        arrow = ggplot2::arrow(length = ggplot2::unit(0.22, "cm"), type = "closed"),
                        linewidth = 0.6) +
  ggrepel::geom_label_repel(aes(x = p1, y = p2, label = gene, color = process),
                            size = 2.0, fontface = "italic", show.legend = FALSE,
                            label.size = 0.12, label.padding = ggplot2::unit(0.06, "lines"),
                            box.padding = 0.45, min.segment.length = 0.15, force = 12,
                            max.overlaps = 20, seed = 12345) +
  bioRtools::scale_color_research() +
  ggplot2::labs(x = "Dim 1 (29.4%)", y = "Dim 2 (10.8%)", color = "Pathway") +
  # 2026-10-19：与投稿面板一致——图例由右侧改为置底 3 行（右侧图例曾把 83 mm 面板的绘图区
  #   压到约 47 mm 宽），并关掉面板 clip 防止标签被裁
  # 2026-10-19：2 列 4 行 + 去掉图例标题（3 列时图例宽过面板，最右列与标题被纸边切掉）
  ggplot2::guides(color = ggplot2::guide_legend(nrow = 4, ncol = 2, byrow = TRUE, title = NULL)) +
  ggplot2::coord_cartesian(clip = "off") +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "bottom", legend.key.size = ggplot2::unit(0.26, "cm"),
                 legend.text = ggplot2::element_text(size = 5.5),
                 legend.title = ggplot2::element_text(size = 6),
                 legend.box.spacing = ggplot2::unit(0.08, "cm"),
                 legend.margin = ggplot2::margin(0, 0, 0, 0))
save_panel(p2e, "Fig2e", 8.3, 9.2)

compose_vectorcraft(
  list(
    list(path = "figures/panels/Fig2a.pdf", rect = c(2, 8, 170, 75), label = "A", label_pos = c(2, 3.5)),
    list(path = "figures/panels/Fig2b.pdf", rect = c(2, 88, 83, 50), label = "B", label_pos = c(2, 83.5)),
    list(path = "figures/panels/Fig2c.pdf", rect = c(91, 88, 83, 50), label = "C", label_pos = c(91, 83.5)),
    list(path = "figures/panels/Fig2d.pdf", rect = c(2, 143, 83, 50), label = "D", label_pos = c(2, 138.5)),
    list(path = "figures/panels/Fig2e.pdf", rect = c(91, 143, 83, 72), label = "E", label_pos = c(91, 138.5))
  ),
  width_mm = 174, height_mm = 237, out_pdf = "figures/Fig2.pdf")

cat("\n[Fig2] sessionInfo()\n")
sessionInfo()
