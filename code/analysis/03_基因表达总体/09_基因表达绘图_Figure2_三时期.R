# ============================================================
# 03_基因表达总体 / 09_基因表达绘图_Figure2_三时期.R
# 作用：绘制**三时期**版 Figure 2A–2E
#       2A 48 个组合（3 时期 × 16 组合）× 92 基因的整体表达热图（行内 0–1 缩放）
#       2B PLS-DA 得分图（Dim 1 × Dim 2，按 时期 × 种植模式 着色，含 0.6 置信椭圆）
#       2C Dim 1 的组间分布与检验 / 2D Dim 2 的组间分布与检验（三时期）
#       2E VIP 前 10 基因在 Dim 1–Dim 2 的载荷图（按通路着色）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_基因表达总体/09_基因表达绘图_Figure2_三时期.R"
#
# 输入：<root>/result/03_基因表达总体/11_表_表达矩阵_三时期.xlsx（combo_matrix，minimal 列名修复）
#       <root>/result/03_基因表达总体/13_表_PLSDA_三时期.xlsx（scores、loadings、vip_top10、variance_explained、axis_tests）
# 输出：<root>/result/03_基因表达总体/14–18_Figure2A–2E_*_三时期.pdf / .png
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 绘图约定：bioRtools::theme_prism()、bioRtools::scale_*_research()、write_data(units = "in", dpi = 500)；
#   图内文字英文；P 值统一定点写法。
# 与两时期版本的差异：① 时期由 2 变 3（热图 6 个分面、得分/分布图 6 组）；② 主轴检验的时期整体比较
#   改为单因素检验（one-way ANOVA / Kruskal–Wallis，见 13 表 axis_tests 的 test_group 列），
#   分布图标注改为"整体检验 P + 各期 Mix vs Pure P"。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
out_dir <- file.path(project_root, "result", "03_基因表达总体")

path_mat <- file.path(out_dir, "11_表_表达矩阵_三时期.xlsx")
path_plsda <- file.path(out_dir, "13_表_PLSDA_三时期.xlsx")

stopifnot(file.exists(path_mat), file.exists(path_plsda))

format_p <- function(p, digits = 3) {
  dplyr::if_else(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, format = "f", digits = digits)))
}

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_mat, "combo_matrix") -> mat.combo
bioRtools::read_data(path_plsda, sheet = "scores") -> scores_tab
bioRtools::read_data(path_plsda, sheet = "loadings") -> loadings_tab
bioRtools::read_data(path_plsda, sheet = "vip_top10") -> vip_top10
bioRtools::read_data(path_plsda, sheet = "variance_explained") -> variance_tab
bioRtools::read_data(path_plsda, sheet = "axis_tests") -> axis_tests
read_sheet_raw(path_mat, "gene_annotation") -> df.gene

meta_cols <- c("combo", "stage", "focal", "neighbor", "group")
gene_cols <- setdiff(colnames(mat.combo), meta_cols)

stage_levels <- c("tillering", "flowering", "filling")
plotgroup_levels <- as.vector(outer(stage_levels, c("Pure", "Mix"), paste))

# ---- Figure 2A：整体表达热图 ----
mat.heat <- as.matrix(mat.combo[, gene_cols])
rownames(mat.heat) <- mat.combo$combo

scale01_by_row <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

mat.heat.scaled <- t(apply(mat.heat, 1, scale01_by_row))
stopifnot(!anyNA(mat.heat.scaled))

gene_order <- loadings_tab %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene") %>%
  dplyr::arrange(process, p1) %>%
  dplyr::pull(gene) %>%
  unique()

df.gene %>%
  dplyr::mutate(process = factor(process, levels = sort(unique(process)))) -> df.gene

as.data.frame(mat.heat.scaled) %>%
  tibble::rownames_to_column(var = "combo") %>%
  tidyr::pivot_longer(cols = -combo, names_to = "gene", values_to = "scaled_expression") %>%
  dplyr::left_join(mat.combo[, meta_cols], by = "combo") %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene") %>%
  dplyr::mutate(
    stage = factor(stage, levels = stage_levels),
    group = factor(group, levels = c("Pure", "Mix")),
    gene = factor(gene, levels = rev(gene_order)),
    combo_label = paste0(focal, "-", neighbor)
  ) -> df.heat.long

df.heat.long %>%
  ggplot2::ggplot(aes(x = combo_label, y = gene, fill = scaled_expression)) +
  ggplot2::geom_tile() +
  ggplot2::facet_grid(stage ~ group, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_gradientn(colours = c("#f2f2f0", "#a8bedb", "#1e4e8c"),
                                name = "Scaled expression (0-1)") +
  ggplot2::labs(x = "Focal-neighbour combination", y = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = ggplot2::element_text(size = 5),
    strip.text = ggplot2::element_text(size = 9, face = "bold"),
    panel.grid = ggplot2::element_blank(),
    legend.position = "right"
  ) -> p.2a

bioRtools::write_data(p.2a, file = file.path(out_dir, "14_Figure2A_整体表达热图_三时期.pdf"),
                      width = 11, height = 13, units = "in")
bioRtools::write_data(p.2a, file = file.path(out_dir, "14_Figure2A_整体表达热图_三时期.png"),
                      width = 11, height = 13, units = "in", dpi = 500)
cat("[09] 输出：14_Figure2A_整体表达热图_三时期.pdf / .png\n")

# ---- 轴标签（实算 R2X）----
pc_lab <- function(i) {
  v <- variance_tab$R2X_single[i]
  paste0("Dim ", i, " (", formatC(100 * v, format = "f", digits = 2), "%)")
}

scores_tab %>%
  dplyr::mutate(
    stage = factor(stage, levels = stage_levels),
    plotgroup = factor(paste0(stage, " ", group), levels = plotgroup_levels)
  ) -> scores_tab

# ---- Figure 2B：PLS-DA 得分图 ----
scores_tab %>%
  ggplot2::ggplot(aes(x = p1, y = p2, color = plotgroup)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
# 椭圆用 type = "norm"（正态近似）而非默认的 MASS::cov.trob（稳健估计）：
#   三时期下 Pure 组每个时期只有 4 个组合均值点，稳健估计会报 "Probable convergence failure"。
  ggplot2::stat_ellipse(level = 0.6, type = "norm", linewidth = 0.5) +
  ggplot2::geom_point(size = 3.2) +
  bioRtools::scale_color_research() +
  ggplot2::labs(x = pc_lab(1), y = pc_lab(2), color = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "right") -> p.2b

bioRtools::write_data(p.2b, file = file.path(out_dir, "15_Figure2B_PLSDA得分图_三时期.pdf"),
                      width = 7.5, height = 5.5, units = "in")
bioRtools::write_data(p.2b, file = file.path(out_dir, "15_Figure2B_PLSDA得分图_三时期.png"),
                      width = 7.5, height = 5.5, units = "in", dpi = 500)
cat("[09] 输出：15_Figure2B_PLSDA得分图_三时期.pdf / .png\n")

# ---- Figure 2C / 2D：主轴分布与检验（三时期）----
plot_axis <- function(axis_name, pid, plabel) {
  axis_col <- dplyr::if_else(axis_name == "Dim 1", "p1", "p2")

  scores_tab %>%
    dplyr::mutate(value = .data[[axis_col]]) -> d

  overall_p <- axis_tests %>%
    dplyr::filter(axis == axis_name, test_group == "stage_overall", test == "kruskal_wallis") %>%
    dplyr::pull(p_value)

  group_p <- axis_tests %>%
    dplyr::filter(axis == axis_name, test_group == "within_stage", test == "wilcoxon") %>%
    dplyr::mutate(
      stage = stringr::str_remove(comparison, ": Mix vs Pure"),
      stage = factor(stage, levels = stage_levels),
      plotgroup = factor(paste0(stage, " Mix"), levels = plotgroup_levels),
      y = -Inf
    )

  d %>%
    ggplot2::ggplot(aes(x = plotgroup, y = value, fill = plotgroup)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.55, size = 1.2) +
    ggplot2::geom_text(
      data = group_p,
      aes(x = plotgroup, y = y, label = format_p(p_value)),
      vjust = -0.8, size = 3.1, inherit.aes = FALSE
    ) +
    ggplot2::annotate("text", x = 3.5, y = Inf,
                      label = paste0("stage (Kruskal-Wallis): ", format_p(overall_p)),
                      vjust = 1.6, size = 3.1) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.14, 0.16))) +
    ggplot2::scale_x_discrete(limits = plotgroup_levels) +
    bioRtools::scale_fill_research() +
    ggplot2::labs(x = "Stage and planting mode", y = paste0(axis_name, " score")) +
    bioRtools::theme_prism() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1)
    ) -> p

  bioRtools::write_data(p, file = file.path(out_dir, paste0(pid, "_", plabel, ".pdf")),
                        width = 8, height = 6, units = "in")
  bioRtools::write_data(p, file = file.path(out_dir, paste0(pid, "_", plabel, ".png")),
                        width = 8, height = 6, units = "in", dpi = 500)
  cat("[09] 输出：", pid, "_", plabel, ".pdf / .png\n", sep = "")
}

plot_axis("Dim 1", "16", "Figure2C_Dim1分布_三时期")
plot_axis("Dim 2", "17", "Figure2D_Dim2分布_三时期")

# ---- Figure 2E：VIP 前 10 基因载荷图 ----
loadings_tab %>%
  dplyr::filter(gene %in% vip_top10$gene) %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene") %>%
  dplyr::left_join(vip_top10, by = "gene") -> load_top

stopifnot(nrow(load_top) == 10)

load_top %>%
  ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_segment(aes(x = 0, y = 0, xend = p1, yend = p2, color = process),
                        arrow = ggplot2::arrow(length = ggplot2::unit(0.22, "cm"), type = "closed"),
                        linewidth = 0.6) +
  ggrepel::geom_label_repel(aes(x = p1, y = p2, label = gene, color = process),
                            size = 3.4, fontface = "italic", max.overlaps = 30,
                            show.legend = FALSE) +
  bioRtools::scale_color_research() +
  ggplot2::labs(x = pc_lab(1), y = pc_lab(2), color = "Pathway") +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "right") -> p.2e

bioRtools::write_data(p.2e, file = file.path(out_dir, "18_Figure2E_VIP载荷图_三时期.pdf"),
                      width = 9, height = 6, units = "in")
bioRtools::write_data(p.2e, file = file.path(out_dir, "18_Figure2E_VIP载荷图_三时期.png"),
                      width = 9, height = 6, units = "in", dpi = 500)
cat("[09] 输出：18_Figure2E_VIP载荷图_三时期.pdf / .png\n")

cat("\n[09] sessionInfo()\n")
sessionInfo()
