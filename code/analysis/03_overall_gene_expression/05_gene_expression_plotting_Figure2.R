# ============================================================
# 03_overall_gene_expression / 05_gene_expression_plotting_Figure2.R
# 作用：绘制 Figure 2A–2E
#       2A 32 个组合（时期 × 纯种/混种组合）× 92 gene的整体expression热图（行内 0–1 缩放）
#       2B PLS-DA 得分图（Dim 1 × Dim 2，按 时期 × 种植模式 着色，含 0.6 置信椭圆）
#       2C Dim 1 的组间分布andtest / 2D Dim 2 的组间分布andtest
#       2E VIP 前 10 gene在 Dim 1–Dim 2 的载荷图（按pathway着色）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_overall_gene_expression/05_gene_expression_plotting_Figure2.R"
#
# 输入：<root>/result/03_overall_gene_expression/01_表_expression矩阵.xlsx（combo_matrix，minimal 列名修复）
#       <root>/result/03_overall_gene_expression/03_表_PLSDA.xlsx（scores、loadings、vip_top10、variance_explained、axis_tests）
# 输出：<root>/result/03_overall_gene_expression/05–09_Figure2A–2E_*.pdf / .png
# 其中 <root> = <repository_root>
#
# plotting约定：bioRtools::theme_prism()；write_data(units = "in", dpi = 500)；图内文字English；
#   P 值统一 format_p() 定点写法。热图用 pheatmap + ggplotify::as.ggplot（旧版用的 pheatmapColor 本机未安装，
#   热图内容and行内缩放一致，仅配色细节略有差异）。
# notes：旧版 Figure 2B/2E 的轴标签（Dim 1 26.00%、Dim 2 8.00%）and 2C/2D 的 P 值是硬编码文本，
#   本脚本一律改为按模型实算（见 03_表_PLSDA.xlsx）。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/03_overall_gene_expression/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
out_dir <- file.path(project_root, "result", "03_overall_gene_expression")

path_mat <- file.path(out_dir, "01_表_expression矩阵.xlsx")
path_plsda <- file.path(out_dir, "03_表_PLSDA.xlsx")

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

# ---- Figure 2A：整体expression热图 ----
# 实现notes：旧版用 pheatmapColor::pheatmap_bio()（本机未安装）；改用 pheatmap + ggplotify 转换时
#   会报 "'gpar' 成分 'fill' 的长度不能为零"（version兼容问题），故直接用 ggplot2::geom_tile 绘制热图：
#   内容一致（行内 0–1 缩放、gene按pathway排序、列按 时期 × 种植模式 分面），且不依赖热图包的渲染转换。
mat.heat <- as.matrix(mat.combo[, gene_cols])
rownames(mat.heat) <- mat.combo$combo

# 行内 0–1 缩放（等价于旧版的 bioRtools::scale01_rows；矩阵逐行处理需自己实现，
# 且某行完全无变异时返回全 0，避免出现 NaN 导致 pheatmap 报 gpar 长度为零）
scale01_by_row <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

mat.heat.scaled <- t(apply(mat.heat, 1, scale01_by_row))
stopifnot(!anyNA(mat.heat.scaled))

# gene排序：先按pathway（process），再按载荷 p1，便于同pathwaygene相邻
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
    stage = factor(stage, levels = c("tillering", "filling")),
    group = factor(group, levels = c("Pure", "Mix")),
    gene = factor(gene, levels = rev(gene_order)),
    combo_label = paste0(focal, "-", neighbor)
  ) -> df.heat.long

df.heat.long %>%
  ggplot2::ggplot(aes(x = combo_label, y = gene, fill = scaled_expression)) +
  ggplot2::geom_tile() +
  ggplot2::facet_grid(stage ~ group, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_gradientn(
    colours = c("#f2f2f0", "#a8bedb", "#1e4e8c"),
    name = "Scaled expression (0-1)"
  ) +
  ggplot2::labs(x = "Focal-neighbor combination", y = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = ggplot2::element_text(size = 5),
    strip.text = ggplot2::element_text(size = 10, face = "bold"),
    panel.grid = ggplot2::element_blank(),
    legend.position = "right"
  ) -> p.2a

bioRtools::write_data(p.2a, file = file.path(out_dir, "05_Figure2A_整体expression热图.pdf"),
                      width = 10, height = 12, units = "in")
bioRtools::write_data(p.2a, file = file.path(out_dir, "05_Figure2A_整体expression热图.png"),
                      width = 10, height = 12, units = "in", dpi = 500)
cat("[05] 输出：05_Figure2A_整体expression热图.pdf / .png\n")

# ---- 轴标签：按模型实算的 R2X ----
pc_lab <- function(i) {
  v <- variance_tab$R2X_single[i]
  paste0("Dim ", i, " (", formatC(100 * v, format = "f", digits = 2), "%)")
}

scores_tab %>%
  dplyr::mutate(
    stage = factor(stage, levels = c("tillering", "filling")),
    plotgroup = factor(paste0(stage, " ", group),
                       levels = c("tillering Pure", "tillering Mix", "filling Pure", "filling Mix"))
  ) -> scores_tab

# ---- Figure 2B：PLS-DA 得分图 ----
scores_tab %>%
  ggplot2::ggplot(aes(x = p1, y = p2, color = plotgroup)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::stat_ellipse(level = 0.6, linewidth = 0.5) +
  ggplot2::geom_point(size = 3.4) +
  bioRtools::scale_color_research() +
  ggplot2::labs(x = pc_lab(1), y = pc_lab(2), color = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "right") -> p.2b

bioRtools::write_data(p.2b, file = file.path(out_dir, "06_Figure2B_PLSDA得分图.pdf"),
                      width = 7, height = 5.5, units = "in")
bioRtools::write_data(p.2b, file = file.path(out_dir, "06_Figure2B_PLSDA得分图.png"),
                      width = 7, height = 5.5, units = "in", dpi = 500)
cat("[05] 输出：06_Figure2B_PLSDA得分图.pdf / .png\n")

# ---- Figure 2C / 2D：主轴分布andtest ----
plot_axis <- function(axis_name, pid, plabel) {
  axis_col <- dplyr::if_else(axis_name == "Dim 1", "p1", "p2")

  scores_tab %>%
    dplyr::mutate(value = .data[[axis_col]]) -> d

  stage_p <- axis_tests %>%
    dplyr::filter(axis == axis_name, comparison == "tillering vs filling", test == "wilcoxon") %>%
    dplyr::pull(p_value)

  group_p <- axis_tests %>%
    dplyr::filter(axis == axis_name, grepl("Mix vs Pure", comparison), test == "wilcoxon") %>%
    dplyr::mutate(
      stage = stringr::str_remove(comparison, ": Mix vs Pure"),
      stage = factor(stage, levels = c("tillering", "filling")),
      plotgroup = factor(paste0(stage, " ", "Mix"),
                         levels = c("tillering Pure", "tillering Mix", "filling Pure", "filling Mix"))
    )

  # 每个时期的 Mix 箱线中心（用于放 P 值标注）
  group_p %>%
    dplyr::mutate(y = -Inf) -> group_p

  d %>%
    ggplot2::ggplot(aes(x = plotgroup, y = value, fill = plotgroup)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.3) +
    ggplot2::geom_text(
      data = group_p,
      aes(x = plotgroup, y = y, label = format_p(p_value)),
      vjust = -0.8, size = 3.4, inherit.aes = FALSE
    ) +
    ggplot2::annotate(
      "text", x = 2.5, y = Inf,
      label = paste0("tillering vs filling: ", format_p(stage_p)),
      vjust = 1.6, size = 3.4
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.14, 0.16))) +
    ggplot2::scale_x_discrete(limits = c("tillering Pure", "tillering Mix", "filling Pure", "filling Mix")) +
    bioRtools::scale_fill_research() +
    ggplot2::labs(x = "Stage and planting mode", y = paste0(axis_name, " score")) +
    bioRtools::theme_prism() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1)
    ) -> p

  bioRtools::write_data(p, file = file.path(out_dir, paste0(pid, "_", plabel, ".pdf")),
                        width = 7, height = 6, units = "in")
  bioRtools::write_data(p, file = file.path(out_dir, paste0(pid, "_", plabel, ".png")),
                        width = 7, height = 6, units = "in", dpi = 500)
  cat("[05] 输出：", pid, "_", plabel, ".pdf / .png\n", sep = "")
}

plot_axis("Dim 1", "07", "Figure2C_Dim1分布")
plot_axis("Dim 2", "08", "Figure2D_Dim2分布")

# ---- Figure 2E：VIP 前 10 gene载荷图 ----
loadings_tab %>%
  dplyr::filter(gene %in% vip_top10$gene) %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene") %>%
  dplyr::left_join(vip_top10, by = "gene") -> load_top

stopifnot(nrow(load_top) == 10)

load_top %>%
  ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_segment(
    aes(x = 0, y = 0, xend = p1, yend = p2, color = process),
    arrow = ggplot2::arrow(length = ggplot2::unit(0.22, "cm"), type = "closed"), linewidth = 0.6
  ) +
  ggrepel::geom_label_repel(
    aes(x = p1, y = p2, label = gene, color = process),
    size = 3.4, fontface = "italic", max.overlaps = 30, show.legend = FALSE
  ) +
  bioRtools::scale_color_research() +
  ggplot2::labs(x = pc_lab(1), y = pc_lab(2), color = "Pathway") +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "right") -> p.2e

bioRtools::write_data(p.2e, file = file.path(out_dir, "09_Figure2E_VIP载荷图.pdf"),
                      width = 9, height = 6, units = "in")
bioRtools::write_data(p.2e, file = file.path(out_dir, "09_Figure2E_VIP载荷图.png"),
                      width = 9, height = 6, units = "in", dpi = 500)
cat("[05] 输出：09_Figure2E_VIP载荷图.pdf / .png\n")

cat("\n[05] sessionInfo()\n")
sessionInfo()
