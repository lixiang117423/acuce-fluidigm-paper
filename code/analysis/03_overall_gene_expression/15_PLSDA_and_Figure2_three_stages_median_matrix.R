# ============================================================
# 03_overall_gene_expression / 15_PLSDA_and_Figure2_three_stages_median_matrix.R
# 作用：按 2026-10-18 决定的**图侧口径**（组合均值改用**中位数**、不剔除任何data）重做：
#       ① PLS-DA（组合均值矩阵）② Figure 2A–2E
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_overall_gene_expression/15_PLSDA_and_Figure2_three_stages_median_matrix.R"
#
# 输入：<root>/result/03_overall_gene_expression/23_表_expression矩阵_three_stages_without_outlier_removal.xlsx（combo_matrix_median、gene_annotation）
# 输出：<root>/result/03_overall_gene_expression/24_表_PLSDA_three_stages_median_matrix.xlsx
#       <root>/result/03_overall_gene_expression/25–29_Figure2A–2E_*_median_matrix.pdf / .png
# 其中 <root> = <repository_root>
#
# 为什么要换：`11_表` 的组合均值是"剔除离群值"后算的（29.55% 的格子and不剔除version不同、个别格子差 4 个数量级），
#   andsample级statistics（limma，未剔除）口径不一致；改用"不剔除 + 中位数"后图and_statistics同源，且热图不被极端sample拉成二值图。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
out_dir <- file.path(project_root, "result", "03_overall_gene_expression")
path_mat <- file.path(out_dir, "23_表_expression矩阵_three_stages_without_outlier_removal.xlsx")
path_out <- file.path(out_dir, "24_表_PLSDA_three_stages_median_matrix.xlsx")

stopifnot(file.exists(path_mat))

stage_levels <- c("tillering", "flowering", "filling")
plotgroup_levels <- as.vector(outer(stage_levels, c("Pure", "Mix"), paste))
sig_levels <- c("FDR < 0.05", "P < 0.05", "n.s.")
sig_colors <- c("FDR < 0.05" = "#a63446", "P < 0.05" = "#e0a458", "n.s." = "grey62")

format_p <- function(p, digits = 3) {
  dplyr::if_else(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, format = "f", digits = digits)))
}

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_mat, "combo_matrix_median") -> mat.combo
bioRtools::read_data(path_mat, sheet = "gene_annotation") -> df.gene

meta_cols <- c("combo", "stage", "focal", "neighbor", "group")
gene_cols <- setdiff(colnames(mat.combo), meta_cols)
stopifnot(length(gene_cols) == 92, nrow(mat.combo) == 48)

mat.combo %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels)) -> mat.combo

# ---------------- PLS-DA ----------------
x <- as.matrix(mat.combo[, gene_cols])
rownames(x) <- mat.combo$combo
y <- as.factor(mat.combo$group)

ropls::opls(x, y, predI = 3) -> model
stopifnot(length(ropls::getScoreMN(model)) > 0)

as.data.frame(model@summaryDF) %>% tibble::rownames_to_column(var = "model_type") -> cv_tab
cat("[15] PLS-DA（中位数组合矩阵）模型摘要：\n")
print(as.data.frame(cv_tab))

ropls::getScoreMN(model) %>%
  as.data.frame() %>% tibble::rownames_to_column(var = "combo") %>%
  dplyr::left_join(mat.combo[, meta_cols], by = "combo") -> scores_tab

ropls::getLoadingMN(model) %>%
  as.data.frame() %>% tibble::rownames_to_column(var = "gene") %>%
  dplyr::select(gene, p1, p2, p3) -> loadings_tab

ropls::getVipVn(model) %>%
  tibble::enframe(name = "gene", value = "vip") %>%
  dplyr::arrange(dplyr::desc(vip)) -> vip_tab
vip_top10 <- dplyr::slice_head(vip_tab, n = 10)
cat("\n[15] VIP 前 10：\n"); print(as.data.frame(vip_top10))

as.data.frame(model@modelDF) %>%
  tibble::rownames_to_column(var = "component") %>%
  dplyr::mutate(R2X_single = R2X, R2Y_single = R2Y) -> variance_tab

score_test <- function(axis_col, label) {
  d <- scores_tab; d$axis_value <- d[[axis_col]]
  aov_p <- summary(stats::aov(axis_value ~ stage, data = d))[[1]][["Pr(>F)"]][1]
  kw <- stats::kruskal.test(axis_value ~ stage, data = d)
  pair_w <- stats::pairwise.wilcox.test(d$axis_value, d$stage, p.adjust.method = "BH")
  pairs_tab <- as.data.frame(as.table(pair_w$p.value)) %>%
    dplyr::filter(!is.na(Freq)) %>%
    dplyr::transmute(contrast = paste(Var1, "vs", Var2), p_value = Freq)
  by_stage <- d %>% dplyr::group_by(stage) %>% dplyr::group_modify(~ {
    tt <- stats::t.test(.x$axis_value ~ .x$group)
    ww <- suppressWarnings(stats::wilcox.test(.x$axis_value ~ .x$group))
    dd <- (mean(.x$axis_value[.x$group == "Mix"]) - mean(.x$axis_value[.x$group == "Pure"])) /
      sqrt((stats::var(.x$axis_value[.x$group == "Pure"]) + stats::var(.x$axis_value[.x$group == "Mix"])) / 2)
    tibble::tibble(n_pure = sum(.x$group == "Pure"), n_mix = sum(.x$group == "Mix"),
                   diff_mix_minus_pure = mean(.x$axis_value[.x$group == "Mix"]) - mean(.x$axis_value[.x$group == "Pure"]),
                   t_p_value = tt$p.value, wilcoxon_p_value = ww$p.value, cohens_d = dd)
  }) %>% dplyr::ungroup()
  dplyr::bind_rows(
    tibble::tibble(axis = label, test_group = "stage_overall", comparison = "all three stages",
                   test = "one_way_ANOVA", estimate = NA_real_, p_value = aov_p, statistic = NA_real_, n = nrow(d)),
    tibble::tibble(axis = label, test_group = "stage_overall", comparison = "all three stages",
                   test = "kruskal_wallis", estimate = NA_real_, p_value = kw$p.value,
                   statistic = unname(kw$statistic), n = nrow(d)),
    pairs_tab %>% dplyr::mutate(axis = label, test_group = "stage_pairwise", test = "wilcoxon_BH",
                                estimate = NA_real_, statistic = NA_real_, n = nrow(d)) %>%
      dplyr::select(axis, test_group, comparison = contrast, test, estimate, p_value, statistic, n),
    by_stage %>%
      tidyr::pivot_longer(cols = c(t_p_value, wilcoxon_p_value), names_to = "test", values_to = "p_value") %>%
      dplyr::mutate(axis = label, test_group = "within_stage",
                    comparison = paste0(stage, ": Mix vs Pure"),
                    test = dplyr::if_else(test == "t_p_value", "welch_t", "wilcoxon")) %>%
      dplyr::select(axis, test_group, comparison, test, estimate = diff_mix_minus_pure, p_value, cohens_d, n_pure, n_mix)
  )
}

test_tab <- dplyr::bind_rows(score_test("p1", "Dim 1"), score_test("p2", "Dim 2"))

bioRtools::write_data(
  list(scores = scores_tab, loadings = loadings_tab, vip = vip_tab, vip_top10 = vip_top10,
       variance_explained = variance_tab, axis_tests = test_tab, cv_summary = cv_tab),
  file = path_out
)
cat("\n[15] 输出：", path_out, "\n")

# ---------------- Figure 2A–2E ----------------
save_panel <- function(p, pid, label, width_in, height_in) {
  bioRtools::write_data(p, file = file.path(out_dir, paste0(pid, "_", label, ".pdf")),
                        width = width_in, height = height_in, units = "in")
  bioRtools::write_data(p, file = file.path(out_dir, paste0(pid, "_", label, ".png")),
                        width = width_in, height = height_in, units = "in", dpi = 500)
  cat("[15] 输出：", pid, "_", label, ".pdf / .png\n", sep = "")
}

# 2A 热图（行内 0–1 缩放）
mat.heat <- as.matrix(mat.combo[, gene_cols])
rownames(mat.heat) <- mat.combo$combo
scale01_by_row <- function(v) {
  rng <- range(v, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0, length(v)))
  (v - rng[1]) / diff(rng)
}
mat.heat.scaled <- t(apply(mat.heat, 1, scale01_by_row))
stopifnot(!anyNA(mat.heat.scaled))

gene_order <- loadings_tab %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene") %>%
  dplyr::arrange(process, p1) %>% dplyr::pull(gene) %>% unique()

as.data.frame(mat.heat.scaled) %>%
  tibble::rownames_to_column(var = "combo") %>%
  tidyr::pivot_longer(cols = -combo, names_to = "gene", values_to = "scaled_expression") %>%
  dplyr::left_join(mat.combo[, meta_cols], by = "combo") %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels),
                group = factor(group, levels = c("Pure", "Mix")),
                gene = factor(gene, levels = rev(gene_order)),
                combo_label = paste0(focal, "-", neighbor)) %>%
  ggplot2::ggplot(aes(x = combo_label, y = gene, fill = scaled_expression)) +
  ggplot2::geom_tile() +
  ggplot2::facet_grid(stage ~ group, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_gradientn(colours = c("#f2f2f0", "#a8bedb", "#1e4e8c"),
                                name = "Scaled median expression (0-1)") +
  ggplot2::labs(x = "Focal-neighbour combination", y = NULL,
                caption = "Combo values are medians across replicate plants (no outlier removal); each row scaled to 0-1") +
  bioRtools::theme_prism() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
                 axis.text.y = ggplot2::element_text(size = 5),
                 strip.text = ggplot2::element_text(size = 9, face = "bold"),
                 panel.grid = ggplot2::element_blank(),
                 legend.position = "right",
                 plot.caption = ggplot2::element_text(size = 7, hjust = 0)) -> p.2a
save_panel(p.2a, "25", "Figure2A_整体expression热图_three_stages_median_matrix", 11, 13)

pc_lab <- function(i) {
  v <- variance_tab$R2X_single[i]
  paste0("Dim ", i, " (", formatC(100 * v, format = "f", digits = 2), "%)")
}

scores_tab %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels),
                plotgroup = factor(paste0(stage, " ", group), levels = plotgroup_levels)) -> scores_plot

scores_plot %>%
  ggplot2::ggplot(aes(x = p1, y = p2, color = plotgroup)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::stat_ellipse(level = 0.6, type = "norm", linewidth = 0.5) +
  ggplot2::geom_point(size = 3.2) +
  bioRtools::scale_color_research() +
  ggplot2::labs(x = pc_lab(1), y = pc_lab(2), color = NULL,
                caption = paste0("Q2(cum) = ", formatC(cv_tab$`Q2(cum)`[1], format = "f", digits = 3),
                                 "; permutation pQ2 = ", formatC(cv_tab$pQ2[1], format = "f", digits = 2),
                                 " -> exploratory only")) +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "right", plot.caption = ggplot2::element_text(size = 7, hjust = 0)) -> p.2b
save_panel(p.2b, "26", "Figure2B_PLSDA得分图_three_stages_median_matrix", 7.5, 5.5)

plot_axis <- function(axis_name, pid, label) {
  axis_col <- dplyr::if_else(axis_name == "Dim 1", "p1", "p2")
  d <- scores_plot %>% dplyr::mutate(value = .data[[axis_col]])
  overall_p <- test_tab %>% dplyr::filter(axis == axis_name, test_group == "stage_overall",
                                          test == "kruskal_wallis") %>% dplyr::pull(p_value)
  group_p <- test_tab %>%
    dplyr::filter(axis == axis_name, test_group == "within_stage", test == "wilcoxon") %>%
    dplyr::mutate(stage = stringr::str_remove(comparison, ": Mix vs Pure"),
                  stage = factor(stage, levels = stage_levels),
                  plotgroup = factor(paste0(stage, " Mix"), levels = plotgroup_levels), y = -Inf)
  d %>%
    ggplot2::ggplot(aes(x = plotgroup, y = value, fill = plotgroup)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.55, size = 1.2) +
    ggplot2::geom_text(data = group_p, aes(x = plotgroup, y = y, label = format_p(p_value)),
                       vjust = -0.8, size = 3.1, inherit.aes = FALSE) +
    ggplot2::annotate("text", x = 3.5, y = Inf,
                      label = paste0("stage (Kruskal-Wallis): ", format_p(overall_p)),
                      vjust = 1.6, size = 3.1) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.14, 0.16))) +
    ggplot2::scale_x_discrete(limits = plotgroup_levels) +
    bioRtools::scale_fill_research() +
    ggplot2::labs(x = "Stage and planting mode", y = paste0(axis_name, " score")) +
    bioRtools::theme_prism() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1)) -> p
  save_panel(p, pid, label, 8, 6)
}
plot_axis("Dim 1", "27", "Figure2C_Dim1分布_three_stages_median_matrix")
plot_axis("Dim 2", "28", "Figure2D_Dim2分布_three_stages_median_matrix")

loadings_tab %>%
  dplyr::filter(gene %in% vip_top10$gene) %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene") %>%
  dplyr::left_join(vip_top10, by = "gene") -> load_top

load_top %>%
  ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_segment(aes(x = 0, y = 0, xend = p1, yend = p2, color = process),
                        arrow = ggplot2::arrow(length = ggplot2::unit(0.22, "cm"), type = "closed"),
                        linewidth = 0.6) +
  ggrepel::geom_label_repel(aes(x = p1, y = p2, label = gene, color = process),
                            size = 3.4, fontface = "italic", max.overlaps = 30, show.legend = FALSE) +
  bioRtools::scale_color_research() +
  ggplot2::labs(x = pc_lab(1), y = pc_lab(2), color = "Pathway") +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "right") -> p.2e
save_panel(p.2e, "29", "Figure2E_VIP载荷图_three_stages_median_matrix", 9, 6)

cat("\n[15] sessionInfo()\n")
sessionInfo()
