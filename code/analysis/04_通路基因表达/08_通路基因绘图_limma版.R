# ============================================================
# 04_通路基因表达 / 08_通路基因绘图_limma版.R
# 作用：用统一口径（limma：focal 区组 + 经验贝叶斯）重绘 Figure 3B / 3C / 3D 与 Figure 4A / 4B / 4C
#       —— 与 03/04/05 脚本的 Wilcoxon 版图并列保留，便于方法学对照
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_通路基因表达/08_通路基因绘图_limma版.R"
#
# 输入：<root>/result/04_通路基因表达/14_表_limma统一口径统计.xlsx
#       （gene_stage_level、gene_stage_focal_level、defense_score、pathway_definition）
#       <root>/result/03_基因表达总体/01_表_表达矩阵.xlsx（sample_matrix，用于画防御基因整体分数）
# 输出：<root>/result/04_通路基因表达/15_Figure3B_防御基因整体分数.pdf / .png
#       16_Figure3C_69分蘖期_limma版、17_Figure3D_87灌浆期_limma版
#       18_Figure4A_防御基因_limma版、19_Figure4B_独脚金内酯_limma版、20_Figure4C_高糖通路_limma版
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 绘图约定：x 轴为 limma 的 moderated logFC（log2(x+1) 尺度、经 focal 校正），
#   颜色表示显著性层级（FDR < 0.05 / P < 0.05 / n.s.），点旁给出 logFC 数值；
#   不画误差棒（limma 的 moderated t 已把方差收缩，误差棒会与 logFC 尺度混淆），显著性以颜色与 FDR 标注表达。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录见 code/04_通路基因表达/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_stat <- file.path(project_root, "result", "04_通路基因表达", "14_表_limma统一口径统计.xlsx")
path_mat <- file.path(project_root, "result", "03_基因表达总体", "01_表_表达矩阵.xlsx")
out_dir <- file.path(project_root, "result", "04_通路基因表达")

stopifnot(file.exists(path_stat), file.exists(path_mat))

sig_levels <- c("FDR < 0.05", "P < 0.05", "n.s.")
sig_colors <- c("FDR < 0.05" = "#a63446", "P < 0.05" = "#e0a458", "n.s." = "grey62")

bioRtools::read_data(path_stat, sheet = "gene_stage_level") -> gene_stage
bioRtools::read_data(path_stat, sheet = "gene_stage_focal_level") -> gene_focal
bioRtools::read_data(path_stat, sheet = "pathway_definition") -> pathway_def
bioRtools::read_data(path_stat, sheet = "defense_score") -> defense_score

# 显著性分级工具
add_sig <- function(d, fdr_col, p_col) {
  d %>%
    dplyr::mutate(
      sig = dplyr::case_when(
        .data[[fdr_col]] < 0.05 ~ "FDR < 0.05",
        .data[[p_col]] < 0.05 ~ "P < 0.05",
        TRUE ~ "n.s."
      ),
      sig = factor(sig, levels = sig_levels)
    )
}

# ---- 通用：logFC 点图 ----
plot_logfc <- function(d, x_col, y_col = "gene", facet_col = NULL, title = NULL) {
  p <- d %>%
    ggplot2::ggplot(aes(x = .data[[x_col]], y = forcats::fct_reorder(.data[[y_col]], .data[[x_col]]))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
    ggplot2::geom_segment(aes(x = 0, xend = .data[[x_col]], yend = forcats::fct_reorder(.data[[y_col]], .data[[x_col]])),
                          color = "grey80", linewidth = 0.4, show.legend = FALSE) +
    ggplot2::geom_point(aes(color = sig), size = 3) +
    ggplot2::geom_text(aes(label = formatC(.data[[x_col]], format = "f", digits = 2)),
                       hjust = -0.25, size = 2.7, color = "grey25") +
    ggplot2::scale_color_manual(values = sig_colors, name = NULL) +
    ggplot2::labs(x = "Moderated log2 fold change (Mix vs Pure)", y = NULL, title = title) +
    bioRtools::theme_prism() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8, face = "italic"),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(size = 10, face = "bold")
    )

  if (!is.null(facet_col)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_col)), scales = "free_y")
  }
  p
}

save_fig <- function(p, pid, label, width_in, height_in) {
  bioRtools::write_data(p, file = file.path(out_dir, paste0(pid, "_", label, ".pdf")),
                        width = width_in, height = height_in, units = "in")
  bioRtools::write_data(p, file = file.path(out_dir, paste0(pid, "_", label, ".png")),
                        width = width_in, height = height_in, units = "in", dpi = 500)
  cat("[08] 输出：", pid, "_", label, ".pdf / .png\n", sep = "")
}

# ---- Figure 3B：防御基因整体分数（每样本一个分数，focal 校正）----
read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}
read_sheet_raw(path_mat, "sample_matrix") -> mat.sample

meta_cols <- c("stage", "sample", "focal", "neighbor")
gene_cols <- setdiff(colnames(mat.sample), meta_cols)
def_genes <- pathway_def %>% dplyr::filter(process == "Defense") %>% dplyr::pull(genes)

mat.sample %>%
  dplyr::mutate(group = factor(dplyr::if_else(focal == neighbor, "Pure", "Mix"),
                              levels = c("Pure", "Mix"))) %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(gene_cols), ~ as.vector(scale(.x)))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(defense_score = rowMeans(dplyr::across(dplyr::all_of(def_genes)))) %>%
  dplyr::mutate(stage = factor(stage, levels = c("tillering", "filling"))) -> d.score

d.score %>%
  dplyr::left_join(defense_score %>% dplyr::select(stage, p_value), by = "stage") %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(y_pos = max(defense_score) + 0.06 * diff(range(defense_score))) %>%
  dplyr::ungroup() -> d.score

p.3b <- d.score %>%
  ggplot2::ggplot(aes(x = group, y = defense_score, fill = group)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  ggplot2::geom_text(aes(x = 1.5, y = y_pos,
                         label = paste0("P = ", formatC(p_value, format = "f", digits = 3))),
                     inherit.aes = FALSE, size = 3.6) +
  ggplot2::facet_wrap(~stage) +
  ggplot2::scale_x_discrete(limits = c("Pure", "Mix")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.14))) +
  bioRtools::scale_fill_research() +
  ggplot2::labs(
    x = NULL,
    y = "Defense pathway score (mean Z of 12 genes)",
    caption = "One value per sample (n = 80 tillering, 88 filling); P from lm(score ~ group + focal)"
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_text(face = "bold"),
    plot.caption = ggplot2::element_text(size = 7, hjust = 0)
  )

save_fig(p.3b, "15", "Figure3B_防御基因整体分数", 8, 6)

# ---- Figure 3C / 3D：焦点纯系（limma，focal 内）----
gene_focal_sig <- gene_focal %>% add_sig("focal_fdr", "focal_p")

p.3c <- gene_focal_sig %>%
  dplyr::filter(stage == "tillering", focal == "69",
                process %in% c("Brassinosteroid", "Auxin"), gene != "OsXTR1") %>%
  plot_logfc("focal_logFC", title = "Focal 69, tillering: BR + auxin (limma, focal-stratified)")

save_fig(p.3c, "16", "Figure3C_69分蘖期_limma版", 7, 4.5)

p.3d <- gene_focal_sig %>%
  dplyr::filter(stage == "filling", focal == "87",
                process %in% c("Brassinosteroid", "Salicylic acid", "Drought"), gene != "OsXTR1") %>%
  plot_logfc("focal_logFC", title = "Focal 87, filling: BR + SA + drought (limma, focal-stratified)")

save_fig(p.3d, "17", "Figure3D_87灌浆期_limma版", 7, 6)

# ---- Figure 4A / 4B / 4C：跨阶段（limma，stage 层面）----
gene_stage_sig <- gene_stage %>% add_sig("limma_fdr", "limma_p")

plot_cross <- function(processes, pid, label, width_in, height_in, title) {
  d <- gene_stage_sig %>%
    dplyr::filter(process %in% processes) %>%
    dplyr::mutate(stage = factor(stage, levels = c("tillering", "filling")))

  p <- plot_logfc(d, "limma_logFC", facet_col = "stage", title = title)
  save_fig(p, pid, label, width_in, height_in)
}

plot_cross("Defense", "18", "Figure4A_防御基因_limma版", 13, 7,
           "Defense genes (12), stage-level limma (focal as block)")
plot_cross("Strigolactone", "19", "Figure4B_独脚金内酯_limma版", 8, 4,
           "Strigolactone genes (3), stage-level limma")
plot_cross("High sugar", "20", "Figure4C_高糖通路_limma版", 9, 5,
           "High sugar genes (5), stage-level limma")

cat("\n[08] sessionInfo()\n")
sessionInfo()
