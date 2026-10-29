# ============================================================
# 04_通路基因表达 / 17_通路基因绘图_limma版_三时期_投稿规格.R
# 作用：**投稿规格版** Figure 3B 与 Figure 4A / 4C（另含 24 号两时期热图、27 号独脚金内酯=Figure S17）
#       与 11 号脚本的差别只在版面：输出尺寸改为最终 rect 尺寸（cm，scale = 1）、
#       去掉面板内 caption、面板标题缩短；输出目录改为 result/08_投稿规格图（原目录只读不改）。
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_通路基因表达/17_通路基因绘图_limma版_三时期_投稿规格.R"
#
# 输入：<root>/result/04_通路基因表达/23_表_limma统一口径统计_三时期.xlsx
#       <root>/result/03_基因表达总体/11_表_表达矩阵_三时期.xlsx（combo_matrix、sample_matrix、gene_annotation）
# 输出：<root>/result/04_通路基因表达/
#       24_Figure3A_防御基因热图_三时期、25_Figure3B_防御基因整体分数_三时期、
#       26_Figure4A_防御基因_limma版_三时期、27_Figure4B_独脚金内酯_limma版_三时期、
#       28_Figure4C_高糖通路_limma版_三时期（各 PDF + PNG）
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 绘图约定：bioRtools::theme_prism() / scale_*_research() / write_data(units = "in", dpi = 500)；
#   图内文字英文；显著性用颜色分层（FDR < 0.05 / P < 0.05 / n.s.）。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_stat <- file.path(project_root, "result", "04_通路基因表达", "23_表_limma统一口径统计_三时期.xlsx")
path_mat <- file.path(project_root, "result", "03_基因表达总体", "11_表_表达矩阵_三时期.xlsx")
out_dir <- file.path(project_root, "result", "04_通路基因表达")
out_pub_dir <- file.path(project_root, "result", "08_投稿规格图")  # 投稿规格面板输出目录
dir.create(out_pub_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(path_stat), file.exists(path_mat))

stage_levels <- c("tillering", "flowering", "filling")
sig_levels <- c("FDR < 0.05", "P < 0.05", "n.s.")
sig_colors <- c("FDR < 0.05" = "#a63446", "P < 0.05" = "#e0a458", "n.s." = "grey62")

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

bioRtools::read_data(path_stat, sheet = "gene_stage_level") -> gene_stage
bioRtools::read_data(path_stat, sheet = "pathway_definition") -> pathway_def
bioRtools::read_data(path_stat, sheet = "defense_score") -> defense_score
read_sheet_raw(path_mat, "combo_matrix") -> mat.combo
read_sheet_raw(path_mat, "sample_matrix") -> mat.sample
read_sheet_raw(path_mat, "gene_annotation") -> df.gene

meta_cols_combo <- c("combo", "stage", "focal", "neighbor", "group")
gene_cols <- setdiff(colnames(mat.combo), meta_cols_combo)

add_sig <- function(d, fdr_col, p_col) {
  d %>%
    dplyr::mutate(
      sig = dplyr::case_when(.data[[fdr_col]] < 0.05 ~ "FDR < 0.05",
                             .data[[p_col]] < 0.05 ~ "P < 0.05",
                             TRUE ~ "n.s."),
      sig = factor(sig, levels = sig_levels)
    )
}

save_fig <- function(p, pid, label, width_cm, height_cm) {
  bioRtools::write_data(p, file = file.path(out_pub_dir, paste0(pid, "_", label, ".pdf")),
                        width = width_cm, height = height_cm, units = "cm")
  bioRtools::write_data(p, file = file.path(out_pub_dir, paste0(pid, "_", label, ".png")),
                        width = width_cm, height = height_cm, units = "cm", dpi = 600)
  cat("[17] 输出：", pid, "_", label, ".pdf / .png\n", sep = "")
}

# ---- Figure 3A：防御相关基因热图（Defense + JA + SA）× 48 组合 ----
def_processes <- c("Defense", "Jasmonic acid", "Salicylic acid")
def_genes <- pathway_def %>%
  dplyr::filter(process %in% def_processes) %>%
  dplyr::distinct(gene) %>%
  dplyr::pull(gene)

cat("[11] Figure 3A 基因数：", length(def_genes), "\n")

mat.def <- as.matrix(mat.combo[, intersect(def_genes, gene_cols)])
rownames(mat.def) <- mat.combo$combo

scale01_by_row <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

mat.def.scaled <- t(apply(mat.def, 1, scale01_by_row))

as.data.frame(mat.def.scaled) %>%
  tibble::rownames_to_column(var = "combo") %>%
  tidyr::pivot_longer(cols = -combo, names_to = "gene", values_to = "scaled_expression") %>%
  dplyr::left_join(mat.combo[, meta_cols_combo], by = "combo") %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene") %>%
  dplyr::mutate(
    stage = factor(stage, levels = stage_levels),
    group = factor(group, levels = c("Pure", "Mix")),
    process = factor(process, levels = def_processes),
    gene = factor(gene, levels = rev(sort(unique(gene)))),
    combo_label = paste0(focal, "-", neighbor)
  ) -> df.def.long

df.def.long %>%
  ggplot2::ggplot(aes(x = combo_label, y = gene, fill = scaled_expression)) +
  ggplot2::geom_tile() +
  ggplot2::facet_grid(process + stage ~ group, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_gradientn(colours = c("#f2f2f0", "#a8bedb", "#1e4e8c"),
                                name = "Scaled expression (0-1)") +
  ggplot2::labs(x = "Focal-neighbour combination", y = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
    axis.text.y = ggplot2::element_text(size = 7, face = "italic"),
    strip.text = ggplot2::element_text(size = 7, face = "bold"),
    panel.grid = ggplot2::element_blank(),
    legend.position = "right"
  ) -> p.3a

save_fig(p.3a, "24", "Figure3A_防御基因热图_三时期", 8.3, 11.0)

# ---- Figure 3B：防御基因整体分数（每样本一个值）----
meta_cols_sample <- c("stage", "sample", "focal", "neighbor")
gene_cols_sample <- setdiff(colnames(mat.sample), meta_cols_sample)
def_genes_score <- pathway_def %>% dplyr::filter(process == "Defense") %>% dplyr::pull(gene)

mat.sample %>%
  dplyr::mutate(group = factor(dplyr::if_else(focal == neighbor, "Pure", "Mix"),
                               levels = c("Pure", "Mix")),
                stage = factor(stage, levels = stage_levels)) %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(gene_cols_sample), ~ as.vector(scale(.x)))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(defense_score_value = rowMeans(dplyr::across(dplyr::all_of(def_genes_score)))) %>%
  dplyr::left_join(defense_score %>% dplyr::select(stage, p_value, n_sample), by = "stage") %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(y_pos = max(defense_score_value) + 0.06 * diff(range(defense_score_value))) %>%
  dplyr::ungroup() -> d.score

d.score %>%
  ggplot2::ggplot(aes(x = group, y = defense_score_value, fill = group)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  ggplot2::geom_text(aes(x = 1.5, y = y_pos,
                         label = paste0("P = ", formatC(p_value, format = "f", digits = 3))),
                     inherit.aes = FALSE, size = 3.3) +
  ggplot2::facet_wrap(~stage) +
  ggplot2::scale_x_discrete(limits = c("Pure", "Mix")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.14))) +
  bioRtools::scale_fill_research() +
  ggplot2::labs(
    # 2026-10-19：旋转 y 轴标题在面板里两头被切 → 缩短 + 缩小字号 + 上下留边距
    x = NULL, y = "Defence pathway score (Z)") +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "none",
                 strip.text = ggplot2::element_text(face = "bold", size = 8),
                 axis.title.y = ggplot2::element_text(size = 6.5),
                 plot.margin = ggplot2::margin(3, 2, 3, 2)) -> p.3b

save_fig(p.3b, "25", "Figure3B_防御基因整体分数_三时期", 17.0, 5.2)

# ---- Figure 4A / 4B / 4C：跨阶段（limma，stage 层面）----
plot_logfc <- function(d, x_col, facet_col = "stage", title = NULL) {
  d %>%
    ggplot2::ggplot(aes(x = .data[[x_col]],
                        y = forcats::fct_reorder(.data[["gene"]], .data[[x_col]]))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
    ggplot2::geom_point(aes(color = sig), size = 3) +
    ggplot2::geom_text(aes(label = formatC(.data[[x_col]], format = "f", digits = 2)),
                       hjust = -0.25, size = 2.2, color = "grey25") +
    ggplot2::scale_color_manual(values = sig_colors, name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.22)) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", facet_col)), scales = "free_y") +
    ggplot2::labs(x = "Moderated log2 fold change (Mix vs Pure)", y = NULL, title = title) +
    bioRtools::theme_prism() +
    # 2026-10-19 修：8 个基因/分面 + 60 mm 高面板时行距小于字高（实测叠 3–12%），降到 6.5 pt
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6.5, face = "italic"),
                   legend.position = "bottom",
                   panel.spacing = ggplot2::unit(0.5, "cm"),
                   plot.title = ggplot2::element_text(size = 8, face = "bold"))
}

gene_stage_sig <- gene_stage %>%
  add_sig("limma_fdr", "limma_p") %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels))

plot_cross <- function(processes, pid, label, width_in, height_in, title) {
  d <- gene_stage_sig %>% dplyr::filter(process %in% processes)
  stopifnot(nrow(d) > 0)
  save_fig(plot_logfc(d, "limma_logFC", title = title), pid, label, width_in, height_in)
}

plot_cross("Defense", "26", "Figure4A_防御基因_limma版_三时期", 17.0, 9.0,
           "Defence genes (12)")
plot_cross("Strigolactone", "27", "Figure4B_独脚金内酯_limma版_三时期", 17.0, 6.5,
           "Strigolactone genes (3)")
plot_cross("High sugar", "28", "Figure4C_高糖通路_limma版_三时期", 17.0, 6.0,
           "High-sugar genes (5)")

cat("\n[11] sessionInfo()\n")
sessionInfo()
