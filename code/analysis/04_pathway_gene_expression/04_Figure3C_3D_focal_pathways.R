# ============================================================
# 04_pathway_gene_expression / 04_Figure3C_3D_focal_pathways.R
# 作用：Figure 3C（69 分蘖期 BR + AUX）、Figure 3D（87 灌浆期 BR + SA + Drought）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_pathway_gene_expression/04_Figure3C_3D_focal_pathways.R"
#
# 输入：<root>/result/04_pathway_gene_expression/01_表_pathwaygene差异.xlsx（clean_long、stage_focal_level）
# 输出：<root>/result/04_pathway_gene_expression/06_Figure3C_69分蘖期_BR_AUX.pdf / .png
#       <root>/result/04_pathway_gene_expression/07_Figure3D_87灌浆期_BR_SA_Drought.pdf / .png
# 其中 <root> = <repository_root>
#
# gene集（and旧版一致）：3C = 69 × tillering × {Brassinosteroid, Auxin}；
#   3D = 87 × filling × {Brassinosteroid, Salicylic acid, Drought}。
# OsXTR1 的处理：该gene在 Mix 条件下expression近似 0（多数sample为 0，log2FC = −6.5 / −Inf），
#   旧版用 `filter(gene != "OsXTR1")` 直接剔除但未notes原因。本次**同样不纳入图内**，但在 README and该gene的statistics行中
#   明确写出原因（expression全为 0 导致比值and可视化无意义），并保留在statistics表里，便于溯源。
# 缩放：expression在gene内做 0–1 缩放；P 值取 01 脚本 stage × focal 层面的 Wilcoxon results。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/04_pathway_gene_expression/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
out_dir <- file.path(project_root, "result", "04_pathway_gene_expression")
path_in <- file.path(out_dir, "01_表_pathwaygene差异.xlsx")

format_p <- function(p, digits = 3) {
  dplyr::if_else(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, format = "f", digits = digits)))
}

bioRtools::read_data(path_in, sheet = "clean_long") -> df.clean
bioRtools::read_data(path_in, sheet = "stage_focal_level") -> focal_tab

plot_focal_pathway <- function(stage_sel, focal_sel, processes, pid, plabel, ncol_facet) {
  df.clean %>%
    dplyr::filter(stage == stage_sel, focal == focal_sel, process %in% processes) %>%
    dplyr::mutate(
      group = factor(group, levels = c("Pure", "Mix")),
      process = factor(process, levels = processes)
    ) %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(rel_expression = as.vector(bioRtools::scale01(expression))) %>%
    dplyr::ungroup() -> d

  focal_tab %>%
    dplyr::filter(stage == stage_sel, focal == focal_sel, process %in% processes) %>%
    dplyr::select(gene, process, wilcoxon_p, log2FC) -> st

  d %>%
    dplyr::left_join(st, by = c("gene", "process")) %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(y_pos = max(rel_expression) + 0.05 * diff(range(rel_expression))) %>%
    dplyr::ungroup() -> d

  cat("\n[04]", plabel, " gene水平statistics：\n")
  print(as.data.frame(st %>% dplyr::arrange(process, gene)))

  d %>%
    dplyr::filter(gene != "OsXTR1") %>%
    ggplot2::ggplot(aes(x = group, y = rel_expression, fill = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.4, size = 0.9) +
    ggplot2::geom_text(aes(x = 1.5, y = y_pos, label = format_p(wilcoxon_p)),
                       inherit.aes = FALSE, size = 3.2) +
    ggplot2::facet_wrap(process ~ gene, scales = "free", ncol = ncol_facet) +
    ggplot2::scale_x_discrete(limits = c("Pure", "Mix")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.14))) +
    bioRtools::scale_fill_research() +
    ggplot2::labs(x = NULL, y = "Relative expression (0-1, per gene)") +
    bioRtools::theme_prism() +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(size = 8, face = "bold")
    ) -> p

  bioRtools::write_data(p, file = file.path(out_dir, paste0(pid, "_", plabel, ".pdf")),
                        width = 8, height = 8, units = "in")
  bioRtools::write_data(p, file = file.path(out_dir, paste0(pid, "_", plabel, ".png")),
                        width = 8, height = 8, units = "in", dpi = 500)
  cat("[04] 输出：", pid, "_", plabel, ".pdf / .png\n", sep = "")
}

plot_focal_pathway("tillering", 69, c("Brassinosteroid", "Auxin"),
                   "06", "Figure3C_69分蘖期_BR_AUX", 3)

plot_focal_pathway("filling", 87, c("Brassinosteroid", "Salicylic acid", "Drought"),
                   "07", "Figure3D_87灌浆期_BR_SA_Drought", 3)

cat("\n[04] sessionInfo()\n")
sessionInfo()
