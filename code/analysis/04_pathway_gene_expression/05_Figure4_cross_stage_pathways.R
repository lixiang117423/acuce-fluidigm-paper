# ============================================================
# 04_pathway_gene_expression / 05_Figure4_cross_stage_pathways.R
# 作用：Figure 4A（12 个防御genecross_stage）、Figure 4B（独脚金内酯 3 genecross_stage，**旧版缺失**）、
#       Figure 4C（高糖pathway 5 genecross_stage）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_pathway_gene_expression/05_Figure4_cross_stage_pathways.R"
#
# 输入：<root>/result/04_pathway_gene_expression/01_表_pathwaygene差异.xlsx（clean_long、stage_level）
# 输出：<root>/result/04_pathway_gene_expression/08_Figure4A_防御genecross_stage.pdf / .png
#       <root>/result/04_pathway_gene_expression/09_Figure4B_独脚金内酯cross_stage.pdf / .png
#       <root>/result/04_pathway_gene_expression/10_Figure4C_高糖pathwaycross_stage.pdf / .png
#       <root>/result/04_pathway_gene_expression/11_表_cross_stagepathway_gene_statistics.xlsx
# 其中 <root> = <repository_root>
#
# 编号notes（2026-09-22 确认）：3A 防御热图 / 3B 防御整体差异 / 3C 69 分蘖期 BR+AUX / 3D 87 灌浆期 BR+SA+Drought；
#   4A defence_pathways分gene / 4B 独脚金内酯 / 4C 高糖pathway。旧版把"defence_pathways分gene"放在 4A、把"独脚金内酯"漏掉，
#   本次按框架complete 4B。
# gene取舍：**不再沿用旧版对 OsGH9B9 / OsPR1b / OsPR5 / OsSPL7 的无理由剔除**，Figure 4A 展示全部 12 个防御gene
#   （OsPR1b / OsPR5 的expression本身较低，但仍如实展示，statisticsresults见 01_表_pathwaygene差异.xlsx）。
# 缩放andtest：expression在gene内做 0–1 缩放；P 值取 01 脚本 stage 层面（4 个 focal 合并）的 Wilcoxon results。
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
bioRtools::read_data(path_in, sheet = "stage_level") -> stage_tab

plot_cross_stage <- function(processes, pid, plabel, ncol_facet, width_in, height_in) {
  df.clean %>%
    dplyr::filter(process %in% processes) %>%
    dplyr::mutate(
      stage = factor(stage, levels = c("tillering", "filling")),
      group = factor(group, levels = c("Pure", "Mix"))
    ) %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(rel_expression = as.vector(bioRtools::scale01(expression))) %>%
    dplyr::ungroup() -> d

  stage_tab %>%
    dplyr::filter(process %in% processes) %>%
    dplyr::select(stage, gene, process, mean_mix, mean_pure, log2FC, wilcoxon_p, wilcoxon_fdr) -> st

  d %>%
    dplyr::left_join(st, by = c("stage", "gene", "process")) %>%
    dplyr::group_by(stage, gene) %>%
    dplyr::mutate(y_pos = max(rel_expression) + 0.05 * diff(range(rel_expression))) %>%
    dplyr::ungroup() -> d

  cat("\n[05]", plabel, " statistics（stage 层面）：\n")
  print(as.data.frame(st %>% dplyr::arrange(gene, stage)))

  d %>%
    ggplot2::ggplot(aes(x = group, y = rel_expression, fill = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.35, size = 0.8) +
    ggplot2::geom_text(aes(x = 1.5, y = y_pos, label = format_p(wilcoxon_p)),
                       inherit.aes = FALSE, size = 3) +
    ggplot2::facet_wrap(stage ~ gene, scales = "free", ncol = ncol_facet) +
    ggplot2::scale_x_discrete(limits = c("Pure", "Mix")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.14))) +
    bioRtools::scale_fill_research() +
    ggplot2::labs(x = NULL, y = "Relative expression (0-1, per gene)") +
    bioRtools::theme_prism() +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(size = 7, face = "bold")
    ) -> p

  bioRtools::write_data(p, file = file.path(out_dir, paste0(pid, "_", plabel, ".pdf")),
                        width = width_in, height = height_in, units = "in")
  bioRtools::write_data(p, file = file.path(out_dir, paste0(pid, "_", plabel, ".png")),
                        width = width_in, height = height_in, units = "in", dpi = 500)
  cat("[05] 输出：", pid, "_", plabel, ".pdf / .png\n", sep = "")

  st
}

st_4a <- plot_cross_stage("Defense", "08", "Figure4A_防御genecross_stage", 6, 15, 6)
st_4b <- plot_cross_stage("Strigolactone", "09", "Figure4B_独脚金内酯cross_stage", 3, 6, 6)
st_4c <- plot_cross_stage("High sugar", "10", "Figure4C_高糖pathwaycross_stage", 5, 13.3, 5)

bioRtools::write_data(
  list(defense = st_4a, strigolactone = st_4b, high_sugar = st_4c),
  file = file.path(out_dir, "11_表_cross_stagepathway_gene_statistics.xlsx")
)
cat("[05] 输出：11_表_cross_stagepathway_gene_statistics.xlsx\n")

cat("\n[05] sessionInfo()\n")
sessionInfo()
