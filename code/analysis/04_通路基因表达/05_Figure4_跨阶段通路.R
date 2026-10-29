# ============================================================
# 04_通路基因表达 / 05_Figure4_跨阶段通路.R
# 作用：Figure 4A（12 个防御基因跨阶段）、Figure 4B（独脚金内酯 3 基因跨阶段，**旧版缺失**）、
#       Figure 4C（高糖通路 5 基因跨阶段）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_通路基因表达/05_Figure4_跨阶段通路.R"
#
# 输入：<root>/result/04_通路基因表达/01_表_通路基因差异.xlsx（clean_long、stage_level）
# 输出：<root>/result/04_通路基因表达/08_Figure4A_防御基因跨阶段.pdf / .png
#       <root>/result/04_通路基因表达/09_Figure4B_独脚金内酯跨阶段.pdf / .png
#       <root>/result/04_通路基因表达/10_Figure4C_高糖通路跨阶段.pdf / .png
#       <root>/result/04_通路基因表达/11_表_跨阶段通路基因统计.xlsx
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 编号说明（2026-09-22 确认）：3A 防御热图 / 3B 防御整体差异 / 3C 69 分蘖期 BR+AUX / 3D 87 灌浆期 BR+SA+Drought；
#   4A 防御通路分基因 / 4B 独脚金内酯 / 4C 高糖通路。旧版把"防御通路分基因"放在 4A、把"独脚金内酯"漏掉，
#   本次按框架补齐 4B。
# 基因取舍：**不再沿用旧版对 OsGH9B9 / OsPR1b / OsPR5 / OsSPL7 的无理由剔除**，Figure 4A 展示全部 12 个防御基因
#   （OsPR1b / OsPR5 的表达量本身较低，但仍如实展示，统计结果见 01_表_通路基因差异.xlsx）。
# 缩放与检验：表达量在基因内做 0–1 缩放；P 值取 01 脚本 stage 层面（4 个 focal 合并）的 Wilcoxon 结果。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录见 code/04_通路基因表达/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
out_dir <- file.path(project_root, "result", "04_通路基因表达")
path_in <- file.path(out_dir, "01_表_通路基因差异.xlsx")

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

  cat("\n[05]", plabel, " 统计（stage 层面）：\n")
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

st_4a <- plot_cross_stage("Defense", "08", "Figure4A_防御基因跨阶段", 6, 15, 6)
st_4b <- plot_cross_stage("Strigolactone", "09", "Figure4B_独脚金内酯跨阶段", 3, 6, 6)
st_4c <- plot_cross_stage("High sugar", "10", "Figure4C_高糖通路跨阶段", 5, 13.3, 5)

bioRtools::write_data(
  list(defense = st_4a, strigolactone = st_4b, high_sugar = st_4c),
  file = file.path(out_dir, "11_表_跨阶段通路基因统计.xlsx")
)
cat("[05] 输出：11_表_跨阶段通路基因统计.xlsx\n")

cat("\n[05] sessionInfo()\n")
sessionInfo()
