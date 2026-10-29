# ============================================================
# 04_pathway_gene_expression / 03_Figure3A_3B_defence_pathways.R
# 作用：Figure 3A（防御相关gene热图）、Figure 3B（防御相关gene整体expression差异）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_pathway_gene_expression/03_Figure3A_3B_defence_pathways.R"
#
# 输入：<root>/result/04_pathway_gene_expression/01_表_pathwaygene差异.xlsx（clean_long、stage_level）
# 输出：<root>/result/04_pathway_gene_expression/03_Figure3A_防御gene热图.pdf / .png
#       <root>/result/04_pathway_gene_expression/04_Figure3B_防御gene整体差异.pdf / .png
# 其中 <root> = <repository_root>
#
# gene集：process ∈ {Defense, Jasmonic acid, Salicylic acid}（and旧版一致，共 17 个gene）。
# 缩放：expression在gene内做 0–1 缩放（and旧版 bioRtools::scale01() 一致），跨gene可比性由该缩放保证。
# 热图实现：and主题 03 相同，用 ggplot2::geom_tile 直接绘制（pheatmapColor 未安装、ggplotify 转换报错）。
# statistics：Figure 3B 的 P 值取 01 脚本 stage 层面、对"每个gene 0–1 缩放值"再做一次 Wilcoxon（and旧版口径一致），
#   并把results写入 result/04_pathway_gene_expression/05_表_防御整体test.xlsx。
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
stopifnot(nrow(df.clean) > 0)

defense_processes <- c("Defense", "Jasmonic acid", "Salicylic acid")

df.clean %>%
  dplyr::filter(process %in% defense_processes) %>%
  dplyr::mutate(
    stage = factor(stage, levels = c("tillering", "filling")),
    group = factor(group, levels = c("Pure", "Mix")),
    process = factor(process, levels = defense_processes)
  ) -> df.def

cat("[03] 防御相关gene：", dplyr::n_distinct(df.def$gene), "个；sample行数：", nrow(df.def), "\n")

# gene内 0–1 缩放
df.def %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(rel_expression = as.vector(bioRtools::scale01(expression))) %>%
  dplyr::ungroup() -> df.def

# ---- Figure 3A：热图（gene × 32 组合）----
df.def %>%
  dplyr::group_by(stage, focal, neighbor, gene, process) %>%
  dplyr::summarise(mean_rel = mean(rel_expression), .groups = "drop") %>%
  dplyr::mutate(
    combo_label = paste0(focal, "-", neighbor),
    group = dplyr::if_else(focal == neighbor, "Pure", "Mix"),
    stage = factor(stage, levels = c("tillering", "filling")),
    group = factor(group, levels = c("Pure", "Mix"))
  ) -> df.def.combo

# gene排序：按pathway、再按 Mix–Pure 的平均差
gene_order <- df.def.combo %>%
  dplyr::group_by(gene, process) %>%
  dplyr::summarise(
    diff_mix_pure = mean(mean_rel[group == "Mix"]) - mean(mean_rel[group == "Pure"]),
    .groups = "drop"
  ) %>%
  dplyr::arrange(process, diff_mix_pure) %>%
  dplyr::pull(gene)

df.def.combo %>%
  dplyr::mutate(gene = factor(gene, levels = gene_order)) %>%
  ggplot2::ggplot(aes(x = combo_label, y = gene, fill = mean_rel)) +
  ggplot2::geom_tile() +
  ggplot2::facet_grid(stage ~ group, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_gradientn(
    colours = c("#f2f2f0", "#f2b8a0", "#a63446"),
    name = "Relative expression (0-1)"
  ) +
  ggplot2::labs(x = "Focal-neighbor combination", y = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
    axis.text.y = ggplot2::element_text(size = 7, face = "italic"),
    strip.text = ggplot2::element_text(size = 10, face = "bold"),
    panel.grid = ggplot2::element_blank(),
    legend.position = "right"
  ) -> p.3a

bioRtools::write_data(p.3a, file = file.path(out_dir, "03_Figure3A_防御gene热图.pdf"),
                      width = 10, height = 6, units = "in")
bioRtools::write_data(p.3a, file = file.path(out_dir, "03_Figure3A_防御gene热图.png"),
                      width = 10, height = 6, units = "in", dpi = 500)
cat("[03] 输出：03_Figure3A_防御gene热图.pdf / .png\n")

# ---- Figure 3B：整体expression差异（stage内 Mix vs Pure）----
stat_tab <- df.def %>%
  dplyr::group_by(stage) %>%
  dplyr::group_modify(~ {
    ww <- suppressWarnings(stats::wilcox.test(rel_expression ~ group, data = .x))
    tt <- stats::t.test(rel_expression ~ group, data = .x)
    tibble::tibble(
      n_pure = sum(.x$group == "Pure"),
      n_mix = sum(.x$group == "Mix"),
      mean_pure = mean(.x$rel_expression[.x$group == "Pure"]),
      mean_mix = mean(.x$rel_expression[.x$group == "Mix"]),
      wilcoxon_p = ww$p.value,
      welch_p = tt$p.value
    )
  }) %>%
  dplyr::ungroup()

cat("\n[03] Figure 3B statistics（gene内 0–1 缩放后的整体comparison）：\n")
print(as.data.frame(stat_tab))

df.def %>%
  dplyr::left_join(stat_tab %>% dplyr::select(stage, wilcoxon_p), by = "stage") %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(y_pos = max(rel_expression) + 0.06 * diff(range(rel_expression))) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot(aes(x = group, y = rel_expression, fill = group)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.35, size = 0.8) +
  ggplot2::geom_text(
    aes(x = 1.5, y = y_pos, label = format_p(wilcoxon_p)),
    inherit.aes = FALSE, size = 3.8
  ) +
  ggplot2::facet_wrap(~stage, scales = "free") +
  ggplot2::scale_x_discrete(limits = c("Pure", "Mix")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.14))) +
  bioRtools::scale_fill_research() +
  ggplot2::labs(x = NULL, y = "Relative expression (0-1, per gene)") +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(face = "bold")) -> p.3b

bioRtools::write_data(p.3b, file = file.path(out_dir, "04_Figure3B_防御gene整体差异.pdf"),
                      width = 8, height = 6, units = "in")
bioRtools::write_data(p.3b, file = file.path(out_dir, "04_Figure3B_防御gene整体差异.png"),
                      width = 8, height = 6, units = "in", dpi = 500)
cat("[03] 输出：04_Figure3B_防御gene整体差异.pdf / .png\n")

bioRtools::write_data(stat_tab, file = file.path(out_dir, "05_表_防御整体test.xlsx"))
cat("[03] 输出：05_表_防御整体test.xlsx\n")

cat("\n[03] sessionInfo()\n")
sessionInfo()
