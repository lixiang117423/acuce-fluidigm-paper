# ============================================================
# 04_pathway_gene_expression / 06_TableS5_high_sugar_defence_correlation.R
# 作用：高糖pathwaygeneand防御/JA/SA pathway_gene_expression的相关性 —— Table S5
#       （旧版脚本最后一段只做 cor_analysis() 打印、未落盘，本次补成正式表格and图）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_pathway_gene_expression/06_TableS5_high_sugar_defence_correlation.R"
#
# 输入：<root>/result/04_pathway_gene_expression/01_表_pathwaygene差异.xlsx（clean_long）
# 输出：<root>/result/04_pathway_gene_expression/12_表_TableS5_high_sugar_defence_correlation.xlsx
#       <root>/result/04_pathway_gene_expression/13_FigureS5_高糖and防御相关性热图.pdf / .png
# 其中 <root> = <repository_root>
#
# 方法：sample级（stage × sample）expression矩阵，取 High sugar 5 个gene × {Defense, Jasmonic acid, Salicylic acid} 17 个gene，
#   计算 Pearson 相关（每个时期单独算），并做 BH-FDR correction；相关性热图展示 r 值。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/04_pathway_gene_expression/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
out_dir <- file.path(project_root, "result", "04_pathway_gene_expression")
path_in <- file.path(out_dir, "01_表_pathwaygene差异.xlsx")

bioRtools::read_data(path_in, sheet = "clean_long") -> df.clean

sugar_genes <- df.clean %>% dplyr::filter(process == "High sugar") %>% dplyr::distinct(gene) %>% dplyr::pull(gene)
def_genes <- df.clean %>%
  dplyr::filter(process %in% c("Defense", "Jasmonic acid", "Salicylic acid")) %>%
  dplyr::distinct(gene) %>%
  dplyr::pull(gene)

cat("[06] High sugar gene：", paste(sugar_genes, collapse = ", "), "\n")
cat("[06] 防御相关gene数：", length(def_genes), "\n")

# sample级宽表（每个sample一个值：sample内先对重复取均值）
df.clean %>%
  dplyr::mutate(sample_id = paste0(stage, "_", sample)) %>%
  dplyr::group_by(sample_id, stage, gene) %>%
  dplyr::summarise(expression = mean(expression), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = gene, values_from = expression) -> mat.sample

cat("[06] sample级矩阵：", nrow(mat.sample), "行 ×", ncol(mat.sample) - 2, "gene\n")

cor_tab <- purrr::map_dfr(c("tillering", "filling"), function(st) {
  d <- mat.sample %>% dplyr::filter(stage == st)
  purrr::map_dfr(sugar_genes, function(g1) {
    purrr::map_dfr(def_genes, function(g2) {
      ct <- stats::cor.test(d[[g1]], d[[g2]], method = "pearson")
      tibble::tibble(
        stage = st, sugar_gene = g1, defense_gene = g2,
        n = nrow(d), r = unname(ct$estimate), p_value = ct$p.value
      )
    })
  })
}) %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(stage, p_value)

cat("\n[06] 相关性最强的前 15 个组合：\n")
cor_tab %>%
  dplyr::select(stage, sugar_gene, defense_gene, r, p_value, p_fdr) %>%
  dplyr::slice_head(n = 15) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n[06] 可计算相关的组合数：", sum(!is.na(cor_tab$p_value)), " / ", nrow(cor_tab),
    "（NA 来自expression无变异的gene，如 OsGH9B9 / OsPR1b / OsSPL7 等大量为 0）\n")
cat("[06] FDR < 0.05 的组合数：", sum(cor_tab$p_fdr < 0.05, na.rm = TRUE),
    "；未correction P < 0.05：", sum(cor_tab$p_value < 0.05, na.rm = TRUE), "\n")

bioRtools::write_data(cor_tab, file = file.path(out_dir, "12_表_TableS5_high_sugar_defence_correlation.xlsx"))
cat("[06] 输出：12_表_TableS5_high_sugar_defence_correlation.xlsx\n")

# ---- 相关性热图 ----
cor_tab %>%
  dplyr::mutate(
    stage = factor(stage, levels = c("tillering", "filling")),
    sugar_gene = factor(sugar_gene, levels = sugar_genes),
    defense_gene = factor(defense_gene, levels = def_genes)
  ) %>%
  ggplot2::ggplot(aes(x = sugar_gene, y = defense_gene, fill = r)) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(aes(label = formatC(r, format = "f", digits = 1)), size = 2.6) +
  ggplot2::facet_wrap(~stage) +
  ggplot2::scale_fill_gradient2(low = "#a63446", mid = "#f7f7f5", high = "#1e4e8c",
                                midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  ggplot2::labs(x = "High sugar pathway gene", y = "Defense-related gene") +
  bioRtools::theme_prism() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 8, face = "italic"),
    axis.text.y = ggplot2::element_text(size = 7, face = "italic"),
    strip.text = ggplot2::element_text(face = "bold"),
    panel.grid = ggplot2::element_blank()
  ) -> p.s5

bioRtools::write_data(p.s5, file = file.path(out_dir, "13_FigureS5_高糖and防御相关性热图.pdf"),
                      width = 8, height = 8, units = "in")
bioRtools::write_data(p.s5, file = file.path(out_dir, "13_FigureS5_高糖and防御相关性热图.png"),
                      width = 8, height = 8, units = "in", dpi = 500)
cat("[06] 输出：13_FigureS5_高糖and防御相关性热图.pdf / .png\n")

cat("\n[06] sessionInfo()\n")
sessionInfo()
