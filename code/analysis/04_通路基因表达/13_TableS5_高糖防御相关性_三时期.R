# ============================================================
# 04_通路基因表达 / 13_TableS5_高糖防御相关性_三时期.R
# 作用：高糖通路基因 × 防御/JA/SA 通路基因的表达相关性（样本级 Pearson），扩展到三个时期 —— Table S5
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_通路基因表达/13_TableS5_高糖防御相关性_三时期.R"
#
# 输入：<root>/result/03_基因表达总体/11_表_表达矩阵_三时期.xlsx（sample_matrix、gene_annotation）
# 输出：<root>/result/04_通路基因表达/31_表_TableS5_高糖防御相关性_三时期.xlsx
#       <root>/result/04_通路基因表达/32_FigureS5_高糖与防御相关性热图_三时期.pdf / .png
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 口径（与两时期脚本 06 一致，仅时期数 2 → 3）：样本级表达矩阵；High sugar 5 个基因 ×
#   {Defense, Jasmonic acid, Salicylic acid} 17 个基因；每个时期单独算 Pearson 相关，
#   **BH-FDR 在每个时期内 85 个组合上校正**（与两时期版"按时期分组校正"一致）。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_mat <- file.path(project_root, "result", "03_基因表达总体", "11_表_表达矩阵_三时期.xlsx")
out_dir <- file.path(project_root, "result", "04_通路基因表达")
path_out <- file.path(out_dir, "31_表_TableS5_高糖防御相关性_三时期.xlsx")

stopifnot(file.exists(path_mat))

stage_levels <- c("tillering", "flowering", "filling")

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_mat, "sample_matrix") -> mat.sample
bioRtools::read_data(path_mat, sheet = "gene_annotation") -> df.gene

sugar_genes <- df.gene %>% dplyr::filter(process == "High sugar") %>% dplyr::pull(gene)
def_genes <- df.gene %>%
  dplyr::filter(process %in% c("Defense", "Jasmonic acid", "Salicylic acid")) %>%
  dplyr::pull(gene)

cat("[13] High sugar 基因（", length(sugar_genes), "）：", paste(sugar_genes, collapse = ", "), "\n")
cat("[13] 防御相关基因数：", length(def_genes), "\n")

mat.sample %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels)) -> d

cor_tab <- purrr::map_dfr(stage_levels, function(st) {
  dd <- d %>% dplyr::filter(stage == st)
  purrr::map_dfr(sugar_genes, function(g1) {
    purrr::map_dfr(def_genes, function(g2) {
      ct <- suppressWarnings(stats::cor.test(dd[[g1]], dd[[g2]], method = "pearson"))
      tibble::tibble(stage = st, sugar_gene = g1, defense_gene = g2,
                     n = sum(stats::complete.cases(dd[[g1]], dd[[g2]])),
                     r = unname(ct$estimate), p_value = ct$p.value)
    })
  })
}) %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(p_fdr = stats::p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(stage, p_value)

cat("\n[13] 每期可计算组合数 / FDR < 0.05 数 / 名义 P < 0.05 数：\n")
print(as.data.frame(cor_tab %>%
                      dplyr::group_by(stage) %>%
                      dplyr::summarise(n_total = dplyr::n(),
                                       n_calc = sum(!is.na(p_value)),
                                       n_fdr005 = sum(p_fdr < 0.05, na.rm = TRUE),
                                       n_p005 = sum(p_value < 0.05, na.rm = TRUE),
                                       min_r_abs_max = max(abs(r), na.rm = TRUE),
                                       .groups = "drop")))

cat("\n[13] 各期相关性最强的 5 个组合：\n")
print(as.data.frame(cor_tab %>% dplyr::group_by(stage) %>%
                      dplyr::slice_min(p_value, n = 5, with_ties = FALSE) %>%
                      dplyr::select(stage, sugar_gene, defense_gene, r, p_value, p_fdr) %>%
                      dplyr::ungroup()))

# ---- 相关性热图（按时期分面）----
cor_tab %>%
  dplyr::mutate(
    stage = factor(stage, levels = stage_levels),
    sig = dplyr::case_when(p_fdr < 0.05 ~ "FDR < 0.05",
                           p_value < 0.05 ~ "P < 0.05",
                           TRUE ~ "n.s."),
    sig = factor(sig, levels = c("FDR < 0.05", "P < 0.05", "n.s.")),
    sugar_gene = factor(sugar_gene, levels = sugar_genes),
    defense_gene = factor(defense_gene, levels = def_genes)
  ) -> d.cor

d.cor %>%
  ggplot2::ggplot(aes(x = defense_gene, y = sugar_gene, fill = r)) +
  ggplot2::geom_tile() +
  ggplot2::geom_point(data = d.cor %>% dplyr::filter(sig != "n.s."),
                      aes(size = sig), shape = 1, color = "grey15", show.legend = TRUE) +
  ggplot2::facet_wrap(~stage) +
  ggplot2::scale_fill_gradient2(low = "#1e4e8c", mid = "#f2f2f0", high = "#a63446",
                                midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  ggplot2::scale_size_manual(values = c("FDR < 0.05" = 2.6, "P < 0.05" = 1.8), name = NULL) +
  ggplot2::labs(x = "Defence-related genes", y = "High-sugar pathway genes",
                caption = paste0("Sample-level Pearson correlations; BH-FDR within each stage (",
                                 length(sugar_genes) * length(def_genes), " combinations per stage)")) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, face = "italic"),
    axis.text.y = ggplot2::element_text(size = 8, face = "italic"),
    strip.text = ggplot2::element_text(face = "bold"),
    panel.grid = ggplot2::element_blank(),
    plot.caption = ggplot2::element_text(size = 7, hjust = 0)
  ) -> p.s5

bioRtools::write_data(p.s5, file = file.path(out_dir, "32_FigureS5_高糖与防御相关性热图_三时期.pdf"),
                      width = 12, height = 6, units = "in")
bioRtools::write_data(p.s5, file = file.path(out_dir, "32_FigureS5_高糖与防御相关性热图_三时期.png"),
                      width = 12, height = 6, units = "in", dpi = 500)
cat("\n[13] 输出：32_FigureS5_高糖与防御相关性热图_三时期.pdf / .png\n")

bioRtools::write_data(list(correlation = cor_tab), file = path_out)
cat("[13] 输出：", path_out, "\n")

cat("\n[13] sessionInfo()\n")
sessionInfo()
