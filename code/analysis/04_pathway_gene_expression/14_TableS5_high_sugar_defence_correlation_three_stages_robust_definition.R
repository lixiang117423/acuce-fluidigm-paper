# ============================================================
# 04_pathway_gene_expression / 14_TableS5_high_sugar_defence_correlation_three_stages_robust_definition.R
# 作用：correction Table S5 的相关口径——**检出率过滤 + Spearman 秩相关**（替代未过滤的 Pearson）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_pathway_gene_expression/14_TableS5_high_sugar_defence_correlation_three_stages_robust_definition.R"
#
# 【为什么改口径】（2026-10-18 实测）
#   脚本 13（未过滤 Pearson）在three_stagesdata上出现 r ≈ 0.9999 的组合，查证为**假相关**：
#   例如 OsPR1b 在分蘖期 80 个sample中有 72 个expression为 0（最大 0.88），
#   OsCRCT and OsPR1b 的 Pearson r = 0.997 完全由 1 个非零点 + 79 个零点决定。
#   旧版两时期 Table S5 同样使用未过滤 Pearson（当时最大 r = 0.79），存在同类风险。
#   因此本脚本：① 只保留"在该时期检出率（expression > 0）≥ 50%"的gene参and相关analysis；
#   ② 用 **Spearman 秩相关**（对非正态、含零的expression更稳健）；③ BH-FDR 在每个时期内correction。
#
# 输入：<root>/result/03_overall_gene_expression/11_表_expression矩阵_three_stages.xlsx（sample_matrix、gene_annotation）
# 输出：<root>/result/04_pathway_gene_expression/33_表_TableS5_high_sugar_defence_correlation_three_stages_robust_definition.xlsx
#       <root>/result/04_pathway_gene_expression/34_FigureS5_高糖and防御相关性热图_three_stages_robust_definition.pdf / .png
# 其中 <root> = <repository_root>
#
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_mat <- file.path(project_root, "result", "03_overall_gene_expression", "11_表_expression矩阵_three_stages.xlsx")
out_dir <- file.path(project_root, "result", "04_pathway_gene_expression")
path_out <- file.path(out_dir, "33_表_TableS5_high_sugar_defence_correlation_three_stages_robust_definition.xlsx")

stopifnot(file.exists(path_mat))

stage_levels <- c("tillering", "flowering", "filling")
detect_cutoff <- 0.5

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_mat, "sample_matrix") -> mat.sample
bioRtools::read_data(path_mat, sheet = "gene_annotation") -> df.gene

sugar_all <- df.gene %>% dplyr::filter(process == "High sugar") %>% dplyr::pull(gene)
def_all <- df.gene %>%
  dplyr::filter(process %in% c("Defense", "Jasmonic acid", "Salicylic acid")) %>%
  dplyr::pull(gene)

mat.sample %>% dplyr::mutate(stage = factor(stage, levels = stage_levels)) -> d

# 每个时期各自计算检出率，并筛选参检gene
detect_tab <- purrr::map_dfr(stage_levels, function(st) {
  dd <- d %>% dplyr::filter(stage == st)
  purrr::map_dfr(c(sugar_all, def_all), function(g) {
    v <- dd[[g]]
    tibble::tibble(stage = st, gene = g, n = length(v),
                   n_detected = sum(v > 0, na.rm = TRUE),
                   detection_rate = mean(v > 0, na.rm = TRUE),
                   role = dplyr::if_else(g %in% sugar_all, "sugar", "defence"))
  })
}) %>%
  dplyr::mutate(passed = detection_rate >= detect_cutoff)

cat("[14] 检出率过滤（>= ", detect_cutoff, "）后各期可用gene：\n", sep = "")
print(as.data.frame(detect_tab %>% dplyr::filter(passed) %>%
                      dplyr::count(stage, role) %>%
                      tidyr::pivot_wider(names_from = role, values_from = n, values_fill = 0)))

cat("\n[14] 被过滤掉的gene（各期）：\n")
print(as.data.frame(detect_tab %>% dplyr::filter(!passed) %>%
                      dplyr::select(stage, gene, role, detection_rate) %>%
                      dplyr::arrange(stage, detection_rate)))

cor_tab <- purrr::map_dfr(stage_levels, function(st) {
  dd <- d %>% dplyr::filter(stage == st)
  ok_sugar <- detect_tab %>% dplyr::filter(stage == st, role == "sugar", passed) %>% dplyr::pull(gene)
  ok_def <- detect_tab %>% dplyr::filter(stage == st, role == "defence", passed) %>% dplyr::pull(gene)
  if (length(ok_sugar) == 0 || length(ok_def) == 0) return(tibble::tibble())
  purrr::map_dfr(ok_sugar, function(g1) {
    purrr::map_dfr(ok_def, function(g2) {
      ct <- suppressWarnings(stats::cor.test(dd[[g1]], dd[[g2]], method = "spearman",
                                             exact = FALSE))
      tibble::tibble(stage = st, sugar_gene = g1, defense_gene = g2, n = nrow(dd),
                     rho = unname(ct$estimate), p_value = ct$p.value)
    })
  })
}) %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(p_fdr = stats::p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(stage, p_value)

cat("\n[14] 每期组合数 / FDR < 0.05 / 名义 P < 0.05 / 最大 |rho|：\n")
print(as.data.frame(cor_tab %>%
                      dplyr::group_by(stage) %>%
                      dplyr::summarise(n_total = dplyr::n(),
                                       n_fdr005 = sum(p_fdr < 0.05, na.rm = TRUE),
                                       n_p005 = sum(p_value < 0.05, na.rm = TRUE),
                                       max_abs_rho = max(abs(rho), na.rm = TRUE),
                                       .groups = "drop")))

cat("\n[14] 各期最强 5 个组合：\n")
print(as.data.frame(cor_tab %>% dplyr::group_by(stage) %>%
                      dplyr::slice_min(p_value, n = 5, with_ties = FALSE) %>%
                      dplyr::select(stage, sugar_gene, defense_gene, rho, p_value, p_fdr) %>%
                      dplyr::ungroup()))

# ---- 热图（仅展示通过检出率过滤的gene）----
cor_tab %>%
  dplyr::mutate(
    stage = factor(stage, levels = stage_levels),
    sig = dplyr::case_when(p_fdr < 0.05 ~ "FDR < 0.05",
                           p_value < 0.05 ~ "P < 0.05",
                           TRUE ~ "n.s."),
    sig = factor(sig, levels = c("FDR < 0.05", "P < 0.05", "n.s."))
  ) -> d.cor

d.cor %>%
  ggplot2::ggplot(aes(x = defense_gene, y = sugar_gene, fill = rho)) +
  ggplot2::geom_tile() +
  ggplot2::geom_point(data = d.cor %>% dplyr::filter(sig != "n.s."),
                      aes(size = sig), shape = 1, color = "grey15") +
  ggplot2::facet_wrap(~stage, scales = "free") +
  ggplot2::scale_fill_gradient2(low = "#1e4e8c", mid = "#f2f2f0", high = "#a63446",
                                midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
  ggplot2::scale_size_manual(values = c("FDR < 0.05" = 2.6, "P < 0.05" = 1.8), name = NULL) +
  ggplot2::labs(x = "Defence-related genes", y = "High-sugar pathway genes",
                caption = paste0("Spearman correlations on sample-level expression; genes with detection rate >= ",
                                 detect_cutoff, " per stage; BH-FDR within each stage")) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, face = "italic"),
    axis.text.y = ggplot2::element_text(size = 8, face = "italic"),
    strip.text = ggplot2::element_text(face = "bold"),
    panel.grid = ggplot2::element_blank(),
    plot.caption = ggplot2::element_text(size = 7, hjust = 0)
  ) -> p.s5

bioRtools::write_data(p.s5, file = file.path(out_dir, "34_FigureS5_高糖and防御相关性热图_three_stages_robust_definition.pdf"),
                      width = 14, height = 5.5, units = "in")
bioRtools::write_data(p.s5, file = file.path(out_dir, "34_FigureS5_高糖and防御相关性热图_three_stages_robust_definition.png"),
                      width = 14, height = 5.5, units = "in", dpi = 500)
cat("\n[14] 输出：34_FigureS5_高糖and防御相关性热图_three_stages_robust_definition.pdf / .png\n")

bioRtools::write_data(list(correlation = cor_tab, detection = detect_tab), file = path_out)
cat("[14] 输出：", path_out, "\n")

cat("\n[14] sessionInfo()\n")
sessionInfo()
