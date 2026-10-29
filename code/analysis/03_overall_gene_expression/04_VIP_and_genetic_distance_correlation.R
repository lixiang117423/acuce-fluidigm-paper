# ============================================================
# 03_overall_gene_expression / 04_VIP_and_genetic_distance_correlation.R
# 作用：VIP 前 10 gene的expression水平and遗传距离的相关性（阴性results）—— Figure S3
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_overall_gene_expression/04_VIP_and_genetic_distance_correlation.R"
#
# 输入：<root>/result/03_overall_gene_expression/01_表_expression矩阵.xlsx（combo_matrix，minimal 列名修复）
#       <root>/result/03_overall_gene_expression/03_表_PLSDA.xlsx（vip_top10）
#       <root>/data/02_raw_data/05_遗传距离/月亮谷纯系的遗传距离.xlsx
# 输出：<root>/result/03_overall_gene_expression/04_表_VIP_and_genetic_distance_correlation.xlsx
#       <root>/result/03_overall_gene_expression/10_FigureS3_VIPand遗传距离.pdf / .png
# 其中 <root> = <repository_root>
#
# 编号notes（2026-09-22 确认的方案）：原paper框架 Figure S3（遗传距离 × 表型）and Figure S4（遗传距离 × VIP gene）
#   内容重复，现统一为：Figure 1F = 遗传距离 × 表型（正文），**Figure S3 = 遗传距离 × VIP 前 10 gene**（本脚本）。
# 方法：仅用 12 个混种组合（Mix）；每个 VIP gene在每个时期对遗传距离做 Pearson 相关（lm），
#   并按时期做 FDR（Benjamini–Hochberg）correction。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：单步耗时 < 60 s，每次运行都重新generate。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/03_overall_gene_expression/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_mat <- file.path(project_root, "result", "03_overall_gene_expression", "01_表_expression矩阵.xlsx")
path_plsda <- file.path(project_root, "result", "03_overall_gene_expression", "03_表_PLSDA.xlsx")
path_gd <- file.path(project_root, "data", "02_raw_data", "05_遗传距离", "月亮谷纯系的遗传距离.xlsx")
out_dir <- file.path(project_root, "result", "03_overall_gene_expression")
path_out <- file.path(out_dir, "04_表_VIP_and_genetic_distance_correlation.xlsx")

format_p <- function(p, digits = 3) {
  dplyr::if_else(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, format = "f", digits = digits)))
}

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_mat, "combo_matrix") -> mat.combo
bioRtools::read_data(path_plsda, sheet = "vip_top10") -> vip_top10
bioRtools::read_data(path_gd) -> gd

focus <- c("15", "49", "69", "87")

# 遗传距离子矩阵 → 长表
gd %>%
  dplyr::select(id, dplyr::all_of(paste0("A", focus))) %>%
  dplyr::filter(id %in% paste0("A", focus)) %>%
  tibble::column_to_rownames(var = "id") -> gd.sub

gd.sub %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "focal") %>%
  tidyr::pivot_longer(cols = -focal, names_to = "neighbor", values_to = "genetic_distance") %>%
  dplyr::mutate(
    focal = stringr::str_remove(focal, "^A"),
    neighbor = stringr::str_remove(neighbor, "^A")
  ) -> gd.long

# VIP gene在 12 个混种组合中的expression
mat.combo %>%
  dplyr::filter(group == "Mix") %>%
  dplyr::select(combo, stage, focal, neighbor, dplyr::all_of(vip_top10$gene)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(vip_top10$gene), names_to = "gene", values_to = "expression") %>%
  dplyr::left_join(gd.long, by = c("focal", "neighbor")) -> df.vip

stopifnot(!anyNA(df.vip$genetic_distance))
cat("[04] 输入：", nrow(df.vip), " 行（", dplyr::n_distinct(df.vip$gene), " gene × ",
    dplyr::n_distinct(df.vip$stage), " 时期 × 12 混种组合）\n")

# 相关analysis（每个gene × 时期）
cor_tab <- df.vip %>%
  dplyr::group_by(stage, gene) %>%
  dplyr::group_modify(~ {
    m <- stats::lm(expression ~ genetic_distance, data = .x)
    s <- summary(m)
    tibble::tibble(
      n = nrow(.x),
      slope = unname(stats::coef(m)[2]),
      r = unname(sqrt(s$r.squared)) * sign(unname(stats::coef(m)[2])),
      r_squared = s$r.squared,
      p_value = s$coefficients[2, 4]
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(p_value)

cat("\n[04] 相关analysisresults（按 P 排序）：\n")
print(as.data.frame(cor_tab))

# ---- Figure S3 ----
cor_lab <- cor_tab %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::mutate(
    stage = factor(stage, levels = c("tillering", "filling")),
    label = paste0(gene, ": ", format_p(p_value))
  )

cat("\n[04] 未correction P < 0.05 的组合数：", nrow(cor_lab), "；FDR < 0.05 的组合数：",
    sum(cor_tab$p_fdr < 0.05), "\n")

df.vip %>%
  dplyr::mutate(
    stage = factor(stage, levels = c("tillering", "filling")),
    gene = factor(gene, levels = vip_top10$gene)
  ) %>%
  ggplot2::ggplot(aes(x = genetic_distance, y = expression)) +
  ggplot2::geom_point(size = 1.6, color = "#1e4e8c", alpha = 0.85) +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "grey35", linewidth = 0.6) +
  ggplot2::facet_grid(stage ~ gene, scales = "free_y") +
  ggplot2::labs(x = "Genetic distance", y = "Relative expression (2^-ΔCt)") +
  bioRtools::theme_prism() +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 7, face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y = ggplot2::element_text(size = 7)
  ) -> p.s3

bioRtools::write_data(p.s3, file = file.path(out_dir, "10_FigureS3_VIPand遗传距离.pdf"),
                      width = 14, height = 6, units = "in")
bioRtools::write_data(p.s3, file = file.path(out_dir, "10_FigureS3_VIPand遗传距离.png"),
                      width = 14, height = 6, units = "in", dpi = 500)
cat("[04] 输出：10_FigureS3_VIPand遗传距离.pdf / .png\n")

bioRtools::write_data(
  list(vip_gene_expression = df.vip, correlation = cor_tab),
  file = path_out
)
cat("[04] 输出：", path_out, "\n")

cat("\n[04] sessionInfo()\n")
sessionInfo()
