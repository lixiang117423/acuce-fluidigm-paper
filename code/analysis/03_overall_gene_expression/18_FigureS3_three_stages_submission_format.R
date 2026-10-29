# ============================================================
# 03_overall_gene_expression / 18_FigureS3_three_stages_submission_format.R
# 作用：Figure S3（VIP 前 10 geneexpression × 遗传距离的相关性，阴性results）**统一成three_stages口径**，
#       并按final版面尺寸（submission_format）出图
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_overall_gene_expression/18_FigureS3_three_stages_submission_format.R"
#
# and原两时期版（04_VIP_and_genetic_distance_correlation.R）的差别：
#   ① expression矩阵改用**three_stages、without_outlier_removal**的组合median_matrix（and Figure 2A/3A 同源）；
#   ② VIP gene改用**three_stagesmedian_matrix的 PLS-DA**（and Figure 2E 同一模型）；
#   ③ 相关analysis从 2 个时期 × 10 gene = 20 个test → **3 个时期 × 10 gene = 30 个test**，FDR 按时期correction；
#   ④ 出图尺寸改为final rect 尺寸（17.4 × 22.0 cm，scale = 1），分面 5 列 × 6 行保证字号；
#   ⑤ 输出写到 result/08_submission_format图，原 result/03_overall_gene_expression 只读不改。
#
# 输入：<root>/result/03_overall_gene_expression/23_表_expression矩阵_three_stages_without_outlier_removal.xlsx（combo_matrix_median）
#       <root>/result/03_overall_gene_expression/24_表_PLSDA_three_stages_median_matrix.xlsx（vip_top10）
#       <root>/data/02_raw_data/05_遗传距离/月亮谷纯系的遗传距离.xlsx
# 输出：<root>/result/08_submission_format图/04b_表_VIP_and_genetic_distance_correlation_three_stages.xlsx
#       <root>/result/08_submission_format图/10b_FigureS3_VIPand遗传距离_three_stages_submission_format.pdf / .png
# 其中 <root> = <repository_root>
#
# 方法：仅用 12 个混种组合（Mix）；每个 VIP gene在每个时期对遗传距离做 Pearson 相关（lm），
#   并按时期做 FDR（Benjamini–Hochberg）correction。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
out_dir <- file.path(project_root, "result", "03_overall_gene_expression")
out_pub_dir <- file.path(project_root, "result", "08_submission_format图")  # submission_format输出目录
dir.create(out_pub_dir, showWarnings = FALSE, recursive = TRUE)
path_mat <- file.path(out_dir, "23_表_expression矩阵_three_stages_without_outlier_removal.xlsx")
path_plsda <- file.path(out_dir, "24_表_PLSDA_three_stages_median_matrix.xlsx")
path_gd <- file.path(project_root, "data", "02_raw_data", "05_遗传距离", "月亮谷纯系的遗传距离.xlsx")
path_out <- file.path(out_pub_dir, "04b_表_VIP_and_genetic_distance_correlation_three_stages.xlsx")

stopifnot(file.exists(path_mat), file.exists(path_plsda), file.exists(path_gd))

stage_levels <- c("tillering", "flowering", "filling")

format_p <- function(p, digits = 3) {
  dplyr::if_else(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, format = "f", digits = digits)))
}

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_mat, "combo_matrix_median") -> mat.combo
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

# VIP gene在 12 个混种组合中的expression（three_stages）
mat.combo %>%
  dplyr::filter(group == "Mix") %>%
  dplyr::select(combo, stage, focal, neighbor, dplyr::all_of(vip_top10$gene)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(vip_top10$gene), names_to = "gene", values_to = "expression") %>%
  dplyr::left_join(gd.long, by = c("focal", "neighbor")) -> df.vip

stopifnot(!anyNA(df.vip$genetic_distance), !anyNA(df.vip$expression))
cat("[18] 输入：", nrow(df.vip), " 行（", dplyr::n_distinct(df.vip$gene), " gene × ",
    dplyr::n_distinct(df.vip$stage), " 时期 × 12 混种组合）\n", sep = "")

# 相关analysis（每个gene × 时期），FDR 按时期correction
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

cat("\n[18] 相关analysisresults（按 P 排序，前 12 行）：\n")
print(utils::head(as.data.frame(cor_tab), 12), row.names = FALSE)

n_tests <- nrow(cor_tab)
n_nom <- sum(cor_tab$p_value < 0.05)
n_fdr <- sum(cor_tab$p_fdr < 0.05)
cat("\n[18] test总数 ", n_tests, "；未correction P < 0.05：", n_nom, "；FDR < 0.05：", n_fdr, "\n", sep = "")

# ---- Figure S3（three_stages，5 列 × 6 行分面，submission_format）----
df.vip %>%
  dplyr::mutate(
    stage = factor(stage, levels = stage_levels),
    gene = factor(gene, levels = vip_top10$gene)
  ) %>%
  ggplot2::ggplot(aes(x = genetic_distance, y = expression)) +
  ggplot2::geom_point(size = 1.5, color = "#1e4e8c", alpha = 0.85) +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "grey35", linewidth = 0.5) +
  ggplot2::facet_wrap(~ stage + gene, ncol = 5, nrow = 6, scales = "free_y") +
  ggplot2::labs(x = "Genetic distance", y = "Relative expression (2^-dCt)") +
  bioRtools::theme_prism() +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 6.5, face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 6.5),
    axis.text.y = ggplot2::element_text(size = 6.5),
    axis.title = ggplot2::element_text(size = 7.5)
  ) -> p.s3

bioRtools::write_data(p.s3, file = file.path(out_pub_dir, "10b_FigureS3_VIPand遗传距离_three_stages_submission_format.pdf"),
                      width = 17.4, height = 22.0, units = "cm")
bioRtools::write_data(p.s3, file = file.path(out_pub_dir, "10b_FigureS3_VIPand遗传距离_three_stages_submission_format.png"),
                      width = 17.4, height = 22.0, units = "cm", dpi = 600)
cat("[18] 输出：10b_FigureS3_VIPand遗传距离_three_stages_submission_format.pdf / .png\n")

bioRtools::write_data(
  list(correlation = cor_tab, vip_top10 = vip_top10, data_joined = df.vip,
       summary = tibble::tibble(n_tests = n_tests, n_nominal_p05 = n_nom, n_fdr05 = n_fdr,
                                stages = paste(stage_levels, collapse = ", "))),
  file = path_out
)
cat("[18] 输出：", path_out, "\n")

cat("\n[18] sessionInfo()\n")
sessionInfo()
