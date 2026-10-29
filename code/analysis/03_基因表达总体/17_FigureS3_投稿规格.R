# ============================================================
# 03_基因表达总体 / 04_VIP与遗传距离相关.R
# 作用：VIP 前 10 基因的表达水平与遗传距离的相关性（阴性结果）—— Figure S3
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_基因表达总体/04_VIP与遗传距离相关.R"
#
# 输入：<root>/result/03_基因表达总体/01_表_表达矩阵.xlsx（combo_matrix，minimal 列名修复）
#       <root>/result/03_基因表达总体/03_表_PLSDA.xlsx（vip_top10）
#       <root>/data/02_原始数据/05_遗传距离/月亮谷纯系的遗传距离.xlsx
# 输出：<root>/result/03_基因表达总体/04_表_VIP与遗传距离相关.xlsx
#       <root>/result/03_基因表达总体/10_FigureS3_VIP与遗传距离.pdf / .png
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 编号说明（2026-09-22 确认的方案）：原论文框架 Figure S3（遗传距离 × 表型）与 Figure S4（遗传距离 × VIP 基因）
#   内容重复，现统一为：Figure 1F = 遗传距离 × 表型（正文），**Figure S3 = 遗传距离 × VIP 前 10 基因**（本脚本）。
# 方法：仅用 12 个混种组合（Mix）；每个 VIP 基因在每个时期对遗传距离做 Pearson 相关（lm），
#   并按时期做 FDR（Benjamini–Hochberg）校正。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：单步耗时 < 60 s，每次运行都重新生成。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录见 code/03_基因表达总体/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_mat <- file.path(project_root, "result", "03_基因表达总体", "01_表_表达矩阵.xlsx")
path_plsda <- file.path(project_root, "result", "03_基因表达总体", "03_表_PLSDA.xlsx")
path_gd <- file.path(project_root, "data", "02_原始数据", "05_遗传距离", "月亮谷纯系的遗传距离.xlsx")
out_dir <- file.path(project_root, "result", "03_基因表达总体")
out_pub_dir <- file.path(project_root, "result", "08_投稿规格图")  # 投稿规格面板输出目录
dir.create(out_pub_dir, showWarnings = FALSE, recursive = TRUE)
path_out <- file.path(out_pub_dir, "98_重出副本_04_表_VIP与遗传距离相关.xlsx")  # 投稿规格版：表写到这里，不动原目录

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

# VIP 基因在 12 个混种组合中的表达
mat.combo %>%
  dplyr::filter(group == "Mix") %>%
  dplyr::select(combo, stage, focal, neighbor, dplyr::all_of(vip_top10$gene)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(vip_top10$gene), names_to = "gene", values_to = "expression") %>%
  dplyr::left_join(gd.long, by = c("focal", "neighbor")) -> df.vip

stopifnot(!anyNA(df.vip$genetic_distance))
cat("[04] 输入：", nrow(df.vip), " 行（", dplyr::n_distinct(df.vip$gene), " 基因 × ",
    dplyr::n_distinct(df.vip$stage), " 时期 × 12 混种组合）\n")

# 相关分析（每个基因 × 时期）
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

cat("\n[04] 相关分析结果（按 P 排序）：\n")
print(as.data.frame(cor_tab))

# ---- Figure S3 ----
cor_lab <- cor_tab %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::mutate(
    stage = factor(stage, levels = c("tillering", "filling")),
    label = paste0(gene, ": ", format_p(p_value))
  )

cat("\n[04] 未校正 P < 0.05 的组合数：", nrow(cor_lab), "；FDR < 0.05 的组合数：",
    sum(cor_tab$p_fdr < 0.05), "\n")

df.vip %>%
  dplyr::mutate(
    stage = factor(stage, levels = c("tillering", "filling")),
    gene = factor(gene, levels = vip_top10$gene)
  ) %>%
  ggplot2::ggplot(aes(x = genetic_distance, y = expression)) +
  ggplot2::geom_point(size = 1.6, color = "#1e4e8c", alpha = 0.85) +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "grey35", linewidth = 0.6) +
  # 投稿规格：原为 facet_grid(stage ~ gene)（20 个小分面，174 mm 宽时每个分面仅 ~17 mm，
  # 轴标签只剩 2.8-4.1 pt）；改为 5 列 × 4 行分面，每格 ~35 mm，字号恢复到 7 pt 以上
  ggplot2::facet_wrap(~ stage + gene, ncol = 5, nrow = 4, scales = "free_y") +
  ggplot2::labs(x = "Genetic distance", y = "Relative expression (2^-ΔCt)") +
  bioRtools::theme_prism() +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 7, face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 6.5),
    axis.text.y = ggplot2::element_text(size = 6.5)
  ) -> p.s3

bioRtools::write_data(p.s3, file = file.path(out_pub_dir, "10_FigureS3_VIP与遗传距离.pdf"),
                      width = 17.4, height = 15.0, units = "cm")
bioRtools::write_data(p.s3, file = file.path(out_pub_dir, "10_FigureS3_VIP与遗传距离.png"),
                      width = 17.4, height = 15.0, units = "cm", dpi = 600)
cat("[04] 输出：10_FigureS3_VIP与遗传距离.pdf / .png\n")

bioRtools::write_data(
  list(vip_gene_expression = df.vip, correlation = cor_tab),
  file = path_out
)
cat("[04] 输出：", path_out, "\n")

cat("\n[04] sessionInfo()\n")
sessionInfo()
