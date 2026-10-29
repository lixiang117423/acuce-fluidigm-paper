# ============================================================
# 03_基因表达总体 / 18_FigureS3_三时期_投稿规格.R
# 作用：Figure S3（VIP 前 10 基因表达 × 遗传距离的相关性，阴性结果）**统一成三时期口径**，
#       并按最终版面尺寸（投稿规格）出图
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_基因表达总体/18_FigureS3_三时期_投稿规格.R"
#
# 与原两时期版（04_VIP与遗传距离相关.R）的差别：
#   ① 表达矩阵改用**三时期、不剔离群值**的组合中位数矩阵（与 Figure 2A/3A 同源）；
#   ② VIP 基因改用**三时期中位数矩阵的 PLS-DA**（与 Figure 2E 同一模型）；
#   ③ 相关分析从 2 个时期 × 10 基因 = 20 个检验 → **3 个时期 × 10 基因 = 30 个检验**，FDR 按时期校正；
#   ④ 出图尺寸改为最终 rect 尺寸（17.4 × 22.0 cm，scale = 1），分面 5 列 × 6 行保证字号；
#   ⑤ 输出写到 result/08_投稿规格图，原 result/03_基因表达总体 只读不改。
#
# 输入：<root>/result/03_基因表达总体/23_表_表达矩阵_三时期_不剔离群值.xlsx（combo_matrix_median）
#       <root>/result/03_基因表达总体/24_表_PLSDA_三时期_中位数矩阵.xlsx（vip_top10）
#       <root>/data/02_原始数据/05_遗传距离/月亮谷纯系的遗传距离.xlsx
# 输出：<root>/result/08_投稿规格图/04b_表_VIP与遗传距离相关_三时期.xlsx
#       <root>/result/08_投稿规格图/10b_FigureS3_VIP与遗传距离_三时期_投稿规格.pdf / .png
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 方法：仅用 12 个混种组合（Mix）；每个 VIP 基因在每个时期对遗传距离做 Pearson 相关（lm），
#   并按时期做 FDR（Benjamini–Hochberg）校正。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
out_dir <- file.path(project_root, "result", "03_基因表达总体")
out_pub_dir <- file.path(project_root, "result", "08_投稿规格图")  # 投稿规格输出目录
dir.create(out_pub_dir, showWarnings = FALSE, recursive = TRUE)
path_mat <- file.path(out_dir, "23_表_表达矩阵_三时期_不剔离群值.xlsx")
path_plsda <- file.path(out_dir, "24_表_PLSDA_三时期_中位数矩阵.xlsx")
path_gd <- file.path(project_root, "data", "02_原始数据", "05_遗传距离", "月亮谷纯系的遗传距离.xlsx")
path_out <- file.path(out_pub_dir, "04b_表_VIP与遗传距离相关_三时期.xlsx")

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

# VIP 基因在 12 个混种组合中的表达（三时期）
mat.combo %>%
  dplyr::filter(group == "Mix") %>%
  dplyr::select(combo, stage, focal, neighbor, dplyr::all_of(vip_top10$gene)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(vip_top10$gene), names_to = "gene", values_to = "expression") %>%
  dplyr::left_join(gd.long, by = c("focal", "neighbor")) -> df.vip

stopifnot(!anyNA(df.vip$genetic_distance), !anyNA(df.vip$expression))
cat("[18] 输入：", nrow(df.vip), " 行（", dplyr::n_distinct(df.vip$gene), " 基因 × ",
    dplyr::n_distinct(df.vip$stage), " 时期 × 12 混种组合）\n", sep = "")

# 相关分析（每个基因 × 时期），FDR 按时期校正
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

cat("\n[18] 相关分析结果（按 P 排序，前 12 行）：\n")
print(utils::head(as.data.frame(cor_tab), 12), row.names = FALSE)

n_tests <- nrow(cor_tab)
n_nom <- sum(cor_tab$p_value < 0.05)
n_fdr <- sum(cor_tab$p_fdr < 0.05)
cat("\n[18] 检验总数 ", n_tests, "；未校正 P < 0.05：", n_nom, "；FDR < 0.05：", n_fdr, "\n", sep = "")

# ---- Figure S3（三时期，5 列 × 6 行分面，投稿规格）----
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

bioRtools::write_data(p.s3, file = file.path(out_pub_dir, "10b_FigureS3_VIP与遗传距离_三时期_投稿规格.pdf"),
                      width = 17.4, height = 22.0, units = "cm")
bioRtools::write_data(p.s3, file = file.path(out_pub_dir, "10b_FigureS3_VIP与遗传距离_三时期_投稿规格.png"),
                      width = 17.4, height = 22.0, units = "cm", dpi = 600)
cat("[18] 输出：10b_FigureS3_VIP与遗传距离_三时期_投稿规格.pdf / .png\n")

bioRtools::write_data(
  list(correlation = cor_tab, vip_top10 = vip_top10, data_joined = df.vip,
       summary = tibble::tibble(n_tests = n_tests, n_nominal_p05 = n_nom, n_fdr05 = n_fdr,
                                stages = paste(stage_levels, collapse = ", "))),
  file = path_out
)
cat("[18] 输出：", path_out, "\n")

cat("\n[18] sessionInfo()\n")
sessionInfo()
