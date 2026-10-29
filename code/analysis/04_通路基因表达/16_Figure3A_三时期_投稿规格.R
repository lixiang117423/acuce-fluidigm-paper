# ============================================================
# 04_通路基因表达 / 16_Figure3A_三时期_投稿规格.R
# 作用：**投稿规格版** Figure 3A（防御相关基因热图；统计口径与 15 号脚本一致）
#       与 15 号的差别在版面与结构：
#         ① 输出尺寸改为最终 rect 尺寸（17.0 × 13.0 cm，scale = 1）；
#         ② **生育期从行分面改为列分面**（原 process + stage ~ group 的结构会把行数
#            乘以 3 个生育期，实际 153 行，任何期刊尺寸都叠字）；
#         ③ 去掉面板内 caption（写进正文 Figure legend）；
#       输出目录改为 result/08_投稿规格图，原 result/04_通路基因表达 只读不改。
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_通路基因表达/16_Figure3A_三时期_投稿规格.R"
#
# 输入：<root>/result/03_基因表达总体/23_表_表达矩阵_三时期_不剔离群值.xlsx（combo_matrix_median、gene_annotation）
# 输出：<root>/result/08_投稿规格图/35_Figure3A_防御基因热图_三时期_中位数矩阵.pdf / .png
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 说明：3B（防御通路分数）与 4A/4B/4C（stage 层面 limma）都用**样本级**数据，不受组合均值口径影响，故沿用 25/26/27/28 号文件。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
out_dir <- file.path(project_root, "result", "04_通路基因表达")
out_pub_dir <- file.path(project_root, "result", "08_投稿规格图")  # 投稿规格面板输出目录
dir.create(out_pub_dir, showWarnings = FALSE, recursive = TRUE)
path_mat <- file.path(project_root, "result", "03_基因表达总体", "23_表_表达矩阵_三时期_不剔离群值.xlsx")

stopifnot(file.exists(path_mat))

stage_levels <- c("tillering", "flowering", "filling")
def_processes <- c("Defense", "Jasmonic acid", "Salicylic acid")

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_mat, "combo_matrix_median") -> mat.combo
bioRtools::read_data(path_mat, sheet = "gene_annotation") -> df.gene

meta_cols <- c("combo", "stage", "focal", "neighbor", "group")
gene_cols <- setdiff(colnames(mat.combo), meta_cols)

def_genes <- df.gene %>% dplyr::filter(process %in% def_processes) %>% dplyr::pull(gene)
cat("[15] 防御相关基因数：", length(def_genes), "\n")

mat.def <- as.matrix(mat.combo[, intersect(def_genes, gene_cols)])
rownames(mat.def) <- mat.combo$combo

scale01_by_row <- function(v) {
  rng <- range(v, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0, length(v)))
  (v - rng[1]) / diff(rng)
}
mat.def.scaled <- t(apply(mat.def, 1, scale01_by_row))

as.data.frame(mat.def.scaled) %>%
  tibble::rownames_to_column(var = "combo") %>%
  tidyr::pivot_longer(cols = -combo, names_to = "gene", values_to = "scaled_expression") %>%
  dplyr::left_join(mat.combo[, meta_cols], by = "combo") %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene") %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels),
                group = factor(group, levels = c("Pure", "Mix")),
                process = factor(process, levels = def_processes),
                gene = factor(gene, levels = rev(sort(unique(gene)))),
                combo_label = paste0(focal, "-", neighbor)) %>%
  ggplot2::ggplot(aes(x = combo_label, y = gene, fill = scaled_expression)) +
  ggplot2::geom_tile() +
  ggplot2::facet_grid(process ~ stage, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_gradientn(colours = c("#f2f2f0", "#a8bedb", "#1e4e8c"),
                                name = "Scaled median expression (0-1)") +
  ggplot2::labs(x = "Focal-neighbour combination", y = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
                 axis.text.y = ggplot2::element_text(size = 5.2, face = "italic"),
                 # 2026-10-19 修：x 轴标题 "Focal-neighbour combination" 左端越出纸边 2.6 pt（"Fo" 被切）→
                 #   字号 7.5 → 7，并给左/下各留 8/3 pt 边距
                 axis.title.x = ggplot2::element_text(size = 7),
                 plot.margin = ggplot2::margin(2, 2, 3, 8),
                 strip.text = ggplot2::element_text(size = 6, face = "bold"),
                 panel.grid = ggplot2::element_blank(),
                 legend.position = "right",
                 legend.key.width = ggplot2::unit(0.3, "cm"),
                 legend.title = ggplot2::element_text(size = 6, angle = 90),
                 legend.text = ggplot2::element_text(size = 5.5)) -> p.3a

bioRtools::write_data(p.3a, file = file.path(out_pub_dir, "35_Figure3A_防御基因热图_三时期_中位数矩阵.pdf"),
                      width = 17.0, height = 14.5, units = "cm")
bioRtools::write_data(p.3a, file = file.path(out_pub_dir, "35_Figure3A_防御基因热图_三时期_中位数矩阵.png"),
                      width = 17.0, height = 14.5, units = "cm", dpi = 600)
cat("[15] 输出：35_Figure3A_防御基因热图_三时期_中位数矩阵.pdf / .png\n")

cat("\n[15] sessionInfo()\n")
sessionInfo()
