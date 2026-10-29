# ============================================================
# 03_overall_gene_expression / 06_expression_matrix_build_three_stages.R
# 作用：按主题 03 现行口径（脚本 01）**重建three_stages**expression矩阵，供下游（PERMANOVA / PLS-DA / Figure 2）使用
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_overall_gene_expression/06_expression_matrix_build_three_stages.R"
#
# 输入：<root>/data/02_raw_data/04_Fluidigmexpression/00_引物and共享表/fluidigm.normaled.expression.3 times.xlsx
#       （three_stagesnormalized长表；time = 1st/2nd/3rd → tillering/flowering/filling）
# 输出：<root>/result/03_overall_gene_expression/11_表_expression矩阵_three_stages.xlsx
#       6 sheet：gene_annotation / sample_info / combo_matrix / sample_matrix / combo_long / outlier_summary
# 其中 <root> = <repository_root>
#
# 口径（and脚本 01 完全一致，仅时期数由 2 变 3）：
#   - 离群值：bioRtools::find_outliers()，按 stage × focal × neighbor × gene 分组判定；组合均值矩阵用剔除后data，
#     sample级矩阵用未剔除data（保证complete，供 limma 使用）；
#   - gene注释取自长表的 id / response / process / reference 列。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_in <- file.path(project_root, "data", "02_raw_data", "04_Fluidigmexpression",
                     "00_引物and共享表", "fluidigm.normaled.expression.3 times.xlsx")
out_dir <- file.path(project_root, "result", "03_overall_gene_expression")
path_out <- file.path(out_dir, "11_表_expression矩阵_three_stages.xlsx")

stopifnot(file.exists(path_in))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("[06] 输入：", path_in, "\n")

stage_map <- c("1st" = "tillering", "2nd" = "flowering", "3rd" = "filling")
stage_levels <- c("tillering", "flowering", "filling")

readxl::read_excel(path_in, .name_repair = "minimal") %>%
  as.data.frame() %>%
  dplyr::mutate(stage = factor(unname(stage_map[time]), levels = stage_levels),
                sample = as.character(sample)) -> d3

stopifnot(
  !anyNA(d3$stage),
  !anyNA(d3$sample), !anyNA(d3$focal), !anyNA(d3$neighbor),
  !anyNA(d3$gene), !anyNA(d3$expression)
)
if (anyNA(d3)) {
  cat("[06] 提示：非关键列存在缺失（各列缺失数）：\n")
  print(colSums(is.na(d3)))
}

cat("[06] 行数 ", nrow(d3), "；sample ", dplyr::n_distinct(d3$sample),
    "；gene ", dplyr::n_distinct(d3$gene), "；pathway ", dplyr::n_distinct(d3$process), "\n")
print(table(d3$stage))

# ---- gene注释表（一gene一行）----
# notes：three_stages长表里 OsWRKY82 有两行注释，唯一差别是 id 字段带一个换行残留 `_x000D_`，
# 其余字段（process / response / reference）完全相同；这里先去残留、再按 gene 去重，
# 保证下游 join 不产生 many-to-many 重复行（主题 08 的 gene_limma 曾因此多出 1 行，数值结论不受影响）。
d3 %>%
  dplyr::mutate(id = stringr::str_remove_all(id, "_x000D_")) %>%
  dplyr::select(gene, id, response, process, reference) %>%
  dplyr::distinct() %>%
  dplyr::arrange(gene, id, process) %>%
  dplyr::group_by(gene) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() -> df.gene
stopifnot(nrow(df.gene) == dplyr::n_distinct(d3$gene))
cat("[06] gene注释表：", nrow(df.gene), "行（一gene一行）\n")

# ---- 剔除离群值（stage × focal × neighbor × gene 内判定）----
d3 %>%
  dplyr::group_by(stage, focal, neighbor, gene) %>%
  dplyr::mutate(is.outlier = suppressWarnings(bioRtools::find_outliers(expression))) %>%
  dplyr::ungroup() -> df.mark

n_out <- sum(df.mark$is.outlier != "no")
cat("[06] 离群值标记数：", n_out, "（占 ", round(100 * n_out / nrow(df.mark), 2), "%）\n")

df.mark %>%
  dplyr::group_by(stage) %>%
  dplyr::summarise(n_total = dplyr::n(), n_outlier = sum(is.outlier != "no"),
                   pct_outlier = round(100 * sum(is.outlier != "no") / dplyr::n(), 2),
                   .groups = "drop") -> outlier_summary
print(as.data.frame(outlier_summary))

df.mark %>%
  dplyr::filter(is.outlier == "no") %>%
  dplyr::select(-is.outlier) -> df.clean

# ---- sample注释 ----
d3 %>%
  dplyr::distinct(stage, sample, focal, neighbor) %>%
  dplyr::mutate(group = dplyr::if_else(focal == neighbor, "Pure", "Mix")) -> df.sample.info

cat("[06] sample级：", nrow(df.sample.info), " 个sample×时期\n")
print(table(df.sample.info$stage, df.sample.info$group))

# ---- (a) 组合均值矩阵（48 行 = 3 期 × 16 组合）----
df.clean %>%
  dplyr::group_by(stage, focal, neighbor, gene) %>%
  dplyr::summarise(mean_expression = mean(expression), n_obs = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(combo = paste0(stage, "_", focal, "_", neighbor),
                group = dplyr::if_else(focal == neighbor, "Pure", "Mix")) %>%
  dplyr::select(combo, stage, focal, neighbor, group, gene, mean_expression, n_obs) -> df.combo.long

df.combo.long %>%
  dplyr::select(combo, stage, focal, neighbor, group, gene, mean_expression) %>%
  tidyr::pivot_wider(names_from = gene, values_from = mean_expression) -> mat.combo

gene_cols_combo <- setdiff(colnames(mat.combo), c("combo", "stage", "focal", "neighbor", "group"))
stopifnot(nrow(mat.combo) == 48, length(gene_cols_combo) == 92, !anyNA(mat.combo))

# ---- (b) sample级矩阵（254 行；未剔离群值，保证complete）----
d3 %>%
  dplyr::select(stage, sample, focal, neighbor, gene, expression) %>%
  tidyr::pivot_wider(names_from = gene, values_from = expression) -> mat.sample

gene_cols_sample <- setdiff(colnames(mat.sample), c("stage", "sample", "focal", "neighbor"))
stopifnot(nrow(mat.sample) == 254, length(gene_cols_sample) == 92,
          sum(is.na(mat.sample)) == 0)

cat("[06] 组合均值矩阵：", nrow(mat.combo), "行 ×", length(gene_cols_combo), "gene\n")
cat("[06] sample级矩阵：", nrow(mat.sample), "行 ×", length(gene_cols_sample), "gene；缺失 0\n")

bioRtools::write_data(
  list(
    gene_annotation = df.gene,
    sample_info = df.sample.info,
    combo_matrix = mat.combo,
    sample_matrix = mat.sample,
    combo_long = df.combo.long,
    outlier_summary = outlier_summary
  ),
  file = path_out
)

cat("\n[06] 输出：", path_out, "\n")

cat("\n[06] sessionInfo()\n")
sessionInfo()
