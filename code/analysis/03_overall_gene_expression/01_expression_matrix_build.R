# ============================================================
# 03_overall_gene_expression / 01_expression_matrix_build.R
# 作用：从 Fluidigm normalizedexpression表build两套analysis矩阵，供 PLS-DA、PERMANOVA and下游plotting复用
#       （a）组合均值矩阵：2 时期 × 16 组合 = 32 个sample（4 个纯种 + 12 个混种组合；and旧版 PLS-DA 输入一致）
#       （b）sample级矩阵：168 个sample（用于 PERMANOVA andgene水平test）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_overall_gene_expression/01_expression_matrix_build.R"
#
# 输入：<root>/data/processed_inputs/04_fluidigm_normalized_expression.xlsx（15456 行 × 10 列）
# 输出：<root>/result/03_overall_gene_expression/01_表_expression矩阵.xlsx（多 sheet）
# 其中 <root> = <repository_root>
#
# data口径（实测核对 2026-09-22）：
#   1) expression 为线性尺度 2^(−ΔCt)（内参 = OseEF-1 / OseIF-4a / OsUBQ5 的 Ct 几何均值），
#      **不是** z-score；需要 0–1 缩放或 Z-score 时在各图/各test里现做；
#   2) 表内 92 个gene / 20 条pathway；引物面板共 96 个位点 = 92 个目标gene + 4 个内参
#      （OsAct / OseEF-1α / OseIF-4a / OsUBQ5，其中 OsAct 在preparestage已剔除），
#      即paper里"96 个gene / 21 条pathway"应表述为"96 个检测位点（92 个gene + 4 个内参）/ 20 条pathway（+ HKG）"；
#   3) Fluidigm 覆盖 4 个焦点纯系的全部 16 个 focal × neighbor 组合（含 4 个纯种），每时期 16 个组合，合计 32 个；
#   4) 每个 (stage, sample, gene) 只有 1 个值，故sample级矩阵在剔除离群值后会出现 NA，
#      sample级矩阵一律保存"未剔离群值"的completeversion；组合均值矩阵保存"已剔离群值"version（and旧版一致）。
# 离群值处理：bioRtools::find_outliers()，按 stage × focal × neighbor × gene 分组判定。
#
# 【重要】矩阵 sheet 的读取方式：readxl/rbioRtools 的默认列名修复（.name_repair = "universal"）会把
#   gene名中的连字符改成点号（OsD14-1 → OsD14.1、OsCPS4-1 → OsCPS4.1、OsNADP-isocitrate → OsNADP.isocitrate），
#   因此**读取本file的两个矩阵 sheet 时必须用 readxl::read_excel(..., .name_repair = "minimal")**，
#   不要用 bioRtools::read_data()（会破坏gene名，导致下游按gene名匹配失败）。
#
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：单步耗时 < 30 s，每次运行都重新generate。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/03_overall_gene_expression/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_in <- file.path(project_root, "data", "processed_inputs", "04_fluidigm_normalized_expression.xlsx")
out_dir <- file.path(project_root, "result", "03_overall_gene_expression")
path_out <- file.path(out_dir, "01_表_expression矩阵.xlsx")

stopifnot(file.exists(path_in))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("[01] 输入：", path_in, "\n")

bioRtools::read_data(path_in) -> d4
stopifnot(nrow(d4) == 15456, ncol(d4) == 10)

cat("[01] 行数 ", nrow(d4), "；gene数 ", dplyr::n_distinct(d4$gene),
    "；pathway数 ", dplyr::n_distinct(d4$process), "\n")

# gene注释表（gene → id / process / response / reference），供下游复用
d4 %>%
  dplyr::select(gene, id, process, response, reference) %>%
  dplyr::distinct() -> df.gene

stopifnot(nrow(df.gene) == dplyr::n_distinct(d4$gene))

# ---- 剔除离群值（stage × focal × neighbor × gene 内判定）----
d4 %>%
  dplyr::group_by(stage, focal, neighbor, gene) %>%
  dplyr::mutate(is.outlier = suppressWarnings(bioRtools::find_outliers(expression))) %>%
  dplyr::ungroup() -> df.mark

n_out <- sum(df.mark$is.outlier != "no")
cat("[01] 离群值标记数：", n_out, "（占 ", round(100 * n_out / nrow(df.mark), 2), "%）\n")

df.mark %>%
  dplyr::filter(is.outlier == "no") %>%
  dplyr::select(-is.outlier) -> df.clean

# ---- sample注释 ----
df.clean %>%
  dplyr::distinct(stage, sample, focal, neighbor) %>%
  dplyr::mutate(group = dplyr::if_else(focal == neighbor, "Pure", "Mix")) -> df.sample.info

cat("[01] sample级：", nrow(df.sample.info), " 个sample×时期\n")
print(table(df.sample.info$stage, df.sample.info$group))

# ---- (a) 组合均值矩阵（40 行）----
df.clean %>%
  dplyr::group_by(stage, focal, neighbor, gene) %>%
  dplyr::summarise(mean_expression = mean(expression), n_obs = dplyr::n(), .groups = "drop") -> df.combo.long

df.combo.long %>%
  dplyr::mutate(
    combo = paste0(stage, "_", focal, "_", neighbor),
    group = dplyr::if_else(focal == neighbor, "Pure", "Mix")
  ) %>%
  dplyr::select(combo, stage, focal, neighbor, group, gene, mean_expression) -> df.combo

df.combo %>%
  dplyr::select(combo, stage, focal, neighbor, group, gene, mean_expression) %>%
  tidyr::pivot_wider(names_from = gene, values_from = mean_expression) -> mat.combo

gene_cols_combo <- setdiff(colnames(mat.combo), c("combo", "stage", "focal", "neighbor", "group"))
stopifnot(nrow(mat.combo) == 32, length(gene_cols_combo) == 92, !anyNA(mat.combo))

# ---- (b) sample级矩阵（168 行，未剔离群值，保证complete）----
d4 %>%
  dplyr::select(stage, sample, focal, neighbor, gene, expression) %>%
  tidyr::pivot_wider(names_from = gene, values_from = expression) -> mat.sample

gene_cols_sample <- setdiff(colnames(mat.sample), c("stage", "sample", "focal", "neighbor"))
stopifnot(nrow(mat.sample) == 168, length(gene_cols_sample) == 92)
cat("[01] sample级矩阵缺失单元格（未剔离群值时应为 0）：", sum(is.na(mat.sample)), "\n")

cat("\n[01] 组合均值矩阵：", nrow(mat.combo), "行 ×", length(gene_cols_combo), "gene\n")
cat("[01] sample级矩阵：", nrow(mat.sample), "行 ×", length(gene_cols_sample), "gene\n")

bioRtools::write_data(
  list(
    gene_annotation = df.gene,
    sample_info = df.sample.info,
    combo_matrix = mat.combo,
    sample_matrix = mat.sample,
    combo_long = df.combo
  ),
  file = path_out
)
cat("[01] 输出：", path_out, "\n")

cat("\n[01] sessionInfo()\n")
sessionInfo()
