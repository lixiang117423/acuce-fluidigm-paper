# ============================================================
# 04_通路基因表达 / 01_通路基因统计.R
# 作用：基因水平的 Mix vs Pure 差异（log2FC + Wilcoxon/Welch + FDR），按 时期 与 时期 × focal 两种粒度输出
#       —— Figure 3A–3D、Figure 4A–4C 与 Table S5 的统计输入
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_通路基因表达/01_通路基因统计.R"
#
# 输入：<root>/data/01_分析输入表/04_fluidigm标准化表达量.xlsx
# 输出：<root>/result/04_通路基因表达/01_表_通路基因差异.xlsx
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 口径说明：
#   1) expression 为 2^(−ΔCt) 线性尺度；Fold change 用均值比（Mix/Pure），log2FC = log2(mean_Mix / mean_Pure)；
#   2) 离群值按 stage × focal × neighbor × gene 剔除（与主题 03 一致）；
#   3) 两种粒度：① stage 层面（把 4 个 focal 的样本合并，与旧版 Figure 4A/4C 的做法一致）；
#      ② stage × focal 层面（用于 Figure 3C/3D 的焦点纯系比较）；
#   4) 检验用 Wilcoxon 秩和（主，旧版亦用）与 Welch t（对照），FDR 用 Benjamini–Hochberg，按检验族分别校正
#      （stage 层面按 stage 校正；stage × focal 层面按 stage × focal 校正）。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：单步耗时 < 60 s，每次运行都重新生成。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录见 code/04_通路基因表达/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_in <- file.path(project_root, "data", "01_分析输入表", "04_fluidigm标准化表达量.xlsx")
out_dir <- file.path(project_root, "result", "04_通路基因表达")
path_out <- file.path(out_dir, "01_表_通路基因差异.xlsx")

stopifnot(file.exists(path_in))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

bioRtools::read_data(path_in) -> d4
stopifnot(nrow(d4) == 15456)

# 基因注释（去重）
d4 %>%
  dplyr::select(gene, id, process) %>%
  dplyr::distinct() -> df.gene

# 剔除离群值
d4 %>%
  dplyr::group_by(stage, focal, neighbor, gene) %>%
  dplyr::mutate(is.outlier = suppressWarnings(bioRtools::find_outliers(expression))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(is.outlier == "no") %>%
  dplyr::select(-is.outlier) %>%
  dplyr::mutate(group = dplyr::if_else(focal == neighbor, "Pure", "Mix")) -> df.clean

cat("[01] 清洗后数据：", nrow(df.clean), "行\n")

# ---- 通用：按给定分组列做基因水平检验 ----
test_by <- function(data, group_cols) {
  data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, "gene")))) %>%
    dplyr::group_modify(~ {
      mix <- .x$expression[.x$group == "Mix"]
      pur <- .x$expression[.x$group == "Pure"]
      if (length(mix) < 3 || length(pur) < 3) {
        return(tibble::tibble(
          n_mix = length(mix), n_pure = length(pur),
          mean_mix = NA_real_, mean_pure = NA_real_, log2FC = NA_real_,
          wilcoxon_p = NA_real_, welch_p = NA_real_, cohens_d = NA_real_
        ))
      }
      ww <- suppressWarnings(stats::wilcox.test(mix, pur))
      tt <- stats::t.test(mix, pur)
      dd <- (mean(mix) - mean(pur)) /
        sqrt((stats::var(pur) + stats::var(mix)) / 2)
      tibble::tibble(
        n_mix = length(mix),
        n_pure = length(pur),
        mean_mix = mean(mix),
        mean_pure = mean(pur),
        log2FC = log2(mean(mix) / mean(pur)),
        wilcoxon_p = ww$p.value,
        welch_p = tt$p.value,
        cohens_d = dd
      )
    }) %>%
    dplyr::ungroup()
}

stage_tab <- test_by(df.clean, "stage") %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(
    wilcoxon_fdr = p.adjust(wilcoxon_p, method = "BH"),
    welch_fdr = p.adjust(welch_p, method = "BH")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(df.gene, by = "gene") %>%
  dplyr::select(stage, gene, id, process, dplyr::everything())

focal_tab <- test_by(df.clean, c("stage", "focal")) %>%
  dplyr::group_by(stage, focal) %>%
  dplyr::mutate(
    wilcoxon_fdr = p.adjust(wilcoxon_p, method = "BH"),
    welch_fdr = p.adjust(welch_p, method = "BH")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(df.gene, by = "gene") %>%
  dplyr::select(stage, focal, gene, id, process, dplyr::everything())

cat("\n[01] stage 层面显著基因数（Wilcoxon FDR < 0.05）：\n")
print(as.data.frame(stage_tab %>% dplyr::filter(wilcoxon_fdr < 0.05) %>%
                      dplyr::count(stage, process, sort = TRUE)))

cat("\n[01] stage × focal 层面（69 分蘖期 / 87 灌浆期）BR 与 AUX/SA/Drought 基因：\n")
focal_tab %>%
  dplyr::filter((stage == "tillering" & focal == 69 & process %in% c("Brassinosteroid", "Auxin")) |
                  (stage == "filling" & focal == 87 &
                     process %in% c("Brassinosteroid", "Salicylic acid", "Drought"))) %>%
  dplyr::select(stage, focal, process, gene, mean_mix, mean_pure, log2FC, wilcoxon_p, wilcoxon_fdr) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ---- 供下游复用的宽表（组合均值，含元信息）----
df.clean %>%
  dplyr::group_by(stage, focal, neighbor, gene) %>%
  dplyr::summarise(mean_expression = mean(expression), .groups = "drop") %>%
  dplyr::mutate(
    group = dplyr::if_else(focal == neighbor, "Pure", "Mix"),
    combo = paste0(stage, "_", focal, "_", neighbor)
  ) -> combo_tab

bioRtools::write_data(
  list(
    stage_level = stage_tab,
    stage_focal_level = focal_tab,
    combo_means = combo_tab,
    clean_long = df.clean,
    gene_annotation = df.gene
  ),
  file = path_out
)
cat("[01] 输出：", path_out, "\n")

cat("\n[01] sessionInfo()\n")
sessionInfo()
