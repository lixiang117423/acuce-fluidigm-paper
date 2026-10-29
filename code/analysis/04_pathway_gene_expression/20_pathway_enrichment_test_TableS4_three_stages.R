# ============================================================
# 04_pathway_gene_expression / 20_pathway_enrichment_test_TableS4_three_stages.R
# 作用：pathway水平富集test（Fisher 精确test）**统一成three_stages口径** —— 供 Table S4 使用
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_pathway_gene_expression/20_pathway_enrichment_test_TableS4_three_stages.R"
#
# and原两时期版（02_pathway_enrichment_test_TableS3S4.R）的差别：
#   ① 输入改用**three_stages limma（统一口径）**的 focal × 时期 results（and Figure 3A/4A/4B、Table S8 同源）；
#   ② 判据 A（主）＝pathway内上调gene（focal_logFC > 0）是否富集；
#      判据 B（辅）＝pathway内显著gene（focal_fdr < 0.05）是否富集；
#   ③ 覆盖 **3 个时期 × 4 个 focal × 2 个判据 × 20 条pathway**，FDR 按 (时期 × focal × 判据) correction；
#   ④ 输出写到 result/08_submission_format图（提交用表），原 result/04_pathway_gene_expression 只读不改。
#
# 输入：<root>/result/04_pathway_gene_expression/23_表_limma_unified_statistics_three_stages.xlsx（sheet gene_stage_focal_level）
# 输出：<root>/result/08_submission_format图/32_表_pathway富集_Fisher_three_stages.xlsx
#       （sheet：TableS4_focal87_filling、TableS4_focal69_tillering、all_stage_focal、summary）
# 其中 <root> = <repository_root>
#
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_in <- file.path(project_root, "result", "04_pathway_gene_expression", "23_表_limma_unified_statistics_three_stages.xlsx")
out_pub_dir <- file.path(project_root, "result", "08_submission_format图")
dir.create(out_pub_dir, showWarnings = FALSE, recursive = TRUE)
path_out <- file.path(out_pub_dir, "32_表_pathway富集_Fisher_three_stages.xlsx")

stopifnot(file.exists(path_in))

stage_levels <- c("tillering", "flowering", "filling")

bioRtools::read_data(path_in, sheet = "gene_stage_focal_level") -> focal_tab
stopifnot(nrow(focal_tab) > 0)

cat("[20] 输入：", nrow(focal_tab), "行；focal × 时期 组合数：",
    dplyr::n_distinct(paste(focal_tab$stage, focal_tab$focal)),
    "；时期：", paste(sort(unique(focal_tab$stage)), collapse = ", "), "\n")

enrich_one <- function(d, criterion) {
  d %>%
    dplyr::mutate(is_hit = if (criterion == "up") focal_logFC > 0 else focal_fdr < 0.05) -> dd

  dd %>%
    dplyr::group_by(process) %>%
    dplyr::summarise(in_hit = sum(is_hit, na.rm = TRUE), in_total = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      out_hit = sum(in_hit) - in_hit,
      out_total = sum(in_total) - in_total,
      in_miss = in_total - in_hit,
      out_miss = out_total - out_hit
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      fisher_p = stats::fisher.test(matrix(c(in_hit, in_miss, out_hit, out_miss), nrow = 2))$p.value,
      odds_ratio = stats::fisher.test(matrix(c(in_hit, in_miss, out_hit, out_miss), nrow = 2))$estimate
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(criterion = criterion)
}

res_tab <- focal_tab %>%
  dplyr::group_by(stage, focal) %>%
  dplyr::group_modify(~ dplyr::bind_rows(enrich_one(.x, "up"), enrich_one(.x, "sig"))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(stage, focal, criterion) %>%
  dplyr::mutate(fisher_fdr = p.adjust(fisher_p, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels)) %>%
  dplyr::arrange(stage, focal, criterion, fisher_p)

cat("\n[20] FDR < 0.05 的pathway（全部 focal × 时期）：\n")
res_tab %>%
  dplyr::filter(fisher_fdr < 0.05) %>%
  dplyr::select(stage, focal, criterion, process, in_hit, in_total, fisher_p, fisher_fdr) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n[20] 未correction P < 0.05 的pathway数：", sum(res_tab$fisher_p < 0.05),
    "；FDR < 0.05：", sum(res_tab$fisher_fdr < 0.05), "\n")

summary_tab <- res_tab %>%
  dplyr::group_by(stage, focal, criterion) %>%
  dplyr::summarise(n_pathway = dplyr::n(),
                   n_p05 = sum(fisher_p < 0.05),
                   n_fdr05 = sum(fisher_fdr < 0.05),
                   min_p = min(fisher_p), .groups = "drop")

bioRtools::write_data(
  list(
    TableS4_focal87_filling = res_tab %>% dplyr::filter(stage == "filling", focal == 87),
    TableS4_focal69_tillering = res_tab %>% dplyr::filter(stage == "tillering", focal == 69),
    all_stage_focal = as.data.frame(res_tab),
    summary = as.data.frame(summary_tab)
  ),
  file = path_out
)
cat("[20] 输出：", path_out, "\n")

cat("\n[20] sessionInfo()\n")
sessionInfo()
