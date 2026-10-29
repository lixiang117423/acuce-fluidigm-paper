# ============================================================
# 04_pathway_gene_expression / 02_pathway_enrichment_test_TableS3S4.R
# 作用：pathway水平富集test（Fisher 精确test）—— Table S3（69 分蘖期）and Table S4（87 灌浆期）
#       并为所有 focal × 时期 组合提供同一套results
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_pathway_gene_expression/02_pathway_enrichment_test_TableS3S4.R"
#
# 输入：<root>/result/04_pathway_gene_expression/01_表_pathwaygene差异.xlsx（sheet stage_focal_level）
# 输出：<root>/result/04_pathway_gene_expression/02_表_pathway富集_TableS3S4.xlsx
# 其中 <root> = <repository_root>
#
# 重要notes：旧版脚本**从未做过卡方/富集test**（两份 qmd 全文都没有 chisq/fisher/p.adjust 调用），
#   paper框架里"69 分蘖期 BR P = 0.041、AUX P = 0.005；87 灌浆期 BR P = 0.023、drought P = 0.010、SA P = 0.022"
#   在代码中没有出处。本脚本按 2026-09-22 确认的方案首次实现。
# 方法（两种判据并列，均用 Fisher 精确test，因为每条pathway的gene数只有 1–12 个，卡方test不适用）：
#   判据 A（主）：上调gene（log2FC > 0）是否在pathway内富集 —— 2 × 2 表 = pathway内上调/未上调 × pathway外上调/未上调；
#   判据 B（辅）：显著gene（Wilcoxon FDR < 0.05）是否在pathway内富集。
#   每条pathway在每个 focal × 时期 内单独test，再按 focal × 时期 做 BH-FDR correction（20 条pathway为一族）。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：单步耗时 < 60 s，每次运行都重新generate。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/04_pathway_gene_expression/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_in <- file.path(project_root, "result", "04_pathway_gene_expression", "01_表_pathwaygene差异.xlsx")
out_dir <- file.path(project_root, "result", "04_pathway_gene_expression")
path_out <- file.path(out_dir, "02_表_pathway富集_TableS3S4.xlsx")

stopifnot(file.exists(path_in))

bioRtools::read_data(path_in, sheet = "stage_focal_level") -> focal_tab

stopifnot(nrow(focal_tab) > 0)

cat("[02] 输入：", nrow(focal_tab), "行；focal × 时期 组合数：",
    dplyr::n_distinct(paste(focal_tab$stage, focal_tab$focal)), "\n")

enrich_one <- function(d, criterion) {
  d %>%
    dplyr::mutate(is_hit = dplyr::case_when(
      criterion == "up" ~ log2FC > 0,
      criterion == "sig" ~ wilcoxon_fdr < 0.05,
      TRUE ~ NA
    )) -> dd

  # 每行一个 (pathway内/pathway外) × (命中/未命中) 的四格计数
  dd %>%
    dplyr::group_by(process) %>%
    dplyr::summarise(
      in_hit = sum(is_hit, na.rm = TRUE),
      in_total = dplyr::n(),
      .groups = "drop"
    ) %>%
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
  dplyr::group_modify(~ {
    dplyr::bind_rows(
      enrich_one(.x, "up"),
      enrich_one(.x, "sig")
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(stage, focal, criterion) %>%
  dplyr::mutate(fisher_fdr = p.adjust(fisher_p, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(stage, focal, criterion, fisher_p)

cat("\n[02] Table S3 —— 69 分蘖期（判据 A：上调gene富集）：\n")
res_tab %>%
  dplyr::filter(stage == "tillering", focal == 69, criterion == "up") %>%
  dplyr::select(process, in_hit, in_total, fisher_p, fisher_fdr, odds_ratio) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n[02] Table S4 —— 87 灌浆期（判据 A：上调gene富集）：\n")
res_tab %>%
  dplyr::filter(stage == "filling", focal == 87, criterion == "up") %>%
  dplyr::select(process, in_hit, in_total, fisher_p, fisher_fdr, odds_ratio) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n[02] 全部 focal × 时期 中 FDR < 0.05 的pathway：\n")
res_tab %>%
  dplyr::filter(fisher_fdr < 0.05) %>%
  dplyr::select(stage, focal, criterion, process, in_hit, in_total, fisher_p, fisher_fdr) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

bioRtools::write_data(
  list(
    TableS3_focal69_tillering = res_tab %>% dplyr::filter(stage == "tillering", focal == 69),
    TableS4_focal87_filling = res_tab %>% dplyr::filter(stage == "filling", focal == 87),
    all_focal_stage = res_tab
  ),
  file = path_out
)
cat("[02] 输出：", path_out, "\n")

cat("\n[02] sessionInfo()\n")
sessionInfo()
