# ============================================================
# 04_通路基因表达 / 20_通路富集检验_TableS4_三时期.R
# 作用：通路水平富集检验（Fisher 精确检验）**统一成三时期口径** —— 供 Table S4 使用
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_通路基因表达/20_通路富集检验_TableS4_三时期.R"
#
# 与原两时期版（02_通路富集检验_TableS3S4.R）的差别：
#   ① 输入改用**三时期 limma（统一口径）**的 focal × 时期 结果（与 Figure 3A/4A/4B、Table S8 同源）；
#   ② 判据 A（主）＝通路内上调基因（focal_logFC > 0）是否富集；
#      判据 B（辅）＝通路内显著基因（focal_fdr < 0.05）是否富集；
#   ③ 覆盖 **3 个时期 × 4 个 focal × 2 个判据 × 20 条通路**，FDR 按 (时期 × focal × 判据) 校正；
#   ④ 输出写到 result/08_投稿规格图（提交用表），原 result/04_通路基因表达 只读不改。
#
# 输入：<root>/result/04_通路基因表达/23_表_limma统一口径统计_三时期.xlsx（sheet gene_stage_focal_level）
# 输出：<root>/result/08_投稿规格图/32_表_通路富集_Fisher_三时期.xlsx
#       （sheet：TableS4_focal87_filling、TableS4_focal69_tillering、all_stage_focal、summary）
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_in <- file.path(project_root, "result", "04_通路基因表达", "23_表_limma统一口径统计_三时期.xlsx")
out_pub_dir <- file.path(project_root, "result", "08_投稿规格图")
dir.create(out_pub_dir, showWarnings = FALSE, recursive = TRUE)
path_out <- file.path(out_pub_dir, "32_表_通路富集_Fisher_三时期.xlsx")

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

cat("\n[20] FDR < 0.05 的通路（全部 focal × 时期）：\n")
res_tab %>%
  dplyr::filter(fisher_fdr < 0.05) %>%
  dplyr::select(stage, focal, criterion, process, in_hit, in_total, fisher_p, fisher_fdr) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n[20] 未校正 P < 0.05 的通路数：", sum(res_tab$fisher_p < 0.05),
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
