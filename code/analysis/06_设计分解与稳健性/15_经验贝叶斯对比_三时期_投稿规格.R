# ============================================================
# 06_设计分解与稳健性 / 06_表达_经验贝叶斯差异检验.R
# 独立分析（续 06 主题，Q8）：用经验贝叶斯收缩（limma 类）替代 Wilcoxon 做基因水平差异检验
# 为什么换方法：主题 04 用 Wilcoxon 逐基因检验，每个基因的方差独立估计、且对 92 个基因做
#   多重校正，检出力很低（分蘖期只有少数基因 FDR < 0.05）。微阵列/低通量 qPCR 面板的常规做法是
#   用经验贝叶斯把各基因的方差向共同方差收缩（limma::eBayes），并用 focal 作为区组因子
#   控制品系效应，从而提高检出力并给出 moderated t 统计量。
# 方法：
#   1) 响应 = log2(expression + 1)；基因 × 样本矩阵按时期分别建模；
#   2) 设计矩阵 ~ 0 + group + focal（focal 作为固定区组因子，Pure/Mix 在每个 focal 内可比）；
#   3) limma::lmFit + contrasts.fit(显式 Mix − Pure 对比) + limma::eBayes(robust = TRUE)；
#      输出 logFC、moderated t、P、BH-FDR；
#   4) 与主题 04 的 Wilcoxon 结果逐基因比较（FDR < 0.05 的重合度、方向一致性）。
#
# 【2026-09-22 更正（重要）】第一版直接用 `topTable(fit, coef = "groupMix")`，但设计是 `~ 0 + group + focal`，
#   该写法下系数是**组均值**而非 Mix − Pure 差异，导致误报"分蘖期 53/92、灌浆期 59/92 个基因 FDR<0.05、
#   OsBLE1 FDR ≈ 1e-11"。现改为显式对比矩阵（groupMix − groupPure），focal 内分析同样改为显式对比。
#   更正后：分蘖期 14 个基因 P<0.05（无 FDR<0.05，最强 OsGNS1 logFC = 1.15、P = 0.0024），灌浆期仅 1 个。
# 注意：本面板每个 gene × sample 只有 1 个技术值，故不设 block/duplicateCorrelation；
#   focal 固定效应已吸收品系差异。
#
# 运行：/opt/homebrew/bin/Rscript "<root>/code/06_设计分解与稳健性/06_表达_经验贝叶斯差异检验.R"
# 输入：<root>/result/03_基因表达总体/11_表_表达矩阵_三时期.xlsx（sample_matrix、gene_annotation，254 样本）
#       <root>/result/04_通路基因表达/23_表_limma统一口径统计_三时期.xlsx（gene_stage_level，官方三时期 limma 结果）
# 输出：<root>/result/08_投稿规格图/06_表_表达_经验贝叶斯差异检验_三时期.xlsx
#       <root>/result/08_投稿规格图/06b_FigureS13_经验贝叶斯火山图_三时期_投稿规格.pdf / .png
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录见 code/06_设计分解与稳健性/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_mat <- file.path(project_root, "result", "03_基因表达总体", "11_表_表达矩阵_三时期.xlsx")
path_limma <- file.path(project_root, "result", "04_通路基因表达", "23_表_limma统一口径统计_三时期.xlsx")
out_dir <- file.path(project_root, "result", "06_设计分解与稳健性")
out_pub_dir <- file.path(project_root, "result", "08_投稿规格图")  # 投稿规格输出目录
dir.create(out_pub_dir, showWarnings = FALSE, recursive = TRUE)
path_out <- file.path(out_pub_dir, "06_表_表达_经验贝叶斯差异检验_三时期.xlsx")
stage_levels <- c("tillering", "flowering", "filling")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stopifnot(file.exists(path_mat), file.exists(path_limma))

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_mat, "sample_matrix") -> mat.sample
bioRtools::read_data(path_mat, sheet = "gene_annotation") -> df.gene
stopifnot(nrow(mat.sample) == 254)

meta_cols <- c("stage", "sample", "focal", "neighbor")
gene_cols <- setdiff(colnames(mat.sample), meta_cols)
stopifnot(length(gene_cols) == 92)

mat.sample %>%
  dplyr::mutate(group = factor(dplyr::if_else(focal == neighbor, "Pure", "Mix"),
                               levels = c("Pure", "Mix"))) -> mat.sample

# ---- 三时期 limma（经验贝叶斯收缩 + focal 区组）：直接读官方结果 ----
# 说明：同一件事不做两套实现（项目规则）——三时期 limma 由
#   code/04_通路基因表达/10_limma统一口径统计_三时期.R 生成（输入同为 11_表 254 样本矩阵），
#   本脚本只读取，保证与 Table S8 / Figure 4A/4B 完全一致。
bioRtools::read_data(path_limma, sheet = "gene_stage_level") -> limma_tab
limma_tab <- limma_tab %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels)) %>%
  dplyr::arrange(stage, limma_p)

# ---- 同口径的 Wilcoxon（每时期 × 每基因，Mix vs Pure；BH 校正 92 个基因）----
wilcox_one <- function(st) {
  d <- mat.sample %>% dplyr::filter(stage == st)
  purrr::map_dfr(gene_cols, function(g) {
    x <- d[[g]][d$group == "Mix"]
    y <- d[[g]][d$group == "Pure"]
    wt <- suppressWarnings(stats::wilcox.test(x, y))
    tibble::tibble(stage = st, gene = g,
                   log2FC = log2(mean(x) + 1) - log2(mean(y) + 1),
                   wilcoxon_p = wt$p.value)
  })
}
wilcox_tab <- purrr::map_dfr(stage_levels, wilcox_one) %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(wilcoxon_fdr = p.adjust(wilcoxon_p, method = "BH")) %>%
  dplyr::ungroup()
cat("[15] 三时期 Wilcoxon（对照口径）生成：", nrow(wilcox_tab), "行\n")

cat("[06] limma（focal 为区组、经验贝叶斯收缩）结果：\n")
print(as.data.frame(limma_tab %>%
                      dplyr::group_by(stage) %>%
                      dplyr::summarise(n_gene = dplyr::n(),
                                       n_fdr_005 = sum(limma_fdr < 0.05),
                                       n_p_005 = sum(limma_p < 0.05),
                                       min_fdr = min(limma_fdr), .groups = "drop")))

cat("\n[06] 分蘖期 FDR < 0.05 的基因（前 12 个）：\n")
print(as.data.frame(limma_tab %>% dplyr::filter(stage == "tillering", limma_fdr < 0.05) %>%
                      dplyr::select(gene, process, limma_logFC, limma_t, limma_fdr) %>%
                      dplyr::slice_head(n = 12)), row.names = FALSE)

cat("\n[06] 灌浆期 FDR < 0.05 的基因（前 12 个）：\n")
print(as.data.frame(limma_tab %>% dplyr::filter(stage == "filling", limma_fdr < 0.05) %>%
                      dplyr::select(gene, process, limma_logFC, limma_t, limma_fdr) %>%
                      dplyr::slice_head(n = 12)), row.names = FALSE)

# ---- 与 Wilcoxon 方法对比 ----
compare_tab <- limma_tab %>%
  dplyr::left_join(wilcox_tab %>% dplyr::select(stage, gene, wilcoxon_p, wilcoxon_fdr, log2FC),
                   by = c("stage", "gene")) %>%
  dplyr::mutate(
    method_agree_direction = sign(limma_logFC) == sign(log2FC),
    sig_limma = limma_fdr < 0.05,
    sig_wilcox = wilcoxon_fdr < 0.05
  )

cat("\n[06] 两种方法的比较（按时期）：\n")
print(as.data.frame(compare_tab %>%
                      dplyr::group_by(stage) %>%
                      dplyr::summarise(
                        n_gene = dplyr::n(),
                        n_sig_limma = sum(sig_limma, na.rm = TRUE),
                        n_sig_wilcoxon = sum(sig_wilcox, na.rm = TRUE),
                        n_sig_both = sum(sig_limma & sig_wilcox, na.rm = TRUE),
                        direction_agree_pct = 100 * mean(method_agree_direction, na.rm = TRUE),
                        min_p_limma = min(limma_p), min_p_wilcoxon = min(wilcoxon_p, na.rm = TRUE),
                        .groups = "drop")))

# 关键基因核对（69 分蘖期 / 87 灌浆期）
key_check <- compare_tab %>%
  dplyr::filter((stage == "tillering" & gene %in% c("OsBLE1", "OsIAA4", "OsGH3.2", "OsGBSSII", "OsGNS1")) |
                  (stage == "filling" & gene %in% c("OsBLE1", "OsWRKY45", "OsWRKY7", "OsGNS1"))) %>%
  dplyr::select(stage, gene, process, limma_logFC, limma_p, limma_fdr,
                wilcoxon_p, wilcoxon_fdr, log2FC)

cat("\n[06] 关键基因在两种方法下的结果：\n")
print(as.data.frame(key_check), row.names = FALSE)

# ---- 稳健性检查：显著基因的效应是否在每个 focal 内方向一致 ----
# 目的：limma 的 group 效应是在 focal 区组内估计的，但表中没有芯片/批次变量；
#   若 Pure 与 Mix 落在不同批次，group 效应可能混入技术批次。这里的检查利用
#   "4 个 focal 分别独立可比"这一结构：效应若在 4 个 focal 内同向，则批次混杂的可能性较低。
limma_by_focal <- function(st, f) {
  d <- mat.sample %>% dplyr::filter(stage == st, focal == f)
  if (dplyr::n_distinct(d$group) < 2) return(NULL)
  expr <- t(log2(as.matrix(d[, gene_cols]) + 1))
  design <- stats::model.matrix(~ 0 + group, data = d)
  colnames(design) <- make.names(colnames(design))
  con <- matrix(0, nrow = ncol(design), ncol = 1,
                dimnames = list(colnames(design), "mix_vs_pure"))
  con["groupMix", "mix_vs_pure"] <- 1
  con["groupPure", "mix_vs_pure"] <- -1
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr, design), con), robust = TRUE)
  limma::topTable(fit, coef = "mix_vs_pure", number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::transmute(stage = st, focal = f, gene = gene, logFC_focal = logFC, p_focal = P.Value)
}

focal_lfc <- tidyr::expand_grid(stage = stage_levels,
                                focal = unique(mat.sample$focal)) %>%
  purrr::pmap_dfr(function(stage, focal) limma_by_focal(stage, focal))

consistency_tab <- focal_lfc %>%
  dplyr::group_by(stage, gene) %>%
  dplyr::summarise(
    n_focal = dplyr::n(),
    n_same_sign = sum(sign(logFC_focal) == sign(logFC_focal[1])),
    min_logFC = min(logFC_focal), max_logFC = max(logFC_focal),
    .groups = "drop"
  ) %>%
  dplyr::left_join(limma_tab %>% dplyr::select(stage, gene, limma_logFC, limma_fdr), by = c("stage", "gene")) %>%
  dplyr::mutate(all_focal_consistent = n_same_sign == n_focal)

cat("\n[06] 稳健性检查（显著基因在 4 个 focal 内的方向一致性）：\n")
print(as.data.frame(consistency_tab %>%
                      dplyr::filter(limma_fdr < 0.05) %>%
                      dplyr::group_by(stage) %>%
                      dplyr::summarise(n_sig = dplyr::n(),
                                       n_all_consistent = sum(all_focal_consistent),
                                       pct_all_consistent = 100 * mean(all_focal_consistent),
                                       .groups = "drop")))

bioRtools::write_data(
  list(limma_results = limma_tab, method_comparison = compare_tab, key_genes = key_check,
       focal_consistency = consistency_tab),
  file = path_out
)
cat("[06] 输出：", path_out, "\n")

# ---- 火山图 ----
p.volcano <- limma_tab %>%
  dplyr::mutate(
    stage = factor(stage, levels = c("tillering", "flowering", "filling")),
    sig = dplyr::case_when(limma_fdr < 0.05 ~ "FDR < 0.05",
                           limma_p < 0.05 ~ "P < 0.05",
                           TRUE ~ "n.s."),
    sig = factor(sig, levels = c("FDR < 0.05", "P < 0.05", "n.s."))
  ) %>%
  ggplot2::ggplot(aes(x = limma_logFC, y = -log10(limma_p))) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  ggplot2::geom_point(aes(color = sig), size = 1.8, alpha = 0.85) +
  ggrepel::geom_text_repel(
    data = . %>% dplyr::filter(limma_fdr < 0.05 | limma_p < 0.01),
    aes(label = gene), size = 2.6, fontface = "italic", max.overlaps = 25
  ) +
  ggplot2::facet_wrap(~stage) +
  ggplot2::scale_color_manual(values = c("FDR < 0.05" = "#a63446", "P < 0.05" = "#e0a458", "n.s." = "grey60")) +
  ggplot2::labs(
    x = "log2 fold change (Mix vs Pure, focal-adjusted)",
    y = "-log10 P (moderated t)",
    color = NULL
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 8),
                 legend.position = "bottom",
                 legend.key.size = ggplot2::unit(0.35, "cm"),
                 legend.text = ggplot2::element_text(size = 7),
                 axis.title = ggplot2::element_text(size = 7.5)) -> p.volcano

bioRtools::write_data(p.volcano,
                      file = file.path(out_pub_dir, "06b_FigureS13_经验贝叶斯火山图_三时期_投稿规格.pdf"),
                      width = 17.0, height = 8.5, units = "cm")
bioRtools::write_data(p.volcano,
                      file = file.path(out_pub_dir, "06b_FigureS13_经验贝叶斯火山图_三时期_投稿规格.png"),
                      width = 17.0, height = 8.5, units = "cm", dpi = 600)
cat("[15] 输出：06b_FigureS13_经验贝叶斯火山图_三时期_投稿规格.pdf / .png\n")

cat("\n[06] sessionInfo()\n")
sessionInfo()
