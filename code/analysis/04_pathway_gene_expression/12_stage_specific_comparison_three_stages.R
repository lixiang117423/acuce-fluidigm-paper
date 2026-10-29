# ============================================================
# 04_pathway_gene_expression / 12_stage_specific_comparison_three_stages.R
# 作用：把"stage特异性"analysis扩展到三个时期——一个 limma 模型同时估计三期各自的 Mix vs Pure effect、
#       两两交互（花期−分蘖、灌浆−分蘖、灌浆−花期），再做重叠结构、方向一致性、effect量comparisonand特异gene归属
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_pathway_gene_expression/12_stage_specific_comparison_three_stages.R"
#
# 输入：<root>/result/03_overall_gene_expression/11_表_expression矩阵_three_stages.xlsx（sample_matrix、gene_annotation）
# 输出：<root>/result/04_pathway_gene_expression/29_表_stage_specific_comparison_three_stages.xlsx
#       <root>/result/04_pathway_gene_expression/30_FigureS15_stage_specific_comparison_three_stages.pdf / .png
# 其中 <root> = <repository_root>
#
# 模型：grp = stage × group（6 个水平），design `~ 0 + grp + focal`（focal 固定区组）；
#   contrasts：tillering / flowering / filling 各自为 (Mix − Pure)；
#              flowering_vs_tillering、filling_vs_tillering、filling_vs_flowering 为两个对比之差。
#   响应 = log2(expression + 1)，limma::eBayes(robust = TRUE)。
# and两时期version的差异：① stage由 2 变 3；② 交互对比由 1 个变 3 个；③ 图改为"跨三期的 logFC 平行坐标图 +
#   各期/各交互的 FDR 计数条形图"两个面板（两时期版是散点＋重叠条形＋特异gene条形）。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_mat <- file.path(project_root, "result", "03_overall_gene_expression", "11_表_expression矩阵_three_stages.xlsx")
out_dir <- file.path(project_root, "result", "04_pathway_gene_expression")
path_out <- file.path(out_dir, "29_表_stage_specific_comparison_three_stages.xlsx")

stopifnot(file.exists(path_mat))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

stage_levels <- c("tillering", "flowering", "filling")

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_mat, "sample_matrix") -> mat.sample
bioRtools::read_data(path_mat, sheet = "gene_annotation") -> df.gene

meta_cols <- c("stage", "sample", "focal", "neighbor")
gene_cols <- setdiff(colnames(mat.sample), meta_cols)
stopifnot(length(gene_cols) == 92)

mat.sample %>%
  dplyr::mutate(
    group = ifelse(focal == neighbor, "Pure", "Mix"),
    stage = factor(stage, levels = stage_levels),
    grp = factor(paste(stage, group, sep = "_"),
                 levels = as.vector(outer(stage_levels, c("Pure", "Mix"), paste, sep = "_"))),
    focal = factor(focal)
  ) -> d

cat("[12] sample：", nrow(d), "；六组sample量：\n")
print(table(d$grp, useNA = "ifany"))

expr <- t(log2(as.matrix(d[, gene_cols]) + 1))
design <- stats::model.matrix(~ 0 + grp + focal, data = d)
colnames(design) <- make.names(colnames(design))

cont <- limma::makeContrasts(
  tillering = grptillering_Mix - grptillering_Pure,
  flowering = grpflowering_Mix - grpflowering_Pure,
  filling = grpfilling_Mix - grpfilling_Pure,
  flowering_vs_tillering = (grpflowering_Mix - grpflowering_Pure) - (grptillering_Mix - grptillering_Pure),
  filling_vs_tillering = (grpfilling_Mix - grpfilling_Pure) - (grptillering_Mix - grptillering_Pure),
  filling_vs_flowering = (grpfilling_Mix - grpfilling_Pure) - (grpflowering_Mix - grpflowering_Pure),
  levels = design
)

fit <- limma::eBayes(limma::lmFit(expr, design), robust = TRUE)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont), robust = TRUE)

coefs <- colnames(cont)

extract_contrast <- function(coef_name) {
  limma::topTable(fit2, coef = coef_name, number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::transmute(gene = gene, logFC = logFC, t = t, p = P.Value, fdr = adj.P.Val, b = B,
                     contrast = coef_name)
}

long_tab <- purrr::map_dfr(coefs, extract_contrast)

wide_tab <- long_tab %>%
  tidyr::pivot_wider(names_from = contrast, values_from = c(logFC, t, p, fdr, b), names_sep = "_") %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene")

stage_tab <- wide_tab %>%
  dplyr::mutate(
    sig_tillering = fdr_tillering < 0.05,
    sig_flowering = fdr_flowering < 0.05,
    sig_filling = fdr_filling < 0.05,
    n_stage_sig = sig_tillering + sig_flowering + sig_filling,
    category = dplyr::case_when(
      n_stage_sig == 0 ~ "Neither",
      n_stage_sig == 3 ~ "All three stages",
      sig_tillering & sig_flowering ~ "Tillering + flowering",
      sig_tillering & sig_filling ~ "Tillering + filling",
      sig_flowering & sig_filling ~ "Flowering + filling",
      sig_tillering ~ "Tillering only",
      sig_flowering ~ "Flowering only",
      TRUE ~ "Filling only"
    ),
    sig_interaction = fdr_flowering_vs_tillering < 0.05 |
      fdr_filling_vs_tillering < 0.05 | fdr_filling_vs_flowering < 0.05,
    max_abs_logFC = pmax(abs(logFC_tillering), abs(logFC_flowering), abs(logFC_filling)),
    stage_of_max = dplyr::case_when(
      abs(logFC_tillering) == max_abs_logFC ~ "tillering",
      abs(logFC_flowering) == max_abs_logFC ~ "flowering",
      TRUE ~ "filling"
    ),
    range_logFC = max_abs_logFC - pmin(abs(logFC_tillering), abs(logFC_flowering), abs(logFC_filling)),
    step_specificity = range_logFC / (abs(logFC_tillering) + abs(logFC_flowering) + abs(logFC_filling) + 1e-9)
  )

overlap_tab <- stage_tab %>% dplyr::count(category) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  dplyr::arrange(dplyr::desc(n))

cat("\n[12] ① 各期 FDR < 0.05 gene数：\n")
print(as.data.frame(stage_tab %>%
                      dplyr::summarise(tillering = sum(sig_tillering),
                                       flowering = sum(sig_flowering),
                                       filling = sum(sig_filling),
                                       interaction_any = sum(sig_interaction))))

cat("\n[12] ② 重叠结构：\n")
print(as.data.frame(overlap_tab), row.names = FALSE)

eff_compare <- purrr::map_dfr(stage_levels, function(st) {
  other <- setdiff(stage_levels, st)
  d1 <- stage_tab[[paste0("logFC_", st)]]
  out <- tibble::tibble(stage = st, n_gene = length(d1),
                        median_abs_logFC = median(abs(d1)),
                        mean_abs_logFC = mean(abs(d1)),
                        max_abs_logFC = max(abs(d1)))
  for (o in other) {
    d2 <- stage_tab[[paste0("logFC_", o)]]
    out[[paste0("wilcoxon_p_vs_", o)]] <-
      suppressWarnings(stats::wilcox.test(abs(d1), abs(d2), paired = TRUE))$p.value
    out[[paste0("pearson_r_vs_", o)]] <- unname(stats::cor.test(d1, d2)$estimate)
  }
  out
})

cat("\n[12] ③ 三期effect量comparisonand logFC 相关性：\n")
print(as.data.frame(eff_compare), row.names = FALSE)

top_specific <- stage_tab %>%
  dplyr::filter(category != "Neither") %>%
  dplyr::arrange(dplyr::desc(step_specificity)) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::select(gene, process, logFC_tillering, logFC_flowering, logFC_filling,
                fdr_tillering, fdr_flowering, fdr_filling, stage_of_max, step_specificity, category)

if (nrow(top_specific) == 0) {
  top_specific <- stage_tab %>%
    dplyr::arrange(dplyr::desc(step_specificity)) %>%
    dplyr::slice_head(n = 15) %>%
    dplyr::select(gene, process, logFC_tillering, logFC_flowering, logFC_filling,
                  fdr_tillering, fdr_flowering, fdr_filling, stage_of_max, step_specificity, category)
}

# ---- 图：跨三期 logFC 平行坐标（effect最强的 20 个gene）+ 各对比 FDR 计数 ----
top20 <- stage_tab %>% dplyr::arrange(dplyr::desc(max_abs_logFC)) %>% dplyr::slice_head(n = 20) %>%
  dplyr::pull(gene)

p.trace <- long_tab %>%
  dplyr::filter(gene %in% top20, contrast %in% stage_levels) %>%
  dplyr::mutate(stage = factor(contrast, levels = stage_levels)) %>%
  ggplot2::ggplot(aes(x = stage, y = logFC, group = gene)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_line(alpha = 0.55, color = "grey35") +
  ggplot2::geom_point(size = 2, color = "#1e4e8c") +
  ggplot2::labs(x = NULL, y = "Moderated log2 fold change (Mix vs Pure)",
                title = "A  Effects of the 20 largest-magnitude genes across three stages") +
  bioRtools::theme_prism() +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 9, face = "bold"),
                 axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1))

count_tab <- tibble::tibble(
  contrast = coefs,
  n_fdr005 = purrr::map_int(coefs, ~ sum(long_tab$fdr[long_tab$contrast == .x] < 0.05)),
  n_p005 = purrr::map_int(coefs, ~ sum(long_tab$p[long_tab$contrast == .x] < 0.05))
)

p.bar <- count_tab %>%
  dplyr::mutate(contrast = factor(contrast, levels = rev(coefs))) %>%
  tidyr::pivot_longer(cols = c(n_p005, n_fdr005), names_to = "criterion", values_to = "n") %>%
  dplyr::mutate(criterion = dplyr::recode(criterion, n_p005 = "P < 0.05", n_fdr005 = "FDR < 0.05")) %>%
  ggplot2::ggplot(aes(x = n, y = contrast, fill = criterion)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
  ggplot2::geom_text(aes(label = n), position = ggplot2::position_dodge(width = 0.7),
                     hjust = -0.25, size = 3) +
  ggplot2::scale_fill_manual(values = c("P < 0.05" = "#e0a458", "FDR < 0.05" = "#a63446"), name = NULL) +
  ggplot2::labs(x = "Number of genes", y = NULL,
                title = "B  Candidate-responsive genes per contrast (92 genes)") +
  bioRtools::theme_prism() +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 9, face = "bold"),
                 legend.position = "bottom")

bioRtools::write_data(p.trace, file = file.path(out_dir, "30_FigureS15_stage_specific_comparison_three_stages.pdf"),
                      width = 8, height = 6, units = "in")
bioRtools::write_data(p.trace, file = file.path(out_dir, "30_FigureS15_stage_specific_comparison_three_stages.png"),
                      width = 8, height = 6, units = "in", dpi = 500)
bioRtools::write_data(p.bar, file = file.path(out_dir, "30b_FigureS15_各对比计数_three_stages.pdf"),
                      width = 8, height = 6, units = "in")
bioRtools::write_data(p.bar, file = file.path(out_dir, "30b_FigureS15_各对比计数_three_stages.png"),
                      width = 8, height = 6, units = "in", dpi = 500)
cat("\n[12] 输出：30_FigureS15_stage_specific_comparison_three_stages.pdf / .png\n")
cat("[12] 输出：30b_FigureS15_各对比计数_three_stages.pdf / .png\n")

bioRtools::write_data(
  list(stage_contrasts = stage_tab, overlap = overlap_tab, effect_compare = eff_compare,
       counts = count_tab, top_specific = top_specific, design_levels = tibble::tibble(grp = levels(d$grp))),
  file = path_out
)
cat("\n[12] 输出：", path_out, "\n")

cat("\n[12] sessionInfo()\n")
sessionInfo()
