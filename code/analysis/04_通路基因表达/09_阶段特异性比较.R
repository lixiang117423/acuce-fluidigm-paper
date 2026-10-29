# ============================================================
# 04_通路基因表达 / 09_阶段特异性比较.R
# 作用：回答标题级问题——"分蘖期是转录响应的主窗口"是否成立，还是"两阶段各有一套部分重叠的响应"
#       做法：一个 limma 模型同时估计两个阶段的 Mix vs Pure 效应及其差值（交互），再做重叠度、
#       方向一致性、效应量比较与阶段特异基因/pathway 归属分析
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_通路基因表达/09_阶段特异性比较.R"
#
# 输入：<root>/result/03_基因表达总体/01_表_表达矩阵.xlsx（sample_matrix、gene_annotation）
# 输出：<root>/result/04_通路基因表达/21_表_阶段特异性比较.xlsx
#       <root>/result/04_通路基因表达/22_FigureS15_阶段特异性比较.pdf / .png
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 模型：grp = stage × group（4 个水平），设计 `~ 0 + grp + focal`（focal 固定区组），
#   contrasts：tillering = grp.tillering_Mix − grp.tillering_Pure；
#              filling   = grp.filling_Mix   − grp.filling_Pure；
#              interaction = tillering − filling（检验两阶段效应是否不同）。
#   响应 = log2(expression + 1)，limma::eBayes(robust = TRUE)。
# 判读标准（本脚本要回答的问题）：
#   ① 显著基因的重叠结构（两阶段共有 / 各自特异）；② 共有基因的方向一致性；
#   ③ 两阶段效应量（|logFC|）是否系统性不同；④ 交互检验显著的基因数。
#   若"分蘖期特异"与"灌浆期特异"数量相当、且共有基因方向一致 → 应写"阶段特异的两套响应"，
#   而不是"分蘖期是唯一主窗口"。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录见 code/04_通路基因表达/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_mat <- file.path(project_root, "result", "03_基因表达总体", "01_表_表达矩阵.xlsx")
out_dir <- file.path(project_root, "result", "04_通路基因表达")
path_out <- file.path(out_dir, "21_表_阶段特异性比较.xlsx")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stopifnot(file.exists(path_mat))

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
    grp = factor(paste(stage, group, sep = "_"),
                 levels = c("tillering_Pure", "tillering_Mix", "filling_Pure", "filling_Mix")),
    focal = factor(focal)
  ) -> d

cat("[09] 样本：", nrow(d), "；四组样本量：\n")
print(table(d$grp, useNA = "ifany"))

expr <- t(log2(as.matrix(d[, gene_cols]) + 1))
design <- stats::model.matrix(~ 0 + grp + focal, data = d)
colnames(design) <- make.names(colnames(design))

need <- c("grptillering_Pure", "grptillering_Mix", "grpfilling_Pure", "grpfilling_Mix")
stopifnot(all(need %in% colnames(design)))

cont <- limma::makeContrasts(
  tillering = grptillering_Mix - grptillering_Pure,
  filling = grpfilling_Mix - grpfilling_Pure,
  interaction = (grptillering_Mix - grptillering_Pure) - (grpfilling_Mix - grpfilling_Pure),
  levels = design
)

fit <- limma::eBayes(limma::lmFit(expr, design), robust = TRUE)
fit2 <- limma::contrasts.fit(fit, cont)
fit2 <- limma::eBayes(fit2, robust = TRUE)

extract_contrast <- function(coef_name) {
  limma::topTable(fit2, coef = coef_name, number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::transmute(gene = gene,
                     logFC = logFC, t = t, p = P.Value, fdr = adj.P.Val, b = B,
                     contrast = coef_name)
}

stage_tab <- dplyr::bind_rows(purrr::map(c("tillering", "filling", "interaction"), extract_contrast)) %>%
  tidyr::pivot_wider(names_from = contrast, values_from = c(logFC, t, p, fdr, b),
                     names_sep = "_") %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene") %>%
  dplyr::mutate(
    sig_tillering = fdr_tillering < 0.05,
    sig_filling = fdr_filling < 0.05,
    sig_interaction = fdr_interaction < 0.05,
    category = dplyr::case_when(
      sig_tillering & sig_filling ~ "Both stages",
      sig_tillering ~ "Tillering only",
      sig_filling ~ "Filling only",
      TRUE ~ "Neither"
    ),
    category = factor(category, levels = c("Both stages", "Tillering only", "Filling only", "Neither")),
    direction_concordant = dplyr::if_else(sign(logFC_tillering) == sign(logFC_filling), TRUE, FALSE),
    abs_diff = abs(logFC_tillering - logFC_filling),
    stage_specificity = abs_diff / (abs(logFC_tillering) + abs(logFC_filling) + 1e-9),
    stronger_stage = dplyr::if_else(abs(logFC_tillering) >= abs(logFC_filling), "tillering", "filling")
  )

# ---- ① 重叠结构 ----
overlap_tab <- stage_tab %>%
  dplyr::count(category, .drop = FALSE) %>%
  dplyr::mutate(pct = 100 * n / sum(n))

# 共有基因可能为 0（本次即是），必须显式处理空集，否则 summarise 返回 0 行、下游会崩
both_tab <- stage_tab %>% dplyr::filter(category == "Both stages")

both_dir <- if (nrow(both_tab) == 0) {
  tibble::tibble(
    n_both = 0L, n_same_direction = 0L, pct_same_direction = NA_real_,
    n_sig_interaction = 0L,
    median_abs_logFC_tillering = NA_real_, median_abs_logFC_filling = NA_real_
  )
} else {
  both_tab %>%
    dplyr::summarise(
      n_both = dplyr::n(),
      n_same_direction = sum(direction_concordant),
      pct_same_direction = 100 * mean(direction_concordant),
      n_sig_interaction = sum(sig_interaction),
      median_abs_logFC_tillering = median(abs(logFC_tillering)),
      median_abs_logFC_filling = median(abs(logFC_filling))
    )
}

cat("\n[09] ① 显著基因的重叠结构：\n")
print(as.data.frame(overlap_tab), row.names = FALSE)
cat("\n[09] 共有基因的方向一致性与交互检验：\n")
print(as.data.frame(both_dir), row.names = FALSE)

# ---- ② 效应量比较（全基因与显著基因）----
eff_compare <- purrr::map_dfr(c("all genes", "Both stages"), function(scope) {
  dd <- if (scope == "all genes") stage_tab else dplyr::filter(stage_tab, category == "Both stages")
  if (nrow(dd) < 3) {
    return(tibble::tibble(
      scope = scope, n_gene = nrow(dd),
      median_abs_tillering = if (nrow(dd) == 0) NA_real_ else median(abs(dd$logFC_tillering)),
      median_abs_filling = if (nrow(dd) == 0) NA_real_ else median(abs(dd$logFC_filling)),
      wilcoxon_p_abs = NA_real_, pearson_r_logFC = NA_real_, pearson_p = NA_real_,
      n_interaction_fdr005 = sum(dd$sig_interaction)
    ))
  }
  wt <- suppressWarnings(stats::wilcox.test(abs(dd$logFC_tillering), abs(dd$logFC_filling), paired = TRUE))
  ct <- stats::cor.test(dd$logFC_tillering, dd$logFC_filling, method = "pearson")
  tibble::tibble(
    scope = scope, n_gene = nrow(dd),
    median_abs_tillering = median(abs(dd$logFC_tillering)),
    median_abs_filling = median(abs(dd$logFC_filling)),
    wilcoxon_p_abs = wt$p.value,
    pearson_r_logFC = unname(ct$estimate), pearson_p = ct$p.value,
    n_interaction_fdr005 = sum(dd$sig_interaction)
  )
})

cat("\n[09] ② 两阶段效应量比较与 logFC 相关性：\n")
print(as.data.frame(eff_compare), row.names = FALSE)

# ---- ③ 阶段特异最强的基因 ----
# 若没有任何显著基因（本次即为 0 个共有基因、仅 8 个分蘖期基因），退化为"效应量最大的 15 个基因"，
# 并在表里保留 category 列以便区分
top_specific <- stage_tab %>%
  dplyr::filter(category != "Neither") %>%
  dplyr::arrange(dplyr::desc(stage_specificity)) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::select(gene, process, logFC_tillering, logFC_filling, fdr_tillering, fdr_filling,
                stage_specificity, stronger_stage, category)

if (nrow(top_specific) == 0) {
  top_specific <- stage_tab %>%
    dplyr::mutate(max_abs = pmax(abs(logFC_tillering), abs(logFC_filling))) %>%
    dplyr::arrange(dplyr::desc(max_abs)) %>%
    dplyr::slice_head(n = 15) %>%
    dplyr::select(gene, process, logFC_tillering, logFC_filling, fdr_tillering, fdr_filling,
                  stage_specificity, stronger_stage, category)
}

cat("\n[09] ③ 阶段特异性最强的 15 个基因：\n")
print(as.data.frame(top_specific), row.names = FALSE)

# ---- ④ 通路层面 ----
pathway_tab <- stage_tab %>%
  dplyr::group_by(process) %>%
  dplyr::summarise(
    n_gene = dplyr::n(),
    n_sig_tillering = sum(sig_tillering),
    n_sig_filling = sum(sig_filling),
    n_both = sum(category == "Both stages"),
    n_tillering_only = sum(category == "Tillering only"),
    n_filling_only = sum(category == "Filling only"),
    mean_abs_tillering = mean(abs(logFC_tillering)),
    mean_abs_filling = mean(abs(logFC_filling)),
    mean_stage_specificity = mean(stage_specificity),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_sig_filling + n_sig_tillering))

cat("\n[09] ④ 通路层面（按显著基因总数排序）：\n")
print(as.data.frame(pathway_tab), row.names = FALSE)

# camera：两阶段分别做通路检验（看方向是否不同）
camera_stage <- purrr::map_dfr(c("tillering", "filling"), function(cn) {
  limma::camera(expr, index = split(df.gene$gene[df.gene$gene %in% gene_cols],
                                    df.gene$process[df.gene$gene %in% gene_cols]),
                design = design, contrast = cont[, cn], inter.gene.cor = 0.01) %>%
    tibble::rownames_to_column(var = "process") %>%
    dplyr::mutate(stage = cn, .before = 1) %>%
    dplyr::select(stage, process, n_gene = NGenes, direction = Direction,
                  camera_p = PValue, camera_fdr = FDR)
})

cat("\n[09] camera 对比（两阶段的通路方向）：\n")
print(as.data.frame(camera_stage %>%
                      dplyr::select(stage, process, n_gene, direction, camera_p) %>%
                      tidyr::pivot_wider(names_from = stage, values_from = c(direction, camera_p)) %>%
                      dplyr::arrange(camera_p_tillering)), row.names = FALSE)

bioRtools::write_data(
  list(
    gene_stage_effects = stage_tab,
    overlap_structure = overlap_tab,
    both_stage_summary = both_dir,
    effect_size_comparison = eff_compare,
    top_stage_specific = top_specific,
    pathway_summary = pathway_tab,
    camera_by_stage = camera_stage
  ),
  file = path_out
)
cat("[09] 输出：", path_out, "\n")

# ---- Figure S15 ----
cat_colors <- c("Both stages" = "#1e4e8c", "Tillering only" = "#e0a458",
                "Filling only" = "#a63446", "Neither" = "grey75")

p.scatter <- stage_tab %>%
  ggplot2::ggplot(aes(x = logFC_tillering, y = logFC_filling, color = category)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey45") +
  ggplot2::geom_point(size = 2, alpha = 0.9) +
  ggplot2::scale_color_manual(values = cat_colors, name = NULL) +
  ggplot2::labs(
    x = "log2 fold change, tillering (Mix vs Pure)",
    y = "log2 fold change, filling (Mix vs Pure)",
    subtitle = paste0("Pearson r = ", formatC(eff_compare$pearson_r_logFC[1], format = "f", digits = 2),
                      " (all 92 genes); interaction FDR < 0.05: ",
                      sum(stage_tab$sig_interaction), " genes")
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    legend.position = "bottom",
    plot.subtitle = ggplot2::element_text(size = 8)
  ) -> p.scatter

p.bar <- overlap_tab %>%
  ggplot2::ggplot(aes(x = category, y = n, fill = category)) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::geom_text(aes(label = paste0(n, " (", formatC(pct, format = "f", digits = 1), "%)")),
                     vjust = -0.4, size = 3.2) +
  ggplot2::scale_fill_manual(values = cat_colors, guide = "none") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
  ggplot2::labs(x = NULL, y = "Number of genes",
                subtitle = if (both_dir$n_both == 0) {
                  "No gene is significant in both stages"
                } else {
                  paste0("Among the ", both_dir$n_both, " shared genes: ",
                         formatC(both_dir$pct_same_direction, format = "f", digits = 0),
                         "% same direction")
                }) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 20, hjust = 1, vjust = 1),
    plot.subtitle = ggplot2::element_text(size = 8)
  )

p.spec <- top_specific %>%
  dplyr::mutate(label = paste0(gene, " (", stronger_stage, ")")) %>%
  ggplot2::ggplot(aes(x = stage_specificity, y = forcats::fct_reorder(label, stage_specificity),
                      fill = stronger_stage)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::scale_fill_manual(values = c(tillering = "#e0a458", filling = "#a63446"), name = "Stronger stage") +
  ggplot2::labs(x = "Stage-specificity index", y = NULL,
                subtitle = "|logFC_t - logFC_f| / (|logFC_t| + |logFC_f|)") +
  bioRtools::theme_prism() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 7, face = "italic"),
    legend.position = "bottom",
    plot.subtitle = ggplot2::element_text(size = 8)
  )

p.s15 <- patchwork::wrap_plots(p.scatter, patchwork::wrap_plots(p.bar, p.spec, ncol = 1), ncol = 2, widths = c(1.1, 1))

bioRtools::write_data(p.s15, file = file.path(out_dir, "22_FigureS15_阶段特异性比较.pdf"),
                      width = 13, height = 8, units = "in")
bioRtools::write_data(p.s15, file = file.path(out_dir, "22_FigureS15_阶段特异性比较.png"),
                      width = 13, height = 8, units = "in", dpi = 500)
cat("[09] 输出：22_FigureS15_阶段特异性比较.pdf / .png\n")

cat("\n[09] sessionInfo()\n")
sessionInfo()
