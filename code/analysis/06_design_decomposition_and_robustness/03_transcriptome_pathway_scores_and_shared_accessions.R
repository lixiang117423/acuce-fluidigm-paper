# ============================================================
# 06_design_decomposition_and_robustness / 03_transcriptome_pathway_scores_and_shared_accessions.R
# 独立analysis（不属于旧版流程）：用"pathway整体分数"替代"逐gene卡方/Fisher test"，并把shared_accessions造成的非独立性纳入模型
#
# 【为什么要换方法】
#   1) 旧版（and我在主题 04 的实现）都按"每条pathway 1–12 个gene、显著/上调gene是否富集"做卡方或 Fisher test。
#      这个思路在pathway只有 2–3 个gene时几乎没有检出力（主题 04 的 Table S3/S4 因此全为阴性），
#      而且"显著gene富集"本身依赖gene水平的显著性判定，属于双重筛选。
#   2) 更合适的做法是把 1 条pathway当作 1 个变量：pathway内gene先按sample Z-score，再取均值得到"pathway_scores"，
#      直接test分数在 Pure / Mix 间的差异。这样使用了全部gene信息，且不依赖单个gene的显著性。
#   3) 12 个混种组合并非独立sample——15-69 and 69-15 共享同一对品系，69-15 and 69-87 shared_accessions 69。
#      因此在模型中加入品系随机截距 (1 | focal)，避免把shared_accessions当作独立重复而高估显著性。
#
# 运行：/opt/homebrew/bin/Rscript "<root>/code/06_design_decomposition_and_robustness/03_transcriptome_pathway_scores_and_shared_accessions.R"
# 输入：<root>/result/03_overall_gene_expression/01_表_expression矩阵.xlsx（sample_matrix、gene_annotation）
# 输出：<root>/result/06_design_decomposition_and_robustness/03_表_pathway_scorestest.xlsx
#       <root>/result/06_design_decomposition_and_robustness/03_FigureS10_pathway_scores森林图.pdf / .png
# 其中 <root> = <repository_root>
#
# 方法细节：① gene Z-score 在**每个时期内部**计算（用期内sample），消除时期间的整体尺度差异；
#   ② pathway_scores = pathway内gene z 的算术均值；③ 主模型 `lmer(score ~ group + (1 | focal))`（Satterthwaite df）；
#   ④ 另做 999 次标签置换（在每个时期内打乱 group 标签）得到经验 P 值，作为不依赖分布假设的对照；
#   ⑤ FDR 按时期在 20 条pathway上做 BH correction。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/06_design_decomposition_and_robustness/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_mat <- file.path(project_root, "result", "03_overall_gene_expression", "01_表_expression矩阵.xlsx")
out_dir <- file.path(project_root, "result", "06_design_decomposition_and_robustness")
path_out <- file.path(out_dir, "03_表_pathway_scorestest.xlsx")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stopifnot(file.exists(path_mat))

# 矩阵 sheet 必须保留gene名原样（bioRtools::read_data 会把 OsD14-1 改成 OsD14.1）
read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_mat, "sample_matrix") -> mat.sample
bioRtools::read_data(path_mat, sheet = "gene_annotation") -> df.gene

meta_cols <- c("stage", "sample", "focal", "neighbor")
gene_cols <- setdiff(colnames(mat.sample), meta_cols)
stopifnot(length(gene_cols) == 92, nrow(mat.sample) == 168)

df.gene.sub <- df.gene %>% dplyr::filter(gene %in% gene_cols)
stopifnot(nrow(df.gene.sub) == 92)

cat("[03] sample级矩阵：", nrow(mat.sample), " 行 ×", length(gene_cols), " gene；pathway ",
    dplyr::n_distinct(df.gene.sub$process), " 条\n")

# ---- 期内gene Z-score ----
mat.sample %>%
  # 显式设定 factor 水平：否则字母序会让 Mix 成为参照水平，系数名变成 groupPure，
  # 后续按 co["groupMix", ...] 取值会下标越界
  dplyr::mutate(group = factor(dplyr::if_else(focal == neighbor, "Pure", "Mix"),
                               levels = c("Pure", "Mix"))) %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(gene_cols), ~ as.vector(scale(.x)))) %>%
  dplyr::ungroup() -> df.z

# ---- pathway_scores ----
pathway_score <- function(st) {
  d <- df.z %>% dplyr::filter(stage == st)
  purrr::map_dfr(unique(df.gene.sub$process), function(pr) {
    genes <- df.gene.sub$gene[df.gene.sub$process == pr]
    genes <- intersect(genes, gene_cols)
    score <- rowMeans(d[, genes, drop = FALSE])
    tibble::tibble(
      stage = st,
      process = pr,
      n_gene = length(genes),
      sample = d$sample,
      focal = d$focal,
      neighbor = d$neighbor,
      group = d$group,
      score = score
    )
  })
}

score_tab <- purrr::map_dfr(c("tillering", "filling"), pathway_score)

cat("[03] pathway_scores表：", nrow(score_tab), " 行\n")

# ---- test：LMM（品系随机截距）+ 置换test ----
test_one <- function(d, pr) {
  m <- lmerTest::lmer(score ~ group + (1 | focal), data = d)
  co <- summary(m)$coefficients

  obs_diff <- mean(d$score[d$group == "Mix"]) - mean(d$score[d$group == "Pure"])
  n_gene_val <- length(intersect(df.gene.sub$gene[df.gene.sub$process == pr], gene_cols))

  perm <- replicate(999, {
    g <- sample(d$group)
    mean(d$score[g == "Mix"]) - mean(d$score[g == "Pure"])
  })

  tibble::tibble(
    process = pr,
    n_gene = n_gene_val,
    n_pure = sum(d$group == "Pure"),
    n_mix = sum(d$group == "Mix"),
    mean_pure = mean(d$score[d$group == "Pure"]),
    mean_mix = mean(d$score[d$group == "Mix"]),
    diff_mix_pure = obs_diff,
    lmm_estimate = co["groupMix", "Estimate"],
    lmm_se = co["groupMix", "Std. Error"],
    lmm_df = co["groupMix", "df"],
    lmm_t = co["groupMix", "t value"],
    lmm_p = co["groupMix", "Pr(>|t|)"],
    perm_p = (sum(abs(perm) >= abs(obs_diff)) + 1) / (length(perm) + 1)
  )
}

score_test <- score_tab %>%
  dplyr::group_by(stage) %>%
  dplyr::group_modify(~ {
    purrr::map_dfr(unique(.x$process), function(pr) {
      test_one(.x %>% dplyr::filter(process == pr), pr)
    })
  }) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(
    lmm_fdr = p.adjust(lmm_p, method = "BH"),
    perm_fdr = p.adjust(perm_p, method = "BH")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(stage, lmm_p)

cat("\n[03] pathway_scorestest（按 LMM P 排序，前 15 条）：\n")
print(as.data.frame(score_test %>%
                      dplyr::select(stage, process, n_gene, diff_mix_pure, lmm_estimate,
                                    lmm_p, lmm_fdr, perm_p, perm_fdr) %>%
                      dplyr::slice_head(n = 15)), row.names = FALSE)

cat("\n[03] LMM FDR < 0.05 的pathway数：", sum(score_test$lmm_fdr < 0.05),
    "；置换test FDR < 0.05 的pathway数：", sum(score_test$perm_fdr < 0.05), "\n")

cat("\n[03] 关键pathway（Auxin / Brassinosteroid / High sugar / Defense / Strigolactone / Jasmonic acid）：\n")
print(as.data.frame(score_test %>%
                      dplyr::filter(process %in% c("Auxin", "Brassinosteroid", "High sugar", "Defense",
                                                   "Strigolactone", "Jasmonic acid")) %>%
                      dplyr::select(stage, process, n_gene, diff_mix_pure, lmm_p, lmm_fdr, perm_p) %>%
                      dplyr::arrange(stage, process)), row.names = FALSE)

cat("[03] 分时期statistics：\n")
print(as.data.frame(score_test %>%
                      dplyr::group_by(stage) %>%
                      dplyr::summarise(n_sig_lmm = sum(lmm_fdr < 0.05),
                                       n_sig_perm = sum(perm_fdr < 0.05),
                                       min_p = min(lmm_p), .groups = "drop")))

bioRtools::write_data(
  list(pathway_scores = score_tab, pathway_tests = score_test),
  file = path_out
)
cat("[03] 输出：", path_out, "\n")

# ---- 图：pathway_scores差异森林图 ----
score_test %>%
  dplyr::mutate(
    stage = factor(stage, levels = c("tillering", "filling")),
    sig = dplyr::case_when(
      lmm_fdr < 0.05 ~ "FDR < 0.05",
      lmm_p < 0.05 ~ "P < 0.05",
      TRUE ~ "n.s."
    ),
    sig = factor(sig, levels = c("FDR < 0.05", "P < 0.05", "n.s."))
  ) %>%
  ggplot2::ggplot(aes(x = lmm_estimate, y = forcats::fct_reorder(process, lmm_estimate), color = sig)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  ggplot2::geom_errorbarh(aes(xmin = lmm_estimate - 1.96 * lmm_se,
                              xmax = lmm_estimate + 1.96 * lmm_se),
                          height = 0.25, linewidth = 0.5) +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::facet_wrap(~stage, scales = "free_x") +
  ggplot2::scale_color_manual(values = c("FDR < 0.05" = "#a63446", "P < 0.05" = "#e0a458", "n.s." = "grey55")) +
  ggplot2::labs(
    x = "Pathway score difference (Mix - Pure, mean Z)",
    y = NULL,
    color = NULL,
    caption = "Pathway score = mean Z-score of genes in the pathway; LMM with line random intercept; BH-FDR across 20 pathways per stage"
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 8),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.caption = ggplot2::element_text(size = 7, hjust = 0)
  ) -> p.s10

bioRtools::write_data(p.s10, file = file.path(out_dir, "03_FigureS10_pathway_scores森林图.pdf"),
                      width = 10, height = 7, units = "in")
bioRtools::write_data(p.s10, file = file.path(out_dir, "03_FigureS10_pathway_scores森林图.png"),
                      width = 10, height = 7, units = "in", dpi = 500)
cat("[03] 输出：03_FigureS10_pathway_scores森林图.pdf / .png\n")

cat("\n[03] sessionInfo()\n")
sessionInfo()
