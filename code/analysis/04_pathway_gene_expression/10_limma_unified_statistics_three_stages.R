# ============================================================
# 04_pathway_gene_expression / 10_limma_unified_statistics_three_stages.R
# 作用：把主题 05 的统一口径（limma + camera + 防御分数）扩展到**三个发育时期**
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_pathway_gene_expression/10_limma_unified_statistics_three_stages.R"
#
# 输入：<root>/result/03_overall_gene_expression/11_表_expression矩阵_three_stages.xlsx（sample_matrix、gene_annotation）
# 输出：<root>/result/04_pathway_gene_expression/23_表_limma_unified_statistics_three_stages.xlsx
#       6 sheet：gene_stage_level、gene_stage_focal_level、pathway_camera、defense_score、
#                key_genes、pathway_definition
# 其中 <root> = <repository_root>
#
# 口径（and脚本 07 完全一致，仅时期数 2 → 3）：
#   - 响应 = log2(expression + 1)；stage 层面design `~ 0 + group + focal`（focal 固定区组），
#     **显式构造 Mix − Pure 对比矩阵**（不可用 topTable(coef = "groupMix")，那是组均值）；
#   - stage × focal 层面在该 focal 内design `~ 0 + group`；
#   - BH-FDR 在每个test族内correction（stage 层面按时期 92 个gene；stage × focal 层面按 时期 × focal）；
#   - pathwaytest用 limma::camera（inter.gene.cor = 0.01）；
#   - 防御pathway_scores = Defense pathwaygene z 均值（每sample一个值）+ `lm(score ~ group + focal)`。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_mat <- file.path(project_root, "result", "03_overall_gene_expression", "11_表_expression矩阵_three_stages.xlsx")
out_dir <- file.path(project_root, "result", "04_pathway_gene_expression")
path_out <- file.path(out_dir, "23_表_limma_unified_statistics_three_stages.xlsx")

stopifnot(file.exists(path_mat))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_mat, "sample_matrix") -> mat.sample
bioRtools::read_data(path_mat, sheet = "gene_annotation") -> df.gene

meta_cols <- c("stage", "sample", "focal", "neighbor")
gene_cols <- setdiff(colnames(mat.sample), meta_cols)
stopifnot(length(gene_cols) == 92)

stage_levels <- c("tillering", "flowering", "filling")

mat.sample %>%
  dplyr::mutate(
    group = factor(dplyr::if_else(focal == neighbor, "Pure", "Mix"), levels = c("Pure", "Mix")),
    focal = factor(focal),
    stage = factor(stage, levels = stage_levels)
  ) -> mat.sample

pathway_list <- df.gene %>%
  dplyr::filter(gene %in% gene_cols) %>%
  dplyr::group_by(process) %>%
  dplyr::summarise(genes = list(intersect(gene, gene_cols)), n_gene = dplyr::n(), .groups = "drop")

pathway_index <- pathway_list$genes
names(pathway_index) <- pathway_list$process

cat("[10] sample ", nrow(mat.sample), "；gene ", length(gene_cols), "；pathway ", nrow(pathway_list), "\n")

# ---- 通用：design矩阵 + 显式对比 ----
make_design <- function(d, formula_rhs) {
  design <- stats::model.matrix(stats::as.formula(formula_rhs), data = d)
  colnames(design) <- make.names(colnames(design))
  design
}

mix_contrast <- function(design) {
  con <- matrix(0, nrow = ncol(design), ncol = 1,
                dimnames = list(colnames(design), "mix_vs_pure"))
  con["groupMix", "mix_vs_pure"] <- 1
  con["groupPure", "mix_vs_pure"] <- -1
  con
}

# ---- stage 层面 ----
run_stage <- function(st) {
  d <- mat.sample %>% dplyr::filter(stage == st)
  expr <- t(log2(as.matrix(d[, gene_cols]) + 1))

  design <- make_design(d, "~ 0 + group + focal")
  con <- mix_contrast(design)

  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr, design), con), robust = TRUE)
  top <- limma::topTable(fit, coef = "mix_vs_pure", number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::transmute(stage = st, gene = gene, limma_logFC = logFC, limma_ave_expr = AveExpr,
                     limma_t = t, limma_p = P.Value, limma_fdr = adj.P.Val, limma_b = B)

  cam <- limma::camera(expr, index = pathway_index, design = design,
                       contrast = con[, "mix_vs_pure"], inter.gene.cor = 0.01) %>%
    tibble::rownames_to_column(var = "process") %>%
    dplyr::mutate(stage = st, .before = 1) %>%
    dplyr::select(stage, process, n_gene = NGenes, direction = Direction,
                  camera_p = PValue, camera_fdr = FDR)

  expr_z <- t(scale(t(expr)))
  def_genes <- pathway_list$genes[[which(pathway_list$process == "Defense")]]
  score <- colMeans(expr_z[def_genes, , drop = FALSE])
  m_score <- stats::lm(score ~ group + focal, data = d)
  co <- summary(m_score)$coefficients
  score_tab <- tibble::tibble(
    stage = st,
    score_name = "Defense pathway score (mean Z)",
    n_gene = length(def_genes),
    n_sample = nrow(d),
    diff_mix_minus_pure = unname(stats::coef(m_score)["groupMix"]),
    se = unname(co[["groupMix", "Std. Error"]]),
    t_value = unname(co[["groupMix", "t value"]]),
    p_value = unname(co[["groupMix", "Pr(>|t|)"]])
  )

  list(gene = top, camera = cam, score = score_tab)
}

res_stage <- purrr::map(stage_levels, run_stage)

gene_stage <- purrr::map_dfr(res_stage, "gene") %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene") %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(limma_fdr = stats::p.adjust(limma_p, method = "BH")) %>%
  dplyr::ungroup()

camera_tab <- purrr::map_dfr(res_stage, "camera") %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(camera_fdr = stats::p.adjust(camera_p, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(stage, camera_p)

score_tab <- purrr::map_dfr(res_stage, "score")

cat("\n[10] stage 层面 limma（每期 92 gene内 BH）：\n")
print(as.data.frame(gene_stage %>%
                      dplyr::group_by(stage) %>%
                      dplyr::summarise(nominal_p005 = sum(limma_p < 0.05),
                                       fdr005 = sum(limma_fdr < 0.05),
                                       min_p = min(limma_p),
                                       max_abs_logFC = max(abs(limma_logFC)), .groups = "drop")))

cat("\n[10] camera（最小 FDR）：\n")
print(as.data.frame(camera_tab %>%
                      dplyr::group_by(stage) %>%
                      dplyr::summarise(min_p = min(camera_p), min_fdr = min(camera_fdr),
                                       n_fdr005 = sum(camera_fdr < 0.05), .groups = "drop")))

cat("\n[10] 防御pathway_scores：\n")
print(as.data.frame(score_tab))

# ---- stage × focal 层面 ----
run_stage_focal <- function(st, fo) {
  d <- mat.sample %>% dplyr::filter(stage == st, focal == fo)
  expr <- t(log2(as.matrix(d[, gene_cols]) + 1))

  design <- make_design(d, "~ 0 + group")
  con <- mix_contrast(design)

  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr, design), con), robust = TRUE)
  limma::topTable(fit, coef = "mix_vs_pure", number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::transmute(stage = st, focal = fo, gene = gene, focal_logFC = logFC,
                     focal_t = t, focal_p = P.Value, focal_fdr = adj.P.Val, n_pure = NA_integer_,
                     n_mix = NA_integer_)
}

gene_focal <- purrr::cross2(stage_levels, as.character(sort(unique(mat.sample$focal)))) %>%
  purrr::map_dfr(~ run_stage_focal(.x[[1]], .x[[2]])) %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene") %>%
  dplyr::group_by(stage, focal) %>%
  dplyr::mutate(focal_fdr = stats::p.adjust(focal_p, method = "BH")) %>%
  dplyr::ungroup()

cat("\n[10] stage × focal 层面（每 时期×focal 内 BH），名义显著计数：\n")
print(as.data.frame(gene_focal %>%
                      dplyr::group_by(stage, focal) %>%
                      dplyr::summarise(nominal_p005 = sum(focal_p < 0.05),
                                       fdr005 = sum(focal_fdr < 0.05), .groups = "drop")))

key_genes <- gene_focal %>%
  dplyr::filter(gene %in% c("OsBLE1", "OsBLE2", "OsGNS1", "OsGBSSII", "OsD14-1", "OsWRKY45",
                            "OsWRKY7", "OsGH3.2", "OsIAA1", "OsIAA4")) %>%
  dplyr::select(stage, focal, gene, process, focal_logFC, focal_p, focal_fdr) %>%
  dplyr::arrange(stage, focal, focal_p)

pathway_def <- pathway_list %>%
  tidyr::unnest_longer(genes, values_to = "gene") %>%
  dplyr::select(process, gene, n_gene)

bioRtools::write_data(
  list(
    gene_stage_level = gene_stage,
    gene_stage_focal_level = gene_focal,
    pathway_camera = camera_tab,
    defense_score = score_tab,
    key_genes = key_genes,
    pathway_definition = pathway_def
  ),
  file = path_out
)

cat("\n[10] 输出：", path_out, "\n")

cat("\n[10] sessionInfo()\n")
sessionInfo()
