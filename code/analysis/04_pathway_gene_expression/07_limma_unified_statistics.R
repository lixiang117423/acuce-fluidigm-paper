# ============================================================
# 04_pathway_gene_expression / 07_limma_unified_statistics.R
# 作用：把主题 04 的statistics口径统一为"focal 区组 + 经验贝叶斯（limma）"，并用
#       limma::camera 做**竞争性gene集test**replace原来的 Fisher 富集test
#       → 新的权威statistics表（Table S3/S4 and Figure 3B/3C/3D、4A/4B/4C 的data来源）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/04_pathway_gene_expression/07_limma_unified_statistics.R"
#
# 【为什么换口径】
#   1) 逐gene Wilcoxon：每个gene独立估方差 + 92 个gene做 BH correction，分蘖期仅 8 个gene FDR<0.05、灌浆期 0 个；
#   2) Fisher 富集：每条pathway只有 1–12 个gene（多数 2–4 个），2×2 表几乎无检出力，全部为阴性；
#   3) limma（focal 作固定区组 + 方差向共同方差收缩）会提升检出力，但**必须显式构造 Mix − Pure 对比**；
#   4) camera 是 limma 体系内的竞争性gene集test，会correctiongene间相关性，适合 2–15 个gene的小pathway。
#
# 【2026-09-22 更正（重要）】本脚本第一版曾直接用 `topTable(fit, coef = "groupMix")`，
#   但design矩阵是 `~ 0 + group + focal`——**在这种"无截距 + 单因子"的写法下，limma 的系数是组均值本身**
#   （groupMix = Mix 组在参考 focal 的对数均值），并不是 Mix − Pure 的差异。第一版因此报出
#   "分蘖期 53/92、灌浆期 59/92 个gene FDR<0.05、OsBLE1 FDR ≈ 1e-11"等results，实际test的是
#   "该gene是否expression"这一平凡假设，属错误结论。现已改为**显式构造对比矩阵**（groupMix − groupPure），
#   camera 亦改为传入该对比向量。更正后的正确results：**分蘖期and灌浆期均无gene达到 FDR < 0.05**
#   （分蘖期 14 个gene P<0.05，最强 OsGNS1 logFC = 1.15、P = 0.0024；灌浆期仅 1 个）。
# 因此本脚本输出的"统一口径"version用于替代 01/02 脚本；原 Wilcoxon/Fisher version保留不动，便于方法学对照。
#
# 输入：<root>/result/03_overall_gene_expression/01_表_expression矩阵.xlsx（sample_matrix、gene_annotation）
# 输出：<root>/result/04_pathway_gene_expression/14_表_limma_unified_statistics.xlsx
# 其中 <root> = <repository_root>
#
# 方法细节：响应 = log2(expression + 1)；stage 层面design `~ 0 + group + focal`（focal 固定区组）；
#   stage × focal 层面在该 focal 内design `~ 0 + group`；均用 limma::eBayes(robust = TRUE)。
#   camera 用 `inter.gene.cor = 0.01`（默认做法；面板gene数少、相关性估计不稳），并在表中record。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/04_pathway_gene_expression/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_mat <- file.path(project_root, "result", "03_overall_gene_expression", "01_表_expression矩阵.xlsx")
out_dir <- file.path(project_root, "result", "04_pathway_gene_expression")
path_out <- file.path(out_dir, "14_表_limma_unified_statistics.xlsx")

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
    group = factor(dplyr::if_else(focal == neighbor, "Pure", "Mix"), levels = c("Pure", "Mix")),
    focal = factor(focal)
  ) -> mat.sample

pathway_list <- df.gene %>%
  dplyr::filter(gene %in% gene_cols) %>%
  dplyr::group_by(process) %>%
  dplyr::summarise(genes = list(intersect(gene, gene_cols)), n_gene = dplyr::n(), .groups = "drop")

# camera 用 index 列表的**名字**作为gene集名称；不设名字会返回行号，故必须显式命名
pathway_index <- pathway_list$genes
names(pathway_index) <- pathway_list$process

cat("[07] gene ", length(gene_cols), " 个；pathway ", nrow(pathway_list), " 条（gene数 ",
    min(pathway_list$n_gene), "-", max(pathway_list$n_gene), "）\n")

# ---- stage 层面：limma（focal 区组）+ camera ----
run_stage <- function(st) {
  d <- mat.sample %>% dplyr::filter(stage == st)
  expr <- t(log2(as.matrix(d[, gene_cols]) + 1))

  design <- stats::model.matrix(~ 0 + group + focal, data = d)
  colnames(design) <- make.names(colnames(design))

  # 显式对比矩阵（必须）：groupMix − groupPure
  con <- matrix(0, nrow = ncol(design), ncol = 1,
                dimnames = list(colnames(design), "mix_vs_pure"))
  con["groupMix", "mix_vs_pure"] <- 1
  con["groupPure", "mix_vs_pure"] <- -1

  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr, design), con), robust = TRUE)
  top <- limma::topTable(fit, coef = "mix_vs_pure", number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::transmute(stage = st, gene = gene,
                     limma_logFC = logFC, limma_ave_expr = AveExpr,
                     limma_t = t, limma_p = P.Value, limma_fdr = adj.P.Val, limma_b = B)

  # camera：竞争性gene集test（correctiongene间相关）
  cam <- limma::camera(expr, index = pathway_index, design = design,
                       contrast = con[, "mix_vs_pure"], inter.gene.cor = 0.01) %>%
    tibble::rownames_to_column(var = "process") %>%
    dplyr::mutate(stage = st, .before = 1) %>%
    dplyr::select(stage, process, n_gene = NGenes, direction = Direction,
                  camera_p = PValue, camera_fdr = FDR)

  # 防御gene整体分数（用于 Figure 3B）：pathway内gene z 均值 → 线性模型（focal 固定区组）
  expr_z <- t(scale(t(expr)))
  def_genes <- pathway_list$genes[[which(pathway_list$process == "Defense")]]
  score <- colMeans(expr_z[def_genes, , drop = FALSE])
  m_score <- stats::lm(score ~ group + focal, data = d)
  score_tab <- tibble::tibble(
    stage = st,
    score_name = "Defense pathway score (mean Z)",
    n_gene = length(def_genes),
    diff_mix_minus_pure = unname(stats::coef(m_score)["groupMix"]),
    se = unname(summary(m_score)$coefficients["groupMix", "Std. Error"]),
    t_value = unname(summary(m_score)$coefficients["groupMix", "t value"]),
    p_value = unname(summary(m_score)$coefficients["groupMix", "Pr(>|t|)"])
  )

  list(gene = top, camera = cam, score = score_tab)
}

res_stage <- purrr::map(c("tillering", "filling"), run_stage)
gene_stage <- purrr::map_dfr(res_stage, "gene") %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene")
camera_tab <- purrr::map_dfr(res_stage, "camera") %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(camera_fdr = p.adjust(camera_p, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(stage, camera_p)
score_tab <- purrr::map_dfr(res_stage, "score")

cat("\n[07] 逐gene limma（FDR < 0.05）数量：\n")
print(as.data.frame(gene_stage %>% dplyr::group_by(stage) %>%
                      dplyr::summarise(n_gene = dplyr::n(),
                                       n_fdr_005 = sum(limma_fdr < 0.05), .groups = "drop")))

cat("\n[07] camera pathwaytest —— 分蘖期：\n")
print(as.data.frame(camera_tab %>% dplyr::filter(stage == "tillering") %>%
                      dplyr::select(process, n_gene, direction, camera_p, camera_fdr)), row.names = FALSE)

cat("\n[07] camera pathwaytest —— 灌浆期：\n")
print(as.data.frame(camera_tab %>% dplyr::filter(stage == "filling") %>%
                      dplyr::select(process, n_gene, direction, camera_p, camera_fdr)), row.names = FALSE)

cat("\n[07] 防御gene整体分数test（Figure 3B 的new_definition）：\n")
print(as.data.frame(score_tab), row.names = FALSE)

# ---- stage × focal 层面：limma（该 focal 内）----
run_focal <- function(st, f) {
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
  limma::topTable(fit, coef = "mix_vs_pure", number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::transmute(stage = st, focal = as.character(f),
                     gene = gene, focal_logFC = logFC, focal_t = t,
                     focal_p = P.Value, focal_fdr = adj.P.Val)
}

focal_grid <- tidyr::expand_grid(stage = c("tillering", "filling"), focal = c("15", "49", "69", "87"))
gene_focal <- purrr::pmap_dfr(list(focal_grid$stage, focal_grid$focal), run_focal) %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene")

cat("\n[07] stage × focal 层面（limma，FDR<0.05 gene数）：\n")
print(as.data.frame(gene_focal %>% dplyr::group_by(stage, focal) %>%
                      dplyr::summarise(n_fdr_005 = sum(focal_fdr < 0.05), .groups = "drop")) %>%
        tidyr::pivot_wider(names_from = focal, values_from = n_fdr_005))

# ---- 关键gene对照（69 分蘖期 / 87 灌浆期）----
key_genes <- c("OsBLE1", "OsBLE2", "OsGH3.2", "OsIAA1", "OsIAA4",
               "OsWRKY45", "OsWRKY7", "OsGBSSII", "OsHXK5", "OsHXK6", "OsGNS1", "OsD14-1")

key_tab <- dplyr::bind_rows(
  gene_focal %>% dplyr::filter(stage == "tillering", focal == "69", gene %in% key_genes) %>%
    dplyr::transmute(scope = "69 tillering", gene, process, logFC = focal_logFC,
                     p_value = focal_p, fdr = focal_fdr),
  gene_focal %>% dplyr::filter(stage == "filling", focal == "87", gene %in% key_genes) %>%
    dplyr::transmute(scope = "87 filling", gene, process, logFC = focal_logFC,
                     p_value = focal_p, fdr = focal_fdr)
) %>%
  dplyr::arrange(scope, p_value)

cat("\n[07] 关键gene（limma，按 focal 分层）：\n")
print(as.data.frame(key_tab), row.names = FALSE)

bioRtools::write_data(
  list(gene_stage_level = gene_stage,
       gene_stage_focal_level = gene_focal,
       pathway_camera = camera_tab,
       defense_score = score_tab,
       key_genes = key_tab,
       pathway_definition = pathway_list %>% tidyr::unnest(genes)),
  file = path_out
)
cat("[07] 输出：", path_out, "\n")

cat("\n[07] sessionInfo()\n")
sessionInfo()
