# ============================================================
# manuscript / source_data / build_source_data.R
# 作用：把正文图（Fig. 1–5）与补充图（Fig. S1–S19）**每个子图**的绘图底层数据
#       汇总到一本工作簿 manuscript/source_data/Source_Data.xlsx —— **一个子图一张 sheet**。
#
# sheet 命名：多面板图 = 图号 + 面板小写字母（Fig1a、Fig1b…、FigS2a…）；
#            单面板图     = 图号（FigS3、FigS4…，不带字母）。
# sheet 布局：第 1 行 = 元信息（figure_panel / code_file / generated / software / note）；
#            第 2 行 = 空行；第 3 行 = 表头；第 4 行起 = 数据。
#            （manuscript/code/*.R 一律用 readxl::read_excel(path, sheet = "Fig1a", skip = 2) 读入）
# 运行：/opt/homebrew/bin/Rscript "<root>/manuscript/source_data/build_source_data.R"
# 说明：本脚本只**读取** result/ 与 data/ 下既有产物并整理成表，不重跑任何分析、不改动上游文件。
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
dir_out <- file.path(project_root, "manuscript", "source_data")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
xlsx_out <- file.path(dir_out, "Source_Data.xlsx")

rd <- function(path, sheet = NULL) {
  p <- file.path(project_root, path)
  if (!file.exists(p)) stop("缺文件：", p)
  if (is.null(sheet)) bioRtools::read_data(p) else bioRtools::read_data(p, sheet = sheet)
}
raw <- function(path, sheet) {
  readxl::read_excel(file.path(project_root, path), sheet = sheet, .name_repair = "minimal") %>%
    as.data.frame()
}

# ---------------- 基础表 ----------------
traits <- raw("data/01_分析输入表/03_分蘖数_株高_单株产量.xlsx", 1)
planned <- rd("result/02_农艺性状/13_表_计划对比与原始尺度效应.xlsx", "planned_contrasts_traits")
idx_con <- rd("result/02_农艺性状/13_表_计划对比与原始尺度效应.xlsx", "height_index_contrasts")
idx_long <- rd("result/02_农艺性状/01_表_株高变化指数.xlsx", "index_long")
pca_scores <- rd("result/02_农艺性状/03_表_表型PCA.xlsx", "scores")
pca_var <- rd("result/02_农艺性状/03_表_表型PCA.xlsx", "variance_explained")
gd_tab <- rd("result/02_农艺性状/05_表_遗传距离与表型相关.xlsx", "data_joined")
gd_cor <- rd("result/02_农艺性状/05_表_遗传距离与表型相关.xlsx", "correlation")
kgw_emm <- rd("result/02_农艺性状/12_表_田间株高与千粒重.xlsx", "emmeans")
kgw_mod <- rd("result/02_农艺性状/12_表_田间株高与千粒重.xlsx", "overall_model")

mat_combo <- raw("result/03_基因表达总体/23_表_表达矩阵_三时期_不剔离群值.xlsx", "combo_matrix_median")
gene_anno <- rd("result/03_基因表达总体/23_表_表达矩阵_三时期_不剔离群值.xlsx", "gene_annotation")
pls_scores <- rd("result/03_基因表达总体/24_表_PLSDA_三时期_中位数矩阵.xlsx", "scores")
pls_axis <- rd("result/03_基因表达总体/24_表_PLSDA_三时期_中位数矩阵.xlsx", "axis_tests")
pls_vip <- rd("result/03_基因表达总体/24_表_PLSDA_三时期_中位数矩阵.xlsx", "vip_top10")
pls_load <- rd("result/03_基因表达总体/24_表_PLSDA_三时期_中位数矩阵.xlsx", "loadings")
pls_var <- rd("result/03_基因表达总体/24_表_PLSDA_三时期_中位数矩阵.xlsx", "variance_explained")
limma_stage <- rd("result/04_通路基因表达/23_表_limma统一口径统计_三时期.xlsx", "gene_stage_level")
limma_focal <- rd("result/04_通路基因表达/23_表_limma统一口径统计_三时期.xlsx", "gene_stage_focal_level")
def_score <- rd("result/04_通路基因表达/23_表_limma统一口径统计_三时期.xlsx", "defense_score")
path_def <- rd("result/04_通路基因表达/23_表_limma统一口径统计_三时期.xlsx", "pathway_definition")
sugar_cor <- rd("result/04_通路基因表达/33_表_TableS5_高糖防御相关性_三时期_稳健口径.xlsx", "correlation")
sugar_det <- rd("result/04_通路基因表达/33_表_TableS5_高糖防御相关性_三时期_稳健口径.xlsx", "detection")
stage_con <- rd("result/04_通路基因表达/29_表_阶段特异性比较_三时期.xlsx", "stage_contrasts")
stage_cnt <- rd("result/04_通路基因表达/29_表_阶段特异性比较_三时期.xlsx", "counts")

inoc_z <- raw("result/01_接种与田间病害/02_表_接种Zscore中间表.xlsx", "Sheet1")
inoc_let <- rd("result/01_接种与田间病害/03_表_接种各focal组合emmeans字母.xlsx", "emmean_letters")
field_main <- rd("result/01_接种与田间病害/04_表_田间病害主要GLMM.xlsx", "emmean_letters_main")
field_all <- rd("result/01_接种与田间病害/05_表_田间病害补充_Smut与36_74.xlsx", "emmean_letters_all")

unit_tab <- rd("result/06_设计分解与稳健性/01_表_接种_盆内聚类与设计分解.xlsx", "inference_unit")
decomp_tab <- rd("result/06_设计分解与稳健性/01_表_接种_盆内聚类与设计分解.xlsx", "design_decomposition")
recip_pairs <- rd("result/06_设计分解与稳健性/01_表_接种_盆内聚类与设计分解.xlsx", "reciprocity_pairs")
recip_sum <- rd("result/06_设计分解与稳健性/01_表_接种_盆内聚类与设计分解.xlsx", "reciprocity_summary")
growth_tab <- rd("result/06_设计分解与稳健性/02_表_株高_生长曲线与起始点.xlsx", "time_point_tests")
cp_grid <- rd("result/06_设计分解与稳健性/05_表_株高_效应形状与变化点.xlsx", "change_point_grid")
shape_tab <- rd("result/06_设计分解与稳健性/05_表_株高_效应形状与变化点.xlsx", "shape_comparison")
boot_runs <- rd("result/06_设计分解与稳健性/05_表_株高_效应形状与变化点.xlsx", "bootstrap_runs")
boot_sum <- rd("result/06_设计分解与稳健性/05_表_株高_效应形状与变化点.xlsx", "bootstrap_summary")
year_tab <- rd("result/06_设计分解与稳健性/07_表_田间_年份与田块稳定性.xlsx", "by_year")
year_main <- rd("result/06_设计分解与稳健性/07_表_田间_年份与田块稳定性.xlsx", "main_model")

horm_class <- rd("result/05_激素_迈维/01_表_检出情况汇总.xlsx", "class_level")
horm_cpd <- rd("result/05_激素_迈维/01_表_检出情况汇总.xlsx", "compound_level")
horm_desc <- rd("result/05_激素_迈维/02_表_全检出物质组间比较.xlsx", "descriptive")
horm_tests <- rd("result/05_激素_迈维/02_表_全检出物质组间比较.xlsx", "tests")

path_tests <- rd("result/08_投稿规格图/03_表_通路分数检验_三时期.xlsx", "pathway_tests")
eta_gene <- rd("result/08_投稿规格图/04_表_伙伴效应方差分解_三时期.xlsx", "gene_level")
eta_top <- rd("result/08_投稿规格图/04_表_伙伴效应方差分解_三时期.xlsx", "top_interaction_genes")
eta_perm <- rd("result/08_投稿规格图/04_表_伙伴效应方差分解_三时期.xlsx", "permutation_reference")
vip_gd <- rd("result/08_投稿规格图/04b_表_VIP与遗传距离相关_三时期.xlsx", "data_joined")
vip_cor <- rd("result/08_投稿规格图/04b_表_VIP与遗传距离相关_三时期.xlsx", "correlation")

stage_levels <- c("tillering", "flowering", "filling")

# ---------------- 派生：热图（与出图脚本同一套计算）----------------
meta_cols <- c("combo", "stage", "focal", "neighbor", "group")
gene_cols <- setdiff(colnames(mat_combo), meta_cols)

scale01_by_row <- function(v) {
  rng <- range(v, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0, length(v)))
  (v - rng[1]) / diff(rng)
}

# Fig2a：通路水平热图（每基因在 stage 内跨 48 组合 Z 标准化 → 通路内取均值）
mat_gene <- as.matrix(mat_combo[, gene_cols])
z_stage <- mat_gene
for (s in stage_levels) {
  idx <- which(mat_combo$stage == s)
  z_stage[idx, ] <- scale(mat_gene[idx, , drop = FALSE])
}
path_levels <- sort(unique(gene_anno$process))
path_mat <- vapply(path_levels, function(pp) {
  gs <- intersect(gene_anno$gene[gene_anno$process == pp], colnames(z_stage))
  rowMeans(z_stage[, gs, drop = FALSE], na.rm = TRUE)
}, numeric(nrow(z_stage)))
colnames(path_mat) <- path_levels
combo_order <- mat_combo %>%
  dplyr::distinct(focal, neighbor, group) %>%
  dplyr::mutate(combo_label = paste0(focal, "-", neighbor)) %>%
  dplyr::arrange(group, focal, neighbor) %>%
  dplyr::pull(combo_label)
fig2a <- as.data.frame(path_mat) %>%
  # 2026-10-19 修：path_mat 无行名，rownames_to_column() 会给出 "1".."48"，
  #   与 mat_combo$combo（"tillering_15_15" 形式）join 不上 → stage/focal/neighbor/group 全 NA。
  dplyr::mutate(combo = mat_combo$combo) %>%
  tidyr::pivot_longer(cols = -combo, names_to = "pathway", values_to = "mean_z") %>%
  dplyr::left_join(mat_combo[, meta_cols], by = "combo") %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels),
                combo_label = paste0(focal, "-", neighbor),
                Panel = "a")

# Fig3a：防御基因热图（与 code/…/16_Figure3A_三时期_投稿规格.R 同口径：按 combo 行做 0-1 缩放）
def_genes <- gene_anno %>% dplyr::filter(process %in% c("Defense", "Jasmonic acid", "Salicylic acid")) %>% dplyr::pull(gene)
mat_def <- as.matrix(mat_combo[, intersect(def_genes, gene_cols)])
rownames(mat_def) <- mat_combo$combo
mat_def_scaled <- t(apply(mat_def, 1, scale01_by_row))
fig3a <- as.data.frame(mat_def_scaled) %>%
  tibble::rownames_to_column(var = "combo") %>%
  tidyr::pivot_longer(cols = -combo, names_to = "gene", values_to = "scaled_expression") %>%
  dplyr::left_join(mat_combo[, meta_cols], by = "combo") %>%
  dplyr::left_join(gene_anno %>% dplyr::select(gene, process), by = "gene") %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels),
                combo_label = paste0(focal, "-", neighbor),
                Panel = "a")

# FigS18：全部 92 基因热图（每行 0-1 缩放）
gene_order <- pls_load %>%
  dplyr::left_join(gene_anno %>% dplyr::select(gene, process), by = "gene") %>%
  dplyr::arrange(process, p1) %>% dplyr::pull(gene) %>% unique()
mat_heat <- as.matrix(mat_combo[, gene_cols])
rownames(mat_heat) <- mat_combo$combo
mat_heat_scaled <- t(apply(mat_heat, 1, scale01_by_row))
figS18 <- as.data.frame(mat_heat_scaled) %>%
  tibble::rownames_to_column(var = "combo") %>%
  tidyr::pivot_longer(cols = -combo, names_to = "gene", values_to = "scaled_expression") %>%
  dplyr::left_join(mat_combo[, meta_cols], by = "combo") %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels),
                gene = factor(gene, levels = rev(gene_order)),
                combo_label = paste0(focal, "-", neighbor))

# ---------------- 组装面板表 ----------------
meta_row <- function(panel, code_file, note = "") {
  c(paste0("figure_panel=", panel), paste0("code_file=", code_file),
    "generated=2026-10-18",
    "software=R 4.6.1; ggplot2 4.0.3; bioRtools",
    paste0("note=", note))
}

sheets <- list()
add <- function(name, df, code_file, note = "") {
  sheets[[name]] <<- list(meta = meta_row(name, code_file, note), data = as.data.frame(df))
}

# ---- Fig. 1 ----
add("Fig1a", idx_con %>% dplyr::mutate(Panel = "a"),
    "Fig1.R", "Height change index contrasts (15 surveys as unit); bars/CI in the figure")
add("Fig1b", traits %>% dplyr::filter(phenotype == "tiller number") %>%
      dplyr::left_join(planned %>% dplyr::filter(phenotype == "tiller number") %>%
                         dplyr::transmute(focal, neighbor = as.character(mix_neighbor),
                                          p_value_mix_vs_pure = p_value, p_fdr_mix_vs_pure = p_fdr),
                       by = c("focal", "neighbor")) %>%
      dplyr::mutate(Panel = "b"),
    "Fig1.R", "Per-plant tiller number + planned contrast P (mix vs its own monoculture, Dunnett)")
add("Fig1c", traits %>% dplyr::filter(phenotype == "yield") %>%
      dplyr::left_join(planned %>% dplyr::filter(phenotype == "yield") %>%
                         dplyr::transmute(focal, neighbor = as.character(mix_neighbor),
                                          p_value_mix_vs_pure = p_value, p_fdr_mix_vs_pure = p_fdr),
                       by = c("focal", "neighbor")) %>%
      dplyr::mutate(Panel = "c"),
    "Fig1.R", "Per-plant yield + planned contrast P (mix vs its own monoculture, Dunnett)")
add("Fig1d", pca_scores %>%
      dplyr::left_join(pca_var %>% dplyr::transmute(component, variance_explained) %>%
                         tidyr::pivot_wider(names_from = component, values_from = variance_explained,
                                            names_prefix = "var_"), by = character(0)) %>%
      dplyr::mutate(Panel = "d"),
    "Fig1.R", "Phenotypic PCA scores of the 16 focal-neighbour combinations (PC1/PC2)")
add("Fig1e", dplyr::bind_rows(
  idx_long %>% dplyr::mutate(focal = as.character(focal), neighbor = as.character(neighbor)) %>%
    dplyr::group_by(focal, neighbor, time) %>%
    dplyr::summarise(value = mean(index, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(trait = "Plant height change index (cm)", Panel = "e"),
  traits %>% dplyr::filter(phenotype %in% c("tiller number", "yield")) %>%
    dplyr::transmute(focal = as.character(focal), neighbor = as.character(neighbor), value,
                     trait = dplyr::if_else(phenotype == "tiller number",
                                            "Tiller number (per plant)", "Yield (g per plant)"),
                     Panel = "e")) %>%
    dplyr::mutate(group = dplyr::if_else(focal == neighbor, "Pure", "Mix")),
  "Fig1.R", "Three traits, pure vs mixture (height index per survey mean; tiller/yield per plant)")

# ---- Fig. 2 ----
add("Fig2a", fig2a, "Fig2.R", "Pathway-level heatmap: per-gene Z within stage across 48 combinations, then pathway mean")
add("Fig2b", pls_scores %>% dplyr::mutate(Panel = "b"), "Fig2.R", "PLS-DA score plot (component 1 vs 2)")
add("Fig2c", pls_scores %>% dplyr::select(combo, stage, focal, neighbor, group, p1) %>%
      dplyr::rename(score = p1) %>%
      dplyr::left_join(pls_axis %>% dplyr::filter(axis == "Dim 1", test_group == "within_stage",
                                                  test == "wilcoxon") %>%
                         dplyr::transmute(stage = stringr::str_remove(comparison, ": Mix vs Pure"),
                                          wilcoxon_p = p_value), by = "stage") %>%
      dplyr::mutate(Panel = "c",
                    kruskal_p = pls_axis$p_value[pls_axis$axis == "Dim 1" &
                                                   pls_axis$test == "kruskal_wallis"][1]),
    "Fig2.R", "Component 1 score by stage and planting mode (+ within-stage Wilcoxon P)")
add("Fig2d", pls_scores %>% dplyr::select(combo, stage, focal, neighbor, group, p2) %>%
      dplyr::rename(score = p2) %>%
      dplyr::left_join(pls_axis %>% dplyr::filter(axis == "Dim 2", test_group == "within_stage",
                                                  test == "wilcoxon") %>%
                         dplyr::transmute(stage = stringr::str_remove(comparison, ": Mix vs Pure"),
                                          wilcoxon_p = p_value), by = "stage") %>%
      dplyr::mutate(Panel = "d",
                    kruskal_p = pls_axis$p_value[pls_axis$axis == "Dim 2" &
                                                   pls_axis$test == "kruskal_wallis"][1]),
    "Fig2.R", "Component 2 score by stage and planting mode (+ within-stage Wilcoxon P)")
add("Fig2e", pls_load %>% dplyr::filter(gene %in% pls_vip$gene) %>%
      dplyr::left_join(gene_anno %>% dplyr::select(gene, process), by = "gene") %>%
      dplyr::left_join(pls_vip, by = "gene") %>%
      dplyr::mutate(Panel = "e"),
    "Fig2.R", "Loading plot of the ten highest-VIP genes")

# ---- Fig. 3 ----
add("Fig3a", fig3a, "Fig3.R", "Defence/JA/SA gene heatmap (values scaled 0-1 within each combination)")
add("Fig3b", {
  sample_mat <- raw("result/03_基因表达总体/11_表_表达矩阵_三时期.xlsx", "sample_matrix")
  def_gene_score <- gene_anno %>% dplyr::filter(process == "Defense") %>% dplyr::pull(gene)
  meta_sample <- c("stage", "sample", "focal", "neighbor")
  gene_sample <- setdiff(colnames(sample_mat), meta_sample)
  sample_mat %>%
    dplyr::mutate(group = factor(dplyr::if_else(focal == neighbor, "Pure", "Mix"),
                                 levels = c("Pure", "Mix"))) %>%
    dplyr::group_by(stage) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(gene_sample), ~ as.vector(scale(.x)))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(defence_score = rowMeans(dplyr::across(dplyr::all_of(intersect(def_gene_score, gene_sample))))) %>%
    dplyr::select(stage, sample, focal, neighbor, group, defence_score) %>%
    dplyr::left_join(def_score %>% dplyr::select(stage, n_sample, diff_mix_minus_pure, se, t_value, p_value),
                     by = "stage") %>%
    dplyr::mutate(Panel = "b")
},
  "Fig3.R", "Defence pathway score per sample (mean Z of the 12 defence genes) + stage-wise test statistics")

# ---- Fig. 4 ----
add("Fig4a", limma_stage %>% dplyr::filter(process == "Defense") %>% dplyr::mutate(Panel = "a"),
    "Fig4.R", "Stage-level limma, defence genes (moderated log2 fold change)")
add("Fig4b", limma_stage %>% dplyr::filter(process == "High sugar") %>% dplyr::mutate(Panel = "b"),
    "Fig4.R", "Stage-level limma, high-sugar genes")

# ---- Fig. 5 ----
add("Fig5a", decomp_tab %>% dplyr::filter(term %in% c("focal", "neighbor", "focal:neighbor")) %>%
      dplyr::mutate(Panel = "a"), "Fig5.R", "Pot-level decomposition of lesion density (F values)")
add("Fig5b", recip_pairs %>% dplyr::mutate(Panel = "b") %>%
      dplyr::left_join(recip_sum %>% dplyr::select(disease, perm_p_symmetry), by = "disease"),
    "Fig5.R", "Reciprocity: effect on focal vs effect on partner (12 unordered pairs)")
add("Fig5c", year_tab %>% dplyr::mutate(Panel = "c") %>%
      dplyr::left_join(year_main %>% dplyr::select(category, p_interaction), by = "category"),
    "Fig5.R", "Year-by-year field effects (diseases: odds ratio; traits: mean difference)")

# ---- Fig. S2 ----
add("FigS2a", inoc_z %>% dplyr::transmute(disease, group, batch = as.character(batch),
                                          focal = as.character(focal), neighbor = as.character(neighbor),
                                          lesion_density_per_cm2 = lesion_number_per_cm2, Panel = "a"),
    "FigS2.R", "Raw lesion density per leaf (greenhouse inoculation)")
add("FigS2b", inoc_z %>% dplyr::group_by(disease, group) %>%
      dplyr::mutate(z = as.vector(scale(lesion_number_per_cm2))) %>% dplyr::ungroup() %>%
      dplyr::transmute(disease, group, batch = as.character(batch), focal = as.character(focal),
                       neighbor = as.character(neighbor), lesion_density_z = z, Panel = "b"),
    "FigS2.R", "Lesion density after Z-score standardisation within disease and batch")
for (fc in c("15", "49", "69", "87")) {
  nm <- paste0("FigS2", c("15" = "c", "49" = "d", "69" = "e", "87" = "f")[[fc]])
  # C–F：每个 focal 的邻体比较——散点（每片叶的 Z 标准化病斑密度）+ 箱线 + Sidak 字母
  add(nm, inoc_z %>%
        dplyr::group_by(disease, group) %>%
        dplyr::mutate(lesion_norm = as.vector(scale(lesion_number_per_cm2))) %>%
        dplyr::ungroup() %>%
        dplyr::filter(as.character(focal) == fc) %>%
        dplyr::transmute(disease, focal = as.character(focal), neighbor = as.character(neighbor),
                         group, batch = as.character(batch), lesion_norm) %>%
        dplyr::left_join(inoc_let %>% dplyr::filter(as.character(focal) == fc) %>%
                           dplyr::transmute(disease, neighbor = as.character(neighbor),
                                            emmean, signif),
                         by = c("disease", "neighbor")) %>%
        dplyr::group_by(disease, neighbor) %>%
        dplyr::mutate(y_letter = max(lesion_norm, na.rm = TRUE) +
                        0.08 * diff(range(lesion_norm, na.rm = TRUE))) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(Panel = substr(nm, 6, 6)),
      "FigS2.R", paste0("Neighbour comparison within focal accession ", fc,
                        ": per-leaf Z-standardised lesion density + Sidak letters"))
}
add("FigS2g", field_main %>% dplyr::mutate(Panel = "g"), "FigS2.R",
    "Field disease incidence, four focal accessions (binomial GLMM)")
add("FigS2h", field_all %>% dplyr::mutate(Panel = "h"), "FigS2.R",
    "Field disease incidence, all six accessions including smut")

# ---- Fig. S3 ----
add("FigS3", vip_gd %>%
      dplyr::left_join(vip_cor %>% dplyr::select(stage, gene, n, slope, r, r_squared, p_value, p_fdr),
                       by = c("stage", "gene")) %>%
      dplyr::mutate(Panel = ""),
    "FigS3.R",
    "Expression of the ten highest-VIP genes vs pairwise genetic distance (12 mixture combinations) + Pearson correlation stats")

# ---- Fig. S4 ----
add("FigS4", kgw_emm %>% dplyr::mutate(Panel = "",
                                       p_overall = kgw_mod$p_value[match(category, kgw_mod$category)]),
    "FigS4.R", "Field plant height and 1000-grain weight (emmeans with 95% CI)")

# ---- Fig. S5 ----
add("FigS5", sugar_cor %>% dplyr::transmute(stage, sugar_gene, defense_gene, n, rho, p_value, p_fdr),
    "FigS5.R", "Spearman correlation between high-sugar and defence genes (detection >= 50%)")

# ---- Fig. S6 ----
add("FigS6a", horm_class %>% dplyr::mutate(Panel = "a"), "FigS6.R",
    "Hormone detection rate by compound class")
add("FigS6b", horm_cpd %>% dplyr::filter(n_detected > 0) %>% dplyr::mutate(Panel = "b"), "FigS6.R",
    "Detection rate of the 38 detected compounds")

# ---- Fig. S7 ----
add("FigS7", {
  cpd_keep <- horm_desc %>% dplyr::distinct(index) %>% dplyr::pull(index)
  horm_raw <- raw("data/01_分析输入表/05_激素_迈维含量表.xlsx", 1)
  horm_raw %>%
    dplyr::filter(index %in% cpd_keep) %>%
    dplyr::transmute(index, class_en, sample, group, focal, replicate, value) %>%
    dplyr::left_join(horm_desc %>% dplyr::select(index, group, n, mean, sd, median),
                     by = c("index", "group")) %>%
    dplyr::left_join(horm_tests %>% dplyr::select(index, kruskal_p), by = "index")
},
  "FigS7.R", "Fully detected compounds: per-sample values by group (+ descriptive stats and Kruskal-Wallis P)")

# ---- Fig. S8 ----
add("FigS8a", unit_tab %>% dplyr::mutate(Panel = "a"), "FigS8.R",
    "Inference-unit comparison: per plant vs pot random effect vs pot aggregation")
add("FigS8b", recip_pairs %>% dplyr::mutate(Panel = "b") %>%
      dplyr::left_join(recip_sum %>% dplyr::select(disease, perm_p_symmetry), by = "disease"),
    "FigS8.R", "Reciprocity of interaction effects between partners")

# ---- Fig. S9 ----
add("FigS9a", growth_tab %>% dplyr::transmute(day, n_pure, n_mix, mean_pure, mean_mix,
                                              diff_mix_pure, SE = NA_real_, Panel = "a"),
    "FigS9.R", "Plant height trajectories: mean height of pure vs mixture by survey day")
add("FigS9b", growth_tab %>% dplyr::transmute(day, t_p_value, t_fdr, Panel = "b"), "FigS9.R",
    "Mix vs Pure contrast P value by survey day")

# ---- Fig. S10 ----
add("FigS10", path_tests, "FigS10.R",
    "Pathway score differences (LMM estimate with 95% CI) at three developmental stages")

# ---- Fig. S11 ----
add("FigS11a", eta_gene %>% dplyr::select(stage, gene, eta2_focal, eta2_neighbor, eta2_interaction) %>%
      dplyr::mutate(Panel = "a"), "FigS11.R", "Partial eta-squared per gene and stage")
add("FigS11b", eta_gene %>% dplyr::select(stage, gene, eta2_focal, eta2_neighbor, eta2_interaction) %>%
      dplyr::mutate(Panel = "b"), "FigS11.R", "Own vs partner genotype effect (scatter)")
add("FigS11c", eta_top %>% dplyr::mutate(Panel = "c"), "FigS11.R",
    "Genes ranked by interaction eta-squared (top 30)")

# ---- Fig. S12 ----
add("FigS12a", {
  tau_hat <- cp_grid$tau[which.max(cp_grid$logLik)]
  ci_set <- cp_grid$tau[cp_grid$in_ci95 %in% c(TRUE, "TRUE")]
  traits %>%
    dplyr::filter(phenotype == "plant height", !is.na(time)) %>%
    dplyr::mutate(group = factor(group, levels = c("Pure", "Mix")),
                  date = as.Date(as.character(time), format = "%Y%m%d"),
                  day = as.numeric(date - min(date))) %>%
    dplyr::group_by(day, group) %>%
    dplyr::summarise(mean_value = mean(value), se = stats::sd(value) / sqrt(dplyr::n()), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = group, values_from = c(mean_value, se)) %>%
    dplyr::mutate(diff_mix_minus_pure = mean_value_Mix - mean_value_Pure,
                  se_diff = sqrt(se_Mix^2 + se_Pure^2),
                  tau_hat = tau_hat,
                  ci95_low = min(ci_set), ci95_high = max(ci_set),
                  Panel = "a")
},
  "FigS12.R", "Plant height difference (Mix - Pure) by survey day with 95% CI and change-point estimate")
add("FigS12b", cp_grid %>% dplyr::select(tau, logLik, aic, b_group, p_group, in_ci95) %>%
      dplyr::mutate(Panel = "b"), "FigS12.R", "Profile likelihood over candidate change points")
add("FigS12c", boot_runs %>% dplyr::mutate(Panel = "c") %>%
      dplyr::left_join(boot_sum %>% dplyr::transmute(pct_step, pct_linear, pct_constant), by = character(0)),
    "FigS12.R", "Bootstrap distribution of change points (199 resamples)")

# ---- Fig. S13 ----
add("FigS13", limma_stage %>% dplyr::transmute(stage, gene, limma_logFC, limma_p, limma_fdr, limma_b),
    "FigS13.R", "Volcano: moderated log2 fold change vs -log10 P for 92 genes at three stages")

# ---- Fig. S14 ----
add("FigS14", year_tab %>% dplyr::mutate(Panel = "") %>%
      dplyr::left_join(year_main %>% dplyr::select(category, p_interaction), by = "category"),
    "FigS14.R", "Year-to-year stability of field effects")

# ---- Fig. S15 ----
add("FigS15a", stage_con %>% dplyr::transmute(gene, stage = "tillering", logFC = logFC_tillering,
                                              Panel = "a") %>%
      dplyr::bind_rows(stage_con %>% dplyr::transmute(gene, stage = "flowering", logFC = logFC_flowering,
                                                      Panel = "a")) %>%
      dplyr::bind_rows(stage_con %>% dplyr::transmute(gene, stage = "filling", logFC = logFC_filling,
                                                      Panel = "a")),
    "FigS15.R", "Stage-specific limma log2 fold changes (20 largest-magnitude genes)")
add("FigS15b", stage_cnt %>% dplyr::mutate(Panel = "b"), "FigS15.R",
    "Number of candidate-responsive genes per contrast")

# ---- Fig. S16 ----
add("FigS16", {
  cor_wide <- gd_cor %>%
    dplyr::filter(dataset == "Mix only (n = 12)") %>%
    dplyr::transmute(trait,
                     r_lab = paste0("r_mean_", stringr::str_remove(trait, "^mean_")),
                     p_lab = paste0("p_mean_", stringr::str_remove(trait, "^mean_"))) %>%
    tidyr::pivot_longer(cols = c(r_lab, p_lab), names_to = "kind", values_to = "colname")
  gd_tab %>%
    dplyr::mutate(
      r_height = gd_cor$r[gd_cor$dataset == "Mix only (n = 12)" & gd_cor$trait == "mean_height"],
      p_height = gd_cor$p_value[gd_cor$dataset == "Mix only (n = 12)" & gd_cor$trait == "mean_height"],
      r_tiller = gd_cor$r[gd_cor$dataset == "Mix only (n = 12)" & gd_cor$trait == "mean_tiller"],
      p_tiller = gd_cor$p_value[gd_cor$dataset == "Mix only (n = 12)" & gd_cor$trait == "mean_tiller"],
      r_yield = gd_cor$r[gd_cor$dataset == "Mix only (n = 12)" & gd_cor$trait == "mean_yield"],
      p_yield = gd_cor$p_value[gd_cor$dataset == "Mix only (n = 12)" & gd_cor$trait == "mean_yield"])
},
  "FigS16.R",
  "Genetic distance vs agronomic traits (16 combinations; the figure uses the 12 mixtures) + Pearson r/P")

# ---- Fig. S17 ----
add("FigS17", limma_stage %>% dplyr::filter(process == "Strigolactone") %>%
      dplyr::transmute(stage, gene, limma_logFC, limma_p, limma_fdr),
    "FigS17.R", "Strigolactone pathway genes across three developmental stages")

# ---- Fig. S18 ----
add("FigS18", figS18, "FigS18.R", "Gene-level heatmap: 92 genes x 48 combinations, scaled 0-1 per combination")

# ---- Fig. S19 ----
add("FigS19a", limma_focal %>% dplyr::filter(stage == "tillering", as.character(focal) == "69",
                                             process %in% c("Brassinosteroid", "Auxin"), gene != "OsXTR1") %>%
      dplyr::mutate(Panel = "a"), "FigS19.R",
    "Accession 69 at tillering: BR and auxin genes (limma within accession)")
add("FigS19b", limma_focal %>% dplyr::filter(stage == "filling", as.character(focal) == "87",
                                             process %in% c("Brassinosteroid", "Salicylic acid", "Drought"),
                                             gene != "OsXTR1") %>% dplyr::mutate(Panel = "b"),
    "FigS19.R", "Accession 87 at grain filling: BR, SA and drought genes")

# ---------------- 写出工作簿 ----------------
wb <- openxlsx::createWorkbook()
for (nm in names(sheets)) {
  openxlsx::addWorksheet(wb, nm)
  openxlsx::writeData(wb, nm, x = sheets[[nm]]$meta, startRow = 1, startCol = 1, colNames = FALSE)
  openxlsx::writeData(wb, nm, x = sheets[[nm]]$data, startRow = 3, startCol = 1, colNames = TRUE)
  openxlsx::freezePane(wb, nm, firstRow = TRUE)
}
openxlsx::saveWorkbook(wb, xlsx_out, overwrite = TRUE)
cat("[source_data] sheet 数：", length(sheets), "\n")
cat("[source_data] 输出：", xlsx_out, "\n")

# ---------------- 回读自检 ----------------
sheets_back <- readxl::excel_sheets(xlsx_out)
cat("[source_data] 回读 sheet 数：", length(sheets_back), "\n")
stopifnot(identical(sheets_back, names(sheets)))
for (nm in sheets_back) {
  d <- readxl::read_excel(xlsx_out, sheet = nm, skip = 2)
  cat(sprintf("  %-8s %5d 行 × %2d 列\n", nm, nrow(d), ncol(d)))
  stopifnot(nrow(d) >= 1)
}
cat("[source_data] 全部 sheet 至少一行数据 ✓\n")

cat("\n[source_data] sessionInfo()\n")
sessionInfo()
