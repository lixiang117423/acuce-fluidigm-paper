# ============================================================
# 03_基因表达总体 / 11_盆号效应敏感性_三时期.R
# 作用：把**盆号接进模型**（审稿意见 M3：同一 Mix 盆两株共享盆环境），量化盆效应并做敏感性对照
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_基因表达总体/11_盆号效应敏感性_三时期.R"
#
# 输入：<root>/result/03_基因表达总体/11_表_表达矩阵_三时期.xlsx（sample_matrix）
#       <root>/result/03_基因表达总体/19_表_样品盆号映射_三时期.xlsx（sample_level：sample → pot）
# 输出：<root>/result/03_基因表达总体/20_表_盆号效应敏感性_三时期.xlsx
#       sheet：model_compare（逐基因两种模型对照）、stage_summary（每期汇总）、pot_structure
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 方法（每期分别做，92 个基因各拟合一次）：
#   ① 现行口径：`lm(log2(expr + 1) ~ group + focal)`（不含盆，等价于 limma 里的 focal 区组思路）；
#   ② 加盆随机截距：`lmerTest::lmer(log2(expr + 1) ~ group + focal + (1 | pot))`，
#      报告 group 的 Satterthwaite p 值与盆方差占比（ICC_pot）。
#   BH-FDR 在每个模型、每期内 92 个基因上校正；比较两种口径下 FDR < 0.05 的基因数（收紧/放宽/不变）。
# 说明：盆随机截距在盆内样本数 = 1–2 时信息有限，可能出现奇异拟合（singular fit），脚本会把
#   "lmer 是否报奇异/收敛警告"一并记录，便于判断结论稳健性。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_mat <- file.path(project_root, "result", "03_基因表达总体", "11_表_表达矩阵_三时期.xlsx")
path_pot <- file.path(project_root, "result", "03_基因表达总体", "19_表_样品盆号映射_三时期.xlsx")
out_dir <- file.path(project_root, "result", "03_基因表达总体")
path_out <- file.path(out_dir, "20_表_盆号效应敏感性_三时期.xlsx")

stopifnot(file.exists(path_mat), file.exists(path_pot))

stage_levels <- c("tillering", "flowering", "filling")

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_mat, "sample_matrix") -> mat.sample
read_sheet_raw(path_pot, "sample_level") -> pot.map

meta_cols <- c("stage", "sample", "focal", "neighbor")
gene_cols <- setdiff(colnames(mat.sample), meta_cols)
stopifnot(length(gene_cols) == 92)

d <- mat.sample %>%
  dplyr::mutate(stage = as.character(stage), sample = as.character(sample),
                group = factor(dplyr::if_else(focal == neighbor, "Pure", "Mix"),
                               levels = c("Pure", "Mix")),
                focal = factor(focal)) %>%
  dplyr::left_join(pot.map %>% dplyr::select(stage, sample, pot, partner_sample),
                   by = c("stage", "sample"))

stopifnot(sum(is.na(d$pot)) == 0)
d$pot <- factor(d$pot)

cat("[11] 样本 ", nrow(d), "；盆 ", dplyr::n_distinct(d$pot), "\n")

pot_structure <- d %>%
  dplyr::count(stage, pot, name = "n_sample") %>%
  dplyr::group_by(stage) %>%
  dplyr::summarise(n_pot = dplyr::n(),
                   n_pot_2_samples = sum(n_sample == 2),
                   n_pot_1_sample = sum(n_sample == 1), .groups = "drop")
print(as.data.frame(pot_structure))

fit_one_stage <- function(st) {
  dd <- d %>% dplyr::filter(stage == st)
  cat("\n[11] === ", st, "：n = ", nrow(dd), "，盆 = ", dplyr::n_distinct(dd$pot), "\n", sep = "")

  purrr::map_dfr(gene_cols, function(g) {
    y <- log2(dd[[g]] + 1)
    dat <- data.frame(y = y, group = dd$group, focal = dd$focal, pot = dd$pot)

    m_lm <- stats::lm(y ~ group + focal, data = dat)
    p_lm <- summary(m_lm)$coefficients["groupMix", "Pr(>|t|)"]

    # 奇异拟合在盆内 n = 1–2 时很常见：把警告静音但把 isSingular 标记记下来
    res <- tryCatch(
      withCallingHandlers({
        m_lmer <- lmerTest::lmer(y ~ group + focal + (1 | pot), data = dat)
        co <- summary(m_lmer)$coefficients
        vc <- as.data.frame(lme4::VarCorr(m_lmer))
        v_pot <- vc$vcov[vc$grp == "pot"][1]
        v_res <- vc$vcov[vc$grp == "Residual"][1]
        list(ok = TRUE, p = co["groupMix", "Pr(>|t|)"],
             icc = v_pot / (v_pot + v_res), singular = lme4::isSingular(m_lmer))
      }, warning = function(w) invokeRestart("muffleWarning")),
      error = function(e) list(ok = FALSE, p = NA_real_, icc = NA_real_, singular = NA)
    )

    tibble::tibble(
      stage = st, gene = g,
      p_lm_no_pot = p_lm,
      p_lmer_pot = res$p,
      icc_pot = res$icc,
      lmer_singular = res$singular,
      fit_ok = res$ok,
      logFC_lm = unname(stats::coef(m_lm)["groupMix"])
    )
  }) %>%
    dplyr::mutate(fdr_lm_no_pot = stats::p.adjust(p_lm_no_pot, method = "BH"),
                  fdr_lmer_pot = stats::p.adjust(p_lmer_pot, method = "BH"))
}

model_compare <- purrr::map_dfr(stage_levels, fit_one_stage)

# 每个基因在每个时期的检出率（表达量 > 0），用于把 "ICC = 1" 这类假象与真实盆效应区分开：
# 大量为 0 的基因其残差方差趋近 0，ICC 会被推到 1，但不代表存在盆效应。
detect_tab <- purrr::map_dfr(stage_levels, function(st) {
  dd <- d %>% dplyr::filter(stage == st)
  purrr::map_dfr(gene_cols, function(g) {
    v <- dd[[g]]
    tibble::tibble(stage = st, gene = g,
                   detection_rate = mean(v > 0, na.rm = TRUE))
  })
})

model_compare <- model_compare %>%
  dplyr::left_join(detect_tab, by = c("stage", "gene")) %>%
  dplyr::mutate(icc_reliable = detection_rate >= 0.5)

stage_summary <- model_compare %>%
  dplyr::group_by(stage) %>%
  dplyr::summarise(
    n_gene = dplyr::n(),
    n_singular = sum(lmer_singular %in% TRUE),
    n_fit_failed = sum(!fit_ok),
    median_icc_pot = median(icc_pot, na.rm = TRUE),
    max_icc_pot = max(icc_pot, na.rm = TRUE),
    n_gene_detect_ge50 = sum(icc_reliable),
    median_icc_pot_detect_ge50 = median(icc_pot[icc_reliable], na.rm = TRUE),
    max_icc_pot_detect_ge50 = max(icc_pot[icc_reliable], na.rm = TRUE),
    n_p005_lm = sum(p_lm_no_pot < 0.05),
    n_p005_lmer = sum(p_lmer_pot < 0.05, na.rm = TRUE),
    n_fdr005_lm = sum(fdr_lm_no_pot < 0.05),
    n_fdr005_lmer = sum(fdr_lmer_pot < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n[11] 每期汇总（含盆 vs 不含盆；ICC 另给 检出率>=50% 子集）：\n")
print(as.data.frame(stage_summary))

status_change <- model_compare %>%
  dplyr::mutate(
    sig_lm = fdr_lm_no_pot < 0.05,
    sig_lmer = fdr_lmer_pot < 0.05,
    change = dplyr::case_when(
      !fit_ok ~ "lmer failed",
      sig_lm & !sig_lmer ~ "lost FDR significance",
      !sig_lm & sig_lmer ~ "gained FDR significance",
      TRUE ~ "unchanged"
    )
  )

cat("\n[11] FDR 显著性状态变化（含盆模型 vs 不含盆模型）：\n")
print(as.data.frame(status_change %>% dplyr::count(stage, change)))

cat("\n[11] 盆效应最大的 10 个基因：\n")
print(as.data.frame(model_compare %>% dplyr::arrange(dplyr::desc(icc_pot)) %>%
                      dplyr::select(stage, gene, icc_pot, p_lm_no_pot, p_lmer_pot,
                                    fdr_lm_no_pot, fdr_lmer_pot, lmer_singular) %>%
                      dplyr::slice_head(n = 10)))

bioRtools::write_data(
  list(model_compare = model_compare, stage_summary = stage_summary,
       pot_structure = pot_structure, status_change = status_change),
  file = path_out
)

cat("\n[11] 输出：", path_out, "\n")

cat("\n[11] sessionInfo()\n")
sessionInfo()
