# ============================================================
# 03_基因表达总体 / 08_PLSDA_三时期.R
# 作用：三时期组合均值矩阵的 PLS-DA（ropls::opls）—— 得分、载荷、VIP、交叉验证统计与主轴检验
#       Figure 2B–2E 的三时期数据来源
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_基因表达总体/08_PLSDA_三时期.R"
#
# 输入：<root>/result/03_基因表达总体/11_表_表达矩阵_三时期.xlsx（combo_matrix，48 行 = 3 期 × 16 组合）
# 输出：<root>/result/03_基因表达总体/13_表_PLSDA_三时期.xlsx
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 口径：与脚本 03（两时期）一致（ropls::opls、predI = 3、orthoI = 0、R2X 取 modelDF 单成分比例），
#   仅两处适配三时期：
#   ① 时期之间的比较由"两组 Welch t + Wilcoxon"改为**单因素检验**（p ~ stage 的 one-way ANOVA、
#      Kruskal–Wallis，以及三对两两 Wilcoxon + BH 校正），因为 t.test/wilcox.test 只接受两组；
#   ② 每个时期内的 Pure vs Mix 检验保持不变（Welch t + Wilcoxon + Cohen's d）。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_in <- file.path(project_root, "result", "03_基因表达总体", "11_表_表达矩阵_三时期.xlsx")
out_dir <- file.path(project_root, "result", "03_基因表达总体")
path_out <- file.path(out_dir, "13_表_PLSDA_三时期.xlsx")

stopifnot(file.exists(path_in))

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_in, "combo_matrix") -> mat.combo
bioRtools::read_data(path_in, sheet = "gene_annotation") -> df.gene

meta_cols <- c("combo", "stage", "focal", "neighbor", "group")
gene_cols <- setdiff(colnames(mat.combo), meta_cols)
stopifnot(length(gene_cols) == 92, nrow(mat.combo) == 48)

stage_levels <- c("tillering", "flowering", "filling")
mat.combo %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels)) -> mat.combo

cat("[08] 输入矩阵：", nrow(mat.combo), "行 ×", length(gene_cols), "基因\n")
print(table(mat.combo$stage, mat.combo$group))

mat.combo %>%
  dplyr::select(dplyr::all_of(meta_cols)) %>%
  dplyr::mutate(plotgroup = paste0(stage, " ", group)) -> df.meta

x <- as.matrix(mat.combo[, gene_cols])
rownames(x) <- mat.combo$combo
y <- as.factor(df.meta$group)

# predI = 3、orthoI = 0（与两时期脚本一致；orthoI = NA 会因 CV 拒绝而返回空模型）
ropls::opls(x, y, predI = 3) -> model

cat("\n[08] 模型摘要：\n")
print(model)
stopifnot(length(ropls::getScoreMN(model)) > 0)

as.data.frame(model@summaryDF) %>%
  tibble::rownames_to_column(var = "model_type") -> cv_tab

cat("\n[08] 交叉验证统计：\n")
print(as.data.frame(cv_tab))

ropls::getScoreMN(model) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "combo") %>%
  dplyr::left_join(df.meta, by = "combo") -> scores_tab

ropls::getLoadingMN(model) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>%
  dplyr::select(gene, p1, p2, p3) -> loadings_tab

ropls::getVipVn(model) %>%
  tibble::enframe(name = "gene", value = "vip") %>%
  dplyr::arrange(dplyr::desc(vip)) -> vip_tab

vip_top10 <- dplyr::slice_head(vip_tab, n = 10)
cat("\n[08] VIP 前 10 基因：\n")
print(as.data.frame(vip_top10))

as.data.frame(model@modelDF) %>%
  tibble::rownames_to_column(var = "component") %>%
  dplyr::mutate(R2X_single = R2X, R2Y_single = R2Y) -> variance_tab

cat("\n[08] 方差解释率：\n")
print(as.data.frame(variance_tab))

# ---- 主轴检验：时期整体差异（单因素）+ 各期内 Pure vs Mix ----
score_test <- function(axis_col, label) {
  d <- scores_tab
  d$axis_value <- d[[axis_col]]

  aov_fit <- stats::aov(axis_value ~ stage, data = d)
  aov_p <- summary(aov_fit)[[1]][["Pr(>F)"]][1]
  kw <- stats::kruskal.test(axis_value ~ stage, data = d)

  pair_w <- stats::pairwise.wilcox.test(d$axis_value, d$stage, p.adjust.method = "BH")
  pairs_tab <- as.data.frame(as.table(pair_w$p.value)) %>%
    dplyr::filter(!is.na(Freq)) %>%
    dplyr::transmute(contrast = paste(Var1, "vs", Var2, sep = " "), p_value = Freq)

  by_stage <- d %>%
    dplyr::group_by(stage) %>%
    dplyr::group_modify(~ {
      tt <- stats::t.test(.x$axis_value ~ .x$group)
      ww <- suppressWarnings(stats::wilcox.test(.x$axis_value ~ .x$group))
      dd <- (mean(.x$axis_value[.x$group == "Mix"]) - mean(.x$axis_value[.x$group == "Pure"])) /
        sqrt((stats::var(.x$axis_value[.x$group == "Pure"]) +
                stats::var(.x$axis_value[.x$group == "Mix"])) / 2)
      tibble::tibble(
        n_pure = sum(.x$group == "Pure"), n_mix = sum(.x$group == "Mix"),
        mean_pure = mean(.x$axis_value[.x$group == "Pure"]),
        mean_mix = mean(.x$axis_value[.x$group == "Mix"]),
        diff_mix_minus_pure = mean_mix - mean_pure,
        t_p_value = tt$p.value, wilcoxon_p_value = ww$p.value, cohens_d = dd
      )
    }) %>%
    dplyr::ungroup()

  dplyr::bind_rows(
    tibble::tibble(axis = label, test_group = "stage_overall", comparison = "all three stages",
                   test = "one_way_ANOVA", estimate = NA_real_, p_value = aov_p,
                   statistic = NA_real_, n = nrow(d)),
    tibble::tibble(axis = label, test_group = "stage_overall", comparison = "all three stages",
                   test = "kruskal_wallis", estimate = NA_real_, p_value = kw$p.value,
                   statistic = unname(kw$statistic), n = nrow(d)),
    pairs_tab %>%
      dplyr::mutate(axis = label, test_group = "stage_pairwise", test = "wilcoxon_BH",
                    estimate = NA_real_, statistic = NA_real_, n = nrow(d)) %>%
      dplyr::select(axis, test_group, comparison = contrast, test, estimate, p_value, statistic, n),
    by_stage %>%
      tidyr::pivot_longer(cols = c(t_p_value, wilcoxon_p_value), names_to = "test", values_to = "p_value") %>%
      dplyr::mutate(axis = label, test_group = "within_stage",
                    comparison = paste0(stage, ": Mix vs Pure"),
                    test = dplyr::if_else(test == "t_p_value", "welch_t", "wilcoxon")) %>%
      dplyr::select(axis, test_group, comparison, test, estimate = diff_mix_minus_pure, p_value,
                    cohens_d, n_pure, n_mix)
  )
}

test_tab <- dplyr::bind_rows(score_test("p1", "Dim 1"), score_test("p2", "Dim 2"))

cat("\n[08] 主轴检验：\n")
print(as.data.frame(test_tab))

bioRtools::write_data(
  list(
    scores = scores_tab,
    loadings = loadings_tab,
    vip = vip_tab,
    vip_top10 = vip_top10,
    variance_explained = variance_tab,
    axis_tests = test_tab,
    cv_summary = cv_tab
  ),
  file = path_out
)
cat("\n[08] 输出：", path_out, "\n")

cat("\n[08] sessionInfo()\n")
sessionInfo()
