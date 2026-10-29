# ============================================================
# 03_overall_gene_expression / 03_PLSDA.R
# 作用：组合均值矩阵的 PLS-DA（ropls::opls）—— 得分、载荷、VIP，以及主轴 1/2 的组间test
#       Figure 2B–2E 的dataand_statistics输入
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_overall_gene_expression/03_PLSDA.R"
#
# 输入：<root>/result/03_overall_gene_expression/01_表_expression矩阵.xlsx
#       （combo_matrix 用 readxl::read_excel(.name_repair = "minimal") 读取以保留gene名）
# 输出：<root>/result/03_overall_gene_expression/03_表_PLSDA.xlsx
# 其中 <root> = <repository_root>
#
# 方法：ropls::opls(x = geneexpression矩阵（32 个组合均值）, y = group（Pure/Mix）, predI = 3)。
#   and旧版一致（旧版亦用 ropls::opls；项目 AGENTS.md 写的 mixOmics and之同源，本次沿用 ropls 以对齐旧图）。
#   方差解释率取模型 @modelDF 的 R2X(cum) 相邻差分（旧版图 2B/2E 轴标签的 26.00% / 8.00% 是硬编码值，本次改为实算）。
# test口径：丢弃旧版 "filter(p1 < 10)" 的临时剔除做法；改为
#   ① 主轴得分在时期之间的comparison（整体）；② 每个时期内 Pure vs Mix 的 Welch t test + Wilcoxon（含effect量 Cohen's d）。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：单步耗时 < 60 s，每次运行都重新generate。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/03_overall_gene_expression/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_in <- file.path(project_root, "result", "03_overall_gene_expression", "01_表_expression矩阵.xlsx")
out_dir <- file.path(project_root, "result", "03_overall_gene_expression")
path_out <- file.path(out_dir, "03_表_PLSDA.xlsx")

# 矩阵 sheet 必须保留gene名原样（read_data 的 universal 修复会把 OsD14-1 改成 OsD14.1）
read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_in, "combo_matrix") -> mat.combo
bioRtools::read_data(path_in, sheet = "gene_annotation") -> df.gene

meta_cols <- c("combo", "stage", "focal", "neighbor", "group")
gene_cols <- setdiff(colnames(mat.combo), meta_cols)
stopifnot(length(gene_cols) == 92, nrow(mat.combo) == 32)

cat("[03] 输入矩阵：", nrow(mat.combo), "行 ×", length(gene_cols), "gene\n")

mat.combo %>%
  dplyr::select(dplyr::all_of(meta_cols)) %>%
  dplyr::mutate(plotgroup = paste0(stage, " ", group)) -> df.meta

x <- as.matrix(mat.combo[, gene_cols])
# 必须设置行名，否则 opls 得分矩阵的行名会退化成 1..32，and组合标签 join 时全部变成 NA
rownames(x) <- mat.combo$combo
y <- as.factor(df.meta$group)

# and旧版一致：predI = 3、默认 orthoI = 0（不加正交成分）。
# 注意：若写 orthoI = NA（自动决定正交成分），ropls 会先做交叉验证，本次data在 CV 下"第一预测成分即不显著"，
#   会直接返回空模型（无 p1/p2 得分），因此不能沿用该写法。
ropls::opls(x, y, predI = 3) -> model

cat("\n[03] 模型摘要：\n")
print(model)

stopifnot(length(ropls::getScoreMN(model)) > 0)

# 交叉验证statistics（R2X / R2Y / Q2 and置换test P 值），随表输出
as.data.frame(model@summaryDF) %>%
  tibble::rownames_to_column(var = "model_type") -> cv_tab

cat("\n[03] 交叉验证statistics：\n")
print(as.data.frame(cv_tab))

# ---- 得分 ----
ropls::getScoreMN(model) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "combo") %>%
  dplyr::left_join(df.meta, by = "combo") -> scores_tab

# ---- 载荷 ----
ropls::getLoadingMN(model) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>%
  dplyr::select(gene, p1, p2, p3) -> loadings_tab

# ---- VIP ----
ropls::getVipVn(model) %>%
  tibble::enframe(name = "gene", value = "vip") %>%
  dplyr::arrange(dplyr::desc(vip)) -> vip_tab

vip_top10 <- dplyr::slice_head(vip_tab, n = 10)
cat("\n[03] VIP 前 10 gene：\n")
print(as.data.frame(vip_top10))

# ---- 方差解释率（R2X 差分）----
model_df <- as.data.frame(model@modelDF) %>%
  tibble::rownames_to_column(var = "component")

model_df %>%
  dplyr::mutate(
    # modelDF 中 R2X / R2Y 本身就是"该成分单独解释的比例"，R2X(cum) / R2Y(cum) 是累计值
    R2X_single = R2X,
    R2Y_single = R2Y
  ) -> variance_tab

cat("\n[03] 方差解释率：\n")
print(as.data.frame(variance_tab))

# ---- 主轴test ----
score_test <- function(axis_col, label) {
  # ① 时期之间的comparison（全部sample）
  t_stage <- stats::t.test(scores_tab[[axis_col]] ~ scores_tab$stage)
  w_stage <- suppressWarnings(stats::wilcox.test(scores_tab[[axis_col]] ~ scores_tab$stage))

  # ② 每个时期内 Pure vs Mix
  by_stage <- scores_tab %>%
    dplyr::group_by(stage) %>%
    dplyr::group_modify(~ {
      tt <- stats::t.test(.x[[axis_col]] ~ .x$group)
      ww <- suppressWarnings(stats::wilcox.test(.x[[axis_col]] ~ .x$group))
      d <- (mean(.x[[axis_col]][.x$group == "Mix"]) - mean(.x[[axis_col]][.x$group == "Pure"])) /
        sqrt((stats::var(.x[[axis_col]][.x$group == "Pure"]) + stats::var(.x[[axis_col]][.x$group == "Mix"])) / 2)
      tibble::tibble(
        n_pure = sum(.x$group == "Pure"),
        n_mix = sum(.x$group == "Mix"),
        mean_pure = mean(.x[[axis_col]][.x$group == "Pure"]),
        mean_mix = mean(.x[[axis_col]][.x$group == "Mix"]),
        diff_mix_minus_pure = mean_mix - mean_pure,
        t_p_value = tt$p.value,
        wilcoxon_p_value = ww$p.value,
        cohens_d = d
      )
    }) %>%
    dplyr::ungroup()

  dplyr::bind_rows(
    tibble::tibble(
      axis = label, comparison = "tillering vs filling", test = "welch_t",
      estimate = unname(diff(t_stage$estimate)), p_value = t_stage$p.value,
      statistic = unname(t_stage$statistic), n = nrow(scores_tab)
    ),
    tibble::tibble(
      axis = label, comparison = "tillering vs filling", test = "wilcoxon",
      estimate = NA_real_, p_value = w_stage$p.value,
      statistic = unname(w_stage$statistic), n = nrow(scores_tab)
    ),
    by_stage %>%
      tidyr::pivot_longer(cols = c(t_p_value, wilcoxon_p_value), names_to = "test", values_to = "p_value") %>%
      dplyr::mutate(
        axis = label,
        comparison = paste0(stage, ": Mix vs Pure"),
        test = dplyr::if_else(test == "t_p_value", "welch_t", "wilcoxon")
      ) %>%
      dplyr::select(axis, comparison, test, estimate = diff_mix_minus_pure, p_value,
                    cohens_d, n_pure, n_mix)
  )
}

test_tab <- dplyr::bind_rows(
  score_test("p1", "Dim 1"),
  score_test("p2", "Dim 2")
)

cat("\n[03] 主轴test（Welch t / Wilcoxon；statistic 列对 t test为 t 值、对 Mix vs Pure 为 Cohen's d）：\n")
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
cat("[03] 输出：", path_out, "\n")

cat("\n[03] sessionInfo()\n")
sessionInfo()
