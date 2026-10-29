# ============================================================
# 02_agronomic_traits / 04_genetic_distance_and_phenotype_correlation.R
# 作用：遗传距离and三性状（height_change_index、tiller_number、单株yield）的相关性 —— Figure 1F（阴性results）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/02_agronomic_traits/04_genetic_distance_and_phenotype_correlation.R"
#
# 输入：<root>/data/02_raw_data/05_遗传距离/月亮谷纯系的遗传距离.xlsx（99 × 100 距离矩阵）
#       <root>/result/02_agronomic_traits/03_表_表型PCA.xlsx（sheet combo_means）
# 输出：<root>/result/02_agronomic_traits/05_表_genetic_distance_and_phenotype_correlation.xlsx
#
# 口径notes：
#   1) 从 99 个纯系的距离矩阵中取出 4 个焦点纯系（A15 / A49 / A69 / A87）构成 4 × 4 子矩阵；
#   2) 主analysis只用 12 个混种组合（Mix），即"互作effect是否由亲缘关系决定"这一问题；
#      同表另存含 4 个纯种组合（距离 = 0）的 16 组合version，便于and旧版 Figure 1F 对照（旧版含纯种点）；
#   3) 相关analysis用 Pearson（lm 斜率 + r + P），并对 3 个性状做 Bonferroni correction（α = 0.05 / 3）。
# 断点续传：单步耗时 < 30 s，每次运行都重新generate。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/02_agronomic_traits/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_gd <- file.path(project_root, "data", "02_raw_data", "05_遗传距离", "月亮谷纯系的遗传距离.xlsx")
path_combo <- file.path(project_root, "result", "02_agronomic_traits", "03_表_表型PCA.xlsx")
out_dir <- file.path(project_root, "result", "02_agronomic_traits")
path_out <- file.path(out_dir, "05_表_genetic_distance_and_phenotype_correlation.xlsx")

stopifnot(file.exists(path_gd), file.exists(path_combo))

bioRtools::read_data(path_gd) -> gd
bioRtools::read_data(path_combo, sheet = "combo_means") -> df.combo

focus <- c("15", "49", "69", "87")

# ---- 取 4 个焦点纯系的 4 × 4 子矩阵 ----
gd %>%
  dplyr::select(id, dplyr::all_of(paste0("A", focus))) %>%
  dplyr::filter(id %in% paste0("A", focus)) %>%
  tibble::column_to_rownames(var = "id") -> gd.sub

stopifnot(nrow(gd.sub) == 4, ncol(gd.sub) == 4)
cat("[04] 遗传距离子矩阵：\n")
print(gd.sub)

gd.sub %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "focal") %>%
  tidyr::pivot_longer(cols = -focal, names_to = "neighbor", values_to = "genetic_distance") %>%
  dplyr::mutate(
    focal = stringr::str_remove(focal, "^A"),
    neighbor = stringr::str_remove(neighbor, "^A")
  ) -> gd.long

# ---- and性状组合均值合并 ----
df.combo %>%
  dplyr::left_join(gd.long, by = c("focal", "neighbor")) -> df.join

stopifnot(!anyNA(df.join$genetic_distance))

# ---- 相关analysis：Mix only（主）and含纯种（对照）----
trait_cols <- c("mean_height", "mean_tiller", "mean_yield")

cor_one <- function(d, label) {
  purrr::map_dfr(trait_cols, function(tr) {
    m <- stats::lm(d[[tr]] ~ d$genetic_distance)
    s <- summary(m)
    tibble::tibble(
      dataset = label,
      trait = tr,
      n = nrow(d),
      slope = unname(stats::coef(m)[2]),
      r = unname(sqrt(s$r.squared)) * sign(unname(stats::coef(m)[2])),
      r_squared = s$r.squared,
      p_value = s$coefficients[2, 4]
    )
  })
}

cor_mix <- cor_one(dplyr::filter(df.join, group == "Mix"), "Mix only (n = 12)")
cor_all <- cor_one(df.join, "All combos (n = 16, 含纯种)")

cor_tab <- dplyr::bind_rows(cor_mix, cor_all) %>%
  dplyr::group_by(dataset) %>%
  dplyr::mutate(p_bonferroni = p.adjust(p_value, method = "bonferroni")) %>%
  dplyr::ungroup()

cat("\n[04] 遗传距离and三性状的相关（主analysis = Mix only）：\n")
print(as.data.frame(dplyr::filter(cor_tab, dataset == "Mix only (n = 12)")))
cat("\n[04] 对照（含纯种组合）：\n")
print(as.data.frame(dplyr::filter(cor_tab, dataset == "All combos (n = 16, 含纯种)")))

bioRtools::write_data(
  list(
    genetic_distance_submatrix = gd.sub %>% as.data.frame() %>% tibble::rownames_to_column(var = "id"),
    data_joined = df.join,
    correlation = cor_tab
  ),
  file = path_out
)
cat("[04] 输出：", path_out, "\n")

cat("\n[04] sessionInfo()\n")
sessionInfo()
