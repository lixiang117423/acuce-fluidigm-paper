# ============================================================
# 02_agronomic_traits / 03_phenotype_PCA_and_euclidean_distance.R
# 作用：三性状（height_change_index、tiller_number、单株yield）的表型 PCA and互作影响强度（欧氏距离）
#       —— Figure 1D 的dataand_statistics输入
# 运行：/opt/homebrew/bin/Rscript "<root>/code/02_agronomic_traits/03_phenotype_PCA_and_euclidean_distance.R"
#
# 输入：<root>/result/02_agronomic_traits/01_表_height_change_index.xlsx（sheet index_long）
#       <root>/data/processed_inputs/03_agronomic_traits_height_tiller_yield.xlsx（tiller_number、yield）
# 输出：<root>/result/02_agronomic_traits/03_表_表型PCA.xlsx
#       <root>/result/02_agronomic_traits/04_表_表型欧氏距离.xlsx
#
# 口径notes：仅用 4 个焦点纯系 15 / 49 / 69 / 87 的 16 个 focal × neighbor 组合（含 4 个纯种对照），
#   36 / 74 不在本主题范围（其plant_height无time信息，tiller_numberandyield也没有record）。
#   三性状先各自 Z-score normalized（prcomp 再中心化并 scale），PCA 用 16 × 3 矩阵；
#   欧氏距离在normalized后的性状空间计算：以每个 focal 的纯种点为原点，Mix 组合到纯种点的距离即互作影响强度。
# 断点续传：单步耗时 < 30 s，每次运行都重新generate。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/02_agronomic_traits/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_idx <- file.path(project_root, "result", "02_agronomic_traits", "01_表_height_change_index.xlsx")
path_in <- file.path(project_root, "data", "processed_inputs", "03_agronomic_traits_height_tiller_yield.xlsx")
out_dir <- file.path(project_root, "result", "02_agronomic_traits")
path_pca <- file.path(out_dir, "03_表_表型PCA.xlsx")
path_dist <- file.path(out_dir, "04_表_表型欧氏距离.xlsx")

stopifnot(file.exists(path_idx), file.exists(path_in))

bioRtools::read_data(path_idx, sheet = "index_long") -> df.idx
bioRtools::read_data(path_in) -> d3

focus <- c("15", "49", "69", "87")

# ---- 三性状的组合均值 ----
df.idx %>%
  dplyr::filter(as.character(focal) %in% focus) %>%
  dplyr::group_by(focal = as.character(focal), neighbor = as.character(neighbor)) %>%
  dplyr::summarise(mean_height = mean(index), .groups = "drop") -> tr.height

d3 %>%
  dplyr::filter(phenotype %in% c("tiller number", "yield"), focal %in% focus) %>%
  dplyr::mutate(phenotype = dplyr::recode(phenotype, "tiller number" = "mean_tiller", "yield" = "mean_yield")) %>%
  dplyr::group_by(focal, neighbor, phenotype) %>%
  dplyr::summarise(value = mean(value), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = phenotype, values_from = value) -> tr.other

tr.height %>%
  dplyr::left_join(tr.other, by = c("focal", "neighbor")) %>%
  dplyr::mutate(group = dplyr::if_else(focal == neighbor, "Pure", "Mix")) %>%
  dplyr::arrange(focal, neighbor) -> df.combo

df.combo <- dplyr::mutate(df.combo, combo = paste0(focal, "-", neighbor))

stopifnot(nrow(df.combo) == 16, !anyNA(df.combo))
cat("[03] 组合均值表：", nrow(df.combo), "行\n")
print(as.data.frame(df.combo))

# ---- PCA ----
df.combo %>%
  dplyr::select(mean_height, mean_tiller, mean_yield) %>%
  as.matrix() -> mat.raw

# 保留组合名作为行名，否则 scale()/prcomp() 的 scores 会丢掉组合标识（无法and表型信息对齐）
rownames(mat.raw) <- df.combo$combo

mat.scaled <- scale(mat.raw)

stats::prcomp(mat.scaled, center = TRUE, scale. = FALSE) -> pca

variance_tab <- tibble::tibble(
  component = paste0("PC", seq_along(pca$sdev)),
  eigenvalue = pca$sdev^2,
  variance_explained = (pca$sdev^2) / sum(pca$sdev^2)
)

scores_tab <- as.data.frame(pca$x) %>%
  tibble::rownames_to_column(var = "combo") %>%
  dplyr::left_join(df.combo, by = "combo")

loadings_tab <- as.data.frame(pca$rotation) %>%
  tibble::rownames_to_column(var = "trait")

cat("\n[03] 方差解释率：\n")
print(as.data.frame(variance_tab))
cat("\n[03] 载荷：\n")
print(as.data.frame(loadings_tab))

# ---- 欧氏距离（normalized性状空间，以纯种点为原点）----
score_mat <- mat.scaled
rownames(score_mat) <- paste0(df.combo$focal, "-", df.combo$neighbor)

pur_idx <- which(df.combo$focal == df.combo$neighbor)

pur_map <- tibble::tibble(
  focal = df.combo$focal[pur_idx],
  pure_row = pur_idx
)

df.combo %>%
  dplyr::mutate(row_id = dplyr::row_number()) %>%
  dplyr::left_join(pur_map, by = "focal") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    euclidean_to_pure = sqrt(sum((score_mat[row_id, ] - score_mat[pure_row, ])^2))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(focal, neighbor, group, mean_height, mean_tiller, mean_yield, euclidean_to_pure) -> dist_tab

# 焦点纯系整体受影响强度（混种组合的平均距离）
dist_tab %>%
  dplyr::filter(group == "Mix") %>%
  dplyr::group_by(focal) %>%
  dplyr::summarise(
    mean_distance = mean(euclidean_to_pure),
    max_distance = max(euclidean_to_pure),
    strongest_neighbor = neighbor[which.max(euclidean_to_pure)],
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(mean_distance)) -> focal_summary

cat("\n[03] 各 focal 受互作影响强度（混种到纯种的欧氏距离）：\n")
print(as.data.frame(focal_summary))

bioRtools::write_data(
  list(combo_means = df.combo, variance_explained = variance_tab,
       scores = scores_tab, loadings = loadings_tab),
  file = path_pca
)
cat("[03] 输出：", path_pca, "\n")

bioRtools::write_data(
  list(distance_long = dist_tab, focal_summary = focal_summary),
  file = path_dist
)
cat("[03] 输出：", path_dist, "\n")

cat("\n[03] sessionInfo()\n")
sessionInfo()
