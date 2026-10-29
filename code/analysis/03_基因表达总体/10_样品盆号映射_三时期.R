# ============================================================
# 03_基因表达总体 / 10_样品盆号映射_三时期.R
# 作用：从三批 Fluidigm 的样品标签表恢复**盆号**，生成"样品 ID → 盆号 → 同盆伙伴样品"映射
#       （审稿意见 M3：同一 Mix 盆的两株共享盆环境，应在模型中体现）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_基因表达总体/10_样品盆号映射_三时期.R"
#
# 输入（项目内只读副本）：
#   <root>/data/02_原始数据/04_Fluidigm表达/01_第一批/第一批96个样品选择表.xlsx
#        sheet `info`：focal / neighbor / label（形如 `2(15_15)` = 盆号(focal_neighbor)）/ tag（= 分析用样品 ID）
#   <root>/data/02_原始数据/04_Fluidigm表达/02_第二批/XL1第二批RNA样品编号.xlsx
#        sheet `Sheet1`：首行表头，其后每 3 行为一组（盆标签行 / `Focal` 样品 ID 行 / `Neighbor` 样品 ID 行）
#   <root>/data/02_原始数据/04_Fluidigm表达/03_第三批/samples.xlsx
#        sheet `samples.tag`：结构同上（每 3 行一组）
#   <root>/result/03_基因表达总体/11_表_表达矩阵_三时期.xlsx（sample_matrix，用于覆盖度核对）
# 输出：<root>/result/03_基因表达总体/19_表_样品盆号映射_三时期.xlsx
#       sheet：map_long（sample → pot / focal / neighbor / 伙伴样品）、pot_summary、coverage
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 上游原始路径（副本来源，未修改；上游标签表 2026-10-18 复制入项目）：
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/田间试验/fluidigm/第二批/XL1第二批RNA样品编号.xlsx
#   /Users/lixiang/NutstoreFiles/01.博士相关/01.data/02.元阳梯田/田间试验/fluidigm/第三批/samples.xlsx
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
dir_flu <- file.path(project_root, "data", "02_原始数据", "04_Fluidigm表达")
path_b1 <- file.path(dir_flu, "01_第一批", "第一批96个样品选择表.xlsx")
path_b2 <- file.path(dir_flu, "02_第二批", "XL1第二批RNA样品编号.xlsx")
path_b3 <- file.path(dir_flu, "03_第三批", "samples.xlsx")
path_mat <- file.path(project_root, "result", "03_基因表达总体", "11_表_表达矩阵_三时期.xlsx")
out_dir <- file.path(project_root, "result", "03_基因表达总体")
path_out <- file.path(out_dir, "19_表_样品盆号映射_三时期.xlsx")

stopifnot(file.exists(path_b1), file.exists(path_b2), file.exists(path_b3), file.exists(path_mat))

read_sheet_raw <- function(path, sheet, col_names = TRUE) {
  readxl::read_excel(path, sheet = sheet, col_names = col_names, .name_repair = "minimal") %>%
    as.data.frame()
}

# 从 `2(15_15)` 形式解析盆号与组合
parse_label <- function(x) {
  x <- as.character(x)
  pot <- suppressWarnings(as.integer(stringr::str_extract(x, "^\\s*[0-9]+")))
  combo <- stringr::str_match(x, "\\(([^)]*)\\)")[, 2]
  part <- stringr::str_split_fixed(combo, "_", 2)
  tibble::tibble(pot = pot, pot_focal = part[, 1], pot_neighbor = part[, 2])
}

# ---- 第一批：一行一个样品 ----
read_sheet_raw(path_b1, "info") %>%
  dplyr::transmute(sample = as.character(tag), focal = as.character(focal),
                   neighbor = as.character(neighbor), label = as.character(label)) -> b1.raw

b1.raw %>% dplyr::bind_cols(parse_label(b1.raw$label)) %>%
  dplyr::mutate(stage = "tillering", role = "focal_plant") -> b1.map

# ---- 第二批 / 第三批：每 3 行为一组（标签行 / Focal / Neighbor）----
parse_layout <- function(path, sheet, stage) {
  raw <- read_sheet_raw(path, sheet, col_names = FALSE)
  # 第 1 行为表头；其后每 3 行为一组（盆标签行 / Focal 行 / Neighbor 行）。
  # 注意：两个批次文件的列偏移不同（第二批 12 列、第三批 13 列），
  # 因此取值列**由标签行里含 "(" 的列动态判定**，不写死列号。
  body <- raw[-1, , drop = FALSE]
  n_group <- floor(nrow(body) / 3)
  out <- vector("list", n_group)
  for (g in seq_len(n_group)) {
    lab_row <- as.character(unlist(body[3 * (g - 1) + 1, , drop = TRUE]))
    value_cols <- which(grepl("\\(", lab_row))
    lab <- lab_row[value_cols]
    foc <- as.character(unlist(body[3 * (g - 1) + 2, value_cols, drop = TRUE]))
    nei <- as.character(unlist(body[3 * (g - 1) + 3, value_cols, drop = TRUE]))
    keep <- !is.na(lab) & nzchar(lab)
    parsed <- parse_label(lab[keep])
    out[[g]] <- tibble::tibble(
      label = lab[keep],
      pot = parsed$pot,
      pot_focal = parsed$pot_focal,
      pot_neighbor = parsed$pot_neighbor,
      sample_focal = foc[keep],
      sample_neighbor = nei[keep],
      stage = stage
    )
  }
  dplyr::bind_rows(out)
}

b2.layout <- parse_layout(path_b2, "Sheet1", "flowering")
b3.layout <- parse_layout(path_b3, "samples.tag", "filling")

cat("[10] 第二批标签行数：", nrow(b2.layout), "；第三批：", nrow(b3.layout), "\n")

bind_layout <- function(d) {
  f <- d %>% dplyr::transmute(stage, pot, label, sample = as.character(sample_focal),
                              focal = pot_focal, neighbor = pot_neighbor,
                              role = "focal_plant")
  n <- d %>% dplyr::transmute(stage, pot, label, sample = as.character(sample_neighbor),
                              focal = pot_neighbor, neighbor = pot_focal,
                              role = "neighbour_plant")
  dplyr::bind_rows(f, n)
}

b23.map <- dplyr::bind_rows(bind_layout(b2.layout), bind_layout(b3.layout)) %>%
  dplyr::filter(!is.na(sample), nzchar(sample))

# ---- 合并三批映射，并标注同盆伙伴样品 ----
dplyr::bind_rows(b1.map, b23.map) %>%
  dplyr::select(stage, pot, sample, focal, neighbor, role, label) %>%
  dplyr::distinct() -> map.long

map.long %>%
  dplyr::group_by(stage, pot) %>%
  dplyr::mutate(n_sample_in_pot = dplyr::n(),
                partner_sample = dplyr::if_else(
                  n_sample_in_pot > 1,
                  vapply(seq_len(dplyr::n()), function(i) {
                    others <- setdiff(sample, sample[i])
                    if (length(others) == 0) NA_character_ else others[1]
                  }, character(1)),
                  NA_character_)) %>%
  dplyr::ungroup() -> map.long

# ---- 覆盖度核对：与三时期样本级矩阵对齐 ----
read_sheet_raw(path_mat, "sample_matrix") %>%
  dplyr::transmute(stage = as.character(stage), sample = as.character(sample),
                   focal, neighbor) -> mat.meta

mat.meta %>%
  dplyr::left_join(map.long, by = c("stage", "sample")) -> coverage_tab

coverage <- coverage_tab %>%
  dplyr::group_by(stage) %>%
  dplyr::summarise(n_sample = dplyr::n(),
                   n_with_pot = sum(!is.na(pot)),
                   pct_with_pot = round(100 * mean(!is.na(pot)), 1),
                   n_unique_pot = dplyr::n_distinct(pot, na.rm = TRUE),
                   .groups = "drop")

# 分析样本里"同一盆出现 ≥2 个样品"的盆数（即真正共享盆环境的 Mix 盆）
pot_in_analysis <- coverage_tab %>%
  dplyr::filter(!is.na(pot)) %>%
  dplyr::count(stage, pot, name = "n_sample_in_analysis")

coverage <- coverage %>%
  dplyr::left_join(
    pot_in_analysis %>%
      dplyr::group_by(stage) %>%
      dplyr::summarise(n_pot_in_analysis = dplyr::n(),
                       n_pot_shared = sum(n_sample_in_analysis >= 2),
                       .groups = "drop"),
    by = "stage"
  )

cat("\n[10] 覆盖度（分析样本中能匹配到盆号的比例）：\n")
print(as.data.frame(coverage))

pot_summary <- map.long %>%
  dplyr::count(stage, pot, n_sample_in_pot, name = "n_rows") %>%
  dplyr::group_by(stage) %>%
  dplyr::summarise(n_pot = dplyr::n(),
                   n_pot_two_samples = sum(n_sample_in_pot == 2),
                   n_pot_one_sample = sum(n_sample_in_pot == 1),
                   .groups = "drop")

cat("\n[10] 盆结构（来自标签表）：\n")
print(as.data.frame(pot_summary))

cat("\n[10] 分析样本中「同盆两株都进入分析」的情况：\n")
print(as.data.frame(coverage_tab %>%
                      dplyr::filter(!is.na(pot), !is.na(partner_sample)) %>%
                      dplyr::distinct(stage, pot) %>%
                      dplyr::count(stage, name = "n_pot_both_plants")))

bioRtools::write_data(
  list(map_long = map.long, pot_summary = pot_summary,
       coverage = coverage, sample_level = coverage_tab),
  file = path_out
)

cat("\n[10] 输出：", path_out, "\n")

cat("\n[10] sessionInfo()\n")
sessionInfo()
