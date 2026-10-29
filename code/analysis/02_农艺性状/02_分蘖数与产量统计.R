# ============================================================
# 02_农艺性状 / 02_分蘖数与产量统计.R
# 作用：分蘖数与单株产量的各 focal 邻体比较（Figure 1B / 1C 的统计输入）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/02_农艺性状/02_分蘖数与产量统计.R"
#
# 输入：<root>/data/01_分析输入表/03_分蘖数_株高_单株产量.xlsx
# 输出：<root>/result/02_农艺性状/02_表_分蘖数与产量统计.xlsx（多 sheet）
#
# 数据口径（实测核对 2026-09-22）：分蘖数为单时点（20220927）、产量为单时点（20221227）；
#   4 个焦点纯系 × 4 邻体（纯种 20 株、各混种组合 10 株，产量个别组合为 7–9 株）。
#   由于是单时点、无批次/区组等随机结构，主分析用单因素方差分析（Tukey HSD 字母），与旧版一致；
#   另用 `value ~ focal * neighbor` 双因素模型给出 focal、neighbor 及交互项的整体检验（正文"品系自身是主要决定因素"的依据）。
# 统计口径：字母用 Tukey HSD；P 值统一定点写法（P < 0.001 等），避免科学计数法误读。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：单步耗时 < 60 s，每次运行都重新生成。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录见 code/02_农艺性状/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_in <- file.path(project_root, "data", "01_分析输入表", "03_分蘖数_株高_单株产量.xlsx")
out_dir <- file.path(project_root, "result", "02_农艺性状")
path_out <- file.path(out_dir, "02_表_分蘖数与产量统计.xlsx")

stopifnot(file.exists(path_in))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

bioRtools::read_data(path_in) -> d3
stopifnot(nrow(d3) == 8095)

phenos <- c("tiller number", "yield")

d3 %>%
  dplyr::filter(phenotype %in% phenos) %>%
  dplyr::mutate(
    focal = as.factor(focal),
    neighbor = as.factor(neighbor),
    group = factor(group, levels = c("Pure", "Mix"))
  ) -> df

cat("[02] 输入行数：", nrow(df), "；phenotype：", paste(unique(df$phenotype), collapse = ", "), "\n")
print(table(df$phenotype, df$focal, df$group))

# ---- 每个 phenotype × focal：单因素 ANOVA + Tukey 字母 ----
df %>%
  dplyr::group_by(phenotype, focal) %>%
  bioRtools::anova_posthoc(value = "value", group = "neighbor") %>%
  dplyr::ungroup() -> posthoc_tab

cat("\n[02] 各 focal 邻体比较（ANOVA + Tukey 字母）：\n")
print(as.data.frame(posthoc_tab))

# ---- 描述性统计 ----
df %>%
  dplyr::group_by(phenotype, focal, neighbor, group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(value),
    sd = sd(value),
    se = sd(value) / sqrt(dplyr::n()),
    .groups = "drop"
  ) -> desc_tab

# ---- 双因素模型整体检验 ----
fit_two_way <- function(ph) {
  d <- df %>% dplyr::filter(phenotype == ph) %>% droplevels()
  m <- stats::lm(value ~ focal * neighbor, data = d)
  as.data.frame(stats::anova(m)) %>%
    tibble::rownames_to_column(var = "term") %>%
    dplyr::rename(df_num = Df, f_value = `F value`, p_value = `Pr(>F)`) %>%
    dplyr::mutate(phenotype = ph, n_obs = nrow(d)) %>%
    dplyr::select(phenotype, term, df_num, f_value, p_value, n_obs)
}

two_way_tab <- purrr::map_dfr(phenos, fit_two_way)

cat("\n[02] 双因素模型（value ~ focal * neighbor）整体检验：\n")
print(as.data.frame(two_way_tab))

bioRtools::write_data(
  list(
    anova_posthoc = posthoc_tab,
    descriptive = desc_tab,
    two_way_anova = two_way_tab
  ),
  file = path_out
)
cat("[02] 输出：", path_out, "\n")

cat("\n[02] sessionInfo()\n")
sessionInfo()
