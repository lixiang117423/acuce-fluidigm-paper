# ============================================================
# 01_接种与田间病害 / 02_接种_各focal组合统计.R
# 作用：按 病害 × focal 纯系，比较各邻体组合（combination = focal-neighbor）的病斑密度差异，
#       输出 emmeans 估计值 + 字母标记（Figure S2C–S2F 的统计输入）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/01_接种与田间病害/02_接种_各focal组合统计.R"
#
# 输入：<root>/data/01_分析输入表/01_接种结果_合并表.xlsx
# 输出：<root>/result/01_接种与田间病害/03_表_接种各focal组合emmeans字母.xlsx（多 sheet）
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 模型：lesion_number_per_cm2 ~ combination + (1 | batch)，逐 病害 × focal 单独拟合
# 多重比较：emmeans 估计组合均值，multcomp::cld() 出字母标记。
#   注意：emmeans 对单一 factor 的全部两两比较会自动把 adjust = "tukey" 改成 "sidak"
#   （旧版脚本同样如此，故旧版字母实际是 Sidak 校正）；本脚本显式写 "sidak"，避免歧义。
#   两两比较的 df 为 136–261（不是组均值的 df ≈ 1–2），可放心用于比较。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：单步耗时 < 30 s，每次运行都重新生成。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录见 code/01_接种与田间病害/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_in <- file.path(project_root, "data", "01_分析输入表", "01_接种结果_合并表.xlsx")
out_dir <- file.path(project_root, "result", "01_接种与田间病害")
path_out <- file.path(out_dir, "03_表_接种各focal组合emmeans字母.xlsx")

stopifnot(file.exists(path_in))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("[02] 输入：", path_in, "\n")

bioRtools::read_data(path_in) -> df.raw

stopifnot(nrow(df.raw) == 1659, ncol(df.raw) == 13)

df.raw %>%
  dplyr::mutate(
    batch = as.factor(batch),
    group = factor(group, levels = c("Pure", "Mix")),
    disease = factor(disease, levels = c("blast", "brown spot"))
  ) -> df

diseases <- levels(df$disease)
focals <- sort(unique(as.character(df$focal)))

grid <- tidyr::expand_grid(disease = diseases, focal = focals)

cat("[02] 待拟合模型数：", nrow(grid), "（", length(diseases), " 病害 x ", length(focals), " focal）\n")

fit_one <- function(ds, f) {
  d <- df %>%
    dplyr::filter(disease == ds, focal == f) %>%
    dplyr::mutate(combination = factor(paste0(focal, "-", neighbor))) %>%
    droplevels()

  m <- lmerTest::lmer(lesion_number_per_cm2 ~ combination + (1 | batch), data = d)
  emm <- emmeans::emmeans(m, "combination")

  emm %>%
    multcomp::cld(alpha = 0.05, Letters = letters, adjust = "sidak", decreasing = TRUE) %>%
    as.data.frame() %>%
    dplyr::mutate(disease = ds, focal = f) %>%
    dplyr::select(disease, focal, combination, emmean, SE, df, lower.CL, upper.CL, signif = .group) -> letters_tab

  # 两两比较（link 尺度）同时落盘，便于正文引用具体 P 值
  emmeans::contrast(emm, "pairwise") %>%
    as.data.frame() %>%
    dplyr::mutate(disease = ds, focal = f) %>%
    dplyr::select(disease, focal, contrast, estimate, SE, df, t.ratio, p.value) -> contrast_tab

  # 与单作对照（focal-focal）的两两比较：正文描述"某邻体是否显著改变抗性"应以这张表为准
  emmeans::contrast(emm, "trt.vs.ctrl",
                    ref = which(levels(d$combination) == paste0(f, "-", f)),
                    adjust = "sidak") %>%
    as.data.frame() %>%
    dplyr::mutate(
      disease = ds,
      focal = f,
      neighbor = stringr::str_remove(contrast, " - .*$"),
      neighbor = stringr::str_remove_all(neighbor, "[()]"),
      neighbor = stringr::str_remove(neighbor, "^[0-9]+-"),
      control = paste0(f, "-", f)
    ) %>%
    dplyr::select(disease, focal, control, neighbor, contrast, estimate, SE, df, t.ratio, p.value) -> ctrl_tab

  list(letters = letters_tab, contrasts = contrast_tab, vs_ctrl = ctrl_tab)
}

res <- purrr::pmap(list(grid$disease, grid$focal), fit_one)

letters_all <- purrr::map_dfr(res, "letters") %>%
  dplyr::mutate(neighbor = stringr::str_remove(combination, "^[0-9]+-")) %>%
  dplyr::select(disease, focal, combination, neighbor, emmean, SE, df, lower.CL, upper.CL, signif)

contrasts_all <- purrr::map_dfr(res, "contrasts")

vs_ctrl_all <- purrr::map_dfr(res, "vs_ctrl") %>%
  dplyr::mutate(
    # 与对照比较时注意方向：estimate < 0 表示混种后病斑数低于单作（抗性提高）
    direction = dplyr::if_else(estimate < 0, "Mix < Pure（病斑更少）", "Mix > Pure（病斑更多）"),
    significant_sidak = dplyr::if_else(p.value < 0.05, "yes", "no")
  ) %>%
  dplyr::select(disease, focal, control, neighbor, estimate, SE, df, t.ratio,
                p.value, significant_sidak, direction, contrast)

cat("\n[02] 字母标记表：", nrow(letters_all), "行\n")
print(as.data.frame(letters_all))

bioRtools::write_data(
  list(emmean_letters = letters_all,
       pure_vs_mix_contrasts = vs_ctrl_all,
       pairwise_contrasts = contrasts_all),
  file = path_out
)
cat("[02] 输出：", path_out, "\n")

cat("\n[02] sessionInfo()\n")
sessionInfo()
