# ============================================================
# 01_inoculation_and_field_disease / 02_inoculation_focal_combination_statistics.R
# 作用：按 病害 × focal 纯系，comparison各邻体组合（combination = focal-neighbor）的病斑密度差异，
#       输出 emmeans 估计值 + 字母标记（Figure S2C–S2F 的statistics输入）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/01_inoculation_and_field_disease/02_inoculation_focal_combination_statistics.R"
#
# 输入：<root>/data/processed_inputs/01_inoculation_merged.xlsx
# 输出：<root>/result/01_inoculation_and_field_disease/03_表_inoculation各focal组合emmeans字母.xlsx（多 sheet）
# 其中 <root> = <repository_root>
#
# 模型：lesion_number_per_cm2 ~ combination + (1 | batch)，逐 病害 × focal 单独拟合
# 多重comparison：emmeans 估计组合均值，multcomp::cld() 出字母标记。
#   注意：emmeans 对单一 factor 的全部两两comparison会自动把 adjust = "tukey" 改成 "sidak"
#   （旧版脚本同样如此，故旧版字母实际是 Sidak correction）；本脚本显式写 "sidak"，避免歧义。
#   两两comparison的 df 为 136–261（不是组均值的 df ≈ 1–2），可放心用于comparison。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：单步耗时 < 30 s，每次运行都重新generate。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/01_inoculation_and_field_disease/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_in <- file.path(project_root, "data", "processed_inputs", "01_inoculation_merged.xlsx")
out_dir <- file.path(project_root, "result", "01_inoculation_and_field_disease")
path_out <- file.path(out_dir, "03_表_inoculation各focal组合emmeans字母.xlsx")

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

  # 两两comparison（link 尺度）同时落盘，便于正文引用具体 P 值
  emmeans::contrast(emm, "pairwise") %>%
    as.data.frame() %>%
    dplyr::mutate(disease = ds, focal = f) %>%
    dplyr::select(disease, focal, contrast, estimate, SE, df, t.ratio, p.value) -> contrast_tab

  # and单作对照（focal-focal）的两两comparison：正文描述"某邻体是否显著改变抗性"应以这张表为准
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
    # and对照comparison时注意方向：estimate < 0 表示混种后病斑数低于单作（抗性提高）
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
