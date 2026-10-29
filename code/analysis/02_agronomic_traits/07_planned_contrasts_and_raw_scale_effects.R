# ============================================================
# 02_agronomic_traits / 07_planned_contrasts_and_raw_scale_effects.R
# 作用：按审稿意见 M4 统一agronomic_traits的statistics口径并补"原始尺度effect量 + 百分比变化"
#       ① tiller_number / 单株yield / height_change_index：改为**计划对比**（每个 Mix 邻体 vs 该 focal 自身的 Pure），
#          取代原来的"各邻体两两 Tukey"；同时给出原始尺度差值and百分比变化及 95% CI；
#       ② 温室inoculation（稻瘟病 / 胡麻叶斑病）：盆级混合模型给出原始尺度effectand百分比变化。
# 运行：/opt/homebrew/bin/Rscript "<root>/code/02_agronomic_traits/07_planned_contrasts_and_raw_scale_effects.R"
#
# 输入：<root>/data/processed_inputs/03_agronomic_traits_height_tiller_yield.xlsx
#       <root>/result/02_agronomic_traits/01_表_height_change_index.xlsx（index_long）
#       <root>result/01_inoculation_and_field_disease/01_表_inoculationstatistics.xlsx（clean_long，若不存在则回退用analysis输入表）
# 输出：<root>/result/02_agronomic_traits/13_表_planned_contrasts_and_raw_scale_effects.xlsx
#       sheet：planned_contrasts_traits / height_index_contrasts / disease_effects / raw_scale_summary
# 其中 <root> = <repository_root>
#
# 口径notes：
#   - "计划对比"= 每个 Mix 邻体组合 vs 该 focal 自己的单作对照（emmeans + Dunnett correction）；
#     这正是审稿意见要求的comparison对象（原版做的是全部邻体两两 Tukey，属于非计划comparison）。
#   - 原始尺度：tiller_number为"个/株"、yield为"g/株"、plant_height为"cm"，差值即原始单位差；
#     百分比变化 = 差值 / 该 focal 单作均值 × 100%。
#   - inoculationdata用盆级混合模型 `lesion ~ planting mode + (1|batch) + (1|pot)`（and主题 06 一致），
#     百分比变化以单作（Pure）的总体均值为分母。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
path_traits <- file.path(project_root, "data", "processed_inputs", "03_agronomic_traits_height_tiller_yield.xlsx")
path_index <- file.path(project_root, "result", "02_agronomic_traits", "01_表_height_change_index.xlsx")
path_inoc <- file.path(project_root, "data", "processed_inputs", "01_inoculation_merged.xlsx")
out_dir <- file.path(project_root, "result", "02_agronomic_traits")
path_out <- file.path(out_dir, "13_表_planned_contrasts_and_raw_scale_effects.xlsx")

stopifnot(file.exists(path_traits), file.exists(path_index), file.exists(path_inoc))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

focus_levels <- c("15", "49", "69", "87")

# ============================================================
# 1) tiller_number / 单株yield：计划对比（每个 Mix 邻体 vs 自身单作）
# ============================================================
bioRtools::read_data(path_traits) %>%
  dplyr::filter(phenotype %in% c("tiller number", "yield")) %>%
  dplyr::mutate(
    focal = factor(focal, levels = focus_levels),
    neighbor = factor(neighbor, levels = focus_levels),
    group = factor(group, levels = c("Pure", "Mix"))
  ) -> df.traits

cat("[07] agronomic_traits行数：", nrow(df.traits), "\n")

planned_contrast_one <- function(d, pheno, fo) {
  dd <- d %>% dplyr::filter(phenotype == pheno, focal == fo) %>% droplevels()
  if (nrow(dd) < 5) return(tibble::tibble())
  m <- stats::lm(value ~ neighbor, data = dd)
  emm <- emmeans::emmeans(m, ~ neighbor)
  # 直接用对比向量构造"每个 Mix 邻体 - 该 focal 单作"，避免 emmeans::trt.vs.ctrl 对 ref 的解析问题
  lev <- levels(dd$neighbor)
  ref_i <- which(lev == fo)
  others <- setdiff(seq_along(lev), ref_i)
  con_list <- lapply(others, function(j) {
    v <- numeric(length(lev)); v[j] <- 1; v[ref_i] <- -1; v
  })
  names(con_list) <- paste0(lev[others], " - ", fo)
  co <- as.data.frame(summary(
    emmeans::contrast(emm, method = con_list, adjust = "dunnett"),
    infer = c(TRUE, TRUE)
  ))
  # 每个对比对应的 Mix sample量（逐对比计算，避免长度回收错误）
  co$n_mix <- vapply(as.character(co$contrast), function(lbl) {
    nb <- stringr::str_trim(stringr::str_remove(lbl, "\\s*-\\s*[^-]+$"))
    sum(dd$group == "Mix" & as.character(dd$neighbor) == nb)
  }, numeric(1))
  pure_mean <- mean(dd$value[dd$group == "Pure"])
  co %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      phenotype = pheno, focal = fo,
      contrast = as.character(contrast),
      # 标签形如 "49 - 15"（对比 - 参照），取前半段作为邻体
      mix_neighbor = stringr::str_trim(stringr::str_remove(as.character(contrast), "\\s*-\\s*[^-]+$")),
      estimate_raw = estimate,
      se = SE, df = df,
      ci_lower = lower.CL, ci_upper = upper.CL,
      p_value = p.value,
      pure_mean = pure_mean,
      pct_change = 100 * estimate / pure_mean,
      pct_ci_lower = 100 * lower.CL / pure_mean,
      pct_ci_upper = 100 * upper.CL / pure_mean,
      n_pure = sum(dd$group == "Pure"),
      n_mix = n_mix
    )
}

planned_traits <- tidyr::expand_grid(
  pheno = c("tiller number", "yield"),
  fo = focus_levels
) %>%
  purrr::pmap_dfr(function(pheno, fo) planned_contrast_one(df.traits, pheno, fo)) %>%
  dplyr::mutate(p_fdr = stats::p.adjust(p_value, method = "BH"),
                .after = p_value)

cat("\n[07] tiller_number / yield的计划对比（每个 Mix 邻体 vs 该 focal 单作；Dunnett correction）：\n")
print(as.data.frame(planned_traits %>%
                      dplyr::select(phenotype, focal, mix_neighbor, estimate_raw, ci_lower, ci_upper,
                                    p_value, p_fdr, pct_change, pct_ci_lower, pct_ci_upper)))

# ============================================================
# 2) height_change_index：每个组合 vs 0（含 95% CI and相对单作百分比）
# ============================================================
bioRtools::read_data(path_index, sheet = "index_long") %>%
  dplyr::mutate(focal = factor(focal, levels = focus_levels),
                neighbor = factor(neighbor, levels = focus_levels),
                group = factor(group, levels = c("Pure", "Mix"))) -> df.index

# 注意：index_long 里每个 (focal, neighbor) 有"15 次调查 × 多株"的行，直接按行数当 n 会重复计数
# （这正是审稿意见 M3/M4 指出的伪重复）。这里先按"每次调查取一个均值"聚合到 15 个time点，
# 再以这 15 个time点为analysis单元做 vs 0 的test；纵向模型的推断结论见主题 06（Figure S9/S12）。
df.index %>%
  dplyr::group_by(focal, neighbor, group, time) %>%
  dplyr::summarise(index_by_survey = mean(index, na.rm = TRUE),
                   n_plant = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(focal, neighbor, group) %>%
  dplyr::summarise(index_mean = mean(index_by_survey, na.rm = TRUE),
                   index_sd = stats::sd(index_by_survey, na.rm = TRUE),
                   n_time = dplyr::n(), .groups = "drop") -> idx_sum

height_index_contrasts <- idx_sum %>%
  dplyr::filter(group == "Mix") %>%
  dplyr::group_by(focal) %>%
  dplyr::mutate(pure_index = idx_sum$index_mean[idx_sum$focal == dplyr::cur_group()$focal &
                                                 idx_sum$group == "Pure"][1]) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    se = index_sd / sqrt(n_time),
    ci_lower = index_mean - stats::qt(0.975, df = max(n_time - 1, 1)) * se,
    ci_upper = index_mean + stats::qt(0.975, df = max(n_time - 1, 1)) * se,
    t_value = index_mean / se,
    p_value = 2 * stats::pt(-abs(t_value), df = max(n_time - 1, 1))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_fdr = stats::p.adjust(p_value, method = "BH"), .after = p_value) %>%
  dplyr::select(focal, neighbor, index_mean, se, ci_lower, ci_upper, t_value, p_value, p_fdr, n_time)

cat("\n[07] height_change_index（Mix 组合 vs 0，按time平均）：\n")
print(as.data.frame(height_index_contrasts))

# ============================================================
# 3) 温室inoculation：盆级模型的原始尺度effectand百分比变化
# ============================================================
bioRtools::read_data(path_inoc) -> df.inoc

cat("\n[07] inoculationdata列名：", paste(colnames(df.inoc), collapse = ", "), "\n")

disease_effect <- purrr::map_dfr(c("blast", "brown spot"), function(dz) {
  dd <- df.inoc %>% dplyr::filter(disease == dz) %>% droplevels()
  if (nrow(dd) == 0) return(tibble::tibble())
  dd <- dd %>%
    dplyr::mutate(mode = factor(dplyr::if_else(focal == neighbor, "Pure", "Mix"),
                                levels = c("Pure", "Mix")),
                  batch = factor(batch), pot = factor(image))
  m <- lmerTest::lmer(lesion_number_per_cm2 ~ mode + (1 | batch) + (1 | pot), data = dd)
  co <- summary(m)$coefficients
  est <- co["modeMix", "Estimate"]
  se <- co["modeMix", "Std. Error"]
  df_res <- co["modeMix", "df"]
  p <- co["modeMix", "Pr(>|t|)"]
  pure_mean <- mean(dd$lesion_number_per_cm2[dd$mode == "Pure"])
  tibble::tibble(
    disease = dz,
    model = "lmer(lesion ~ mode + (1|batch) + (1|pot))",
    n_leaf = nrow(dd), n_pot = dplyr::n_distinct(dd$pot),
    estimate_raw = est, se = se, df = df_res,
    ci_lower = est - stats::qt(0.975, df_res) * se,
    ci_upper = est + stats::qt(0.975, df_res) * se,
    p_value = p,
    pure_mean = pure_mean,
    pct_change = 100 * est / pure_mean,
    pct_ci_lower = 100 * (est - stats::qt(0.975, df_res) * se) / pure_mean,
    pct_ci_upper = 100 * (est + stats::qt(0.975, df_res) * se) / pure_mean
  )
})

cat("\n[07] inoculation病害的盆级effect（原始尺度 + 百分比变化）：\n")
print(as.data.frame(disease_effect))

raw_scale_summary <- dplyr::bind_rows(
  planned_traits %>%
    dplyr::transmute(trait = paste0(phenotype, " (", focal, " x ", mix_neighbor, ")"),
                     unit = dplyr::if_else(phenotype == "yield", "g per plant", "tillers per plant"),
                     estimate_raw, ci_lower, ci_upper, p_value, p_fdr,
                     pct_change, pct_ci_lower, pct_ci_upper, source = "greenhouse / field phenotyping"),
  height_index_contrasts %>%
    dplyr::transmute(trait = paste0("plant height index (", focal, " x ", neighbor, ")"),
                     unit = "cm", estimate_raw = index_mean, ci_lower, ci_upper,
                     p_value, p_fdr,
                     pct_change = NA_real_, pct_ci_lower = NA_real_, pct_ci_upper = NA_real_,
                     source = "field, averaged over 15 surveys"),
  disease_effect %>%
    dplyr::transmute(trait = paste0("lesion density (", disease, ")"),
                     unit = "lesions per cm2", estimate_raw, ci_lower, ci_upper,
                     p_value, p_fdr = NA_real_,
                     pct_change, pct_ci_lower, pct_ci_upper, source = "greenhouse, pot-level model")
)

bioRtools::write_data(
  list(
    planned_contrasts_traits = planned_traits,
    height_index_contrasts = height_index_contrasts,
    disease_effects = disease_effect,
    raw_scale_summary = raw_scale_summary
  ),
  file = path_out
)

cat("\n[07] 输出：", path_out, "\n")

cat("\n[07] sessionInfo()\n")
sessionInfo()
