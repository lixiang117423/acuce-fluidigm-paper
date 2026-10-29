# ============================================================
# 06_设计分解与稳健性 / 16_第二轮审稿_正确实验单元重算.R
# 作用：按第二轮审稿意见重做三个核心统计口径：
#   1) 接种试验的 focal / neighbour / focal x neighbour 分解改为 pot-level；
#   2) 田间病害与农艺性状改为 field x year x focal-neighbour plot-level；
#   3) 株高主分析继续以 combination 为纵向单元，并给出“盆/株随机斜率无法可靠识别”的
#      数据证据，避免用不可靠的 pseudoreplication 修补掩盖原始记录缺口。
# 输入：
#   <root>/data/01_分析输入表/01_接种结果_合并表.xlsx
#   <root>/data/01_分析输入表/02_元阳梯田田间病害调查表.xlsx
#   <root>/data/01_分析输入表/03_分蘖数_株高_单株产量.xlsx
# 输出：
#   <root>/result/06_设计分解与稳健性/16_表_第二轮审稿_实验单元重算.xlsx
#   <root>/result/06_设计分解与稳健性/16_Figure5_实验单元重算.pdf / .png
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 关键口径：
#   - 接种：每个 image = 一盆，盆内有 6-14 株、每株一片叶；主检验先用盆内均值;
#     `lesion_pot_mean ~ focal * neighbor + batch`，batch 为固定效应（仅 2-3 个水平）。
#   - 田间：原始表每行 = 一株评分；Field = 3 个田块，Year = 1999/2020/2021。
#     原始表没有 plot 内重复号，故把 Year x Field x focal x neighbour x condition 作为
#     可观测的最小独立单元；同一个单元内的多株先取均值（病害按 10 株计算二项计数），
#     绝不再把单株当独立重复。
#   - 田间病害：quasibinomial GLM，`cbind(diseased, total-diseased) ~ condition * Year + Field`。
#   - 田间性状：`value ~ condition * Year + Field`，并给 Field 随机截距敏感性。
#   - 株高：沿用组合作为纵向单元，报告 group x day；同时输出对“每盆、每株”随机斜率的
#     可识别性检查；原始表缺少稳定单株 ID，不能伪造 `(day|plant)`。
#
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录见 code/06_设计分解与稳健性/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_inoc <- file.path(project_root, "data", "01_分析输入表", "01_接种结果_合并表.xlsx")
path_field <- file.path(project_root, "data", "01_分析输入表", "02_元阳梯田田间病害调查表.xlsx")
path_trait <- file.path(project_root, "data", "01_分析输入表", "03_分蘖数_株高_单株产量.xlsx")
path_height_raw <- file.path(project_root, "data", "02_原始数据", "03_农艺性状", "株高调查结果.xlsx")
out_dir <- file.path(project_root, "result", "06_设计分解与稳健性")
path_out <- file.path(out_dir, "16_表_第二轮审稿_实验单元重算.xlsx")
path_fig <- file.path(out_dir, "16_Figure5_实验单元重算")

stopifnot(file.exists(path_inoc), file.exists(path_field), file.exists(path_trait),
          file.exists(path_height_raw))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

format_p <- function(p, digits = 3) {
  dplyr::if_else(
    is.na(p), "NA",
    dplyr::if_else(p < 0.001, "P < 0.001",
                   paste0("P = ", formatC(p, format = "f", digits = digits)))
  )
}

# ============================================================
# 一、接种试验：pot-level 结构分解
# ============================================================

bioRtools::read_data(path_inoc) -> inoc_raw

inoc_raw %>%
  dplyr::mutate(
    batch = factor(batch),
    group = factor(group, levels = c("Pure", "Mix")),
    focal = factor(focal),
    neighbor = factor(neighbor),
    image = factor(image),
    disease = factor(disease)
  ) -> inoc

# 每盆一个均值：这是主分析单元，而不是每片叶。
inoc %>%
  dplyr::group_by(disease, batch, image, group, focal, neighbor) %>%
  dplyr::summarise(
    lesion_pot_mean = mean(lesion_number_per_cm2, na.rm = TRUE),
    n_leaf = dplyr::n(),
    .groups = "drop"
  ) -> inoc_pot

cat("[16] 接种：叶片 ", nrow(inoc), "，盆 ", dplyr::n_distinct(inoc$image), "\n", sep = "")

inoc_model_tab <- purrr::map_dfr(levels(inoc$disease), function(ds) {
  d <- inoc_pot %>% dplyr::filter(disease == ds) %>% droplevels()

  m0 <- stats::lm(lesion_pot_mean ~ group + batch, data = d)
  m1 <- stats::lm(lesion_pot_mean ~ group * batch, data = d)

  d_mix <- d %>% dplyr::filter(group == "Mix") %>% droplevels()

  # 不完整 4 x 4 设计（4 个自配对单元不存在），car::Anova(type = 2) 会产生别名项。
  # 为避免把无法估计的项误写成显著，这里用显式嵌套模型比较：
  #   focal:        y ~ focal + neighbor + batch  vs  y ~ neighbor + batch
  #   neighbor:     y ~ focal + neighbor + batch  vs  y ~ focal + batch
  #   interaction:  y ~ focal * neighbor + batch  vs  y ~ focal + neighbor + batch
  # 对每个对比取 F 值、分子/分母自由度和 P 值。
  m_base <- stats::lm(lesion_pot_mean ~ batch, data = d_mix)
  m_focal_only <- stats::lm(lesion_pot_mean ~ focal + batch, data = d_mix)
  m_neighbor_only <- stats::lm(lesion_pot_mean ~ neighbor + batch, data = d_mix)
  m_add <- stats::lm(lesion_pot_mean ~ focal + neighbor + batch, data = d_mix)
  m_full <- stats::lm(lesion_pot_mean ~ focal * neighbor + batch, data = d_mix)

  cmp_focal <- stats::anova(m_neighbor_only, m_add, test = "F")
  cmp_neighbor <- stats::anova(m_focal_only, m_add, test = "F")
  cmp_inter <- stats::anova(m_add, m_full, test = "F")

  ss_res <- sum(stats::residuals(m_full)^2)
  df_res <- stats::df.residual(m_full)
  ss_focal <- cmp_focal$`Sum of Sq`[2]
  ss_neighbor <- cmp_neighbor$`Sum of Sq`[2]
  ss_inter <- cmp_inter$`Sum of Sq`[2]

  tibble::tibble(
    disease = ds,
    n_leaf = nrow(inoc %>% dplyr::filter(disease == ds)),
    n_pot = nrow(d),
    n_pot_mix = nrow(d_mix),
    model = "lesion_pot_mean ~ focal * neighbor + batch",
    focal_F = cmp_focal$`F`[2],
    focal_df1 = cmp_focal$Df[2],
    focal_p = cmp_focal$`Pr(>F)`[2],
    neighbor_F = cmp_neighbor$`F`[2],
    neighbor_df1 = cmp_neighbor$Df[2],
    neighbor_p = cmp_neighbor$`Pr(>F)`[2],
    interaction_F = cmp_inter$`F`[2],
    interaction_df1 = cmp_inter$Df[2],
    interaction_df2 = cmp_inter$Res.Df[2],
    interaction_p = cmp_inter$`Pr(>F)`[2],
    partial_eta2_focal = ss_focal / (ss_focal + ss_res),
    partial_eta2_neighbor = ss_neighbor / (ss_neighbor + ss_res),
    partial_eta2_interaction = ss_inter / (ss_inter + ss_res),
    batch_test = "fixed effect",
    note = "Incomplete 4x4 mixture design; interaction df = 5; pot means used as observations"
  )
})

cat("\n[16] 接种 pot-level 结构分解：\n")
print(as.data.frame(inoc_model_tab), row.names = FALSE)

# 保留总体 Pure vs Mix 的 pot-level 估计，便于在正文和图注里统一。
inoc_overall <- purrr::map_dfr(levels(inoc$disease), function(ds) {
  d <- inoc_pot %>% dplyr::filter(disease == ds) %>% droplevels()
  m <- stats::lm(lesion_pot_mean ~ group + batch, data = d)
  co <- summary(m)$coefficients
  est <- co["groupMix", "Estimate"]
  se <- co["groupMix", "Std. Error"]
  tibble::tibble(
    disease = ds,
    estimate = est,
    se = se,
    ci_lower = est - 1.96 * se,
    ci_upper = est + 1.96 * se,
    p_value = co["groupMix", "Pr(>|t|)"],
    n_pot = nrow(d),
    n_leaf = nrow(inoc %>% dplyr::filter(disease == ds))
  )
})

# ============================================================
# 二、田间：聚合到 Year x Field x focal-neighbour 的最小独立单元
# ============================================================

bioRtools::read_data(path_field) -> field_raw

focus_lines <- c("15", "49", "69", "87")

field_raw %>%
  dplyr::select(-1) %>%
  dplyr::mutate(
    FocalPlant = as.character(FocalPlant),
    Neighbor = as.character(Neighbor),
    condition = factor(stringr::str_to_title(condition), levels = c("Pure", "Mix")),
    Field = factor(Field),
    Year = factor(Year),
    Category = as.character(Category)
  ) %>%
  dplyr::filter(FocalPlant %in% focus_lines, !is.na(Value)) %>%
  dplyr::mutate(combo = paste0(FocalPlant, "-", Neighbor)) -> field_focus

# 原始记录是每个 unit 内 10 株/3 行? 只依据现有列，严禁编造 plot 内重复号。
# 这里的最小可观测单元 = Year x Field x combo x condition；同一单元的所有 Value 取均值。
field_plot <- field_focus %>%
  dplyr::group_by(Category, Year, Field, FocalPlant, Neighbor, condition, combo) %>%
  dplyr::summarise(
    value = mean(Value, na.rm = TRUE),
    diseased = sum(Value, na.rm = TRUE),
    total = 10 * dplyr::n(),
    n_scored = dplyr::n(),
    .groups = "drop"
  )

disease_cats <- c("Blast", "BrownSpot", "Smut")
trait_cats <- c("PH", "KGW")

field_disease_tab <- purrr::map_dfr(disease_cats, function(cat) {
  d <- field_plot %>% dplyr::filter(Category == cat) %>% droplevels()
  d$prop <- d$diseased / d$total
  m0 <- stats::glm(cbind(diseased, total - diseased) ~ condition + Year + Field,
                   family = quasibinomial, data = d)
  m1 <- stats::glm(cbind(diseased, total - diseased) ~ condition * Year + Field,
                   family = quasibinomial, data = d)
  a <- stats::anova(m0, m1, test = "F")
  co <- summary(m1)$coefficients

  # marginal / sequential tables for both Year and interaction (not overclaiming causality)
  a_seq <- stats::anova(m1, test = "F")
  tibble::tibble(
    category = cat,
    n_plot = nrow(d),
    model = "quasibinomial: diseased/total ~ condition * Year + Field",
  year_F = a_seq["Year", "F"],
  year_p = a_seq["Year", "Pr(>F)"],
    interaction_F = a$F[2],
    interaction_df1 = a$Df[2],
    interaction_df2 = a$Res.Df[2],
    interaction_p = a$`Pr(>F)`[2],
    condition_baseline_OR = exp(co["conditionMix", "Estimate"]),
    condition_baseline_p = co["conditionMix", "Pr(>|t|)"],
    note = "Field and Year are fixed; plot-level aggregates are the unit"
  )
})

field_trait_tab <- purrr::map_dfr(trait_cats, function(cat) {
  d <- field_plot %>% dplyr::filter(Category == cat) %>% droplevels()
  m0 <- stats::lm(value ~ condition + Year + Field, data = d)
  m1 <- stats::lm(value ~ condition * Year + Field, data = d)
  a <- stats::anova(m0, m1)
  co <- summary(m1)$coefficients
  # Random-intercept sensitivity for Field
  m2 <- lmerTest::lmer(value ~ condition * Year + (1 | Field), data = d)
  a2 <- anova(m2, type = 3)
  tibble::tibble(
    category = cat,
    n_plot = nrow(d),
    model = "Gaussian: value ~ condition * Year + Field",
    interaction_F = a$F[2],
    interaction_df1 = a$Df[2],
    interaction_df2 = a$Res.Df[2],
    interaction_p = a$`Pr(>F)`[2],
    condition_baseline_estimate = co["conditionMix", "Estimate"],
    condition_baseline_p = co["conditionMix", "Pr(>|t|)"],
    random_field_interaction_F = a2["condition:Year", "F value"],
    random_field_interaction_p = a2["condition:Year", "Pr(>F)"],
    random_field_singular = lme4::isSingular(m2),
    note = "Field fixed in primary model; Field random-intercept sensitivity also reported"
  )
})

field_interaction_tab <- dplyr::bind_rows(field_disease_tab, field_trait_tab)

# 逐年效应：每个 Year 内以 plot 为单位，疾病给 OR，性状给均值差。
field_by_year <- purrr::map_dfr(c(disease_cats, trait_cats), function(cat) {
  yrs <- field_plot %>% dplyr::filter(Category == cat) %>% dplyr::pull(Year) %>% droplevels() %>% levels()
  purrr::map_dfr(yrs, function(yr) {
    d <- field_plot %>% dplyr::filter(Category == cat, Year == yr) %>% droplevels()
    if (nrow(d) == 0 || dplyr::n_distinct(d$condition) < 2) return(tibble::tibble())
    if (cat %in% disease_cats) {
      m <- stats::glm(cbind(diseased, total - diseased) ~ condition + Field,
                      family = quasibinomial, data = d)
      co <- summary(m)$coefficients
      est <- co["conditionMix", "Estimate"]
      se <- co["conditionMix", "Std. Error"]
      tibble::tibble(
        category = cat, year = yr, scale = "OR",
        estimate = exp(est), ci_lower = exp(est - 1.96 * se),
        ci_upper = exp(est + 1.96 * se),
        p_value = co["conditionMix", "Pr(>|t|)"], n_plot = nrow(d)
      )
    } else {
      m <- stats::lm(value ~ condition + Field, data = d)
      co <- summary(m)$coefficients
      est <- co["conditionMix", "Estimate"]
      se <- co["conditionMix", "Std. Error"]
      tibble::tibble(
        category = cat, year = yr, scale = "mean difference",
        estimate = est, ci_lower = est - 1.96 * se, ci_upper = est + 1.96 * se,
        p_value = co["conditionMix", "Pr(>|t|)"], n_plot = nrow(d)
      )
    }
  })
})

# ============================================================
# 三、株高：恢复原始盆号，pot-level 纵向模型 + 首次调查基线
# ============================================================
#
# 重要数据事实（由原始表核对，不靠猜测）：
#   - `株高调查结果.xlsx` 的 greenhouse 子集有 position / label，label 形如 "94(87_87)"，
#     可提取物理盆号 1-100；合并后的分析输入表把这两个字段丢掉了。
#   - Pure 盆的两个植株在原始表中按不同 focal 方向重复出现，同一盆、同一天的两个高度值
#     有少量日期只剩 1 个值（可能自然衰老/漏测）；这里先按盆取均值，作为“盆级纵向观测”，
#     不把同一个体重复当独立样本。
#   - Mix 盆：同一物理盆有两个 focal 方向；每个方向取均值，并以 `pot_orientation`
#     （物理盆 × focal）作为纵向随机斜率单元；同时保留物理盆随机截距，检验方向相关性。
#   - 首次调查（day = 0）为 baseline 协变量，因为审稿人指出移栽后第一次调查混种已高 1.8 cm。
#
# 主模型：`value ~ group * day100 + baseline_c + (1 + day100 | pot)`（物理盆随机斜率；
# day100 = day / 100，避免随机斜率数值病态）。首次调查基线按盆方向中心化。
# 敏感性：`(1 + day100 | pot_orientation)`；混合盆两方向共享盆物理环境，需要一并报告。
# 结论只按能稳定收敛的主模型报告；奇异拟合明确标注，不能把奇异模型当确证。

readxl::read_excel(path_height_raw, .name_repair = "minimal") -> height_raw

height_raw %>%
  dplyr::filter(
    `group.2` == "greenhouse",
    !is.na(label), !is.na(time), !is.na(value)
  ) %>%
  dplyr::mutate(
    pot_num = as.integer(stringr::str_extract(label, "^[0-9]+")),
    focal = as.factor(focal),
    neighbor = as.factor(neighbor),
    group = factor(group, levels = c("Pure", "Mix")),
    date = as.Date(as.character(time), format = "%Y%m%d"),
    day = as.numeric(date - min(date)),
    combo = paste(focal, neighbor, sep = "-"),
    pot_orientation = paste(pot_num, focal, sep = "_")
  ) -> height_raw_clean

# Pure 组：原始合并表把同一 Pure 盆的记录复制到了不同 neighbour 标签下；
# 先去掉只由 neighbour 标签造成的重复，再按物理盆/时间取均值。
height_raw_clean %>%
  dplyr::filter(group == "Pure") %>%
  dplyr::distinct(pot_num, focal, time, value, .keep_all = TRUE) %>%
  dplyr::group_by(pot_num, focal, time, day, date) %>%
  dplyr::summarise(value = mean(value), n_plants = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(
    group = factor("Pure", levels = c("Pure", "Mix")),
    neighbor = focal,
    combo = paste(focal, focal, sep = "-"),
    pot_orientation = paste(pot_num, focal, sep = "_")
  ) -> height_pure

height_raw_clean %>%
  dplyr::filter(group == "Mix") %>%
  dplyr::group_by(pot_num, focal, neighbor, time, day, date) %>%
  dplyr::summarise(value = mean(value), n_plants = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(
    group = factor("Mix", levels = c("Pure", "Mix")),
    combo = paste(focal, neighbor, sep = "-"),
    pot_orientation = paste(pot_num, focal, sep = "_")
  ) -> height_mix

height_df <- dplyr::bind_rows(height_pure, height_mix) %>%
  dplyr::mutate(
    pot = factor(pot_num),
    pot_orientation = factor(pot_orientation)
  )

height_baseline <- height_df %>%
  dplyr::filter(day == min(day)) %>%
  dplyr::group_by(pot_orientation) %>%
  dplyr::summarise(baseline = mean(value), .groups = "drop")

height_df <- height_df %>%
  dplyr::left_join(height_baseline, by = "pot_orientation") %>%
  dplyr::filter(!is.na(baseline)) %>%
  dplyr::mutate(
    day100 = day / 100,
    baseline_c = baseline - mean(baseline)
  )

fit_height_model <- function(random_formula, label) {
  f <- stats::as.formula(
    paste("value ~ group * day100 + baseline_c +", random_formula)
  )
  warnings <- character()
  m <- withCallingHandlers(
    lmerTest::lmer(f, data = height_df,
                   control = lme4::lmerControl(
                     optimizer = "bobyqa",
                     optCtrl = list(maxfun = 200000)
                   )),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  co <- summary(m)$coefficients
  tibble::tibble(
    model = label,
    random_effect = random_formula,
    term = rownames(co),
    estimate = co[, "Estimate"],
    se = co[, "Std. Error"],
    df = co[, "df"],
    t_value = co[, "t value"],
    p_value = co[, "Pr(>|t|)"],
    n_obs = nrow(height_df),
    n_pot = dplyr::n_distinct(height_df$pot),
    n_pot_orientation = dplyr::n_distinct(height_df$pot_orientation),
    singular = lme4::isSingular(m),
    convergence_warning = paste(unique(warnings), collapse = " | ")
  )
}

height_rate <- dplyr::bind_rows(
  fit_height_model("(1 + day100 | pot)", "primary: physical pot random slope"),
  fit_height_model("(1 + day100 | pot_orientation)", "sensitivity: pot x focal random slope"),
  fit_height_model("(1 + day100 | pot) + (1 | pot_orientation)", "sensitivity: pot slope + orientation intercept")
)

height_identifiability <- tibble::tibble(
  table = "02_原始数据/03_农艺性状/株高调查结果.xlsx",
  subset = "greenhouse",
  pots_recovered = dplyr::n_distinct(height_df$pot),
  pot_orientations = dplyr::n_distinct(height_df$pot_orientation),
  times = dplyr::n_distinct(height_df$day),
  pure_records_after_dedup = nrow(height_pure),
  mix_records_after_pot_mean = nrow(height_mix),
  baseline_available = all(!is.na(height_df$baseline)),
  n_rows_used = nrow(height_df),
  pots_used = dplyr::n_distinct(height_df$pot),
  pot_orientations_used = dplyr::n_distinct(height_df$pot_orientation),
  note = "Pot number recovered from label prefix; Pure rows deduplicated across neighbour labels before pot-time averaging; rows without day-0 baseline excluded"
)

# ============================================================
# 四、汇总输出
# ============================================================

out_list <- list(
  inoculation_pot_decomposition = inoc_model_tab,
  inoculation_overall = inoc_overall,
  field_plot_data = field_plot,
  field_interaction = field_interaction_tab,
  field_by_year = field_by_year,
  height_rate = height_rate,
  height_identifiability = height_identifiability,
  methods_note = tibble::tibble(
    item = c("inoculation_unit", "field_unit", "height_unit", "batch", "field"),
    definition = c(
      "Pot (image); plants within a pot averaged before inference",
      "Year x Field x focal-neighbour x condition; 10 scored plants per unit averaged/aggregated",
      "Focal-neighbour combination; repeated surveys as longitudinal observations",
      "Inoculation batch included as fixed effect (2-3 levels)",
      "Field included as fixed effect; random-intercept sensitivity reported"
    )
  )
)

bioRtools::write_data(out_list, file = path_out)
cat("\n[16] 输出：", path_out, "\n")

# ============================================================
# 五、Figure 5：pot-level 分解 + 田间逐年效应（无旧 reciprocity 过强结论）
# ============================================================

inoc_plot_tab <- inoc_model_tab %>%
  dplyr::select(disease, focal_F, neighbor_F, interaction_F,
                focal_p, neighbor_p, interaction_p) %>%
  tidyr::pivot_longer(
    cols = c(focal_F, neighbor_F, interaction_F),
    names_to = "term", values_to = "statistic"
  ) %>%
  dplyr::mutate(
    term = factor(
      term,
      levels = c("focal_F", "neighbor_F", "interaction_F"),
      labels = c("Focal identity", "Neighbour identity", "Focal x neighbour")
    ),
    p_value = dplyr::case_when(
      term == "Focal identity" ~ focal_p,
      term == "Neighbour identity" ~ neighbor_p,
      TRUE ~ interaction_p
    ),
    disease = factor(disease, levels = c("blast", "brown spot"),
                     labels = c("Rice blast", "Brown spot")),
    label = format_p(p_value, 2)
  )

p_inoc <- inoc_plot_tab %>%
  ggplot2::ggplot(ggplot2::aes(x = term, y = statistic, fill = term)) +
  ggplot2::geom_col(width = 0.68, show.legend = FALSE) +
  ggplot2::geom_text(
    ggplot2::aes(label = label),
    vjust = -0.25, size = 2.7, lineheight = 0.95
  ) +
  ggplot2::facet_wrap(~disease) +
  bioRtools::scale_fill_research() +
  ggplot2::scale_y_continuous(
    limits = c(0, NA),
    expand = ggplot2::expansion(mult = c(0, 0.24))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "F value (pot means)",
    title = "Pot-level decomposition",
    subtitle = "Incomplete 4 x 4 mixture design; interaction df = 5"
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1, size = 7),
    strip.text = ggplot2::element_text(face = "bold", size = 8),
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6.5)
  )

plot_order <- c("Blast", "BrownSpot", "Smut", "PH", "KGW")
ref_lines <- tibble::tibble(
  category = factor(plot_order, levels = plot_order),
  ref = c(1, 1, 1, 0, 0)
)

field_plot_tab <- field_by_year %>%
  dplyr::mutate(
    category = factor(category, levels = plot_order),
    year = factor(year)
  )

inter_lab <- field_interaction_tab %>%
  dplyr::mutate(
    category = factor(category, levels = plot_order),
    label = paste0("condition x Year\n", format_p(interaction_p, 3))
  )

p_field <- field_plot_tab %>%
  ggplot2::ggplot(ggplot2::aes(x = estimate, y = year)) +
  ggplot2::geom_vline(
    data = ref_lines,
    ggplot2::aes(xintercept = ref),
    linetype = "dashed", color = "grey40"
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
    height = 0.2, color = "grey35"
  ) +
  ggplot2::geom_point(ggplot2::aes(color = year), size = 2.6) +
  ggplot2::geom_text(
    data = inter_lab,
    ggplot2::aes(x = -Inf, y = Inf, label = label),
    hjust = -0.05, vjust = 1.4, size = 2.5, lineheight = 0.95,
    inherit.aes = FALSE
  ) +
  ggplot2::facet_wrap(~category, scales = "free_x", nrow = 1) +
  ggplot2::scale_y_discrete(limits = rev(levels(field_plot_tab$year))) +
  bioRtools::scale_color_research() +
  ggplot2::labs(
    x = "Effect per year (disease: OR; trait: mean difference)",
    y = "Year",
    title = "Field effects by year",
    subtitle = "Plot-level aggregates"
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  bioRtools::theme_prism() +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 8, face = "bold"),
    strip.clip = "off",
    legend.position = "none",
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6.5),
    axis.text.x = ggplot2::element_text(size = 7)
  )

p_fig <- patchwork::wrap_plots(p_inoc, p_field, ncol = 1, heights = c(1, 1.2))

bioRtools::write_data(p_fig, file = paste0(path_fig, ".pdf"),
                      width = 17, height = 16, units = "cm")
bioRtools::write_data(p_fig, file = paste0(path_fig, ".png"),
                      width = 17, height = 16, units = "cm", dpi = 600)
cat("[16] 输出：", paste0(path_fig, ".pdf / .png"), "\n")

cat("\n[16] sessionInfo()\n")
sessionInfo()
