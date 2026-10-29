# ============================================================
# 06_设计分解与稳健性 / 01_接种_盆内聚类与设计分解.R
# 独立分析（不属于旧版流程）：温室接种实验的盆内聚类校正 + focal/neighbor 双向分解 + reciprocity 检验
# 运行：/opt/homebrew/bin/Rscript "<root>/code/06_设计分解与稳健性/01_接种_盆内聚类与设计分解.R"
# 命名沿革：本脚本原名 01_接种_伪重复与设计分解.R（2026-09-22 创建），因"伪重复"表述不准确，
#   经用户确认于同日改名为 01_接种_盆内聚类与设计分解.R，输出表同步改名；新旧名映射登记在 项目文件记录.md。
#
# 【为什么要做这一步】
#   采样结构（2026-09-22 由用户确认，并已用数据复核）：**每一株就是一个个体，一株只取一片叶**，
#   因此表中一行 = 一株植物（不是"同一株的重复叶片"），叶级并不存在"同一株重复采样"式的伪重复。
#   但同一张图像（同一盆）内的多株植物生长在同一个盆里、共享土壤与邻体环境，彼此并非独立观测：
#   实测每个 image 含 6–14 株、leaf_number 在 image 内唯一、且一个 image 只对应一个 focal × neighbor 组合，
#   盆间方差显著（胡麻叶斑病 ICC = 0.12）。因此这里做的是**盆内聚类校正**：
#   把 pot/image 作为随机截距后，有效独立单位从"株数"降到接近"盆数"（1043 株 → 96 盆），
#   旧版按株数给出 df ≈ 1039 会高估精度。这正是本脚本 Q1 要量化的东西。
#   同时旧版按 focal 逐个切片比较，无法回答"邻体身份本身贡献了多少"、"同一对品种的双向响应是否对称"。
#   本脚本回答三个新问题：
#   Q1 聚类校正：按株（旧版）vs 加盆随机截距 vs 按盆聚合，P 值与结论是否改变？
#   Q2 设计分解：在 Mix 组合中，focal 身份、neighbor 身份分别解释多少病斑密度变异？交互项是否显著？
#   Q3 reciprocity：对每个无序品种对 (i, j)，i 受 j 影响与 j 受 i 影响是否对称（幅度与方向）？
#
# 输入：<root>/data/01_分析输入表/01_接种结果_合并表.xlsx
# 输出：<root>/result/06_设计分解与稳健性/01_表_接种_盆内聚类与设计分解.xlsx
#       <root>/result/06_设计分解与稳健性/01_FigureS8_接种设计分解.pdf / .png
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录见 code/06_设计分解与稳健性/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_in <- file.path(project_root, "data", "01_分析输入表", "01_接种结果_合并表.xlsx")
out_dir <- file.path(project_root, "result", "06_设计分解与稳健性")
path_out <- file.path(out_dir, "01_表_接种_盆内聚类与设计分解.xlsx")

stopifnot(file.exists(path_in))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

bioRtools::read_data(path_in) -> d1
stopifnot(nrow(d1) == 1659)

d1 %>%
  dplyr::mutate(
    batch = as.factor(batch),
    group = factor(group, levels = c("Pure", "Mix")),
    focal = as.factor(focal),
    neighbor = as.factor(neighbor),
    image = as.factor(image),
    disease = as.factor(disease)
  ) -> df

cat("[01] 数据：", nrow(df), " 片叶；", dplyr::n_distinct(df$image), " 个 pot 图像\n")
cat("[01] 每个 pot 的叶片数：", min(dplyr::count(df, image)$n), "–", max(dplyr::count(df, image)$n), "\n")

# ---------- Q1：聚类校正（按株 vs 加盆随机截距 vs 按盆聚合）----------
unit_tab <- purrr::map_dfr(levels(df$disease), function(ds) {
  d <- df %>% dplyr::filter(disease == ds) %>% droplevels()

  m_leaf <- lmerTest::lmer(lesion_number_per_cm2 ~ group + (1 | batch), data = d)
  m_pot <- lmerTest::lmer(lesion_number_per_cm2 ~ group + (1 | batch) + (1 | image), data = d)

  # 盆级聚合：同一 image（盆）内的多株植物取均值，再以盆为单位检验（有效 n = 盆数）
  d_pot <- d %>%
    dplyr::group_by(image, batch, group) %>%
    dplyr::summarise(lesion_pot = mean(lesion_number_per_cm2), n_leaf = dplyr::n(), .groups = "drop")
  m_agg <- lmerTest::lmer(lesion_pot ~ group + (1 | batch), data = d_pot)

  vc_pot <- lme4::VarCorr(m_pot)
  tibble::tibble(
    disease = ds,
    n_leaf = nrow(d),
    n_pot = dplyr::n_distinct(d$image),
    estimate_leaf = summary(m_leaf)$coefficients["groupMix", "Estimate"],
    p_leaf = summary(m_leaf)$coefficients["groupMix", "Pr(>|t|)"],
    df_leaf = summary(m_leaf)$coefficients["groupMix", "df"],
    estimate_pot_random = summary(m_pot)$coefficients["groupMix", "Estimate"],
    p_pot_random = summary(m_pot)$coefficients["groupMix", "Pr(>|t|)"],
    df_pot_random = summary(m_pot)$coefficients["groupMix", "df"],
    var_image = as.numeric(vc_pot$image),
    var_resid = attr(vc_pot, "sc")^2,
    icc_pot = var_image / (var_image + attr(vc_pot, "sc")^2),
    estimate_pot_agg = summary(m_agg)$coefficients["groupMix", "Estimate"],
    p_pot_agg = summary(m_agg)$coefficients["groupMix", "Pr(>|t|)"],
    df_pot_agg = summary(m_agg)$coefficients["groupMix", "df"]
  )
})

cat("\n[01] Q1 聚类校正对比（按株 / 加盆随机截距 / 按盆聚合；一行 = 一株、一株一片叶）：\n")
print(as.data.frame(unit_tab %>% dplyr::select(disease, n_leaf, n_pot, estimate_leaf, p_leaf, df_leaf,
                                              p_pot_random, df_pot_random, icc_pot, p_pot_agg, df_pot_agg)),
      row.names = FALSE)

# ---------- Q2：Mix 组合内的 focal / neighbor 分解 ----------
decomp_tab <- purrr::map_dfr(levels(df$disease), function(ds) {
  d <- df %>%
    dplyr::filter(disease == ds, group == "Mix") %>%
    droplevels()

  m_full <- lmerTest::lmer(lesion_number_per_cm2 ~ focal * neighbor + (1 | batch), data = d)
  aov_full <- as.data.frame(anova(m_full)) %>% tibble::rownames_to_column(var = "term")

  # 仅主效应的部分 R²（用两两模型的残差平方和之差衡量）
  ss_total <- sum((d$lesion_number_per_cm2 - mean(d$lesion_number_per_cm2))^2)
  m_focal <- stats::lm(lesion_number_per_cm2 ~ focal, data = d)
  m_neigh <- stats::lm(lesion_number_per_cm2 ~ neighbor, data = d)

  tibble::tibble(
    disease = ds,
    term = aov_full$term,
    df_num = aov_full$NumDF,
    den_df = aov_full$DenDF,
    f_value = aov_full$`F value`,
    p_value = aov_full$`Pr(>F)`,
    r2_partial_focal = 1 - sum(stats::residuals(m_focal)^2) / ss_total,
    r2_partial_neighbor = 1 - sum(stats::residuals(m_neigh)^2) / ss_total
  )
})

cat("\n[01] Q2 Mix 组合内的 focal / neighbor / 交互项检验：\n")
print(as.data.frame(decomp_tab %>% dplyr::select(disease, term, df_num, den_df, f_value, p_value)),
      row.names = FALSE)

# ---------- Q3：reciprocity ----------
# 每个无序品种对 (i, j)：effect_ij = mean(Mix i-j) − mean(Pure i)
# 效应量在**按 病害 × 批次 Z-score 后**计算，以消除批次间感染压力差异（与主题 01 的标准化口径一致）
df %>%
  dplyr::group_by(disease, batch) %>%
  dplyr::mutate(lesion_norm = as.vector(scale(lesion_number_per_cm2))) %>%
  dplyr::ungroup() -> df.z

df.z %>%
  dplyr::filter(group == "Pure") %>%
  dplyr::group_by(disease, focal) %>%
  dplyr::summarise(pure_mean = mean(lesion_norm), .groups = "drop") -> pure_mean

df.z %>%
  dplyr::filter(group == "Mix") %>%
  dplyr::group_by(disease, focal, neighbor) %>%
  dplyr::summarise(mix_mean = mean(lesion_norm), .groups = "drop") %>%
  dplyr::left_join(pure_mean, by = c("disease", "focal")) %>%
  dplyr::mutate(effect = mix_mean - pure_mean) -> pair_long

# 配对：把 (i,j) 与 (j,i) 拼在一起
pair_long %>%
  dplyr::rename(focal_ab = focal, neighbor_ab = neighbor, effect_ab = effect) -> a_side

pair_long %>%
  dplyr::rename(neighbor_ba = focal, focal_ba = neighbor, effect_ba = effect) -> b_side

recip_tab <- a_side %>%
  dplyr::inner_join(b_side, by = c("disease", "focal_ab" = "focal_ba", "neighbor_ab" = "neighbor_ba")) %>%
  dplyr::filter(as.integer(focal_ab) < as.integer(neighbor_ab)) %>%
  dplyr::mutate(pair = paste0(focal_ab, "-", neighbor_ab)) %>%
  dplyr::select(disease, pair, focal = focal_ab, neighbor = neighbor_ab,
                effect_focal = effect_ab, effect_neighbor = effect_ba) %>%
  dplyr::arrange(disease, pair)

recip_summary <- recip_tab %>%
  dplyr::group_by(disease) %>%
  dplyr::group_modify(~ {
    ct <- stats::cor.test(.x$effect_focal, .x$effect_neighbor)
    # 对称性检验：观察到的平均绝对差 vs 随机翻转符号的置换分布
    obs <- mean(abs(.x$effect_focal - .x$effect_neighbor))
    n <- nrow(.x)
    perm <- replicate(999, {
      signs <- sample(c(-1, 1), n, replace = TRUE)
      mean(abs(.x$effect_focal - signs * .x$effect_neighbor))
    })
    tibble::tibble(
      n_pair = n,
      r = unname(ct$estimate),
      cor_p = ct$p.value,
      mean_abs_diff = obs,
      perm_p_symmetry = (sum(perm >= obs) + 1) / (length(perm) + 1),
      mean_effect_focal = mean(.x$effect_focal),
      mean_effect_neighbor = mean(.x$effect_neighbor)
    )
  }) %>%
  dplyr::ungroup()

cat("\n[01] Q3 reciprocity（12 个品种对）：\n")
print(as.data.frame(recip_tab), row.names = FALSE)
cat("\n")
print(as.data.frame(recip_summary), row.names = FALSE)
cat("\n[01] 绘图数据检查：行数 =", nrow(recip_tab),
    "；x/y 缺失数 =", sum(is.na(recip_tab$effect_focal)) + sum(is.na(recip_tab$effect_neighbor)), "\n")

bioRtools::write_data(
  list(
    inference_unit = unit_tab,
    design_decomposition = decomp_tab,
    reciprocity_pairs = recip_tab,
    reciprocity_summary = recip_summary
  ),
  file = path_out
)
cat("[01] 输出：", path_out, "\n")

# ---------- 图：Q1 P 值对比 + Q3 reciprocity ----------
p.unit <- unit_tab %>%
  dplyr::select(disease, p_leaf, p_pot_random, p_pot_agg) %>%
  tidyr::pivot_longer(cols = -disease, names_to = "model", values_to = "p_value") %>%
  dplyr::mutate(
    model = factor(model,
                   levels = c("p_leaf", "p_pot_random", "p_pot_agg"),
                   labels = c("Per plant\n(old inference)", "Per plant +\npot random effect", "Pot-level\naggregation")),
    p_label = paste0("P = ", formatC(p_value, format = "f", digits = 3))
  ) %>%
  ggplot2::ggplot(aes(x = model, y = p_value, fill = disease)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.65) +
  ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30") +
  ggplot2::geom_text(aes(label = p_label),
                     position = ggplot2::position_dodge(width = 0.7), vjust = -0.4, size = 3) +
  bioRtools::scale_fill_research() +
  ggplot2::labs(x = NULL, y = "P value of Pure vs Mix effect", fill = "Disease") +
  bioRtools::theme_prism() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8)) -> p.unit

p.recip <- recip_tab %>%
  dplyr::filter(is.finite(effect_focal), is.finite(effect_neighbor)) %>%
  ggplot2::ggplot(aes(x = effect_focal, y = effect_neighbor)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  ggplot2::geom_point(size = 2.6, color = "#1e4e8c") +
  ggrepel::geom_text_repel(aes(label = pair), size = 2.8, max.overlaps = 20) +
  ggplot2::facet_wrap(~disease, scales = "free") +
  ggplot2::labs(
    x = "Effect on focal (lesions/cm2, Mix - Pure)",
    y = "Effect on partner (lesions/cm2, Mix - Pure)"
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold")) -> p.recip

# 注意：不能用 aplot::insert_right()——它会强制对齐两个子图的 y 轴，
# 导致 reciprocity 面板中 y 值超出左图 y 范围的点被当作超出 scale 范围而静默删除
# （实测 12 个点中 7 个被删）。改用 patchwork 拼图，不跨面板对齐坐标轴。
p.s8 <- patchwork::wrap_plots(p.unit, p.recip, widths = c(1, 1.6))

bioRtools::write_data(p.s8, file = file.path(out_dir, "01_FigureS8_接种设计分解.pdf"),
                      width = 14, height = 6, units = "in")
bioRtools::write_data(p.s8, file = file.path(out_dir, "01_FigureS8_接种设计分解.png"),
                      width = 14, height = 6, units = "in", dpi = 500)
cat("[01] 输出：01_FigureS8_接种设计分解.pdf / .png\n")

cat("\n[01] sessionInfo()\n")
sessionInfo()
