# ============================================================
# 02_agronomic_traits / 09_agronomic_traits_plotting_Figure1_submission_format.R
# 作用：**submission_format_version** Figure 1 面板（2026-10-18 新增；plotting口径anddata同 08 号脚本）
#       and 08 号脚本的差别只在版面：
#         ① 输出尺寸 = 面板在组图中的**final rect 尺寸**（scale = 1），保证图内字号
#            落版后仍是design字号（08 号在大画布上导出、拼图时整体缩到 174 mm，
#            图内文字final只剩 1.9–4.6 pt，低于期刊下限）；
#         ② 去掉面板内 caption（"P 值是计划对比""误差棒定义"等统一写进正文
#            Figure legend），避免小面板被 caption 挤空；
#         ③ 1E 分面标签缩写，保证 83 mm 宽面板不越界；
#         ④ 输出写到 result/08_submission_format图，原 result/02_agronomic_traits 只读不改。
#       1A height_change_index（cm；以 15 次调查为单元的 95% CI）
#       1B tiller_number / 1C 单株yield（原始尺度箱线 + 抖动散点，**标注and自身单作的计划对比 P**）
#       1D 三性状表型 PCA（16 个组合；Pure → Mix 连线）
#       1E 三性状 Pure vs Mix 整体comparison（每个分面自带单位）
#       1F 遗传距离and表型 → **按审稿意见 M12 移出main_figures**，另存为supplementary_figures（19）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/02_agronomic_traits/09_agronomic_traits_plotting_Figure1_submission_format.R"
#
# 输入：<root>/data/processed_inputs/03_agronomic_traits_height_tiller_yield.xlsx
#       <root>/result/02_agronomic_traits/13_表_planned_contrasts_and_raw_scale_effects.xlsx（planned_contrasts_traits、height_index_contrasts）
#       <root>/result/02_agronomic_traits/03_表_表型PCA.xlsx（scores、variance_explained）
#       <root>/result/02_agronomic_traits/05_表_genetic_distance_and_phenotype_correlation.xlsx（data_joined、correlation）
# 输出：<root>/result/08_submission_format图/14–18_Figure1A–1E_*_new_definition.pdf / .png
#       <root>/result/08_submission_format图/19_FigureS_遗传距离and表型_new_definition.pdf / .png
# 其中 <root> = <repository_root>
#
# and旧版图的差异（对应审稿意见 M4/M10/M11/M12 and"图内单位/配色"条目）：
#   ① 1A 的误差棒and显著性改为"15 次调查为单元"，不再用行数当 n；
#   ② 1B/1C 的标注由"两两 Tukey 字母"改为"每个 Mix 邻体 vs 该 focal 自身单作"的计划对比 P（Dunnett）；
#   ③ 每个分面/面板自带单位（tillers per plant、g per plant、cm）；
#   ④ 配色and图例由 ggplot 直接generate（不再用文字描述颜色），避免图例and图内颜色不符；
#   ⑤ 遗传距离—表型analysis移出main_figures。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
out_dir <- file.path(project_root, "result", "02_agronomic_traits")
out_pub_dir <- file.path(project_root, "result", "08_submission_format图")  # submission_format面板输出目录
dir.create(out_pub_dir, showWarnings = FALSE, recursive = TRUE)
path_traits <- file.path(project_root, "data", "processed_inputs", "03_agronomic_traits_height_tiller_yield.xlsx")
path_planned <- file.path(out_dir, "13_表_planned_contrasts_and_raw_scale_effects.xlsx")
path_pca <- file.path(out_dir, "03_表_表型PCA.xlsx")
path_cor <- file.path(out_dir, "05_表_genetic_distance_and_phenotype_correlation.xlsx")

stopifnot(file.exists(path_traits), file.exists(path_planned), file.exists(path_pca), file.exists(path_cor))

focus_levels <- c("15", "49", "69", "87")

format_p <- function(p, digits = 3) {
  dplyr::if_else(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, format = "f", digits = digits)))
}

# submission_format：紧凑写法（无空格、2 位小数），保证 170 mm 宽面板里 4 个邻体的 P 值不叠字；
#          精确值见正文 Figure legend 指向的 Table S1
format_p_compact <- function(p, digits = 2) {
  dplyr::if_else(p < 0.001, "P<0.001", paste0("P=", formatC(p, format = "f", digits = digits)))
}

bioRtools::read_data(path_planned, sheet = "planned_contrasts_traits") -> planned
bioRtools::read_data(path_planned, sheet = "height_index_contrasts") -> idx_con
bioRtools::read_data(path_pca, sheet = "scores") -> pca_scores
bioRtools::read_data(path_pca, sheet = "variance_explained") -> pca_var
bioRtools::read_data(path_cor, sheet = "data_joined") -> cor_join
bioRtools::read_data(path_cor, sheet = "correlation") -> cor_tab
bioRtools::read_data(path_traits) -> df.raw

# ---- Figure 1A：height_change_index（15 次调查为单元）----
idx_con %>%
  dplyr::mutate(focal = factor(focal, levels = focus_levels),
                neighbor = factor(neighbor, levels = focus_levels)) -> idx_plot

p.1a <- idx_plot %>%
  ggplot2::ggplot(aes(x = neighbor, y = index_mean)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  ggplot2::geom_col(fill = "grey88", width = 0.62) +
  ggplot2::geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.18, color = "grey25") +
  # P 值标注：用 ggrepel 自动避让（同一 facet 内相邻柱的 P 标签曾互相压住 12%）
  ggrepel::geom_text_repel(aes(y = ci_upper, label = format_p(p_value)),
                           vjust = -0.4, size = 2.4, color = "grey20",
                           direction = "y", box.padding = 0.18, min.segment.length = 0.2,
                           segment.size = 0.2, seed = 12345) +
  ggplot2::facet_wrap(~focal, nrow = 1, scales = "free_y") +
  ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 4),
                              expand = ggplot2::expansion(mult = c(0.12, 0.30))) +
  ggplot2::labs(x = "Neighbour accession", y = "Height index (cm)") +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 8),
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))

# submission_format：尺寸一律按**final rect 尺寸（cm）**给，落版即design字号（scale = 1）
save_panel <- function(p, pid, label, width_cm, height_cm) {
  bioRtools::write_data(p, file = file.path(out_pub_dir, paste0(pid, "_", label, ".pdf")),
                        width = width_cm, height = height_cm, units = "cm")
  bioRtools::write_data(p, file = file.path(out_pub_dir, paste0(pid, "_", label, ".png")),
                        width = width_cm, height = height_cm, units = "cm", dpi = 600)
  cat("[09] 输出：", pid, "_", label, ".pdf / .png\n", sep = "")
}

save_panel(p.1a, "14", "Figure1A_height_change_index_new_definition", 17.0, 5.0)

# ---- Figure 1B / 1C：原始尺度箱线 + 计划对比 P ----
plot_trait <- function(pheno, pid, label, unit_label, width_in, height_in) {
  df.raw %>%
    dplyr::filter(phenotype == pheno) %>%
    dplyr::mutate(focal = factor(focal, levels = focus_levels),
                  neighbor = factor(neighbor, levels = focus_levels),
                  group = factor(group, levels = c("Pure", "Mix"))) -> d

  ann <- planned %>%
    dplyr::filter(phenotype == pheno) %>%
    dplyr::mutate(focal = factor(focal, levels = focus_levels),
                  neighbor = factor(mix_neighbor, levels = focus_levels))

  ymax <- d %>% dplyr::group_by(focal) %>%
    dplyr::summarise(y = max(value, na.rm = TRUE), .groups = "drop")

  # submission_format（2026-10-19 修）：同一 focal 内 3 个 Mix 邻体（邻体 = 该 focal 自身的那个是单作对照，P 为 NA，不参and）
  #   的 P 值**合并成一行**放在分面顶部——
  #   逐柱标记时每个标注宽约 26 pt、柱间距仅约 21 pt，必然互相压住（实测重叠 5.4 pt）；
  #   合并成单个字符串后总宽约 90 pt < 分面宽约 116 pt，既保留数值又不可能叠字；
  #   顺序 = 按邻体 accession 升序（= 除该 focal 自身之外的三个 accession），显著者（BH-FDR < 0.05）后加 "*"。
  rng_tab <- d %>% dplyr::group_by(focal) %>%
    dplyr::summarise(ymax = max(value, na.rm = TRUE),
                     rng = diff(range(value, na.rm = TRUE)), .groups = "drop")
  ann_line <- ann %>% dplyr::left_join(rng_tab, by = "focal") %>%
    dplyr::arrange(focal, neighbor) %>%
    dplyr::group_by(focal) %>%
    dplyr::summarise(ymax = dplyr::first(ymax), rng = dplyr::first(rng),
                     lab = paste0("P: ", paste0(
                       formatC(p_value, format = "f", digits = 2),
                       dplyr::if_else(!is.na(p_fdr) & p_fdr < 0.05, "*", ""),
                       collapse = " / ")),
                     .groups = "drop") %>%
    dplyr::mutate(y = ymax + 0.08 * rng)

  d %>%
    ggplot2::ggplot(aes(x = neighbor, y = value)) +
    ggplot2::geom_boxplot(aes(fill = group), width = 0.62, outlier.shape = NA, alpha = 0.85) +
    ggplot2::geom_jitter(aes(color = group), width = 0.16, alpha = 0.45, size = 1) +
    ggplot2::geom_text(data = ann_line, aes(x = 2.5, y = y, label = lab),
                       inherit.aes = FALSE, size = 2.4, color = "grey20", vjust = 0) +
    ggplot2::facet_wrap(~focal, nrow = 1, scales = "free_y") +
    ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 3),
                                expand = ggplot2::expansion(mult = c(0.05, 0.22))) +
    ggplot2::scale_fill_manual(values = c(Pure = "grey85", Mix = "#1e4e8c"), name = NULL) +
    ggplot2::scale_color_manual(values = c(Pure = "grey45", Mix = "#1e4e8c"), name = NULL) +
  ggplot2::labs(x = "Neighbour accession", y = unit_label) +
    bioRtools::theme_prism() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"),
                 legend.position = "top",
                 legend.key.size = ggplot2::unit(0.35, "cm"),
                 legend.margin = ggplot2::margin(0, 0, 0, 0),
                 legend.box.spacing = ggplot2::unit(0.08, "cm")) -> p

  save_panel(p, pid, label, width_in, height_in)
}

plot_trait("tiller number", "15", "Figure1B_tiller_number_new_definition", "Tiller number per plant", 17.0, 4.6)
plot_trait("yield", "16", "Figure1C_单株yield_new_definition", "Yield per plant (g)", 17.0, 4.6)

# ---- Figure 1D：表型 PCA ----
pc_lab <- function(i) {
  v <- pca_var$variance_explained[pca_var$component == paste0("PC", i)]
  paste0("PC", i, " (", formatC(100 * v, format = "f", digits = 2), "%)")
}

pca_scores %>%
  dplyr::mutate(focal = factor(focal, levels = focus_levels),
                group = factor(group, levels = c("Pure", "Mix"))) %>%
  dplyr::group_by(focal) %>%
  dplyr::mutate(center_x = mean(PC1), center_y = mean(PC2)) %>%
  dplyr::ungroup() -> d.pca

p.1d <- d.pca %>%
  ggplot2::ggplot(aes(x = PC1, y = PC2)) +
  ggplot2::geom_segment(aes(xend = center_x, yend = center_y, color = focal),
                        linetype = "dashed", linewidth = 0.6, alpha = 0.75) +
  ggplot2::geom_point(aes(color = focal, shape = group), size = 3.6) +
  ggplot2::scale_shape_manual(values = c(Pure = 16, Mix = 17), name = "Planting mode") +
  bioRtools::scale_color_research() +
  ggplot2::labs(x = pc_lab(1), y = pc_lab(2), color = "Focal accession") +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "right",
                 legend.key.size = ggplot2::unit(0.32, "cm"),
                 legend.text = ggplot2::element_text(size = 7),
                 legend.title = ggplot2::element_text(size = 7.5),
                 legend.spacing.y = ggplot2::unit(0.05, "cm"),
                 legend.margin = ggplot2::margin(0, 0, 0, 0))

save_panel(p.1d, "17", "Figure1D_表型PCA_new_definition", 8.3, 5.2)

# ---- Figure 1E：三性状 Pure vs Mix（每个分面自带单位）----
# 注意：13 表的 height_index_contrasts 只包含混种组合（供计划对比用），做 1E 需要把单作也带上，
# 因此这里直接从 01 表的 index_long reanalysis"每次调查一个均值"的指数（and 1A 的test单元一致）。
bioRtools::read_data(file.path(out_dir, "01_表_height_change_index.xlsx"), sheet = "index_long") %>%
  dplyr::mutate(focal = factor(focal, levels = focus_levels),
                neighbor = factor(neighbor, levels = focus_levels)) %>%
  dplyr::group_by(focal, neighbor, time) %>%
  dplyr::summarise(index_by_survey = mean(index, na.rm = TRUE), .groups = "drop") %>%
  dplyr::transmute(combo = paste0(focal, "-", neighbor), focal, neighbor,
                   value = index_by_survey,
                   trait = "Plant height change index (cm)") -> height_index

df.raw %>%
  dplyr::filter(phenotype %in% c("tiller number", "yield")) %>%
  dplyr::transmute(combo = paste0(focal, "-", neighbor), focal, neighbor, value,
                   trait = dplyr::if_else(phenotype == "tiller number",
                                          "Tiller number (per plant)", "Yield (g per plant)")) -> other_traits

d.1e <- dplyr::bind_rows(height_index, other_traits) %>%
  dplyr::mutate(focal = factor(focal, levels = focus_levels),
                neighbor = factor(neighbor, levels = focus_levels),
                group = factor(dplyr::if_else(focal == neighbor, "Pure", "Mix"),
                               levels = c("Pure", "Mix")),
                trait = factor(trait, levels = c("Plant height change index (cm)",
                                                 "Tiller number (per plant)",
                                                 "Yield (g per plant)"),
                               # 2026-10-19 修（硬伤）：分面title原为 "Height index (cm)" / "Tillers per plant" /
                               #   "Yield (g per plant)"，在 83 mm 宽 3 分面面板里**比 strip 框还宽**，
                               #   被 ggplot2 的 strip.clip（默认 "on"）**两侧切掉**（实测渲染成
                               #   "eight index (cm" / "Tillers per plan" / "ield (g per plan"）。
                               #   改为短标签（单位在正文 Figure 1 图注里给全称），并同时关掉 strip.clip。
                               labels = c("Height (cm)", "Tillers", "Yield (g)")))

p_tab <- d.1e %>%
  dplyr::group_by(trait) %>%
  rstatix::t_test(value ~ group) %>%
  dplyr::ungroup() %>%
  dplyr::select(trait, p)

d.1e %>%
  dplyr::left_join(p_tab, by = "trait") %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(y_pos = max(value, na.rm = TRUE) + 0.08 * diff(range(value, na.rm = TRUE))) %>%
  dplyr::ungroup() -> d.1e.ann

p.1e <- d.1e.ann %>%
  ggplot2::ggplot(aes(x = group, y = value)) +
  ggplot2::geom_violin(aes(fill = trait), alpha = 0.45, color = NA, show.legend = FALSE) +
  ggplot2::geom_boxplot(aes(fill = trait), width = 0.28, outlier.shape = NA,
                        alpha = 0.85, show.legend = FALSE) +
  ggplot2::geom_jitter(width = 0.12, alpha = 0.25, size = 0.6, color = "grey30") +
  ggplot2::geom_text(aes(x = 1.5, y = y_pos, label = format_p(p)), size = 3.2, check_overlap = TRUE) +
  ggplot2::facet_wrap(~trait, scales = "free_y") +
  ggplot2::scale_x_discrete(limits = c("Pure", "Mix")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.16))) +
  bioRtools::scale_fill_research() +
  ggplot2::labs(x = NULL, y = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 6.5),
                 strip.clip = "off",   # 2026-10-19：确保分面title不会被 strip 框裁掉
                 legend.position = "none")

save_panel(p.1e, "18", "Figure1E_三性状整体comparison_new_definition", 8.3, 5.2)

# ---- 遗传距离and表型（按 M12 移出main_figures，作supplementary_figures）----
cor_lab <- cor_tab %>%
  dplyr::filter(dataset == "Mix only (n = 12)") %>%
  dplyr::mutate(
    trait = factor(trait, levels = c("mean_height", "mean_tiller", "mean_yield"),
                   labels = c("Plant height index (cm)", "Tiller number (per plant)", "Yield (g per plant)")),
    label = paste0("italic(r) == ", formatC(r, format = "f", digits = 2), "~~", format_p(p_value))
  )

cor_join %>%
  dplyr::filter(group == "Mix") %>%
  tidyr::pivot_longer(cols = c(mean_height, mean_tiller, mean_yield),
                      names_to = "trait", values_to = "value") %>%
  dplyr::mutate(trait = factor(trait, levels = c("mean_height", "mean_tiller", "mean_yield"),
                               labels = c("Plant height index (cm)", "Tiller number (per plant)",
                                          "Yield (g per plant)"))) %>%
  ggplot2::ggplot(aes(x = genetic_distance, y = value)) +
  ggplot2::geom_point(size = 2.6, color = "#1e4e8c") +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "grey30", fill = "grey85") +
  ggplot2::geom_text(data = cor_lab, aes(x = -Inf, y = Inf, label = label),
                     hjust = -0.08, vjust = 1.6, size = 3.2, parse = TRUE, inherit.aes = FALSE) +
  ggplot2::facet_wrap(~trait, scales = "free_y") +
  ggplot2::labs(x = "Pairwise genetic distance", y = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 8)) -> p.s

save_panel(p.s, "19", "FigureS_遗传距离and表型_new_definition", 17.0, 6.95)

cat("\n[08] sessionInfo()\n")
sessionInfo()
