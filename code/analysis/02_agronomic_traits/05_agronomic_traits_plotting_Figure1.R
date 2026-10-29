# ============================================================
# 02_agronomic_traits / 05_agronomic_traits_plotting_Figure1.R
# 作用：绘制 Figure 1A–1F
#       1A height_change_index（按 focal 分组、neighbor 并列，虚线 = 0）
#       1B tiller_number / 1C 单株yield（各 focal 的邻体comparison）
#       1D 三性状表型 PCA（Pure → Mix 连线）
#       1E 三性状 Pure vs Mix 整体comparison
#       1F 遗传距离and三性状（阴性results，Mix 组合 n = 12）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/02_agronomic_traits/05_agronomic_traits_plotting_Figure1.R"
#
# 输入：<root>/result/02_agronomic_traits/01_表_height_change_index.xlsx（index_long、anova_posthoc）
#       <root>/result/02_agronomic_traits/02_表_tiller_and_yield_statistics.xlsx（anova_posthoc）
#       <root>/result/02_agronomic_traits/03_表_表型PCA.xlsx（scores、variance_explained）
#       <root>/result/02_agronomic_traits/05_表_genetic_distance_and_phenotype_correlation.xlsx（data_joined、correlation）
# 输出：<root>/result/02_agronomic_traits/06–11_Figure1A–1F_*.pdf / .png
# 其中 <root> = <repository_root>
#
# plotting约定：bioRtools::theme_prism()；scale_fill_research() / scale_color_research()；
#          bioRtools::write_data(..., units = "in", dpi = 500)；图内文字English。
#          P 值统一 format_p() 定点写法（P < 0.001 / P = 0.012），不用科学计数法。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/02_agronomic_traits/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
out_dir <- file.path(project_root, "result", "02_agronomic_traits")

path_idx <- file.path(out_dir, "01_表_height_change_index.xlsx")
path_tiller <- file.path(out_dir, "02_表_tiller_and_yield_statistics.xlsx")
path_pca <- file.path(out_dir, "03_表_表型PCA.xlsx")
path_cor <- file.path(out_dir, "05_表_genetic_distance_and_phenotype_correlation.xlsx")

stopifnot(file.exists(path_idx), file.exists(path_tiller), file.exists(path_pca), file.exists(path_cor))

format_p <- function(p, digits = 3) {
  dplyr::if_else(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, format = "f", digits = digits)))
}

bioRtools::read_data(path_idx, sheet = "index_long") -> df.idx
bioRtools::read_data(path_idx, sheet = "anova_posthoc") -> posthoc.height
bioRtools::read_data(path_tiller, sheet = "anova_posthoc") -> posthoc.tiller
bioRtools::read_data(path_pca, sheet = "scores") -> df.score
bioRtools::read_data(path_pca, sheet = "variance_explained") -> var.exp
bioRtools::read_data(path_cor, sheet = "data_joined") -> df.join
bioRtools::read_data(path_cor, sheet = "correlation") -> cor.tab

focus_levels <- c("15", "49", "69", "87")

df.idx %>%
  dplyr::mutate(
    focal = factor(focal, levels = focus_levels),
    neighbor = factor(neighbor, levels = focus_levels),
    group = factor(group, levels = c("Pure", "Mix"))
  ) -> df.idx

# 字母位置：每个 focal 面板内、箱线图下限下方
letter_pos <- function(d, facet_col, y_col, dodge = 0.75) {
  d %>%
    dplyr::group_by(.data[[facet_col]]) %>%
    dplyr::summarise(
      y_letter = min(.data[[y_col]]) - 0.07 * diff(range(.data[[y_col]])),
      .groups = "drop"
    )
}

# ---- Figure 1A：height_change_index ----
posthoc.height %>%
  dplyr::mutate(
    focal = factor(focal, levels = focus_levels),
    neighbor = factor(neighbor, levels = focus_levels),
    signif = trimws(tukey_signif)
  ) -> ph.height

df.idx %>%
  dplyr::group_by(focal) %>%
  dplyr::summarise(y_letter = min(index) - 0.07 * diff(range(index)), .groups = "drop") -> pos.height

ph.height %>%
  dplyr::left_join(pos.height, by = "focal") -> ph.height

df.idx %>%
  ggplot2::ggplot(aes(x = focal, y = index, color = neighbor)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.7) +
  ggplot2::geom_point(
    position = ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
    alpha = 0.5, size = 0.9
  ) +
  ggplot2::geom_text(
    data = ph.height,
    aes(x = focal, y = y_letter, label = signif, group = neighbor),
    color = "black", size = 3.6,
    position = ggplot2::position_dodge(width = 0.7), check_overlap = TRUE
  ) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  ggplot2::geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dashed", color = "grey80") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.05))) +
  bioRtools::scale_color_research() +
  ggplot2::labs(
    x = "Focal accession",
    y = "Plant height change index (cm)",
    color = "Neighbor"
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "right") -> p.1a

bioRtools::write_data(p.1a, file = file.path(out_dir, "06_Figure1A_height_change_index.pdf"),
                      width = 8, height = 6, units = "in")
bioRtools::write_data(p.1a, file = file.path(out_dir, "06_Figure1A_height_change_index.png"),
                      width = 8, height = 6, units = "in", dpi = 500)
cat("[05] 输出：06_Figure1A_height_change_index.pdf / .png\n")

# ---- Figure 1B / 1C：tiller_numberandyield ----
plot_trait <- function(ph, pid, plabel, ylab) {
  d <- bioRtools::read_data(
    file.path(project_root, "data", "processed_inputs", "03_agronomic_traits_height_tiller_yield.xlsx")
  ) %>%
    dplyr::filter(phenotype == ph) %>%
    dplyr::mutate(
      focal = factor(focal, levels = focus_levels),
      neighbor = factor(neighbor, levels = focus_levels)
    )

  posthoc.tiller %>%
    dplyr::filter(phenotype == ph) %>%
    dplyr::mutate(
      focal = factor(focal, levels = focus_levels),
      neighbor = factor(neighbor, levels = focus_levels),
      signif = trimws(tukey_signif)
    ) -> ph_tab

  d %>%
    dplyr::group_by(focal) %>%
    dplyr::summarise(y_letter = min(value) - 0.07 * diff(range(value)), .groups = "drop") %>%
    dplyr::left_join(ph_tab, by = "focal") -> ph_pos

  d %>%
    ggplot2::ggplot(aes(x = focal, y = value, color = neighbor)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.7) +
    ggplot2::geom_point(
      position = ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
      alpha = 0.5, size = 0.9
    ) +
    ggplot2::geom_text(
      data = ph_pos,
      aes(x = focal, y = y_letter, label = signif, group = neighbor),
      color = "black", size = 3.6,
      position = ggplot2::position_dodge(width = 0.7), check_overlap = TRUE
    ) +
    ggplot2::geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dashed", color = "grey80") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.05))) +
    bioRtools::scale_color_research() +
    ggplot2::labs(x = "Focal accession", y = ylab, color = "Neighbor") +
    bioRtools::theme_prism() +
    ggplot2::theme(legend.position = "right") -> p

  bioRtools::write_data(p, file = file.path(out_dir, paste0(pid, "_", plabel, ".pdf")),
                        width = 8, height = 6, units = "in")
  bioRtools::write_data(p, file = file.path(out_dir, paste0(pid, "_", plabel, ".png")),
                        width = 8, height = 6, units = "in", dpi = 500)
  cat("[05] 输出：", pid, "_", plabel, ".pdf / .png\n", sep = "")
}

plot_trait("tiller number", "07", "Figure1B_tiller_number", "Tiller number")
plot_trait("yield", "08", "Figure1C_单株yield", "Grain yield per plant (g)")

# ---- Figure 1D：表型 PCA ----
pc1_lab <- paste0("PC1 (", formatC(100 * var.exp$variance_explained[1], format = "f", digits = 2), "%)")
pc2_lab <- paste0("PC2 (", formatC(100 * var.exp$variance_explained[2], format = "f", digits = 2), "%)")

df.score %>%
  dplyr::mutate(
    focal = factor(focal, levels = focus_levels),
    group = factor(group, levels = c("Pure", "Mix"))
  ) -> df.score

df.score %>%
  dplyr::filter(group == "Pure") %>%
  dplyr::select(focal, center_x = PC1, center_y = PC2) -> df.center

df.score %>%
  dplyr::left_join(df.center, by = "focal") %>%
  ggplot2::ggplot(aes(x = PC1, y = PC2, color = focal)) +
  ggplot2::geom_segment(
    aes(xend = center_x, yend = center_y),
    linetype = "dashed", linewidth = 0.7, alpha = 0.8
  ) +
  ggplot2::geom_point(aes(shape = group), size = 4.5) +
  ggplot2::scale_shape_manual(values = c("Pure" = 16, "Mix" = 17)) +
  bioRtools::scale_color_research() +
  ggplot2::labs(x = pc1_lab, y = pc2_lab, color = "Focal", shape = "Planting mode") +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "right") -> p.1d

bioRtools::write_data(p.1d, file = file.path(out_dir, "09_Figure1D_表型PCA.pdf"),
                      width = 8, height = 6, units = "in")
bioRtools::write_data(p.1d, file = file.path(out_dir, "09_Figure1D_表型PCA.png"),
                      width = 8, height = 6, units = "in", dpi = 500)
cat("[05] 输出：09_Figure1D_表型PCA.pdf / .png\n")

# ---- Figure 1E：三性状 Pure vs Mix 整体comparison ----
bioRtools::read_data(file.path(project_root, "data", "processed_inputs", "03_agronomic_traits_height_tiller_yield.xlsx")) %>%
  dplyr::filter(phenotype %in% c("tiller number", "yield")) %>%
  dplyr::select(focal, neighbor, group, value, phenotype) -> df.other

df.idx %>%
  dplyr::select(focal, neighbor, group, value = index, phenotype) %>%
  dplyr::bind_rows(df.other) %>%
  dplyr::mutate(
    phenotype = factor(phenotype,
                       levels = c("plant height", "tiller number", "yield"),
                       labels = c("Plant height index", "Tiller number", "Yield per plant")),
    group = factor(group, levels = c("Pure", "Mix"))
  ) %>%
  dplyr::group_by(phenotype) %>%
  rstatix::t_test(value ~ group) %>%
  dplyr::ungroup() %>%
  dplyr::select(phenotype, p) -> p_tab

cat("\n[05] Figure 1E 的 t test：\n")
print(as.data.frame(p_tab))

df.idx %>%
  dplyr::select(focal, neighbor, group, value = index, phenotype) %>%
  dplyr::bind_rows(df.other) %>%
  dplyr::mutate(
    phenotype = factor(phenotype,
                       levels = c("plant height", "tiller number", "yield"),
                       labels = c("Plant height index", "Tiller number", "Yield per plant")),
    group = factor(group, levels = c("Pure", "Mix"))
  ) %>%
  dplyr::left_join(p_tab, by = "phenotype") %>%
  dplyr::group_by(phenotype) %>%
  dplyr::mutate(y_pos = max(value) + 0.08 * diff(range(value))) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot(aes(x = group, y = value)) +
  ggplot2::geom_violin(aes(fill = phenotype), alpha = 0.5, color = NA, show.legend = FALSE) +
  ggplot2::geom_boxplot(aes(fill = phenotype), width = 0.3, outlier.shape = NA,
                        alpha = 0.8, show.legend = FALSE) +
  ggplot2::geom_jitter(width = 0.12, alpha = 0.25, size = 0.6, color = "grey30") +
  ggplot2::geom_text(aes(x = 1.5, y = y_pos, label = format_p(p)), size = 4, check_overlap = TRUE) +
  ggplot2::facet_wrap(~phenotype, scales = "free") +
  ggplot2::scale_x_discrete(limits = c("Pure", "Mix")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
  bioRtools::scale_fill_research() +
  ggplot2::labs(x = NULL, y = "Value") +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(face = "bold")) -> p.1e

bioRtools::write_data(p.1e, file = file.path(out_dir, "10_Figure1E_三性状整体comparison.pdf"),
                      width = 10, height = 5, units = "in")
bioRtools::write_data(p.1e, file = file.path(out_dir, "10_Figure1E_三性状整体comparison.png"),
                      width = 10, height = 5, units = "in", dpi = 500)
cat("[05] 输出：10_Figure1E_三性状整体comparison.pdf / .png\n")

# ---- Figure 1F：遗传距离and三性状（Mix only）----
df.join %>%
  dplyr::filter(group == "Mix") %>%
  tidyr::pivot_longer(cols = c(mean_height, mean_tiller, mean_yield),
                      names_to = "trait", values_to = "value") %>%
  dplyr::mutate(trait = factor(trait,
                               levels = c("mean_height", "mean_tiller", "mean_yield"),
                               labels = c("Plant height index", "Tiller number", "Yield per plant"))) -> df.f

cor.tab %>%
  dplyr::filter(dataset == "Mix only (n = 12)") %>%
  dplyr::mutate(
    trait = factor(trait,
                   levels = c("mean_height", "mean_tiller", "mean_yield"),
                   labels = c("Plant height index", "Tiller number", "Yield per plant")),
    label = paste0("italic(r) == ", formatC(r, format = "f", digits = 2), "~~", format_p(p_value))
  ) -> cor.f

df.f %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(y_pos = max(value) - 0.05 * diff(range(value))) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot(aes(x = genetic_distance, y = value)) +
  ggplot2::geom_point(size = 2.6, color = "#1e4e8c") +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "grey30", fill = "grey85") +
  ggplot2::geom_text(
    data = cor.f,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.08, vjust = 1.6, size = 3.6, parse = TRUE, inherit.aes = FALSE
  ) +
  ggplot2::facet_wrap(~trait, scales = "free_y") +
  ggplot2::labs(x = "Genetic distance", y = "Value") +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold")) -> p.1f

bioRtools::write_data(p.1f, file = file.path(out_dir, "11_Figure1F_遗传距离and表型.pdf"),
                      width = 10, height = 4.5, units = "in")
bioRtools::write_data(p.1f, file = file.path(out_dir, "11_Figure1F_遗传距离and表型.png"),
                      width = 10, height = 4.5, units = "in", dpi = 500)
cat("[05] 输出：11_Figure1F_遗传距离and表型.pdf / .png\n")

cat("\n[05] sessionInfo()\n")
sessionInfo()
