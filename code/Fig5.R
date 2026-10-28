# ============================================================
# Fig5.R —— Figure 5（Interaction architecture: pathogen-specific and environment-dependent）
# Source Data：只读 source_data/Source_Data.xlsx（sheet：Fig5a–Fig5c）
# 运行：cd manuscript && Rscript code/Fig5.R
# 输出：figures/Fig5.pdf / .png；面板 figures/panels/Fig5a–c.pdf
# 软件：R 4.6.1；ggplot2；bioRtools::theme_prism()；set.seed(12345)
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/Fig5.R
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/_vectorcraft.R")
sd_file <- "source_data/Source_Data.xlsx"
rd <- function(sheet) readxl::read_excel(sd_file, sheet = sheet, skip = 2)
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)
save_panel <- function(p, name, w_cm, h_cm) {
  bioRtools::write_data(p, file = file.path("figures/panels", paste0(name, ".pdf")),
                        width = w_cm, height = h_cm, units = "cm")
}
fmt_p <- function(p, digits = 2) dplyr::if_else(p < 0.001, "P<0.001", paste0("P=", formatC(p, format = "f", digits = digits)))

# 5A：结构分解（F 值，按病原分面）
d5a <- rd("Fig5a") %>%
  dplyr::mutate(term = factor(term, levels = c("focal", "neighbor", "focal:neighbor"),
                              labels = c("Focal identity", "Neighbor identity",
                                         "Focal x neighbor (specific pairing)")),
                disease = factor(disease, levels = c("blast", "brown spot"),
                                 labels = c("Rice blast", "Brown spot")),
                label = paste0("F = ", formatC(f_value, format = "f", digits = 1), "\n",
                               dplyr::if_else(p_value < 0.001, "P<0.001",
                                              paste0("P=", formatC(p_value, format = "f", digits = 2)))))
p5a <- d5a %>%
  ggplot2::ggplot(aes(x = term, y = f_value, fill = term)) +
  ggplot2::geom_col(width = 0.68, show.legend = FALSE) +
  ggplot2::geom_text(aes(label = label), vjust = -0.25, size = 2.9, lineheight = 0.95) +
  ggplot2::facet_wrap(~disease) +
  bioRtools::scale_fill_research() +
  ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.28))) +
  ggplot2::labs(x = NULL, y = "F value (batch as random effect)",
                title = "Interaction structure differs between pathogens",
                subtitle = "Mix combinations form an incomplete 4 x 4 design; interaction df = 5") +
  bioRtools::theme_prism() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1, vjust = 1, size = 8),
                 strip.text = ggplot2::element_text(face = "bold"),
                 plot.title = ggplot2::element_text(size = 8, face = "bold"),
                 plot.subtitle = ggplot2::element_text(size = 6.5))
save_panel(p5a, "Fig5a", 17.0, 7.0)

# 5B：reciprocity（focal 效应 vs 伙伴效应）
d5b <- rd("Fig5b") %>%
  dplyr::mutate(disease = factor(disease, levels = c("blast", "brown spot"),
                                 labels = c("Rice blast", "Brown spot"))) %>%
  dplyr::mutate(lab = paste0("symmetry permutation ", fmt_p(perm_p_symmetry, 3)))
p5b <- d5b %>%
  ggplot2::ggplot(aes(x = effect_focal, y = effect_neighbor)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey45") +
  ggplot2::geom_point(size = 2.6, color = "#1e4e8c") +
  ggrepel::geom_text_repel(aes(label = pair), size = 2.7, max.overlaps = 20, seed = 12345) +
  # 面板内不再重复标注置换检验 P（会与点标签互压）；该信息见正文 Figure legend
  ggplot2::facet_wrap(~disease, scales = "free") +
  # 2026-10-19：旋转 y 轴标题比绘图区还长、上下两头被切 → 缩短 + 缩小字号 + 留边距（单位写进图注）
  ggplot2::labs(x = "Effect on focal (Z)",
                y = "Effect on the partner (Z)",
                title = "Interaction effects are asymmetric between partners") +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"),
                 axis.title = ggplot2::element_text(size = 6.5),
                 plot.margin = ggplot2::margin(3, 2, 3, 2),
                 plot.title = ggplot2::element_text(size = 8, face = "bold"))
save_panel(p5b, "Fig5b", 17.0, 8.0)

# 5C：逐年效应
d5c <- rd("Fig5c") %>%
  dplyr::mutate(category = factor(category, levels = c("Blast", "BrownSpot", "Smut", "PH", "KGW")),
                year = factor(year)) %>%
  # 2026-10-19：与投稿面板一致——P 标注改两行（原单行比所在分面还宽，被面板 clip 切掉尾部）
  dplyr::mutate(lab = paste0("condition x Year\n", fmt_p(p_interaction, 2)))
ref_lines <- tibble::tibble(category = factor(c("Blast", "BrownSpot", "Smut", "PH", "KGW"),
                                              levels = c("Blast", "BrownSpot", "Smut", "PH", "KGW")),
                            ref = c(1, 1, 1, 0, 0))
p5c <- d5c %>%
  ggplot2::ggplot(aes(x = estimate, y = year)) +
  ggplot2::geom_vline(data = ref_lines, aes(xintercept = ref), linetype = "dashed", color = "grey40") +
  ggplot2::geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, color = "grey35") +
  ggplot2::geom_point(aes(color = year), size = 2.6) +
  ggplot2::geom_text(aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.4,
                     size = 2.2, lineheight = 0.95) +
  ggplot2::facet_wrap(~category, scales = "free_x", nrow = 1) +
  ggplot2::scale_y_discrete(limits = rev(levels(d5c$year))) +
  bioRtools::scale_color_research() +
  ggplot2::labs(x = "Effect estimate per year (diseases: odds ratio; traits: mean difference)", y = "Year",
                title = "Field effects are not stable across years") +
  ggplot2::coord_cartesian(clip = "off") +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 8),
                 strip.clip = "off",
                 legend.position = "none", plot.title = ggplot2::element_text(size = 8, face = "bold"),
                 axis.text.x = ggplot2::element_text(size = 7))
save_panel(p5c, "Fig5c", 17.0, 6.4)

compose_vectorcraft(
  list(
    list(path = "figures/panels/Fig5a.pdf", rect = c(2, 8, 170, 70), label = "A", label_pos = c(2, 3.5)),
    list(path = "figures/panels/Fig5b.pdf", rect = c(2, 83, 170, 80), label = "B", label_pos = c(2, 78.5)),
    list(path = "figures/panels/Fig5c.pdf", rect = c(2, 168, 170, 60), label = "C", label_pos = c(2, 163.5))
  ),
  width_mm = 174, height_mm = 230, out_pdf = "figures/Fig5.pdf")

cat("\n[Fig5] sessionInfo()\n")
sessionInfo()
