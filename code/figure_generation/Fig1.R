# ============================================================
# Fig1_v32.R —— Figure 1（second-round revision）
# Source Data：**只读** source_data/Source_Data.xlsx 的 sheet Fig1a–Fig1e
# 运行：cd manuscript && Rscript code/Fig1.R
# 输出：figures/Fig1.pdf（+ figures/Fig1.png）；面板 figures/panels/Fig1a–e.pdf
# 软件：R 4.6.1；ggplot2；bioRtools::theme_prism()；随机种子 set.seed(12345)
# 组图：由 code/_vectorcraft.R 调 VectorCraft CLI 完成（版面 rect 见文件末）
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/Fig1.R
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/_vectorcraft.R")

sd_file <- "source_data/Source_Data_v32.xlsx"
rd <- function(sheet) readxl::read_excel(sd_file, sheet = sheet, skip = 2)
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)
save_panel <- function(p, name, w_cm, h_cm) {
  bioRtools::write_data(p, file = file.path("figures/panels", paste0(name, ".pdf")),
                        width = w_cm, height = h_cm, units = "cm")
}

fmt_p <- function(p, digits = 2) dplyr::if_else(p < 0.001, "P<0.001", paste0("P=", formatC(p, format = "f", digits = digits)))

d1a <- rd("Fig1a") %>%
  dplyr::mutate(focal = factor(as.character(focal), levels = c("15", "49", "69", "87")),
                neighbor = factor(as.character(neighbor), levels = c("15", "49", "69", "87")))

p1a <- d1a %>%
  ggplot2::ggplot(aes(x = neighbor, y = index_mean)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  ggplot2::geom_col(fill = "grey88", width = 0.62) +
  ggplot2::geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.18, color = "grey25") +
  ggplot2::facet_wrap(~focal, nrow = 1, scales = "free_y") +
  ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 4),
                              expand = ggplot2::expansion(mult = c(0.12, 0.30))) +
  ggplot2::labs(x = "Neighbour accession", y = "Height index (cm)") +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 8),
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
save_panel(p1a, "Fig1a", 17.0, 5.0)

panel_trait <- function(sheet, ylab, name) {
  d <- rd(sheet) %>%
    dplyr::mutate(focal = factor(as.character(focal), levels = c("15", "49", "69", "87")),
                  neighbor = factor(as.character(neighbor), levels = c("15", "49", "69", "87")),
                  group = factor(group, levels = c("Pure", "Mix")))
  p <- d %>%
    ggplot2::ggplot(aes(x = neighbor, y = value)) +
    ggplot2::geom_boxplot(aes(fill = group), width = 0.62, outlier.shape = NA, alpha = 0.85) +
    ggplot2::geom_jitter(aes(color = group), width = 0.16, alpha = 0.45, size = 1) +
    ggplot2::facet_wrap(~focal, nrow = 1, scales = "free_y") +
    ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 3),
                                expand = ggplot2::expansion(mult = c(0.05, 0.22))) +
    ggplot2::scale_fill_manual(values = c(Pure = "grey85", Mix = "#1e4e8c"), name = NULL) +
    ggplot2::scale_color_manual(values = c(Pure = "grey45", Mix = "#1e4e8c"), name = NULL) +
    ggplot2::labs(x = "Neighbour accession", y = ylab) +
    bioRtools::theme_prism() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"),
                   legend.position = "top",
                   legend.key.size = ggplot2::unit(0.35, "cm"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   axis.text.y = ggplot2::element_text(size = 6.5))
  save_panel(p, name, 17.0, 4.6)
}
panel_trait("Fig1b", "Tiller number per plant", "Fig1b")
panel_trait("Fig1c", "Yield per plant (g)", "Fig1c")

d1d <- rd("Fig1d") %>%
  dplyr::mutate(focal = factor(as.character(focal), levels = c("15", "49", "69", "87")),
                group = factor(group, levels = c("Pure", "Mix"))) %>%
  dplyr::group_by(focal) %>%
  dplyr::mutate(center_x = mean(PC1), center_y = mean(PC2)) %>%
  dplyr::ungroup()
p1d <- d1d %>%
  ggplot2::ggplot(aes(x = PC1, y = PC2)) +
  ggplot2::geom_segment(aes(xend = center_x, yend = center_y, color = focal),
                        linetype = "dashed", linewidth = 0.6, alpha = 0.75) +
  ggplot2::geom_point(aes(color = focal, shape = group), size = 3.2) +
  ggplot2::scale_shape_manual(values = c(Pure = 16, Mix = 17), name = "Planting mode") +
  bioRtools::scale_color_research() +
  ggplot2::labs(x = "PC1 (49.31%)", y = "PC2 (38.85%)", color = "Focal accession") +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "right",
                 legend.key.size = ggplot2::unit(0.32, "cm"),
                 legend.text = ggplot2::element_text(size = 7),
                 legend.title = ggplot2::element_text(size = 7.5))
save_panel(p1d, "Fig1d", 8.3, 5.2)

d1e <- rd("Fig1e") %>%
  dplyr::mutate(group = factor(group, levels = c("Pure", "Mix")),
                trait = factor("Final plant height (cm)",
                               levels = "Final plant height (cm)"))
p1e <- d1e %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(y_pos = max(value, na.rm = TRUE) + 0.08 * diff(range(value, na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot(aes(x = group, y = value)) +
  ggplot2::geom_violin(aes(fill = trait), alpha = 0.45, color = NA) +
  ggplot2::geom_boxplot(aes(fill = trait), width = 0.28, outlier.shape = NA, alpha = 0.85) +
  ggplot2::geom_jitter(width = 0.12, alpha = 0.25, size = 0.6, color = "grey30") +
  ggplot2::facet_wrap(~trait, scales = "free_y") +
  ggplot2::scale_x_discrete(limits = c("Pure", "Mix")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.16))) +
  bioRtools::scale_fill_research() +
  ggplot2::labs(x = NULL, y = "Final plant height (cm)") +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 6.5),
                 strip.clip = "off",
                 legend.position = "none")
save_panel(p1e, "Fig1e", 8.3, 5.2)

compose_vectorcraft(
  list(
    list(path = "figures/panels/Fig1a.pdf", rect = c(2, 8, 170, 50), label = "A", label_pos = c(2, 3.5)),
    list(path = "figures/panels/Fig1b.pdf", rect = c(2, 63, 170, 46), label = "B", label_pos = c(2, 58.5)),
    list(path = "figures/panels/Fig1c.pdf", rect = c(2, 114, 170, 46), label = "C", label_pos = c(2, 109.5)),
    list(path = "figures/panels/Fig1d.pdf", rect = c(2, 165, 83, 52), label = "D", label_pos = c(2, 160.5)),
    list(path = "figures/panels/Fig1e.pdf", rect = c(91, 165, 83, 52), label = "E", label_pos = c(91, 160.5))
  ),
  width_mm = 174, height_mm = 219, out_pdf = "figures/Fig1_v32.pdf")

cat("\n[Fig1] sessionInfo()\n")
sessionInfo()
