# ============================================================
# FigS2.R —— Figure S2（Greenhouse inoculation results and field disease survey；8 个子图 a–h）
# Source Data：只读 source_data/Source_Data.xlsx（sheet：FigS2a–FigS2h）
# 运行：cd manuscript && Rscript code/FigS2.R
# 输出：figures/FigS2.pdf / .png；面板 figures/panels/FigS2a–h.pdf
# 软件：R 4.6.1；ggplot2；bioRtools::theme_prism()；set.seed(12345)
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/FigS2.R
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

inoc_box <- function(d, col, ylab) {
  d %>% dplyr::mutate(group = factor(group, levels = c("Pure", "Mix"))) %>%
    ggplot2::ggplot(aes(x = group, y = .data[[col]])) +
    ggplot2::geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.6, width = 0.6) +
    ggplot2::geom_jitter(aes(color = batch), width = 0.2, alpha = 0.8, size = 1.2) +
    ggplot2::facet_wrap(~disease, scales = "free_y") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.05))) +
    bioRtools::scale_fill_research() + bioRtools::scale_color_research() +
    ggplot2::labs(x = NULL, y = ylab) +
    bioRtools::theme_prism() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8, face = "bold"),
                   strip.clip = "off",   # 2026-10-19：防止比 strip 宽的标题被两侧切掉
                   legend.position = "right", legend.key.size = ggplot2::unit(0.3, "cm"),
                   legend.text = ggplot2::element_text(size = 6))
}
d2a <- rd("FigS2a")
# 2026-10-19：y 轴单位用真上标 cm^2（原为纯文本 "cm2"）
save_panel(inoc_box(d2a, "lesion_density_per_cm2",
                    expression("Lesion number per cm"^2)), "FigS2a", 8.3, 6.0)
d2b <- rd("FigS2b")
save_panel(inoc_box(d2b, "lesion_density_z", "Relative lesion density (Z-score)"), "FigS2b", 8.3, 6.0)

for (nm in c("FigS2c", "FigS2d", "FigS2e", "FigS2f")) {
  d <- rd(nm) %>% dplyr::mutate(group = factor(group, levels = c("Pure", "Mix")),
                                neighbor = factor(as.character(neighbor), levels = c("15", "49", "69", "87")))
  l <- d %>% dplyr::group_by(disease, neighbor) %>%
    dplyr::summarise(signif = dplyr::first(signif), y_letter = dplyr::first(y_letter), .groups = "drop")
  p <- d %>%
    ggplot2::ggplot(aes(x = neighbor, y = lesion_norm)) +
    ggplot2::geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.6, width = 0.6) +
    ggplot2::geom_jitter(aes(color = batch), width = 0.2, alpha = 0.8, size = 1.2) +
    ggplot2::geom_text(data = l, aes(x = neighbor, y = y_letter, label = signif),
                       inherit.aes = FALSE, size = 3) +
    ggplot2::facet_wrap(~disease, scales = "free_y") +
    # 2026-10-19：与投稿面板一致——图例置底（原右侧图例把两分面各挤到约 50 pt 宽），
    #   分面标题 8 → 6.5 pt 并关掉 strip.clip（默认 "on" 会把比 strip 宽的标题两侧切掉）
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
                    color = ggplot2::guide_legend(nrow = 1, byrow = TRUE)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.05))) +
    bioRtools::scale_fill_research() + bioRtools::scale_color_research() +
    # 2026-10-19：旋转 y 轴标题两头被切 → 缩短 + 缩小字号 + 留边距；图例改上下两行堆叠（并排会超宽）
    ggplot2::labs(x = "Neighbour accession", y = "Relative lesion density (Z)") +
    bioRtools::theme_prism() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 6.5, face = "bold"),
                   strip.clip = "off",
                   legend.position = "bottom", legend.box = "vertical",
                   legend.key.size = ggplot2::unit(0.26, "cm"),
                   legend.text = ggplot2::element_text(size = 5.8),
                   legend.title = ggplot2::element_text(size = 6),
                   legend.box.spacing = ggplot2::unit(0.08, "cm"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   axis.text.x = ggplot2::element_text(size = 6),
                   axis.title.y = ggplot2::element_text(size = 6.5),
                   plot.margin = ggplot2::margin(3, 2, 2, 2))
  save_panel(p, nm, 8.3, 6.0)   # 2026-10-19：保持 6.0 cm（S2 整页已 365 mm，不再加高）
}

field_bar <- function(d, name) {
  d <- d %>% dplyr::mutate(focal = factor(as.character(focal)),
                           condition = factor(condition, levels = c("Pure", "Mix")))
  p <- d %>%
    ggplot2::ggplot(aes(x = focal, y = prob, fill = condition)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.62) +
    ggplot2::geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                           position = ggplot2::position_dodge(width = 0.7), width = 0.15) +
    ggplot2::facet_wrap(~disease, scales = "free_y") +
    bioRtools::scale_fill_research() +
    # 2026-10-19：旋转 y 轴标题被切 → 缩短 + 缩小字号 + 留边距
    ggplot2::labs(x = "Focal accession", y = "Disease incidence (prob.)", fill = NULL) +
    bioRtools::theme_prism() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8, face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 6.5),
                   plot.margin = ggplot2::margin(3, 2, 2, 2),
                   legend.position = "bottom", legend.key.size = ggplot2::unit(0.3, "cm"))
  save_panel(p, name, 17.0, 7.0)
}
field_bar(rd("FigS2g"), "FigS2g")
save_panel({
  d <- rd("FigS2h") %>% dplyr::mutate(focal = factor(as.character(focal)),
                                      condition = factor(condition, levels = c("Pure", "Mix")))
  d %>% ggplot2::ggplot(aes(x = focal, y = prob, fill = condition)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.62) +
    ggplot2::geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                           position = ggplot2::position_dodge(width = 0.7), width = 0.15) +
    ggplot2::facet_wrap(~disease, scales = "free_y", nrow = 1) +
    bioRtools::scale_fill_research() +
    ggplot2::labs(x = "Focal accession", y = "Disease incidence (probability)", fill = NULL) +
    bioRtools::theme_prism() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8, face = "bold"),
                   legend.position = "bottom", legend.key.size = ggplot2::unit(0.3, "cm"))
}, "FigS2h", 17.0, 8.5)

compose_vectorcraft(
  list(
    list(path = "figures/panels/FigS2a.pdf", rect = c(2, 8, 83, 60), label = "A", label_pos = c(2, 3.5)),
    list(path = "figures/panels/FigS2b.pdf", rect = c(91, 8, 83, 60), label = "B", label_pos = c(91, 3.5)),
    list(path = "figures/panels/FigS2c.pdf", rect = c(2, 73, 83, 60), label = "C", label_pos = c(2, 68.5)),
    list(path = "figures/panels/FigS2d.pdf", rect = c(91, 73, 83, 60), label = "D", label_pos = c(91, 68.5)),
    list(path = "figures/panels/FigS2e.pdf", rect = c(2, 138, 83, 60), label = "E", label_pos = c(2, 133.5)),
    list(path = "figures/panels/FigS2f.pdf", rect = c(91, 138, 83, 60), label = "F", label_pos = c(91, 133.5)),
    list(path = "figures/panels/FigS2g.pdf", rect = c(2, 203, 170, 70), label = "G", label_pos = c(2, 198.5)),
    list(path = "figures/panels/FigS2h.pdf", rect = c(2, 278, 170, 85), label = "H", label_pos = c(2, 273.5))
  ),
  width_mm = 174, height_mm = 365, out_pdf = "figures/FigS2.pdf")

cat("\n[FigS2] sessionInfo()\n")
sessionInfo()
