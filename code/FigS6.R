# ============================================================
# FigS6.R —— Figure S6（Detection summary of the root-exudate hormone panel；2 个子图 a、b）
# Source Data：只读 source_data/Source_Data.xlsx（sheet：FigS6a、FigS6b）
# 运行：cd manuscript && Rscript code/FigS6.R
# 输出：figures/FigS6.pdf / .png；面板 figures/panels/FigS6a–b.pdf
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/FigS6.R
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

d6a <- rd("FigS6a") %>% dplyr::mutate(class_en = factor(class_en, levels = class_en))
p6a <- d6a %>%
  ggplot2::ggplot(aes(x = class_en, y = detect_rate)) +
  ggplot2::geom_col(fill = "#1e4e8c", alpha = 0.85, width = 0.65) +
  # 2026-10-19：与投稿面板一致——标签由 "n/N compounds" 缩短为 "n/N"（原标签比柱尾到轴上界
  #   剩下的空间还长，被面板 clip 切掉尾部），字号 2.9 → 2.6，coord_flip 里关掉 clip
  ggplot2::geom_text(aes(label = paste0(n_compound_detected, "/", n_compound)),
                     hjust = -0.08, size = 2.6) +
  ggplot2::coord_flip(clip = "off") +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                              expand = ggplot2::expansion(mult = c(0.02, 0.30))) +
  ggplot2::labs(x = NULL, y = "Detection rate (fraction of 12 samples)",
                title = "Hormone detection rate by compound class",
                # 2026-10-19：subtitle 加长后越出纸边 2 pt → 回到短版（标签含义写进正文图注）
                subtitle = "109 compounds x 12 samples; no brassinosteroid in this panel") +
  bioRtools::theme_prism() +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold"),
                 plot.subtitle = ggplot2::element_text(size = 6.5))
save_panel(p6a, "FigS6a", 10.4, 10.0)

d6b <- rd("FigS6b") %>% dplyr::mutate(index = factor(index, levels = rev(index)),
                                      class_en = factor(class_en))
p6b <- d6b %>%
  ggplot2::ggplot(aes(x = index, y = detect_rate, fill = class_en)) +
  ggplot2::geom_col(width = 0.75) +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                              expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
  bioRtools::scale_fill_research() +
  ggplot2::labs(x = NULL, y = "Detection rate", fill = "Class",
                title = "38 of 109 compounds detected") +
  bioRtools::theme_prism() +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 4.5),
                 plot.title = ggplot2::element_text(size = 8, face = "bold"),
                 legend.key.size = ggplot2::unit(0.3, "cm"),
                 legend.text = ggplot2::element_text(size = 6))
save_panel(p6b, "FigS6b", 6.6, 10.0)

compose_vectorcraft(
  list(
    list(path = "figures/panels/FigS6a.pdf", rect = c(2, 8, 104, 100), label = "A", label_pos = c(2, 3.5)),
    list(path = "figures/panels/FigS6b.pdf", rect = c(108, 8, 66, 100), label = "B", label_pos = c(108, 3.5))
  ),
  width_mm = 174, height_mm = 110, out_pdf = "figures/FigS6.pdf")

cat("\n[FigS6] sessionInfo()\n")
sessionInfo()
