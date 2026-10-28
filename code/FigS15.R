# ============================================================
# FigS15.R —— Figure S15（Stage specificity of the transcriptional response；2 个子图 a、b）
# Source Data：只读 source_data/Source_Data.xlsx（sheet：FigS15a、FigS15b）
# 运行：cd manuscript && Rscript code/FigS15.R
# 输出：figures/FigS15.pdf / .png；面板 figures/panels/FigS15a–b.pdf
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/FigS15.R
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
stage_levels <- c("tillering", "flowering", "filling")

d15a <- rd("FigS15a") %>% dplyr::mutate(stage = factor(stage, levels = stage_levels))
top20 <- d15a %>% dplyr::group_by(gene) %>% dplyr::summarise(m = max(abs(logFC)), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(m)) %>% dplyr::slice_head(n = 20) %>% dplyr::pull(gene)
p15a <- d15a %>%
  dplyr::filter(gene %in% top20) %>%
  ggplot2::ggplot(aes(x = stage, y = logFC, group = gene)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_line(alpha = 0.55, color = "grey35") +
  ggplot2::geom_point(size = 2, color = "#1e4e8c") +
  # 2026-10-19：旋转 y 轴标题两头被切 → 缩短 + 缩小字号 + 上下留边距
  ggplot2::labs(x = NULL, y = "Moderated log2FC (Mix vs Pure)",
                title = "Largest-magnitude genes across stages") +
  bioRtools::theme_prism() +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold"),
                 axis.title.y = ggplot2::element_text(size = 6.5),
                 plot.margin = ggplot2::margin(3, 2, 3, 2),
                 axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1))
save_panel(p15a, "FigS15a", 8.3, 6.5)

d15b <- rd("FigS15b")
p15b <- d15b %>%
  dplyr::mutate(contrast = factor(contrast, levels = rev(contrast))) %>%
  tidyr::pivot_longer(cols = c(n_p005, n_fdr005), names_to = "criterion", values_to = "n") %>%
  dplyr::mutate(criterion = dplyr::recode(criterion, n_p005 = "P<0.05", n_fdr005 = "FDR<0.05")) %>%
  ggplot2::ggplot(aes(x = n, y = contrast, fill = criterion)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
  # 2026-10-19：与投稿面板一致——柱尾数字超出面板右界被 clip 切掉（"24" 缺 6.5 pt），
  #   字号 3 → 2.6、给 x 轴留右扩展、并关掉面板 clip
  ggplot2::geom_text(aes(label = dplyr::if_else(n == 0, "", as.character(n))),
                     position = ggplot2::position_dodge(width = 0.7), hjust = -0.25, size = 2.6) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.14))) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::scale_fill_manual(values = c("P<0.05" = "#e0a458", "FDR<0.05" = "#a63446"), name = NULL) +
  ggplot2::labs(x = "Number of genes", y = NULL, title = "Response candidates per contrast") +
  bioRtools::theme_prism() +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold"),
                 legend.position = "bottom", legend.key.size = ggplot2::unit(0.3, "cm"),
                 legend.text = ggplot2::element_text(size = 7))
save_panel(p15b, "FigS15b", 8.3, 6.0)

compose_vectorcraft(
  list(
    list(path = "figures/panels/FigS15a.pdf", rect = c(2, 8, 83, 65), label = "A", label_pos = c(2, 3.5)),
    list(path = "figures/panels/FigS15b.pdf", rect = c(91, 8, 83, 60), label = "B", label_pos = c(91, 3.5))
  ),
  width_mm = 174, height_mm = 78, out_pdf = "figures/FigS15.pdf")

cat("\n[FigS15] sessionInfo()\n")
sessionInfo()
