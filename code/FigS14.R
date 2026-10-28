# ============================================================
# FigS14.R —— Figure S14（Year-to-year stability of field effects；单面板）
# Source Data：只读 source_data/Source_Data.xlsx（sheet：FigS14）
# 运行：cd manuscript && Rscript code/FigS14.R
# 输出：figures/FigS14.pdf / .png
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/FigS14.R
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/_vectorcraft.R")
sd_file <- "source_data/Source_Data.xlsx"
d <- readxl::read_excel(sd_file, sheet = "FigS14", skip = 2) %>%
  dplyr::mutate(category = factor(category, levels = unique(category)), year = factor(year),
                lab = paste0("condition x Year ",
                             dplyr::if_else(p_interaction < 0.001, "P<0.001",
                                            paste0("P=", formatC(p_interaction, format = "f", digits = 2)))))
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)

ref_lines <- d %>% dplyr::distinct(category, scale) %>%
  dplyr::mutate(ref = dplyr::if_else(scale == "odds_ratio", 1, 0))
p <- d %>%
  ggplot2::ggplot(aes(x = estimate, y = year)) +
  ggplot2::geom_vline(data = ref_lines, aes(xintercept = ref), linetype = "dashed", color = "grey40") +
  ggplot2::geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, color = "grey35") +
  ggplot2::geom_point(aes(color = year), size = 2.6) +
  ggplot2::geom_text(aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.6, size = 2.4) +
  ggplot2::facet_wrap(~category, scales = "free", nrow = 2) +
  bioRtools::scale_color_research() +
  ggplot2::labs(x = "Effect estimate per year (diseases: odds ratio; traits: mean difference)",
                y = "Year", color = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 8),
                 legend.position = "none", axis.text.x = ggplot2::element_text(size = 7))
bioRtools::write_data(p, file = "figures/panels/FigS14.pdf", width = 17.4, height = 8.7, units = "cm")

compose_vectorcraft(
  list(list(path = "figures/panels/FigS14.pdf", rect = c(2, 8, 170, 87), label = "", label_pos = c(2, 3.5))),
  width_mm = 174, height_mm = 95, out_pdf = "figures/FigS14.pdf")

cat("\n[FigS14] sessionInfo()\n")
sessionInfo()
