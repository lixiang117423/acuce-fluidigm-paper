# ============================================================
# FigS16.R —— Figure S16（Relationship between genetic distance and agronomic traits；单面板）
# Source Data：只读 data/source_data/Source_Data_v45.xlsx（sheet：FigS16；图用 12 个混种组合）
# 运行：Rscript code/figure_generation/FigS16.R
# 输出：figures/FigS16.pdf / .png
# 本机绝对路径：<repository_root>/manuscript/code/FigS16.R
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/figure_generation/_vectorcraft.R")
sd_file <- "data/source_data/Source_Data_v45.xlsx"
d <- readxl::read_excel(sd_file, sheet = "FigS16", skip = 2)
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)

d_mix <- d %>% dplyr::filter(group == "Mix")
plot_d <- d_mix %>%
  dplyr::select(focal, neighbor, genetic_distance, mean_height, mean_tiller, mean_yield) %>%
  tidyr::pivot_longer(cols = c(mean_height, mean_tiller, mean_yield),
                      names_to = "trait", values_to = "value") %>%
  dplyr::mutate(trait = factor(trait, levels = c("mean_height", "mean_tiller", "mean_yield"),
                               labels = c("Plant height index (cm)", "Tiller number (per plant)",
                                          "Yield (g per plant)")),
                lab = dplyr::case_when(
                  trait == "Plant height index (cm)" ~ paste0("r = ", formatC(d$r_height[1], format = "f", digits = 2),
                                                              ", P = ", formatC(d$p_height[1], format = "f", digits = 2)),
                  trait == "Tiller number (per plant)" ~ paste0("r = ", formatC(d$r_tiller[1], format = "f", digits = 2),
                                                                ", P = ", formatC(d$p_tiller[1], format = "f", digits = 2)),
                  TRUE ~ paste0("r = ", formatC(d$r_yield[1], format = "f", digits = 2),
                                ", P = ", formatC(d$p_yield[1], format = "f", digits = 2))))

p <- plot_d %>%
  ggplot2::ggplot(aes(x = genetic_distance, y = value)) +
  ggplot2::geom_point(size = 2.6, color = "#1e4e8c") +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "grey30", fill = "grey85") +
  ggplot2::geom_text(data = plot_d %>% dplyr::distinct(trait, lab),
                     aes(x = -Inf, y = Inf, label = lab), hjust = -0.08, vjust = 1.6,
                     size = 3, inherit.aes = FALSE) +
  ggplot2::facet_wrap(~trait, scales = "free_y") +
  ggplot2::labs(x = "Pairwise genetic distance", y = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 8))
bioRtools::write_data(p, file = "figures/panels/FigS16.pdf", width = 17.0, height = 6.95, units = "cm")

compose_vectorcraft(
  list(list(path = "figures/panels/FigS16.pdf", rect = c(2, 8, 170, 69.5), label = "", label_pos = c(2, 3.5))),
  width_mm = 174, height_mm = 77, out_pdf = "figures/FigS16.pdf")

cat("\n[FigS16] sessionInfo()\n")
sessionInfo()
