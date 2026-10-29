# ============================================================
# 01_inoculation_and_field_disease / 05_field_disease_plotting.R
# 作用：绘制field_disease调查图
#       Figure S2G：焦点纯系 15 / 49 / 69 / 87 × Blast / BrownSpot（主analysis，正文引用）
#       Figure S2H：全部 6 个纯系 × Blast / BrownSpot / Smut（补充，含 36 / 74 and Smut）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/01_inoculation_and_field_disease/05_field_disease_plotting.R"
#
# 输入：<root>/result/01_inoculation_and_field_disease/04_表_field_disease_main_GLMM.xlsx（sheet emmean_letters_main）
#       <root>/result/01_inoculation_and_field_disease/05_表_field_disease_supplement_Smutand36_74.xlsx（sheet emmean_letters_all）
# 输出：<root>/result/01_inoculation_and_field_disease/13_FigureS2G_field_disease调查.pdf / .png
#       <root>/result/01_inoculation_and_field_disease/14_FigureS2H_field_disease调查_全部纯系andSmut.pdf / .png
# 其中 <root> = <repository_root>
#
# plotting约定：bioRtools::theme_prism()；bioRtools::scale_color_research()；
#          bioRtools::write_data() 输出，尺寸单位显式 units = "in"。
# notes：点为二项 GLMM 的 emmean（Prob 尺度），误差棒为 95% 渐近置信区间，字母为 Sidak correction的两两comparison。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/01_inoculation_and_field_disease/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
out_dir <- file.path(project_root, "result", "01_inoculation_and_field_disease")
path_main <- file.path(out_dir, "04_表_field_disease_main_GLMM.xlsx")
path_supp <- file.path(out_dir, "05_表_field_disease_supplement_Smutand36_74.xlsx")

stopifnot(file.exists(path_main), file.exists(path_supp))

bioRtools::read_data(path_main, sheet = "emmean_letters_main") -> df.main
bioRtools::read_data(path_supp, sheet = "emmean_letters_all") -> df.all

stopifnot(nrow(df.main) > 0, nrow(df.all) > 0)

# 统一因子水平，保证面板and小平面标签顺序稳定
prep_tab <- function(d) {
  d %>%
    dplyr::mutate(
      disease = factor(disease, levels = c("Blast", "BrownSpot", "Smut")),
      focal = factor(focal, levels = c("15", "36", "49", "69", "74", "87")),
      neighbor = factor(neighbor, levels = c("15", "36", "49", "69", "74", "87")),
      condition = factor(condition, levels = c("Pure", "Mix"))
    )
}

df.main <- prep_tab(df.main)
df.all <- prep_tab(df.all)

plot_field <- function(d, size) {
  d %>%
    ggplot2::ggplot(aes(x = neighbor, y = prob, color = condition)) +
    ggplot2::geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::geom_text(
      aes(y = asymp.UCL, label = signif),
      vjust = -0.6, size = 3.4, show.legend = FALSE
    ) +
    ggplot2::facet_grid(disease ~ focal, scales = "free_x", space = "free_x") +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      expand = ggplot2::expansion(mult = c(0.02, 0.10))
    ) +
    bioRtools::scale_color_research() +
    ggplot2::labs(
      x = "Neighbor",
      y = "Disease incidence (probability)",
      color = "Planting mode"
    ) +
    bioRtools::theme_prism() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 11, face = "bold"),
      legend.position = "bottom"
    ) -> p

  p
}

p.s2g <- plot_field(df.main)

bioRtools::write_data(p.s2g, file = file.path(out_dir, "13_FigureS2G_field_disease调查.pdf"),
                      width = 12, height = 6, units = "in")
bioRtools::write_data(p.s2g, file = file.path(out_dir, "13_FigureS2G_field_disease调查.png"),
                      width = 12, height = 6, units = "in", dpi = 500)
cat("[05] 输出：13_FigureS2G_field_disease调查.pdf / .png\n")

p.s2h <- plot_field(df.all)

bioRtools::write_data(p.s2h, file = file.path(out_dir, "14_FigureS2H_field_disease调查_全部纯系andSmut.pdf"),
                      width = 16, height = 9, units = "in")
bioRtools::write_data(p.s2h, file = file.path(out_dir, "14_FigureS2H_field_disease调查_全部纯系andSmut.png"),
                      width = 16, height = 9, units = "in", dpi = 500)
cat("[05] 输出：14_FigureS2H_field_disease调查_全部纯系andSmut.pdf / .png\n")

cat("\n[05] sessionInfo()\n")
sessionInfo()
