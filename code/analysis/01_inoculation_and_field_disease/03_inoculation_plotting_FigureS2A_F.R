# ============================================================
# 01_inoculation_and_field_disease / 03_inoculation_plotting_FigureS2A_F.R
# 作用：绘制 Figure S2A–S2F
#       S2A inoculationdata整体 Pure vs Mix（normalized前）
#       S2B inoculationdata整体 Pure vs Mix（按 病害 × 批次 Z-score normalized后）
#       S2C–S2F focal 15 / 49 / 69 / 87 各邻体组合comparison（Z-score，字母标记）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/01_inoculation_and_field_disease/03_inoculation_plotting_FigureS2A_F.R"
#
# 输入：<root>/result/01_inoculation_and_field_disease/01_表_inoculationLMM总体test.xlsx（sheet fixed_effect）
#       <root>/result/01_inoculation_and_field_disease/02_表_inoculationZscore中间表.xlsx
#       <root>/result/01_inoculation_and_field_disease/03_表_inoculation各focal组合emmeans字母.xlsx（sheet emmean_letters）
# 输出：<root>/result/01_inoculation_and_field_disease/07–12_FigureS2A–S2F_*.pdf / .png
# 其中 <root> = <repository_root>
#
# plotting约定：bioRtools::theme_prism() 主题；bioRtools::scale_fill_research() / scale_color_research() 配色；
#          bioRtools::write_data() 输出（PDF 矢量 + PNG 500 dpi），尺寸以英寸为单位（write_data 默认 cm，必须显式 units = "in"）。
# statistics标注notes：图内 P 值来自 01 脚本的固定effect contrast（Satterthwaite）；字母来自 02 脚本的
#          两两comparison Sidak correction（and旧版 emmeans 自动把 tukey 改为 sidak 的实际行为一致）。
#          P 值统一用 format_p() 输出为定点写法（P < 0.001 / P = 0.012），
#          避免 R 默认对小于 0.001 的数值输出科学计数法（如 2.1e-04）造成误读。
#          各组合的 emmean 因批次方差主导、df ≈ 1–2，其 CI 不可用于comparison，故本图不作误差棒。
# 坐标notes：本表不含gene组坐标，不涉及 0-based / 1-based 换算。
# versionrecord：脚本末尾打印 sessionInfo()；运行record见 code/01_inoculation_and_field_disease/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- normalizePath(getwd())
out_dir <- file.path(project_root, "result", "01_inoculation_and_field_disease")
path_lmm <- file.path(out_dir, "01_表_inoculationLMM总体test.xlsx")
path_z <- file.path(out_dir, "02_表_inoculationZscore中间表.xlsx")
path_letters <- file.path(out_dir, "03_表_inoculation各focal组合emmeans字母.xlsx")

stopifnot(file.exists(path_lmm), file.exists(path_z), file.exists(path_letters))

bioRtools::read_data(path_z) -> df.z
bioRtools::read_data(path_lmm, sheet = "fixed_effect") -> df.fixed
bioRtools::read_data(path_letters, sheet = "emmean_letters") -> df.letters

stopifnot(nrow(df.z) == 1659, nrow(df.letters) == 32)

# 统一的 P 值写法：避免科学计数法（e 记号）
format_p <- function(p, digits = 3) {
  dplyr::if_else(
    p < 0.001,
    "P < 0.001",
    paste0("P = ", formatC(p, format = "f", digits = digits))
  )
}

# 病害顺序and分组因子（图例顺序）
df.z %>%
  dplyr::mutate(
    disease = factor(disease, levels = c("blast", "brown spot")),
    group = factor(group, levels = c("Pure", "Mix")),
    batch = as.factor(batch),
    focal = factor(focal, levels = c("15", "49", "69", "87")),
    neighbor = factor(neighbor, levels = c("15", "49", "69", "87"))
  ) -> df.z

# 总体 Pure vs Mix 的 P 值（固定effect contrast）
df.fixed %>%
  dplyr::filter(term == "groupMix") %>%
  dplyr::select(disease, p_value) %>%
  dplyr::mutate(
    disease = factor(disease, levels = c("blast", "brown spot")),
    label = format_p(p_value)
  ) -> df.p

# 字母标记在 Z-score 面板下方的位置（按 facet 计算，避免被裁切）
df.z %>%
  dplyr::group_by(disease) %>%
  dplyr::summarise(y_letter = min(lesion_norm) - 0.08 * diff(range(lesion_norm)), .groups = "drop") -> df.letter.pos

df.letters %>%
  dplyr::mutate(
    disease = factor(disease, levels = c("blast", "brown spot")),
    neighbor = factor(neighbor, levels = c("15", "49", "69", "87"))
  ) %>%
  dplyr::left_join(df.letter.pos, by = "disease") -> df.letters

# P 值标注在面板内的位置
df.z %>%
  dplyr::group_by(disease) %>%
  dplyr::summarise(y_p = max(lesion_norm) + 0.05 * diff(range(lesion_norm)), .groups = "drop") %>%
  dplyr::left_join(df.p, by = "disease") -> df.p.pos

# ---- Figure S2A：normalized前 ----
df.z %>%
  ggplot2::ggplot(aes(x = group, y = lesion_number_per_cm2)) +
  ggplot2::geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.6, width = 0.6) +
  ggplot2::geom_jitter(aes(color = batch), width = 0.2, alpha = 0.8, size = 1.2) +
  ggplot2::facet_wrap(~disease, scales = "free_y") +
  bioRtools::scale_fill_research() +
  bioRtools::scale_color_research() +
  ggplot2::labs(
    x = "Planting mode",
    y = expression("Lesion number per cm"^2),
    fill = "Planting mode",
    color = "Batch"
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) -> p.s2a

bioRtools::write_data(p.s2a, file = file.path(out_dir, "07_FigureS2A_inoculation总体_normalized前.pdf"),
                      width = 8, height = 6, units = "in")
bioRtools::write_data(p.s2a, file = file.path(out_dir, "07_FigureS2A_inoculation总体_normalized前.png"),
                      width = 8, height = 6, units = "in", dpi = 500)
cat("[03] 输出：07_FigureS2A_inoculation总体_normalized前.pdf / .png\n")

# ---- Figure S2B：Z-score normalized后 ----
df.z %>%
  ggplot2::ggplot(aes(x = group, y = lesion_norm)) +
  ggplot2::geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.6, width = 0.6) +
  ggplot2::geom_jitter(aes(color = batch), width = 0.2, alpha = 0.8, size = 1.2) +
  ggplot2::geom_text(
    data = df.p.pos,
    aes(x = 1.5, y = y_p, label = label),
    inherit.aes = FALSE, size = 4
  ) +
  ggplot2::facet_wrap(~disease, scales = "free_y") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.10))) +
  bioRtools::scale_fill_research() +
  bioRtools::scale_color_research() +
  ggplot2::labs(
    x = "Planting mode",
    y = "Relative lesion density (Z-score)",
    fill = "Planting mode",
    color = "Batch"
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) -> p.s2b

bioRtools::write_data(p.s2b, file = file.path(out_dir, "08_FigureS2B_inoculation总体_normalized后.pdf"),
                      width = 8, height = 6, units = "in")
bioRtools::write_data(p.s2b, file = file.path(out_dir, "08_FigureS2B_inoculation总体_normalized后.png"),
                      width = 8, height = 6, units = "in", dpi = 500)
cat("[03] 输出：08_FigureS2B_inoculation总体_normalized后.pdf / .png\n")

# ---- Figure S2C–S2F：各 focal 的邻体组合comparison ----
plot_by_focal <- function(f, panel_id, panel_label) {
  df.z %>%
    dplyr::filter(as.character(focal) == f) %>%
    droplevels() -> d

  df.letters %>%
    dplyr::filter(as.character(focal) == f) %>%
    droplevels() -> l

  d %>%
    ggplot2::ggplot(aes(x = neighbor, y = lesion_norm)) +
    ggplot2::geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.6, width = 0.6) +
    ggplot2::geom_jitter(aes(color = batch), width = 0.2, alpha = 0.8, size = 1.2) +
    ggplot2::geom_text(
      data = l,
      aes(x = neighbor, y = y_letter, label = signif),
      inherit.aes = FALSE, size = 4, color = "black"
    ) +
    ggplot2::facet_wrap(~disease, scales = "free_y") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.05))) +
    bioRtools::scale_fill_research() +
    bioRtools::scale_color_research() +
    ggplot2::labs(
      x = "Neighbor",
      y = "Relative lesion density (Z-score)",
      fill = "Planting mode",
      color = "Batch"
    ) +
    bioRtools::theme_prism() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      legend.position = "right"
    ) -> p

  bioRtools::write_data(p, file = file.path(out_dir, paste0(panel_id, "_", panel_label, ".pdf")),
                        width = 8, height = 6, units = "in")
  bioRtools::write_data(p, file = file.path(out_dir, paste0(panel_id, "_", panel_label, ".png")),
                        width = 8, height = 6, units = "in", dpi = 500)
  cat("[03] 输出：", panel_id, "_", panel_label, ".pdf / .png\n", sep = "")
}

plot_by_focal("15", "09", "FigureS2C_inoculation各focalcomparison_15")
plot_by_focal("49", "10", "FigureS2D_inoculation各focalcomparison_49")
plot_by_focal("69", "11", "FigureS2E_inoculation各focalcomparison_69")
plot_by_focal("87", "12", "FigureS2F_inoculation各focalcomparison_87")

cat("\n[03] sessionInfo()\n")
sessionInfo()
