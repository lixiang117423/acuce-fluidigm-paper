# ============================================================
# 01_接种与田间病害 / 03_接种绘图_FigureS2A_F.R
# 作用：绘制 Figure S2A–S2F
#       S2A 接种数据整体 Pure vs Mix（标准化前）
#       S2B 接种数据整体 Pure vs Mix（按 病害 × 批次 Z-score 标准化后）
#       S2C–S2F focal 15 / 49 / 69 / 87 各邻体组合比较（Z-score，字母标记）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/01_接种与田间病害/03_接种绘图_FigureS2A_F.R"
#
# 输入：<root>/result/01_接种与田间病害/01_表_接种LMM总体检验.xlsx（sheet fixed_effect）
#       <root>/result/01_接种与田间病害/02_表_接种Zscore中间表.xlsx
#       <root>/result/01_接种与田间病害/03_表_接种各focal组合emmeans字母.xlsx（sheet emmean_letters）
# 输出：<root>/result/01_接种与田间病害/07–12_FigureS2A–S2F_*.pdf / .png
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 绘图约定：bioRtools::theme_prism() 主题；bioRtools::scale_fill_research() / scale_color_research() 配色；
#          bioRtools::write_data() 输出（PDF 矢量 + PNG 500 dpi），尺寸以英寸为单位（write_data 默认 cm，必须显式 units = "in"）。
# 统计标注说明：图内 P 值来自 01 脚本的固定效应 contrast（Satterthwaite）；字母来自 02 脚本的
#          两两比较 Sidak 校正（与旧版 emmeans 自动把 tukey 改为 sidak 的实际行为一致）。
#          P 值统一用 format_p() 输出为定点写法（P < 0.001 / P = 0.012），
#          避免 R 默认对小于 0.001 的数值输出科学计数法（如 2.1e-04）造成误读。
#          各组合的 emmean 因批次方差主导、df ≈ 1–2，其 CI 不可用于比较，故本图不作误差棒。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录见 code/01_接种与田间病害/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
out_dir <- file.path(project_root, "result", "01_接种与田间病害")
out_pub_dir <- file.path(project_root, "result", "08_投稿规格图")  # 投稿规格面板输出目录
dir.create(out_pub_dir, showWarnings = FALSE, recursive = TRUE)
path_lmm <- file.path(out_dir, "01_表_接种LMM总体检验.xlsx")
path_z <- file.path(out_dir, "02_表_接种Zscore中间表.xlsx")
path_letters <- file.path(out_dir, "03_表_接种各focal组合emmeans字母.xlsx")

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

# 病害顺序与分组因子（图例顺序）
df.z %>%
  dplyr::mutate(
    disease = factor(disease, levels = c("blast", "brown spot")),
    group = factor(group, levels = c("Pure", "Mix")),
    batch = as.factor(batch),
    focal = factor(focal, levels = c("15", "49", "69", "87")),
    neighbor = factor(neighbor, levels = c("15", "49", "69", "87"))
  ) -> df.z

# 总体 Pure vs Mix 的 P 值（固定效应 contrast）
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

# ---- Figure S2A：标准化前 ----
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
    strip.text = ggplot2::element_text(size = 8, face = "bold"),
    strip.clip = "off",
    legend.position = "right"
  ) -> p.s2a

bioRtools::write_data(p.s2a, file = file.path(out_pub_dir, "07_FigureS2A_接种总体_标准化前.pdf"),
                      width = 8.3, height = 6.0, units = "cm")
bioRtools::write_data(p.s2a, file = file.path(out_pub_dir, "07_FigureS2A_接种总体_标准化前.png"),
                      width = 8.3, height = 6.0, units = "cm", dpi = 600)
cat("[03] 输出：07_FigureS2A_接种总体_标准化前.pdf / .png\n")

# ---- Figure S2B：Z-score 标准化后 ----
df.z %>%
  ggplot2::ggplot(aes(x = group, y = lesion_norm)) +
  ggplot2::geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.6, width = 0.6) +
  ggplot2::geom_jitter(aes(color = batch), width = 0.2, alpha = 0.8, size = 1.2) +
  ggplot2::geom_text(
    data = df.p.pos,
    # 2026-10-19 修：原先标注在图内居中（x = 1.5）、字号 4（约 11 pt），会压住 y 轴刻度；
    #   改为贴面板左缘（x = -Inf）左对齐、字号 2.6（约 7.4 pt）
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, size = 2.6, hjust = -0.05, vjust = 1.6
  ) +
  ggplot2::facet_wrap(~disease, scales = "free_y") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.10))) +
  bioRtools::scale_fill_research() +
  bioRtools::scale_color_research() +
  ggplot2::labs(
    x = "Planting mode",
    # 2026-10-19：旋转的 y 轴标题在 60 mm 高面板里两头被切（缺 2.4 pt）→ 缩短 + 缩小字号
    y = "Relative lesion density (Z)",
    fill = "Planting mode",
    color = "Batch"
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 8, face = "bold"),
    strip.clip = "off",
    axis.title.y = ggplot2::element_text(size = 6.5),
    plot.margin = ggplot2::margin(3, 2, 2, 2),
    legend.position = "right"
  ) -> p.s2b

bioRtools::write_data(p.s2b, file = file.path(out_pub_dir, "08_FigureS2B_接种总体_标准化后.pdf"),
                      width = 8.3, height = 6.0, units = "cm")
bioRtools::write_data(p.s2b, file = file.path(out_pub_dir, "08_FigureS2B_接种总体_标准化后.png"),
                      width = 8.3, height = 6.0, units = "cm", dpi = 600)
cat("[03] 输出：08_FigureS2B_接种总体_标准化后.pdf / .png\n")

# ---- Figure S2C–S2F：各 focal 的邻体组合比较 ----
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
      inherit.aes = FALSE, size = 3, color = "black"
    ) +
    ggplot2::facet_wrap(~disease, scales = "free_y") +
    # 2026-10-19 修：图例置底（原先右侧的 8 条图例把两分面各挤到约 50 pt 宽），
    #   分面标题 12 pt → 6.5 pt，并关掉 strip.clip（默认 "on" 会把比 strip 宽的标题两侧切掉，
    #   实测 S2D 的 "brown spot" 被渲染成 "own s"）
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
                    color = ggplot2::guide_legend(nrow = 1, byrow = TRUE)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.05))) +
    bioRtools::scale_fill_research() +
    bioRtools::scale_color_research() +
    ggplot2::labs(
      x = "Neighbor",
      # 2026-10-19：旋转的 y 轴标题在 60 mm 高面板里两头被切 → 缩短 + 缩小字号
      y = "Relative lesion density (Z)",
      fill = "Planting mode",
      color = "Batch"
    ) +
    bioRtools::theme_prism() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 6.5, face = "bold"),
      strip.clip = "off",
      # 二次修：置底图例两列并排时比面板宽（左切 "Pl"、右切 "3"）→ 改为上下两行堆叠
      legend.position = "bottom",
      legend.box = "vertical",
      legend.key.size = ggplot2::unit(0.26, "cm"),
      legend.text = ggplot2::element_text(size = 5.5),
      legend.title = ggplot2::element_text(size = 6),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = ggplot2::unit(0.08, "cm"),
      axis.text.x = ggplot2::element_text(size = 6),
      axis.title.y = ggplot2::element_text(size = 6.5),
      plot.margin = ggplot2::margin(3, 2, 2, 2)
    ) -> p

  # 投稿规格：面板按最终 rect 尺寸（8.3 × 6.0 cm）导出，落版即设计字号
  # （2026-10-19：图例移到面板底部后绘图区略矮，但**不改面板高度**——S2 整页已 365 mm，
  #   再长就要拆页；组图脚本里 C–F 的 h 仍为 60 mm）
  bioRtools::write_data(p, file = file.path(out_pub_dir, paste0(panel_id, "_", panel_label, ".pdf")),
                        width = 8.3, height = 6.0, units = "cm")
  bioRtools::write_data(p, file = file.path(out_pub_dir, paste0(panel_id, "_", panel_label, ".png")),
                        width = 8.3, height = 6.0, units = "cm", dpi = 600)
  cat("[03] 输出：", panel_id, "_", panel_label, ".pdf / .png\n", sep = "")
}

plot_by_focal("15", "09", "FigureS2C_接种各focal比较_15")
plot_by_focal("49", "10", "FigureS2D_接种各focal比较_49")
plot_by_focal("69", "11", "FigureS2E_接种各focal比较_69")
plot_by_focal("87", "12", "FigureS2F_接种各focal比较_87")

cat("\n[03] sessionInfo()\n")
sessionInfo()
