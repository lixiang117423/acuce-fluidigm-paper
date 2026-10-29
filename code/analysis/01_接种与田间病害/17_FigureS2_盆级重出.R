# ============================================================
# 01_接种与田间病害 / 17_FigureS2_盆级重出.R
# 作用：按审稿意见，以盆均值重做 Figure S2A–S2F
# 输入：data/01_分析输入表/01_接种结果_合并表.xlsx
# 输出：
#   result/01_接种与田间病害/17_表_接种FigureS2_盆级统计.xlsx
#   result/08_投稿规格图/FigureS2_v43/07–12_FigureS2A–S2F_盆级*.pdf / .png
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_in <- file.path(project_root, "data", "01_分析输入表", "01_接种结果_合并表.xlsx")
out_dir <- file.path(project_root, "result", "01_接种与田间病害")
out_pub <- file.path(project_root, "result", "08_投稿规格图", "FigureS2_v43r2")
path_out <- file.path(out_dir, "17b_表_接种FigureS2_盆级统计.xlsx")

stopifnot(file.exists(path_in))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_pub, showWarnings = FALSE, recursive = TRUE)

bioRtools::read_data(path_in) -> raw

raw %>%
  dplyr::mutate(
    batch = factor(batch),
    group = factor(group, levels = c("Pure", "Mix")),
    focal = factor(focal, levels = c("15", "49", "69", "87")),
    neighbor = factor(neighbor, levels = c("15", "49", "69", "87")),
    disease = factor(disease, levels = c("blast", "brown spot")),
    image = factor(image)
  ) -> d

# image 是盆；同一盆内叶片先取均值。
d %>%
  dplyr::group_by(disease, batch, image, group, focal, neighbor) %>%
  dplyr::summarise(
    lesion_pot_mean = mean(lesion_number_per_cm2, na.rm = TRUE),
    n_leaf = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::group_by(disease, batch) %>%
  dplyr::mutate(
    lesion_pot_norm = as.numeric(scale(lesion_pot_mean))
  ) %>%
  dplyr::ungroup() -> pot

stopifnot(nrow(pot) == 160)

overall <- purrr::map_dfr(levels(pot$disease), function(ds) {
  x <- pot %>% dplyr::filter(disease == ds) %>% droplevels()
  m_raw <- stats::lm(lesion_pot_mean ~ group + batch, data = x)
  m_norm <- stats::lm(lesion_pot_norm ~ group + batch, data = x)
  tibble::tibble(
    disease = ds,
    n_pot = nrow(x),
    raw_estimate = unname(stats::coef(m_raw)["groupMix"]),
    raw_p = summary(m_raw)$coefficients["groupMix", "Pr(>|t|)"],
    norm_estimate = unname(stats::coef(m_norm)["groupMix"]),
    norm_p = summary(m_norm)$coefficients["groupMix", "Pr(>|t|)"]
  )
})

grid <- tidyr::expand_grid(
  disease = levels(pot$disease),
  focal = levels(pot$focal)
) %>%
  dplyr::mutate(
    disease = factor(disease, levels = levels(pot$disease)),
    focal = factor(focal, levels = levels(pot$focal))
  ) %>%
  dplyr::arrange(disease, focal)

letters_tab <- purrr::pmap_dfr(grid, function(disease, focal) {
  x <- pot %>%
    dplyr::filter(.data$disease == .env$disease, .data$focal == .env$focal) %>%
    droplevels()
  m <- stats::lm(lesion_pot_mean ~ neighbor + batch, data = x)
  emm <- emmeans::emmeans(m, "neighbor")
  emmeans::contrast(emm, "pairwise", adjust = "sidak") %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(disease = .env$disease, focal = .env$focal)
})

letters_cld <- purrr::pmap_dfr(grid, function(disease, focal) {
  x <- pot %>%
    dplyr::filter(.data$disease == .env$disease, .data$focal == .env$focal) %>%
    droplevels()
  m <- stats::lm(lesion_pot_mean ~ neighbor + batch, data = x)
  emmeans::emmeans(m, "neighbor") %>%
    multcomp::cld(
      alpha = 0.05,
      Letters = letters,
      adjust = "sidak",
      decreasing = TRUE
    ) %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(disease = .env$disease, focal = .env$focal)
})

format_p <- function(p, digits = 3) {
  dplyr::if_else(
    p < 0.001,
    "P < 0.001",
    paste0("P = ", formatC(p, format = "f", digits = digits))
  )
}

overall <- overall %>%
  dplyr::mutate(
    raw_label = format_p(raw_p),
    norm_label = format_p(norm_p)
  )

# A：原始盆均值
pot %>%
  ggplot2::ggplot(aes(x = group, y = lesion_pot_mean)) +
  ggplot2::geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.6, width = 0.6) +
  ggplot2::geom_jitter(aes(color = batch), width = 0.15, alpha = 0.8, size = 1.4) +
  ggplot2::geom_text(
    data = overall,
    aes(x = -Inf, y = Inf, label = raw_label),
    inherit.aes = FALSE, hjust = -0.05, vjust = 1.7, size = 2.6
  ) +
  ggplot2::facet_wrap(~disease, scales = "free_y") +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.5)) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.12))) +
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
    plot.margin = ggplot2::margin(3, 2, 2, 2),
    legend.position = "right"
  ) -> p_a

bioRtools::write_data(
  p_a,
  file = file.path(out_pub, "07_FigureS2A_接种总体_盆均值.pdf"),
  width = 8.3, height = 6.0
)
bioRtools::write_data(
  p_a,
  file = file.path(out_pub, "07_FigureS2A_接种总体_盆均值.png"),
  width = 8.3, height = 6.0, dpi = 600
)

# B：盆均值按 病害 × batch Z-score
pot %>%
  ggplot2::ggplot(aes(x = group, y = lesion_pot_norm)) +
  ggplot2::geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.6, width = 0.6) +
  ggplot2::geom_jitter(aes(color = batch), width = 0.15, alpha = 0.8, size = 1.4) +
  ggplot2::facet_wrap(~disease, scales = "free_y") +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.5)) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.12))) +
  bioRtools::scale_fill_research() +
  bioRtools::scale_color_research() +
  ggplot2::labs(
    x = "Planting mode",
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
  ) -> p_b

bioRtools::write_data(
  p_b,
  file = file.path(out_pub, "08_FigureS2B_接种总体_盆均值标准化.pdf"),
  width = 8.3, height = 6.0
)
bioRtools::write_data(
  p_b,
  file = file.path(out_pub, "08_FigureS2B_接种总体_盆均值标准化.png"),
  width = 8.3, height = 6.0, dpi = 600
)

plot_focal <- function(f, panel_id, label) {
  dd <- pot %>% dplyr::filter(focal == f) %>% droplevels()
  ll <- letters_cld %>% dplyr::filter(focal == f) %>% droplevels()
  dd %>%
    ggplot2::ggplot(aes(x = neighbor, y = lesion_pot_mean)) +
    ggplot2::geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.6, width = 0.6) +
    ggplot2::geom_jitter(aes(color = batch), width = 0.12, alpha = 0.8, size = 1.5) +
    ggplot2::geom_text(
      data = ll,
      aes(x = neighbor, y = emmean, label = .group),
      inherit.aes = FALSE, vjust = -1.2, size = 2.8
    ) +
    ggplot2::facet_wrap(~disease, scales = "free_y") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.14))) +
    bioRtools::scale_fill_research() +
    bioRtools::scale_color_research() +
    ggplot2::labs(
      x = "Neighbour accession",
      y = expression("Lesion number per cm"^2),
      fill = "Planting mode",
      color = "Batch",
      title = paste("Focal accession", f)
    ) +
    bioRtools::theme_prism() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 8, face = "bold", hjust = 0),
      strip.text = ggplot2::element_text(size = 7, face = "bold"),
      strip.clip = "off",
      legend.position = "bottom",
      legend.box = "vertical",
      legend.key.size = ggplot2::unit(0.24, "cm"),
      legend.text = ggplot2::element_text(size = 5.4),
      legend.title = ggplot2::element_text(size = 5.8),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = ggplot2::unit(0.08, "cm"),
      axis.text.x = ggplot2::element_text(size = 6),
      axis.title.y = ggplot2::element_text(size = 6.5),
      plot.margin = ggplot2::margin(3, 2, 2, 2)
    ) -> p

  bioRtools::write_data(
    p,
    file = file.path(out_pub, paste0(panel_id, "_", label, ".pdf")),
    width = 8.3, height = 6.0
  )
  bioRtools::write_data(
    p,
    file = file.path(out_pub, paste0(panel_id, "_", label, ".png")),
    width = 8.3, height = 6.0, dpi = 600
  )
  cat("[17] 输出：", panel_id, "_", label, "\n", sep = "")
}

plot_focal("15", "09", "FigureS2C_盆级_focal15")
plot_focal("49", "10", "FigureS2D_盆级_focal49")
plot_focal("69", "11", "FigureS2E_盆级_focal69")
plot_focal("87", "12", "FigureS2F_盆级_focal87")

writexl::write_xlsx(
  list(
    pot_data = pot,
    overall = overall,
    pairwise_sidak = letters_tab,
    emmeans_letters = letters_cld,
    methods = tibble::tibble(
      unit = "pot (image)",
      model = "pot mean lesion density ~ group + batch; focal-wise neighbour models use lm",
      adjustment = "Sidak",
      note = "Figure S2C-F no longer use leaf-level observations or leaf-level significance letters"
    )
  ),
  path = path_out
)

cat("[17] 输出：", path_out, "\n")
cat("\n[17] sessionInfo()\n")
sessionInfo()
