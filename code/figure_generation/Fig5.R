# ============================================================
# Fig5_v32.R —— Figure 5（second-round revision）
# Source Data：source_data/Source_Data_v32.xlsx（Fig5a–Fig5c）
# A: pot-level decomposition; B: plot-level yearly field effects;
# C: plot-level condition x year interaction tests.
# 运行：cd manuscript && /opt/homebrew/bin/Rscript code/Fig5_v32.R
# 输出：figures/Fig5_v32.pdf / .png；面板 figures/panels/Fig5a_v32.pdf 等
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

source("code/_vectorcraft.R")
sd_file <- "source_data/Source_Data_v32.xlsx"
rd <- function(sheet) readxl::read_excel(sd_file, sheet = sheet, skip = 2)
dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)
save_panel <- function(p, name, w_cm, h_cm) {
  bioRtools::write_data(
    p,
    file = file.path("figures/panels", paste0(name, ".pdf")),
    width = w_cm, height = h_cm, units = "cm"
  )
}
fmt_p <- function(p, digits = 2) {
  dplyr::if_else(
    p < 0.001, "P < 0.001",
    paste0("P = ", formatC(p, format = "f", digits = digits))
  )
}

# A: pot-level decomposition
d5a <- rd("Fig5a") %>%
  dplyr::mutate(
    term = factor(
      term,
      levels = c("focal", "neighbor", "focal:neighbor"),
      labels = c("Focal identity", "Neighbour identity", "Focal x neighbour")
    ),
    disease = factor(
      disease,
      levels = c("blast", "brown spot"),
      labels = c("Rice blast", "Brown spot")
    ),
    lab = fmt_p(p_value)
  )

p5a <- d5a %>%
  ggplot2::ggplot(ggplot2::aes(x = term, y = f_value, fill = term)) +
  ggplot2::geom_col(width = 0.68, show.legend = FALSE) +
  ggplot2::geom_text(ggplot2::aes(label = lab), vjust = -0.3, size = 2.8) +
  ggplot2::facet_wrap(~disease) +
  bioRtools::scale_fill_research() +
  ggplot2::scale_y_continuous(
    limits = c(0, NA),
    expand = ggplot2::expansion(mult = c(0, 0.25))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "F value (pot means)",
    title = "Pot-level decomposition",
    subtitle = "Batch fixed; incomplete mixture design"
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 8),
    axis.text.x = ggplot2::element_text(angle = 20, hjust = 1, vjust = 1, size = 7.5),
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6.5)
  )
save_panel(p5a, "Fig5a_v32", 17.0, 7.0)

# B: yearly field effects
d5b <- rd("Fig5b") %>%
  dplyr::mutate(
    category = factor(category, levels = c("Blast", "BrownSpot", "Smut", "PH", "KGW")),
    year = factor(year)
  )

ref_lines <- tibble::tibble(
  category = factor(c("Blast", "BrownSpot", "Smut", "PH", "KGW"),
                    levels = c("Blast", "BrownSpot", "Smut", "PH", "KGW")),
  ref = c(1, 1, 1, 0, 0)
)

inter_lab <- d5b %>%
  dplyr::distinct(category, p_interaction) %>%
  dplyr::mutate(lab = paste0("condition x Year\n", fmt_p(p_interaction)))

p5b <- d5b %>%
  ggplot2::ggplot(ggplot2::aes(x = estimate, y = year)) +
  ggplot2::geom_vline(
    data = ref_lines,
    ggplot2::aes(xintercept = ref),
    linetype = "dashed", color = "grey40"
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
    height = 0.2, color = "grey35"
  ) +
  ggplot2::geom_point(ggplot2::aes(color = year), size = 2.6) +
  ggplot2::geom_text(
    data = inter_lab,
    ggplot2::aes(x = -Inf, y = Inf, label = lab),
    hjust = -0.05, vjust = 1.35, size = 2.4, lineheight = 0.95,
    inherit.aes = FALSE
  ) +
  ggplot2::facet_wrap(~category, scales = "free_x", nrow = 1) +
  ggplot2::scale_y_discrete(limits = rev(levels(d5b$year))) +
  bioRtools::scale_color_research() +
  ggplot2::labs(
    x = "Effect per year (disease: OR; trait: mean difference)",
    y = "Year",
    title = "Year-by-year field effects",
    subtitle = "Plot-level aggregates"
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  bioRtools::theme_prism() +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 8),
    strip.clip = "off",
    legend.position = "none",
    plot.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6.5),
    axis.text.x = ggplot2::element_text(size = 7)
  )
save_panel(p5b, "Fig5b_v32", 17.0, 8.5)

# C: plot-level condition x year interaction tests
d5c <- rd("Fig5c") %>%
  dplyr::mutate(
    category = factor(category, levels = c("Blast", "BrownSpot", "Smut", "PH", "KGW")),
    lab = fmt_p(interaction_p)
  )

p5c <- d5c %>%
  ggplot2::ggplot(ggplot2::aes(x = category, y = interaction_F, fill = category)) +
  ggplot2::geom_col(width = 0.65, show.legend = FALSE) +
  ggplot2::geom_text(ggplot2::aes(label = lab), vjust = -0.3, size = 2.5) +
  ggplot2::scale_y_continuous(
    limits = c(0, NA),
    expand = ggplot2::expansion(mult = c(0, 0.28))
  ) +
  bioRtools::scale_fill_research() +
  ggplot2::labs(
    x = NULL,
    y = "Condition x Year F",
    title = "No significant year x treatment interaction"
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1, size = 7),
    plot.title = ggplot2::element_text(size = 8, face = "bold")
  )
save_panel(p5c, "Fig5c_v32", 17.0, 5.0)

compose_vectorcraft(
  list(
    list(path = "figures/panels/Fig5a_v32.pdf", rect = c(2, 8, 170, 70),
         label = "A", label_pos = c(2, 3.5)),
    list(path = "figures/panels/Fig5b_v32.pdf", rect = c(2, 88, 170, 85),
         label = "B", label_pos = c(2, 83.5)),
    list(path = "figures/panels/Fig5c_v32.pdf", rect = c(2, 178, 170, 50),
         label = "C", label_pos = c(2, 173.5))
  ),
  width_mm = 174, height_mm = 232, out_pdf = "figures/Fig5_v32.pdf"
)

cat("\n[Fig5_v32] sessionInfo()\n")
sessionInfo()
