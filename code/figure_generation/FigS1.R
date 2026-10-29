# Figure S1: experimental design.
# This is a schematic; no Source Data sheet is required.
# Input: values transcribed from the manuscript Methods and verified analysis metadata.
# Output: figures/FigS1.pdf and figures/FigS1.png.

library(ggplot2)

boxes <- data.frame(
  panel = c(rep("A", 3), rep("B", 6)),
  xmin = c(0.4, 5.4, 13.1, 0.7, 5.0, 9.7, 14.0, 0.7, 9.0),
  xmax = c(5.0, 12.7, 17.6, 4.3, 9.0, 13.3, 17.3, 8.3, 17.3),
  ymin = c(rep(6.7, 3), rep(3.0, 4), 0.5, 0.5),
  ymax = c(rep(8.4, 3), rep(4.7, 4), 2.2, 2.2),
  label = c(
    "Plant material<br>4 Acuce lines<br>15 / 49 / 69 / 87",
    "16 focal-neighbour combinations<br>Pure or Mix<br>2-3 batches; 4-6 pots per combination<br>64 blast / 96 brown spot",
    "Inoculation<br>*M. oryzae* CH1189<br>*B. oryzae* CH9015",
    "Seedlings<br>one plant per pot",
    "Outdoor PHIM<br>2 plants per pot<br>10 replicate pots",
    "Weekly height<br>15 surveys<br>Yield at 153 d",
    "Three sampling stages<br>Tillering / flowering<br>Grain filling",
    "RNA extraction and Fluidigm<br>92 target genes, 20 pathways<br>80 / 86 / 88 samples",
    "Field survey in the Yuanyang terraces<br>3 years x 3 fields x 6 accessions<br>Blast, brown spot, smut, height, grain weight"
  )
)

arrows <- data.frame(
  x = c(5.0, 12.7, 4.3, 9.0, 13.3, 15.65, 8.3),
  y = c(7.55, 7.55, 3.85, 3.85, 3.85, 3.0, 1.35),
  xend = c(5.4, 13.1, 5.0, 9.7, 14.0, 15.65, 9.0),
  yend = c(7.55, 7.55, 3.85, 3.85, 3.85, 2.2, 1.35)
)

label_lines <- boxes %>%
  tidyr::separate_longer_delim(label, delim = "<br>") %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(
    line_n = dplyr::n(),
    line_id = dplyr::row_number(),
    x = (xmin + xmax) / 2,
    y = (ymin + ymax) / 2 + ((line_n + 1) / 2 - line_id) * 0.19
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    label_plot = dplyr::case_when(
      label == "*M. oryzae* CH1189" ~ "italic('M. oryzae')~CH1189",
      label == "*B. oryzae* CH9015" ~ "italic('B. oryzae')~CH9015",
      TRUE ~ label
    ),
    is_math = stringr::str_detect(label, "^\\*")
  )

p <- ggplot() +
  geom_rect(data = boxes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = panel),
            color = "grey30", linewidth = 0.35) +
  geom_text(data = dplyr::filter(label_lines, !is_math),
            aes(x = x, y = y, label = label),
            family = "Arial", size = 1.85) +
  geom_text(data = dplyr::filter(label_lines, is_math),
            aes(x = x, y = y, label = label_plot),
            family = "Arial", size = 1.85, parse = TRUE) +
  geom_segment(data = arrows,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = grid::arrow(length = grid::unit(0.14, "cm"), type = "closed"),
               color = "grey30", linewidth = 0.38) +
  annotate("text", x = 0.7, y = 8.8,
           label = "A  Greenhouse inoculation experiment",
           family = "Arial", size = 2.8, fontface = "bold", hjust = 0) +
  annotate("text", x = 0.7, y = 5.1,
           label = "B  Outdoor full-life-cycle experiment and field survey",
           family = "Arial", size = 2.8, fontface = "bold", hjust = 0) +
  scale_fill_manual(values = c(A = "#EEF2F7", B = "#F7F2EE"), guide = "none") +
  coord_cartesian(xlim = c(0, 18), ylim = c(0, 9.2), expand = FALSE) +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(1, 1, 1, 1, unit = "mm"))

out_pdf <- "figures/FigS1.pdf"
out_png <- "figures/FigS1.png"
if (file.exists(out_pdf) || file.exists(out_png)) {
  stop(
    "Refusing to overwrite an existing Figure S1. Archive the current PDF and PNG first.",
    call. = FALSE
  )
}

grDevices::cairo_pdf(out_pdf, width = 180 / 25.4, height = 120 / 25.4,
                     family = "Arial")
print(p)
dev.off()

ragg::agg_png(out_png, width = 180, height = 120,
              units = "mm", res = 600, background = "white")
print(p)
dev.off()
