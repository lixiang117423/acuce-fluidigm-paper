# ============================================================
# check_source_data_v45.R
# Source Data v45 自检：sheet and重编号后的supplementary_figures一一对应
# 运行：Rscript code/figure_generation/check_source_data_v45.R
# ============================================================

rm(list = ls())

library(tidyverse)

sd_file <- "data/source_data/Source_Data_v45.xlsx"
stopifnot(file.exists(sd_file))

panels <- c(
  paste0("Fig1", letters[1:5]),
  paste0("Fig2", letters[1:5]),
  paste0("Fig3", letters[1:2]),
  paste0("Fig4", letters[1:2]),
  paste0("Fig5", letters[1:3]),
  paste0("FigS2", letters[1:8]),
  "FigS3", "FigS4", "FigS5",
  paste0("FigS6", letters[1:2]),
  paste0("FigS7", letters[1:2]),
  "FigS8",
  paste0("FigS9", letters[1:3]),
  paste0("FigS10", letters[1:3]),
  "FigS11",
  paste0("FigS12", letters[1:2]),
  "FigS13", "FigS14", "FigS15",
  paste0("FigS16", letters[1:2])
)

code_map <- c(
  FigS2 = "code/figure_generation/FigS2.R",
  FigS3 = "code/figure_generation/FigS3.R",
  FigS4 = "code/figure_generation/FigS4.R",
  FigS5 = "code/figure_generation/FigS5.R",
  FigS6 = "code/figure_generation/FigS6.R",
  FigS7 = "code/figure_generation/FigS7.R",
  FigS8 = "code/figure_generation/FigS8.R",
  FigS9 = "code/figure_generation/FigS9.R",
  FigS10 = "code/figure_generation/FigS10.R",
  FigS11 = "code/figure_generation/FigS11.R",
  FigS12 = "code/figure_generation/FigS12.R",
  FigS13 = "code/figure_generation/FigS13.R",
  FigS14 = "code/figure_generation/FigS14.R",
  FigS15 = "code/figure_generation/FigS15.R",
  FigS16 = "code/figure_generation/FigS16.R"
)

sheets <- readxl::excel_sheets(sd_file)

ok <- TRUE
report <- function(cond, msg) {
  cat(if (cond) "  [OK]   " else "  [FAIL] ", msg, "\n", sep = "")
  if (!cond) ok <<- FALSE
}

cat("=== Source Data v45 自检 ===\n")
report(identical(sheets, panels), paste0("sheet 名单一一对应（sheet = ", length(sheets),
                                         "，panel = ", length(panels), "）"))

code_files <- unique(c("code/figure_generation/Fig1.R", "code/figure_generation/Fig2.R", "code/figure_generation/Fig3.R", "code/figure_generation/Fig4.R",
                       "code/figure_generation/Fig5.R", unname(code_map)))
report(all(file.exists(code_files)), paste0("plotting/analysis脚本均存在（", length(code_files), " 个）"))

bad_panel <- c()
empty <- c()
for (s in sheets) {
  d <- readxl::read_excel(sd_file, sheet = s, skip = 2)
  if (nrow(d) < 1) empty <- c(empty, s)
  suffix <- stringr::str_extract(s, "[a-z]$")
  if (!is.na(suffix)) {
    if (!"Panel" %in% colnames(d) || !all(as.character(d$Panel) == suffix)) {
      bad_panel <- c(bad_panel, s)
    }
  }
}
report(length(empty) == 0, paste0("无空 sheet（空：", length(empty), "）"))
report(length(bad_panel) == 0, paste0("多面板 Panel 列and sheet 后缀一致（异常：",
                                      length(bad_panel), "）"))

report(all(c("FigS6a", "FigS6b", "FigS7a", "FigS7b", "FigS8", "FigS9a",
             "FigS10a", "FigS11", "FigS12a", "FigS13", "FigS14", "FigS15",
             "FigS16a", "FigS16b") %in% sheets),
       "重编号后的supplementary_figures sheet 均存在")
report(!any(c("FigS6a_old", "FigS7", "FigS14_old") %in% sheets),
       "已移除的根系分泌物and重复图未残留")

cat("\n结论：", if (ok) "全部通过" else "存在未通过项", "\n", sep = "")
sessionInfo()
