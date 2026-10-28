# ============================================================
# check_source_data.R —— Source Data 自检（子图清单 ↔ sheet 名单 ↔ 代码）
# 运行：cd manuscript && Rscript code/check_source_data.R
# 依据：AGENTS「Source Data 交付前自检（三者对账，缺一不算完成）」
# 本机绝对路径：/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/check_source_data.R
# ============================================================

rm(list = ls())

library(tidyverse)

sd_file <- "source_data/Source_Data.xlsx"
stopifnot(file.exists(sd_file))

# 子图清单（正文 Fig. 1–5 + 补充 Fig. S1–S19；Fig. S1 为示意图，按 README 登记为例外）
panels <- c(
  paste0("Fig1", letters[1:5]),
  paste0("Fig2", letters[1:5]),
  paste0("Fig3", letters[1:2]),
  paste0("Fig4", letters[1:2]),
  paste0("Fig5", letters[1:3]),
  paste0("FigS2", letters[1:8]),
  "FigS3", "FigS4", "FigS5",
  paste0("FigS6", letters[1:2]),
  "FigS7",
  paste0("FigS8", letters[1:2]),
  paste0("FigS9", letters[1:2]),
  "FigS10",
  paste0("FigS11", letters[1:3]),
  paste0("FigS12", letters[1:3]),
  "FigS13", "FigS14",
  paste0("FigS15", letters[1:2]),
  "FigS16", "FigS17", "FigS18",
  paste0("FigS19", letters[1:2])
)
exempt <- c("FigS1")  # 示意图：无 Source Data（README 已登记）
code_files <- paste0("code/", unique(stringr::str_remove(panels, "[a-z]$")), ".R")

sheets <- readxl::excel_sheets(sd_file)

ok <- TRUE
report <- function(cond, msg) {
  cat(if (cond) "  [OK]   " else "  [FAIL] ", msg, "\n", sep = "")
  if (!cond) ok <<- FALSE
}

cat("=== Source Data 自检 ===\n")
report(identical(sheets, panels), paste0("sheet 名单与子图清单一一对应（sheet = ", length(sheets),
                                         "，子图 = ", length(panels), "）"))
report(all(!stringr::str_detect(sheets, paste0("^(", paste0(stringr::str_remove(panels, "[a-z]$"),
                                                            collapse = "|"), ")$")) |
             sheets %in% c("FigS3", "FigS4", "FigS5", "FigS7", "FigS10", "FigS13",
                           "FigS14", "FigS16", "FigS17", "FigS18")),
       "多面板图不得出现裸 FigN、单面板图不得带字母")
report(all(file.exists(code_files)), paste0("每个图都有对应脚本（", length(code_files), " 个）"))

bad_panel <- c()
empty <- c()
for (s in sheets) {
  d <- readxl::read_excel(sd_file, sheet = s, skip = 2)
  if (nrow(d) < 1) empty <- c(empty, s)
  suffix <- stringr::str_extract(s, "[a-z]$")
  if (!is.na(suffix)) {
    if (!"Panel" %in% colnames(d) || !all(as.character(d$Panel) == suffix)) bad_panel <- c(bad_panel, s)
  }
}
report(length(empty) == 0, paste0("每个 sheet 至少一行数据（空 sheet：", length(empty), " 个）"))
report(length(bad_panel) == 0, paste0("多面板 sheet 的 Panel 列与后缀一致（异常：", length(bad_panel), " 个）"))

cat("\n=== 例外登记（README）===\n")
cat("  ", paste(exempt, collapse = ", "), "：示意图，无 Source Data\n", sep = "")
cat("\n结论：", if (ok) "全部通过 ✓" else "存在未通过项 ✗", "\n", sep = "")

cat("\n[check_source_data] sessionInfo()\n")
sessionInfo()
