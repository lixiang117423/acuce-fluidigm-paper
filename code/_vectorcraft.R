# ============================================================
# _vectorcraft.R —— 组图薄封装（供 Fig*.R 调用）
# 规则依据：AGENTS「组图（多面板拼 Figure）一律用 VectorCraft」——本文件只做
#   "读布局 → 发 file.new / file.place / text.create → --export"，版面计算交给 VectorCraft。
# 输入：panels_mm = list(list(path=, rect=c(x,y,w,h), label="A", label_pos=c(x,y)))，单位 mm（原点左上）
# 输出：out_pdf（矢量）＋ 同名 .png（600 dpi，由同一 PDF 转，保证同源）
# 本机绝对路径（登记用）：
#   /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理/manuscript/code/_vectorcraft.R
# ============================================================

MM_TO_PT <- 72 / 25.4

compose_vectorcraft <- function(panels_mm, width_mm, height_mm, out_pdf,
                                label_size = 10, font = "Arial",
                                cli = Sys.getenv(
                                  "VECTORCRAFT_CLI",
                                  "/Users/lixiang/tools/vectorcraft/cli-0.12.0/vectorcraft-cli-0.12.0-macos-universal/vectorcraft-cli")) {
  if (!file.exists(cli)) {
    stop("VectorCraft CLI 未找到：", cli,
         "；请设置环境变量 VECTORCRAFT_CLI 指向 vectorcraft-cli 可执行文件")
  }
  layout <- list(
    width_mm = width_mm, height_mm = height_mm, font = font, label_size = label_size,
    panels = lapply(panels_mm, function(p) {
      list(path = normalizePath(p$path, mustWork = TRUE),
           rect = round(p$rect * MM_TO_PT, 2),
           label = p$label,
           label_pos = round(p$label_pos * MM_TO_PT, 2))
    })
  )
  layout_file <- tempfile(fileext = "_layout.json")
  jsonlite::write_json(layout, layout_file, auto_unbox = TRUE, pretty = TRUE)

  args <- c("run", "--cmd", "file.new", "--params",
            jsonlite::toJSON(list(width = paste0(width_mm, " mm"),
                                  height = paste0(height_mm, " mm"),
                                  units = "Millimeters"), auto_unbox = TRUE))
  for (p in layout$panels) {
    args <- c(args, "--cmd", "file.place", "--params",
              jsonlite::toJSON(list(path = p$path, rect = p$rect), auto_unbox = TRUE))
  }
  for (p in layout$panels) {
    if (!is.null(p$label) && nzchar(p$label)) {
      args <- c(args, "--cmd", "text.create", "--params",
                jsonlite::toJSON(list(x = p$label_pos[1], y = p$label_pos[2], text = p$label,
                                      size = label_size, font = font, style = "bold"),
                                 auto_unbox = TRUE))
    }
  }
  args <- c(args, "--export", out_pdf)
  # 注意：system2() 不对参数做引号处理，而 --params 的 JSON 里含空格（如 "174 mm"）会被拆开，
  # 因此优先用 processx::run()（逐参数传递）；没有 processx 时退回 system() + shQuote()。
  if (requireNamespace("processx", quietly = TRUE)) {
    res <- processx::run(cli, args, error_on_status = FALSE)
    status <- res$status
    if (status != 0) {
      stop("VectorCraft 导出失败（退出码 ", status, "）：\n", res$stderr)
    }
  } else {
    cmd <- paste(shQuote(cli), paste(shQuote(args), collapse = " "))
    status <- system(cmd)
    if (status != 0) stop("VectorCraft 导出失败（退出码 ", status, "）")
  }
  if (!file.exists(out_pdf)) stop("VectorCraft 未生成输出文件：", out_pdf)

  png_out <- sub("\\.pdf$", "", out_pdf)
  system2("pdftoppm", c("-r", "600", "-png", "-singlefile", out_pdf, png_out))
  cat("[vectorcraft] 输出：", out_pdf, " / ", png_out, ".png\n", sep = "")
  invisible(out_pdf)
}
