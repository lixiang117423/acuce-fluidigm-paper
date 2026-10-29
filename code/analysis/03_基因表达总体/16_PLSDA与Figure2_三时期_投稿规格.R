# ============================================================
# 03_基因表达总体 / 16_PLSDA与Figure2_三时期_投稿规格.R
# 作用：**投稿规格版** Figure 2 面板（2026-10-18 新增；统计口径与 15 号脚本完全一致）
#       与 15 号脚本的差别只在版面：
#         ① 输出尺寸 = 面板在组图中的**最终 rect 尺寸**（cm，scale = 1），保证图内字号落版不缩水；
#         ② 2A 热图 92 个基因拆成**左右两块**（各 46 个基因），每块 83 mm 宽，行高 2.4 mm，
#            行标签 5-6 pt 仍不叠字（15 号是 11×13 in 大画布，缩到 83 mm 后行标签只剩 1.5 pt）；
#         ③ 去掉面板内 caption（Q2、检验口径等写进正文 Figure legend）；
#         ④ 输出目录改为 result/08_投稿规格图，原 result/03_基因表达总体 只读不改。
# 运行：/opt/homebrew/bin/Rscript "<root>/code/03_基因表达总体/16_PLSDA与Figure2_三时期_投稿规格.R"
#
# 输入：<root>/result/03_基因表达总体/23_表_表达矩阵_三时期_不剔离群值.xlsx（combo_matrix_median、gene_annotation）
# 输出：<root>/result/08_投稿规格图/24_表_PLSDA_三时期_中位数矩阵.xlsx（同 15 号脚本，仅换目录）
#       <root>/result/08_投稿规格图/25a/25b_Figure2A_*、26–29_Figure2B–2E_*.pdf / .png
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 为什么要换：`11_表` 的组合均值是"剔除离群值"后算的（29.55% 的格子与不剔除版本不同、个别格子差 4 个数量级），
#   与样本级统计（limma，未剔除）口径不一致；改用"不剔除 + 中位数"后图与统计同源，且热图不被极端样本拉成二值图。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()。
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
out_dir <- file.path(project_root, "result", "03_基因表达总体")
out_pub_dir <- file.path(project_root, "result", "08_投稿规格图")  # 投稿规格面板输出目录
dir.create(out_pub_dir, showWarnings = FALSE, recursive = TRUE)
path_mat <- file.path(out_dir, "23_表_表达矩阵_三时期_不剔离群值.xlsx")
path_out <- file.path(out_pub_dir, "24_表_PLSDA_三时期_中位数矩阵.xlsx")

stopifnot(file.exists(path_mat))

stage_levels <- c("tillering", "flowering", "filling")
plotgroup_levels <- as.vector(outer(stage_levels, c("Pure", "Mix"), paste))
sig_levels <- c("FDR < 0.05", "P < 0.05", "n.s.")
sig_colors <- c("FDR < 0.05" = "#a63446", "P < 0.05" = "#e0a458", "n.s." = "grey62")

format_p <- function(p, digits = 3) {
  dplyr::if_else(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, format = "f", digits = digits)))
}

# 投稿规格：紧凑写法（无空格、2 位小数），避免小面板里检验 P 值叠字；精确值见正文 Figure legend 指向的附表
format_p_compact <- function(p, digits = 2) {
  dplyr::if_else(p < 0.001, "P<0.001", paste0("P=", formatC(p, format = "f", digits = digits)))
}

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_mat, "combo_matrix_median") -> mat.combo
bioRtools::read_data(path_mat, sheet = "gene_annotation") -> df.gene

meta_cols <- c("combo", "stage", "focal", "neighbor", "group")
gene_cols <- setdiff(colnames(mat.combo), meta_cols)
stopifnot(length(gene_cols) == 92, nrow(mat.combo) == 48)

mat.combo %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels)) -> mat.combo

# ---------------- PLS-DA ----------------
x <- as.matrix(mat.combo[, gene_cols])
rownames(x) <- mat.combo$combo
y <- as.factor(mat.combo$group)

ropls::opls(x, y, predI = 3) -> model
stopifnot(length(ropls::getScoreMN(model)) > 0)

as.data.frame(model@summaryDF) %>% tibble::rownames_to_column(var = "model_type") -> cv_tab
cat("[15] PLS-DA（中位数组合矩阵）模型摘要：\n")
print(as.data.frame(cv_tab))

ropls::getScoreMN(model) %>%
  as.data.frame() %>% tibble::rownames_to_column(var = "combo") %>%
  dplyr::left_join(mat.combo[, meta_cols], by = "combo") -> scores_tab

ropls::getLoadingMN(model) %>%
  as.data.frame() %>% tibble::rownames_to_column(var = "gene") %>%
  dplyr::select(gene, p1, p2, p3) -> loadings_tab

ropls::getVipVn(model) %>%
  tibble::enframe(name = "gene", value = "vip") %>%
  dplyr::arrange(dplyr::desc(vip)) -> vip_tab
vip_top10 <- dplyr::slice_head(vip_tab, n = 10)
cat("\n[15] VIP 前 10：\n"); print(as.data.frame(vip_top10))

as.data.frame(model@modelDF) %>%
  tibble::rownames_to_column(var = "component") %>%
  dplyr::mutate(R2X_single = R2X, R2Y_single = R2Y) -> variance_tab

score_test <- function(axis_col, label) {
  d <- scores_tab; d$axis_value <- d[[axis_col]]
  aov_p <- summary(stats::aov(axis_value ~ stage, data = d))[[1]][["Pr(>F)"]][1]
  kw <- stats::kruskal.test(axis_value ~ stage, data = d)
  pair_w <- stats::pairwise.wilcox.test(d$axis_value, d$stage, p.adjust.method = "BH")
  pairs_tab <- as.data.frame(as.table(pair_w$p.value)) %>%
    dplyr::filter(!is.na(Freq)) %>%
    dplyr::transmute(contrast = paste(Var1, "vs", Var2), p_value = Freq)
  by_stage <- d %>% dplyr::group_by(stage) %>% dplyr::group_modify(~ {
    tt <- stats::t.test(.x$axis_value ~ .x$group)
    ww <- suppressWarnings(stats::wilcox.test(.x$axis_value ~ .x$group))
    dd <- (mean(.x$axis_value[.x$group == "Mix"]) - mean(.x$axis_value[.x$group == "Pure"])) /
      sqrt((stats::var(.x$axis_value[.x$group == "Pure"]) + stats::var(.x$axis_value[.x$group == "Mix"])) / 2)
    tibble::tibble(n_pure = sum(.x$group == "Pure"), n_mix = sum(.x$group == "Mix"),
                   diff_mix_minus_pure = mean(.x$axis_value[.x$group == "Mix"]) - mean(.x$axis_value[.x$group == "Pure"]),
                   t_p_value = tt$p.value, wilcoxon_p_value = ww$p.value, cohens_d = dd)
  }) %>% dplyr::ungroup()
  dplyr::bind_rows(
    tibble::tibble(axis = label, test_group = "stage_overall", comparison = "all three stages",
                   test = "one_way_ANOVA", estimate = NA_real_, p_value = aov_p, statistic = NA_real_, n = nrow(d)),
    tibble::tibble(axis = label, test_group = "stage_overall", comparison = "all three stages",
                   test = "kruskal_wallis", estimate = NA_real_, p_value = kw$p.value,
                   statistic = unname(kw$statistic), n = nrow(d)),
    pairs_tab %>% dplyr::mutate(axis = label, test_group = "stage_pairwise", test = "wilcoxon_BH",
                                estimate = NA_real_, statistic = NA_real_, n = nrow(d)) %>%
      dplyr::select(axis, test_group, comparison = contrast, test, estimate, p_value, statistic, n),
    by_stage %>%
      tidyr::pivot_longer(cols = c(t_p_value, wilcoxon_p_value), names_to = "test", values_to = "p_value") %>%
      dplyr::mutate(axis = label, test_group = "within_stage",
                    comparison = paste0(stage, ": Mix vs Pure"),
                    test = dplyr::if_else(test == "t_p_value", "welch_t", "wilcoxon")) %>%
      dplyr::select(axis, test_group, comparison, test, estimate = diff_mix_minus_pure, p_value, cohens_d, n_pure, n_mix)
  )
}

test_tab <- dplyr::bind_rows(score_test("p1", "Dim 1"), score_test("p2", "Dim 2"))

bioRtools::write_data(
  list(scores = scores_tab, loadings = loadings_tab, vip = vip_tab, vip_top10 = vip_top10,
       variance_explained = variance_tab, axis_tests = test_tab, cv_summary = cv_tab),
  file = path_out
)
cat("\n[15] 输出：", path_out, "\n")

# ---------------- Figure 2A–2E ----------------
save_panel <- function(p, pid, label, width_cm, height_cm) {
  bioRtools::write_data(p, file = file.path(out_pub_dir, paste0(pid, "_", label, ".pdf")),
                        width = width_cm, height = height_cm, units = "cm")
  bioRtools::write_data(p, file = file.path(out_pub_dir, paste0(pid, "_", label, ".png")),
                        width = width_cm, height = height_cm, units = "cm", dpi = 600)
  cat("[16] 输出：", pid, "_", label, ".pdf / .png\n", sep = "")
}

# 2A 热图（行内 0–1 缩放）
mat.heat <- as.matrix(mat.combo[, gene_cols])
rownames(mat.heat) <- mat.combo$combo
scale01_by_row <- function(v) {
  rng <- range(v, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0, length(v)))
  (v - rng[1]) / diff(rng)
}
mat.heat.scaled <- t(apply(mat.heat, 1, scale01_by_row))
stopifnot(!anyNA(mat.heat.scaled))

gene_order <- loadings_tab %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene") %>%
  dplyr::arrange(process, p1) %>% dplyr::pull(gene) %>% unique()

as.data.frame(mat.heat.scaled) %>%
  tibble::rownames_to_column(var = "combo") %>%
  tidyr::pivot_longer(cols = -combo, names_to = "gene", values_to = "scaled_expression") %>%
  dplyr::left_join(mat.combo[, meta_cols], by = "combo") %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels),
                group = factor(group, levels = c("Pure", "Mix")),
                combo_label = paste0(focal, "-", neighbor)) -> d.heat

# 投稿规格 2A：**通路水平热图**（20 条通路 × 48 个组合），生育期放在横轴分面
#   原因：15 号脚本把生育期做成行分面后，行数 = 3 个生育期 × 92 个基因 = 276 行，
#   即使整页排版也无法保证 ≥5 pt 且不叠字；改把生育期放横轴后行数 = 20 条通路，
#   在 170 × 65 mm 内即可做到 6-7 pt、零叠字。基因级热图移到附图（25c，整页）。
pathway_def_tab <- df.gene %>% dplyr::select(gene, process)
mat_gene <- as.matrix(mat.combo[, intersect(pathway_def_tab$gene, colnames(mat.combo))])
z_stage <- mat_gene
for (s in stage_levels) {
  idx <- which(mat.combo$stage == s)
  z_stage[idx, ] <- scale(mat_gene[idx, , drop = FALSE])
}
path_levels <- sort(unique(pathway_def_tab$process))
path_mat <- vapply(path_levels, function(pp) {
  gs <- intersect(pathway_def_tab$gene[pathway_def_tab$process == pp], colnames(z_stage))
  rowMeans(z_stage[, gs, drop = FALSE], na.rm = TRUE)
}, numeric(nrow(z_stage)))
colnames(path_mat) <- path_levels
stopifnot(!anyNA(path_mat))

# x 轴顺序：每个生育期内先 Pure（4 个）后 Mix（12 个）
mat.combo %>%
  dplyr::distinct(focal, neighbor, group) %>%
  dplyr::mutate(combo_label = paste0(focal, "-", neighbor)) %>%
  dplyr::arrange(group, focal, neighbor) %>%
  dplyr::pull(combo_label) -> combo_order

as.data.frame(path_mat) %>%
  # 2026-10-19 修（硬伤）：path_mat 由矩阵运算得到、**没有行名**，用 rownames_to_column() 会得到 "1".."48"，
  #   与 mat.combo$combo（如 "tillering_15_15"）类型不一致 → left_join 全部落空 → stage/focal/neighbor/group
  #   全变 NA，48 个组合被塌缩成一列（x 轴只显示 "NA"）。改为直接按行序取 mat.combo$combo。
  dplyr::mutate(combo = mat.combo$combo) %>%
  tidyr::pivot_longer(cols = -combo, names_to = "pathway", values_to = "mean_z") %>%
  dplyr::left_join(mat.combo[, meta_cols], by = "combo") %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels),
                group = factor(group, levels = c("Pure", "Mix")),
                pathway = factor(pathway, levels = rev(path_levels)),
                combo_label = factor(paste0(focal, "-", neighbor), levels = combo_order)) %>%
  ggplot2::ggplot(aes(x = combo_label, y = pathway, fill = mean_z)) +
  ggplot2::geom_tile() +
  ggplot2::facet_grid(. ~ stage, scales = "free_x", space = "free_x") +
  # 配色：本项目 five（深紫 #3B1F4E）/ hundred（深红 #780018）两端 + 白色中点，
  #      避免红绿并列（AGENTS「配色规范」要求从预定义配色取值）
  ggplot2::scale_fill_gradient2(low = "#3B1F4E", mid = "#FFFFFF", high = "#780018",
                                # 2026-10-19：图例标题旋转 90°，原 "Mean scaled expression (Z)" 超出图例
                                #   键高、顶端被切 → 缩短为 "Scaled expression (Z)"，字号 6.5 → 5.5
                                midpoint = 0, name = "Scaled expression (Z)") +
  ggplot2::labs(x = "Focal-neighbour combination (orange = mixture)", y = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
                 axis.text.y = ggplot2::element_text(size = 5.5),
                 axis.title.x = ggplot2::element_text(size = 7.5),
                 strip.text = ggplot2::element_text(size = 8, face = "bold"),
                 panel.grid = ggplot2::element_blank(),
                 legend.position = "right",
                 legend.key.width = ggplot2::unit(0.3, "cm"),
                 legend.key.height = ggplot2::unit(1.0, "cm"),
                 # 2026-10-19：图例标题改水平（原来旋转 90°，比图例键高还长，顶端被纸边切掉）
                 legend.title = ggplot2::element_text(size = 5.5, angle = 0),
                 legend.text = ggplot2::element_text(size = 6)) -> p.2a
save_panel(p.2a, "25", "Figure2A_通路水平热图_三时期_中位数矩阵", 17.0, 7.5)

# 25c：**基因级热图（92 基因 × 48 组合）**，改把生育期放横轴后行数 = 92，
#      整页（170 × 220 mm）可保证 5.5 pt 行标签不叠字 → 用作附图（Figure S18）
d.heat %>%
  dplyr::mutate(gene = factor(gene, levels = rev(gene_order))) %>%
  ggplot2::ggplot(aes(x = combo_label, y = gene, fill = scaled_expression)) +
  ggplot2::geom_tile() +
  ggplot2::facet_grid(. ~ stage, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_gradientn(colours = c("#f2f2f0", "#a8bedb", "#1e4e8c"),
                                name = "Scaled median expression (0-1)") +
  ggplot2::labs(x = "Focal-neighbour combination", y = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
                 axis.text.y = ggplot2::element_text(size = 5.2, face = "italic"),
                 axis.title.x = ggplot2::element_text(size = 7.5),
                 strip.text = ggplot2::element_text(size = 8, face = "bold"),
                 panel.grid = ggplot2::element_blank(),
                 legend.position = "right",
                 legend.key.width = ggplot2::unit(0.3, "cm"),
                 legend.key.height = ggplot2::unit(1.2, "cm"),
                 legend.title = ggplot2::element_text(size = 5.5, angle = 0),
                 legend.text = ggplot2::element_text(size = 6)) -> p.2a_full
save_panel(p.2a_full, "25c", "FigureS18_基因级热图_92基因_三时期_中位数矩阵", 17.0, 22.0)

pc_lab <- function(i) {
  v <- variance_tab$R2X_single[i]
  paste0("Dim ", i, " (", formatC(100 * v, format = "f", digits = 2), "%)")
}

scores_tab %>%
  dplyr::mutate(stage = factor(stage, levels = stage_levels),
                plotgroup = factor(paste0(stage, " ", group), levels = plotgroup_levels)) -> scores_plot

scores_plot %>%
  ggplot2::ggplot(aes(x = p1, y = p2, color = plotgroup)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::stat_ellipse(level = 0.6, type = "norm", linewidth = 0.5) +
  ggplot2::geom_point(size = 3.2) +
  bioRtools::scale_color_research() +
  ggplot2::labs(x = pc_lab(1), y = pc_lab(2), color = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "right",
                 legend.key.size = ggplot2::unit(0.3, "cm"),
                 legend.text = ggplot2::element_text(size = 6),
                 legend.margin = ggplot2::margin(0, 0, 0, 0)) -> p.2b
save_panel(p.2b, "26", "Figure2B_PLSDA得分图_三时期_中位数矩阵", 8.3, 5.0)

plot_axis <- function(axis_name, pid, label) {
  axis_col <- dplyr::if_else(axis_name == "Dim 1", "p1", "p2")
  d <- scores_plot %>% dplyr::mutate(value = .data[[axis_col]])
  overall_p <- test_tab %>% dplyr::filter(axis == axis_name, test_group == "stage_overall",
                                          test == "kruskal_wallis") %>% dplyr::pull(p_value)
  group_p <- test_tab %>%
    dplyr::filter(axis == axis_name, test_group == "within_stage", test == "wilcoxon") %>%
    dplyr::mutate(stage = stringr::str_remove(comparison, ": Mix vs Pure"),
                  stage = factor(stage, levels = stage_levels),
                  plotgroup = factor(paste0(stage, " Mix"), levels = plotgroup_levels), y = -Inf)
  d %>%
    ggplot2::ggplot(aes(x = plotgroup, y = value, fill = plotgroup)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.55, size = 1.2) +
    ggplot2::geom_text(data = group_p, aes(x = plotgroup, y = y, label = format_p_compact(p_value)),
                       vjust = -0.8, size = 2.6, inherit.aes = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.14, 0.16))) +
    ggplot2::scale_x_discrete(limits = plotgroup_levels) +
    bioRtools::scale_fill_research() +
    ggplot2::labs(x = NULL, y = paste0(axis_name, " score"),
                  subtitle = paste0("Stage effect (Kruskal-Wallis): ", format_p_compact(overall_p))) +
    bioRtools::theme_prism() +
    ggplot2::theme(legend.position = "none",
                   plot.subtitle = ggplot2::element_text(size = 7, hjust = 0),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1, size = 7)) -> p
  save_panel(p, pid, label, 8.3, 5.0)
}
plot_axis("Dim 1", "27", "Figure2C_Dim1分布_三时期_中位数矩阵")
plot_axis("Dim 2", "28", "Figure2D_Dim2分布_三时期_中位数矩阵")

loadings_tab %>%
  dplyr::filter(gene %in% vip_top10$gene) %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene") %>%
  dplyr::left_join(vip_top10, by = "gene") -> load_top

load_top %>%
  ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_segment(aes(x = 0, y = 0, xend = p1, yend = p2, color = process),
                        arrow = ggplot2::arrow(length = ggplot2::unit(0.22, "cm"), type = "closed"),
                        linewidth = 0.6) +
  ggrepel::geom_label_repel(aes(x = p1, y = p2, label = gene, color = process),
                            size = 2.4, fontface = "italic", show.legend = FALSE,
                            label.size = 0.12, label.padding = ggplot2::unit(0.06, "lines"),
                            box.padding = 0.30, min.segment.length = 0.15, force = 12, max.overlaps = 15,
                            seed = 12345) +
  bioRtools::scale_color_research() +
  ggplot2::labs(x = pc_lab(1), y = pc_lab(2), color = "Pathway") +
  # 2026-10-19 修（用户反馈"2E 的 X 轴有点窄，应该是图例的问题"）：
  #   原图例放右侧，8 条通路名把 83 mm 宽的面板挤掉约 23 mm（绘图区只剩约 47 mm 宽）；
  #   改为图例置底、3 行铺开，绘图区拿回整幅宽度；同时关掉面板 clip，避免标签被裁。
  # 2026-10-19 二次修：图例置底但只给"2 列"——3 列时图例比 83 mm 面板还宽，
  #   最右一列（Drought / Primary metabolism）与图例标题 Pathway 被纸边切掉（实测越界 0.7–0.9 pt）；
  #   改为 2 列 4 行、去掉图例标题（"colour = pathway" 写进正文图注），文字 5.8 → 5.5 pt
  ggplot2::guides(color = ggplot2::guide_legend(nrow = 4, ncol = 2, byrow = TRUE, title = NULL)) +
  ggplot2::coord_cartesian(clip = "off") +
  bioRtools::theme_prism() +
  ggplot2::theme(legend.position = "bottom",
                 legend.key.size = ggplot2::unit(0.26, "cm"),
                 legend.text = ggplot2::element_text(size = 5.5),
                 legend.title = ggplot2::element_text(size = 6),
                 legend.box.spacing = ggplot2::unit(0.08, "cm"),
                 legend.margin = ggplot2::margin(0, 0, 0, 0)) -> p.2e
# 2026-10-19 修：10 个基因标签在 76 mm 高的面板里仍有 22% 叠字 → 面板加高到 92 mm
save_panel(p.2e, "29", "Figure2E_VIP载荷图_三时期_中位数矩阵", 8.3, 9.2)

cat("\n[15] sessionInfo()\n")
sessionInfo()
