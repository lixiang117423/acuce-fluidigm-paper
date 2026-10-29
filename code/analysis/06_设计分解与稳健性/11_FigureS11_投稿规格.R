# ============================================================
# 06_设计分解与稳健性 / 04_转录组_伙伴效应方差分解.R
# 独立分析（续 06 主题，Q6）：把每个基因的表达变异拆成三部分——
#   自身基因型效应（focal）、伙伴基因型效应（neighbor）、特定组合效应（focal x neighbor）
# 为什么值得单独做：旧版与主题 04 都只比较 Mix vs Pure，无法回答"邻体的遗传身份本身贡献了多少表达变异"，
#   也无法区分"某品系自身表达高"与"跟谁种在一起才表达高"（后者即特定组合效应 / SCA）。
#
# 【方法选择说明】先用交叉随机效应模型做方差成分分解（(1|focal) + (1|neighbor) + (1|combo)），
#   结果 184 个模型中 174 个为奇异拟合（成分被压到 0），残差占比 96%–97%，
#   说明在本样本量下（每时期 80–88 个样本、16 个组合、每个因子 4 个水平）随机效应方差成分不可识别。
#   因此改用**两因素方差分析分解**（响应为 log2(expression+1)）：
#     y ~ focal * neighbor（Type III SS，sum-to-zero 对照），
#   逐基因输出三项的 partial eta^2 与 F/P，再按时期对 92 个基因做 BH-FDR。
#   该模型的残差自由度 = 样本数 − 16（tillering 64、filling 72），可给出可靠检验。
#
# 运行：/opt/homebrew/bin/Rscript "<root>/code/06_设计分解与稳健性/04_转录组_伙伴效应方差分解.R"
# 输入：<root>/result/03_基因表达总体/01_表_表达矩阵.xlsx（sample_matrix、gene_annotation）
# 输出：<root>/result/06_设计分解与稳健性/04_表_伙伴效应方差分解.xlsx
#       <root>/result/06_设计分解与稳健性/04_FigureS11_伙伴效应方差分解.pdf / .png
# 其中 <root> = /Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理
#
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录见 code/06_设计分解与稳健性/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_mat <- file.path(project_root, "result", "03_基因表达总体", "01_表_表达矩阵.xlsx")
out_dir <- file.path(project_root, "result", "06_设计分解与稳健性")
out_pub_dir <- file.path(project_root, "result", "08_投稿规格图")  # 投稿规格面板输出目录
dir.create(out_pub_dir, showWarnings = FALSE, recursive = TRUE)
path_out <- file.path(out_pub_dir, "98_重出副本_04_表_伙伴效应方差分解.xlsx")  # 投稿规格版：表写到这里，不动原目录

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stopifnot(file.exists(path_mat))

read_sheet_raw <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") %>% as.data.frame()
}

read_sheet_raw(path_mat, "sample_matrix") -> mat.sample
bioRtools::read_data(path_mat, sheet = "gene_annotation") -> df.gene

meta_cols <- c("stage", "sample", "focal", "neighbor")
gene_cols <- setdiff(colnames(mat.sample), meta_cols)
stopifnot(length(gene_cols) == 92, nrow(mat.sample) == 168)

cat("[04] 输入：", nrow(mat.sample), " 个样本；", length(gene_cols), " 个基因；",
    dplyr::n_distinct(paste0(mat.sample$focal, "-", mat.sample$neighbor)), " 个组合\n")

# ---- 单基因 × 单时期：两因素 ANOVA 分解（Type III SS）----
anova_one <- function(st, gene) {
  d <- mat.sample %>%
    dplyr::filter(stage == st) %>%
    dplyr::transmute(
      y = log2(.data[[gene]] + 1),
      focal = factor(focal),
      neighbor = factor(neighbor)
    )

  m <- stats::lm(y ~ focal * neighbor, data = d)
  a <- car::Anova(m, type = 3,
                  contrasts = list(focal = stats::contr.sum, neighbor = stats::contr.sum)) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "term")

  ss_resid <- a$`Sum Sq`[a$term == "Residuals"]
  ss_term <- function(nm) {
    v <- a$`Sum Sq`[a$term == nm]
    if (length(v) == 0) NA_real_ else v
  }

  tibble::tibble(
    stage = st,
    gene = gene,
    n_obs = nrow(d),
    ss_focal = ss_term("focal"),
    ss_neighbor = ss_term("neighbor"),
    ss_interaction = ss_term("focal:neighbor"),
    ss_residual = ss_resid,
    f_focal = a$`F value`[a$term == "focal"],
    f_neighbor = a$`F value`[a$term == "neighbor"],
    f_interaction = a$`F value`[a$term == "focal:neighbor"],
    p_focal = a$`Pr(>F)`[a$term == "focal"],
    p_neighbor = a$`Pr(>F)`[a$term == "neighbor"],
    p_interaction = a$`Pr(>F)`[a$term == "focal:neighbor"]
  ) %>%
    dplyr::mutate(
      # partial eta^2 = SS_term / (SS_term + SS_residual)
      eta2_focal = ss_focal / (ss_focal + ss_residual),
      eta2_neighbor = ss_neighbor / (ss_neighbor + ss_residual),
      eta2_interaction = ss_interaction / (ss_interaction + ss_residual),
      eta2_total_explained = 1 - ss_residual / (ss_focal + ss_neighbor + ss_interaction + ss_residual)
    )
}

grid_task <- tidyr::expand_grid(stage = c("tillering", "filling"), gene = gene_cols)
cat("[04] 待拟合模型数：", nrow(grid_task), "（2 时期 × 92 基因）\n")

eta_tab <- purrr::pmap_dfr(list(grid_task$stage, grid_task$gene), anova_one) %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(
    fdr_focal = p.adjust(p_focal, method = "BH"),
    fdr_neighbor = p.adjust(p_neighbor, method = "BH"),
    fdr_interaction = p.adjust(p_interaction, method = "BH")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(df.gene %>% dplyr::select(gene, process), by = "gene")

cat("[04] 逐基因 FDR < 0.05 的基因数：\n")
print(as.data.frame(eta_tab %>%
                      dplyr::group_by(stage) %>%
                      dplyr::summarise(
                        n_sig_focal = sum(fdr_focal < 0.05),
                        n_sig_neighbor = sum(fdr_neighbor < 0.05),
                        n_sig_interaction = sum(fdr_interaction < 0.05),
                        .groups = "drop"
                      )))

summary_tab <- eta_tab %>%
  dplyr::group_by(stage) %>%
  dplyr::summarise(
    n_gene = dplyr::n(),
    median_eta2_focal = median(eta2_focal),
    median_eta2_neighbor = median(eta2_neighbor),
    median_eta2_interaction = median(eta2_interaction),
    median_eta2_total = median(eta2_total_explained),
    q75_eta2_interaction = quantile(eta2_interaction, 0.75),
    n_neighbor_gt_focal = sum(eta2_neighbor > eta2_focal),
    .groups = "drop"
  )

cat("\n[04] partial eta^2 汇总：\n")
print(as.data.frame(summary_tab), row.names = FALSE)

pair_test <- eta_tab %>%
  dplyr::group_by(stage) %>%
  dplyr::group_modify(~ {
    wt <- suppressWarnings(stats::wilcox.test(.x$eta2_neighbor, .x$eta2_focal, paired = TRUE))
    tibble::tibble(
      median_diff_neighbor_minus_focal = median(.x$eta2_neighbor - .x$eta2_focal),
      wilcoxon_p = wt$p.value
    )
  }) %>%
  dplyr::ungroup()

cat("\n[04] 伙伴效应 vs 自身效应（配对 Wilcoxon，92 个基因）：\n")
print(as.data.frame(pair_test), row.names = FALSE)

# ---- 置换参照：打破"配对关系"后的交互项 eta^2 期望值 ----
# 说明：交互项自由度为 9，而两个主效应各为 3，直接比较 eta^2 会偏向交互项。
#   因此把 neighbor 标签在 focal 组内随机打乱（保留每个 focal 的伙伴组成、破坏具体配对），
#   重算全部基因的交互项 eta^2，得到"纯加性下的期望水平"作为参照。
perm_reference <- function(st, n_perm = 99) {
  d0 <- mat.sample %>% dplyr::filter(stage == st)

  one_perm <- function(seed) {
    set.seed(seed)
    d <- d0 %>%
      dplyr::group_by(focal) %>%
      dplyr::mutate(neighbor_perm = sample(neighbor)) %>%
      dplyr::ungroup()

    purrr::map_dbl(gene_cols, function(g) {
      dd <- d %>% dplyr::transmute(y = log2(.data[[g]] + 1),
                                   focal = factor(focal), neighbor = factor(neighbor_perm))
      a <- car::Anova(stats::lm(y ~ focal * neighbor, data = dd), type = 3,
                      contrasts = list(focal = stats::contr.sum, neighbor = stats::contr.sum)) %>%
        as.data.frame()
      ss_int <- a$`Sum Sq`[rownames(a) == "focal:neighbor"]
      ss_res <- a$`Sum Sq`[rownames(a) == "Residuals"]
      ss_int / (ss_int + ss_res)
    }) %>% stats::median()
  }

  perm_med <- purrr::map_dbl(seq_len(n_perm), one_perm)

  obs <- eta_tab %>% dplyr::filter(stage == st) %>% dplyr::pull(eta2_interaction) %>% stats::median()

  tibble::tibble(
    stage = st,
    observed_median_eta2_interaction = obs,
    null_median_mean = mean(perm_med),
    null_median_q025 = quantile(perm_med, 0.025),
    null_median_q975 = quantile(perm_med, 0.975),
    perm_p = (sum(perm_med >= obs) + 1) / (n_perm + 1)
  )
}

cat("\n[04] 置换参照（打乱配对关系后的交互项 eta^2 中位数，99 次置换）：\n")
perm_tab <- purrr::map_dfr(c("tillering", "filling"), perm_reference)
print(as.data.frame(perm_tab), row.names = FALSE)

top_interaction <- eta_tab %>%
  dplyr::group_by(stage) %>%
  dplyr::arrange(dplyr::desc(eta2_interaction), .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup() %>%
  dplyr::select(stage, gene, process, eta2_focal, eta2_neighbor, eta2_interaction,
                f_interaction, fdr_interaction)

cat("\n[04] 组合特异（focal x neighbor 交互）效应最强的 10 个基因（每时期）：\n")
print(as.data.frame(top_interaction), row.names = FALSE)

pathway_tab <- eta_tab %>%
  dplyr::group_by(stage, process) %>%
  dplyr::summarise(
    n_gene = dplyr::n(),
    mean_eta2_focal = mean(eta2_focal),
    mean_eta2_neighbor = mean(eta2_neighbor),
    mean_eta2_interaction = mean(eta2_interaction),
    n_sig_interaction = sum(fdr_interaction < 0.05),
    .groups = "drop"
  ) %>%
  dplyr::arrange(stage, dplyr::desc(mean_eta2_interaction))

cat("\n[04] 通路层面（按交互项 eta^2 排序，前 12 条）：\n")
print(as.data.frame(pathway_tab %>% dplyr::slice_head(n = 12)), row.names = FALSE)

bioRtools::write_data(
  list(gene_level = eta_tab, stage_summary = summary_tab, partner_vs_own_test = pair_test,
       permutation_reference = perm_tab,
       top_interaction_genes = top_interaction, pathway_level = pathway_tab),
  file = path_out
)
cat("[04] 输出：", path_out, "\n")

# ---- 图 ----
# 投稿规格：分面标签缩短（"Own genotype" / "Partner genotype" / "Specific pairing"），
# x 轴标签用 30° 倾斜，避免 85 mm 窄面板里文字越界/叠字
p.eta <- eta_tab %>%
  dplyr::select(stage, gene, eta2_focal, eta2_neighbor, eta2_interaction) %>%
  tidyr::pivot_longer(cols = dplyr::starts_with("eta2_"), names_to = "term", values_to = "eta2") %>%
  dplyr::mutate(
    term = factor(term,
                  levels = c("eta2_focal", "eta2_neighbor", "eta2_interaction"),
                  labels = c("Own genotype", "Partner genotype", "Specific pairing")),
    stage = factor(stage, levels = c("tillering", "filling"))
  ) %>%
  ggplot2::ggplot(aes(x = term, y = eta2, fill = term)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  ggplot2::facet_wrap(~stage) +
  bioRtools::scale_fill_bioinfo() +
  ggplot2::labs(x = NULL, y = expression("Partial " * eta^2 * " per gene"), fill = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1, size = 7),
    strip.text = ggplot2::element_text(face = "bold")
  ) -> p.eta

p.scatter <- eta_tab %>%
  dplyr::mutate(stage = factor(stage, levels = c("tillering", "filling"))) %>%
  ggplot2::ggplot(aes(x = eta2_focal, y = eta2_neighbor)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::geom_point(aes(color = eta2_interaction), size = 2, alpha = 0.85) +
  ggplot2::facet_wrap(~stage) +
  ggplot2::scale_color_gradientn(colours = c("#d9d9d6", "#a63446"),
                                 name = expression("Interaction " * eta^2)) +
  ggplot2::labs(x = expression("Own genotype " * eta^2),
                y = expression("Partner genotype " * eta^2)) +
  bioRtools::theme_prism() +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"),
                 axis.text = ggplot2::element_text(size = 6.5),
                 axis.title = ggplot2::element_text(size = 7)) -> p.scatter

p.top <- top_interaction %>%
  dplyr::mutate(stage = factor(stage, levels = c("tillering", "filling"))) %>%
  ggplot2::ggplot(aes(x = eta2_interaction, y = forcats::fct_reorder(gene, eta2_interaction),
                      fill = eta2_interaction)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::facet_wrap(~stage, scales = "free_y") +
  ggplot2::scale_fill_gradientn(colours = c("#a8bedb", "#1e4e8c"),
                                name = expression("Interaction " * eta^2)) +
  ggplot2::labs(x = expression("Interaction " * eta^2 * " (focal x neighbor)"), y = NULL) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 7, face = "italic"),
    strip.text = ggplot2::element_text(face = "bold")
  ) -> p.top

# 投稿规格：三个面板分别导出（组图交 VectorCraft）
bioRtools::write_data(p.eta, file = file.path(out_pub_dir, "04a_FigureS11_eta分布.pdf"),
                      width = 8.5, height = 7.0, units = "cm")
bioRtools::write_data(p.eta, file = file.path(out_pub_dir, "04a_FigureS11_eta分布.png"),
                      width = 8.5, height = 7.0, units = "cm", dpi = 600)
bioRtools::write_data(p.scatter, file = file.path(out_pub_dir, "04b_FigureS11_伙伴效应散点.pdf"),
                      width = 8.5, height = 7.0, units = "cm")
bioRtools::write_data(p.scatter, file = file.path(out_pub_dir, "04b_FigureS11_伙伴效应散点.png"),
                      width = 8.5, height = 7.0, units = "cm", dpi = 600)
bioRtools::write_data(p.top, file = file.path(out_pub_dir, "04c_FigureS11_基因排序.pdf"),
                      width = 17.0, height = 7.5, units = "cm")
bioRtools::write_data(p.top, file = file.path(out_pub_dir, "04c_FigureS11_基因排序.png"),
                      width = 17.0, height = 7.5, units = "cm", dpi = 600)
cat("[11] 输出：04a / 04b / 04c_FigureS11 面板（PDF + PNG）\n")

cat("\n[04] sessionInfo()\n")
sessionInfo()
