# ============================================================
# 02_农艺性状 / 06_田间株高与千粒重补充.R
# 作用：田间实验（2019–2021、3 田块、6 纯系）的株高（PH）与千粒重（KGW）补充分析
#       —— 作为主题 02 农艺性状的独立年份/独立性状验证（Figure S4）
# 运行：/opt/homebrew/bin/Rscript "<root>/code/02_农艺性状/06_田间株高与千粒重补充.R"
#
# 输入：<root>/data/01_分析输入表/02_元阳梯田田间病害调查表.xlsx
#       （Category 中 PH 3564 行、KGW 1140 行，此前旧版未分析）
# 输出：<root>/result/02_农艺性状/12_表_田间株高与千粒重.xlsx
#       <root>/result/02_农艺性状/12_FigureS4_田间株高与千粒重.pdf / .png
#
# 模型：value ~ condition + FocalPlant + (1 | Year)，分 PH / KGW 各拟合一次（condition = Pure / Mix）；
#       同时给出每个 focal 的 condition 效应（emmeans 估计值 + 95% 渐近 CI），用于判断混种是否改变田间株高与千粒重。
# 口径说明：36 / 74 两个纯系在田间仅与固定伙伴配对（36-49、69-74），故其 condition 对比同样可比；
#       本图按 6 个纯系全部展示，主结论以 4 个焦点纯系（15/49/69/87）为准。
# 坐标说明：本表不含基因组坐标，不涉及 0-based / 1-based 换算。
# 断点续传：单步耗时 < 60 s，每次运行都重新生成。
# 版本记录：脚本末尾打印 sessionInfo()；运行记录见 code/02_农艺性状/README.md
# ============================================================

rm(list = ls())

library(tidyverse)

set.seed(12345)

project_root <- "/Users/lixiang/NutstoreFiles/02.博后相关/02.data/90.fluidigm文章重新整理"
path_in <- file.path(project_root, "data", "01_分析输入表", "02_元阳梯田田间病害调查表.xlsx")
out_dir <- file.path(project_root, "result", "02_农艺性状")
path_out <- file.path(out_dir, "12_表_田间株高与千粒重.xlsx")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stopifnot(file.exists(path_in))

format_p <- function(p, digits = 3) {
  dplyr::if_else(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, format = "f", digits = digits)))
}

# emmeans 在观测数 > 3000 时会关闭自由度计算，此时置信区间列名为 asymp.LCL / asymp.UCL，
# 统一改成 lower.CL / upper.CL，便于下游绘图与选列
norm_ci <- function(x) {
  if (!"lower.CL" %in% names(x) && "asymp.LCL" %in% names(x)) {
    x <- dplyr::rename(x, lower.CL = asymp.LCL, upper.CL = asymp.UCL)
  }
  x
}

bioRtools::read_data(path_in) -> df.raw

df.raw %>%
  dplyr::filter(Category %in% c("PH", "KGW")) %>%
  dplyr::mutate(
    FocalPlant = as.character(FocalPlant),
    Neighbor = as.character(Neighbor),
    Year = as.factor(Year),
    Field = as.factor(Field),
    condition = factor(stringr::str_to_title(condition), levels = c("Pure", "Mix"))
  ) -> df

cat("[06] 输入：PH", sum(df$Category == "PH"), "行；KGW", sum(df$Category == "KGW"), "行\n")
cat("[06] PH/KGW 数据覆盖（纯系 x 年份）：\n")
print(table(df$Category, df$FocalPlant, df$Year))

# ---- 总体 condition 效应 ----
fit_overall <- function(ct) {
  d <- df %>% dplyr::filter(Category == ct, !is.na(Value)) %>% droplevels()
  m <- lmerTest::lmer(Value ~ condition + FocalPlant + (1 | Year), data = d)
  s <- summary(m)$coefficients
  tibble::tibble(
    category = ct,
    n_obs = nrow(d),
    estimate = s["conditionMix", "Estimate"],
    std_error = s["conditionMix", "Std. Error"],
    df = s["conditionMix", "df"],
    t_value = s["conditionMix", "t value"],
    p_value = s["conditionMix", "Pr(>|t|)"],
    mean_pure = mean(d$Value[d$condition == "Pure"]),
    mean_mix = mean(d$Value[d$condition == "Mix"])
  )
}

overall_tab <- purrr::map_dfr(c("PH", "KGW"), fit_overall)

cat("\n[06] 总体 condition 效应（LMM）：\n")
print(as.data.frame(overall_tab))

# ---- 每个纯系的 condition 均值与 CI ----
emm_tab <- purrr::map_dfr(c("PH", "KGW"), function(ct) {
  d <- df %>% dplyr::filter(Category == ct, !is.na(Value)) %>% droplevels()
  m <- lmerTest::lmer(Value ~ condition + FocalPlant + (1 | Year), data = d)
  emmeans::emmeans(m, ~ condition | FocalPlant, lmer.df = "satterthwaite") %>%
    as.data.frame() %>%
    norm_ci() %>%
    dplyr::mutate(category = ct) %>%
    dplyr::select(category, focal = FocalPlant, condition, emmean, SE, df, lower.CL, upper.CL)
})

cat("\n[06] 各纯系 condition 均值：\n")
print(as.data.frame(emm_tab))

# ---- 描述性统计 ----
df %>%
  dplyr::filter(!is.na(Value)) %>%
  dplyr::group_by(Category, FocalPlant, condition, Year) %>%
  dplyr::summarise(n = dplyr::n(), mean = mean(Value), sd = sd(Value), .groups = "drop") -> desc_tab

bioRtools::write_data(
  list(overall_model = overall_tab, emmeans = emm_tab, descriptive = desc_tab),
  file = path_out
)
cat("[06] 输出：", path_out, "\n")

# ---- Figure S4 ----
emm_tab %>%
  dplyr::mutate(
    category = factor(category, levels = c("PH", "KGW"),
                      labels = c("Plant height (cm)", "1000-grain weight (g)")),
    focal = factor(focal, levels = c("15", "36", "49", "69", "74", "87"))
  ) %>%
  ggplot2::ggplot(aes(x = focal, y = emmean, color = condition)) +
  ggplot2::geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.18,
                         position = ggplot2::position_dodge(width = 0.5)) +
  ggplot2::geom_point(size = 2.6, position = ggplot2::position_dodge(width = 0.5)) +
  ggplot2::geom_text(
    aes(y = upper.CL, label = focal), alpha = 0, show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~category, scales = "free_y") +
  bioRtools::scale_color_research() +
  ggplot2::labs(
    x = "Focal accession",
    y = "Estimated mean (95% CI)",
    color = "Planting mode",
    caption = paste0(
      "Field trial 2019-2021 (3 years x 3 fields). Overall condition effect: ",
      paste0(overall_tab$category, " ", format_p(overall_tab$p_value), collapse = "; ")
    )
  ) +
  bioRtools::theme_prism() +
  ggplot2::theme(
    legend.position = "bottom",
    strip.text = ggplot2::element_text(face = "bold"),
    plot.caption = ggplot2::element_text(size = 8, hjust = 0)
  ) -> p.s4

bioRtools::write_data(p.s4, file = file.path(out_dir, "12_FigureS4_田间株高与千粒重.pdf"),
                      width = 10, height = 5, units = "in")
bioRtools::write_data(p.s4, file = file.path(out_dir, "12_FigureS4_田间株高与千粒重.png"),
                      width = 10, height = 5, units = "in", dpi = 500)
cat("[06] 输出：12_FigureS4_田间株高与千粒重.pdf / .png\n")

cat("\n[06] sessionInfo()\n")
sessionInfo()
