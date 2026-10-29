# Source Data v45

文件：`Source_Data_v45.xlsx`

本版按第二轮审稿最终意见调整：

- Figure S2A–S2F 改为盆均值数据；S2C–S2F 的字母标记由盆均值模型计算。
- 删除原根系分泌物激素 Figure S6、S7 对应 sheet。
- 删除与主图 Figure 5 重复的原 Figure S14 sheet。
- 删除上述图后，后续补充图连续重编号。
- 共 47 个 sheet：正文 Fig1–Fig5 的 17 个子图 + 补充图 30 个 sheet。

## 补充图重编号

| 新编号 | 旧编号 | sheet | 绘图代码（保留原文件名） |
| --- | --- | --- | --- |
| Figure S2 | Figure S2 | FigS2a–FigS2h | `17_FigureS2_盆级重出.R`；G/H 沿用 `FigS2.R` |
| Figure S3 | Figure S3 | FigS3 | `FigS3.R` |
| Figure S4 | Figure S4 | FigS4 | `FigS4.R` |
| Figure S5 | Figure S5 | FigS5 | `FigS5.R` |
| Figure S6 | Figure S8 | FigS6a、FigS6b | `FigS8.R` |
| Figure S7 | Figure S9 | FigS7a、FigS7b | `FigS9.R` |
| Figure S8 | Figure S10 | FigS8 | `FigS10.R` |
| Figure S9 | Figure S11 | FigS9a–FigS9c | `FigS11.R` |
| Figure S10 | Figure S12 | FigS10a–FigS10c | `FigS12.R` |
| Figure S11 | Figure S13 | FigS11 | `FigS13.R` |
| Figure S12 | Figure S15 | FigS12a、FigS12b | `FigS15.R` |
| Figure S13 | Figure S16 | FigS13 | `FigS16.R` |
| Figure S14 | Figure S17 | FigS14 | `FigS17.R` |
| Figure S15 | Figure S18 | FigS15 | `FigS18.R` |
| Figure S16 | Figure S19 | FigS16a、FigS16b | `FigS19.R` |

Figure S1 为示意图，无 Source Data。

## Figure S2 的统计口径

- 推断单元：盆，`image` 为盆标识。
- S2A/S2B：`lm(pot mean lesion density ~ group + batch)`。
- S2C–S2F：逐 disease × focal 拟合
  `lm(pot mean lesion density ~ neighbour + batch)`。
- 多重比较：Sidak。
- `FigS2a`、`FigS2b` 的元信息中记录 raw 与 standardised model 的 P 值和盆数。
- `FigS2c`–`FigS2f` 均含 `emmean` 与 `signif` 列。

## 移除的数据

- 旧 FigS6a、FigS6b、FigS7：根系分泌物激素面板，另行成文。
- 旧 FigS14：与主图 Figure 5 的逐年/交互结果重复。

旧 `Source_Data_v32.xlsx` 与旧图文件均保留，不回写。
