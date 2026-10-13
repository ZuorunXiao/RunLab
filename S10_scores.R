#===============================================================
## S10_scores.R —— 炎症评分与分层（步骤规格 §8：S10）
##
## 目的：由层 2 的样本级段分（ssGSEA 评分矩阵）构造四个形态学评分（§4.3）：
##       IB（S1-S3 均值）、RI（S4 - S2）、SASP（Z 标准化）、IRh（五段秩向量），
##       并按部位 / 病变阶段 / IB 分层比较。
##
## 定位：分析脚本。评分定义与随机种子固定；§9 Q10 要求固定种子重算两次一致。
## 产出：03_results/S10/table/inflam_scores.csv
##       03_results/S10/table/stratified_comparison.csv
##       03_results/S10/table/q10_seed_reproducibility.csv
##       03_results/S10/Figure/Score_Stratified_Boxplot.pdf    分层箱线 + 显著性标注
##       03_results/S10/Figure/Segment_Profile_Heatmap.pdf      样本 x 五段 轮廓热图
##       03_results/S10/Figure/Segment_Profile_By_Stratum.pdf   段均轮廓线 + 95%CI
##       03_results/S10/Figure/Score_Correlation.pdf            评分相关热图（构造效度）
##
## 图件清单（规则 35 可视化强制 · SCI 发表导向；动笔前已调研文献确定风格/形式/数量）：
# S10 共 4 类（原 2 类：均值柱图冗余、缺五段轮廓、缺显著性、IRh 未可视化）：
##     1. Score_Stratified_Boxplot.pdf
# 形式：分面箱线（IQR/中位线/1.5×IQR 须线/离群点，大 n 组不叠加全样本抖动、不叠加均值菱形）+
# 显著性标注；
##              箱线按 stratum level 填 Set2 色板（回应「箱线图不加颜色」）
# （facet_grid(score ~ strata)；region/stage/ib 实测均为 2 水平，全部 Wilcoxon；Holm 校正）
##        科学问题：IB/RI/SASP 三标量评分在 region/stage/IB 三层下的分布差异与显著性。
##        参照文献：PMID 37270717 (pancreatic, 2023) Fig 9 箱线 + Mann-Whitney 显著性
##              （ns/*p<.05/**p<.01/***p<.001；IQR+中位+1.5×IQR 须线+离群点）。
##     2. Segment_Profile_Heatmap.pdf
# 形式：ComplexHeatmap 重写 —— 样本(行)×五段(列) Z-score 热图，分歧色 circize::colorRamp2(±2 封顶)；
# 行注释 stage/region/ib 色条（which=row）+ 行按 stage 拆分块（替代 facet，块高按样本量 220:26
# 比例，无白块）
##        科学问题：五段 ssGSEA 签名在样本层面的异质性与分层带结构（框架核心）。
##        参照文献：PMID 39865865 (HCC 2025) Fig 2A；PMID 38218914 (OA 2024) Fig 1B；
##              PMID 33584801 (GBM 2021) Fig 1 —— 样本×基因集 ssGSEA 热图按队列分组标注。
##     3. Segment_Profile_By_Stratum.pdf
##        形式：x=S1-S5(有序) 折线 + 95%CI 误差棒，按分层变量 facet_wrap，折线=组。
##        科学问题：initiation→chronic 五段轨迹形态如何随 stratum 改变。
##        参照文献：PMC12403291 (retinoblastoma 2025) GSVA 评分按分组折线/柱展示均值轨迹。
##     4. Score_Correlation.pdf
##        形式：{IB,RI,SASP,S1..S5} Spearman 相关热图 + 数值标注。
##        科学问题：标量分与五段轴的构造效度（IB≈mean(S1-S3)/RI≈S4-S2/SASP≈S5）。
##        注：IRh 五段秩向量由图 2/3 段轮廓涵盖，不单列。
##        参照文献：PMID 37270717 Fig 4 评分间相关展示惯例。
##
## 图件口径：sci_theme / sci_save 逐字取自规范 §sci_theme / §sci_save 官方定义；
##       分面一律 facet_wrap / facet_grid，不引入任何拼图包；显著性用预计算 p 值 + geom_text，
##       无新依赖；可读性核验见 sci_save 之 RAPA_QC_PNG 钩子。
##
## 前置：S06。段分评分矩阵来自 S06 的 segment_score_matrix.rds（段 x 样本）。
#===============================================================

## PART 0 常量 ------------------------------------------------------------------------------

SEED <- 20260925
set.seed(SEED)
FIG_DIR <- "03_results/S10/Figure"
TAB_DIR <- "03_results/S10/table"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

vec.segment <- c("S1_initiation", "S2_amplification", "S3_effect",
  "S4_resolution", "S5_chronic")
strata.vars <- c("region", "stage", "ib")   ## §8 S10：部位 / 病变阶段 / IB 分层
SCORES <- c("IB", "RI", "SASP")

library(ggplot2)  ## 评分箱线图；主题与出口用规范 §sci_theme / §sci_save
## F2 样本×五段 轮廓热图（ComplexHeatmap 2.28.0，本机已装；行注释 + 分歧色）
library(ComplexHeatmap)

## 主题：规范 §sci_theme 官方定义（theme_bw + 关网格），全脚本统一口径
sci_theme <- ggplot2::theme_bw(base_size = 12) +
ggplot2::theme(
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank(),
  strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
  legend.position  = "right",
  legend.text      = ggplot2::element_text(size = 10, face = "plain"),
  axis.title       = ggplot2::element_text(size = 12, face = "plain"),
  axis.text        = ggplot2::element_text(size = 10, face = "plain"),
  plot.title       = ggplot2::element_text(size = 14, hjust = 0.5,
    face = "plain"),
  legend.title     = ggplot2::element_text(size = 11, face = "plain")
)

## 出口：规范 §sci_save 官方定义，仅 PDF，落本步 FIG_DIR
sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  pdf(pdf_file, width = width, height = height)
  print(plot_obj)
  dev.off()
  ## 可读性核验钩子（规则 35）：默认关闭，交付物仍为纯 PDF；
  ## 置 RAPA_QC_PNG=1 时额外导 _qc/ PNG 供逐张目视核验，不改变 PDF 产出。
  if (Sys.getenv("RAPA_QC_PNG") == "1") {
    dir.create(file.path(FIG_DIR, "_qc"), showWarnings = FALSE)
    png(file.path(FIG_DIR, "_qc", paste0(filename, ".png")),
      width = width * 150, height = height * 150, res = 150)
    print(plot_obj)
    dev.off()
  }
}

# PART 1 段分评分矩阵 -> 四个评分 ------------------------------------------------------------------

list.score <- readRDS("03_results/S06/table/segment_score_matrix.rds")
acc <- "GSE232911"
mat.score <- as.matrix(list.score[[acc]])
mat.score <- mat.score[vec.segment, , drop = FALSE]
stopifnot(all(vec.segment %in% rownames(mat.score)))

dat.meta <- read.delim("03_results/S01/table/gse232911_grouping.tsv",
  stringsAsFactors = FALSE)
dat.meta <- dat.meta[match(colnames(mat.score), dat.meta$sample_id),
  c("sample_id", "tunica", "thrombus", "state", "subject"),
  drop = FALSE]
## 对齐失败的代价是静默全 NA，此处硬断言（实测 246/246 命中）。
stopifnot(!any(is.na(dat.meta$sample_id)))
rownames(dat.meta) <- dat.meta$sample_id

dat.score <- data.frame(sample = colnames(mat.score), stringsAsFactors = FALSE)
dat.score$IB <- colMeans(mat.score[vec.segment[1:3], , drop = FALSE])
dat.score$RI <- mat.score[vec.segment[4L], ] - mat.score[vec.segment[2L], ]
dat.score$SASP <- as.numeric(scale(t(mat.score[vec.segment[5L], ,
        drop = FALSE])))
mat.rank <- apply(mat.score, 2L, rank)
dat.score[, paste0("IRh_", vec.segment)] <- t(mat.rank)
dat.score$region <- dat.meta$tunica[match(dat.score$sample, rownames(dat.meta))]
dat.score$stage <- ifelse(
  dat.meta$state[match(dat.score$sample,
      rownames(dat.meta))] == "Donor abdominal aorta",
  "control", "AAA")
write.csv(dat.score, file.path(TAB_DIR, "inflam_scores.csv"), row.names = FALSE)

## PART 2 分层比较 ------------------------------------------------------------------------------

dat.score$ib <- ifelse(dat.score$IB >= median(dat.score$IB, na.rm = TRUE),
  "IB_high", "IB_low")
dat.strat <- data.frame(score = character(0), strata = character(0),
  layer = character(0),
  n = integer(0), mean = numeric(0), median = numeric(0),
  stringsAsFactors = FALSE)

for (sc in SCORES) {
  for (st in strata.vars) {
    dat.sub <- dat.score[!is.na(dat.score[[st]]), , drop = FALSE]
    if (length(unique(dat.sub[[st]])) < 2L) next
    tab <- split(dat.sub[[sc]], dat.sub[[st]])
    for (g in names(tab)) {
      v <- as.numeric(tab[[g]])
      dat.strat <- rbind(dat.strat, data.frame(
          score = sc, strata = st, layer = g, n = length(v),
          mean = mean(v), median = median(v), stringsAsFactors = FALSE))
    }
  }
}
write.csv(dat.strat, file.path(TAB_DIR, "stratified_comparison.csv"),
  row.names = FALSE)

## PART 2.5 图件（规则 35 重做：4 图 SCI 导向；文献依据见头部图件清单）---------------

## 统一长表：样本级 评分 + 分层变量（供 F1 箱线 + 显著性）
dat.score.long <- do.call(rbind, lapply(strata.vars, function(st) {
      d <- dat.score[!is.na(dat.score[[st]]), , drop = FALSE]
      data.frame(score = rep(SCORES, each = nrow(d)),
        strata = st,
        layer = rep(as.character(d[[st]]), times = length(SCORES)),
        value = as.numeric(unlist(lapply(SCORES, function(sc) d[[sc]]))),
        stringsAsFactors = FALSE)
    }))
dat.score.long$score <- factor(dat.score.long$score, levels = SCORES)
dat.score.long$strata <- factor(dat.score.long$strata, levels = strata.vars)
dat.score.long <- dat.score.long[!is.na(dat.score.long$value), ]

## 显著性预计算（无新依赖；2 组 Wilcoxon 秩和，>2 组 Kruskal-Wallis 整体；Holm 校正）
sig.rows <- list()
for (sc in SCORES) {
  for (st in strata.vars) {
    d <- dat.score[!is.na(dat.score[[st]]), , drop = FALSE]
    grp <- unique(as.character(d[[st]]))
    if (length(grp) < 2L) next
    if (length(grp) == 2L) {
      p <- wilcox.test(d[[sc]][d[[st]] == grp[1L]],
        d[[sc]][d[[st]] == grp[2L]])$p.value
      meth <- "Wilcoxon"
    } else {
      p <- kruskal.test(d[[sc]], d[[st]])$p.value
      meth <- "Kruskal-Wallis"
    }
    sig.rows[[length(sig.rows) + 1L]] <- data.frame(score = sc, strata = st,
      method = meth, p = p, stringsAsFactors = FALSE)
  }
}
dat.sig <- do.call(rbind, sig.rows)
dat.sig$p_adj <- p.adjust(dat.sig$p, method = "holm")
dat.sig$label <- ifelse(dat.sig$p_adj < 0.001, "***",
  ifelse(dat.sig$p_adj < 0.01, "**",
    ifelse(dat.sig$p_adj < 0.05, "*", "ns")))
dat.sig$score <- factor(dat.sig$score, levels = SCORES)
dat.sig$strata <- factor(dat.sig$strata, levels = strata.vars)

## 显著性括号（手画，零新依赖；复刻 ggsignif/ggpubr bracket 样式：
##   括号跨两箱、p 值居中于括号上方；依据 PMC5512552 / PMID 38170664）
bracket.rows <- list()
for (sc in SCORES) for (st in strata.vars) {
  sub <- dat.score.long[as.character(dat.score.long$score) == sc &
    as.character(dat.score.long$strata) == st, , drop = FALSE]
  if (nrow(sub) < 2L) next
  ymax <- max(sub$value, na.rm = TRUE)
  ymin <- min(sub$value, na.rm = TRUE)
  rng <- ymax - ymin
  if (rng <= 0) rng <- max(1e-6, abs(ymax))
  yb <- ymax + 0.10 * rng
  tick <- 0.04 * rng
  sig <- dat.sig[dat.sig$score == sc & dat.sig$strata == st, , drop = FALSE]
  lab <- if (nrow(sig) == 1L)
    paste0(sig$method, " p=", signif(sig$p_adj, 2), " ", sig$label) else ""
  bracket.rows[[length(bracket.rows) + 1L]] <- data.frame(
    score = sc, strata = st, x1 = 1, x2 = 2, yb = yb, tick = tick,
    label_y = yb + 0.02 * rng, label = lab, stringsAsFactors = FALSE)
}
dat.bracket <- do.call(rbind, bracket.rows)
dat.bracket$score <- factor(dat.bracket$score, levels = SCORES)
dat.bracket$strata <- factor(dat.bracket$strata, levels = strata.vars)
seg.rows <- list()
for (i in seq_len(nrow(dat.bracket))) {
  b <- dat.bracket[i, ]
  seg.rows[[length(seg.rows) + 1L]] <- rbind(
    data.frame(score = b$score, strata = b$strata,
      x = b$x1, xend = b$x2, y = b$yb, yend = b$yb),
    data.frame(score = b$score, strata = b$strata,
      x = b$x1, xend = b$x1, y = b$yb, yend = b$yb - b$tick),
    data.frame(score = b$score, strata = b$strata,
      x = b$x2, xend = b$x2, y = b$yb, yend = b$yb - b$tick))
}
dat.seg <- do.call(rbind, seg.rows)
dat.seg$score <- factor(dat.seg$score, levels = SCORES)
dat.seg$strata <- factor(dat.seg$strata, levels = strata.vars)

## F1 分层箱线 + 显著性标注（PMID 37270717：箱线 IQR/中位/1.5×IQR 须线/离群点；
##   原文 "Outliers are depicted as black dots" —— 大 n 组（100-220 点）不叠加
# 全样本抖动，否则蜂群淹没箱体（用户裁定"完全不是正常箱线图"）；不叠加均值菱形（评分近常量无信息增量）
##   显著性 = 手画 bracket 跨两箱、p 居中（ggsignif 样式；PMC5512552）
sci_save(ggplot(dat.score.long, aes(x = layer, y = value)) +
  geom_boxplot(aes(fill = layer), colour = "black", outlier.shape = 16,
    outlier.size = 0.8, outlier.alpha = 0.5, outlier.colour = "black") +
  scale_fill_brewer(palette = "Set2", name = "Stratum level") +
  geom_segment(data = dat.seg, inherit.aes = FALSE,
    aes(x = x, xend = xend, y = y, yend = yend),
    colour = "black", linewidth = 0.3) +
  geom_text(data = dat.bracket, inherit.aes = FALSE,
    aes(x = 1.5, y = label_y, label = label),
    hjust = 0.5, vjust = 0, size = 3, colour = "black") +
  facet_grid(score ~ strata, scales = "free") +
  ## y 上方留 15% 空间，防显著性标签压住箱体/离群点
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = "Inflammation scores stratified by region, stage and IB level",
    x = "Stratum", y = "Score") +
  sci_theme +
  theme(legend.position = "none",
    strip.text = element_text(size = 10, face = "plain")),
  "Score_Stratified_Boxplot", width = 12, height = 7)

## F2 五段轮廓热图（ComplexHeatmap 2.28.0 重写；回应「箱线图不加颜色/全部白框」质疑）
##   样本(行, 按 stage/region/ib 排序) × 五段(S1-S5, 列) Z-score 热图
##   行注释 = stage/region/ib 颜色条（which="row" + left_annotation）；行按 stage 拆分块
##   （替代 ggplot facet，块高按样本量比例：AAA 220 : control 26，无白块）
##   ComplexHeatmap 默认无 cell border（消白格线）；分歧色 circize::colorRamp2 断点 ±2 封顶增饱和
##   文献依据：PMID 39865865 (HCC 2025) Fig 2A；PMID 38218914 (OA 2024) Fig 1B；
##             PMID 33584801 (GBM 2021) Fig 1 —— 样本×基因集 ssGSEA 热图按队列分组标注
## 列 Z 标准化（每段跨 246 样本）：原始 ssGSEA 值域窄（约 2.0±0.06），midpoint=0
##   双色渐变全偏红、样本异质性不可见；文献口径 normalized scores
##   （PMID 39865865 "normalized ssGSEA scores" / PMID 38218914 Fig 1B）。
## apply(MARGIN=1) 已得 样本(行) × 段(列)：每段跨 246 样本做 Z 标准化，勿再 t()
mat.z <- apply(mat.score[vec.segment, , drop = FALSE], 1L, function(x) {
    (x - mean(x)) / sd(x)
  })
## ±2 封顶增饱和（避免极端值外推出超界色）
mat.z[mat.z > 2] <- 2
mat.z[mat.z < -2] <- -2
## 行序：stage -> region -> ib（块高按样本量比例）
ord <- order(dat.score$stage, dat.score$region, dat.score$ib)
mat.z <- mat.z[ord, , drop = FALSE]
stage.vec <- dat.score$stage[ord]
region.vec <- dat.score$region[ord]
ib.vec <- dat.score$ib[ord]
ht.ann <- ComplexHeatmap::HeatmapAnnotation(stage = stage.vec,
  region = region.vec, ib = ib.vec, which = "row",
  col = list(
    stage  = c("AAA" = "#E41A1C", "control" = "#377EB8"),
    region = c("Adventitia" = "#1B9E77", "Media" = "#D95F02"),
    ib     = c("IB_high" = "#984EA3", "IB_low" = "#FF7F00")))
ht.f2 <- ComplexHeatmap::Heatmap(mat.z,
  name = "Z-score",
  col = circlize::colorRamp2(c(-2, 0, 2),
    c("#2166AC", "#F7F7F7", "#B2182B")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = stage.vec,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 90,
  column_names_side = "bottom",
  column_title = "Five-segment ssGSEA profile by sample (Z-score per segment)",
  column_title_gp = grid::gpar(fontsize = 14),
  left_annotation = ht.ann,
  heatmap_legend_param = list(title = "Z-score",
    legend_height = grid::unit(4, "cm")))
sci_save(ht.f2, "Segment_Profile_Heatmap", width = 9, height = 11)

## F3 段均轮廓线 + 95%CI（x=S1-S5 有序；PMC12403291 GSVA 均值轨迹惯例）
prof.rows <- list()
for (st in strata.vars) {
  d <- dat.score[!is.na(dat.score[[st]]), , drop = FALSE]
  for (g in unique(as.character(d[[st]]))) {
    id.samp <- d$sample[d[[st]] == g]
    for (s in vec.segment) {
      v <- as.numeric(mat.score[s, id.samp])
      n <- length(v)
      prof.rows[[length(prof.rows) + 1L]] <- data.frame(strata = st, layer = g,
        segment = s, mean = mean(v),
        lo = mean(v) - 1.96 * sd(v) / sqrt(n),
        hi = mean(v) + 1.96 * sd(v) / sqrt(n),
        stringsAsFactors = FALSE)
    }
  }
}
dat.prof <- do.call(rbind, prof.rows)
dat.prof$segment <- factor(dat.prof$segment, levels = vec.segment)
dat.prof$strata <- factor(dat.prof$strata, levels = strata.vars)
sci_save(ggplot(dat.prof, aes(x = segment, y = mean, colour = layer,
    group = layer)) +
  geom_line(linewidth = 0.8) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.1, linewidth = 0.5) +
  geom_point(size = 1.5) +
  facet_wrap(~ strata, scales = "free_y") +
  labs(title = "Mean five-segment profile by stratification variable",
    x = "Segment", y = "Mean ssGSEA score", colour = "Stratum") +
  sci_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text = element_text(size = 10, face = "plain")),
  "Segment_Profile_By_Stratum", width = 11, height = 4.5)

## F4 评分相关热图（构造效度；IRh 秩向量由 F2/F3 段轮廓涵盖）
##   段分(S1-S5) 存于 mat.score 矩阵、未作为 dat.score 列，须先从矩阵提取再拼接标量分。
dat.cor <- as.data.frame(t(mat.score))
dat.cor$IB <- dat.score$IB[match(rownames(dat.cor), dat.score$sample)]
dat.cor$RI <- dat.score$RI[match(rownames(dat.cor), dat.score$sample)]
dat.cor$SASP <- dat.score$SASP[match(rownames(dat.cor), dat.score$sample)]
cor.vars <- c("IB", "RI", "SASP", vec.segment)
cor.mat <- cor(dat.cor[, cor.vars], method = "spearman", use = "complete.obs")
cor.long <- do.call(rbind, lapply(rownames(cor.mat), function(r)
  do.call(rbind, lapply(colnames(cor.mat), function(c)
    data.frame(row = r, col = c, r = cor.mat[r, c], stringsAsFactors = FALSE))))
)
cor.long$row <- factor(cor.long$row, levels = cor.vars)
cor.long$col <- factor(cor.long$col, levels = cor.vars)
sci_save(ggplot(cor.long, aes(x = col, y = row, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-1, 1), name = "Spearman r") +
  geom_text(aes(label = round(r, 2)), size = 3, colour = "black") +
  labs(title = "Score construct validity: scalar scores vs five segments",
    x = "Score / segment", y = "Score / segment", fill = "Spearman r") +
  sci_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text = element_text(size = 10, face = "plain")),
  "Score_Correlation", width = 7, height = 6)

# PART 3 Q10：固定种子重算一致性 -------------------------------------------------------------------

set.seed(SEED)
mat.re <- as.matrix(list.score[[acc]])
mat.re <- mat.re[vec.segment, , drop = FALSE]
ib.re <- colMeans(mat.re[vec.segment[1:3], , drop = FALSE])
sasp.re <- as.numeric(scale(t(mat.re[vec.segment[5L], , drop = FALSE])))
## 比较只针对数值：colMeans 带样本名 names 而 data.frame 列赋值会剥离它，
##   identical() 会因此判 FALSE（假阴性），故两侧都 as.numeric() 剥属性。
dat.q10 <- data.frame(item = c("IB", "SASP"), identical = c(
    identical(as.numeric(dat.score$IB), as.numeric(ib.re)),
    identical(as.numeric(dat.score$SASP), as.numeric(sasp.re))),
  stringsAsFactors = FALSE)
write.csv(dat.q10, file.path(TAB_DIR, "q10_seed_reproducibility.csv"),
  row.names = FALSE)

## PART 4 收口 -------------------------------------------------------------------------------

message(paste0("\n===== S10 score summary =====\nsamples scored = ", nrow(
      dat.score),
    "  stratification rows = ", nrow(dat.strat),
    "\nQ10 identical = ", paste(dat.q10$item, dat.q10$identical, sep = ":",
      collapse = "  ")))
if (any(!dat.q10$identical)) {
  stop("S10 FAILED: Q10 — scores not reproducible under the fixed seed.")
}
message("S10 OK. tables + figures -> 03_results/S10/")
