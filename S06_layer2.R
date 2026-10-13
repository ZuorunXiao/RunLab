#===============================================================
## S06_layer2.R —— 层 2：样本级五段评分（步骤规格 §8：S06）
##
## 目的：把 S02 锁定的五段基因集按方案 §4.2「**段内并集先合并**」各自并成一个段基因集，
##       用 ssGSEA 打到每个样本，得到 样本 x 段（5 行）的评分矩阵，
##       再按 §5 的设计做组间对比，输出方向 + 效应量（Cohen's d）+ BH 校正。
##
## 定位：分析脚本。层 2 = 样本级段分，层 3（跨集方向）在 S07。
##
##
## 产出：03_results/S06/table/layer2_segment_<acc>.csv   单系列段分 + 组间判定（5 行）
##       03_results/S06/table/layer2_segment_all.csv     全系列汇总（供 S07 元分析 / S13 对账）
##       03_results/S06/table/segment_score_<acc>.rds    原始 ssGSEA 段分矩阵（段 x 样本）
##       03_results/S06/table/segment_score_matrix.rds   全系列 段分矩阵（供 S08/S10）
##       03_results/S06/table/segment_union_genes.csv    五段并集成员基因（段定义单一事实源）
##       03_results/S06/Figure/Layer2_Segment_FoldChange_By_Series.pdf 热图：段 x 系列 log2FC
##       03_results/S06/Figure/Layer2_CohensD_By_Series.pdf            分组柱：段 x 系列 Cohen's d
##
## 图件清单（规则 35 可视化强制 · SCI 发表导向；动笔前已按方案 §4.2 与已发表惯例定风格）：
##   可视化: S06 共 2 幅
##     1. Layer2_Segment_FoldChange_By_Series.pdf
##        —— 形式：单面板热图（x = 五段，y = 系列，fill = log2FC 双向色；格内标数值；
##           显著格加黑色描边，未达判读门槛者数值用常规字重）
##        —— 科学问题：五段段分在各系列的组间方向与强度（层 2 样本级段分读数）。
# —— 参照：fold-change 热图为差异表达矩阵标准画法；与本项目 S05 图 1 同形，保证同项目风格一致。
##        —— 为何不用分面柱：4 系列 x 5 段 = 20 格，分面柱会拆成 4 个面板且无法一眼横向比较。
##     2. Layer2_CohensD_By_Series.pdf
##        —— 形式：单面板分组柱（x = 五段，y = Cohen's d，fill = 系列；d = 0 参考虚线）
##        —— 科学问题：各段效应的效应量（不只看 p 值；§4.2 要求并列报告 Cohen's d）。
##        —— 参照：效应量与 p 值并列报告为差异分析规范（Cohen 1988）。
##
## 图件口径：sci_theme / sci_save 逐字取自规范 §sci_theme / §sci_save 官方定义；
##       不引入任何拼图包；可读性核验见 sci_save 之 RAPA_QC_PNG 钩子。
##
## 前置：S02（基因集表）+ S05（contrast_spec.csv = 设计公式 + 显式对比；design_used_<acc>.csv
##       含 sample_id 的实际设计矩阵）。**本脚本不臆造分组**：公式与对比一律取自 S05 落盘产物。
#===============================================================

## PART 0 常量 ------------------------------------------------------------------------------

SEED <- 20260925
set.seed(SEED)
## 常量（§10 预注册）
## 注：层 1 的 |log2FC| >= 0.585 倍数阈**不在本步**。层 2 读数是 ssGSEA 评分（normalize=TRUE
## 后按数据集内极差缩放），量纲不是表达倍数；§4.2 对层 2 只预注册「显著且方向一致（同号）」。
## 依据实测：01_data/verify/S06_run_2026-09-29b.log（20/20 被 0.585 抹平，其中 12/20 BH<0.05）。
PVAL_CUTOFF <- 0.05
MIN_GS_SIZE <- 20L  ## §0.2：官方默认 minSize=1，本项目按 §10 预注册改为 20
ALPHA_SSGSEA <- 0.25  ## 官方 ssgseaParam 的 alpha 默认 0.25，与 ssGSEA 原文一致
## §4.2：段间 Jaccard > 0.5 时该段对联合贡献 ×0.5（S10 应用）
JACCARD_CUT <- 0.5

vec.segment <- c("S1_initiation", "S2_amplification", "S3_effect",
  "S4_resolution", "S5_chronic")

## 层 2 的适用数据集 = S05 的对比注册表登记的系列（单一事实源，不在本脚本另立清单）
F_SPEC <- file.path("03_results", "S05", "table", "contrast_spec.csv")
F_L1 <- file.path("03_results", "S05", "table", "layer1_gene_dir.csv")
F_GS <- file.path("03_results", "S02", "table",
  "inflammation_geneset_final.csv")
F_OV <- file.path("03_results", "S02", "table", "segment_overlap.csv")

TAB_DIR <- "03_results/S06/table"
FIG_DIR <- "03_results/S06/Figure"
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

library(GSVA)  ## 层 2 样本级段分 ssGSEA（官方 ssgseaParam 实现）
library(limma)                            ## 层 2 设计矩阵与对比，口径对齐 S05
library(ggplot2)  ## 层 2 方向图；主题与出口用规范 §sci_theme / §sci_save

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

## 出口：规范 §sci_save 官方定义，仅 PDF，落本步 FIG_DIR。
## 可读性核验钩子（规则 22 · 交付前须核验可读性）：仅当 RAPA_QC_PNG=1 时额外导 PNG
## 到 FIG_DIR/_qc/ 供核验；默认不产出，交付物仍为纯 PDF（符合 §sci_save）。
sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  pdf(pdf_file, width = width, height = height)
  print(plot_obj)
  dev.off()
  if (Sys.getenv("RAPA_QC_PNG") == "1") {
    dir.create(file.path(FIG_DIR, "_qc"), showWarnings = FALSE)
    png(file.path(FIG_DIR, "_qc", paste0(filename, ".png")),
      width = width * 150, height = height * 150, res = 150)
    print(plot_obj)
    dev.off()
  }
}

## PART 1 段基因集：段内并集先合并（方案 §4.2） --------------------------------------------

if (!file.exists(F_GS)) stop(paste0("S06: 需要 ", F_GS, "（S02 产物）。"))
dat.final <- read.csv(F_GS, stringsAsFactors = FALSE)

vec.req <- c("segment", "gs_name", "members")
if (!all(vec.req %in% colnames(dat.final))) {
  stop(paste0("S06: ", F_GS, " 缺列 ", paste(setdiff(vec.req, colnames(
            dat.final)),
        collapse = ", "),
      "；实测列: ", paste(colnames(dat.final), collapse = ", ")))
}

## 每段 = 该段全部基因集成员基因的并集（§4.2「段内并集先合并」）
list.union <- lapply(vec.segment, function(seg) {
    d <- dat.final[dat.final$segment == seg, ]
    g <- unique(unlist(strsplit(d$members, ";", fixed = TRUE)))
    sort(g[nzchar(g)])
  })
names(list.union) <- vec.segment

dat.seg <- data.frame(
  segment = vec.segment,
  n_gene_set = vapply(vec.segment, function(s) {
      length(unique(dat.final$gs_name[dat.final$segment == s]))
    }, integer(1L)),
  n_member_nominal = vapply(vec.segment, function(s) length(list.union[[s]]),
    integer(1L)),
  stringsAsFactors = FALSE)
for (s in vec.segment) {
  if (dat.seg$n_member_nominal[dat.seg$segment == s] < MIN_GS_SIZE) {
    stop(paste0("S06 FAILED: 段 ", s, " 的并集名义成员数 ",
        dat.seg$n_member_nominal[dat.seg$segment == s], " < minSize ",
        MIN_GS_SIZE))
  }
  if (dat.seg$n_gene_set[dat.seg$segment == s] < 1L) {
    stop(paste0("S06 FAILED: 段 ", s, " 在 S02 表中没有任何基因集。"))
  }
}

## 段间重叠（§4.2 规则）：本步只供出该段与非自身段的**最大** Jaccard，
## 「×0.5」的折扣在复合评分（S10）落地 —— 本脚本不改写评分规则。
dat.ov <- read.csv(F_OV, stringsAsFactors = FALSE)
dat.seg$max_jaccard <- vapply(vec.segment, function(s) {
    v <- dat.ov$jaccard[dat.ov$row.segment == s & dat.ov$col.segment != s]
    if (length(v) == 0L) NA_real_ else max(v)
  }, numeric(1L))
dat.seg$overlap_discounted <- dat.seg$max_jaccard > JACCARD_CUT

cat("===== S06 段基因集（段内并集先合并，方案 §4.2）=====\n")
print(dat.seg, row.names = FALSE)

## 段并集的成员基因落盘，作为下游（S07 富集检验 / S09 / S10 / S12）的**单一事实源**：
## 「段 = 该段全部基因集成员的并集」这条规则只在本脚本实现一次，
## 否则 S07 等脚本各写一遍 union 逻辑，日后 S02 口径一变就会静默漂移。
dat.union.out <- data.frame(
  segment = vec.segment,
  n_gene = dat.seg$n_member_nominal,
  members = vapply(list.union, paste, character(1L), collapse = ";"),
  stringsAsFactors = FALSE)
write.csv(dat.union.out, file.path(TAB_DIR, "segment_union_genes.csv"),
  row.names = FALSE)

## PART 2 逐系列算段分 -----------------------------------------------------------------------

## 2026-09-29 起读**基因级**矩阵 <acc>_gene.rds（行 = 基因符号、列 = 真实样本 ID）。
## GSVA 以基因符号为键；探针级 _expr.rds 的行名是 p1..pN，对不上任何基因集且不报错
## （属「静默产出全零」型错误，比崩溃更危险）。
dat.spec <- read.csv(F_SPEC, stringsAsFactors = FALSE)
## [user ruling 2026-09-29 GSE6002] 层 2 的系列集 = **有 primary 对比**的系列。
## GSE6002 是双色补充臂（primary = FALSE、n = 4 功效下限，方案 §487「非主聚合」），
## 若纳入：① 下方 dat.spec[acc & primary] 得 0 行 -> stop「primary 应恰为 1 个」；
##   ② 即便绕过也会改层 3 通过率分母，推翻 2821/8224 正典。故在此显式排除，
## 补充臂的方向读数由 S11 单独消费（S05 已落 layer1_gene_dir.csv 的该系列行）。
vec.acc <- unique(dat.spec$acc[dat.spec$primary])

## 臂长与「是否探索性」的单一事实源 = S05 的层 1 表（避免在此重复 MIN_ARM_N 判定逻辑）
dat.l1 <- read.csv(F_L1, stringsAsFactors = FALSE)
dat.arm <- unique(dat.l1[dat.l1$primary, c("acc", "contrast", "n_num", "n_den",
      "exploratory")])

dat.layer2 <- data.frame()
list.score <- list()

for (acc in vec.acc) {
  f.expr <- file.path("03_results", "S03", "table", paste0(acc, "_gene.rds"))
  if (!file.exists(f.expr)) {
    stop(paste0("S06: 缺基因级矩阵 ", f.expr,
        "（请先跑 S03）；不跳过、不用探针级矩阵顶替。"))
  }
  expr <- readRDS(f.expr)$expr
  cat(sprintf("\n[S06] %s  (%d symbols x %d samples)\n", acc, nrow(expr),
      ncol(expr)))

  list.union.acc <- lapply(list.union, function(g) g[g %in% rownames(expr)])
  vec.eff <- vapply(list.union.acc, length, integer(1L))
  if (any(vec.eff < MIN_GS_SIZE)) {
    stop(paste0("S06: ", acc, " 的段并集在该平台的有效成员数不足 minSize=",
        MIN_GS_SIZE,
        "：", paste(sprintf("%s=%d", names(vec.eff), vec.eff),
          collapse = " ")))
  }

  obj.ss <- GSVA::gsva(GSVA::ssgseaParam(expr, list.union.acc,
      minSize = MIN_GS_SIZE,
      alpha = ALPHA_SSGSEA, normalize = TRUE),
    verbose = FALSE)
  mat.score <- as.matrix(obj.ss)
  if (!identical(sort(rownames(mat.score)), sort(vec.segment))) {
    stop(paste0("S06: ", acc, " 的 ssGSEA 返回行名与五段不符：",
        paste(rownames(mat.score), collapse = ", ")))
  }
  if (!identical(colnames(mat.score), colnames(expr))) {
    stop(paste0("S06: ", acc, " 的 ssGSEA 返回列名与表达矩阵不一致。"))
  }
  mat.score <- mat.score[vec.segment, , drop = FALSE]
  list.score[[acc]] <- mat.score
  saveRDS(mat.score, file.path(TAB_DIR, paste0("segment_score_", acc, ".rds")))

  ## 设计：单一事实源 = contrast_spec.csv（formula / term / primary 对比）
  ##       + design_used_<acc>.csv（含 sample_id 的实际设计矩阵，由 S05 落盘）
  dat.p <- dat.spec[dat.spec$acc == acc & dat.spec$primary, , drop = FALSE]
  if (nrow(dat.p) != 1L) {
    stop(paste0("S06: ", acc,
        " 在 contrast_spec.csv 中的 primary 对比应恰为 1 个，实得 ",
        nrow(dat.p), " 个。"))
  }
  
  if (!"factor.vars" %in% colnames(dat.spec)) {
    stop(paste0("S06: contrast_spec.csv 缺 factor.vars",
      "列（分类协变量类型契约），请先重跑 S05。"))
  }
  f.file <- file.path("03_results", "S05", "table", paste0("design_used_", acc,
      ".csv"))
  if (!file.exists(f.file)) {
    stop(paste0("S06: 缺设计表 ", f.file, "（请先跑 S05）。"))
  }
  dat.design <- read.csv(f.file, stringsAsFactors = FALSE)
  dat.design <- dat.design[match(colnames(mat.score), dat.design$sample_id), ,
    drop = FALSE]
  rownames(dat.design) <- colnames(mat.score)

  dat.f <- stats::model.frame(as.formula(dat.p$formula), data = dat.design)
  tt.design <- attr(dat.f, "terms")
  dat.f[[dat.p$term]] <- droplevels(as.factor(dat.f[[dat.p$term]]))
  vec.fv <- strsplit(dat.p$factor.vars, ";", fixed = TRUE)[[1L]]
  vec.fv <- vec.fv[nzchar(vec.fv)]
  for (v in vec.fv) {
    if (!v %in% colnames(dat.f)) {
      stop(paste0("S06: ", acc, " 声明的分类协变量 '", v,
          "' 不在设计表内（实测列: ",
          paste(colnames(dat.f), collapse = ", "), "）。"))
    }
    dat.f[[v]] <- droplevels(as.factor(dat.f[[v]]))
  }
  dat.f <- dat.f[, vapply(dat.f, function(x) !is.numeric(x) || length(unique(
          x)) > 1L, TRUE),
    drop = FALSE]
  mat.design <- stats::model.matrix(tt.design, data = dat.f)

  vec.lv <- levels(dat.f[[dat.p$term]])
  vec.col <- paste0(dat.p$term, vec.lv)
  i.num <- match(dat.p$num, vec.lv)
  i.den <- match(dat.p$den, vec.lv)
  if (anyNA(c(i.num, i.den))) {
    stop(paste0("S06: ", acc, " 的 primary 对比水平不存在: ", dat.p$num, " vs ",
        dat.p$den,
        " ; 实测水平: ", paste(vec.lv, collapse = " | ")))
  }
  lgl.has <- vec.col[c(i.num, i.den)] %in% colnames(mat.design)
  if (!any(lgl.has)) {
    stop(paste0("S06: ", acc, " 的 primary 对比两侧都被吸收为基准，无法估计。"))
  }
  mat.con <- matrix(0, nrow = ncol(mat.design), ncol = 1L,
    dimnames = list(colnames(mat.design), dat.p$contrast))
  mat.con[vec.col[c(i.num, i.den)][lgl.has], 1L] <- c(1, -1)[lgl.has]

  fit <- limma::lmFit(mat.score, design = mat.design)
  fit <- limma::contrasts.fit(fit, mat.con)
  fit <- limma::eBayes(fit)

  ## 直接从 fit 取（5 行 x 1 对比），不经 topTable 做行名匹配 —— 少一次可能错配的查表。
  ## 后验标准误 = stdev.unscaled * sqrt(s2.post)，供 S07 的随机效应元分析（rma 需要 sei）。
  vec.lfc <- fit$coefficients[, 1L]
  vec.se <- fit$stdev.unscaled[, 1L] * sqrt(fit$s2.post)
  vec.t <- fit$t[, 1L]
  vec.p <- fit$p.value[, 1L]
  vec.padj <- p.adjust(vec.p, method = "BH")

  ## 臂长/探索性：取自 S05 实测（同一设计、同一矩阵列），并做一致性断言防设计漂移
  dat.a <- dat.arm[dat.arm$acc == acc & dat.arm$contrast == dat.p$contrast, ,
    drop = FALSE]
  if (nrow(dat.a) != 1L) {
    stop(paste0("S06: ", acc, " / ", dat.p$contrast,
        " 在 layer1_gene_dir.csv 中无唯一臂长记录（实得 ", nrow(dat.a),
        " 行）。"))
  }
  tab.arm <- table(dat.f[[dat.p$term]])
  if (as.integer(tab.arm[dat.p$num]) != dat.a$n_num ||
    as.integer(tab.arm[dat.p$den]) != dat.a$n_den) {
    stop(paste0("S06: ", acc, " 的臂长与 S05 不一致（本步 ", dat.p$num, "=",
        as.integer(tab.arm[dat.p$num]), " / ", dat.p$den, "=",
        as.integer(tab.arm[dat.p$den]), "；S05 ", dat.a$n_num, " / ",
        dat.a$n_den,
        "）—— 设计与层 1 漂移。"))
  }

  for (seg in vec.segment) {
    v.num <- mat.score[seg, dat.f[[dat.p$term]] == vec.lv[i.num]]
    v.den <- mat.score[seg, dat.f[[dat.p$term]] == vec.lv[i.den]]
    d.cohen <- (mean(v.num) - mean(v.den)) /
    sqrt((var(v.num) * (length(v.num) - 1L) + var(v.den) * (length(
            v.den) - 1L)) /
      (length(v.num) + length(v.den) - 2L))
    ## 两臂的均值/标准差一并落盘，供 S07 的随机效应元分析按官方 metafor::escalc 计算
    ## 标准化效应量（ssGSEA 评分经 normalize=TRUE 后量纲随数据集而变，跨集合并须用无量纲效应）。
    direction.seg <- ifelse(vec.padj[seg] >= PVAL_CUTOFF, "indeterminate",
      ifelse(vec.lfc[seg] > 0, "up", "down"))
    dat.layer2 <- rbind(dat.layer2, data.frame(
        acc = acc, segment = seg,
        n_gene_set = dat.seg$n_gene_set[dat.seg$segment == seg],
        n_member_nominal = dat.seg$n_member_nominal[dat.seg$segment == seg],
        n_member_eff = vec.eff[[seg]],
        max_jaccard = dat.seg$max_jaccard[dat.seg$segment == seg],
        overlap_discounted = dat.seg$overlap_discounted[dat.seg$segment == seg],
        contrast = dat.p$contrast, num = dat.p$num, den = dat.p$den,
      n_num = dat.a$n_num, n_den = dat.a$n_den, exploratory = dat.a$exploratory,
        mean_num = mean(v.num), sd_num = stats::sd(v.num),
        mean_den = mean(v.den), sd_den = stats::sd(v.den),
        log2FC = vec.lfc[seg], se = vec.se[seg], t = vec.t[seg],
        P.Value = vec.p[seg], adj.P.Val = vec.padj[seg], cohens_d = d.cohen,
        direction = direction.seg, stringsAsFactors = FALSE))
  }
  cat(sprintf(paste0("segments scored = 5 ; contrast = %s ; n = %d vs %d ;",
    "exploratory = %s\n"),
      dat.p$contrast, dat.a$n_num, dat.a$n_den, dat.a$exploratory))
}

saveRDS(list.score, file.path(TAB_DIR, "segment_score_matrix.rds"))
write.csv(dat.layer2, file.path(TAB_DIR, "layer2_segment_all.csv"),
  row.names = FALSE)
for (acc in vec.acc) {
  write.csv(dat.layer2[dat.layer2$acc == acc, , drop = FALSE],
    file.path(TAB_DIR, paste0("layer2_segment_", acc,
        ".csv")), row.names = FALSE)
}

## PART 2.5 图件 -------------------------------------------------------------------------------
## 命名遵守规范 16（描述性英文，禁 Fig 前缀）；图内文字全英文（规范 17）；一律单面板。

dat.l2.plot <- dat.layer2
dat.l2.plot$segment <- factor(dat.l2.plot$segment, levels = vec.segment)
## 探索性系列（臂 n < 10，§9 Q12）在轴标签上直接标出，避免与预注册强度混同（规则 35 可判读性）
vec.lab <- vec.acc
vec.lab[vec.acc %in% unique(dat.layer2$acc[dat.layer2$exploratory])] <-
paste0(vec.lab[vec.acc %in% unique(dat.layer2$acc[dat.layer2$exploratory])],
  " (exploratory)")
dat.l2.plot$acc <- factor(as.character(dat.l2.plot$acc),
  levels = rev(vec.acc), labels = rev(vec.lab))
dat.l2.plot$sig <- !is.na(dat.l2.plot$direction) &
dat.l2.plot$direction %in% c("up", "down")

## 图 1：段 x 系列 的 log2FC 热图。格内标数值，显著格以 `*` 标注
## （user 2026-09-29 要求去除原黑色描边；`*` 为 SCI 热图惯例，且不遮挡双向色）。
## 数值保留 3 位：ssGSEA 评分经 normalize=TRUE 后尺度很小（实测 |log2FC| <= 0.13），
## 2 位小数会把几乎所有格显示成 0.01/-0.00，读者无法判读（规则 35）。
p.fc <- ggplot(dat.l2.plot, aes(x = segment, y = acc, fill = log2FC)) +
geom_tile(color = "white", linewidth = 0.5) +
geom_text(aes(label = sprintf("%.3f%s", log2FC, ifelse(sig, "*", ""))),
  size = 3.4, colour = "grey15") +
scale_fill_gradient2(low = "#3B6BA5", mid = "white", high = "#C44E52",
  midpoint = 0,
  name = expression(log[2]~FC)) +
labs(title = "Segment fold change by series",
  subtitle = "ssGSEA score of the segment gene-set union; asterisk = BH < 0.05",
  x = "Segment", y = "Series") +
sci_theme
sci_save(p.fc, "Layer2_Segment_FoldChange_By_Series", width = 9, height = 4.5)

## 图 2：段 x 系列 的 Cohen's d 分组柱（§4.2 要求效应量与显著性并列报告）
sci_save(ggplot(dat.l2.plot, aes(x = segment, y = cohens_d, fill = acc)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(title = "Cohen's d per segment and series",
    subtitle = paste0("Positive d = higher score in the numerator arm of the ",
      "declared primary contrast"),
    x = "Segment", y = "Cohen's d") +
  scale_fill_discrete(name = "Series") +
  sci_theme,
  "Layer2_CohensD_By_Series", width = 9, height = 4.8)

## PART 3 收口 -------------------------------------------------------------------------------

n.indet <- sum(dat.layer2$direction == "indeterminate")
message(paste0("\n===== S06 layer2 summary =====\nseries scored = ", length(
      vec.acc),
    "  rows = ", nrow(dat.layer2),
    "\nindeterminate = ", n.indet, " / ", nrow(dat.layer2),
    " (", round(100 * n.indet / max(nrow(dat.layer2),
        1L)), "%)  <- §4.2 须随结果报告",
    paste0("\nexploratory contrasts (arms n < 10, must be excluded",
      "from C_k weights): "),
    paste(unique(dat.layer2$acc[dat.layer2$exploratory]),
      collapse = ", "),
    "\nsegments with max_jaccard > ", JACCARD_CUT,
    " (discounted in S10 composite score): ",
    paste(unique(dat.layer2$segment[dat.layer2$overlap_discounted]),
      collapse = ", ")))
if (nrow(dat.layer2) == 0L) {
  stop("S06 FAILED: no segment score produced; check S02/S03/S05.")
}
message("S06 OK. tables -> 03_results/S06/table/")
