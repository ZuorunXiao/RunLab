#===============================================================
## S11_perturbation.R —— 扰动层与方向仲裁（步骤规格 §8：S11）
##
## 目的：在 L-P 层做雷帕霉素的方向仲裁：
##       (1) GSE124020：人巨噬 M(IFN-g)/M(IL-10) x +/- 雷帕霉素 2x2 析因；
##       (2) GSE323983：HCAEC 四个臂，两个 counts 文件为块结构，须先做块诊断；
##       (3) GSE6002：CD14+ 单核 LPS +/- 雷帕霉素**双色比值**（唯一走双色的系列）；
##       (4) GSE65145：反向警示（静态参考签名，非拟合臂）——Rapa -> IL-6 家族上调；
##       (5) tacrolimus 作为**机制判别对照**（不是"无效应"参照，§5 对照体系重构）。
##
## 定位：分析脚本。产出「臂 x 段」方向矩阵，供 S13 与层 1-3 对齐。
##
## 产出：03_results/S11/table/perturb_arm_matrix.csv     臂 x 段 方向
# 03_results/S11/table/direction_alignment.csv    与 tacrolimus 的机制判别结论 + GSE65145 反向参考
# 03_results/S11/table/gse65145_reverse_signature.csv  GSE65145 静态反向签名（45 基因，作者发表）
# 03_results/S11/Figure/Perturbation_Direction_Forest.pdf      森林图：系列 x 臂 方向指数+95%CI
# 03_results/S11/Figure/Perturbation_Direction_Landscape.pdf   方向景观热图（rapamycin 实测，单列）
##
## 图件清单：
##   可视化: S11 共 2 类重写为定量方向指数图，替代原样本数柱状图/空瓦片）
##     1. Perturbation_Direction_Forest.pdf
##        —— 形式：森林图（y = 系列:臂，x = 有符号方向指数 dir.index=(n.up-n.down)/(n.up+n.down)，
##           点+95%CI 须线，填色=方向；含 GSE65145 静态反向参考点为对照）
##        —— 科学问题：雷帕霉素扰动「往哪走、走多强、多确定」——实测各臂下调主导（抗炎），
##           唯成纤维细胞静态参考上调（反向警示）；indeterminate 臂无显著方向，标于 0。
##        —— 参照文献：meta 分析森林图合成惯例（Borenstein 2009, Introduction to Meta-Analysis,
##           ISBN 978-0-470-05724-7）；方向指数 95% CI 取 Wilson 得分区间（小 n.sig 不致过度自信）。
##     2. Perturbation_Direction_Landscape.pdf
##        —— 形式：方向景观瓦片（y = 系列:臂，单列 rapamycin，fill = dir.index 连续 RdBu，
##           格内写 up/down 计数；indeterminate 臂灰格）
##        —— 科学问题：雷帕霉素各臂方向景观（幅度+计数）。tacrolimus 机制判别缺口不画图、
##           改由正文/Limitations 陈述
##        —— 参照文献：方向景观热图惯例（同森林图）
##
## 图件口径：sci_theme / sci_save 逐字取自规范 §sci_theme / §sci_save 官方定义；
##       分面一律 facet_wrap / facet_grid，不引入任何拼图包；可读性核验见
##       sci_save 之 RAPA_QC_PNG 钩子。
#===============================================================

## PART 0 常量 ------------------------------------------------------------------------------

PVAL_CUTOFF <- 0.05
LOGFC_MIN <- 0.585
vec.two.color.only <- "GSE6002"       ## §9 Q7：双色流程白名单

TAB_DIR <- "03_results/S11/table"
FIG_DIR <- "03_results/S11/Figure"
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

library(limma)                            ## 扰动臂差异检验
library(ggplot2)  ## 扰动图；主题与出口用规范 §sci_theme / §sci_save

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
  ## 可读性核验钩子：默认关闭，交付物仍为纯 PDF；
  ## 置 RAPA_QC_PNG=1 时额外导 _qc/ PNG 供逐张目视核验，不改变 PDF 产出。
  if (Sys.getenv("RAPA_QC_PNG") == "1") {
    dir.create(file.path(FIG_DIR, "_qc"), showWarnings = FALSE)
    png(file.path(FIG_DIR, "_qc", paste0(filename, ".png")),
      width = width * 150, height = height * 150, res = 150)
    print(plot_obj)
    dev.off()
  }
}

## 臂方向汇总：显著基因 up/down 多数表决
summarise_direction <- function(tt, pcut, lfc) {
  ## 返回方向判定 + 显著 up/down 基因数 + 有符号方向指数
  ## （dir.index = (n.up - n.down)/(n.up + n.down)；负=下调主导/抗炎，正=上调主导；
  ##  NA = 无显著方向。供图件表达「方向幅度与置信度」，非仅分类标签）
  sig <- tt$adj.P.Val < pcut
  n.up <- sum(sig & tt$logFC >= lfc, na.rm = TRUE)
  n.down <- sum(sig & tt$logFC <= -lfc, na.rm = TRUE)
  if (n.up == 0L && n.down == 0L) dir <- "indeterminate"
  else if (n.up > n.down) dir <- "up"
  else if (n.down > n.up) dir <- "down"
  else dir <- "mixed"
  n.sig <- n.up + n.down
  dir.index <- if (n.sig > 0L) (n.up - n.down) / n.sig else NA_real_
  list(direction = dir, n.up = n.up, n.down = n.down, dir.index = dir.index)
}

## PART 1 析因模型：GSE124020（2x2 polarization x drug） -------------------------------------------

dat.arm <- data.frame(acc = character(0), arm = character(0), n = integer(0),
  n.up = integer(0), n.down = integer(0), direction = character(0),
  dir.index = numeric(0), stringsAsFactors = FALSE)
dat.align <- data.frame(comparison = character(0), segment = character(0),
  rapa_dir = character(0), tac_dir = character(0),
  verdict = character(0), stringsAsFactors = FALSE)

## PART 2 双色流程白名单校验（Q7） ---------------------------------------------------------------

## 2026-09-29 机械迁移：S03 的消费契约由探针级 _expr.rds 改为基因级 _gene.rds。
## 注意：本段的 Q7 断言存在**语义缺陷**（见 01_data/verify/S03_gene_level_report_2026-09-29.md）——
## 它要求「白名单系列 GSE6002 必须已产出 bulk 矩阵」，但 GSE6002 走的正是双色路径（B2 未实现），
## 永远不会产出 bulk 矩阵，故该断言恒不成立。属既有缺陷，未擅自改写语义，待批准后修。
dat.expr <- list.files(file.path("03_results", "S03", "table"),
  pattern = "_gene.rds$", full.names = FALSE)
vec.expr.acc <- gsub("_gene\\.rds$", "", dat.expr)
vec.in.two.color <- intersect(vec.two.color.only, vec.expr.acc)
if (length(vec.in.two.color) == 0L) {
  stop(paste0("S11 Q7 FAILED: two-colour whitelist series ", vec.two.color.only,
      " has no expression matrix from S03; the whitelist check cannot pass. ",
      "Do not relax it by adding other series."))
}

## PART 3 逐扰动数据集取数并拟合 -------------------------------------------------------------------

for (acc in c("GSE124020", "GSE323983", "GSE6002")) {
  ## 2026-09-29 起改读基因级矩阵（行 = 基因符号、列 = 真实样本 ID）：
  ## 方向签名要与 GSE65145 反向签名按符号比对，探针级 p1..pN 行名对不上任何符号。
  f.expr <- file.path("03_results", "S03", "table", paste0(acc, "_gene.rds"))
  if (!file.exists(f.expr)) {
    stop(paste0("S11 STOP: ", acc, " has no expression matrix from S03; ",
        "S11 will not skip it silently."))
  }
  obj.expr <- readRDS(f.expr)
  cat(sprintf("\n[S11] %s  (%d probes x %d samples)\n", acc, nrow(
        obj.expr$expr),
      ncol(obj.expr$expr)))
  dat.design <- data.frame(sample = colnames(obj.expr$expr),
    stringsAsFactors = FALSE)

  ## 组名按臂推断：优先用 S01 的 characteristics 键，键不可得则不臆造
  if (acc == "GSE124020") {
    # 2x2 析因：极化(IFN-gamma / IL-10) x 药(Rapa / Ctrl)；设计键来自 S01 design_keys.tsv 的 agent
    # 列。
    ## 实测：expr 列名是 GSM 编号（GSM3518827..838），不含 rap/ifn 文本——
    ## 原 grepl(colnames) 推断失效，arm 退化为单水平因子，model.matrix(~ 0 + arm) 报 contrasts 错。
    ## 故改从 design_keys.tsv 按 sample_id = colnames 取 agent，解析极化与药，
    ## 再用 ~ 0 + arm 组均值模型 + makeContrasts 取两个「Rapa vs Ctrl」臂内对比（不预设交互）。
    f.keys <- "03_results/S01/table/design_keys.tsv"
    dk <- read.delim(f.keys, stringsAsFactors = FALSE, check.names = FALSE)
    dk <- dk[dk$acc == "GSE124020" & dk$key == "agent", ]
    agent.vec <- dk$value[match(colnames(obj.expr$expr), dk$sample_id)]
    stopifnot(!any(is.na(agent.vec)))
    polarization <- ifelse(grepl("il-?10", agent.vec, ignore.case = TRUE),
      "IL10", "IFNg")
    drug <- ifelse(grepl("rapamy", agent.vec, ignore.case = TRUE), "Rapa",
      "Ctrl")
    arm <- factor(paste0(drug, "_", polarization))
    dat.design <- data.frame(sample = colnames(obj.expr$expr),
      polarization = polarization, drug = drug, stringsAsFactors = FALSE)
    design <- model.matrix(~ 0 + arm, dat.design)
    colnames(design) <- levels(arm)
    fit <- limma::lmFit(obj.expr$expr, design = design)
    cm <- limma::makeContrasts(
      Rapa_vs_Ctrl_IFNg = Rapa_IFNg - Ctrl_IFNg,
      Rapa_vs_Ctrl_IL10 = Rapa_IL10 - Ctrl_IL10,
      levels = colnames(design))
    fit.c <- limma::contrasts.fit(fit, cm)
    fit.c <- limma::eBayes(fit.c)
    for (cn in c("Rapa_vs_Ctrl_IFNg", "Rapa_vs_Ctrl_IL10")) {
      tt <- limma::topTable(fit.c, coef = cn, number = Inf, adjust = "BH")
      s <- summarise_direction(tt, PVAL_CUTOFF, LOGFC_MIN)
      dat.arm <- rbind(dat.arm, data.frame(acc = acc, arm = cn,
        n = ncol(obj.expr$expr), n.up = s$n.up, n.down = s$n.down,
        direction = s$direction, dir.index = s$dir.index,
        stringsAsFactors = FALSE))
    }
  } else if (acc == "GSE323983") {
    ## 块诊断：donor 为块、arm 为处理（Ctrl/IR/Rapa/IR with Rapa，来自 S01 sample_axis）
    ## 两 counts 文件已并入 expr（13425x23）；按样本名对齐设计，禁静默错配
    f.axis <- "03_results/S01/table/sample_axis.tsv"
    ax <- read.delim(f.axis, stringsAsFactors = FALSE, check.names = FALSE)
    ax <- ax[ax$acc == "GSE323983", ]
    ax <- ax[match(colnames(obj.expr$expr), ax$sample_id),
      c("sample_id", "donor", "arm"), drop = FALSE]
    stopifnot(!any(is.na(ax$sample_id)))
    dat.design <- data.frame(sample = colnames(obj.expr$expr),
      donor = factor(ax$donor), arm = factor(ax$arm),
      stringsAsFactors = FALSE)
    design <- model.matrix(~ donor + arm, dat.design)
    fit <- limma::lmFit(obj.expr$expr, design = design)
    fit <- limma::eBayes(fit)
    for (a in c("Rapa", "IR", "IR with Rapa")) {
      coef.name <- paste0("arm", a)
      if (!coef.name %in% colnames(fit)) next
      tt <- limma::topTable(fit, coef = coef.name, number = Inf,
        adjust = "BH")
      s <- summarise_direction(tt, PVAL_CUTOFF, LOGFC_MIN)
      dat.arm <- rbind(dat.arm, data.frame(acc = acc,
        arm = paste0(a, "_vs_Ctrl"), n = ncol(obj.expr$expr),
        n.up = s$n.up, n.down = s$n.down, direction = s$direction,
        dir.index = s$dir.index, stringsAsFactors = FALSE))
    }
  } else if (acc == "GSE6002") {
    ## 双色：S03 已 loess 归一化并落 M=log2(Cy5/Cy3)=Rapa 在 LPS 背景效应
    ## （expr 1042x4，已为 M 值，勿再 read.maimages）；每阵列即单次对比，
    ## 方向 = M 跨 4 donor 单样本 t 检验（n=4，Holm；S05 同口径）
    M <- obj.expr$expr
    mean.m <- rowMeans(M, na.rm = TRUE)
    p.m <- apply(M, 1L, function(v) {
        v <- v[!is.na(v)]
        if (length(v) < 2L) return(NA_real_)
        t.test(v)$p.value
      })
    p.adj <- p.adjust(p.m, method = "holm")
    tt <- data.frame(logFC = mean.m, adj.P.Val = p.adj,
      stringsAsFactors = FALSE)
    s <- summarise_direction(tt, PVAL_CUTOFF, LOGFC_MIN)
    dat.arm <- rbind(dat.arm, data.frame(acc = acc, arm = "Rapa_vs_LPS",
        n = ncol(M), n.up = s$n.up, n.down = s$n.down,
        direction = s$direction, dir.index = s$dir.index,
        stringsAsFactors = FALSE))
  }
}

write.csv(dat.arm, file.path(TAB_DIR, "perturb_arm_matrix.csv"),
  row.names = FALSE)

# PART 3.5 反向警示静态参考签名（GSE65145，D' 处置）
# ------------------------------------------------
## 来源：作者发表的上调签名（Gillespie 2015, Nucleus, PMID 26652669, Table 3 的 ≥5 倍上调并集），
##       已落盘 01_data/verify/gse65145_rapamycin_up_signature.csv（45 基因 + 命中通路数）。
## 口径：GSE65145 不重定量、不进 S03/S11 拟合流水线（RAW.tar 内为 SeqMonk bedGraph RPKM，非 counts；
# 作者定量本已可用，方案 A 重比对代价不成比例）。仅作方向性诚实披露的静态参考，供 S13 方向性轴引用。
## 权重：Q12 锁死补充面（n=6 < 10，不进主聚合）；非血管细胞（人包皮成纤维 2DD），不可外推至动脉。

file.gse65145.sig <- "01_data/verify/gse65145_rapamycin_up_signature.csv"
if (!file.exists(file.gse65145.sig)) {
  stop("S11 STOP: GSE65145 reverse signature missing at ", file.gse65145.sig,
    "; it must exist as the static reference (D' disposition).")
}
dat.gse65145.sig <- read.csv(file.gse65145.sig, stringsAsFactors = FALSE)
cat(sprintf("\n[S11] GSE65145 静态反向签名：%d 基因（作者发表上调并集）\n",
    nrow(dat.gse65145.sig)))
dat.align <- rbind(dat.align, data.frame(
    comparison = "GSE65145_reverse_reference",
    segment    = "S2_S5_cytokine",
    rapa_dir   = "up",
    tac_dir    = "n/a",
    verdict    = paste0("REVERSE WARNING (static, author-published): rapa up-regulates IL-6 family ",
      "in foreskin fibroblasts; non-vascular, n=6, NOT re-quantified; ",
      "see 01_data/verify/gse65145_rapamycin_up_signature.csv"),
    stringsAsFactors = FALSE))
write.csv(dat.gse65145.sig, file.path(TAB_DIR,
    "gse65145_reverse_signature.csv"),
  row.names = FALSE)

# PART 4 tacrolimus 机制判别对照
# --------------------------------------------------------------------

dat.align <- rbind(dat.align, data.frame(
    comparison = "tacrolimus_vs_rapa", segment = "S1-S5",
    rapa_dir = if (nrow(dat.arm) > 0L) dat.arm$direction[1L] else "n/a",
    tac_dir = "pending",
  verdict = paste0("pending: no tacrolimus expression series in project ",
    "(plan §5 defines it only conceptually); mechanism discrimination ",
    "not data-driven — tacrolimus left as pending cell in the figure"),
    stringsAsFactors = FALSE))
write.csv(dat.align, file.path(TAB_DIR, "direction_alignment.csv"),
  row.names = FALSE)

## PART 4.5 图件 ------------------------------------------------------------------------------
## 分面一律用 facet_wrap，不引入任何拼图包。

## PART 4.5 图件
## 旧图问题（用户裁定无法发表、看不出结果）：F1 以「臂内样本数」为柱高（同系列柱等高，
##   完全不表达扰动结果），仅用颜色标 down/indeterminate 分类；F2 为 2x1 分类瓦片，几乎空白。
## 新图核心：有符号方向指数 dir.index = (n.up - n.down)/(n.up + n.down)
##   （负=下调主导/抗炎，正=上调主导），配 95% CI（显著基因中 down 比例的正态近似）
##   与 up/down 基因计数——直接回答「雷帕霉素往哪走、走多强、多确定」。
## 清理上一版图件（命名已变更，避免 stale PDF 残留误导）
for (f.old in c("Perturbation_Arm_Direction_By_Series",
    "Perturbation_Direction_Alignment")) {
  fp <- file.path(FIG_DIR, paste0(f.old, ".pdf"))
  if (file.exists(fp)) file.remove(fp)
}

## 方向指数与 95% CI（Wilson 得分区间；小 n.sig 时不会过度自信，优于正态近似）
z.ci <- qnorm(0.975)
dat.arm$n.sig <- dat.arm$n.up + dat.arm$n.down
dat.arm$p.down <- dat.arm$n.down / dat.arm$n.sig
center <- (dat.arm$p.down + z.ci^2 / (2 * dat.arm$n.sig)) /
  (1 + z.ci^2 / dat.arm$n.sig)
half <- z.ci * sqrt(dat.arm$p.down * (1 - dat.arm$p.down) / dat.arm$n.sig +
  z.ci^2 / (4 * dat.arm$n.sig^2)) / (1 + z.ci^2 / dat.arm$n.sig)
dat.arm$ci.lo <- 1 - 2 * pmin(1, center + half)
dat.arm$ci.hi <- 1 - 2 * pmax(0, center - half)

## 反向警示静态参考（GSE65145，作者发表上调签名）作为对照点纳入森林图，
## 凸显「实测扰动均下调，唯成纤维细胞静态参考上调（反向）」；CI 取 Wilson 同口径（45 基因全 up）。
n.g <- nrow(dat.gse65145.sig); p.g <- 0
c.g <- (p.g + z.ci^2 / (2 * n.g)) / (1 + z.ci^2 / n.g)
h.g <- z.ci * sqrt(p.g * (1 - p.g) / n.g + z.ci^2 / (4 * n.g^2)) /
  (1 + z.ci^2 / n.g)
ci.lo.g <- 1 - 2 * pmin(1, c.g + h.g)
ci.hi.g <- 1 - 2 * pmax(0, c.g - h.g)
dat.forest <- dat.arm
dat.forest <- rbind(dat.forest, data.frame(acc = "GSE65145",
  arm = "rapa_up_signature (static)", n = n.g, n.up = n.g, n.down = 0L,
  direction = "up", dir.index = 1, n.sig = n.g, p.down = 0,
  ci.lo = ci.lo.g, ci.hi = ci.hi.g, stringsAsFactors = FALSE))
dat.forest$arm.lbl <- paste0(dat.forest$acc, " : ", dat.forest$arm,
  " (n.sig=", dat.forest$n.sig, ")")
lvl <- with(dat.forest, order(acc, -dir.index))
dat.forest$arm.lbl <- factor(dat.forest$arm.lbl,
  levels = rev(dat.forest$arm.lbl[lvl]))

DIR_COL <- c(down = "#2166AC", up = "#B2182B", mixed = "#7F7F7F",
  indeterminate = "#D9D9D9")

sci_save(ggplot(dat.forest, aes(x = dir.index, y = arm.lbl)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  geom_pointrange(data = dat.forest[!is.na(dat.forest$dir.index), ],
    aes(xmin = ci.lo, xmax = ci.hi, colour = direction),
    size = 0.6) +
  geom_text(data = dat.forest[is.na(dat.forest$dir.index), ],
    aes(x = 0, label = "indeterminate"), colour = "#888888", size = 3,
    vjust = -0.8) +
  scale_colour_manual(name = "Direction", values = DIR_COL) +
  coord_cartesian(xlim = c(-1.15, 1.15)) +
  labs(title = "Rapamycin direction index across perturbation arms",
    subtitle = paste0("Negative = down-regulation (anti-inflammatory);",
      "positive = up-regulation"),
    caption = paste0("n.sig = significant genes (padj < 0.05 & |log2FC| >= 0.585), not samples.\n",
      "95% CI: Wilson interval. GSE65145 = static signature. ",
      "n.sig = 1 arms: underpowered, tentative."),
    x = "Direction index = (n.up - n.down) / (n.up + n.down)",
    y = "Dataset : arm") +
  sci_theme,
  "Perturbation_Direction_Forest", width = 10, height = 5)

## 图 2：雷帕霉素方向景观（单系列热图）。tacrolimus 无项目内表达系列，
##   机制判别缺口改由正文/Limitations 陈述，不在图中以「pending」占位格呈现
dat.heat <- data.frame(
  acc = dat.arm$acc,
  arm = dat.arm$arm,
  dir.index = dat.arm$dir.index,
  n.up = dat.arm$n.up,
  n.down = dat.arm$n.down,
  n.sig = dat.arm$n.sig,
  stringsAsFactors = FALSE)
## 臂标签格式必须与森林图 arm.lbl（L365-366）逐字一致，否则 factor 对不上
dat.heat$arm.lbl <- paste0(dat.heat$acc, " : ", dat.heat$arm,
  " (n.sig=", dat.heat$n.sig, ")")
dat.heat$arm.lbl <- factor(dat.heat$arm.lbl,
  levels = levels(dat.forest$arm.lbl))
stopifnot(!any(is.na(dat.heat$arm.lbl)))
dat.heat$label <- ifelse(is.na(dat.heat$dir.index),
  sprintf("indeterminate (n.sig=%d)", dat.heat$n.sig),
  sprintf("%d up / %d down (n.sig=%d)", dat.heat$n.up, dat.heat$n.down,
    dat.heat$n.sig))
sci_save(ggplot(dat.heat, aes(x = factor("rapamycin"), y = arm.lbl)) +
  geom_tile(aes(fill = dir.index), colour = "white") +
  geom_text(aes(label = label), size = 2.8, colour = "black") +
  scale_fill_gradient2(name = "Direction index", low = "#2166AC",
    mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limits = c(-1, 1),
    na.value = "#D9D9D9") +
  labs(title = "Rapamycin direction landscape per perturbation arm",
    subtitle = paste0("Mechanism discrimination vs tacrolimus not feasible:\n",
      "no in-project tacrolimus expression series (plan S5)"),
    x = "Agent", y = "Dataset : arm") +
  sci_theme,
  "Perturbation_Direction_Landscape", width = 9.5, height = 5)

## PART 5 收口 ------------------------------------------------------------------------------

message(paste0("\n===== S11 perturbation summary =====\narm rows = ", nrow(
      dat.arm),
    "  alignment rows = ", nrow(dat.align)))
message(paste0("S11 OK: GSE124020/GSE323983/GSE6002",
  "三系列方向仲裁完成（tacrolimus 项目内无数据，留",
  "pending）。 "),
  "tables -> 03_results/S11/table/")
