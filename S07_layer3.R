#===============================================================
## S07_layer3.R —— 层 3：跨数据集方向一致性（步骤规格 §8：S07）
##
## 目的：① 对五段分别做三路并列的富集检验（cameraPR + roast + fgsea），
##         检验单元 = 该段的**成员并集**（与层 2 同一对象，方案 §4.2「段内并集先合并」）；
##         该三路检验经实测定为**覆盖度依赖性旁证**，不作层 3 主结论；
##       ② 对每个段做随机效应元分析（metafor），按 §4.2 判定
##         `confirmed` / `conflicted` / `indeterminate`，并报告 I² 与 indeterminate 占比。
##         **层 3 主读出 = 本条 SMD 元分析 + 段内逐基因集敏感性轨道的 FDR 通过率。**
##
## 定位：分析脚本。层 3 = 段效应的跨数据集一致性；层 2（样本级段分）在 S06。
## 产出：03_results/S07/table/layer3_evidence_union.csv      主表：段并集的三路检验
##       03_results/S07/table/layer3_sensitivity_sets.csv    敏感性：段内逐基因集的检验留痕
##       03_results/S07/table/layer3_sensitivity_summary.csv 敏感性：段级通过 FDR 的比例
##       03_results/S07/table/layer3_meta_dir.csv            元分析主表（SMD）+ 判定
##       03_results/S07/table/layer3_meta_sensitivity.csv    元分析敏感性表（log2FC + limma SE）
##       03_results/S07/table/layer3_verdict_counts.csv      判定计数
##       03_results/S07/table/checkpoints/S07_<acc>_layer3.rds  PART 级断点（规则 33）
##       03_results/S07/Figure/Layer3_Evidence_Pvalue_Grid.pdf      分面网格：方法 x 数据集
##       03_results/S07/Figure/Layer3_Metaanalysis_Forest.pdf       森林图（合并 SMD + 95% CI）
##       03_results/S07/Figure/Layer3_Heterogeneity_By_Verdict.pdf  I² 分面柱（按判定着色）
##       03_results/S07/Figure/Layer3_Sensitivity_PassRate.pdf      敏感性：段内通过 FDR 的集合比例
##
## 图件清单（规则 35 可视化强制 · SCI 发表导向；动笔前已按方案 §4.2 与已发表惯例定风格/数量）：
##   可视化: S07 共 4 幅
##     1. Layer3_Evidence_Pvalue_Grid.pdf
##        —— 形式：分面网格散点（x = 段，y = -log10 p，facet 方法 x 数据集），虚线标 p = 0.05
##        —— 科学问题：三路独立富集检验（cameraPR / roast / fgsea）在各段的证据强度与一致性。
##        —— 参照：多方法富集 p 值并列网格为证据稳健性标准画法
##           （limma cameraPR, Wu & Smyth 2012, doi:10.1093/nar/gks461, PMID 22684630）。
##     2. Layer3_Metaanalysis_Forest.pdf
##        —— 形式：单面板森林图（x = 合并标准化效应 SMD 与 95% CI，y = 段）
##        —— 科学问题：跨数据集随机效应合并后的段效应方向与区间。
##        —— 参照：森林图为元分析标准呈现（metafor; Viechtbauer 2010, J Stat Softw,
##           doi:10.18637/jss.v036.i03）。
##     3. Layer3_Heterogeneity_By_Verdict.pdf
##        —— 形式：单面板柱（x = 段，y = I²，fill = §4.2 判定）
##        —— 科学问题：各段跨研究异质性 I² 与最终判定的对应。
##        —— 参照：I² 异质性报告为元分析标准（Higgins 2003, BMJ,
##           doi:10.1136/bmj.327.7414.557, PMID 12958120）。
##     4. Layer3_Sensitivity_PassRate.pdf
##        —— 形式：分面柱（x = 段，y = 通过 FDR 的基因集比例，fill = 检验方法，按数据集分面）
##        —— 科学问题：把段当作**单个大集合**时，若该集合在平台上覆盖过广（实测 S3 在
##           GSE323983 覆盖 89% 已测基因），竞争性检验的对比度会被稀释；本图给出
##           「段内逐基因集」轨道下通过 FDR 的比例，用于判读主表结论是否被覆盖度掩盖。
##        —— 参照：逐基因集富集 + 通过率汇总为通路层面的稳健性惯例（同 1 的文献）。
##        —— 为何新增（user ruling 2026-09-29「并集为主 + 段内逐集合作敏感性轨道」）。
##
## 图件口径：sci_theme / sci_save 逐字取自规范 §sci_theme / §sci_save 官方定义；
# 分面用 facet_wrap / facet_grid，不引入任何拼图包；可读性核验见 sci_save 之 RAPA_QC_PNG 钩子。
##
## 前置：S02（基因集）+ S03（基因级矩阵）+ S05（contrast_spec / design_used）+
##       S06（segment_union_genes.csv = 段定义单一事实源；layer2_segment_all.csv = 元分析输入）。
#===============================================================

## PART 0 常量 ------------------------------------------------------------------------------

SEED <- 20260925
set.seed(SEED)
PVAL_CUTOFF <- 0.05
MIN_GS_SIZE <- 20L
MAX_GS_SIZE <- 5000L  ## 逐基因集轨道的上界；段并集轨道不设上界，理由见 PART 2
MIN_N_STUDY <- 2L             ## §4.2：随机效应元分析至少需要 2 个数据集
I2_CUT <- 50  ## I² 异质性判读切点（仅用于图表标注，判定按 §4.2 三段式）
## roast 旋转数（limma 官方默认）；实测提高旋转数只移动下限、
NROT_ROAST <- 1999L
FGSEA_EPS <- 1e-50
## roast 双侧/Mixed 行的最小可能 p（实测 = 5.0e-04）
ROAST_FLOOR <- 2 / (2 * NROT_ROAST + 1)

vec.segment <- c("S1_initiation", "S2_amplification", "S3_effect",
  "S4_resolution", "S5_chronic")

## fgsea 并行后端：显式串行。不写 nproc：官方 setUpBPPARAM 的
## `if (nproc == 1)` 之后是独立 if，Windows 分支会覆盖它，实测 nproc = 1 仍返回 SnowParam。
BP_SERIAL <- BiocParallel::SerialParam()
## 并集轨道三路检验的证据角色：作旁证登记，不作层 3 主结论
ROLE_UNION <- "corroborative_coverage_dependent"

F_SPEC <- file.path("03_results", "S05", "table", "contrast_spec.csv")
F_UNION <- file.path("03_results", "S06", "table", "segment_union_genes.csv")
F_L2 <- file.path("03_results", "S06", "table", "layer2_segment_all.csv")
F_GS <- file.path("03_results", "S02", "table",
  "inflammation_geneset_final.csv")

TAB_DIR <- "03_results/S07/table"
FIG_DIR <- "03_results/S07/Figure"
CKPT_DIR <- file.path(TAB_DIR, "checkpoints")
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)

##   任何一项变化即判缓存未命中，杜绝「上游已重算、断点照旧恢复」的静默错配。
Fingerprint <- function(acc) {
  vec.f <- c(file.path("03_results", "S03", "table", paste0(acc, "_gene.rds")),
    file.path("03_results", "S05", "table", paste0("design_used_",
      acc, ".csv")),
    file.path("03_results", "S05", "table", "contrast_spec.csv"),
    F_GS,
    file.path("03_results", "S06", "table", "segment_union_genes.csv"))
  out <- unname(tools::md5sum(vec.f))
  c(out, SEED = as.character(SEED))
}

library(limma)  ## 层 3 cameraPR 主证据 + roast 自包含检验
library(fgsea)                            ## 层 3 GSEA 置换检验交叉验证
library(BiocParallel)  ## fgsea 并行后端选择
library(metafor)  ## 层 3 随机效应合并（escalc SMD + rma）与 I²
library(ggplot2)  ## 层 3 证据图；主题与出口用规范 §sci_theme / §sci_save

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

## PART 1 输入与设计构造（单一事实源：S05 的 contrast_spec + design_used） -------------------

for (f in c(F_SPEC, F_UNION, F_L2, F_GS)) {
  if (!file.exists(f)) {
    stop(paste0("S07: 缺前置产物 ", f, "（S02 / S05 / S06 须先跑通）。"))
  }
}
dat.spec <- read.csv(F_SPEC, stringsAsFactors = FALSE)
dat.union <- read.csv(F_UNION, stringsAsFactors = FALSE)
dat.gs <- read.csv(F_GS, stringsAsFactors = FALSE)

## 段并集（单一事实源 = S06 落盘）：命名 list<段 = 成员基因向量>
list.union <- lapply(vec.segment, function(seg) {
  strsplit(dat.union$members[dat.union$segment == seg], ";", fixed = TRUE)[[1L]]
  })
names(list.union) <- vec.segment
## 逐基因集的成员（敏感性轨道用）；gs_name 在 S02 表内可能重复出现，取首个非空串
map.mem <- tapply(dat.gs$members, dat.gs$gs_name, function(x) x[1L])
list.setmem <- lapply(vec.segment, function(seg) {
    v <- unique(dat.gs$gs_name[dat.gs$segment == seg])
  out <- lapply(v, function(g) strsplit(map.mem[[g]], ";", fixed = TRUE)[[1L]])
    names(out) <- v
    out
  })
names(list.setmem) <- vec.segment

## PART 2 逐数据集：三路并列检验（主表 = 段并集；敏感性 = 段内逐基因集） ----------------------

## [user ruling 2026-09-29 GSE6002] 层 3 的系列集 = **有 primary 对比**的系列
## （与 S06 同口径）。GSE6002 为双色补充臂（primary = FALSE、n = 4），纳入会
## ① 使 dat.spec[acc & primary] 得 0 行而 stop；② 改变逐集 FDR 通过率分母，
## 推翻 2821/8224 正典。其方向读数由 S11 单独消费，不进主聚合。
vec.acc <- unique(dat.spec$acc[dat.spec$primary])
dat.ev <- data.frame()
dat.sens <- data.frame()

for (acc in vec.acc) {
  f.ck <- file.path(CKPT_DIR, paste0("S07_", acc, "_layer3.rds"))
  if (file.exists(f.ck)) {
    ## 规则 33：读取失败一律报错，不以 tryCatch 吞错充当续传（规则 26）
    obj.ck <- readRDS(f.ck)
    if (!identical(obj.ck$fp, Fingerprint(acc))) {
      message(paste0("[S07] ", acc,
        " checkpoint fingerprint mismatch -> recompute"))
    } else {
      dat.ev <- rbind(dat.ev, obj.ck$ev)
      dat.sens <- rbind(dat.sens, obj.ck$sens)
      cat(sprintf("[S07] %-11s checkpoint restored\n", acc))
      next
    }
  }

  expr <- readRDS(file.path("03_results", "S03", "table", paste0(acc,
        "_gene.rds")))$expr
  dat.p <- dat.spec[dat.spec$acc == acc & dat.spec$primary, , drop = FALSE]
  if (nrow(dat.p) != 1L) {
    stop(paste0("S07: ", acc,
        " 在 contrast_spec.csv 中的 primary 对比应恰为 1 个，实得 ",
        nrow(dat.p), " 个。"))
  }
  if (!"factor.vars" %in% colnames(dat.spec)) {
    stop(paste0("S07: contrast_spec.csv 缺 factor.vars",
      "列（分类协变量类型契约），请先重跑 S05。"))
  }
  dat.design <- read.csv(file.path("03_results", "S05", "table",
      paste0("design_used_", acc, ".csv")),
    stringsAsFactors = FALSE)
  dat.design <- dat.design[match(colnames(expr), dat.design$sample_id), ,
    drop = FALSE]
  rownames(dat.design) <- colnames(expr)
  dat.f <- stats::model.frame(as.formula(dat.p$formula), data = dat.design)
  tt.design <- attr(dat.f, "terms")
  dat.f[[dat.p$term]] <- droplevels(as.factor(dat.f[[dat.p$term]]))

  vec.fv <- strsplit(dat.p$factor.vars, ";", fixed = TRUE)[[1L]]
  vec.fv <- vec.fv[nzchar(vec.fv)]
  for (v in vec.fv) {
    if (!v %in% colnames(dat.f)) {
      stop(paste0("S07: ", acc, " 声明的分类协变量 '", v,
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
    stop(paste0("S07: ", acc, " 的 primary 对比水平不存在: ", dat.p$num, " vs ",
        dat.p$den))
  }
  lgl.has <- vec.col[c(i.num, i.den)] %in% colnames(mat.design)
  if (!any(lgl.has)) {
    stop(paste0("S07: ", acc, " 的 primary 对比两侧都被吸收为基准，无法估计。"))
  }
  mat.con <- matrix(0, nrow = ncol(mat.design), ncol = 1L,
    dimnames = list(colnames(mat.design), dat.p$contrast))
  mat.con[vec.col[c(i.num, i.den)][lgl.has], 1L] <- c(1, -1)[lgl.has]

  fit <- limma::lmFit(expr, design = mat.design)
  fit <- limma::contrasts.fit(fit, mat.con)
  fit <- limma::eBayes(fit)

  ## 全基因 t 统计量（三路检验共用）。cameraPR / roast 的 index 是「对 statistic 位置」的选择，
  ## 排序会破坏位置对应，故这二者用**未排序版**；fgsea 另要求预排序基因列表，单用排序版。
  vec.t <- as.numeric(fit$t[, 1L])
  names(vec.t) <- rownames(expr)
  vec.t <- vec.t[!is.na(vec.t)]
  if (length(vec.t) < 3L) {
    stop(paste0("S07: ", acc, " 可用基因数 ", length(vec.t),
        " < 3，cameraPR 无法运行。"))
  }
  vec.t.sorted <- sort(vec.t, decreasing = TRUE)

  vec.ev <- list()
  vec.sn <- list()

  for (seg in vec.segment) {
    ## ---- 轨道 A（主表）：段并集作为单个基因集 ----
    ## 覆盖度必须一并登记：实测 S3 并集在 GSE323983 覆盖该平台 89% 的已测基因，
    ## 竞争性检验（cameraPR / fgsea）的对比基准只剩余下部分，判读时须知情。
    vec.mem <- intersect(list.union[[seg]], names(vec.t))
    if (length(vec.mem) < MIN_GS_SIZE) {
      stop(paste0("S07: ", acc, "/", seg, " 的并集有效成员数 ", length(vec.mem),
          " < minSize ", MIN_GS_SIZE, "。"))
    }
    cov.seg <- length(vec.mem) / length(vec.t)

    ## (1) cameraPR —— 主检验。官方 \value 返回 NGenes / Direction / PValue / FDR 等列，
    ##     必须**按列名**取 PValue；按位置取会拿到 NGenes（基因集大小）当作 p 值。
    res.pr <- limma::cameraPR(statistic = vec.t, index = list(tmp = vec.mem),
      sort = FALSE)
    vec.miss <- setdiff(c("PValue", "NGenes", "Direction"), names(res.pr))
    if (length(vec.miss) > 0L) {
      stop(paste0("S07: cameraPR 返回列缺 ", paste(vec.miss, collapse = ", "),
          "；实测列: ", paste(names(res.pr), collapse = ", ")))
    }
    vec.ev[[length(vec.ev) + 1L]] <- data.frame(
      acc = acc, segment = seg, tool = "cameraPR", track = "union",
      n_gene_eff = length(vec.mem), coverage = cov.seg,
      statistic = as.character(res.pr$Direction[1L]),
      pval = as.numeric(res.pr$PValue[1L]),
      pval_censored = FALSE, na_reason = "", evidence_role = ROLE_UNION,
      stringsAsFactors = FALSE)

    ## (2) roast —— 官方返回 class "Roast"；$p.value 是数据框，列名为 Active.Prop / P.Value，
    ##     行名为 Down / Up / UpOrDown / Mixed（实测无 direction 字段）。
    res.roast <- limma::roast(expr, index = vec.mem, design = mat.design,
      coef = 1L,
      nrot = NROT_ROAST)
    vec.rn <- rownames(res.roast$p.value)
    vec.want <- c("Down", "Up", "UpOrDown", "Mixed")
    vec.miss <- setdiff(vec.want, vec.rn)
    if (length(vec.miss) > 0L) {
      stop(paste0("S07: roast 返回行缺 ", paste(vec.miss, collapse = ", "),
          "；实测行名: ", paste(vec.rn, collapse = ", ")))
    }
    vec.p.roast <- vapply(vec.want, function(k) {
        as.numeric(res.roast$p.value[k, "P.Value"])
      }, numeric(1L))
    vec.ev[[length(vec.ev) + 1L]] <- data.frame(
      acc = acc, segment = seg, tool = "roast", track = "union",
      n_gene_eff = length(vec.mem), coverage = cov.seg,
      statistic = paste0(paste(vec.want, sprintf("%.3g", vec.p.roast),
          sep = "=",
          collapse = ", "),
        ";nrot=", NROT_ROAST),
      pval = vec.p.roast[["Mixed"]],
      pval_censored = vec.p.roast[["Mixed"]] <= ROAST_FLOOR,
      na_reason = "", evidence_role = ROLE_UNION, stringsAsFactors = FALSE)

    ## (3) fgsea —— 预排序基因列表。maxSize 不设上界（= 排序列表长度）：
    ##     段并集的大小是方案 §4.2 的设计后果，不能被 size 上限静默剔除；
    ##     覆盖度已单独登记，供判读时折价。
    set.seed(SEED)
    res.fg <- fgsea::fgsea(pathways = list(tmp = vec.mem), stats = vec.t.sorted,
      minSize = MIN_GS_SIZE, maxSize = length(vec.t.sorted),
      BPPARAM = BP_SERIAL)
    if (nrow(res.fg) == 0L) {
      stop(paste0("S07: fgsea 未返回结果（", acc, "/", seg,
          "）—— 不得静默跳过。"))
    }
    vec.ev[[length(vec.ev) + 1L]] <- data.frame(
      acc = acc, segment = seg, tool = "fgsea", track = "union",
      n_gene_eff = length(vec.mem), coverage = cov.seg,
      statistic = sprintf("ES=%.4f,NES=%.4f", res.fg$ES[1L], res.fg$NES[1L]),
      pval = as.numeric(res.fg$pval[1L]),
      pval_censored = isTRUE(as.numeric(res.fg$pval[1L]) <= FGSEA_EPS),
      na_reason = if (is.na(res.fg$pval[1L])) {
        paste0("fgsea p not estimable (modeFraction < 10; eps floor = ",
          FGSEA_EPS,
          ")")
      } else {
        ""
      },
      evidence_role = ROLE_UNION, stringsAsFactors = FALSE)

    ## ---- 轨道 B（敏感性，user ruling 2026-09-29）：段内逐基因集 ----
    ## 只纳入在该平台有效成员数 >= minSize 的集合（其余会被官方 minSize 过滤掉）
    list.mem.ok <- lapply(list.setmem[[seg]], function(v) intersect(v,
        names(vec.t)))
    list.mem.ok <- list.mem.ok[vapply(list.mem.ok, length,
        integer(1L)) >= MIN_GS_SIZE]
    if (length(list.mem.ok) == 0L) next

    ## cameraPR 支持一次传 list，一趟出每集合一行（比逐集合调用快得多）
    res.pr.s <- limma::cameraPR(statistic = vec.t, index = list.mem.ok,
      sort = FALSE)
    res.pr.s$gs_name <- rownames(res.pr.s)
    vec.sn[[length(vec.sn) + 1L]] <- data.frame(
      acc = acc, segment = seg, tool = "cameraPR", gs_name = res.pr.s$gs_name,
      n_gene_eff = res.pr.s$NGenes, pval = res.pr.s$PValue,
      padj = if ("FDR" %in% names(res.pr.s)) res.pr.s$FDR else NA_real_,
      direction = res.pr.s$Direction, pval_censored = FALSE,
      stringsAsFactors = FALSE)


    set.seed(SEED)
    res.fg.s <- fgsea::fgsea(pathways = list.mem.ok, stats = vec.t.sorted,
      minSize = MIN_GS_SIZE, maxSize = MAX_GS_SIZE,
      BPPARAM = BP_SERIAL)
    if (nrow(res.fg.s) > 0L) {
      vec.sn[[length(vec.sn) + 1L]] <- data.frame(
        acc = acc, segment = seg, tool = "fgsea", gs_name = res.fg.s$pathway,
        n_gene_eff = res.fg.s$size, pval = res.fg.s$pval, padj = res.fg.s$padj,
        direction = ifelse(res.fg.s$NES > 0, "Up", "Down"),
        pval_censored = !is.na(res.fg.s$pval) & res.fg.s$pval <= FGSEA_EPS,
        stringsAsFactors = FALSE)
    }
  }

  dat.ev.acc <- do.call(rbind, vec.ev)
  dat.sn.acc <- if (length(vec.sn) == 0L) data.frame() else do.call(rbind,
    vec.sn)
  saveRDS(list(ev = dat.ev.acc, sens = dat.sn.acc,
    fp = Fingerprint(acc)), f.ck)
  dat.ev <- rbind(dat.ev, dat.ev.acc)
  dat.sens <- rbind(dat.sens, dat.sn.acc)
  cat(sprintf("[S07] %-11s done (segments=%d, per-set rows=%d)\n",
      acc, length(vec.segment), nrow(dat.sn.acc)))
}

write.csv(dat.ev, file.path(TAB_DIR, "layer3_evidence_union.csv"),
  row.names = FALSE)
write.csv(dat.sens, file.path(TAB_DIR, "layer3_sensitivity_sets.csv"),
  row.names = FALSE)

## 敏感性轨道按段汇总：段级证据 = 「通过 FDR 的集合数 / 该段实际检验的集合数」
dat.sens.sum <- do.call(rbind, lapply(split(dat.sens, list(dat.sens$acc,
        dat.sens$segment,
        dat.sens$tool), drop = TRUE),
    function(d) {
      data.frame(acc = d$acc[1L], segment = d$segment[1L], tool = d$tool[1L],
        n_set = nrow(d), n_pass_fdr = sum(d$padj < PVAL_CUTOFF,
          na.rm = TRUE),
        n_up = sum(d$padj < PVAL_CUTOFF & d$direction == "Up",
          na.rm = TRUE),
        n_down = sum(d$padj < PVAL_CUTOFF & d$direction == "Down",
          na.rm = TRUE),
        n_censored = sum(d$pval_censored, na.rm = TRUE),
        stringsAsFactors = FALSE)
    }))
dat.sens.sum$pass_rate <- dat.sens.sum$n_pass_fdr / dat.sens.sum$n_set
write.csv(dat.sens.sum, file.path(TAB_DIR, "layer3_sensitivity_summary.csv"),
  row.names = FALSE)

## PART 3 随机效应元分析（metafor） ------------------------------------------------------------

dat.l2 <- read.csv(F_L2, stringsAsFactors = FALSE)
vec.req <- c("acc", "segment", "mean_num", "sd_num", "n_num", "mean_den",
  "sd_den", "n_den",
  "log2FC", "se", "adj.P.Val", "exploratory")
vec.miss <- setdiff(vec.req, colnames(dat.l2))
if (length(vec.miss) > 0L) {
  stop(paste0("S07: layer2_segment_all.csv 缺列 ", paste(vec.miss,
        collapse = ", "),
      "；S06 契约变更。"))
}

## §4.2 判定（三段式，显著性为先决条件）：
##   不显著（合并 p >= PVAL_CUTOFF）           -> indeterminate
##   显著 且 各数据集效应**同号**（与合并效应同号）-> confirmed
##   显著 且 方向不一致                        -> conflicted（并列报告，禁止择优）
MetaSeg <- function(dat.study, yi, vi) {
  mod <- metafor::rma(yi = yi, vi = vi, method = "REML", test = "z")
  b.pooled <- as.numeric(mod$b)
  p.pooled <- as.numeric(mod$pval)
  lgl.agree <- all(sign(yi) == sign(b.pooled))
  data.frame(n_study = length(yi), pooled = b.pooled,
    ci_lo = as.numeric(mod$ci.lb), ci_hi = as.numeric(mod$ci.ub),
    pval = p.pooled, I2 = as.numeric(mod$I2), tau2 = as.numeric(mod$tau2),
    sign_agree = lgl.agree,
    verdict = if (p.pooled >= PVAL_CUTOFF) "indeterminate"
    else if (lgl.agree) "confirmed" else "conflicted",
    stringsAsFactors = FALSE)
}

dat.meta <- data.frame()
dat.meta.sens <- data.frame()
vec.skip <- character(0)

for (seg in vec.segment) {
  d <- dat.l2[dat.l2$segment == seg, , drop = FALSE]
  if (nrow(d) < MIN_N_STUDY) {
    vec.skip <- c(vec.skip, seg)
    next
  }
  ## 主表：escalc(measure="SMD") -> Hedges' g + 抽样方差
  es <- metafor::escalc(measure = "SMD",
    m1i = d$mean_num, sd1i = d$sd_num, n1i = d$n_num,
    m2i = d$mean_den, sd2i = d$sd_den, n2i = d$n_den)
  m1 <- MetaSeg(d, es$yi, es$vi)
  m1$segment <- seg
  m1$n_sig_layer2 <- sum(d$adj.P.Val < PVAL_CUTOFF)
  m1$n_exploratory <- sum(d$exploratory)
  m1$contributing <- paste(d$acc, collapse = ";")
  dat.meta <- rbind(dat.meta, m1)

  ## 敏感性表：limma 的 log2FC 与后验 SE（保留配对/区组设计信息，但量纲随数据集而变）
  m2 <- MetaSeg(d, d$log2FC, d$se^2)
  m2$segment <- seg
  m2$contributing <- paste(d$acc, collapse = ";")
  dat.meta.sens <- rbind(dat.meta.sens, m2)
}

write.csv(dat.meta, file.path(TAB_DIR, "layer3_meta_dir.csv"),
  row.names = FALSE)
write.csv(dat.meta.sens, file.path(TAB_DIR, "layer3_meta_sensitivity.csv"),
  row.names = FALSE)

## PART 3.5 图件 -------------------------------------------------------------------------------
## 命名遵守规范 16（描述性英文，禁 Fig 前缀）；图内文字全英文（规范 17）；
## 分面一律用 facet_wrap / facet_grid，不引入任何拼图包。

n.ev.na <- sum(is.na(dat.ev$pval))
n.ev.floor <- sum(dat.ev$pval_censored, na.rm = TRUE)
dat.ev.plot <- dat.ev[!is.na(dat.ev$pval) & dat.ev$pval > 0, ]
n.ev.omit <- nrow(dat.ev) - nrow(dat.ev.plot)
dat.ev.plot$neglog10p <- -log10(dat.ev.plot$pval)
dat.ev.plot$segment <- factor(dat.ev.plot$segment, levels = vec.segment)
dat.ev.plot$tool <- factor(dat.ev.plot$tool, levels = c("cameraPR", "roast",
    "fgsea"))
sci_save(ggplot(dat.ev.plot, aes(x = segment, y = neglog10p)) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_hline(yintercept = -log10(PVAL_CUTOFF), linetype = "dashed",
    colour = "grey40") +
  facet_grid(tool ~ acc, scales = "free_y") +
  labs(title = "Three independent enrichment tests on each segment union",
    subtitle = paste0(
      "Union track = corroborative evidence (coverage-dependent); ",
      "y axis is free per test row"),
    caption = paste0("Dashed line = p = 0.05; n = ", n.ev.omit,
      " p omitted (NA or zero, of which NA = ", n.ev.na, ")",
      "\nPoints at ", -log10(FGSEA_EPS),
      " or the roast floor are censored, not estimates"),
    x = "Segment", y = expression(-log[10]~p)) +
  sci_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    plot.caption = element_text(hjust = 0)),
  "Layer3_Evidence_Pvalue_Grid", width = 12, height = 7)

## 图 2：森林图（合并标准化效应 SMD + 95% CI）
dat.meta.plot <- dat.meta
dat.meta.plot$segment <- factor(dat.meta.plot$segment,
  levels = rev(vec.segment))
sci_save(ggplot(dat.meta.plot, aes(x = pooled, y = segment)) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
    orientation = "y", width = 0.22) +
  geom_point(aes(colour = verdict), size = 3.2) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(title = "Random-effects meta-analysis of segment effects",
    subtitle = paste0("Primary layer-3 readout; pooled Hedges' g across ",
      "datasets, positive = higher in numerator arm"),
    x = "Pooled SMD (95% CI)", y = "Segment") +
  scale_colour_discrete(name = "Verdict") +
  sci_theme,
  "Layer3_Metaanalysis_Forest", width = 8.5, height = 5)

## 图 3：异质性 I² 与判定
dat.meta.plot$segment <- factor(dat.meta.plot$segment, levels = vec.segment)
dat.meta.plot$verdict <- factor(dat.meta.plot$verdict,
  levels = c("confirmed", "conflicted",
    "indeterminate"))
sci_save(ggplot(dat.meta.plot, aes(x = segment, y = I2, fill = verdict)) +
  geom_col() +
  geom_hline(yintercept = I2_CUT, linetype = "dashed",
    colour = "grey40") +
  scale_y_continuous(limits = c(0, 100),
    expand = expansion(mult = c(0, 0.04))) +
  labs(title = "Heterogeneity (I-squared) by meta-analysis verdict",
    subtitle = paste0("Dashed line = I2 = ", I2_CUT,
      "%; verdict follows the pre-registered rule"),
    x = "Segment", y = expression(I^2~"(%)")) +
  scale_fill_discrete(name = "Verdict") +
  sci_theme,
  "Layer3_Heterogeneity_By_Verdict", width = 8.5, height = 4.8)

## 图 4：敏感性轨道 —— 段内逐基因集通过 FDR 的比例
dat.sens.sum$segment <- factor(dat.sens.sum$segment, levels = vec.segment)
dat.sens.sum$tool <- factor(dat.sens.sum$tool, levels = c("cameraPR", "fgsea"))
dat.sens.sum$rate <- dat.sens.sum$pass_rate
sci_save(ggplot(dat.sens.sum, aes(x = segment, y = rate, fill = tool)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_text(aes(label = paste0(n_pass_fdr, "/", n_set)),
    position = position_dodge(width = 0.8), vjust = -0.4,
    size = 2.6, colour = "grey20") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  facet_wrap(~ acc) +
  labs(
    title = "Sensitivity: share of gene sets passing FDR within each segment",
    subtitle = paste0("Labels = sets with FDR < 0.05 / sets tested; ",
      "roast is not run per set (cost)"),
    caption = paste0("Sets whose fgsea p reached the eps floor are ",
      "counted as passing; see n_censored in the summary table"),
    x = "Segment", y = "Share passing FDR") +
  scale_fill_discrete(name = "Test") +
  sci_theme +

  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    plot.caption = element_text(hjust = 0)),
  "Layer3_Sensitivity_PassRate", width = 12, height = 6)

## PART 4 收口 -------------------------------------------------------------------------------

n.verdict <- table(factor(dat.meta$verdict,
    levels = c("confirmed", "conflicted",
      "indeterminate")))
write.csv(data.frame(verdict = names(n.verdict), n = as.integer(n.verdict),
    stringsAsFactors = FALSE),
  file.path(TAB_DIR, "layer3_verdict_counts.csv"), row.names = FALSE)

message(paste0("\n===== S07 layer3 summary =====\n",
    "primary readout: SMD meta-analysis (layer3_meta_dir.csv) + ",
    "per-set FDR pass rate (layer3_sensitivity_summary.csv)\n",
    "union-track (corroborative, coverage-dependent) rows = ", nrow(dat.ev),
    "  of which censored = ", n.ev.floor, "  non-estimable NA = ", n.ev.na,
    "\n  note: fgsea p is censored at its eps floor ", FGSEA_EPS,
    "; roast Mixed floor = ", sprintf("%.3g",
      ROAST_FLOOR), " (nrot = ",
    NROT_ROAST, ")\n",
    "per-set sensitivity rows = ", nrow(dat.sens),
    "  censored = ", sum(dat.sens$pval_censored, na.rm = TRUE),
    "\n  per-set FDR pass = ", sum(dat.sens$padj < PVAL_CUTOFF,
      na.rm = TRUE),
    " / ", nrow(dat.sens),
    "\nmeta rows (SMD) = ", nrow(dat.meta),
    "  skipped (n_study < ", MIN_N_STUDY, ") = ",
    if (length(vec.skip) == 0L) "none" else paste(vec.skip,
      collapse = ", "),
    "\nverdicts = ", paste(names(n.verdict), n.verdict, sep = ":",
      collapse = " "),
    "\nindeterminate share = ",
    round(100 * n.verdict[["indeterminate"]] / max(nrow(dat.meta),
        1L)), "%",
    "  <- §4.2 须随结果报告",
    "\nI2 range = ", sprintf("%.1f%% - %.1f%%", min(dat.meta$I2),
      max(dat.meta$I2)),
    "\nexploratory-only segments = ",
    paste(dat.meta$segment[dat.meta$n_exploratory >= dat.meta$n_study],
      collapse = ", ")))
if (nrow(dat.meta) == 0L) {
  stop(paste0("S07 FAILED: no meta-analysis row produced ; fewer than",
    "2 studies per segment."))
}
message("S07 OK. tables -> 03_results/S07/table/")
