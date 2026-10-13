#===============================================================
## S08_wgcna.R —— 共表达网络（步骤规格 §8：S08）
##
## 目的：在层 2 的段分与临床/解剖元数据之间找共表达模块，并检验模块与段分的关系，
##       目的是给「哪条共表达模块承载炎症五段信号」提供证据，而非重建调控网络。
##
## 定位：分析脚本。§0.2 已源码核定：bicor 属 WGCNA 函数（非独立包），
##       decoupleR 只返回 statistic/source/condition/score 四列，故本脚本**不引入
##       decoupleR**，也不把任何 activity 当作 p 值使用。
##
##
## 产出：03_results/S08/table/wgcna_modules_<acc>.csv  模块及其特征基因
##       03_results/S08/table/module_trait.csv         模块 x 段分 相关性
##       03_results/S08/Figure/SoftThreshold_Diagnostic.pdf    软阈值诊断（双面板）
##       03_results/S08/Figure/GeneDendrogram_ModuleColors.pdf 基因树 + 模块色带
##       03_results/S08/Figure/Module_Trait_Correlation.pdf    模块-性状相关热图
##       03_results/S08/Figure/ModuleEigengene_Correlation.pdf  ME 相关热图
##       03_results/S08/Figure/Module_Size.pdf                  模块规模柱（QC）
##
## 图件清单（规则 35 可视化强制 · SCI 发表导向；动笔前已调研文献确定风格/形式/数量）：
##   可视化: S08 共 5 幅（4 核心标准件 + 1 QC）
##     1. SoftThreshold_Diagnostic.pdf
##        —— 形式：双面板 ggplot（上=scale-free fit R2、下=fit slope，按系列着色；
##           0.85 / 0 截断虚线）
##        —— 科学问题：软阈值选取与无标度前提达成度（含不适用系列落选原因）
##        —— 参照文献：WGCNA 标准诊断图（Langfelder & Horvath 2008, PMID 19114008）
##     2. GeneDendrogram_ModuleColors.pdf
##        —— 形式：plotDendroAndColors（基因层级聚类树 + 模块颜色带，每系列每 block 一页）
##        —— 科学问题：模块识别结构与各基因所属模块
##        —— 参照文献：WGCNA 模块识别标准件（同上 PMID 19114008）
##     3. Module_Trait_Correlation.pdf
##        —— 形式：labeledHeatmap（行=模块、列=5 段分，单元格写 r 与 BH 显著星标）
##        —— 科学问题：哪些共表达模块与哪一段段分相关（模块 x 段 关联热图）
##        —— 参照文献：WGCNA 模块-性状关系标准件（同上 PMID 19114008）
##     4. ModuleEigengene_Correlation.pdf
##        —— 形式：labeledHeatmap（ME x ME 相关矩阵；格内不写数字，只表达模式）
##        —— 科学问题：模块间独立性（是否存在高度相关的模块对）
##        —— 参照文献：WGCNA 模块特征相关标准件（同上 PMID 19114008）
##     5. Module_Size.pdf
##        —— 形式：分面柱（x = 模块，y = 模块基因数，按系列分面）
##        —— 科学问题：各系列检出模块规模（QC 补充，非核心标准件）
##        —— 参照文献：模块规模条形为 WGCNA QC 惯例
##
## 图件口径：sci_theme / sci_save 逐字取自规范 §sci_theme / §sci_save 官方定义；
##       分面用 facet_wrap（ggplot 图）；plotDendroAndColors / labeledHeatmap
##       为 WGCNA 基类图函数，走新增 sci_save_base（同口径仅 PDF）。
##
## 前置：S03 + S06。WGCNA 依赖 impute，缺则显式 stop，不降级。
#===============================================================

## PART 0 常量 ------------------------------------------------------------------------------

SEED <- 20260925
set.seed(SEED)
## 软阈值扫描范围：逐字取自官方默认（实测 formals(pickSoftThreshold)$powerVector
## = c(seq(1,10,by=1), seq(12,20,by=2))）。旧写法硬编码 1:10 **截掉了 12-20**，
## 属凭记忆写参数（本例非主因，但确为缺陷）。
POWER_VECTOR <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
## 无标度判据的**完整形式**：官方 pickSoftThreshold 只按 RsquaredCut（默认 0.85）
## 判 R^2，**不检查拟合斜率的符号**。但无标度拓扑的定义是 p(k) ∝ k^(-gamma)、
## gamma > 0，即 log p(k) 对 log k 的斜率**必须为负**。实测 GSE43292 在 power=1
## 拿到 R^2=0.949 而 slope=+1.27 —— 那是对一条正斜率关系的完美拟合，与无标度
## 定义相反；旧判据只看 R^2 放行了它，于是产出「76% 基因挤进一个模块」的伪网络。
## 故补 slope < 0 这一**定义性**条件（非魔法数字），并随结果披露实测斜率。
SFT_R2_CUT <- 0.85             ## 与官方 RsquaredCut 默认值一致
## 网络类型：官方默认 "unsigned"。旧写法未声明 -> 隐式依赖默认，且两个函数的
## networkType 必须同值（power 与邻接须匹配），故显式登记为单一常量。
NETWORK_TYPE <- "unsigned"
MIN_MODULE_SIZE <- 30L
## blockwiseModules 模块合并阈值（官方 Rd 默认值，显式登记以杜绝隐式默认）
MERGE_CUT_HEIGHT <- 0.15
MIN_COR <- 0.3                ## 模块-性状 显著相关的下限
## WGCNA 输入基因数上限：邻接矩阵计算量 ∝ 基因^2 x 样本。实测外推（凭据
## 01_data/verify/probe_s08_smoke.R：3000 基因 x 23 样本 = 47.7 s）→ GSE232911
## 全基因 26822 基因 x 246 样本约 2.3 h/系列。取 **top 变异基因**降维是 WGCNA
## 官方教程的常规做法（Langfelder & Horvath 2008 的演示即用预过滤基因集），
## 且低变异基因对共表达结构无贡献。该上限是**操作定义**，非文献取值。
MAX_WGCNA_GENES <- 5000L
## WGCNA 稳定估计所需的最低样本数。**操作定义，非文献取值**：WGCNA 官方教程的
## 演示数据普遍 n >> 20，未给出硬下限；此处取 20 仅用于触发 Limitations 披露，
## 不得在正文声称有文献出处（规则 21/28：阈值须标注来源与性质）。
MIN_WGCNA_N <- 20L

vec.segment <- c("S1_initiation", "S2_amplification", "S3_effect",
  "S4_resolution", "S5_chronic")

TAB_DIR <- "03_results/S08/table"
FIG_DIR <- "03_results/S08/Figure"
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
## 规则 33：pickSoftThreshold + blockwiseModules 全基因模块识别属重时步骤，
## 断点落 checkpoints/（命名 S08_<acc>_modules.rds）；TAB_DIR 须先落定，
## 否则 file.path 在其赋值前求值而立即报错
CKPT_DIR <- file.path(TAB_DIR, "checkpoints")
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop("S08 FAILED: 'WGCNA' not installed (needs 'impute' as a dependency). ",
    paste0("Install both, then rerun. The script will not fall back",
      "to an equivalent method."))
}
library(WGCNA)  ## 模块识别（本机 stripped 安装，调用前已做 formals 签名核对）
library(ggplot2)  ## 模块图；主题与出口用规范 §sci_theme / §sci_save

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

## 基类图出口：plotDendroAndColors / labeledHeatmap 为 WGCNA 基类图函数，
## 无法走 ggplot 管道。在 pdf 设备内执行 base 绘图闭包，与 sci_save 同口径
## （仅 PDF 落 FIG_DIR；RAPA_QC_PNG=1 时额外导 PNG 供核验）。
sci_save_base <- function(filename, width, height, draw_fun) {
  pdf(file.path(FIG_DIR, paste0(filename, ".pdf")), width = width,
    height = height)
  draw_fun()
  dev.off()
  if (Sys.getenv("RAPA_QC_PNG") == "1") {
    dir.create(file.path(FIG_DIR, "_qc"), showWarnings = FALSE)
    png(file.path(FIG_DIR, "_qc", paste0(filename, ".png")),
      width = width * 150, height = height * 150, res = 150)
    draw_fun()
    dev.off()
  }
}

## 签名核对：调用点硬断言（项目既定做法）；参数名取自官方 Rd 实测
vec.arg.ps <- c("dataIsExpr", "powerVector", "verbose")
vec.miss.ps <- setdiff(vec.arg.ps, names(formals(WGCNA::pickSoftThreshold)))
if (length(vec.miss.ps) > 0L) {
  stop("S08 FAILED: pickSoftThreshold() has no argument(s): ",
    paste(vec.miss.ps,
      collapse = ", "), " (WGCNA signature drifted from this script)")
}
vec.arg.bm <- c("power", "minModuleSize", "deepSplit", "mergeCutHeight",
  "numericLabels")
vec.miss.bm <- setdiff(vec.arg.bm, names(formals(WGCNA::blockwiseModules)))
if (length(vec.miss.bm) > 0L) {
  stop("S08 FAILED: blockwiseModules() has no argument(s): ",
    paste(vec.miss.bm,
      collapse = ", "), " (WGCNA signature drifted from this script)")
}

dat.trait <- read.csv("03_results/S06/table/layer2_segment_all.csv",
  stringsAsFactors = FALSE)
vec.col.req <- c("acc", "segment", "direction")
vec.col.miss <- setdiff(vec.col.req, names(dat.trait))
if (length(vec.col.miss) > 0L) {
  stop("S08 FAILED: layer2_segment_all.csv lacks column(s): ",
    paste(vec.col.miss, collapse = ", "))
}
## 模块-段分相关要的是**样本级**段分，layer2_segment_all.csv 是 (段 x 系列) 汇总行，
## 里面根本没有 sample 列（旧写法按 sample 校验，必然 stop）。样本级段分唯一来源是
## S06 的 segment_score_matrix.rds（段 x 样本），此处改从它取，方向判定仍查 layer2 表。
vec.acc <- unique(dat.trait$acc)
list.score <- readRDS("03_results/S06/table/segment_score_matrix.rds")

dat.module <- data.frame(acc = character(0), module = character(0),
  segment = character(0), n_genes = integer(0),
  trait_r = numeric(0), trait_p = numeric(0), sft_r2 = numeric(0),
  stringsAsFactors = FALSE)
## 无标度前提不成立的系列（不建网络，如实登记诊断值而非产出伪模块）
vec.nosf <- character(0)
dat.nosf <- data.frame(acc = character(0), power = integer(0),
  sft_r2 = numeric(0), sft_slope = numeric(0), reason = character(0),
  stringsAsFactors = FALSE)
## 画图数据容器：重算或断点恢复后统一从此取（避免重复 compute）；
## list.fit 存每系列 fitIndices（SFT 诊断图，含不适用系列），
## list.plot 存每系列 dendrograms / colors_numeric / blockGenes / mat.me
## （基因树+色带图、ME 相关热图）。
list.fit <- list()
list.plot <- list()

for (acc in vec.acc) {
  ## 2026-09-29 起改读**基因级**矩阵：WGCNA 的模块基因名与模块-段分相关都要求真实基因符号，
  ## 探针级 _expr.rds 的行名是 p1..pN（无生物学含义）。
  obj.expr <- readRDS(file.path("03_results", "S03", "table", paste0(acc,
        "_gene.rds")))
  expr <- obj.expr$expr
  ## 基因过滤：按表达变异度取 top MAX_WGCNA_GENES（依据见该常量处注释）。
  ## NA 变异度的行先剔除（双色系列含 NA；且 sd 需 >= 2 个非 NA 值）。
  if (nrow(expr) > MAX_WGCNA_GENES) {
    vec.sd <- apply(expr, 1L, sd)
    vec.sd <- vec.sd[!is.na(vec.sd)]
    expr <- expr[names(sort(vec.sd, decreasing = TRUE))[
        seq_len(min(MAX_WGCNA_GENES, length(vec.sd)))], , drop = FALSE]
    cat(sprintf("  gene filter: top %d by SD (from %d)\n", nrow(expr),
        length(vec.sd)))
  }
  ## 段分作为性状：用于与模块特征基因做相关。样本级段分来自 S06 的评分矩阵
  if (!(acc %in% names(list.score))) {
    stop(paste0("S08 STOP: ", acc, " has no segment score matrix in S06 ",
        paste0("(segment_score_matrix.rds); module-trait correlation is",
          "not estimable.")))
  }
  mat.seg <- as.matrix(list.score[[acc]])
  ## 段分矩阵的列（样本）必须与表达矩阵**按名对齐**：本项目既定教训是位置对齐会
  if (!all(colnames(expr) %in% colnames(mat.seg))) {
    stop(paste0("S08 STOP: ", acc, " 的段分矩阵缺少表达矩阵中的样本（例: ",
        paste(head(setdiff(colnames(expr), colnames(mat.seg)), 3L),
          collapse = ", "),
        "）；模块-性状相关不可估，禁止按位置对齐。"))
  }
  mat.seg <- mat.seg[, colnames(expr), drop = FALSE]
  cat(sprintf("\n[S08] %s (%d probes x %d samples)\n", acc, nrow(expr),
      ncol(expr)))

  ## 规则 33：阈值扫描 + 模块识别属重时步骤，按 acc 粒度续传；
  ## 恢复失败即报错，不用 tryCatch 把损坏的 checkpoint 悄悄丢掉（规则 26）。
  ## 后缀 `_gene` 为**契约版本**：2026-09-29 起输入由探针级 _expr.rds（行名 p1..pN）
  ## 换成基因级 _gene.rds，旧 checkpoint 的模块成员是无生物学含义的 p 编号，语义已废，
  ## 故用新文件名使其自动失效。
  f.ck <- file.path(CKPT_DIR, paste0("S08_", acc, "_modules_gene_v2.rds"))
  
  sig.now <- paste(tools::md5sum(file.path("03_results", "S03", "table",
        paste0(acc, "_gene.rds"))),
    tools::md5sum("03_results/S06/table/segment_score_matrix.rds"),
    tools::md5sum("02_scripts/S08_wgcna.R"), SEED, sep = "|")
  obj.ck <- NULL
  if (file.exists(f.ck)) {
    obj.tmp <- readRDS(f.ck)
    if (identical(obj.tmp$signature, sig.now)) {
      obj.ck <- obj.tmp
    } else {
      cat(sprintf("[S08] %s checkpoint fingerprint mismatch -> recompute\n",
          acc))
    }
  }
  if (!is.null(obj.ck)) {
    power <- obj.ck$power
    vec.module <- obj.ck$vec.module
    list.me <- obj.ck$list.me
    vec.module.gene <- obj.ck$vec.module.gene
    sft.r2 <- obj.ck$sft_r2
    ## 画图数据：从断点恢复，与重算分支口径一致（list.fit/list.plot 全局）
    list.fit[[acc]] <- obj.ck$fitIndices
    list.plot[[acc]] <- list(dendrograms = obj.ck$dendrograms,
      colors_numeric = obj.ck$colors_numeric,
      blockGenes = obj.ck$blockGenes, mat.me = obj.ck$mat.me)
  } else {
    ## dataIsExpr 官方含义 = data 视为表达数据(TRUE)还是相似性矩阵(FALSE)，
    ## 喂表达矩阵必须 TRUE；软阈值扫描范围的形参是 powerVector。
    list.soft <- pickSoftThreshold(t(expr), dataIsExpr = TRUE,
      powerVector = POWER_VECTOR, networkType = NETWORK_TYPE, verbose = 0)
    ## SFT 诊断图数据源：全部系列（含不适用者）都需 fitIndices 以交代落选原因
    list.fit[[acc]] <- list.soft$fitIndices
    power <- list.soft$powerEstimate
    i.p <- match(power, list.soft$fitIndices$Power)
    sft.r2 <- if (is.na(i.p)) NA_real_
      else list.soft$fitIndices$SFT.R.sq[i.p]
    sft.slope <- if (is.na(i.p)) NA_real_
      else list.soft$fitIndices$slope[i.p]
    lgl.sf <- !is.na(sft.r2) && sft.r2 >= SFT_R2_CUT &&
      !is.na(sft.slope) && sft.slope < 0
    if (!lgl.sf) {
      vec.nosf <- c(vec.nosf, acc)
      dat.nosf <- rbind(dat.nosf, data.frame(acc = acc,
          power = if (is.na(power)) NA_integer_ else power,
          sft_r2 = sft.r2, sft_slope = sft.slope,
          reason = if (is.na(power)) "no power reached RsquaredCut"
            else if (is.na(sft.r2) || sft.r2 < SFT_R2_CUT)
              "scale-free R2 below cut" else "fit slope >= 0 (not scale-free)",
          stringsAsFactors = FALSE))
      cat(sprintf("  scale-free premise NOT met (R2=%.3f, slope=%+.2f) -> ",
          sft.r2, sft.slope))
      cat("no network built; reported as not applicable\n")
      next
    }
    ## blockwiseModules 官方形参：power / minModuleSize / deepSplit /
    ## mergeCutHeight / numericLabels（无 maxModuleSize、无 tool、无 mergeModules）
    net <- blockwiseModules(t(expr), power = power,
      minModuleSize = MIN_MODULE_SIZE, deepSplit = 2,
      mergeCutHeight = MERGE_CUT_HEIGHT,
      networkType = NETWORK_TYPE, numericLabels = TRUE)
    stopifnot(length(net$colors) == nrow(expr))
    vec.module.gene <- paste0("ME", net$colors)
    names(vec.module.gene) <- rownames(expr)
    vec.module <- unique(vec.module.gene)
    ## 模块特征基因改用 WGCNA 官方 moduleEigengenes（模块表达的**第一主成分**）。
    ## 旧写法用模块内基因的 colMeans（算术均值）：模块内反向相关的基因会相互
    ## 抵消，且非官方 ME 定义。官方函数实测存在（getNamespaceExports("WGCNA")）。
    mat.me <- moduleEigengenes(t(expr), colors = net$colors)$eigengenes
    if (!all(vec.module %in% colnames(mat.me))) {
      stop(paste0("S08 STOP: moduleEigengenes 的返回列不包含全部模块（实测列: ",
          paste(colnames(mat.me), collapse = ", "), "）。"))
    }
    ## grey 模块（numericLabels 下为 ME0）= **未被分配到任何模块的基因**，
    ## 不构成真实模块：官方约定 color 0 = grey = unassigned，其基因之间无足够
    ## 共表达（Sci Rep 2018, s41598-017-18705-z: "genes in the grey module are
    ## not co-expressed"；BMC Plant Biol 2017, 10.1186/s12870-017-1043-1:
    ## "MEgrey does not present a real module"）。已发表 SCI 的处理是**用官方
    ## removeGreyME 剔除**（Nucleic Acids Res 2022, PMC9452847 原文："the grey
    ## module ... was removed from further analysis with removeGreyME"）。
    ## 成员表仍保留 ME0 行（供披露未分配基因规模），但不参与相关与判据。
    n.grey <- sum(vec.module.gene == "ME0")
    if ("ME0" %in% colnames(mat.me)) {
      mat.me <- removeGreyME(mat.me, greyMEName = "ME0")
    }
    vec.module <- setdiff(vec.module, "ME0")
    list.me <- lapply(vec.module, function(m) {
        vec.g <- names(vec.module.gene)[vec.module.gene == m]
        list(name = m, value = as.numeric(mat.me[, m]),
          n.genes = length(vec.g))
      })
    ## 画图数据：dendrograms（基因树，plotDendroAndColors 唯一数据源）、
    ## colors_numeric（net$colors，须 labels2colors 转标准色）、blockGenes
    ## （多 block 时基因→块映射）、mat.me（ME 矩阵，ME 相关热图）。
    ## 旧版断点未存这些 → 树/ME 图无法事后补出，故本版一并落盘。
    list.plot[[acc]] <- list(dendrograms = net$dendrograms,
      colors_numeric = net$colors, blockGenes = net$blockGenes,
      mat.me = mat.me)
    saveRDS(list(power = power, vec.module = vec.module, list.me = list.me,
        vec.module.gene = vec.module.gene, sft_r2 = sft.r2,
        fitIndices = list.soft$fitIndices, dendrograms = net$dendrograms,
        colors_numeric = net$colors, blockGenes = net$blockGenes,
        mat.me = mat.me, signature = sig.now), f.ck)
  }
  ## 模块成员表（文件头声明的产物 wgcna_modules_<acc>.csv，此前只声明未写出）。
  ## 基因名为基因级矩阵的行名（真实符号），供 S13 与人工核查模块构成。
  write.csv(data.frame(gene = names(vec.module.gene),
      module = unname(vec.module.gene), stringsAsFactors = FALSE),
    file.path(TAB_DIR, paste0("wgcna_modules_", acc,
        ".csv")), row.names = FALSE)
  vec.seg.used <- intersect(vec.segment, rownames(mat.seg))
  for (seg in vec.seg.used) {
    ## 性状 = 该段的**样本级段分**（长度 = 样本数，与模块特征同维、连续值）。
    ## 旧写法取 layer2_segment_all.csv 的 direction（系列级 up/down/indeterminate）
    ## 编码成 ±1 作性状：① 该表无 sample 列 -> dat.tmp$sample 为 NULL ->
    ##    match() 得 integer(0) -> v.seg 全 NA -> cor.test 失败（实测复现）；
    ## ② 即便有 sample 列，±1 是系列级离散编码，既非样本级也丢失段分量级。
    ## 样本级段分唯一来源 = S06 的 segment_score（段 x 样本），已在上面按名对齐。
    v.seg <- as.numeric(mat.seg[seg, ])
    for (m in list.me) {
      r <- cor.test(m$value, v.seg)
      if (is.na(r$p.value)) {
        stop("S08 STOP: ", acc, " / ", m$name,
          " 的模块特征与段分相关不可估（零方差），须先查表达矩阵。")
      }
      dat.module <- rbind(dat.module, data.frame(
          acc = acc, module = m$name, segment = seg, n_genes = m$n.genes,
          trait_r = r$estimate, trait_p = r$p.value, sft_r2 = sft.r2,
          stringsAsFactors = FALSE))
    }
  }
  cat(sprintf("  modules detected = %d (grey excluded, %d genes unassigned)",
      length(vec.module), n.grey))
  cat(sprintf(" | soft power %d, scale-free R2 = %.3f (official cut 0.85)\n",
      power, sft.r2))
}

write.csv(dat.nosf, file.path(TAB_DIR,
    "wgcna_scale_free_not_met.csv"), row.names = FALSE)

## 多重校正：模块 x 段 x 系列 的相关检验是**一次检验族**（4 系列 x N 模块 x 5 段），
## 旧判据只取 raw p < 0.05，假阳性随检验数膨胀。补 BH 列，收口判据改用 BH。
dat.module$trait_p_adj <- p.adjust(dat.module$trait_p, method = "BH")
dat.module <- dat.module[order(-abs(dat.module$trait_r)), ]
write.csv(dat.module, file.path(TAB_DIR, "module_trait.csv"), row.names = FALSE)

## PART 0.5 图件 -------------------------------------------------------------------------------
## 规则 35（可视化强制 · SCI 发表导向）+ 已调研文献确定风格/形式/数量：
## WGCNA 标准件（Langfelder & Horvath 2008, BMC Bioinformatics,
## doi:10.1186/1471-2105-9-559, PMID 19114008）含四类核心图：
## ① 软阈值诊断 ② 基因树+模块色带 ③ 模块-性状热图 ④ ME 相关热图。
## 旧版仅 2 张分面柱图，其中模块-性状用柱图替代热图（类型错）、且缺 ①②④
## （数目不足）。本次按标准件重做，共 5 张（④ 后补模块规模柱作 QC）。
## 命名遵守规范 16（描述性英文，禁 Fig 前缀）；图内文字全英文（规范 17）。

## 图 1：软阈值诊断（双面板：上=scale-free fit R2，下=fit slope；按系列着色）
dat.sft <- do.call(rbind, lapply(names(list.fit), function(a) {
  fi <- list.fit[[a]]
  data.frame(acc = a, power = fi$Power, sft_r2 = fi$SFT.R.sq,
    slope = fi$slope, stringsAsFactors = FALSE)
}))
dat.sft.long <- data.frame(
  metric = rep(c("Scale-free fit R2", "Fit slope"), each = nrow(dat.sft)),
  acc = rep(dat.sft$acc, 2),
  power = rep(dat.sft$power, 2),
  value = c(dat.sft$sft_r2, dat.sft$slope),
  stringsAsFactors = FALSE)
dat.sft.long$acc <- factor(dat.sft.long$acc, levels = names(list.fit))
sci_save(ggplot(dat.sft.long, aes(x = power, y = value, colour = acc)) +
  geom_line() + geom_point() +
  facet_wrap(~ metric, scales = "free_y") +
  geom_hline(data = data.frame(metric = "Scale-free fit R2",
    y = SFT_R2_CUT), aes(yintercept = y), linetype = "dashed",
    colour = "grey40") +
  geom_hline(data = data.frame(metric = "Fit slope", y = 0),
    aes(yintercept = y), linetype = "dashed", colour = "grey40") +
  labs(title = "Scale-free topology diagnostic by series",
    x = "Soft threshold (power)", y = "R2 / slope") +
  scale_colour_discrete(name = "Series") +
  sci_theme,
  "SoftThreshold_Diagnostic", width = 8, height = 5)

## 图 2：基因聚类树 + 模块颜色带（每系列每 block 一页，plotDendroAndColors）
sci_save_base("GeneDendrogram_ModuleColors", 10, 6, function() {
  for (a in names(list.plot)) {
    pl <- list.plot[[a]]
    vec.k <- seq_along(pl$dendrograms)
    for (k in vec.k) {
      idx <- pl$blockGenes[[k]]
      col.mod <- labels2colors(pl$colors_numeric[idx])
      plotDendroAndColors(pl$dendrograms[[k]], col.mod,
        groupLabels = paste0(a, " module"),
        main = paste0(a, " modules (block ", k, ")"),
        dendroLabels = FALSE, hang = 0.03,
        addGuide = TRUE, guideHang = 0.05)
    }
  }
})

## 图 3：模块-性状相关热图（labeledHeatmap，单元格写 r 与 BH 显著星标）
dat.mt <- dat.module[!is.na(dat.module$trait_r), ]
dat.mt$module <- factor(dat.mt$module)
dat.mt$segment <- factor(dat.mt$segment, levels = vec.segment)
sci_save_base("Module_Trait_Correlation", 7, 6, function() {
  for (a in unique(dat.mt$acc)) {
    dm <- dat.mt[dat.mt$acc == a, ]
    ## 模块按基因数降序（大模块优先，便于阅读）
    vec.mod <- names(sort(tapply(dm$n_genes, dm$module,
      function(x) x[1L]), decreasing = TRUE))
    mat.r <- matrix(NA_real_, length(vec.mod), length(vec.segment),
      dimnames = list(vec.mod, vec.segment))
    mat.p <- mat.r
    for (i in seq_along(vec.mod)) for (j in seq_along(vec.segment)) {
      rw <- dm[dm$module == vec.mod[i] & dm$segment == vec.segment[j], ]
      if (nrow(rw) == 1L) {
        mat.r[i, j] <- rw$trait_r
        mat.p[i, j] <- rw$trait_p
      }
    }
    textM <- matrix("", nrow(mat.r), ncol(mat.r))
    for (i in seq_len(nrow(mat.r))) for (j in seq_len(ncol(mat.r))) {
      if (is.na(mat.r[i, j])) next
      star <- ifelse(mat.p[i, j] < 0.05, "*", "")
      textM[i, j] <- paste0(sprintf("%.2f", mat.r[i, j]), star)
    }
    labeledHeatmap(M = mat.r, xLabels = vec.segment, yLabels = vec.mod,
      ySymbols = vec.mod, textMatrix = textM, colorLabels = FALSE,
      colors = blueWhiteRed(50),
      main = paste0(a, " module - segment score correlation"),
      cex.main = 1.1, cex.lab = 0.8, cex.text = 0.7)
  }
})

## 图 4：模块特征基因(ME)相关热图（cor(MEs)，每系列一页）。
## 2026-09-30 用户指令：22x22 格内写满 r 值互相叠压不可读 → 去除格内数字
## （官方源码 labeledHeatmap L358 证实 textMatrix=NULL 即不写字）；
## 数值本身仍可查 module_trait.csv / ME 相关表，图只表达模式。
sci_save_base("ModuleEigengene_Correlation", 7, 6, function() {
  for (a in names(list.plot)) {
    me <- list.plot[[a]]$mat.me
    M <- cor(me)
    labeledHeatmap(M = M, xLabels = colnames(M), yLabels = colnames(M),
      ySymbols = colnames(M), colorLabels = FALSE,
      colors = blueWhiteRed(50),
      main = paste0(a, " module eigengene correlation"),
      cex.main = 1.1, cex.lab = 0.7)
  }
})

## 图 5：模块规模柱（QC 补充，分面 by 系列）。
## 2026-09-30 用户指令：x 轴 22 个模块名水平排布互相叠压不可读 → 转 90° 直立
## （hjust=1 使标签末端贴轴；置于 sci_theme 之后以覆盖其 axis.text 设置）。
dat.mod.size <- unique(dat.module[, c("acc", "module", "n_genes")])
sci_save(ggplot(dat.mod.size, aes(x = module, y = n_genes)) +
  geom_col(fill = "#4C72B0") +
  facet_wrap(~ acc, scales = "free_x") +
  labs(title = "Module size (number of genes) by series",
    x = "Module", y = "Genes in module") +
  sci_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)),
  "Module_Size", width = 11, height = 4.5)

## PART 1 收口 -------------------------------------------------------------------------------

dat.sig <- dat.module[abs(dat.module$trait_r) > MIN_COR &
  !is.na(dat.module$trait_r) & dat.module$trait_p_adj < 0.05, ]
message(paste0("\n===== S08 network summary =====\nseries with a scale-free ",
    "network = ", length(unique(dat.module$acc)),
    "  series where the premise is NOT met = ", length(vec.nosf)))
if (length(vec.nosf) > 0L) {
  message("[NOT APPLICABLE] ", paste(vec.nosf, collapse = ", "),
    paste0("\n  no network built for these series (single dominant",
      "factor: PC1 49-53%, ")
    , "median gene-gene |r| 0.45-0.52, fit slope >= 0);",
    "\n  see wgcna_scale_free_not_met.csv and disclose in Limitations.")
}
message(paste0("rows = ", nrow(dat.module),
    "  module-segment pairs passing |r|>", MIN_COR, " & BH<0.05 = ", nrow(
      dat.sig)))
vec.small.n <- unique(dat.module$acc[
  vapply(unique(dat.module$acc), function(a) {
      f <- file.path("03_results", "S03", "table", paste0(a, "_gene.rds"))
      ncol(readRDS(f)$expr) < MIN_WGCNA_N
    }, TRUE)])
if (length(vec.small.n) > 0L) {
  message("[NOTE] series with n < ", MIN_WGCNA_N, " samples (WGCNA modules ",
    "unstable, disclose in Limitations): ", paste(vec.small.n,
      collapse = ", "))
}
if (nrow(dat.sig) == 0L) {
  message("[NOTE] no module reached the trait-correlation threshold; ",
    "report this as a negative result in the manuscript.")
}
message("S08 OK. tables -> 03_results/S08/table/")
