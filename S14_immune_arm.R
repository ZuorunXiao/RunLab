#===============================================================
## S14_immune_arm.R —— 独立免疫臂：rapa × 血管炎症 免疫微环境富集
##
## 目的：把免疫读数与 rapa 扰动轴显式对齐，作为主线（rapa × 血管炎症）的独立
##       证据臂，弥补 S06/S07 五段级联中免疫信号被"嵌入"而非独立成章的缺口。
##       三件事：
##         ① 样本级免疫基因集 ssGSEA（C7 IMMUNESIGDB + GO BP 免疫 + KEGG 免疫 精选面板）；
##         ② rapa 轴对齐——GSE323983 用 limma 检验 Rapa / IR+Rapa vs 对照 在免疫评分上的差异；
##         ③ 独立臂：rapa 扰动轴（GSE323983 两对比）× 全 C7 IMMUNESIGDB（4872 集，最纯
##            "免疫特异性基因集"）的免疫特异性基因集富集，fgsea + cameraPR 双路 concordance，
##            独立表 + 独立图，明确对齐 rapa 扰动轴；
##         ④ 血管炎症背景富集：各系列主对比 vs 精选面板（C7 子集 + GO BP 免疫 + KEGG 免疫），
##            作为 context 段与独立臂物理分离（不同表/不同图）。
##
## 定位：分析脚本。读 S03 基因矩阵 + S05 设计/对比；产出独立免疫臂表与图。
##
##
## 产出：03_results/S14/table/immune_genesets.csv              精选免疫基因集（段/来源/成员）
##       03_results/S14/table/immune_ssgsea_scores.rds        样本 × 精选免疫集 ssGSEA 评分
# 03_results/S14/table/immune_score_rapa_limma.csv     GSE323983 rapa 轴 limma（评分级，PART 3）
# 03_results/S14/table/immune_fgsea.csv                 master：fgsea+cameraPR（rapa 全 C7 +
# 各系列精选面板）
# 03_results/S14/table/immune_specific_rapa_enrichment.csv  独立臂：仅 rapa
# 对齐免疫富集（fgsea+cameraPR 双列）
##       03_results/S14/table/immune_scRNA_rapa.csv           单细胞 rapa 方向一致性（S12）
##       03_results/S14/table/immune_arm_summary.csv          供 S13 可选层并入的汇总
##       03_results/S14/Figure/Immune_ssGSEA_Heatmap.pdf         样本×免疫集 热图（按系列+状态分面）
##       03_results/S14/Figure/Immune_Score_By_Rapa.pdf          GSE323983 免疫评分按 arm 分组箱线
# 03_results/S14/Figure/Immune_Specific_Enrichment_By_Rapa.pdf  独立臂：rapa 轴免疫特异性 GS
# 富集火山图
##       03_results/S14/Figure/Immune_Context_Enrichment.pdf      context：血管炎症背景富集火山图
##
## 图件清单：
##   可视化: S14 共 4 幅
##     1. Immune_ssGSEA_Heatmap
##        —— 形式：ggplot2 热图（样本×精选免疫集，Z 标准化，按系列+状态分面注释）
##        —— 科学问题：各样本免疫富集谱的异质性与 rapa/疾病分组结构。
##     2. Immune_Score_By_Rapa
##        —— 形式：分组箱线（x = arm，y = 免疫集 ssGSEA 评分，按精选免疫集分面）
##        —— 科学问题：rapa 是否改变免疫微环境特征（GSE323983 唯一 rapa 系列，评分级）。
##     3. Immune_Specific_Enrichment_By_Rapa  【独立臂】
##        —— 形式：火山图（x = NES，y = -log10 p，颜色 = NES，形状 = BH<0.05，按 rapa 对比分面）
##        —— 科学问题：rapa 转录响应是否富集免疫特异性基因集（C7），呈现免疫抑制方向性。
##     4. Immune_Context_Enrichment  【context 段】
##        —— 形式：火山图（x = NES，y = -log10 p，按系列主对比分面）
##        —— 科学问题：血管炎症系列广谱免疫富集背景，与独立臂物理分离。
##
## Figure legend：
##   Immune_Specific_Enrichment_By_Rapa: Gene-set enrichment (fgsea) of the rapa
##     perturbation transcriptional response against the C7 IMMUNESIGDB collection
##     (4,872 immune cell/perturbation signatures, size 10-500). Panels: Rapa vs
##     Ctrl (control), and rapa + irradiation vs irradiation. Triangles, BH
##     FDR < 0.05; color, normalized enrichment score (NES).
##   Immune_Context_Enrichment: fgsea of each vascular-inflammation cohort's
##     primary disease contrast against the curated immune panel (GO:BP immune
##     terms + KEGG (legacy) immune pathways + C7 immune signatures; 1,882 sets
##     after size filtering 10-500). Triangles, BH FDR < 0.05; color, NES.
##   Immune_ssGSEA_Heatmap: ssGSEA (alpha = 0.25) scores of the focused immune
##     panel (GO:BP/KEGG immune terms prioritized, up to 120 sets), row Z-scored
##     per set; columns grouped by cohort.
##
## 图件口径：sci_theme / sci_save 逐字取自规范 §sci_theme / §sci_save（与 S06/S09 同口径）；
##       单面板 ggplot2 出图；RAPA_QC_PNG=1 时额外导 _qc PNG 供目视核验。
##
## 前置：S03（基因矩阵）、S05（设计/对比）、S12（单细胞 rapa 方向）。
#===============================================================

## PART 0 常量与目录 ---------------------------------------------------------------

TAB_DIR <- "03_results/S14/table"
FIG_DIR <- "03_results/S14/Figure"
CKPT_DIR <- file.path(TAB_DIR, "checkpoints")
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)

library(ggplot2)   ## 免疫臂图；主题与出口用规范 §sci_theme / §sci_save
## 免疫基因集（C7 / GO BP / KEGG）
library(msigdbr)
## 样本级免疫集 ssGSEA（官方 ssgseaParam 实现，与 S06 同口径）
library(GSVA)
library(limma)     ## rapa 轴 limma + fgsea 排序 t 统计量
library(fgsea)     ## rapa 响应 DEG 免疫集富集交叉验证
library(BiocParallel)  ## fgsea 并行后端

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
    legend.title     = ggplot2::element_text(size = 11, face = "plain"))

## 出口：规范 §sci_save 官方定义，仅 PDF，落本步 FIG_DIR；RAPA_QC_PNG=1 额外导 PNG
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

MIN_GS <- 10L        ## 基因集最小成员（与 S02/S06 min.sz 一致）
MAX_GS <- 500L       ## 基因集最大成员（ssGSEA/fgsea 上限，防超大集淹没信号）
ALPHA_SSGSEA <- 0.25 ## 官方 ssgseaParam 默认 alpha（与 S06 同口径）

## rapa 轴唯一系列 + 其余血管炎症系列（与 S06 四系列口径一致）
vec.acc <- c("GSE232911", "GSE57691", "GSE43292", "GSE323983")
## GSE323983 arm 水平显式锁定（避免 factor 参考水平随机）
lvl.arm <- c("Ctrl", "Rapa", "IR", "IR with Rapa")
## 免疫关键词（GO BP / KEGG 名称过滤，大小写不敏感；覆盖细胞型/状态/通路）
kw.immune <- c("immune", "immun", "cytokine", "interleukin", "interferon",
  "t cell", "t-cell", "b cell", "b-cell", "macrophage", "neutrophil",
  "dendritic", "nk cell", "nk-cell", "monocyte", "leukocyte", "antigen",
  "mhc", "tcr", "chemokine", "inflammat", "lymphocyte", "granulocyte",
  "mast cell", "eosinophil", "il-", "tnf", "ifn")
## 聚焦热图面板关键词（更窄，保证可读）与截断上限
kw.heat <- c("b cell", "t cell", "cd4", "cd8", "nk cell", "macrophage",
  "dendritic", "monocyte", "neutrophil", "eosinophil", "mast cell",
  "cytokine", "interferon", "interleukin", "chemokine", "adaptive immune",
  "innate immune", "inflammatory response", "t cell activation",
    "b cell activation")
rep.terms <- c("B cell", "T cell", "CD4", "CD8", "NK", "Macrophage",
  "Dendritic", "Monocyte", "Neutrophil", "Interferon", "Interleukin",
    "Cytokine")
MAX_HEAT <- 120L   ## 热图面板截断上限

## PART 1 免疫基因集锁定--------------------------------

cat(sprintf(paste0("\n[PART 1] fetching immune gene sets via msigdbr",
  "(species = Homo sapiens)\n")))
dat.ms <- msigdbr::msigdbr(species = "Homo sapiens")
## 全 C7 IMMUNESIGDB（fgsea rapa 用，不一并进 ssGSEA 热图以免不可读）
dat.c7 <- dat.ms[dat.ms$gs_collection == "C7" &
    dat.ms$gs_subcollection == "IMMUNESIGDB", ]
nm.norm <- tolower(gsub("_", " ", dat.ms$gs_name))
kw.vec <- paste(tolower(kw.immune), collapse = "|")
is.immune <- grepl(kw.vec, nm.norm)
## 精选面板：GO:BP 免疫 + KEGG(LEGACY) 免疫 + C7 免疫签名（ssGSEA 热图 + context fgsea 用）
dat.cur <- dat.ms[((dat.ms$gs_collection == "C5" &
      dat.ms$gs_subcollection == "GO:BP") |
    (dat.ms$gs_collection == "C2" & dat.ms$gs_subcollection ==
      "CP:KEGG_LEGACY") |
    (dat.ms$gs_collection == "C7")) & is.immune, ]

## 去重基因集（同名可能跨子库重复），保留来源标签
BuildSetList <- function(dat.gs) {
  lst <- split(dat.gs$gene_symbol, dat.gs$gs_name)
  lst <- lst[sapply(lst, function(g) {
    n <- length(unique(g)); n >= MIN_GS && n <= MAX_GS })]
  lst
}
list.c7 <- BuildSetList(dat.c7)                 ## 全 C7（rapa fgsea）
list.cur <-
  BuildSetList(dat.cur)               ## 精选面板（各系列主对比 fgsea）
PickRepSets <- function(terms, names.all) {
  nm.sp <- tolower(gsub("_", " ", names.all)) 
  out <- character(0)
  for (t in terms) {
    hit <- grep(tolower(t), nm.sp, fixed = FALSE)
    if (length(hit) >= 1L) out <- c(out, names.all[hit[1L]])
  }
  out
}
dat.heat <- dat.cur[grepl(paste(tolower(kw.heat), collapse = "|"),
  tolower(gsub("_", " ", dat.cur$gs_name))), ]   
list.heat <- BuildSetList(dat.heat)
is.gk.heat <- !grepl("^GSE[0-9]+_", names(list.heat))
list.heat <- list.heat[order(!is.gk.heat)]
if (length(list.heat) > MAX_HEAT) list.heat <- list.heat[seq_len(MAX_HEAT)]
PrettyLabel <- function(nm) {
  u <- unique(nm)
  s <- sub("^GSE[0-9]+_|^GOBP_|^KEGG_", "", u)
  s <- gsub("_", " ", s)
  dup <- duplicated(s) | duplicated(s, fromLast = TRUE)
  lut <- ifelse(dup, u, s)
  names(lut) <- u
  as.character(lut[nm])
}
cat(sprintf("  C7 IMMUNESIGDB sets (size %d-%d) = %d\n", MIN_GS, MAX_GS,
  length(list.c7)))
cat(sprintf("  curated immune panel (GO/KEGG/C7-cell) = %d\n",
  length(list.cur)))
cat(sprintf("  focused heatmap panel = %d (capped at %d)\n", length(list.heat),
  MAX_HEAT))
rep.sets <- PickRepSets(rep.terms, names(list.heat))
stopifnot(length(rep.sets) >= 5L)   
cat(sprintf("  representative cell-type sets for boxplot = %d\n",
  length(rep.sets)))

## 精选面板落盘（段/来源/成员），供复核与 S13 并入
dat.gset.out <- do.call(rbind, lapply(names(list.cur), function(nm) {
  src <- unique(dat.cur$gs_collection[dat.cur$gs_name == nm])
  data.frame(gs_name = nm, source = paste(src, collapse = ";"),
    n_gene = length(list.cur[[nm]]),
    members = paste(unique(list.cur[[nm]]), collapse = ";"),
    stringsAsFactors = FALSE)
}))
write.csv(dat.gset.out, file.path(TAB_DIR, "immune_genesets.csv"),
  row.names = FALSE)
cat(sprintf("  immune_genesets.csv written: %d sets\n", nrow(dat.gset.out)))

## PART 2 样本级免疫 ssGSEA--------------------

cat(sprintf(paste0("\n[PART 2] sample-level immune ssGSEA (per series,",
  "checkpointed)\n")))
list.score.all <- list()
for (acc in vec.acc) {
  ckpt <- file.path(CKPT_DIR, paste0("immune_ssgsea_", acc, ".rds"))
  if (file.exists(ckpt)) {
    list.score.all[[acc]] <- readRDS(ckpt)
    cat(sprintf("  %-11s checkpoint reused\n", acc))
    next
  }
  file.expr <- file.path("03_results", "S03", "table", paste0(acc, "_gene.rds"))
  if (!file.exists(file.expr)) {
    stop(paste0("S14 STOP: gene matrix for ", acc, " missing (run S03)."))
  }
  mat.expr <- readRDS(file.expr)$expr
  if (!is.matrix(mat.expr) || is.null(rownames(mat.expr))) {
    stop(paste0("S14 STOP: ", acc, " expr matrix lacks gene-symbol rownames."))
  }
  ## 聚焦面板求交表达矩阵（ssGSEA 以基因符号为键）
  vec.keep <- sapply(list.heat,
    function(g) sum(g %in% rownames(mat.expr)) >= MIN_GS)
  list.use <- list.heat[vec.keep]
  if (length(list.use) == 0L) {
    stop(paste0("S14 STOP: ", acc,
      " no focused immune set overlaps expr matrix."))
  }
  obj.ss <- GSVA::gsva(GSVA::ssgseaParam(mat.expr, list.use,
    alpha = ALPHA_SSGSEA, normalize = TRUE))
  stopifnot(ncol(obj.ss) == ncol(mat.expr))
  rownames(obj.ss) <- names(list.use)
  colnames(obj.ss) <- colnames(mat.expr)
  saveRDS(obj.ss, ckpt)
  list.score.all[[acc]] <- obj.ss
  cat(sprintf("  %-11s ssGSEA done: %d sets x %d samples\n", acc,
    nrow(obj.ss), ncol(obj.ss)))
}
## 合并为长表 + 样本注释（系列 + rapa/disease 状态）
dat.score.long <- do.call(rbind, lapply(names(list.score.all), function(a) {
  m <- as.matrix(list.score.all[[a]])
  data.frame(acc = a, immune_set = rep(rownames(m), times = ncol(m)),
    sample = rep(colnames(m), each = nrow(m)),
    score = as.numeric(m), stringsAsFactors = FALSE)
}))
saveRDS(list.score.all, file.path(TAB_DIR, "immune_ssgsea_scores.rds"))
cat(sprintf("  combined ssGSEA long table: %d rows\n", nrow(dat.score.long)))

## PART 3 rapa 轴 limma（GSE323983：免疫评分 ~ arm）-------------------------------

cat(sprintf("\n[PART 3] rapa-axis limma on immune scores (GSE323983)\n"))
file.design <- file.path("03_results", "S05", "table",
  "design_used_GSE323983.csv")
if (!file.exists(file.design)) {
  stop("S14 STOP: GSE323983 design missing (run S05).")
}
dat.design <- read.csv(file.design, stringsAsFactors = FALSE)
dat.design$arm <- factor(gsub(" ", "", dat.design$arm),
  levels = c("Ctrl", "Rapa", "IR", "IRwithRapa"))  ## 去空格保证系数名合法
mat.score <- list.score.all[["GSE323983"]]
stopifnot(all(dat.design$sample_id %in% colnames(mat.score)))  ## 设计样本必须全在评分矩阵
mat.score <- mat.score[, dat.design$sample_id,
  drop = FALSE]  ## 按设计顺序显式重排
des <- stats::model.matrix(~ donor + arm, data = dat.design)
fit <- limma::lmFit(mat.score, des)
ct <- limma::makeContrasts(
  Rapa_vs_Ctrl = armRapa,
  RapaIR_vs_IR = armIRwithRapa - armIR,
  levels = des)
fit.c <- limma::eBayes(limma::contrasts.fit(fit, ct))
dat.rapa <- do.call(rbind, lapply(c("Rapa_vs_Ctrl", "RapaIR_vs_IR"),
  function(cn) {
    tt <- limma::topTable(fit.c, coef = cn, number = Inf, sort.by = "none")
    data.frame(contrast = cn, immune_set = rownames(tt),
      logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
      P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val,
      stringsAsFactors = FALSE)
  }))
write.csv(dat.rapa, file.path(TAB_DIR, "immune_score_rapa_limma.csv"),
  row.names = FALSE)
n.sig <- sum(dat.rapa$adj.P.Val < 0.05)
cat(sprintf(paste0("GSE323983 rapa limma: %d immune sets x 2 contrasts;",
  "BH<0.05 = %d\n"),
  length(unique(dat.rapa$immune_set)), n.sig))

## PART 4 富集：rapa 扰动轴（独立免疫臂，fgsea + cameraPR 双路）+ 血管炎症背景（context）

cat(sprintf("\n[PART 4] enrichment (SerialParam; set.seed before each call)\n"))
bp <-
  BiocParallel::SerialParam()  
dat.fg <- data.frame(acc = character(0), scope = character(0),
  gs_name = character(0), n_gene = integer(0),
  pval = numeric(0), padj = numeric(0), NES = numeric(0),
  cam_n_gene = integer(0), cam_direction = numeric(0), cam_FDR = numeric(0),
  stringsAsFactors = FALSE)

## 富集执行器：fgsea（排序 t 统计自检）+ cameraPR（竞争性基因集检验）双路 concordance
RunEnrich <- function(acc, term, num, den, list.sets, scope) {
  file.expr <- file.path("03_results", "S03", "table", paste0(acc, "_gene.rds"))
  mat.expr <- readRDS(file.expr)$expr
  file.design <- file.path("03_results", "S05", "table",
    paste0("design_used_", acc, ".csv"))
  dat.d <- read.csv(file.design, stringsAsFactors = FALSE)
  if (!(term %in% colnames(dat.d))) {
    stop(paste0("S14 STOP: term '", term, "' absent in ", acc, " design."))
  }
  keep <- dat.d[[term]] %in% c(num, den)
  dat.s <- dat.d[keep, ]
  stopifnot(all(dat.s$sample_id %in% colnames(mat.expr)))  
  ## 去空格保证 grp 系数名合法（如 "large AAA" -> "largeAAA"，与对比名构造一致）
  grp <- factor(gsub(" ", "", dat.s[[term]]),
    levels = gsub(" ", "", c(num, den)))
  mat.s <- mat.expr[, as.character(dat.s$sample_id), drop = FALSE]
  des <- stats::model.matrix(~ 0 + grp)
  ct.name <- paste0("grp", gsub(" ", "", num), " - grp", gsub(" ", "", den))
  fit <- limma::lmFit(mat.s, des)
  fit.c <- limma::eBayes(limma::contrasts.fit(fit,
    limma::makeContrasts(contrasts = ct.name, levels = des)))
  tstat <- setNames(fit.c$t[, 1], rownames(mat.s))
  tstat <- sort(tstat, decreasing = TRUE)
  set.seed(123L)   ## 可复现（fgsea 蒙特卡洛）
  res <- fgsea::fgsea(pathways = list.sets, stats = tstat,
    minSize = MIN_GS, maxSize = MAX_GS, BPPARAM = bp)
  res$gs_name <- names(list.sets)[match(res$pathway, names(list.sets))]
  ## cameraPR：以表达矩阵 symbol 索引匹配基因集，竞争性检验做第二路证据
  idx <- lapply(list.sets, function(s) which(rownames(mat.s) %in% s))
  ov <- lengths(idx) >= 1L
  if (!all(ov)) {
    cat(sprintf(paste0("[warn] %s: %d sets have 0 genes in %s matrix (cameraPR",
      "skipped)\n"),
      scope, sum(!ov), acc))   
  }
  cm <- limma::makeContrasts(contrasts = ct.name, levels = des)
  cam <- limma::cameraPR(mat.s, index = idx[ov], design = des,
    contrast = as.numeric(cm[, 1]))
  cam$gs_name <- names(list.sets)[ov]
  res.m <- merge(res[, c("gs_name", "size", "pval", "padj", "NES")],
    cam[, c("gs_name", "NGenes", "Direction", "FDR")], by = "gs_name",
    all.x = TRUE)
  data.frame(acc = acc, scope = scope, gs_name = res.m$gs_name,
    n_gene = res.m$size, pval = res.m$pval, padj = res.m$padj, NES = res.m$NES,
    cam_n_gene = res.m$NGenes, cam_direction = res.m$Direction,
    cam_FDR = res.m$FDR, stringsAsFactors = FALSE)
}

dat.spec <- read.csv(file.path("03_results", "S05", "table",
  "contrast_spec.csv"), stringsAsFactors = FALSE)
## 4a 独立臂：rapa 扰动轴 × 全 C7（4872 集，最纯免疫特异性基因集富集，对齐 rapa 轴）
dat.fg <- rbind(dat.fg, RunEnrich("GSE323983", "arm", "Rapa", "Ctrl",
  list.c7, "rapa_Rapa_vs_Ctrl_C7"))
dat.fg <- rbind(dat.fg, RunEnrich("GSE323983", "arm", "IR with Rapa", "IR",
  list.c7, "rapa_RapaIR_vs_IR_C7"))
## 4b 血管炎症背景富集（context 段，与独立臂物理分离）：各系列主对比 vs 精选面板
for (i in seq_len(nrow(dat.spec))) {
  if (!dat.spec$primary[i]) next
  acc <- dat.spec$acc[i]; term <- dat.spec$term[i]
  num <- dat.spec$num[i]; den <- dat.spec$den[i]
  if (!acc %in% vec.acc) next
  dat.fg <- rbind(dat.fg, RunEnrich(acc, term, num, den,
    list.cur, paste0("primary_", acc)))
  cat(sprintf("  fgsea+cameraPR %-11s %s_vs_%s done\n", acc, num, den))
}
write.csv(dat.fg, file.path(TAB_DIR, "immune_fgsea.csv"), row.names = FALSE)

## 4c 独立臂物理分离：仅 rapa 对齐的免疫特异性富集（fgsea + cameraPR 双列），独立表交付
dat.arm <- dat.fg[grepl("^rapa_", dat.fg$scope), ]
stopifnot(nrow(dat.arm) > 0)   
write.csv(dat.arm, file.path(TAB_DIR, "immune_specific_rapa_enrichment.csv"),
  row.names = FALSE)
n.fgsea.arm <- sum(dat.arm$padj < 0.05, na.rm = TRUE)
n.cam.arm <- sum(dat.arm$cam_FDR < 0.05, na.rm = TRUE)
cat(sprintf(paste0("[ARM] immune-specific GS enrichment x rapa axis: %d",
  "sets; fgsea BH<0.05=%d; cameraPR FDR<0.05=%d\n"),
  nrow(dat.arm), n.fgsea.arm, n.cam.arm))
cat(sprintf("  immune_fgsea.csv written: %d rows; BH<0.05 = %d\n",
  nrow(dat.fg), sum(dat.fg$padj < 0.05, na.rm = TRUE)))

## PART 5 单细胞 rapa 方向一致性-------------------------

cat(sprintf("\n[PART 5] single-cell rapa direction concordance (S12)\n"))
file.sc <- file.path("03_results", "S12", "table",
  "gse212388_rapa_direction.csv")
file.ct <- file.path("03_results", "S12", "table", "segment_by_celltype.csv")
dat.sc <- if (file.exists(file.sc)) {
  read.csv(file.sc, stringsAsFactors = FALSE)
} else {
  stop("S14 STOP: S12 rapa direction missing (run S12).")
}
dat.ct <- if (file.exists(file.ct)) {
  read.csv(file.ct, stringsAsFactors = FALSE)
} else {
  stop("S14 STOP: S12 segment_by_celltype missing (run S12).")
}
## scRNA rapa 方向（GSE212388 Rapamycin_vs_Control）逐段，全部 indeterminate（已记）
dat.sc.sum <- data.frame(
  acc = unique(dat.sc$acc),
  arm = unique(dat.sc$arm),
  n_segment = nrow(dat.sc),
  n_indeterminate = sum(dat.sc$direction == "indeterminate"),
  min_BH = min(dat.sc$min.BH, na.rm = TRUE),
  stringsAsFactors = FALSE)
## 免疫相关细胞类型（segment_by_celltype）的段评分均值，作为免疫微环境单细胞侧面
vec.immune.ct <- c("B_cell", "T_cell", "NK_cell", "Macrophage",
  "Dendritic", "Monocyte", "Neutrophil", "Endothelial", "Smooth_muscle")
dat.ct.imm <- dat.ct[dat.ct$celltype %in% vec.immune.ct, ]
dat.ct.imm.sum <- aggregate(mean_score ~ celltype, data = dat.ct.imm,
  FUN = mean)
write.csv(dat.sc.sum, file.path(TAB_DIR, "immune_scRNA_rapa.csv"),
  row.names = FALSE)
cat(sprintf("  scRNA rapa: %d segments, %d indeterminate, min BH = %.3f\n",
  dat.sc.sum$n_segment, dat.sc.sum$n_indeterminate, dat.sc.sum$min_BH))
cat(sprintf("  immune cell types summarized: %d\n", nrow(dat.ct.imm.sum)))

MakeShortLabel <- function(x) {
  lab <- gsub("^GSE[0-9]+_", "",
    x)            ## 剥免疫签名来源研究的 GEO 编号前缀
  lab <- gsub("^(GO|KEGG|REACTOME)_", "", lab)
  lab <- gsub("_UP$", " (up)", lab)
  lab <- gsub("_DN$", " (down)", lab)
  lab <- gsub("_VS_", " vs ", lab, fixed = TRUE)
  lab <- gsub("TCELL", " T cell", lab, fixed = TRUE)
  lab <- gsub("BCELL", " B cell", lab, fixed = TRUE)
  lab <- gsub("NKCELL", " NK cell", lab, fixed = TRUE)
  lab <- tolower(gsub("_", " ", lab, fixed = TRUE))
  lab <- gsub("cd([48])", "CD\\1", lab)
  lab <- gsub("\\btcr\\b", "TCR", lab)
  lab <- gsub("\\bnk\\b", "NK", lab)
  lab <- gsub("\\btnf\\b", "TNF", lab)
  lab <- gsub("\\bmhc\\b", "MHC", lab)
  lab <- gsub("ifn([a-z])", "IFN-\\1", lab)
  lab <- gsub("\\bil([0-9])", "IL-\\1", lab)
  lab <- gsub("mcsf", "M-CSF", lab)
  lab <- gsub("\\btreg\\b", "Treg", lab)
  substr(lab, 1, 1) <- toupper(substr(lab, 1, 1))
  lab
}
lvl.cat <- c("CD4 T cell", "CD8 T cell", "T cell", "B cell", "NK cell",
  "Dendritic cell", "Macrophage / Monocyte", "Granulocyte",
  "Cytokine / IFN response", "Other immune")
ClassifyImmuneSet <- function(x) {
  x.l <- tolower(x)
  out <- rep(NA_character_, length(x))
  out[grepl("cd4", x.l)] <- "CD4 T cell"
  hit <- is.na(out) & grepl("cd8", x.l); out[hit] <- "CD8 T cell"
  hit <- is.na(out) & grepl("tcell|t[-_ ]cell|tcr|treg|th[0-9]|tfh", x.l)
  out[hit] <- "T cell"
  hit <- is.na(out) & grepl("bcell|b[-_ ]cell", x.l); out[hit] <- "B cell"
  hit <- is.na(out) & grepl("nk[-_ ]?cell|\\bnk\\b", x.l); out[hit] <- "NK cell"
  hit <- is.na(out) & grepl("dendritic", x.l); out[hit] <- "Dendritic cell"
  hit <- is.na(out) & grepl("macrophage|monocyte|mcsf", x.l)
  out[hit] <- "Macrophage / Monocyte"
  hit <- is.na(out) & grepl("neutrophil|granulocyte|eosinophil|basophil|mast",
    x.l)
  out[hit] <- "Granulocyte"
  hit <- is.na(out) & grepl(paste(c("interferon", "ifn", "interleukin",
    "\\bil[-_ ]?[0-9]", "tnf", "chemokine", "cytokine"), collapse = "|"),
    x.l)
  out[hit] <- "Cytokine / IFN response"
  out[is.na(out)] <- "Other immune"
  out
}

## PART 6 图件 -----------------------------------------------------------------------

cat(sprintf("\n[PART 6] figures\n"))
## 图 1：ssGSEA 热图（样本×精选免疫集，每集跨样本 Z 标准化；按系列+状态注释）
dat.hm <- dat.score.long
dat.hm$acc <- factor(dat.hm$acc, levels = vec.acc)
## 状态注释：GSE323983 按 arm（rapa/ctrl/ir），其余按系列
dat.ann <- data.frame(sample = character(0), state = character(0),
  stringsAsFactors = FALSE)
for (acc in vec.acc) {
  if (acc == "GSE323983") {
    dd <- read.csv(file.path("03_results", "S05", "table",
      "design_used_GSE323983.csv"), stringsAsFactors = FALSE)
    dd$arm <- factor(dd$arm, levels = lvl.arm)
    dd$state <- as.character(dd$arm)
    dd$sample <- dd$sample_id   ## 与设计表对齐：列名 sample_id -> sample
  } else {
    dd <- data.frame(sample = unique(dat.score.long$sample[dat.score.long$acc ==
      acc]),
      state = acc, stringsAsFactors = FALSE)
  }
  dat.ann <- rbind(dat.ann, dd[, c("sample", "state")])
}
dat.hm <- merge(dat.hm, dat.ann, by = "sample", all.x = TRUE)
stopifnot(!any(is.na(dat.hm$state)))   ## 注释对齐断言：缺注释即停
## 样本按 系列->状态->样本 排序并固化为因子水平，保证分面内顺序稳定可读
dat.hm <- dat.hm[order(dat.hm$acc, dat.hm$state, dat.hm$sample), ]
dat.hm$sample <- factor(dat.hm$sample, levels = unique(dat.hm$sample))

dat.hm$z <- ave(dat.hm$score, dat.hm$immune_set,
  FUN = function(x) as.numeric(scale(x)))
p.heat <- ggplot(dat.hm, aes(x = sample, y = immune_set, fill = z)) +
  geom_tile() +
  facet_grid(~ acc, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, name = "Z score") +
  scale_x_discrete(breaks = NULL) +
  scale_y_discrete(labels = function(x) PrettyLabel(x)) +
  labs(title = "Immune ssGSEA profile per sample (Z-score)",
    x = "Sample (grouped by series)", y = "Immune gene set") +
  sci_theme +
  theme(axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 8),
    panel.spacing.x = unit(0.15, "lines"),
    plot.title = element_text(hjust = 0.5))
sci_save(p.heat, "Immune_ssGSEA_Heatmap", width = 18, height = 12)

## 图 2：GSE323983 免疫评分按 arm 分组箱线（rapa 轴核心图，代表性细胞型）
rep.avail <- rep.sets[rep.sets %in% rownames(list.score.all[["GSE323983"]])]
stopifnot(length(rep.avail) >= 5L)   ## 该系列可用代表集过少即停，不静默降级
dat.bx <- merge(dat.score.long[dat.score.long$acc == "GSE323983" &
    dat.score.long$immune_set %in% rep.avail, ], dat.ann, by = "sample",
  all.x = TRUE)
stopifnot(all(rep.avail %in% dat.bx$immune_set))   ## 代表性集必须全在评分矩阵
dat.bx$state <- factor(dat.bx$state, levels = lvl.arm)
dat.bx$immune_set <- factor(PrettyLabel(dat.bx$immune_set),
  levels = PrettyLabel(rep.avail))
p.bx <- ggplot(dat.bx, aes(x = state, y = score, fill = state)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(~ immune_set, scales = "free_y") +
  scale_fill_manual(values = c("Ctrl" = "#4477AA", "Rapa" = "#EE6677",
    "IR" = "#228833", "IR with Rapa" = "#CCBB44"), name = "Arm") +
  labs(title = "Immune ssGSEA score by rapa arm (GSE323983)",
    x = "Arm", y = "ssGSEA score") +
  sci_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 7))
sci_save(p.bx, "Immune_Score_By_Rapa", width = 18, height = 8)

## 图 3（独立臂）：免疫特异性基因集富集 × rapa 扰动轴（fgsea NES 火山图，按 rapa 对比分面）
dat.arm.fig <- dat.fg[grepl("^rapa_", dat.fg$scope), ]
stopifnot(nrow(dat.arm.fig) > 0)
dat.arm.fig$contrast <- gsub("rapa_|_C7$", "", dat.arm.fig$scope)
dat.arm.fig$sig <- dat.arm.fig$padj < 0.05
dat.arm.fig$nlp <- -log10(pmax(dat.arm.fig$pval, 1e-300))
p.arm <- ggplot(dat.arm.fig, aes(x = NES, y = nlp)) +
  geom_point(aes(color = NES, shape = sig), size = 1.6, alpha = 0.75) +
  facet_wrap(~ contrast) +
  scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, name = "NES") +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16),
    breaks = c(TRUE, FALSE), labels = c("TRUE", "FALSE"), name = "BH < 0.05") +
  labs(title = paste0("Immune-specific gene set enrichment along the rapa",
    "perturbation axis"),
    x = "Normalized enrichment score (NES)", y = "-log10(p)") +
  sci_theme
sci_save(p.arm, "Immune_Specific_Enrichment_By_Rapa", width = 10, height = 6)

## 图 4（context 段）：血管炎症背景富集（各系列主对比 vs 精选面板），与独立臂物理分离
dat.ctx <- dat.fg[grepl("^primary_", dat.fg$scope), ]
if (nrow(dat.ctx) > 0) {
  dat.ctx$series <- gsub("^primary_", "", dat.ctx$scope)
  dat.ctx$sig <- dat.ctx$padj < 0.05
  dat.ctx$nlp <- -log10(pmax(dat.ctx$pval, 1e-300))
  p.ctx <- ggplot(dat.ctx, aes(x = NES, y = nlp)) +
    geom_point(aes(color = NES, shape = sig), size = 1.4, alpha = 0.7) +
    facet_wrap(~ series) +
    scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, name = "NES") +
    scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16),
      breaks = c(TRUE, FALSE), labels = c("TRUE", "FALSE"),
        name = "BH < 0.05") +
    labs(title = "Vascular-inflammation immune enrichment (context)",
      x = "Normalized enrichment score (NES)", y = "-log10(p)") +
    sci_theme
  sci_save(p.ctx, "Immune_Context_Enrichment", width = 12, height = 6)
}

## PART 7 汇总（供 S13 可选层并入）---------------------------------------------------

cat(sprintf("\n[PART 7] summary for S13\n"))
dat.sum <- data.frame(
  arm = "immune (independent)",
  n_immune_sets = length(list.cur),
  n_c7_sets = length(list.c7),
  rapa_limma_BH = n.sig,
  rapa_arm_fgsea_BH = n.fgsea.arm,
  rapa_arm_cameraPR_BH = n.cam.arm,
  fgsea_BH = sum(dat.fg$padj < 0.05, na.rm = TRUE),
  scRNA_rapa_n_segment = nrow(dat.sc),
  scRNA_rapa_indeterminate = sum(dat.sc$direction == "indeterminate"),
  scRNA_rapa_min_BH = min(dat.sc$min.BH, na.rm = TRUE),
  stringsAsFactors = FALSE)
stopifnot(nrow(dat.sum) == 1L)
write.csv(dat.sum, file.path(TAB_DIR, "immune_arm_summary.csv"),
  row.names = FALSE)
cat(sprintf("\n===== S14 immune arm summary =====\n"))
print(dat.sum)
cat("S14 OK. tables + figures -> 03_results/S14/\n")
