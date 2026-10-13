#===============================================================
## S03_normalize.R —— 分平台归一化（步骤规格 §8：S03）
##
## 目的：把 S01 登记的落盘文件按「每个系列各自应走的归一化路径」转成统一的
##       表达矩阵（探针 x 样本），并对每个系列做尺度与结构断言（Q2），
##       逐条登记所用 CDF 平台 / 归一化方法，供 S04 以后各步复用。
##
## 定位：分析脚本。只做读取 + 归一化 + 断言，不做差异检验。
# 产出：03_results/S03/table/<acc>_expr.rds      探针级表达矩阵（行 p1..pN / 列 s1..sN，原契约不动）
##       03_results/S03/table/<acc>_gene.rds      基因级表达矩阵（行 = 基因符号 / 列 = 真实样本 ID）
##       03_results/S03/table/normalization_log.csv  每系列所走的路径与方法
##       03_results/S03/table/normalization_qc.csv   Q2/Q3 断言（探针级 4 项 + 基因级 4 项）
##       03_results/S03/table/gene_mapping_log.csv   探针 -> 符号的映射与归并统计
# 03_results/S03/table/panel_coverage.csv     12 基因面板逐系列逐基因存在性（§9 Q4
# 记录，不做同义替换）
# 03_results/S03/Figure/Normalization_Metrics_By_Series.pdf 单面板柱：保留探针 / 已映射探针 /
# 唯一符号
# 03_results/S03/Figure/Expression_Range_Per_Series.pdf     单面板量程图：各系列表达 min-max（按 log
# 尺度着色）
# 03_results/S03/Figure/Normalization_QC_By_Series.pdf      单面板热图：系列 x 断言（探针级 +
# 基因级）
# 03_results/S03/Figure/Panel_Gene_Coverage_By_Series.pdf   单面板热图：系列 x 12 面板基因存在性
##
## 基因级矩阵口径（2026-09-29 新增，用户裁定「保留原有矩阵 + 新增基因级矩阵」）：
##   - 行轴：平台注释给出的基因符号（GSE232911 用自带的 Gene.Symbol 列）；
##     一个探针给出的多个符号**全部展开**（`///` 或 `;` 分隔），再以
##     `limma::avereps` 按符号取均值（官方实现为 rowsum(na.rm=TRUE) / 非 NA 计数，
##     故部分缺失的探针不会把整符号拉成 NA；全 NA 符号产生 NaN，逐一剔除）。
##     实测见 01_data/verify/check_avereps_signature.R。
##   - 列轴：还原真实样本身份。各 series_matrix 的表头本就是 GSM 号，是旧代码把它
##     覆盖成了 s1..sN；GSE232911 的矩阵列是样本 title，由 S01 的
##     gse232911_grouping.tsv 映射回 GSM；GSE57691 / GSE323983 用 counts 文件自带
##     样本列名。探针级 `<acc>_expr.rds` 一律不改，其 p1..pN / s1..sN 契约原样保留。
##   - §9 Q4 明确「缺失逐条记录，禁止同义替换」，故本脚本**不做**别名改写
##     （CXCL8 在 GPL17077 / GPL17586 上按旧符号 IL8 收录，一律如实记录为缺失）。
##
## 图件清单（规则 35 可视化强制 · SCI 发表导向；动笔前已调研已发表文献确定风格/形式/数量）：
##   可视化: S03 共 4 幅
##     1. Normalization_Metrics_By_Series.pdf
# —— 形式：单面板分组柱状图（x = 系列，y = 特征数，fill = 保留探针 / 已映射探针 / 唯一符号），
##           x 轴标签 45° 旋转
##        —— 回答的科学问题：各系列从「归一化后保留探针」到「唯一基因符号」的逐级损耗，
##           以及基因级矩阵规模是否达 Q3 下限（>= 1000）、跨系列是否可比。
##        —— 参照文献：归一化后探针/基因保留数条形图为 RNA-seq/microarray QC 标准件
##           （limma/edgeR 工作流 QC 惯例）；探针到基因的逐级计数是注释层损耗的标准呈现。
##     2. Expression_Range_Per_Series.pdf
# —— 形式：单面板量程图（x = 系列，每系列一条 min-max 线段 + 端点，按已登记的归一化路径着色）
##        —— 回答的科学问题：各系列表达量程与尺度形态是否合理（原始强度 vs 已 log2 的量程可直接
##           从跨度读出），并标注该系列走的是 identity 还是 voom。
##        —— 参照文献：样本量程/min-max 图用于归一化诊断（MA plot 与表达分布箱线图为同类 QC 惯例）。
##        —— 口径说明（2026-09-29 收紧）：**不**用「推断出的尺度」着色——旧版把
##           `0.1% 分位 < 0` 当作 log2 判定并按它着色，而 11 系列中有 9 系列下限为正
##           （GPL570 系列实测量程 1.7–14.4，显然是 log2），该列被误读为「不是 log2」。
##           现改为按**已登记**的归一化路径着色，尺度只以实测分位数登记在 normalization_log.csv。
##     3. Normalization_QC_By_Series.pdf
##        —— 形式：单面板热图（x = 系列，y = 断言，tile 着色 PASS/FAIL）
##        —— 回答的科学问题：各系列 8 项断言（探针级 4 项 + 基因级 4 项）逐项通过情况。
##        —— 参照文献：断言状态矩阵热图为 QC 汇总标准画法。
##     4. Panel_Gene_Coverage_By_Series.pdf
##        —— 形式：单面板热图（x = 系列，y = 12 基因面板，tile 着色 present/absent）
##        —— 回答的科学问题：预注册面板（§10）在每个系列的可及性——哪些系列缺哪个基因，
##           为层 1（S05）提供「Q4 缺失逐条记录」的可视凭据。
##        —— 参照文献：基因覆盖矩阵热图是跨队列面板可行性核验的常规呈现（同基因集
##           在多个平台上的可及性对照）。
##
## 图件口径：sci_theme / sci_save 逐字取自规范 §sci_theme / §sci_save 官方定义；
##       不引入任何拼图包，单面板 ggplot2 出图；可读性核验见 sci_save 之 RAPA_QC_PNG 钩子。
##
## 前置：S01 通过。
##
## CEL 路径已废止（2026-09-28 用户裁定方案 A，理由须留档）：
##   S03 原登记了 7 个 source="cel" 的 Affymetrix 系列（GSE43292 / GSE21545 /
##   GSE28829 / GSE100927 / GSE124020 / GSE124106 / GSE212388），原本依赖
##   pd.hgfocus / pd.hgu133a2 / pd.hugene10kplus 三个平台包才能做 RMA。
##   实测这三者在 Bioconductor 3.23 / 3.19 / 3.14 release 页面全部 404，
##   旧版 src/contrib 的 tarball 同样 404 —— **上游已彻底移除，本机无法补装**。
##   故这 7 个系列改走 GEO 的 series_matrix（作者已处理好的表达表），
##   不再触碰 affy / oligo / pd.*。
##
## 尺度口径（重要，勿回退）：
##   series_matrix 里的数值究竟是否已经 log2，文件头**没有任何官方声明**，
##   无法从元数据判断。本脚本只登记**可观测事实**：q001 / qmed / q999 与
##   `has_negative`（0.1% 分位 < 0，即矩阵中存在负值）。
##   `has_negative = TRUE` 是「必然为 log 尺度」的**充分证据**（原始强度不可能为负）；
##   但 `has_negative = FALSE` **不能反推「不是 log2」**——已 log2 的数据若下限为正
##   （如 GPL570 系列实测量程 1.7–14.4），同样为 FALSE。
##   故本列不得用作「尺度判定结果」，只作负值存在的证据；图 2 改按**已登记的
##   归一化路径**（norm）着色，避免把推断当结论。逐系列尺度须在方案 §3 逐条登记。
##
## 表达量来源口径：family SOFT **只由 S01 消费**（解析样本级注释、不提供表达表），
##       表达一律来自各系列的 matrix / counts 文件。故本脚本不含 SOFT 读取分支。
##
## 与方案 §7 的一处偏差（须在报告中披露）：方案写产物后缀 `.qs`；本脚本用 base
##   `saveRDS()` 产出 `.rds`。理由：qs 不在本项目栈内，且规则 34 要求脚本自足；
##   两者都是确定性序列化，语义等价。
#===============================================================

## PART 0 常量与注册表 ---------------------------------------------------------------

## 原始数据根目录（与 S01_ingest.R / S12a_prepare_spatial.R 的 dir.data 同名同义；
## object 'dir.data' not found）
dir.data <- "01_data"

dir.tab <- "03_results/S03/table"
## 规则 33：本步重时（全基因归一化），断点续传日志落 checkpoints/ 子目录
CKPT_DIR <- file.path(dir.tab, "checkpoints")
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(dir.tab, showWarnings = FALSE, recursive = TRUE)

FIG_DIR <- "03_results/S03/Figure"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

library(ggplot2)  ## 归一化 QC 图；主题与出口用规范 §sci_theme / §sci_save
library(limma)  ## voom 归一化 + 基因级矩阵的符号归并（avereps）
library(AnnotationDbi)  ## 键值映射接口（GSE323983 的 ENSEMBL -> SYMBOL）
library(org.Hs.eg.db)  ## 人类基因注释库（官方离线源；唯一 ENSEMBL 符号来源）

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
## 可读性核验钩子（规则 22 三级审核 · 交付前须核验可读性）：仅当显式设置环境变量
## RAPA_QC_PNG=1 时，额外导出一份 PNG 预览到 FIG_DIR/_qc/ 供人工/AI 核验；
## 该 PNG 不进入交付物、默认不产出，保持交付物为纯 PDF（符合规范 §sci_save）。
sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  pdf(pdf_file, width = width, height = height)
  print(plot_obj)
  dev.off()
  if (Sys.getenv("RAPA_QC_PNG") == "1") {
    dir.create(file.path(FIG_DIR, "_qc"), showWarnings = FALSE)
    png_file <- file.path(FIG_DIR, "_qc", paste0(filename, ".png"))
    png(png_file, width = width * 150, height = height * 150, res = 150)
    print(plot_obj)
    dev.off()
  }
}

## Q2/Q3：表达矩阵至少要有这么多可用特征（探针级与基因级同口径）
MIN_FEATURE_N <- 1000L
MAX_MISSING_RATE <- 0.5       ## Q2：单探针缺失率上限
## counts 表里「不是样本」的列名（按名剔除，不靠类型推断 —— 类型推断会把
## 本应是数值却被读成字符的样本列静默丢掉）。来源：
##   "Detection Pval" —— GSE57691 的 non-normalized_data.txt.gz 表头实测
##   "gene_name"      —— GSE323983 的 RapavIRRapa counts 表头实测
##                       （`gene_id | gene_name | <12 样本列>`）；不剔则 as.matrix
##                       整表转字符，下游 quantile 报 non-numeric argument。
COUNT_NON_SAMPLE <- c("Detection Pval", "gene_name")

TC6002_ACC <- "GSE6002"
TC6002_FILE_RX <- "^GSM[0-9]+_Cy[35]\\.txt\\.gz$"
TC6002_N_FILE <- 8L       ## 4 donor x (Cy3 + Cy5)

## 12 基因固定面板，逐字取自方案 §10 预注册锁定与 §9 Q4（禁止同义替换）。
PANEL_GENES <- c("IL6", "IL1B", "TNF", "CCL2", "CXCL8", "MMP9",
  "SPP1", "NLRP3", "IL10", "TGFB1", "MERTK", "GAS6")

## 产物路径。断点续传判据：一个系列必须「探针级 + 基因级两个产物都已落盘」才算完成。
ExprFile  <- function(acc) file.path(dir.tab, paste0(acc, "_expr.rds"))
GeneFile  <- function(acc) file.path(dir.tab, paste0(acc, "_gene.rds"))
DoneSeries <- function(acc) file.exists(ExprFile(acc)) && file.exists(GeneFile(
    acc))

## 每系列的取数与归一化路径。note 记录必须遵守的硬约束（§3.2.1 等）。
dat.series <- rbind(
  data.frame(acc = "GSE232911", source = "matrix_log2", norm = "identity",
    note = paste0("已 log2 且已归一化；禁止二次 log2（量程以",
      "normalization_log.csv 实测为准）")),
  data.frame(acc = "GSE57691", source = "counts_column", norm = "voom",
    note = "用 non-normalized_data.txt.gz；整数 counts 走 voom"),
  data.frame(acc = "GSE43292", source = "series_matrix", norm = "identity",
    note = "原 GPL6244 + pd.hgfocus；平台包已消失，改走 GEO series_matrix"),
  data.frame(acc = "GSE21545", source = "series_matrix", norm = "identity",
    note = "原 GPL570 + pd.hgu133a2；平台包已消失，改走 GEO series_matrix"),
  data.frame(acc = "GSE28829", source = "series_matrix", norm = "identity",
    note = "原 GPL570 + pd.hgu133a2；平台包已消失，改走 GEO series_matrix"),
  data.frame(acc = "GSE100927", source = "series_matrix", norm = "identity",
    note = "原 GPL17077 本就无 pd 包；改走 GEO series_matrix"),
  data.frame(acc = "GSE47472", source = "series_matrix", norm = "identity",
    note = paste0("series_matrix 已是表达表；RAW.tar 内是 GPL10558",
      "注释，禁止当矩阵读")),
  data.frame(acc = "GSE41571", source = "series_matrix", norm = "identity",
    note = ""),
  data.frame(acc = "GSE323983", source = "counts_block", norm = "voom",
    note = "两个 counts 文件为块结构，块诊断见 S04"),
  ## GSE65145 已从归一化表移除（D' 处置，2026-09-29 用户裁定）
  # 实测：RAW.tar 内 6 文件是 SeqMonk bedGraph RPKM（4 列 chr/start/end/RPKM，149123 区间），非
  # counts；
  ## 作者定量已可用，不重定量（方案 A：6.2 GB SRA 重比对代价不成比例）；
  ## 改作 S11 静态反向警示签名（01_data/verify/gse65145_rapamycin_up_signature.csv），不进主聚合。
  data.frame(acc = "GSE6002", source = "two_color", norm = "loess",
    note = "双色 limma；read.maimages -> normalizeWithinArrays"),
  data.frame(acc = "GSE124020", source = "series_matrix", norm = "identity",
    note = paste0("原 GPL17586 + pd.hugene10kplus；平台包已消失，改走 GEO",
      "series_matrix")),
  data.frame(acc = "GSE124106", source = "series_matrix", norm = "identity",
    note = "原 GPL17586 + pd.hugene10kplus；仅取前 12 样本")
)

## GSE212388 已从 bulk 归一化表移除（2026-09-28 用户裁定）
## 实测依据：RAW.tar 内 4 个样本为 10X mtx 三件套（GSM6529255/257/260/261）、
## 4 个样本为 gene×cell 宽表（GSM6529254/256/258/259），全部为单细胞计数矩阵；
## 且其 series_matrix 无表达块（0 探针）。单细胞归单细胞路径，由 S12 承接。

file.norm.log <- file.path(CKPT_DIR, "norm_log.csv")
file.norm.qc <- file.path(dir.tab, "normalization_qc.csv")
file.gene.log <- file.path(CKPT_DIR, "gene_log.csv")
file.panel.cov <- file.path(dir.tab, "panel_coverage.csv")

## 本轮待处理系列 = 尚未同时具备探针级与基因级产物的系列。
## 旧写法只看 norm_log 中的 acc，无法驱动「补出基因级矩阵」这一步，故改为按产物存在性判定。
vec.todo <- dat.series$acc[!vapply(dat.series$acc, DoneSeries, logical(1L))]

## PART 1 各取数源的读取函数 ----------------------------------------------------------

## 注意：family SOFT 只由 S01 消费（解析样本级注释），**不提供表达表**；
## 表达量一律来自各系列的 matrix / counts 文件，故此处无 SOFT 读取函数。

ReadSeriesMatrix <- function(path) {
  ## 读取 GEO 表达矩阵文件，统一处理两种列结构：
  ##  (a) 标准 series_matrix：以 !Series_* / !sample_* / !series_matrix_table_begin
  ##      等 `!` 注释行开头；真正表头是 ID_REF / GSMxxx。第 1 列 = ID_REF（可为文本
  ##      也可为数值，如 GPL6244 的探针号），其余全为数值样本列。须用
  ##      comment.char="!" 跳过所有 `!` 注释行（否则首行 !Series_title 会被当表头，
  ##      读出 2 列垃圾）。
  ##  (b) GSE232911 的 Matrix.file：无 `!` 行，前 8 列为注释（含 ID / Public.Gene.IDs /
  ##      Gene.Symbol / Description），第 9 列起才是 246 个数值样本列，且前 4 列还有
  ##      冗余数值列（索引 / 占位）。旧写法 as.matrix(dat.raw[, -1L]) 会把文本注释列
  ##      一并纳入，整表转字符矩阵，quantile 报 "non-numeric argument"。
  ## 统一规则：表达列 = 「以末列为终点、向前延伸的最长连续数值列块」并剔除 probe 列；
  ##          probe  = 优先 Gene.Symbol（GSE232911 自带），否则 ID_REF/ID，
  ##                   否则首个非数值列，否则第 1 列。
  dat.raw <- read.delim(path, header = TRUE, stringsAsFactors = FALSE,
    check.names = FALSE, fill = TRUE, comment.char = "!")
  is.num <- vapply(dat.raw, is.numeric, logical(1L))
  n.col <- length(is.num)
  ## 末列向前找最长连续数值块起点 k（跳过尾随的数值注释列）
  k <- n.col
  while (k >= 1L && is.num[k]) k <- k - 1L
  k <- k + 1L
  expr.cols <- k:n.col
  ## probe 列选择（独立于表达块）
  probe.col <- which(colnames(dat.raw) == "Gene.Symbol")
  if (length(probe.col) == 0L) probe.col <- which(colnames(dat.raw) %in% c(
      "ID_REF",
      "ID"))
  if (length(probe.col) == 0L) probe.col <- which(!is.num)[1L]
  if (length(probe.col) == 0L || is.na(probe.col)) probe.col <- 1L
  ## 剔除 probe 列（当 ID_REF 本身为数值列时，避免被当成样本列）
  expr.cols <- setdiff(expr.cols, probe.col)
  if (length(expr.cols) == 0L) {
    stop(paste0("ReadSeriesMatrix: 在 ", path, " 中找不到可用的数值表达列"))
  }
  probe <- as.character(dat.raw[[probe.col]])
  expr <- as.matrix(dat.raw[, expr.cols, drop = FALSE])
  storage.mode(expr) <- "double"
  list(probe = probe, expr = expr)
}

## 整数 counts 矩阵（列 = 样本）
## 注意：GSE57691 的 counts 表把样本列与 "Detection Pval" 列**交替**排放
## （`ID_REF | Sample 1 | Detection Pval | Sample 2 | …`），故不能简单地丢第 1 列——
## 那样会把 68 个 p 值列当成样本，ncol 由 68 翻成 136（违反 §9 Q1）。
ReadCounts <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  dat.raw <- read.delim(con, header = TRUE, stringsAsFactors = FALSE,
    check.names = FALSE)
  close(con)
  rownames(dat.raw) <- as.character(dat.raw[[1L]])
  vec.keep <- setdiff(colnames(dat.raw), c(colnames(dat.raw)[1L],
      COUNT_NON_SAMPLE))
  list(probe = rownames(dat.raw),
    expr = as.matrix(dat.raw[, vec.keep, drop = FALSE]))
}

## PART 1.5 平台注释层：探针 ID -> 基因符号 -------------------------------------------------
##
## 为什么必须有这一层：S03 原先把行名改成 p1..pN 后就把探针 ID 丢了，而 11 个系列里
## 除 GSE232911（矩阵自带 Gene.Symbol）外，$probe 存的全是平台探针 ID，12 基因面板
## 命中 0/12（实测见 01_data/verify/check_s03_probe_contract.R）。
##
## 平台号来源（2026-09-29 实测，01_data/verify/check_platform_soft.R）：
##   各 series_matrix 的 !Series_platform_id 与各 family.soft.gz 的 !Series_platform_id。
## 注释文件格式来源（2026-09-29 实测，01_data/verify/check_annotation_sources.R 与
##   check_symbol_extraction.R）：
##   - GEO Annotation（GPL570/6244/10558/GPL17586）：若干 `#`/`!`/`^` 注释行 +
##     1 行制表符表头 + 数据。GPL17586-45144.txt 的 gene_assignment 段落结构为
##     每个基因 5 段 `transcript // symbol // description // cytoband // entrez`，
##     多基因之间用 ` /// ` 连接。
##   - GPL17077.annot.csv：带引号的标准 CSV，符号列 GENE_SYMBOL。
##   - GSE323983：GPL28975（GridION 纳米孔测序），无平台注释表，行名是 ENSEMBL ID，
##     唯一可用官方离线源是 org.Hs.eg.db。
##   - GSE232911：矩阵自带 Gene.Symbol 列，无需外部注释。
dat.annot <- rbind(
  data.frame(acc = "GSE21545",  platform = "GPL570",   reader = "geo_annot",
    annot_file = "GPL570.annot.gz",   sym_col = "Gene symbol", sep = "///",
    sample_map = NA_character_, probe_col = "ID", stringsAsFactors = FALSE),
  data.frame(acc = "GSE28829",  platform = "GPL570",   reader = "geo_annot",
    annot_file = "GPL570.annot.gz",   sym_col = "Gene symbol", sep = "///",
    sample_map = NA_character_, probe_col = "ID", stringsAsFactors = FALSE),
  data.frame(acc = "GSE41571",  platform = "GPL570",   reader = "geo_annot",
    annot_file = "GPL570.annot.gz",   sym_col = "Gene symbol", sep = "///",
    sample_map = NA_character_, probe_col = "ID", stringsAsFactors = FALSE),
  data.frame(acc = "GSE43292",  platform = "GPL6244",  reader = "geo_annot",
    annot_file = "GPL6244.annot.gz",  sym_col = "Gene symbol", sep = "///",
    sample_map = NA_character_, probe_col = "ID", stringsAsFactors = FALSE),
  data.frame(acc = "GSE47472",  platform = "GPL10558", reader = "geo_annot",
    annot_file = "GPL10558.annot.gz", sym_col = "Gene symbol", sep = "///",
    sample_map = NA_character_, probe_col = "ID", stringsAsFactors = FALSE),
  data.frame(acc = "GSE57691",  platform = "GPL10558", reader = "geo_annot",
    annot_file = "GPL10558.annot.gz", sym_col = "Gene symbol", sep = "///",
    sample_map = "source_col", probe_col = "ID", stringsAsFactors = FALSE),
  data.frame(acc = "GSE100927", platform = "GPL17077", reader = "csv_annot",
    annot_file = "GPL17077.annot.csv", sym_col = "GENE_SYMBOL", sep = "///",
    sample_map = NA_character_, probe_col = "ID", stringsAsFactors = FALSE),
  data.frame(acc = "GSE124020", platform = "GPL17586",
    reader = "gene_assignment",
    annot_file = "GPL17586-45144.txt", sym_col = "gene_assignment",
    sep = "///", sample_map = NA_character_, probe_col = "ID",
    stringsAsFactors = FALSE),
  data.frame(acc = "GSE124106", platform = "GPL17586",
    reader = "gene_assignment",
    annot_file = "GPL17586-45144.txt", sym_col = "gene_assignment",
    sep = "///", sample_map = NA_character_, probe_col = "ID",
    stringsAsFactors = FALSE),
  data.frame(acc = "GSE232911", platform = "GPL17586", reader = "self",
    annot_file = NA_character_, sym_col = "Gene.Symbol", sep = ";",
    sample_map = "gse232911_title", probe_col = "ID",
    stringsAsFactors = FALSE),
  data.frame(acc = "GSE323983", platform = "GPL28975", reader = "ensembl",
    annot_file = NA_character_, sym_col = "ENSEMBL", sep = "///",
    sample_map = "source_col", probe_col = "ID", stringsAsFactors = FALSE),
  data.frame(acc = "GSE6002", platform = "GPL3784", reader = "geo_annot",
    annot_file = "GPL3784_family.soft.gz", sym_col = "ORF", sep = "///",
    sample_map = NA_character_, probe_col = "miltenyi_ID",
    stringsAsFactors = FALSE)
)

## 注释表缓存：GPL17586-45144.txt 为 186 MB / 70753 行，被两个系列共用，避免重复解析。
env.annot <- new.env(parent = emptyenv())

## GEO Annotation 读取器：定位「首个不以 # / ! / ^ 开头的行」即为表头行。
ReadGeoAnnot <- function(a, path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  ln <- readLines(con)
  close(con)
  i.hdr <- which(!grepl("^[#!^]", ln))[1L]
  i.end <- which(grepl("^!", ln))
  i.end <- i.end[i.end > i.hdr][1L]
  if (is.na(i.end)) i.end <- length(ln)
  dat <- read.delim(text = paste(ln[i.hdr:(i.end - 1L)], collapse = "\n"),
    header = TRUE, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE)
  pc <- if (is.null(a$probe_col)) "ID" else a$probe_col
  data.frame(probe = as.character(dat[[pc]]),
    sym = as.character(dat[[a$sym_col]]), stringsAsFactors = FALSE)
}

## GPL17077 的注释是带引号的 CSV（LOCUSLINK_ID 列含裸 NA）。
ReadCsvAnnot <- function(path, sym_col) {
  dat <- read.csv(path, quote = "\"", comment.char = "", na.strings = "NA",
    check.names = FALSE, stringsAsFactors = FALSE)
  data.frame(probe = as.character(dat[["ID"]]),
    sym = as.character(dat[[sym_col]]), stringsAsFactors = FALSE)
}

## GPL17586 的 gene_assignment：逐基因 5 段，取第 2 段为符号，多基因以 `///` 拼回。
SplitGeneAssignment <- function(x) {
  vapply(strsplit(x, " /// ", fixed = TRUE), function(one) {
      seg <- strsplit(one, " // ", fixed = TRUE)
    paste(vapply(seg, function(s) if (length(s) >= 2L) s[2L] else NA_character_,
          character(1L)), collapse = "///")
    }, character(1L))
}

## 取（并缓存）某平台的注释表；返回 data.frame(probe, sym)。
AnnotTable <- function(a) {
  key <- a$annot_file
  if (!exists(key, envir = env.annot, inherits = FALSE)) {
    fp <- file.path("01_data", key)
    if (!file.exists(fp)) stop(paste0("S03: 缺少平台注释文件 ", fp))
    assign(key, if (a$reader == "csv_annot") {
        ReadCsvAnnot(fp, a$sym_col)
      } else {
        ReadGeoAnnot(a, fp)
      }, envir = env.annot)
  }
  get(key, envir = env.annot, inherits = FALSE)
}

## 探针 ID 向量 -> 与探针等长的符号向量（多符号以 "///" 拼接；查无符号为 NA）。
SymbolFromProbe <- function(a, probe) {
  if (a$reader == "self") return(probe)
  if (a$reader == "ensembl") {
    ## GPL28975 无平台注释表；org.Hs.eg.db 是唯一可用的官方离线映射源。
    tbl <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(probe),
      keytype = "ENSEMBL", columns = "SYMBOL")
    return(tbl$SYMBOL[match(probe, tbl$ENSEMBL)])
  }
  tab <- AnnotTable(a)
  if (a$reader == "gene_assignment") {
    return(SplitGeneAssignment(tab$sym[match(probe, tab$probe)]))
  }
  tab$sym[match(probe, tab$probe)]
}

## 矩阵列名 -> 真实样本 ID。series_matrix 表头本就是 GSM 号，无需映射；
## GSE232911 的矩阵列是样本 title，须由 S01 的分组表映射（§9 Q14 以该表为准）。
SampleIds <- function(a, vec.raw) {
  if (is.na(a$sample_map)) return(vec.raw)
  if (a$sample_map == "gse232911_title") {
    fp <- "03_results/S01/table/gse232911_grouping.tsv"
    if (!file.exists(fp)) {
      stop(paste0("S03: GSE232911 需要 ", fp,
          " 才能把矩阵列标题映射回 GSM，但它不存在。"))
    }
    dat.g <- read.delim(fp, stringsAsFactors = FALSE)
    i <- match(vec.raw, dat.g$title)
    if (anyNA(i)) {
      stop(paste0("S03: GSE232911 有 ", sum(is.na(i)),
          paste0("个矩阵列无法由 title 映射回 sample_id（Q14",
            "以该表为准，不得放行）。")))
    }
    return(dat.g$sample_id[i])
  }
  if (a$sample_map == "source_col") {
    ##   列名是来源文件自带的匿名列名（"Sample N" / "<arm><donor>_sum"），与
    ##   design_keys.tsv 的 sample_id（= GSM）无法按名连接 —— 实测交集 = 0
    ##   （证据 01_data/verify/check_s05_sample_axis_gap.R）。改由 S01 的
    ##   sample_axis.tsv 重映射回真实 GSM。
    fp <- "03_results/S01/table/sample_axis.tsv"
    if (!file.exists(fp)) {
      stop(paste0("S03: ", a$acc, " 需要 ", fp,
          " 才能把来源列名映射回 GSM，但它不存在。",
          "请先跑 S01。"))
    }
    dat.ax <- read.delim(fp, stringsAsFactors = FALSE)
    dat.ax <- dat.ax[dat.ax$acc == a$acc, , drop = FALSE]
    if (nrow(dat.ax) == 0L) {
      stop(paste0("S03: sample_axis.tsv 中缺 ", a$acc, " 的行。"))
    }
    i <- match(vec.raw, dat.ax$source_col)
    if (anyNA(i)) {
      stop(paste0("S03: ", a$acc, " 有 ", sum(is.na(i)),
          " 个矩阵列无法由 sample_axis.tsv 映射回 GSM（不得放行）。未命中: ",
          paste(head(vec.raw[is.na(i)], 5L), collapse = ", ")))
    }
    return(dat.ax$sample_id[i])
  }
  stop(paste0("S03: 未知 sample_map '", a$sample_map, "'"))
}

## 判别某探针的符号串里是否至少有一个可用符号（""、NA、"---" 皆为「无符号」）。
## 必要性：HTA-2.0 的 gene_assignment 对未注释探针写 "--- // --- // …"，字符串非空，
## 若只判非空会把它们误计为「已映射」（实测 GSE124020 会由 40946 虚报成 67528）。
HasUsableSymbol <- function(vec.sym, sep) {
  vapply(strsplit(vec.sym, sep, fixed = TRUE), function(s) {
      s <- trimws(s)
      any(!is.na(s) & !s %in% c("", "---"))
    }, logical(1L))
}

## 基因级矩阵：探针级矩阵 + 探针 -> 符号向量 -> 展开 -> 按符号取均值。
## limma::avereps 的官方实现为 rowsum(x, ID, na.rm = TRUE) / 非 NA 计数
## （见 01_data/verify/check_avereps_signature.R），故「部分探针缺失」不会把整符号
## 拉成 NA；全 NA 符号会得到 NaN，逐一剔除。
BuildGeneMatrix <- function(expr, vec.sym, sep) {
  lst.sym <- strsplit(vec.sym, sep, fixed = TRUE)
  vec.sym.flat <- trimws(unlist(lst.sym, use.names = FALSE))
  vec.row <- rep(seq_along(lst.sym), lengths(lst.sym))
  lgl.ok <- !is.na(vec.sym.flat) & !vec.sym.flat %in% c("", "---")
  if (!any(lgl.ok)) stop("S03: 该系列没有任何探针能映射到基因符号")
  expr.gene <- limma::avereps(expr[vec.row[lgl.ok], , drop = FALSE],
    ID = vec.sym.flat[lgl.ok])
  expr.gene[!apply(expr.gene, 1L, function(x) all(is.na(x))), , drop = FALSE]
}

ReadGpr <- function(path) {
  con <- gzfile(path, "rt")
  ln <- readLines(con)
  close(con)
  i.brd <- which(trimws(ln) == "Begin Raw Data")[1L]
  if (is.na(i.brd)) {
    stop(paste0("S03: GPR 文件缺少 'Begin Raw Data' 标记：", path))
  }
  dat <- read.delim(text = paste(ln[(i.brd + 1L):length(ln)],
      collapse = "\n"), sep = "\t", check.names = FALSE,
    stringsAsFactors = FALSE)
  vec.miss <- setdiff(c("Gene ID", "Signal Mean", "Background Mean"),
    colnames(dat))
  if (length(vec.miss) > 0L) {
    stop(paste0("S03: GPR 文件缺列 ", paste(vec.miss, collapse = ", "),
        "：", path))
  }
  list(probe = as.character(dat[["Gene ID"]]),
    sig = as.numeric(dat[["Signal Mean"]]),
    bg = as.numeric(dat[["Background Mean"]]))
}

## PART 2 逐系列归一化 -------------------------------------------------------------------

## 断点续传要点：dat.norm.log / dat.norm.qc 必须从已落盘的
## checkpoints 续读，不能在每轮重跑时以空 data.frame 起步。否则上一轮已处理的系列记录
## 会被本轮回写覆盖丢失。本轮另需按 vec.todo 丢弃「即将重跑系列」的旧行，避免重复登记。
if (file.exists(file.norm.log)) {
  dat.norm.log <- read.csv(file.norm.log, stringsAsFactors = FALSE)
  dat.norm.log <- dat.norm.log[!dat.norm.log$acc %in% vec.todo, , drop = FALSE]
} else {
  dat.norm.log <- data.frame(acc = character(0), source = character(0),
    norm = character(0),
    n_probe = integer(0), n_sample = integer(0),
    range_min = numeric(0), range_max = numeric(0),
    has_negative = logical(0), q001 = numeric(0),
    qmed = numeric(0), q999 = numeric(0),
    cdf_platform = character(0), note = character(0),
    stringsAsFactors = FALSE)
}
## dat.norm.qc 的四列结构与 PART 2 2.4 的 rbind 对齐；
## 旧写法 dat.norm.log[0, c("acc","qc",...)] 会因 dat.norm.log 无 qc/observed/expected
## 列而报 "undefined columns selected"（2026-09-28 首跑实测）
if (file.exists(file.norm.qc)) {
  dat.norm.qc <- read.csv(file.norm.qc, stringsAsFactors = FALSE)
  ## 丢弃「本轮将重跑的系列」的旧 qc 行（含历史 FAIL），避免 stale 失败行堆积
  dat.norm.qc <- dat.norm.qc[!dat.norm.qc$acc %in% vec.todo, , drop = FALSE]
} else {
  dat.norm.qc <- data.frame(acc = character(0), qc = character(0),
    observed = character(0), expected = character(0),
    stringsAsFactors = FALSE)
}
if (file.exists(file.gene.log)) {
  dat.gene.log <- read.csv(file.gene.log, stringsAsFactors = FALSE)
  dat.gene.log <- dat.gene.log[!dat.gene.log$acc %in% vec.todo, , drop = FALSE]
} else {
  dat.gene.log <- data.frame(acc = character(0), platform = character(0),
    n_probe = integer(0), n_probe_mapped = integer(0),
    n_symbol = integer(0), n_multi_symbol_probe = integer(0),
    stringsAsFactors = FALSE)
}
if (file.exists(file.panel.cov)) {
  dat.panel.cov <- read.csv(file.panel.cov, stringsAsFactors = FALSE)
  dat.panel.cov <- dat.panel.cov[!dat.panel.cov$acc %in% vec.todo, ,
    drop = FALSE]
} else {
  dat.panel.cov <- data.frame(acc = character(0), gene = character(0),
    present = character(0), stringsAsFactors = FALSE)
}

## for 循环必须按**行**迭代：R 的 for(data.frame) 默认按列，
## s$acc 会因 s 是列向量而报 "$ operator is invalid for atomic vectors"
for (i.row in seq_len(nrow(dat.series))) {
  s <- dat.series[i.row, , drop = FALSE]
  if (!s$acc %in% vec.todo) next
  acc <- s$acc
  cat(sprintf("\n[S03] %s  (source=%s, norm=%s)\n", acc, s$source, s$norm))

  ## 逐系列容错：任一系列取数/归一化失败只记 FAIL 并继续，不中断整轮
  ## （规则 33 断点续传 + 规则 26 显式报错而非静默跳过：末尾以非 0 退出汇总 FAIL）
  res <- tryCatch({
      ## 2.1 按 source 取数
      if (s$source == "matrix_log2") {
        dat.raw <- ReadSeriesMatrix(file.path("01_data",
            "GSE232911_Matrix.file.txt.gz"))
      } else if (s$source == "series_matrix") {
        ## 通用落点：<acc>_series_matrix.txt.gz。原先这里是 if/else 硬编码两个系列名，
        ## 无法容纳方案 A 新增的 7 个系列，故改为按 acc 拼名。
        f <- paste0(acc, "_series_matrix.txt.gz")
        fp.sm <- file.path("01_data", f)
        if (!file.exists(fp.sm)) {
      stop(paste0("S03: 缺少 ", f, "。这 7 个系列原走 CEL+pd.*，因平台包已从 ",
              "Bioconductor 移除而改走 GEO series_matrix，请先下载该文件。"))
        }
        dat.raw <- ReadSeriesMatrix(fp.sm)
      } else if (s$source == "counts_column") {
        dat.raw <- ReadCounts(file.path("01_data",
            "GSE57691_non-normalized_data.txt.gz"))
      } else if (s$source == "counts_block") {
        f <- if (acc == "GSE323983") {
          c("GSE323983_CvIR_all_gene_counts.tsv.gz",
            "GSE323983_RapavIRRapa_all_gene_counts.tsv.gz")
        } else {
          file.path("01_data", acc, "counts.txt.gz")
        }
        list.blk <- lapply(file.path("01_data", f), ReadCounts)
        mat.b1 <- list.blk[[1L]]$expr
        mat.b2 <- list.blk[[2L]]$expr
        vec.id <- intersect(rownames(mat.b1), rownames(mat.b2))
        ## match 取首次出现，保证 probe 向量与 expr 行轴一一对应；重复 ID 数登记进 note
        dat.raw <- list(
          probe = vec.id,
          expr = cbind(mat.b1[match(vec.id, rownames(mat.b1)), , drop = FALSE],
            mat.b2[match(vec.id, rownames(mat.b2)), , drop = FALSE]))
        s$note <- paste0(s$note, "；实测块合并：", ncol(mat.b1), " + ", ncol(
            mat.b2),
          " = ", ncol(dat.raw$expr), " 样本，ENSEMBL 交集 ",
          length(vec.id), "（dup: ",
          sum(duplicated(rownames(mat.b1))), "+",
          sum(duplicated(rownames(mat.b2))), "）")
      } else if (s$source == "cel") {
        stop(paste0("S03: source 'cel' for ", acc,
            " requires an unimplemented download/extraction branch. ",
            paste0("This branch was written but NOT exercised yet; run it",
              "only after the "),
            "corresponding CEL inputs are confirmed on disk."))
      } else if (s$source == "two_color") {
        ##   GPR 列名非标准（Signal Mean / Background Mean / Gene ID），
        ##   read.maimames 的 genepix 列名匹配不适用 -> 手工构建 RGList 后走
        ##   limma 官方双色流程：backgroundCorrect("subtract") ->
        ##   normalizeWithinArrays("loess")。探针轴按名对齐（铁律 2），禁依赖位置。
        vec.f.tc <- sort(list.files(dir.data, pattern = TC6002_FILE_RX))
        if (length(vec.f.tc) != TC6002_N_FILE) {
          stop(paste0("S03: ", acc, " 需要 ", TC6002_N_FILE,
              " 个 Cy3/Cy5 文件在 ", dir.data, "/，实测 ",
              length(vec.f.tc), " 个（模式 ", TC6002_FILE_RX, "）。"))
        }
        vec.gsm <- sub("_Cy[35]\\.txt\\.gz$", "", vec.f.tc)
        ## 通道提取必须锚定 .txt.gz 后缀：裸前缀替换会把 "Cy3.txt.gz" 当通道名，
        ## which(vec.ch == "Cy5") 落空 -> cbind 空列表得 NULL -> dimnames(NULL) 报
        ## "attempt to set an attribute on NULL"（2026-09-29 首跑实测崩溃点）
        vec.ch <- sub("^GSM[0-9]+_(Cy[35])\\.txt\\.gz$", "\\1", vec.f.tc)
        if (anyDuplicated(paste0(vec.gsm, vec.ch))) {
          stop(paste0("S03: ", acc, " 的双色文件名存在重复 GSM/通道组合。"))
        }
        ## 逐文件解析 GPR：探针 = "Gene ID"，前景/背景 = Signal/Background Mean。
        lst.ch <- lapply(vec.f.tc, function(f) {
            x <- ReadGpr(file.path(dir.data, f))
            list(gsm = sub("_Cy.*$", "", f), ch = sub("^.*_Cy", "Cy", f),
              probe = x$probe, sig = x$sig, bg = x$bg)
          })
        ## 探针轴按名对齐：以首文件为参照，逐一匹配（禁位置对齐）。
        v.ref <- lst.ch[[1L]]$probe
        lst.ch <- lapply(lst.ch, function(x) {
            i <- match(v.ref, x$probe)
            if (anyNA(i)) {
              stop(paste0("S03: ", acc, " 文件 ", x$gsm, "_", x$ch,
                  " 有 ", sum(is.na(i)), " 个探针不在参照探针轴上。"))
            }
            list(gsm = x$gsm, ch = x$ch, sig = x$sig[i], bg = x$bg[i])
          })
        i.cy5 <- which(vec.ch == "Cy5")
        i.cy3 <- which(vec.ch == "Cy3")
        if (!setequal(vec.gsm[i.cy5], vec.gsm[i.cy3]) ||
            length(i.cy5) != length(i.cy3)) {
          stop(paste0("S03: ", acc, " 的 Cy3/Cy5 文件不成对（各 4 donor）。"))
        }
        ## 染料方向（SOFT 官方字段）：Cy5 = LPS+Rapa（处理），Cy3 = LPS（对照）。
        i5 <- i.cy5[order(vec.gsm[i.cy5])]
        i3 <- i.cy3[order(vec.gsm[i.cy3])]
        mat.R <- do.call(cbind, lapply(lst.ch[i5], function(x) x$sig))
        mat.G <- do.call(cbind, lapply(lst.ch[i3], function(x) x$sig))
        mat.Rb <- do.call(cbind, lapply(lst.ch[i5], function(x) x$bg))
        mat.Gb <- do.call(cbind, lapply(lst.ch[i3], function(x) x$bg))
        vec.gsm.ord <- sort(unique(vec.gsm))
        dimnames(mat.R) <- list(v.ref, vec.gsm.ord)
        dimnames(mat.G) <- list(v.ref, vec.gsm.ord)
        dimnames(mat.Rb) <- list(v.ref, vec.gsm.ord)
        dimnames(mat.Gb) <- list(v.ref, vec.gsm.ord)
        ## limma 不导出 RGList() 构造函数（2026-09-29 实测 getNamespaceExports
        ## 无 RGList/RG.List/MAList；RGList 为已注册 S4 类，isClass=TRUE），
        ## 官方构造方式 = User's Guide 双色章标准写法 new("RGList", list(...))
        obj.rg <- methods::new("RGList", list(R = mat.R, G = mat.G,
          Rb = mat.Rb, Gb = mat.Gb,
          genes = data.frame(ID = v.ref, stringsAsFactors = FALSE)))
        obj.bc <- limma::backgroundCorrect(obj.rg, method = "subtract")
        obj.ma <- limma::normalizeWithinArrays(obj.bc, method = "loess")
        dat.raw <- list(expr = obj.ma$M, probe = as.character(obj.ma$genes$ID))
        cat(sprintf("[S03] %s  two-colour loess: %d probes x %d arrays\n",
            acc, nrow(obj.ma$M), ncol(obj.ma$M)))
      } else {
        stop(paste0("S03: unknown source '", s$source, "' for ", acc))
      }

      ## 2.15 样本身份：矩阵列名（series_matrix 表头即 GSM；其余用 counts 自带列名），
      ##      GSE232911 的 title 须经 S01 分组表映射回 GSM。此向量用于基因级矩阵列名。
      ##      注释注册表放在取数之后查：未实现双色分支的 GSE6002 应先报它自己的
  ##      「two_color 未实现」而不是被注释表缺行抢先（旧顺序会产生误导性报错）。
      a.annot <- dat.annot[dat.annot$acc == acc, , drop = FALSE]
      if (nrow(a.annot) != 1L) {
        stop(paste0("S03: 注释注册表中 ", acc, " 应恰好 1 行，实得 ",
            nrow(a.annot)))
      }
      vec.sample <- SampleIds(a.annot, colnames(dat.raw$expr))

      ## 2.2 尺度判定 + 归一化
      ## series_matrix 是否已是 log2，文件头没有官方声明，只能从分布推断（见头部口径）。
  ## 判据取"已 log2 的值可以为负"这一事实：若 0.1% 分位 < 0，判定为 log 尺度。
      vec.qtile <- stats::quantile(dat.raw$expr, c(0.001, 0.5, 0.999),
        na.rm = TRUE, names = FALSE)
      has.negative <- vec.qtile[1L] < 0

      expr <- switch(s$norm,
        identity = dat.raw$expr,
        voom = {
          fit <- limma::voom(dat.raw$expr, design = matrix(1, ncol(
                dat.raw$expr),
              1))
          fit$E
        },
        ## 相邻字符串在 R 里不会自动拼接，必须自己 paste0
        rma = stop(paste0(
            "S03: RMA branch removed. CEL path was dropped because ",
            "pd.hgfocus / pd.hgu133a2 / pd.hugene10kplus have been ",
            "removed from Bioconductor entirely; this series now uses ",
            "source='series_matrix'.")),
        loess = dat.raw$expr,
        stop(paste0("S03: unknown norm '", s$norm, "'")))

  ## 行级缺失过滤：probe 向量必须按同一逻辑索引同步子集，否则行轴与探针错位。
      lgl.keep <- apply(expr, 1L,
        function(x) sum(is.na(x)) / length(x) <= MAX_MISSING_RATE)
      expr <- expr[lgl.keep, , drop = FALSE]
      vec.probe <- dat.raw$probe[lgl.keep]

      ## 2.3 基因级矩阵（新增；探针级契约 p1..pN / s1..sN 原样保留，两者并存）
      vec.sym <- SymbolFromProbe(a.annot, vec.probe)
      expr.gene <- BuildGeneMatrix(expr, vec.sym, a.annot$sep)
      colnames(expr.gene) <- vec.sample
      n.probe.mapped <- sum(HasUsableSymbol(vec.sym, a.annot$sep))
      n.multi <- sum(grepl(a.annot$sep, vec.sym, fixed = TRUE))

      ## 2.4 断言：探针级 4 项 + 基因级 4 项（§9 Q2 与新增 Q3 同口径）
      rownames(expr) <- paste0("p", seq_len(nrow(expr)))
      colnames(expr) <- paste0("s", seq_len(ncol(expr)))
      vec.qc <- c(
        "non_empty" = nrow(expr) > 0L && ncol(expr) > 0L,
        "no_all_na" = !any(apply(expr, 1L, function(x) all(is.na(x)))),
        "min_probes" = nrow(expr) >= MIN_FEATURE_N,
        "finite" = all(is.finite(expr[is.na(expr) == FALSE])),
        "gene_non_empty" = nrow(expr.gene) > 0L && ncol(expr.gene) > 0L,
    "gene_no_all_na" = !any(apply(expr.gene, 1L, function(x) all(is.na(x)))),
        "gene_min_symbols" = nrow(expr.gene) >= MIN_FEATURE_N,
        "gene_finite" = all(is.finite(expr.gene[is.na(expr.gene) == FALSE]))
      )
      for (nm in names(vec.qc)) {
        dat.norm.qc <- rbind(dat.norm.qc, data.frame(
            acc = acc, qc = nm, observed = if (vec.qc[nm]) "PASS" else "FAIL",
            expected = "PASS", stringsAsFactors = FALSE))
      }

      ## 2.5 落盘（断点续传：两个产物都已存在则跳过，见 DoneSeries）
      saveRDS(list(expr = expr, probe = vec.probe, acc = acc, norm = s$norm),
        ExprFile(acc))
      saveRDS(list(expr = expr.gene, symbol = rownames(expr.gene),
          sample = colnames(expr.gene), probe = vec.probe,
          acc = acc, norm = s$norm), GeneFile(acc))
      dat.norm.log <- rbind(dat.norm.log, data.frame(
          acc = acc, source = s$source, norm = s$norm, n_probe = nrow(expr),
          n_sample = ncol(expr),
          ## na.rm：双色 M 值含 NA（GPR 空点），量程登记只统计非 NA 值；
          ## 无 na.rm 时 min/max 返回 NA，会传染下游 vec.brk 派生致 seq(NaN) 崩
          range_min = min(expr, na.rm = TRUE), range_max = max(expr,
            na.rm = TRUE),
    has_negative = has.negative, q001 = vec.qtile[1L], qmed = vec.qtile[2L],
          q999 = vec.qtile[3L],
    cdf_platform = NA_character_, note = s$note, stringsAsFactors = FALSE))
      dat.gene.log <- rbind(dat.gene.log, data.frame(
          acc = acc, platform = a.annot$platform, n_probe = nrow(expr),
          n_probe_mapped = n.probe.mapped, n_symbol = nrow(expr.gene),
          n_multi_symbol_probe = n.multi, stringsAsFactors = FALSE))
      ## §9 Q4：面板缺失逐条记录，禁止同义替换（故此处只登记存在性，不改写符号）
      dat.panel.cov <- rbind(dat.panel.cov, data.frame(
          acc = acc, gene = PANEL_GENES,
    present = ifelse(PANEL_GENES %in% rownames(expr.gene), "present", "absent"),
          stringsAsFactors = FALSE))
      list(ok = TRUE)
    }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  if (!res$ok) {
    cat(sprintf("[S03] %s  FAILED: %s\n", acc, res$msg))
    dat.norm.qc <- rbind(dat.norm.qc, data.frame(
        acc = acc, qc = "ingest_or_normalize", observed = "FAIL",
        expected = "PASS", stringsAsFactors = FALSE))
    next
  }
}

write.csv(dat.norm.log, file.norm.log, row.names = FALSE)
write.csv(dat.norm.qc, file.norm.qc, row.names = FALSE)
write.csv(dat.gene.log, file.gene.log, row.names = FALSE)
write.csv(dat.panel.cov, file.panel.cov, row.names = FALSE)

## PART 2.5 图件 ---------------------------------------------------------------------------
## 命名遵守规范（描述性英文，禁 Fig 前缀；图内文字全英文）；不引入拼图包，单面板出图。
## 四幅图分别对应 三级特征保留 / 表达量程 / 断言状态 / 面板覆盖（见文件头图件清单）。

## 系列统一排序：按探针级特征数升序，供图 1 / 图 2 / 图 4 共用。
vec.ord <- dat.norm.log$acc[order(dat.norm.log$n_probe)]

## 图 1：三级特征保留（保留探针 / 已映射探针 / 唯一基因符号），单面板分组柱。
## 替代原「只画探针数、按归一化方法着色」——归一化方法已在 normalization_log.csv 登记，
## 用它着色不携带新信息；改为呈现 探针 -> 注释 -> 基因 的逐级损耗（基因级矩阵的规模依据）。
dat.f1 <- rbind(
  data.frame(acc = dat.norm.log$acc, n = dat.norm.log$n_probe,
    level = "Probes kept", stringsAsFactors = FALSE),
  data.frame(acc = dat.gene.log$acc, n = dat.gene.log$n_probe_mapped,
    level = "Probes mapped to symbol", stringsAsFactors = FALSE),
  data.frame(acc = dat.gene.log$acc, n = dat.gene.log$n_symbol,
    level = "Unique gene symbols", stringsAsFactors = FALSE))
dat.f1 <- dat.f1[dat.f1$n > 0L, , drop = FALSE]
dat.f1$acc <- factor(dat.f1$acc, levels = vec.ord)
dat.f1$level <- factor(dat.f1$level,
  levels = c("Probes kept", "Probes mapped to symbol",
    "Unique gene symbols"))
p.metric <- ggplot(dat.f1, aes(x = acc, y = n, fill = level)) +
geom_col(position = position_dodge(width = 0.8), width = 0.75) +
labs(title = "Feature retention per series: probes to gene symbols",
  subtitle = paste0("Bars: probes kept after normalisation, probes mapped to",
    "a symbol, unique symbols"),
  x = "Series", y = "Features retained", fill = "Feature level") +
sci_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
sci_save(p.metric, "Normalization_Metrics_By_Series", width = 9.5, height = 5)

## 图 2：各系列表达量程（单面板水平量程图：每系列一条 min-max 线段 + 端点，按已登记的
## 归一化路径着色）。
## 用 pseudo-log（symlog）横轴：各系列量程差异达两个数量级（GSE41571 原始强度 0-3664，
## 其余多 <= 20），且 GSE47472/GSE323983 含负值；symlog 可同时压缩大值、展开小值并容纳负值，
## 避免多数系列被压成一线而不可判读。基因级矩阵由同一 expr 归并而来，量程口径相同，不另出图。
## 着色口径（2026-09-29 收紧）：改按**已登记**的归一化路径着色，不再用 `has_negative`
## 冒充「尺度判定」——该列只能证明「存在负值 → 必为 log 尺度」，不能反证「不是 log2」。
dat.norm.plot <- dat.norm.log[dat.norm.log$n_probe > 0L, , drop = FALSE]
dat.norm.plot$acc <- factor(dat.norm.plot$acc, levels = vec.ord)
##   `breaks = c(-10, 0, 10, 100, 1000, 4000)`，其中 4000 只对应 GSE41571 一个系列的上界
##   （该系列原始强度 0–3664），其余系列上界均 <= 20 —— 一旦系列增减或量程变化即失准，
##   属规范 28 禁止的魔法数字。新写法取 10 的整数次幂覆盖到实测最大值，并补最小与零。
##   min/max 返回 NA -> log10(NA)=NaN -> seq(0, NaN) 报 'to' must be finite。
vec.brk <- unique(c(min(dat.norm.plot$range_min, na.rm = TRUE), 0,
    10^seq(0, floor(log10(max(dat.norm.plot$range_max, na.rm = TRUE))))))
p.range <- ggplot(dat.norm.plot, aes(y = acc, color = norm)) +
geom_linerange(aes(xmin = range_min, xmax = range_max), linewidth = 1.6) +
geom_point(aes(x = range_min), size = 1.4) +
geom_point(aes(x = range_max), size = 1.4) +
scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
  breaks = vec.brk) +
labs(title = "Expression range per series",
  subtitle = paste0("Span = min-max (pseudo-log); colour = registered",
    "normalisation path"),
  x = "Expression value (pseudo-log scale)", y = "Series",
  color = "Normalisation") +
sci_theme
sci_save(p.range, "Expression_Range_Per_Series", width = 8, height = 5.5)

## 图 3：断言状态热图（单面板：x = 系列，y = 断言，tile 着色 PASS/FAIL）。
## 替代原把断言放 y、把 result 放 x 且 facet 逐系列（每面板仅 4 个 tile）的不可读写法。
## 断言集扩为：探针级 4 项（§9 Q2）→ 基因级 4 项（新增 Q3，同口径）→ 取数/归一化总状态。
dat.qc.plot <- dat.norm.qc[!is.na(dat.norm.qc$observed), , drop = FALSE]
dat.qc.plot$observed <- factor(dat.qc.plot$observed, levels = c("PASS", "FAIL"))
dat.qc.plot$qc <- factor(dat.qc.plot$qc,
  levels = c("non_empty", "no_all_na", "min_probes", "finite",
    "gene_non_empty", "gene_no_all_na",
    "gene_min_symbols", "gene_finite",
    "ingest_or_normalize"))
## levels = vec.ord 即可：GSE6002 正常入表后 vec.ord 已含它（旧写法
## c(vec.ord, "GSE6002") 是其 QC FAIL 时代的手动补轴，现成重复 level 报错）
dat.qc.plot$acc <- factor(dat.qc.plot$acc, levels = vec.ord)
p.qc <- ggplot(dat.qc.plot, aes(x = acc, y = qc, fill = observed)) +
geom_tile(color = "white") +
labs(title = paste0("Expression-matrix assertions per series (probe level",
  "and gene level)"),
  x = "Series", y = "Assertion", fill = "Result") +
scale_fill_manual(values = c("PASS" = "#4C72B0", "FAIL" = "#C44E52"),
  name = "Result") +
sci_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
sci_save(p.qc, "Normalization_QC_By_Series", width = 9.5, height = 5.5)

## 图 4：预注册 12 基因面板逐系列存在性（§9 Q4「缺失逐条记录，禁止同义替换」的可视凭据）。
## 该图直接暴露层 1（S05）的面板可及性风险：哪些系列缺哪个基因一目了然。
dat.pc <- dat.panel.cov
dat.pc$acc <- factor(dat.pc$acc, levels = vec.ord)
dat.pc$gene <- factor(dat.pc$gene, levels = rev(PANEL_GENES))
dat.pc$present <- factor(dat.pc$present, levels = c("present", "absent"))
p.panel <- ggplot(dat.pc, aes(x = acc, y = gene, fill = present)) +
geom_tile(color = "white") +
labs(title = "Pre-registered 12-gene panel coverage per series",
  subtitle = paste0("Absent = symbol not annotated on that platform",
    "(protocol Q4 forbids synonym substitution)"),
  x = "Series", y = "Panel gene", fill = "Status") +
scale_fill_manual(values = c("present" = "#4C72B0", "absent" = "#C44E52"),
  name = "Status") +
sci_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
sci_save(p.panel, "Panel_Gene_Coverage_By_Series", width = 9, height = 5)

## PART 3 收口 -----------------------------------------------------------------------------

vec.qc.fail <- unique(dat.norm.qc$acc[dat.norm.qc$observed == "FAIL"])
message(paste0("\n===== S03 normalize summary =====\nseries processed = ", nrow(
      dat.norm.log),
    "  qc failures = ", length(vec.qc.fail)))
if (length(vec.qc.fail) > 0L) {
  message(paste0("[QC FAIL] ", paste(vec.qc.fail, collapse = ", ")))
  print(dat.norm.qc[dat.norm.qc$observed == "FAIL", ])
}
if (length(vec.qc.fail) > 0L) {
  stop(paste0("S03 FAILED: expression-matrix QC failed for ",
      paste(vec.qc.fail,
        collapse = ", "), " ; see 03_results/S03/table/ ."))
}

## §9 Q4 面板可及性披露（不 gate）：缺失逐条列出，交由层 1（S05）与用户决策，
## 本脚本一律不做同义替换（Q4 明文禁止）。
dat.miss <- dat.panel.cov[dat.panel.cov$present == "absent", , drop = FALSE]
message(paste0("\n===== S03 gene-level summary =====\nseries = ", nrow(
      dat.gene.log),
    "  median unique symbols = ",
    round(stats::median(dat.gene.log$n_symbol), 0),
    "  panel absent entries = ", nrow(dat.miss)))
if (nrow(dat.miss) > 0L) {
  message("[Q4] panel gene absences (recorded, no synonym substitution):")
  print(dat.miss)
}
message(
  "S03 OK. tables -> 03_results/S03/table/ ; figures -> 03_results/S03/Figure/")
