################################################################
#  CRLM 泛素相关基因研究 - Step 1: 数据获取与预处理
#
#  数据集概览:
#    GSE49355  (训练集,   n=57)  → GPL96    → raw/GPL96.annot.gz
#    GSE159216 (验证集1, n=283) → GPL17586 → raw/GPL17586.txt
#    GSE39582  (验证集2, n=585) → GPL570   → raw/GPL570.annot.gz
#    GSE72718  (独立集,  n=28)  → GPL17586 → raw/GPL17586.txt
#
#  GPL 注释策略:
#    GPL96 / GPL570  : 读取本地 raw/GPL*.annot.gz  (SOFT annot 格式)
#    GPL17586        : 读取本地 raw/GPL17586.txt
#                      解析 gene_assignment 列提取 Gene Symbol
#                      解析结果缓存至 raw/GPL17586_anno.rds（第一次运行约1分钟）
#
#  输出目录:
#    data/processed/  → 标准化矩阵、临床信息、批次校正结果 (RDS)
#    results/figures/ → QC 可视化图表 (PDF, QC1~QC17)
#    results/tables/  → 汇总统计表格 (CSV)
################################################################

# ============================================================
# Section 0: 加载依赖包
# ============================================================
library(GEOquery)       # getGEO() 读取 GEO series_matrix 文件
library(limma)          # normalizeBetweenArrays(), lmFit(), eBayes()
library(ggplot2)        # 所有 ggplot2 图表
library(ggrepel)        # geom_text_repel() 文字标注防重叠
library(pheatmap)       # pheatmap() 热图
library(dplyr)          # 数据操作管道
library(tibble)         # column_to_rownames()
library(Biobase)        # exprs(), pData(), annotation()
library(WGCNA)          # 备用：goodSamplesGenes() 质控
library(RColorBrewer)   # brewer.pal() 配色
library(scales)         # scales::comma() 数字格式化
library(tidyr)          # pivot_longer()

# ============================================================
# Section 1: 工作目录与路径设置
# ============================================================
setwd("X:/CRLM")

# 创建输出目录（已存在时静默跳过）
dir.create("data/processed",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)

# 原始数据目录：所有 GEO 文件和 GPL 注释文件均位于此处
raw_dir <- "raw"

# ============================================================
# Section 2: GPL 平台注释函数定义
# ============================================================

# ------------------------------------------------------------
# 函数 A: 读取本地 annot.gz 文件 (适用于 GPL96 / GPL570)
# 文件格式: SOFT annot 格式，以 ^/!/# 开头的行为元数据，其余为数据
# ------------------------------------------------------------
read_annot_gz <- function(annot_file) {
  con   <- gzfile(annot_file, "rt")
  lines <- readLines(con, warn = FALSE)
  close(con)

  # 找到第一行不以 ^ ! # 开头的行，即表头行
  data_start <- which(
    substr(lines, 1, 1) != "^" &
    substr(lines, 1, 1) != "!" &
    substr(lines, 1, 1) != "#"
  )[1]
  if (is.na(data_start)) {
    stop("read_annot_gz: 无法找到数据起始行，请检查文件: ", annot_file)
  }

  header     <- strsplit(lines[data_start], "\t")[[1]]
  data_lines <- lines[(data_start + 1):length(lines)]
  data_lines <- data_lines[nchar(trimws(data_lines)) > 0]  # 去除空行

  df_list <- lapply(data_lines, function(l) {
    parts <- strsplit(l, "\t")[[1]]
    # 列数不足时补 NA，列数超出时截断
    if (length(parts) < length(header)) {
      parts <- c(parts, rep(NA_character_, length(header) - length(parts)))
    } else if (length(parts) > length(header)) {
      parts <- parts[seq_len(length(header))]
    }
    parts
  })

  df <- as.data.frame(
    do.call(rbind, df_list),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  colnames(df) <- header

  cat("  [annot.gz] 读取完成:", basename(annot_file),
      "| 行:", nrow(df), "| 列:", ncol(df), "\n")
  return(df)
}

# ------------------------------------------------------------
# 函数 B: 解析 GPL17586.txt 平台注释文件
# 适用于 GSE159216 (GPL17586) 和 GSE72718 (GPL17586)
# gene_assignment 列格式: "NM_xxx // GeneSymbol // desc /// ..."
# 解析结果缓存为 GPL17586_anno.rds，避免重复解析（约1分钟）
# ------------------------------------------------------------
read_gpl17586_family_soft <- function(soft_file) {
  rds_cache <- sub("\\.txt$", "_anno.rds", soft_file)

  # 优先读取缓存
  if (file.exists(rds_cache)) {
    cat("  [GPL17586] 读取缓存注释:", basename(rds_cache), "\n")
    return(readRDS(rds_cache))
  }

  cat("  [GPL17586] 首次解析（约1分钟）:", basename(soft_file), "\n")

  # 内部函数：从 gene_assignment 字段提取第一个 Gene Symbol
  extract_sym <- function(ga) {
    if (is.na(ga) || ga == "" || ga == "---") return(NA_character_)
    # 取第一个 /// 分隔的 entry
    first_entry <- strsplit(ga, " /// ", fixed = TRUE)[[1]][1]
    # 再按 // 分隔，第二个字段是 Symbol
    parts <- strsplit(first_entry, " // ", fixed = TRUE)[[1]]
    if (length(parts) < 2) return(NA_character_)
    sym <- trimws(parts[2])
    if (sym == "---" || sym == "") return(NA_character_)
    return(sym)
  }

  # 读取文件头部，定位 !platform_table_begin 标记
  con      <- file(soft_file, "rt")
  preamble <- readLines(con, n = 30L, warn = FALSE)
  tb_idx <- which(grepl("^!platform_table_begin", preamble, ignore.case = TRUE))
  if (length(tb_idx) == 0L) {
    close(con)
    stop("GPL17586.txt: 未找到 !platform_table_begin 标记（前30行内）")
  }
  tb_pos      <- tb_idx[1L]
  header_line <- preamble[tb_pos + 1L]
  header      <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  id_idx <- which(header == "ID")
  ga_idx <- which(header == "gene_assignment")
  if (length(id_idx) == 0L || length(ga_idx) == 0L) {
    close(con)
    stop("GPL17586 注释表缺少 'ID' 或 'gene_assignment' 列。实际列名: ",
         paste(header, collapse = ", "))
  }
  cat("  [GPL17586] 列索引 | ID:", id_idx, "| gene_assignment:", ga_idx, "\n")
  cat("  [GPL17586] 逐行读取探针数据...\n")

  # 预分配向量（GPL17586 约 67K 条探针）
  probe_vec  <- character(100000L)
  symbol_vec <- character(100000L)
  n <- 0L

  # 逐行处理函数，遇到 !platform_table_end 时返回 FALSE 停止
  process_line <- function(line) {
    if (grepl("^!platform_table_end", line, ignore.case = TRUE)) return(FALSE)
    parts          <- strsplit(line, "\t", fixed = TRUE)[[1]]
    n             <<- n + 1L
    probe_vec[n]  <<- if (length(parts) >= id_idx) parts[id_idx]  else NA_character_
    symbol_vec[n] <<- extract_sym(
                        if (length(parts) >= ga_idx) parts[ga_idx] else NA_character_
                      )
    return(TRUE)
  }

  # 先处理 preamble 中 table_begin 之后的剩余行
  if ((tb_pos + 2L) <= length(preamble)) {
    for (ln in preamble[(tb_pos + 2L):length(preamble)]) {
      if (!process_line(ln)) {
        close(con)
        break
      }
    }
  }

  # 继续从文件流读取剩余行
  if (isOpen(con)) {
    repeat {
      line <- readLines(con, n = 1L, warn = FALSE)
      if (length(line) == 0L) { close(con); break }
      if (!process_line(line)) { close(con); break }
      if (n %% 50000L == 0L) cat("  [GPL17586] 已读取", n, "条探针...\n")
    }
  }

  anno <- data.frame(
    probe  = probe_vec[1:n],
    symbol = symbol_vec[1:n],
    stringsAsFactors = FALSE
  )
  # 去除无法映射的探针
  anno <- anno[!is.na(anno$symbol) & anno$symbol != "", ]

  cat("  [GPL17586] 解析完成 | 总探针:", n,
      "| 有symbol:", nrow(anno),
      "| 唯一基因:", length(unique(anno$symbol)), "\n")
  saveRDS(anno, rds_cache)
  cat("  [GPL17586] 注释缓存已写入:", basename(rds_cache), "\n")
  return(anno)
}

# ------------------------------------------------------------
# 函数 C: 读取 series_matrix.txt.gz，并自动关联本地 annot 注释
# 用于 GPL96 (GSE49355) 和 GPL570 (GSE39582)
# 注意: GPL17586 的数据集单独调用 getGEO，注释在 Part 2 统一加载
# ------------------------------------------------------------
read_geo_no_auto_gpl <- function(gz_file) {
  gse_tmp  <- getGEO(filename = gz_file, getGPL = FALSE)
  gse_obj  <- if (is.list(gse_tmp)) gse_tmp[[1]] else gse_tmp
  gpl_name <- annotation(gse_obj)
  cat("  [GEO] 检测到 GPL 平台:", gpl_name, "\n")

  local_annot <- file.path(raw_dir, paste0(gpl_name, ".annot.gz"))
  anno_result <- NULL

  if (file.exists(local_annot)) {
    cat("  [GEO] 找到本地 annot 文件:", basename(local_annot), "\n")
    anno_df <- read_annot_gz(local_annot)
    # 第一列是探针ID，找 Gene Symbol 列（含 "symbol" 字样）
    id_col  <- colnames(anno_df)[1]
    sym_col <- grep("symbol", colnames(anno_df), value = TRUE, ignore.case = TRUE)[1]
    cat("  [GEO] 探针列:", id_col, "| Symbol列:", sym_col, "\n")
    if (!is.na(sym_col)) {
      anno_sub           <- anno_df[, c(id_col, sym_col), drop = FALSE]
      colnames(anno_sub) <- c("probe", "symbol")
      # 去除空值和多mapping探针（含 /// 的为多基因映射）
      anno_sub <- anno_sub[anno_sub$symbol != "" & !is.na(anno_sub$symbol), ]
      anno_sub <- anno_sub[!grepl("///", anno_sub$symbol), ]
      anno_result <- anno_sub
    }
  } else {
    cat("  [GEO] 警告: 未找到本地 annot 文件，anno 将为 NULL\n")
  }

  return(list(gse = gse_obj, anno = anno_result))
}

# ------------------------------------------------------------
# 函数 D: 探针矩阵 → 基因Symbol矩阵（多探针取均值）
# 参数:
#   expr_probe_matrix : 行为探针ID、列为样本的表达矩阵
#   anno_df           : data.frame，含 probe 和 symbol 两列
# 返回: 行为 Gene Symbol、列为样本的表达矩阵（多探针均值）
# ------------------------------------------------------------
map_probes_to_genes <- function(expr_probe_matrix, anno_df) {
  probes_in_expr <- rownames(expr_probe_matrix)
  # 只保留注释表与表达矩阵共有的探针
  probes_matched <- anno_df$probe[anno_df$probe %in% probes_in_expr]
  cat("  [mapProbes] 探针匹配:", length(probes_matched),
      "/", length(probes_in_expr), "个探针有注释\n")

  expr_sub  <- expr_probe_matrix[probes_matched, , drop = FALSE]
  anno_sub  <- anno_df[match(probes_matched, anno_df$probe), , drop = FALSE]

  # 合并 symbol 列，按 symbol 取多探针均值
  df         <- as.data.frame(expr_sub, stringsAsFactors = FALSE)
  df$symbol  <- anno_sub$symbol

  result <- df %>%
    group_by(symbol) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
    column_to_rownames("symbol") %>%
    as.matrix()

  return(result)
}

# ============================================================
# Section 3: 读取本地 GEO series_matrix 文件
# ============================================================

# ------------------------------------------------------------
# GSE49355: 训练集，n=57，GPL96 芯片（Affymetrix HG-U133A）
# 三组: Normal Colon / Primary Tumor / Liver Metastasis
# ------------------------------------------------------------
cat("\n>>> [Step1] 加载 GSE49355 (GPL96)...\n")
tmp49355    <- read_geo_no_auto_gpl(file.path(raw_dir, "GSE49355_series_matrix.txt.gz"))
gse49355    <- tmp49355$gse
anno49355   <- tmp49355$anno
expr49355_probe <- exprs(gse49355)
pheno49355      <- pData(gse49355)
cat("GSE49355 | 探针数:", nrow(expr49355_probe),
    "| 样本数:", ncol(expr49355_probe), "\n")
cat("GSE49355 | 有效注释探针数:", nrow(anno49355), "\n")

# ------------------------------------------------------------
# GSE159216: 验证集1，n=283，GPL17586 芯片（Affymetrix HTA 2.0）
# 注: GPL17586 注释单独在 Section 4 处理
# ------------------------------------------------------------
cat("\n>>> [Step1] 加载 GSE159216 (GPL17586)...\n")
# 【注意】GPL17586 的注释文件较大，此处只读表达矩阵，注释在 Section 4 处理
gse159216_tmp    <- getGEO(filename = file.path(raw_dir, "GSE159216_series_matrix.txt.gz"),
                           getGPL = FALSE)
gse159216        <- if (is.list(gse159216_tmp)) gse159216_tmp[[1]] else gse159216_tmp
expr159216_probe <- exprs(gse159216)
pheno159216      <- pData(gse159216)
cat("GSE159216 | 探针数:", nrow(expr159216_probe),
    "| 样本数:", ncol(expr159216_probe), "\n")

# ------------------------------------------------------------
# GSE39582: 验证集2，n=585，GPL570 芯片（Affymetrix HG-U133 Plus 2.0）
# 含丰富临床信息: TNM 分期、分子亚型、OS 生存数据
# ------------------------------------------------------------
cat("\n>>> [Step1] 加载 GSE39582 (GPL570)...\n")
tmp39582    <- read_geo_no_auto_gpl(file.path(raw_dir, "GSE39582_series_matrix.txt.gz"))
gse39582    <- tmp39582$gse
anno39582   <- tmp39582$anno
expr39582_probe <- exprs(gse39582)
pheno39582      <- pData(gse39582)
cat("GSE39582 | 探针数:", nrow(expr39582_probe),
    "| 样本数:", ncol(expr39582_probe), "\n")
cat("GSE39582 | 有效注释探针数:", nrow(anno39582), "\n")

# ------------------------------------------------------------
# GSE72718: 独立验证集，n=28，GPL17586 芯片
# 含肝转移 vs 原发灶对比样本
# ------------------------------------------------------------
cat("\n>>> [Step1] 加载 GSE72718 (GPL17586)...\n")
gse72718_tmp    <- getGEO(filename = file.path(raw_dir, "GSE72718_series_matrix.txt.gz"),
                          getGPL = FALSE)
gse72718        <- if (is.list(gse72718_tmp)) gse72718_tmp[[1]] else gse72718_tmp
expr72718_probe <- exprs(gse72718)
pheno72718      <- pData(gse72718)
cat("GSE72718 | 探针数:", nrow(expr72718_probe),
    "| 样本数:", ncol(expr72718_probe), "\n")

# ============================================================
# Section 4: 探针 → Gene Symbol 映射
# ============================================================

# ------------------------------------------------------------
# GPL17586 注释（GSE159216 和 GSE72718 共用）
# 首次运行会解析 GPL17586.txt 并缓存，后续自动读缓存
# ------------------------------------------------------------
cat("\n>>> [Step1] 加载 GPL17586 注释 (raw/GPL17586.txt)...\n")
anno17586 <- read_gpl17586_family_soft(file.path(raw_dir, "GPL17586.txt"))

# GSE49355 (GPL96) 探针 → 基因映射
cat("\n>>> [Step1] GSE49355 探针注释映射 (GPL96)...\n")
expr49355 <- map_probes_to_genes(expr49355_probe, anno49355)
cat("GSE49355 | 映射后基因数:", nrow(expr49355),
    "| 样本数:", ncol(expr49355), "\n")

# GSE159216 (GPL17586) 探针 → 基因映射
cat("\n>>> [Step1] GSE159216 探针注释映射 (GPL17586)...\n")
expr159216 <- map_probes_to_genes(expr159216_probe, anno17586)
cat("GSE159216 | 映射后基因数:", nrow(expr159216),
    "| 样本数:", ncol(expr159216), "\n")

# GSE39582 (GPL570) 探针 → 基因映射
cat("\n>>> [Step1] GSE39582 探针注释映射 (GPL570)...\n")
expr39582 <- map_probes_to_genes(expr39582_probe, anno39582)
cat("GSE39582 | 映射后基因数:", nrow(expr39582),
    "| 样本数:", ncol(expr39582), "\n")

# GSE72718 (GPL17586) 探针 → 基因映射
cat("\n>>> [Step1] GSE72718 探针注释映射 (GPL17586)...\n")
expr72718 <- map_probes_to_genes(expr72718_probe, anno17586)
cat("GSE72718 | 映射后基因数:", nrow(expr72718),
    "| 样本数:", ncol(expr72718), "\n")

# ============================================================
# Section 5: 样本分组信息构建
# ============================================================

# ------------------------------------------------------------
# GSE49355 样本分组
# title 列格式示例:
#   "Liver_Metastasis-HM_016_MV"  → LM
#   "Primary_Tumor-CT_016_MV"     → PT
#   "Normal_Colon-CN_016_MV"      → Normal
# 【核查点】运行后请打印 table(group49355) 确认分组数目
# ------------------------------------------------------------
cat("\n>>> [Step1] GSE49355 样本分组...\n")
cat("  title 前3个示例:", head(pheno49355$title, 3), "\n")

group49355 <- ifelse(
  grepl("^[Ll]iver[_ ]?[Mm]etasta", pheno49355$title), "LM",
  ifelse(
    grepl("^[Pp]rimary[_ ]?[Tt]umor", pheno49355$title), "PT",
    ifelse(
      grepl("^[Nn]ormal[_ ]?[Cc]olon", pheno49355$title), "Normal",
      "Unclassified"   # 原代码所有非LM非PT全归Normal，改为显式标注
    )
  )
)
# names 必须与 expr49355 列名一致，后续图表着色依赖此命名
names(group49355) <- colnames(expr49355)

cat("\nGSE49355 分组统计（应为 Normal/PT/LM 三组，无 Unclassified）:\n")
print(table(group49355, useNA = "ifany"))

# 若出现 Unclassified，请打印 pheno49355$title 核查 title 格式
if (any(group49355 == "Unclassified")) {
  warning("GSE49355 存在无法分类的样本，请人工检查 pheno49355$title")
  print(pheno49355$title[group49355 == "Unclassified"])
}
pheno49355$group <- group49355

# ------------------------------------------------------------
# GSE39582 临床信息提取
# 关键列: dataset:ch1 / tnm.stage:ch1 / cit.molecularsubtype:ch1
#          os.event:ch1 / os.delay (months):ch1
#          rfs.event:ch1 / rfs.delay:ch1
# 【注意】使用 [[ ]] 访问含特殊字符（冒号、空格、括号）的列名
# ------------------------------------------------------------
cat("\n>>> [Step1] GSE39582 临床信息提取...\n")

# 先打印关键列供人工核查（【核查点】）
clin_cols <- grep("characteristics|os|surv|metastas|status|tnm|subtype|dataset",
                  colnames(pheno39582), ignore.case = TRUE, value = TRUE)
cat("  GSE39582 可用临床列:\n")
print(clin_cols)
cat("  前5行样本临床信息:\n")
print(head(pheno39582[, clin_cols[1:min(6, length(clin_cols))]], 5))

# dataset:ch1 → 样本类型 (Normal / Tumor)
ds39582_raw <- pheno39582[["dataset:ch1"]]
ds39582 <- ifelse(
  grepl("Non Tumoral", ds39582_raw, ignore.case = TRUE), "Normal", "Tumor"
)
names(ds39582) <- colnames(expr39582)
cat("  dataset 分布:\n")
print(table(ds39582, useNA = "ifany"))

# tnm.stage:ch1 → TNM 分期
tnm39582_raw <- pheno39582[["tnm.stage:ch1"]]
# 去除列名前缀（如 "tnm.stage: 2" → "2"）
tnm39582 <- gsub("^tnm\\.stage:\\s*", "", tnm39582_raw)
tnm39582 <- ifelse(tnm39582 %in% c("N/A", "", "NA"), NA_character_, tnm39582)
tnm39582 <- factor(tnm39582, levels = c("0", "1", "2", "3", "4"))
names(tnm39582) <- colnames(expr39582)
cat("  TNM stage 分布:\n")
print(table(tnm39582, useNA = "ifany"))

# cit.molecularsubtype:ch1 → 分子亚型 (CIT C1-C6)
ms39582_raw <- pheno39582[["cit.molecularsubtype:ch1"]]
ms39582 <- gsub("^cit\\.molecularsubtype:\\s*", "", ms39582_raw)
ms39582 <- ifelse(ms39582 %in% c("N/A", "", "NA"), NA_character_, ms39582)
ms39582 <- factor(ms39582, levels = c("C1", "C2", "C3", "C4", "C5", "C6"))
names(ms39582) <- colnames(expr39582)
cat("  Molecular subtype 分布:\n")
print(table(ms39582, useNA = "ifany"))

# os.event:ch1 → OS 生存事件 (1=死亡, 0=存活)
os_event39582_raw <- pheno39582[["os.event:ch1"]]
os_event39582 <- gsub("^os\\.event:\\s*", "", os_event39582_raw)
os_event39582 <- ifelse(os_event39582 %in% c("N/A", "", "NA"), NA_character_, os_event39582)
names(os_event39582) <- colnames(expr39582)
cat("  OS event 分布:\n")
print(table(os_event39582, useNA = "ifany"))

# ============================================================
# Section 6: 数据标准化（Log2 转换 + 分位数标准化）
# ============================================================
# 逻辑说明:
#   步骤1: 判断是否需要 log2 转换（max > 100 认为是原始强度）
#   步骤2: 在 log2 转换前保存快照（_pre），用于后续可视化对比
#          原代码先赋值 _pre 再做转换，_pre 存的是未转换数据——
#          实际上 _pre 的目的是存"log2转换后、分位数标准化前"的状态，
#          因此正确顺序是：log2 → 保存为 _pre → quantile normalize
#   步骤3: normalizeBetweenArrays(..., method="quantile") 分位数标准化

cat("\n>>> [Step1] 数据标准化...\n")

# ----- GSE49355 -----
if (max(expr49355, na.rm = TRUE) > 100) {
  expr49355 <- log2(expr49355 + 1)
  cat("GSE49355: 已执行 log2(x+1) 转换\n")
} else {
  cat("GSE49355: 判断已为 log2 scale，跳过转换\n")
}
expr49355_pre  <- expr49355                          # 保存 log2 后、QN 前的快照
expr49355_norm <- normalizeBetweenArrays(expr49355, method = "quantile")
cat("GSE49355: 分位数标准化完成\n")

# ----- GSE159216 -----
if (max(expr159216, na.rm = TRUE) > 100) {
  expr159216 <- log2(expr159216 + 1)
  cat("GSE159216: 已执行 log2(x+1) 转换\n")
} else {
  cat("GSE159216: 判断已为 log2 scale，跳过转换\n")
}
expr159216_pre  <- expr159216
expr159216_norm <- normalizeBetweenArrays(expr159216, method = "quantile")
cat("GSE159216: 分位数标准化完成\n")

# ----- GSE39582 -----
if (max(expr39582, na.rm = TRUE) > 100) {
  expr39582 <- log2(expr39582 + 1)
  cat("GSE39582: 已执行 log2(x+1) 转换\n")
} else {
  cat("GSE39582: 判断已为 log2 scale，跳过转换\n")
}
expr39582_pre  <- expr39582
expr39582_norm <- normalizeBetweenArrays(expr39582, method = "quantile")
cat("GSE39582: 分位数标准化完成\n")

# ----- GSE72718 -----
if (max(expr72718, na.rm = TRUE) > 100) {
  expr72718 <- log2(expr72718 + 1)
  cat("GSE72718: 已执行 log2(x+1) 转换\n")
} else {
  cat("GSE72718: 判断已为 log2 scale，跳过转换\n")
}
expr72718_pre  <- expr72718
expr72718_norm <- normalizeBetweenArrays(expr72718, method = "quantile")
cat("GSE72718: 分位数标准化完成\n")

# ============================================================
# Section 7: 数据质控可视化
# ============================================================

# 预计算 GSE49355 的高变异基因集（后续多处复用）
rv_49355    <- apply(expr49355_norm, 1, var, na.rm = TRUE)
top5k_49355 <- names(sort(rv_49355, decreasing = TRUE))[1:5000]

# 定义统一颜色与形状（整个 Section 7 复用）
sci_colors <- c("Normal" = "#0072B2", "PT" = "#D55E00", "LM" = "#CC79A7")
sci_shapes <- c("Normal" = 16,        "PT" = 17,        "LM" = 15)

# ------------------------------------------------------------
# QC1: GSE49355 样本层次聚类树状图（带底部分组色带）
# 目的: 验证样本按组织类型正确聚类，排查混杂样本
# ------------------------------------------------------------
cat("\n>>> 生成 QC1: GSE49355 层次聚类...\n")
sampleDist_49355 <- dist(t(expr49355_norm[top5k_49355, ]), method = "euclidean")
hc_49355         <- hclust(sampleDist_49355, method = "average")

# dend_colors 应按聚类顺序（hc_49355$order 是整数位置索引）
# group49355 是有names的字符向量，用整数位置索引是合法的
# 但保险起见用 names(group49355)[hc_49355$order] 取对应样本名再查表
ordered_samples <- names(group49355)[hc_49355$order]
dend_colors     <- sci_colors[group49355[ordered_samples]]

pdf("results/figures/QC1_clustering_GSE49355.pdf", width = 14, height = 7)
tryCatch({
  par(mar = c(8, 4, 3, 1), family = "sans")
  plot(hc_49355,
       main  = "Hierarchical Clustering Dendrogram of GSE49355 Samples",
       xlab  = "", sub = "", cex = 0.5, hang = -1,
       col   = dend_colors,   # 叶子颜色按分组
       lwd   = 0.8)

  # 底部分组色带（每个样本位置绘制一个小矩形）
  usr     <- par("usr")
  bar_h   <- 0.02 * diff(usr[3:4])
  for (i in seq_along(hc_49355$order)) {
    rect(i - 0.5, usr[3] - bar_h, i + 0.5, usr[3],
         col = dend_colors[i], border = NA)
  }

  legend("topright",
         legend = c(
           paste0("Normal Colon (n=",   sum(group49355 == "Normal"), ")"),
           paste0("Primary Tumor (n=",  sum(group49355 == "PT"),     ")"),
           paste0("Liver Metastasis (n=",sum(group49355 == "LM"),    ")")
         ),
         fill = sci_colors, border = NA, bty = "n", cex = 0.9)
}, finally = dev.off())

saveRDS(list(hc = hc_49355, distance_matrix = as.matrix(sampleDist_49355)),
        "data/processed/qc_clustering_49355.rds")
cat("  保存: QC1_clustering_GSE49355.pdf\n")

# ------------------------------------------------------------
# QC2: GSE49355 PCA（95%置信椭圆，按三组着色）
# 目的: 评估样本整体分离情况，是否与组织类型对应
# ------------------------------------------------------------
cat(">>> 生成 QC2: GSE49355 PCA...\n")
pca_49355     <- prcomp(t(expr49355_norm[top5k_49355, ]), scale. = TRUE)
var_exp_49355 <- summary(pca_49355)$importance[2, ] * 100
pca_df_49355  <- as.data.frame(pca_49355$x[, 1:2])
pca_df_49355$group <- group49355

p1 <- ggplot(pca_df_49355, aes(PC1, PC2, color = group, shape = group)) +
  geom_point(size = 3, alpha = 0.85, stroke = 0.5) +
  # 【兼容性】stat_ellipse + fill=polygon 在新版ggplot2中分开写以确保边框可见
  stat_ellipse(aes(fill = group), level = 0.95, geom = "polygon",
               alpha = 0.08, color = NA) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.5) +
  scale_color_manual(values = sci_colors,
                     labels = c("LM"     = "Liver Metastasis",
                                "Normal" = "Normal Colon",
                                "PT"     = "Primary Tumor")) +
  scale_fill_manual(values  = sci_colors,
                    labels  = c("LM"     = "Liver Metastasis",
                                "Normal" = "Normal Colon",
                                "PT"     = "Primary Tumor"),
                    guide   = "none") +
  scale_shape_manual(values = sci_shapes,
                     labels = c("LM"     = "Liver Metastasis",
                                "Normal" = "Normal Colon",
                                "PT"     = "Primary Tumor")) +
  labs(
    title    = "Principal Component Analysis of GSE49355",
    subtitle = paste0("n = ", nrow(pca_df_49355),
                      " samples, top 5,000 variable genes, 95% confidence ellipses"),
    x        = paste0("PC1 (", round(var_exp_49355[1], 1), "%)"),
    y        = paste0("PC2 (", round(var_exp_49355[2], 1), "%)"),
    color    = "Tissue Type", shape = "Tissue Type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor    = element_blank(),
    panel.grid.major    = element_line(linewidth = 0.2, colour = "#E0E0E0"),
    plot.title          = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle       = element_text(hjust = 0.5, color = "grey40", size = 10),
    legend.position     = "right",
    legend.title        = element_text(face = "bold", size = 11),
    axis.title          = element_text(size = 11)
  )

ggsave("results/figures/QC2_PCA_GSE49355.pdf", p1, width = 7, height = 5.5)
saveRDS(list(
  pca              = pca_49355,
  variance_explained = var_exp_49355,
  pca_scores       = pca_df_49355,
  top_variable_genes = top5k_49355), "data/processed/qc_pca_49355.rds")
cat("  保存: QC2_PCA_GSE49355.pdf\n")

# ------------------------------------------------------------
# QC3: GSE49355 样本异常值检测（Pearson 相关系数箱线图）
# 方法: 计算样本间两两相关系数 → 每样本均值 → Z-score
#       |Z| > 3 标记为潜在异常值
# 目的: 识别技术错误或样本污染
# ------------------------------------------------------------
cat(">>> 生成 QC3: GSE49355 异常样本检测...\n")
sample_cor_49355 <- cor(expr49355_norm, method = "pearson",
                        use = "pairwise.complete.obs")
mean_cor_49355   <- rowMeans(sample_cor_49355)
z_49355          <- scale(mean_cor_49355)[, 1]   # 取向量（drop matrix维度）
outlier_49355    <- names(z_49355[abs(z_49355) > 3])

cat("GSE49355 潜在异常样本（|Z| > 3）:",
    ifelse(length(outlier_49355) == 0, "无", paste(outlier_49355, collapse = ", ")), "\n")

outlier_df_49355 <- data.frame(
  sample     = names(mean_cor_49355),
  mean_cor   = mean_cor_49355,
  z_score    = z_49355,
  is_outlier = ifelse(names(mean_cor_49355) %in% outlier_49355, "Outlier", "Normal"),
  group      = group49355[names(mean_cor_49355)],
  stringsAsFactors = FALSE
)

p3 <- ggplot(outlier_df_49355, aes(x = group, y = mean_cor, fill = group)) +
  geom_boxplot(outlier.size = 1, outlier.shape = 21, outlier.stroke = 0.5,
               alpha = 0.7, linewidth = 0.4, coef = 1.5) +
  geom_jitter(aes(color = group, shape = is_outlier), width = 0.15,
              size = 2.5, alpha = 0.8) +
  geom_hline(yintercept = mean(mean_cor_49355), linetype = "dashed",
             color = "#333333", linewidth = 0.5) +
  scale_fill_manual(values  = sci_colors) +
  scale_color_manual(values = sci_colors) +
  scale_shape_manual(values = c("Normal" = 16, "Outlier" = 17)) +
  scale_x_discrete(labels = c(
    "LM"     = paste0("Liver Metastasis\n(n=",   sum(group49355 == "LM"),     ")"),
    "Normal" = paste0("Normal Colon\n(n=",        sum(group49355 == "Normal"), ")"),
    "PT"     = paste0("Primary Tumor\n(n=",       sum(group49355 == "PT"),     ")")
  )) +
  labs(
    title    = "Sample Quality Assessment of GSE49355",
    subtitle = "Inter-sample Pearson correlation by tissue type (|Z| > 3 flagged as outlier)",
    x        = "Tissue Type",
    y        = "Mean Pearson Correlation Coefficient",
    fill     = "Tissue Type", color = "Tissue Type", shape = "Sample Status"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor      = element_blank(),
    panel.grid.major.x    = element_blank(),
    panel.grid.major      = element_line(linewidth = 0.2, colour = "#E0E0E0"),
    plot.title            = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle         = element_text(hjust = 0.5, color = "grey40", size = 10),
    legend.position       = "right",
    legend.title          = element_text(face = "bold", size = 11),
    axis.title            = element_text(size = 11)
  )

ggsave("results/figures/QC3_outlier_GSE49355.pdf", p3, width = 8, height = 5.5)
saveRDS(list(
  correlation_matrix = sample_cor_49355,
  mean_correlation   = mean_cor_49355,
  z_scores           = z_49355,
  outlier_samples    = outlier_49355,
  outlier_summary    = outlier_df_49355), "data/processed/qc_outlier_49355.rds")
cat("  保存: QC3_outlier_GSE49355.pdf\n")

# ------------------------------------------------------------
# QC4: GSE49355 表达密度图（标准化前后对比）
# 目的: 展示分位数标准化后各样本分布一致性
# 注: 使用 _pre（log2后、QN前）和 _norm（QN后）对比
# ------------------------------------------------------------
cat(">>> 生成 QC4: GSE49355 表达密度图...\n")

# 基于 _pre 版本重新计算高变异基因（与 _norm 用不同基因集，防止循环依赖）
rv_pre  <- apply(expr49355_pre, 1, var, na.rm = TRUE)
top5k_d <- names(sort(rv_pre, decreasing = TRUE))[1:5000]

dens_pre_list  <- apply(expr49355_pre[top5k_d, ],  2, density, na.rm = TRUE)
dens_norm_list <- apply(expr49355_norm[top5k_d, ], 2, density, na.rm = TRUE)

make_dens_df <- function(dens_list, label) {
  do.call(rbind, lapply(seq_along(dens_list), function(i) {
    data.frame(x = dens_list[[i]]$x, y = dens_list[[i]]$y,
               Sample = i, Stage = label)
  }))
}

dens_all_df <- rbind(
  make_dens_df(dens_pre_list,  "Before Normalization"),
  make_dens_df(dens_norm_list, "After Normalization")
)
dens_all_df$Stage <- factor(dens_all_df$Stage,
                            levels = c("Before Normalization", "After Normalization"))

p_dens <- ggplot(dens_all_df, aes(x = x, y = y, group = Sample)) +
  geom_line(alpha = 0.15, linewidth = 0.4, color = "#0072B2") +
  facet_wrap(~Stage, ncol = 2) +
  labs(
    title    = "Expression Distribution of GSE49355 Samples",
    subtitle = "Top 5,000 most variable genes (log2 scale)",
    x        = "Gene Expression (log2)",
    y        = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, colour = "#E0E0E0"),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 10),
    strip.text       = element_text(face = "bold", size = 11),
    axis.title       = element_text(size = 11)
  )

ggsave("results/figures/QC4_density_GSE49355.pdf", p_dens, width = 12, height = 5)
cat("  保存: QC4_density_GSE49355.pdf\n")

# ------------------------------------------------------------
# QC5: 四个数据集表达箱线图（原始 → log2 → QN 三步对比）
# 目的: 验证标准化步骤有效消除数据集间分布差异
# 注: 随机抽样 500 个基因加速渲染（set.seed(42) 保证可重复）
# ------------------------------------------------------------
cat(">>> 生成 QC5: 表达箱线图（三步对比）...\n")

expr_list_raw  <- list(
  "GSE49355"  = expr49355_probe,
  "GSE159216" = expr159216_probe,
  "GSE39582"  = expr39582_probe,
  "GSE72718"  = expr72718_probe
)
expr_list_pre  <- list(
  "GSE49355"  = expr49355_pre,
  "GSE159216" = expr159216_pre,
  "GSE39582"  = expr39582_pre,
  "GSE72718"  = expr72718_pre
)
expr_list_norm <- list(
  "GSE49355"  = expr49355_norm,
  "GSE159216" = expr159216_norm,
  "GSE39582"  = expr39582_norm,
  "GSE72718"  = expr72718_norm
)

make_box_df <- function(expr_list, stage_label, max_genes = 500) {
  do.call(rbind, lapply(names(expr_list), function(ds) {
    m <- expr_list[[ds]]
    if (nrow(m) > max_genes) {
      set.seed(42)
      m <- m[sample(nrow(m), max_genes), ]
    }
    df <- as.data.frame(t(m))
    df$Dataset <- ds
    df$Stage   <- stage_label
    tidyr::pivot_longer(df, cols = -c(Dataset, Stage),
                        names_to = "Gene", values_to = "Expression")
  }))
}

box_all_df <- rbind(
  make_box_df(expr_list_raw,  "Raw (Original)"),
  make_box_df(expr_list_pre,  "log2 Transformed"),
  make_box_df(expr_list_norm, "Quantile Normalized")
)
box_all_df$Stage <- factor(box_all_df$Stage,
                           levels = c("Raw (Original)", "log2 Transformed",
                                      "Quantile Normalized"))

p_box <- ggplot(box_all_df, aes(x = Dataset, y = Expression, fill = Dataset)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2, alpha = 0.7,
               linewidth = 0.4, coef = 1.5) +
  facet_wrap(~Stage, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(
    "GSE49355"  = "#E69F00",
    "GSE159216" = "#56B4E9",
    "GSE39582"  = "#009E73",
    "GSE72718"  = "#F0E442"
  )) +
  labs(
    title    = "Expression Distribution Across Processing Steps",
    subtitle = "Raw probe intensities, log2 transformation, and quantile normalization (500 genes sampled)",
    x        = "", y = "Expression Value", fill = "Dataset"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor    = element_blank(),
    panel.grid.major.x  = element_blank(),
    panel.grid.major    = element_line(linewidth = 0.2, colour = "#E0E0E0"),
    plot.title          = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle       = element_text(hjust = 0.5, color = "grey40", size = 10),
    strip.text          = element_text(face = "bold", size = 11),
    axis.text.x         = element_text(angle = 30, hjust = 1, size = 10),
    axis.title          = element_text(size = 11),
    legend.position     = "none"
  )

ggsave("results/figures/QC5_boxplot_allDatasets.pdf", p_box, width = 14, height = 5)
cat("  保存: QC5_boxplot_allDatasets.pdf\n")

# ------------------------------------------------------------
# QC6: GSE49355 样本相关性热图（按组排序 + 组间分隔线）
# 目的: 展示组内相关性高于组间，确认数据分组正确
# ------------------------------------------------------------
cat(">>> 生成 QC6: GSE49355 样本相关性热图...\n")
sample_cor <- cor(expr49355_norm[top5k_49355, ], method = "pearson",
                  use = "pairwise.complete.obs")

# 按组排序: Normal → PT → LM，组内按组内均值相关性排列
order_idx <- integer(0)
for (g in c("Normal", "PT", "LM")) {
  g_idx <- which(group49355 == g)
  if (length(g_idx) > 1) {
    g_sub_cor <- sample_cor[g_idx, g_idx]
    g_mean    <- rowMeans(g_sub_cor)
    order_idx <- c(order_idx, g_idx[order(g_mean, decreasing = TRUE)])
  } else {
    order_idx <- c(order_idx, g_idx)
  }
}
sample_cor_ordered <- sample_cor[order_idx, order_idx]

ann_col <- data.frame(
  `Tissue Type` = factor(
    group49355[colnames(sample_cor_ordered)],
    levels = c("Normal", "PT", "LM"),
    labels = c("Normal Colon", "Primary Tumor", "Liver Metastasis")
  ),
  row.names = colnames(sample_cor_ordered),
  check.names = FALSE
)

ann_colors <- list(
  `Tissue Type` = c(
    "Normal Colon"    = "#0072B2",
    "Primary Tumor"   = "#D55E00",
    "Liver Metastasis"= "#CC79A7"
  )
)

# 组间分隔线位置（按排序后的 Normal/PT/LM 数量累计）
n_normal <- sum(group49355 == "Normal")
n_pt     <- sum(group49355 == "PT")
gap_rows <- c(n_normal, n_normal + n_pt)

pdf("results/figures/QC6_correlation_heatmap_GSE49355.pdf", width = 10, height = 9)
tryCatch({
  pheatmap::pheatmap(
    sample_cor_ordered,
    annotation_col      = ann_col,
    annotation_row      = ann_col,
    annotation_colors   = ann_colors,
    color               = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
    border_color        = NA,
    show_rownames       = FALSE,
    show_colnames       = FALSE,
    main                = "Inter-sample Pearson Correlation Heatmap, GSE49355 (Top 5,000 Variable Genes)",
    fontsize_main       = 13,
    annotation_legend   = TRUE,
    legend              = TRUE,
    annotation_names_row = TRUE,
    annotation_names_col = TRUE,
    padding             = 2,
    cluster_rows        = FALSE,
    cluster_cols        = FALSE,
    gaps_row            = gap_rows,
    gaps_col            = gap_rows,
    fontsize            = 10
  )
}, finally = dev.off())

saveRDS(list(
  correlation_matrix = sample_cor,
  annotation         = ann_col,
  top_variable_genes = top5k_49355), "data/processed/qc_correlation_49355.rds")
cat("  保存: QC6_correlation_heatmap_GSE49355.pdf\n")

# ------------------------------------------------------------
# QC7: 分位数标准化前后 PCA 对比
# 目的: 直观展示 QN 对消除技术变异、保留生物学变异的效果
# ------------------------------------------------------------
cat(">>> 生成 QC7: 标准化前后 PCA 对比...\n")
rv_c  <- apply(expr49355_pre, 1, var, na.rm = TRUE)
top3k <- names(sort(rv_c, decreasing = TRUE))[1:3000]

pca_bef     <- prcomp(t(expr49355_pre[top3k,  ]), scale. = TRUE)
pca_aft     <- prcomp(t(expr49355_norm[top3k, ]), scale. = TRUE)
var_exp_bef <- summary(pca_bef)$importance[2, 1:2] * 100
var_exp_aft <- summary(pca_aft)$importance[2, 1:2] * 100

pca_bef_df <- data.frame(
  PC1   = pca_bef$x[, 1], PC2 = pca_bef$x[, 2],
  group = group49355, time = "Before Normalization"
)
pca_aft_df <- data.frame(
  PC1   = pca_aft$x[, 1], PC2 = pca_aft$x[, 2],
  group = group49355, time = "After Normalization"
)

pca_norm_compare       <- rbind(pca_bef_df, pca_aft_df)
pca_norm_compare$time  <- factor(pca_norm_compare$time,
                                 levels = c("Before Normalization", "After Normalization"))

# 每个 facet 的方差解释比注释标签
var_exp_labels <- data.frame(
  time  = c("Before Normalization", "After Normalization"),
  label = c(
    paste0("PC1: ", round(var_exp_bef[1], 1), "% | PC2: ", round(var_exp_bef[2], 1), "%"),
    paste0("PC1: ", round(var_exp_aft[1], 1), "% | PC2: ", round(var_exp_aft[2], 1), "%")
  ),
  stringsAsFactors = FALSE
)

p5 <- ggplot(pca_norm_compare, aes(PC1, PC2, color = group, shape = group)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.5) +
  facet_wrap(~time, ncol = 2) +
  scale_color_manual(values = sci_colors,
                     labels = c("LM"     = "Liver Metastasis",
                                "Normal" = "Normal Colon",
                                "PT"     = "Primary Tumor")) +
  scale_shape_manual(values = sci_shapes,
                     labels = c("LM"     = "Liver Metastasis",
                                "Normal" = "Normal Colon",
                                "PT"     = "Primary Tumor")) +
  geom_text(data = var_exp_labels,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.5, size = 3, color = "grey30",
            inherit.aes = FALSE) +
  labs(
    title    = "Effect of Quantile Normalization on Sample Separation",
    subtitle = "GSE49355, top 3,000 variable genes",
    x        = "PC1", y = "PC2",
    color    = "Tissue Type", shape = "Tissue Type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, colour = "#E0E0E0"),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 10),
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold", size = 11)
  )

ggsave("results/figures/QC7_normalization_effect_GSE49355.pdf", p5, width = 12, height = 5.5)
saveRDS(list(
  pca_before     = pca_bef,
  pca_after      = pca_aft,
  var_exp_before = var_exp_bef,
  var_exp_after  = var_exp_aft,
  top_genes      = top3k), "data/processed/qc_normalization_compare.rds")
cat("  保存: QC7_normalization_effect_GSE49355.pdf\n")

# ------------------------------------------------------------
# QC8: 探针-基因映射汇总柱状图
# 目的: 展示各数据集从原始探针到基因symbol的转换数量
# ------------------------------------------------------------
cat(">>> 生成 QC8: 探针-基因映射汇总...\n")
probe_count <- c(
  "GSE49355"  = nrow(expr49355_probe),
  "GSE159216" = nrow(expr159216_probe),
  "GSE39582"  = nrow(expr39582_probe),
  "GSE72718"  = nrow(expr72718_probe)
)
gene_count <- c(
  "GSE49355"  = nrow(expr49355_norm),
  "GSE159216" = nrow(expr159216_norm),
  "GSE39582"  = nrow(expr39582_norm),
  "GSE72718"  = nrow(expr72718_norm)
)

sum_df <- data.frame(
  Dataset      = factor(names(probe_count), levels = names(probe_count)),
  Raw_Probes   = unname(probe_count),
  Mapped_Genes = unname(gene_count[match(names(probe_count), names(gene_count))]),
  row.names    = NULL
)

sum_long <- tidyr::pivot_longer(sum_df, cols = -Dataset,
                                names_to = "Type", values_to = "Count")

p_count <- ggplot(sum_long, aes(x = Dataset, y = Count / 1000, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.85) +
  geom_text(
    aes(label = ifelse(Count >= 1000,
                       paste0(round(Count / 1000, 1), "K"),
                       as.character(Count))),
    stat = "identity", position = position_dodge(width = 0.7),
    vjust = -0.5, size = 3, color = "#333333"
  ) +
  scale_fill_manual(values = c("Raw_Probes" = "#4393C3", "Mapped_Genes" = "#D6604D"),
                    labels = c("Raw Probes", "Mapped Genes")) +
  labs(
    title    = "Probe-to-Gene Mapping Summary",
    subtitle = "Number of raw probes and annotated unique genes per dataset",
    x        = "", y = "Count (\u00D7 1,000)", fill = ""
  ) +
  coord_cartesian(ylim = c(0, max(sum_long$Count / 1000) * 1.15)) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major   = element_line(linewidth = 0.2, colour = "#E0E0E0"),
    plot.title         = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle      = element_text(hjust = 0.5, color = "grey40", size = 10),
    legend.position    = "top",
    legend.title       = element_text(face = "bold", size = 11),
    axis.text.x        = element_text(size = 10),
    axis.title         = element_text(size = 11)
  )

ggsave("results/figures/QC8_gene_counts.pdf", p_count, width = 7, height = 5)
cat("  保存: QC8_gene_counts.pdf\n")

# ------------------------------------------------------------
# QC9: GSE72718 PCA（按三类组织类型着色）
# 目的: 验证肝转移灶与原发灶的基因表达差异
# tissue type:ch1 列取值: "colorectal cancer liver metastasis" /
#                          "primary tumor with liver metastasis" /
#                          "primary tumor without liver metastasis"
# ------------------------------------------------------------
cat(">>> 生成 QC9: GSE72718 PCA...\n")
rv_72718 <- apply(expr72718_norm, 1, var, na.rm = TRUE)
# 过滤方差为0或NA的基因（避免 prcomp 报错）
expr72718_norm_pca <- expr72718_norm[rv_72718 > 0 & !is.na(rv_72718), ]

pca_72718     <- prcomp(t(expr72718_norm_pca), scale. = TRUE)
var_exp_72718 <- summary(pca_72718)$importance[2, ] * 100
pca_df_72718  <- as.data.frame(pca_72718$x[, 1:2])

# tissue_72718_raw 是从 pData 取出的列向量，其 names 为行名（即样本名）
# rownames(pca_df_72718) 是样本名，可直接作为字符索引
tissue_72718_raw <- pheno72718[["tissue type:ch1"]]
names(tissue_72718_raw) <- rownames(pheno72718)  # 确保 names 与样本名一致

pca_df_72718$tissue <- sapply(rownames(pca_df_72718), function(s) {
  val <- as.character(tissue_72718_raw[s])
  if (is.na(val) || val == "") return("Unknown")
  if (grepl("without liver metastasis", val, ignore.case = TRUE)) {
    return("Primary without LM")
  } else if (grepl("colorectal cancer liver metastasis", val, ignore.case = TRUE)) {
    return("Liver Metastasis")
  } else {
    return("Primary with LM")
  }
})

pca_df_72718$tissue <- factor(
  pca_df_72718$tissue,
  levels = c("Liver Metastasis", "Primary with LM", "Primary without LM")
)

sci_72718_colors <- c("Liver Metastasis"    = "#CC79A7",
                      "Primary with LM"     = "#D55E00",
                      "Primary without LM"  = "#0072B2")
sci_72718_shapes <- c("Liver Metastasis"    = 15,
                      "Primary with LM"     = 17,
                      "Primary without LM"  = 16)

p72718 <- ggplot(pca_df_72718, aes(PC1, PC2, color = tissue, shape = tissue)) +
  geom_point(size = 3, alpha = 0.8, stroke = 0.5) +
  stat_ellipse(aes(fill = tissue), level = 0.95, geom = "polygon",
               alpha = 0.08, color = NA) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.5) +
  scale_color_manual(values = sci_72718_colors) +
  scale_fill_manual(values  = sci_72718_colors, guide = "none") +
  scale_shape_manual(values = sci_72718_shapes) +
  labs(
    title    = "Principal Component Analysis of GSE72718",
    subtitle = paste0("n = ", nrow(pca_df_72718),
                      " samples, all variable genes, 95% confidence ellipses"),
    x        = paste0("PC1 (", round(var_exp_72718[1], 1), "%)"),
    y        = paste0("PC2 (", round(var_exp_72718[2], 1), "%)"),
    color    = "Tissue Type", shape = "Tissue Type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, colour = "#E0E0E0"),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 10),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 11),
    axis.title       = element_text(size = 11)
  )

ggsave("results/figures/QC9_PCA_GSE72718.pdf", p72718, width = 7, height = 5.5)
saveRDS(list(pca = pca_72718, variance_explained = var_exp_72718,
             pca_scores = pca_df_72718),
        "data/processed/qc_pca_72718.rds")
cat("  保存: QC9_PCA_GSE72718.pdf\n")

# ------------------------------------------------------------
# QC10: GSE39582 PCA（按样本类型: 肿瘤 vs 正常）
# ------------------------------------------------------------
cat(">>> 生成 QC10: GSE39582 PCA...\n")
rv_39582    <- apply(expr39582_norm, 1, var, na.rm = TRUE)
top5k_39582 <- names(sort(rv_39582, decreasing = TRUE))[1:5000]

pca_39582     <- prcomp(t(expr39582_norm[top5k_39582, ]), scale. = TRUE)
var_exp_39582 <- summary(pca_39582)$importance[2, ] * 100
pca_df_39582  <- as.data.frame(pca_39582$x[, 1:2])

pca_df_39582$sample_type <- ds39582[rownames(pca_df_39582)]
pca_df_39582$sample_type <- ifelse(pca_df_39582$sample_type == "Normal",
                                   "Normal Tissue", "Tumor Tissue")

sci_39582_colors <- c("Normal Tissue" = "#0072B2", "Tumor Tissue" = "#D55E00")
sci_39582_shapes <- c("Normal Tissue" = 16, "Tumor Tissue" = 17)

p2 <- ggplot(pca_df_39582, aes(PC1, PC2, color = sample_type, shape = sample_type)) +
  geom_point(size = 3, alpha = 0.8, stroke = 0.5) +
  stat_ellipse(aes(fill = sample_type), level = 0.95, geom = "polygon",
               alpha = 0.08, color = NA) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.5) +
  scale_color_manual(values = sci_39582_colors) +
  scale_fill_manual(values  = sci_39582_colors, guide = "none") +
  scale_shape_manual(values = sci_39582_shapes) +
  labs(
    title    = "Principal Component Analysis of GSE39582",
    subtitle = paste0("n = ", nrow(pca_df_39582),
                      " samples, top 5,000 variable genes, 95% confidence ellipses"),
    x        = paste0("PC1 (", round(var_exp_39582[1], 1), "%)"),
    y        = paste0("PC2 (", round(var_exp_39582[2], 1), "%)"),
    color    = "Sample Type", shape = "Sample Type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, colour = "#E0E0E0"),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 10),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 11),
    axis.title       = element_text(size = 11)
  )

ggsave("results/figures/QC10_PCA_GSE39582.pdf", p2, width = 7, height = 5.5)
saveRDS(list(pca = pca_39582, variance_explained = var_exp_39582,
             pca_scores = pca_df_39582,
             top_variable_genes = top5k_39582),
        "data/processed/qc_pca_39582.rds")
cat("  保存: QC10_PCA_GSE39582.pdf\n")

# ------------------------------------------------------------
# QC11: GSE159216 PCA（按样本类型着色）
# 注: pheno159216 中通过 dataset:ch1 列判断 Normal/Tumor
# ------------------------------------------------------------
cat(">>> 生成 QC11: GSE159216 PCA...\n")
rv_159216    <- apply(expr159216_norm, 1, var, na.rm = TRUE)
top5k_159216 <- names(sort(rv_159216, decreasing = TRUE))[1:5000]

pca_159216     <- prcomp(t(expr159216_norm[top5k_159216, ]), scale. = TRUE)
var_exp_159216 <- summary(pca_159216)$importance[2, ] * 100
pca_df_159216  <- as.data.frame(pca_159216$x[, 1:2])

# GSE159216 样本类型判断
# 【核查点】若 dataset:ch1 列不存在或内容与预期不符，请打印 head(pheno159216) 检查
ds159216_raw <- pheno159216[["dataset:ch1"]]
if (!is.null(ds159216_raw) && any(!is.na(ds159216_raw))) {
  ds159216 <- ifelse(
    grepl("Normal|non.tumor|control", ds159216_raw, ignore.case = TRUE),
    "Normal Tissue", "Tumor"
  )
} else {
  cat("  警告: GSE159216 dataset:ch1 列不存在或全为NA，默认全部归为 Tumor\n")
  ds159216 <- rep("Tumor", ncol(expr159216_norm))
}
names(ds159216) <- colnames(expr159216_norm)

pca_df_159216$sample_type <- ds159216[rownames(pca_df_159216)]

sci_159216_colors <- c("Normal Tissue" = "#0072B2", "Tumor" = "#D55E00")
sci_159216_shapes <- c("Normal Tissue" = 16, "Tumor" = 17)

p159216 <- ggplot(pca_df_159216, aes(PC1, PC2, color = sample_type, shape = sample_type)) +
  geom_point(size = 3, alpha = 0.8, stroke = 0.5) +
  stat_ellipse(aes(fill = sample_type), level = 0.95, geom = "polygon",
               alpha = 0.08, color = NA) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.5) +
  scale_color_manual(values = sci_159216_colors) +
  scale_fill_manual(values  = sci_159216_colors, guide = "none") +
  scale_shape_manual(values = sci_159216_shapes) +
  labs(
    title    = "Principal Component Analysis of GSE159216",
    subtitle = paste0("n = ", nrow(pca_df_159216),
                      " samples, top 5,000 variable genes, 95% confidence ellipses"),
    x        = paste0("PC1 (", round(var_exp_159216[1], 1), "%)"),
    y        = paste0("PC2 (", round(var_exp_159216[2], 1), "%)"),
    color    = "Sample Type", shape = "Sample Type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, colour = "#E0E0E0"),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 10),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 11),
    axis.title       = element_text(size = 11)
  )

ggsave("results/figures/QC11_PCA_GSE159216.pdf", p159216, width = 7, height = 5.5)
saveRDS(list(pca = pca_159216, variance_explained = var_exp_159216,
             pca_scores = pca_df_159216,
             top_variable_genes = top5k_159216),
        "data/processed/qc_pca_159216.rds")
cat("  保存: QC11_PCA_GSE159216.pdf\n")

# ------------------------------------------------------------
# QC12: GSE39582 样本相关性热图（带 TNM / 分子亚型临床注释）
# 目的: 评估样本质量，揭示与分子亚型的相关性模式
# ------------------------------------------------------------
cat(">>> 生成 QC12: GSE39582 相关性热图...\n")
rv_39582_hm    <- apply(expr39582_norm, 1, var, na.rm = TRUE)
top5k_39582_hm <- names(sort(rv_39582_hm, decreasing = TRUE))[1:5000]
sample_cor_39582 <- cor(expr39582_norm[top5k_39582_hm, ], method = "pearson",
                        use = "pairwise.complete.obs")

ann_col_39582 <- data.frame(
  `Sample Type`  = factor(ds39582[colnames(sample_cor_39582)],
                          levels = c("Normal", "Tumor"),
                          labels = c("Normal Tissue", "Tumor Tissue")),
  `TNM Stage`    = tnm39582[colnames(sample_cor_39582)],
  `Mol. Subtype` = ms39582[colnames(sample_cor_39582)],
  row.names = colnames(sample_cor_39582),
  check.names = FALSE
)

ann_colors_39582 <- list(
  `Sample Type`  = c("Normal Tissue" = "#0072B2", "Tumor Tissue" = "#D55E00"),
  `TNM Stage`    = c("0" = "#F7F7F7", "1" = "#4393C3", "2" = "#92C5DE",
                     "3" = "#D6604D", "4" = "#B2182B"),
  `Mol. Subtype` = c("C1" = "#E41A1C", "C2" = "#377EB8", "C3" = "#4DAF4A",
                     "C4" = "#984EA3", "C5" = "#FF7F00", "C6" = "#A65628")
)

pdf("results/figures/QC12_correlation_heatmap_GSE39582.pdf", width = 10, height = 9)
tryCatch({
  pheatmap::pheatmap(
    sample_cor_39582,
    annotation_col      = ann_col_39582,
    annotation_row      = ann_col_39582,
    annotation_colors   = ann_colors_39582,
    color               = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
    border_color        = NA,
    show_rownames       = FALSE,
    show_colnames       = FALSE,
    main                = "Inter-sample Pearson Correlation Heatmap, GSE39582 (Top 5,000 Variable Genes)",
    fontsize_main       = 13,
    legend              = TRUE,
    clustering_method   = "average",
    annotation_names_row = TRUE,
    annotation_names_col = TRUE,
    fontsize            = 10
  )
}, finally = dev.off())

saveRDS(list(
  correlation_matrix  = sample_cor_39582,
  top_variable_genes  = top5k_39582_hm,
  clinical_annotation = ann_col_39582,
  ds39582             = ds39582,
  tnm39582            = tnm39582,
  ms39582             = ms39582,
  os_event39582       = os_event39582), "data/processed/qc_correlation_39582.rds")
cat("  保存: QC12_correlation_heatmap_GSE39582.pdf\n")

# ------------------------------------------------------------
# QC13: GSE49355 Top 50 高变异基因热图（Z-score 标准化，按组排序）
# 目的: 直观展示三组间表达差异最显著的基因
# ------------------------------------------------------------
cat(">>> 生成 QC13: Top 50 高变异基因热图...\n")
top50_49355 <- names(sort(apply(expr49355_norm, 1, var, na.rm = TRUE),
                          decreasing = TRUE))[1:50]
mat_top50   <- expr49355_norm[top50_49355, ]

# 对每个基因进行行方向 Z-score 标准化（便于颜色对比）
mat_top50_z <- t(scale(t(mat_top50)))

ann_col_top50 <- data.frame(
  `Tissue Type` = factor(
    group49355[colnames(mat_top50_z)],
    levels = c("Normal", "PT", "LM"),
    labels = c("Normal Colon", "Primary Tumor", "Liver Metastasis")
  ),
  row.names = colnames(mat_top50_z),
  check.names = FALSE
)

ann_colors_top50 <- list(
  `Tissue Type` = c(
    "Normal Colon"    = "#0072B2",
    "Primary Tumor"   = "#D55E00",
    "Liver Metastasis"= "#CC79A7"
  )
)

# 按组排序: Normal → PT → LM
hm_order <- unlist(lapply(c("Normal", "PT", "LM"), function(g) {
  which(group49355[colnames(mat_top50_z)] == g)
}))
mat_top50_z   <- mat_top50_z[, hm_order]
ann_col_top50 <- ann_col_top50[colnames(mat_top50_z), , drop = FALSE]
gap_top50     <- cumsum(table(group49355[colnames(mat_top50_z)])[c("Normal", "PT")])

pdf("results/figures/QC13_top50_variable_genes_heatmap.pdf", width = 8, height = 10)
tryCatch({
  pheatmap::pheatmap(
    mat_top50_z,
    annotation_col    = ann_col_top50,
    annotation_colors = ann_colors_top50,
    color             = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
    border_color      = NA,
    show_rownames     = TRUE,
    show_colnames     = FALSE,
    main              = "Top 50 Most Variable Genes in GSE49355",
    fontsize_main     = 13,
    fontsize_row      = 7,
    annotation_legend = TRUE,
    legend            = TRUE,
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    gaps_col          = gap_top50,
    clustering_method = "average",
    scale             = "none"
  )
}, finally = dev.off())
cat("  保存: QC13_top50_variable_genes_heatmap.pdf\n")

# ------------------------------------------------------------
# QC14: GSE49355 差异表达火山图（LM vs Normal + PT vs Normal）
# 方法: limma lmFit → eBayes → topTable (BH校正)
# 显著性: FDR < 0.05 且 |log2FC| > 1
# 目的: 初步了解三组间的整体差异表达概况
# ------------------------------------------------------------
cat(">>> 生成 QC14: 差异表达火山图...\n")

# 构建 design matrix（~ 0 + factor 格式，各列对应一组）
# 【核查点】factor 水平顺序决定列名顺序，需与下一行 colnames 一致
design_49355       <- model.matrix(~ 0 + factor(group49355,
                                                levels = c("LM", "Normal", "PT")))
colnames(design_49355) <- c("LM", "Normal", "PT")
fit_49355 <- lmFit(expr49355_norm, design_49355)

# LM vs Normal
contrast_lm <- makeContrasts(LM_vs_Normal = LM - Normal, levels = design_49355)
fit_lm      <- eBayes(contrasts.fit(fit_49355, contrast_lm))
res_lm      <- topTable(fit_lm, number = Inf, adjust.method = "BH")

# PT vs Normal
contrast_pt <- makeContrasts(PT_vs_Normal = PT - Normal, levels = design_49355)
fit_pt      <- eBayes(contrasts.fit(fit_49355, contrast_pt))
res_pt      <- topTable(fit_pt, number = Inf, adjust.method = "BH")

volcano_all <- rbind(
  data.frame(logFC      = res_lm$logFC,
             negLog10P  = -log10(res_lm$P.Value),
             FDR        = res_lm$adj.P.Val,
             Comparison = "LM vs Normal",
             Gene       = rownames(res_lm)),
  data.frame(logFC      = res_pt$logFC,
             negLog10P  = -log10(res_pt$P.Value),
             FDR        = res_pt$adj.P.Val,
             Comparison = "PT vs Normal",
             Gene       = rownames(res_pt))
)
volcano_all$Significance <- ifelse(
  volcano_all$FDR < 0.05 & abs(volcano_all$logFC) > 1,
  "Significant", "Not Significant"
)
volcano_all$Comparison <- factor(volcano_all$Comparison,
                                 levels = c("LM vs Normal", "PT vs Normal"))

p_volcano <- ggplot(volcano_all, aes(x = logFC, y = negLog10P, color = Significance)) +
  geom_point(alpha = 0.4, size = 1.2) +
  scale_color_manual(values = c("Significant"     = "#D55E00",
                                "Not Significant" = "#999999")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "#333333", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "#333333", linewidth = 0.3) +
  facet_wrap(~Comparison, ncol = 2) +
  labs(
    title    = "Differential Expression Analysis of GSE49355",
    subtitle = "LM vs Normal and PT vs Normal (|log2FC| > 1 and FDR < 0.05)",
    x        = expression(log[2]~"Fold Change"),
    y        = expression(-log[10]("P value")),
    color    = "Significance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, colour = "#E0E0E0"),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 10),
    strip.text       = element_text(face = "bold", size = 11),
    legend.position  = "top",
    legend.title     = element_text(face = "bold", size = 11),
    axis.title       = element_text(size = 11)
  )

ggsave("results/figures/QC14_volcano_GSE49355.pdf", p_volcano, width = 10, height = 5)
saveRDS(list(
  lm_vs_normal = res_lm,
  pt_vs_normal = res_pt,
  volcano_data = volcano_all), "data/processed/qc_differential_expression.rds")
cat("  保存: QC14_volcano_GSE49355.pdf\n")

# ------------------------------------------------------------
# QC15: 各数据集基因完整性柱状图
# 仅统计在所有样本中均无缺失值的基因（MissingRate == 0）
# 目的: 展示各数据集的数据完整性，排查批量缺失问题
# ------------------------------------------------------------
cat(">>> 生成 QC15: 数据完整性柱状图...\n")

calc_missing <- function(expr_mat, ds_name) {
  missing_rate <- rowMeans(is.na(expr_mat))
  data.frame(
    Gene        = rownames(expr_mat),
    MissingRate = missing_rate,
    Dataset     = ds_name,
    stringsAsFactors = FALSE
  )
}

missing_all <- rbind(
  calc_missing(expr49355_norm,  "GSE49355"),
  calc_missing(expr159216_norm, "GSE159216"),
  calc_missing(expr39582_norm,  "GSE39582"),
  calc_missing(expr72718_norm,  "GSE72718")
)

# 只保留完整基因（MissingRate == 0）用于绘图
complete_genes_df <- missing_all[missing_all$MissingRate == 0, ]

p_missing <- ggplot(complete_genes_df, aes(x = Dataset, fill = Dataset)) +
  geom_bar() +
  scale_fill_manual(values = c("GSE49355"  = "#E69F00",
                               "GSE159216" = "#56B4E9",
                               "GSE39582"  = "#009E73",
                               "GSE72718"  = "#F0E442")) +
  geom_text(aes(label = scales::comma(after_stat(count))),
            stat = "count", vjust = -0.3, size = 3.5, color = "#333333") +
  coord_cartesian(ylim = c(0, max(table(complete_genes_df$Dataset)) * 1.12)) +
  labs(
    title    = "Gene Expression Data Completeness",
    subtitle = "Number of genes with no missing values across all samples per dataset",
    x        = "", y = "Number of Complete Genes (MissingRate = 0)",
    fill     = "Dataset"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major   = element_line(linewidth = 0.2, colour = "#E0E0E0"),
    plot.title         = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle      = element_text(hjust = 0.5, color = "grey40", size = 10),
    legend.position    = "none",
    axis.title         = element_text(size = 11),
    axis.text.x        = element_text(size = 10)
  )

ggsave("results/figures/QC15_data_completeness.pdf", p_missing, width = 6, height = 5)
cat("  保存: QC15_data_completeness.pdf\n")

# ------------------------------------------------------------
# QC16: 数据集间基因重叠柱状图（两两 pairwise overlap）
# 替代4集 Venn（Venn 超过3个集合可读性差）
# 目的: 了解哪些基因跨数据集通用，为后续合并分析奠定基础
# ------------------------------------------------------------
cat(">>> 生成 QC16: 数据集间基因重叠图...\n")

gene_sets <- list(
  GSE49355  = rownames(expr49355_norm),
  GSE159216 = rownames(expr159216_norm),
  GSE39582  = rownames(expr39582_norm),
  GSE72718  = rownames(expr72718_norm)
)

# 安全检查
gene_sets <- lapply(gene_sets, function(x) if (is.null(x)) character(0) else x)
cat("  各数据集基因数:\n")
print(sapply(gene_sets, length))
cat("  四数据集共有基因数:",
    length(Reduce(intersect, gene_sets)), "\n")

set_names <- names(gene_sets)
n_sets    <- length(set_names)

upset_rows <- list()
for (i in 1:(n_sets - 1)) {
  for (j in (i + 1):n_sets) {
    pair_name <- paste(set_names[i], set_names[j], sep = " / ")
    overlap   <- length(intersect(gene_sets[[set_names[i]]], gene_sets[[set_names[j]]]))
    upset_rows[[pair_name]] <- overlap
  }
}

upset_df <- data.frame(
  Pair    = factor(names(upset_rows), levels = names(upset_rows)),
  Overlap = unname(unlist(upset_rows)),
  stringsAsFactors = FALSE
)

p_upset <- ggplot(upset_df, aes(x = Pair, y = Overlap, fill = Pair)) +
  geom_bar(stat = "identity", alpha = 0.85, width = 0.7) +
  geom_text(aes(label = scales::comma(Overlap)),
            vjust = -0.3, size = 3.5, color = "#333333") +
  scale_fill_brewer(palette = "Set2") +
  coord_cartesian(ylim = c(0, max(upset_df$Overlap) * 1.12)) +
  labs(
    title    = "Gene Overlap Between Datasets",
    subtitle = "Number of shared genes for each pairwise comparison",
    x        = "", y = "Shared Gene Count", fill = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major   = element_line(linewidth = 0.2, colour = "#E0E0E0"),
    plot.title         = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle      = element_text(hjust = 0.5, color = "grey40", size = 10),
    legend.position    = "none",
    axis.text.x        = element_text(angle = 30, hjust = 1, size = 10),
    axis.title         = element_text(size = 11)
  )

ggsave("results/figures/QC16_gene_overlap.pdf", p_upset, width = 8, height = 5)
saveRDS(gene_sets, "data/processed/qc_gene_sets.rds")
cat("  保存: QC16_gene_overlap.pdf\n")

# ------------------------------------------------------------
# QC17: GSE39582 分子亚型 PCA（CIT C1-C6 着色）
# 目的: 展示结直肠癌 CIT 分子亚型的 PCA 分布
# 注: ms39582 已在 Section 5 中提取，此处直接复用 pca_df_39582
# ------------------------------------------------------------
cat(">>> 生成 QC17: GSE39582 分子亚型 PCA...\n")

pca_df_39582$mol_subtype <- as.character(ms39582[rownames(pca_df_39582)])
pca_df_39582$mol_subtype[is.na(pca_df_39582$mol_subtype)] <- "Unknown"

subtype_levels <- c("C1", "C2", "C3", "C4", "C5", "C6", "Unknown")
subtype_colors <- c("C1" = "#E41A1C", "C2" = "#377EB8", "C3" = "#4DAF4A",
                    "C4" = "#984EA3", "C5" = "#FF7F00", "C6" = "#A65628",
                    "Unknown" = "#BBBBBB")
# shape 向量长度与 subtype_levels 严格对齐，避免循环截断错误
subtype_shapes <- setNames(c(16, 17, 15, 18, 8, 3, 1), subtype_levels)

pca_df_39582$mol_subtype <- factor(pca_df_39582$mol_subtype, levels = subtype_levels)

p_subtype <- ggplot(pca_df_39582, aes(PC1, PC2, color = mol_subtype, shape = mol_subtype)) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = subtype_colors, drop = FALSE) +
  scale_shape_manual(values = subtype_shapes,  drop = FALSE) +
  labs(
    title    = "Molecular Subtype Separation in GSE39582",
    subtitle = "CIT Consensus Molecular Subtype classification (C1-C6)",
    x        = paste0("PC1 (", round(var_exp_39582[1], 1), "%)"),
    y        = paste0("PC2 (", round(var_exp_39582[2], 1), "%)"),
    color    = "CIT Subtype", shape = "CIT Subtype"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, colour = "#E0E0E0"),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 10),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 10),
    axis.title       = element_text(size = 11)
  )

ggsave("results/figures/QC17_GSE39582_subtype_PCA.pdf", p_subtype, width = 7.5, height = 5.5)
cat("  保存: QC17_GSE39582_subtype_PCA.pdf\n")

# ============================================================
# Section 8: 汇总表格输出
# ============================================================

# ------------------------------------------------------------
# Table 1: 数据集概览表
# ------------------------------------------------------------
cat(">>> 生成 Table 1: 数据集概览...\n")
overview_df <- data.frame(
  Dataset      = c("GSE49355", "GSE159216", "GSE39582", "GSE72718"),
  Platform     = c("GPL96", "GPL17586", "GPL570", "GPL17586"),
  Array_Type   = c("HG-U133A", "HTA 2.0", "HG-U133 Plus 2.0", "HTA 2.0"),
  Role         = c("Training", "Validation 1", "Validation 2", "Independent"),
  n_Samples    = c(ncol(expr49355_norm), ncol(expr159216_norm),
                   ncol(expr39582_norm),  ncol(expr72718_norm)),
  Raw_Probes   = c(nrow(expr49355_probe), nrow(expr159216_probe),
                   nrow(expr39582_probe),  nrow(expr72718_probe)),
  Mapped_Genes = c(nrow(expr49355_norm),  nrow(expr159216_norm),
                   nrow(expr39582_norm),   nrow(expr72718_norm)),
  check.names = FALSE
)
write.csv(overview_df, "results/tables/Table1_dataset_overview.csv", row.names = FALSE)
cat("  保存: Table1_dataset_overview.csv\n")
print(overview_df)
saveRDS(list(
  probe_counts = probe_count,
  gene_counts  = gene_count,
  overview     = overview_df), "data/processed/qc_mapping_summary.rds")

# ------------------------------------------------------------
# Table 2: GSE49355 样本分组详情表
# ------------------------------------------------------------
cat(">>> 生成 Table 2: GSE49355 分组详情...\n")
group_detail_df <- data.frame(
  Sample_ID   = colnames(expr49355_norm),
  Title       = pheno49355$title,
  Group       = group49355,
  Group_Label = dplyr::recode(group49355,
                              "Normal" = "Normal Colon",
                              "PT"     = "Primary Tumor",
                              "LM"     = "Liver Metastasis"),
  stringsAsFactors = FALSE
)
write.csv(group_detail_df, "results/tables/Table2_GSE49355_group_detail.csv",
          row.names = FALSE)
cat("  保存: Table2_GSE49355_group_detail.csv\n")

# ------------------------------------------------------------
# Table 3: GSE39582 临床特征汇总表
# tnm_table["NA"] 用字符串"NA"索引在有些情况下会返回 NA，
#         改用 sum(is.na(...)) 直接统计
# ------------------------------------------------------------
cat(">>> 生成 Table 3: GSE39582 临床特征汇总...\n")

n_normal_39582 <- sum(ds39582 == "Normal", na.rm = TRUE)
n_tumor_39582  <- sum(ds39582 == "Tumor",  na.rm = TRUE)
n_total_39582  <- length(ds39582)

# TNM 分期计数（用 sum 而非 table 索引，避免 NA 命名问题）
tnm_counts <- sapply(c("0", "1", "2", "3", "4"), function(lvl)
  sum(tnm39582 == lvl, na.rm = TRUE))
tnm_na_count <- sum(is.na(tnm39582))

# 分子亚型计数
ms_counts <- sapply(c("C1", "C2", "C3", "C4", "C5", "C6"), function(lvl)
  sum(ms39582 == lvl, na.rm = TRUE))
ms_na_count <- sum(is.na(ms39582))

# OS 事件计数
os_event_clean   <- os_event39582[!is.na(os_event39582)]
os_event1_count  <- sum(os_event_clean == "1")
os_event0_count  <- sum(os_event_clean == "0")
os_event_na_count <- sum(is.na(os_event39582))

pct <- function(n, d) sprintf("%.1f%%", n / d * 100)

clin_summary <- data.frame(
  Feature    = c(
    "Total Samples", "Normal Tissue", "Tumor Tissue", "",
    "TNM Stage 0", "TNM Stage 1", "TNM Stage 2", "TNM Stage 3", "TNM Stage 4", "TNM Stage NA", "",
    "CIT C1", "CIT C2", "CIT C3", "CIT C4", "CIT C5", "CIT C6", "CIT NA", "",
    "OS Event = 1 (deceased)", "OS Event = 0 (alive/censored)", "OS Event NA"
  ),
  N          = c(
    n_total_39582, n_normal_39582, n_tumor_39582, NA,
    tnm_counts, tnm_na_count, NA,
    ms_counts, ms_na_count, NA,
    os_event1_count, os_event0_count, os_event_na_count
  ),
  Percentage = c(
    "100.0%", pct(n_normal_39582, n_total_39582), pct(n_tumor_39582, n_total_39582), NA,
    sapply(tnm_counts, pct, d = n_total_39582), pct(tnm_na_count, n_total_39582), NA,
    sapply(ms_counts,  pct, d = n_total_39582), pct(ms_na_count,  n_total_39582), NA,
    pct(os_event1_count,   length(os_event_clean)),
    pct(os_event0_count,   length(os_event_clean)),
    pct(os_event_na_count, n_total_39582)
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write.csv(clin_summary, "results/tables/Table3_GSE39582_clinical_summary.csv",
          row.names = FALSE)
cat("  保存: Table3_GSE39582_clinical_summary.csv\n")

# ============================================================
# Section 9b: 泛素相关基因库构建（UbRG Gene Set）
#
# 目的: 生成 ubiquitin_genes.rds，供 Step 2 差异分析时筛选 UbRG
#
# 【来源】与 GeneSelection2.R 的原始设计完全一致:
#   主体: org.Hs.eg.db 本地 GO 注释库（6 个泛素核心 GO 条目）
#   补充: 核心手工整理基因（E1/E2/E3/DUB/蛋白酶体/SUMO/NEDD8 等）
#
# GO 条目:
#   GO:0016567  protein ubiquitination
#   GO:0006511  ubiquitin-dependent protein catabolic process
#   GO:0043161  proteasome-mediated ubiquitin-dependent protein catabolic
#   GO:0000209  protein polyubiquitination
#   GO:0016579  protein deubiquitination
#   GO:0031396  regulation of protein ubiquitination
#
# 【核查点】运行后打印基因数，预期约 1221 个（GO ~1202 + 手工补充 ~21）
# ============================================================
cat("\n>>> [Step1] 构建泛素相关基因库 (UbRG)...\n")

# ── 加载 GO 注释包 ─────────────────────────────────────────────
suppressMessages(library(org.Hs.eg.db))

# ── 主体：从 org.Hs.eg.db 本地 GO 注释库提取 ──────────────────
ub_go_terms <- c(
  "GO:0016567",  # protein ubiquitination
  "GO:0006511",  # ubiquitin-dependent protein catabolic process
  "GO:0043161",  # proteasome-mediated ubiquitin-dependent protein catabolic
  "GO:0000209",  # protein polyubiquitination
  "GO:0016579",  # protein deubiquitination
  "GO:0031396"   # regulation of protein ubiquitination
)

go2gene <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = ub_go_terms,
    columns = c("SYMBOL", "GOALL"),
    keytype = "GOALL"
  )
)
genes_go <- unique(as.character(na.omit(go2gene$SYMBOL)))
cat(sprintf("  [GO 注释库] %d 个泛素相关基因（6 个核心 GO 条目）\n",
            length(genes_go)))

# ── 补充：核心手工整理基因（GO 漏网的重要 UPS 成员）──────────────
genes_manual <- c(
  # E1 激活酶
  "UBA1", "UBA2", "UBA3", "UBA5", "UBA6", "UBA7",
  # E2 结合酶
  "UBE2A", "UBE2B", "UBE2C", "UBE2D1", "UBE2D2", "UBE2D3", "UBE2D4",
  "UBE2E1", "UBE2E2", "UBE2F", "UBE2G1", "UBE2G2", "UBE2H", "UBE2I",
  "UBE2J1", "UBE2J2", "UBE2K", "UBE2L3", "UBE2L6", "UBE2M", "UBE2N",
  "UBE2O", "UBE2Q1", "UBE2Q2", "UBE2R1", "UBE2R2", "UBE2S", "UBE2T",
  "UBE2U", "UBE2V1", "UBE2V2", "UBE2W", "UBE2Z",
  # E3 连接酶（RING/HECT/RBR）
  "MDM2", "BRCA1", "TRAF2", "TRAF6", "SIAH1", "SIAH2", "STUB1",
  "FBXW7", "FBXO4", "FBXO11", "FBXL3", "SKP2", "BTRC",
  "CUL1", "CUL2", "CUL3", "CUL4A", "CUL4B", "CUL5", "CUL7",
  "RBX1", "RBX2", "NEDD4", "NEDD4L", "ITCH",
  "HERC1", "HERC2", "HERC3", "HERC4", "HUWE1",
  "WWP1", "WWP2", "TRIP12", "HACE1", "SMURF1", "SMURF2",
  "RNF2", "RNF8", "RNF43", "ZNRF3", "RNF20", "RNF40",
  "MARCH1", "MARCH2", "MARCHF8", "MUL1", "MARCHF5",
  "VHL", "ELOC", "ELOB", "KEAP1", "SPOP",
  "PARK2", "PRKN", "ARIH1", "ARIH2",
  # 去泛素化酶（DUBs）
  "USP1", "USP2", "USP3", "USP4", "USP5", "USP6", "USP7", "USP8", "USP9X",
  "USP10", "USP11", "USP12", "USP13", "USP14", "USP15", "USP16", "USP17L2",
  "USP18", "USP19", "USP20", "USP21", "USP22", "USP24", "USP25", "USP26",
  "USP27X", "USP28", "USP29", "USP30", "USP31", "USP32", "USP33", "USP34",
  "USP36", "USP37", "USP38", "USP39", "USP40", "USP41", "USP42", "USP43",
  "USP44", "USP45", "USP46", "USP47", "USP48", "USP49", "USP50", "USP51",
  "USP52", "USP53", "USP54",
  "OTUB1", "OTUB2", "OTUD1", "OTUD3", "OTUD4", "OTUD5",
  "OTUD6A", "OTUD6B", "OTUD7A", "OTUD7B", "VCPIP1", "BAP1",
  "JOSD1", "JOSD2", "UCHL1", "UCHL3", "UCHL5", "CYLD", "TNFAIP3",
  # 蛋白酶体亚基
  "PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMA8",
  "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7", "PSMB8",
  "PSMB9", "PSMB10", "PSMB11",
  "PSMC1", "PSMC2", "PSMC3", "PSMC4", "PSMC5", "PSMC6",
  "PSMD1", "PSMD2", "PSMD3", "PSMD4", "PSMD5", "PSMD6", "PSMD7",
  "PSMD8", "PSMD9", "PSMD10", "PSMD11", "PSMD12", "PSMD13", "PSMD14",
  # 泛素本身
  "UBB", "UBC", "UBA52", "RPS27A",
  # SUMO/NEDD8 通路
  "SUMO1", "SUMO2", "SUMO3", "SUMO4",
  "SAE1", "UBA2", "UBE2I",
  "NEDD8", "NAE1", "UBA3", "UBE2M", "UBE2F",
  # 关键调控蛋白
  "VCP", "UFD1L", "NPLOC4", "UBQLN1", "UBQLN2", "UBQLN4",
  "RAD23A", "RAD23B", "SQSTM1", "NBR1", "TAX1BP1", "CALCOCO2"
)

# ── 合并去重 ───────────────────────────────────────────────────
ubiquitin_genes <- unique(c(genes_go, genes_manual))
n_supplement    <- length(setdiff(genes_manual, genes_go))
cat(sprintf("  [手工补充] %d 个（GO 未覆盖的核心 UPS 基因）\n", n_supplement))
cat(sprintf("  [合并总数] %d 个 UbRG（目标 ~1221）\n",
            length(ubiquitin_genes)))

# ------------------------------------------------------------
# 统计各数据集的覆盖情况（仅供参考，不过滤）
# Step2 会各自与本数据集表达矩阵取交集，无需在此预先过滤
# ------------------------------------------------------------
all_expr_genes <- unique(c(
  rownames(expr49355_norm),
  rownames(expr39582_norm),
  rownames(expr159216_norm),
  rownames(expr72718_norm)
))

n_total_ubrg    <- length(ubiquitin_genes)
n_in_any        <- length(intersect(ubiquitin_genes, all_expr_genes))
n_in_39582      <- length(intersect(ubiquitin_genes, rownames(expr39582_norm)))
n_in_49355      <- length(intersect(ubiquitin_genes, rownames(expr49355_norm)))

cat(sprintf("  UbRG 总数目:           %d 个\n", n_total_ubrg))
cat(sprintf("  在任一数据集中覆盖:    %d 个\n", n_in_any))
cat(sprintf("  GSE39582 中覆盖:       %d 个（差异分析用）\n", n_in_39582))
cat(sprintf("  GSE49355 中覆盖:       %d 个（差异验证用）\n", n_in_49355))
cat("  【注】保存完整列表，Step2 各自与表达矩阵取交集\n")

# 保存
saveRDS(ubiquitin_genes, "data/processed/ubiquitin_genes.rds")
cat(sprintf("\n  保存: data/processed/ubiquitin_genes.rds (%d 个 UbRG)\n",
            length(ubiquitin_genes)))

# ============================================================
# Section 10: 保存最终处理结果（供后续 Step 2-5 使用）
# ============================================================
cat("\n>>> [Step1] 保存最终结果...\n")

# 标准化表达矩阵（行: Gene Symbol, 列: Sample ID）
saveRDS(expr49355_norm,       "data/processed/expr49355_norm.rds")
saveRDS(expr159216_norm,      "data/processed/expr159216_norm.rds")
saveRDS(expr39582_norm,       "data/processed/expr39582_norm.rds")
saveRDS(expr72718_norm,       "data/processed/expr72718_norm.rds")

# 表型/临床数据（来自 pData）
saveRDS(pheno49355,  "data/processed/pheno49355.rds")
saveRDS(pheno39582,  "data/processed/pheno39582.rds")
saveRDS(pheno159216, "data/processed/pheno159216.rds")
saveRDS(pheno72718,  "data/processed/pheno72718.rds")

# 分组信息（命名字符向量）
saveRDS(group49355,  "data/processed/group49355.rds")
saveRDS(ds39582,     "data/processed/group_39582.rds")  

# GSE39582 完整临床注释（Step 3 预后模型需要）
# 【注意】os.delay (months):ch1 含括号，使用 [[]] 访问
# ── P0 : 数值列统一转为 numeric，清理 N/A → NA ──────────────────────
clean_numeric <- function(x) {
  # 输入 character 向量（如 "10", "N/A", NA），输出 numeric
  x[x == "N/A" | x == "" | x == "NA" | is.na(x)] <- NA_character_
  as.numeric(x)
}

os_event_val <- clean_numeric(os_event39582)
os_delay_val <- clean_numeric(pheno39582[["os.delay (months):ch1"]])
rfs_event_val <- clean_numeric(pheno39582[["rfs.event:ch1"]])
rfs_delay_val <- clean_numeric(pheno39582[["rfs.delay:ch1"]])

# ── P2: dataset 统一转为 factor ────────────────────────────────────
dataset_factor <- factor(ds39582, levels = c("Normal", "Tumor"))
age_val  <- clean_numeric(pheno39582[["age.at.diagnosis (year):ch1"]])
sex_val  <- pheno39582[["Sex:ch1"]]
sex_val  <- ifelse(sex_val %in% c("M", "F", "Male", "Female"), sex_val, NA_character_)

saveRDS(list(
  dataset           = dataset_factor,
  tnm_stage         = tnm39582,
  molecular_subtype = ms39582,
  age               = age_val,
  sex               = sex_val,
  os_event          = os_event_val,
  os_delay          = os_delay_val,
  rfs_event         = rfs_event_val,
  rfs_delay         = rfs_delay_val), "data/processed/clin_39582.rds")

cat("\n")
cat("================================================================\n")
cat("  Step 1 完成！\n")
cat("  标准化数据  → data/processed/\n")
cat("  QC 图表     → results/figures/   (QC1 ~ QC17)\n")
cat("  汇总表格    → results/tables/    (Table 1~3)\n")
cat("  下一步      → 运行 CRLM_UbRG_Step2_GeneSelection.R\n")
cat("================================================================\n")
