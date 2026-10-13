###############################################################################
## S09a_build_reference.R —— 构建 BayesPrism 动脉参考 GEP（S09 前置）
##
## PART 目的：把动脉粥样硬化斑块单细胞 atlas（Traeuble et al. Nat Commun 2025,
##       16:8255, DOI 10.1038/s41467-025-63202-x；同疾病直接先例——该文自身
##       即用 BayesPrism 对 bulk RNA-seq 去卷积）按细胞类型聚合成 type-by-gene
##       GEP，供 S09 的 BayesPrism 通道以 input.type="GEP" 作先验参考。
##
## PART 方法：CELLxGENE Plaque Atlas h5ad 的 /raw/X（CSR 稀疏 counts）经
##       rhdf5 直读；按细胞类型分块累加求 gene-wise 和；再按 biotype 过滤、
##       同名 symbol 合并。全程不把 259k x 25k 矩阵稠密化。
##
## PART 数据：01_data/ad9f52c0-054b-4dc5-944c-46c376a544cf.h5ad
##       （259,116 细胞 x 25,379 基因，CELLxGENE schema 7.1.0，
##         collection db70986c-7d91-49fe-a399-a4730be394ac）
##
## 产出：01_data/references/BayesPrism_reference_atherosclerosis.rds
##         = list(gep               type-by-gene 求和 counts 矩阵
##                cell.type.labels  细胞类型名（= rownames(gep)）
##                cell.state.labels 状态名（state == type 粒度）
##                provenance        实测数据来源、聚合统计与 QA 断言结果)
##       03_results/S09a/table/reference_summary.csv（长表）
##       03_results/S09a/table/reference_celltypes.csv（逐类型统计）
##       03_results/S09a/checkpoints/（PART 级断点）
##
## 方法依据（全部 2026-09-28 实测于本机 R 4.6.1 + rhdf5 2.56.1）：
##   [1] BayesPrism new.prism.Rd 原文：GEP 输入用 "summed raw count for each
##       cell state is preferred"，并明写 "log transformation should be
##       avoided" —— 故本脚本必须用 /raw/X 的 counts，绝不可用 /X。
##   [2] zellkonverter 1.22.0 的 readH5AD(reader="R") 内部 .native_reader
##       只读 X 与 layers，**忽略 anndata 的 /raw 组**。实测 /X 首值
##       1.899254 / 1.345590 / 2.514612（已变换），/raw/X 首值 2/1/2/1…
##       （整数 counts，首行 802 个非零，范围 1~70）—— 二者不可互换。
##   [3] 实测 rhdf5 2.56.1 无 h5open（仅 H5Fopen/H5Fclose），H5Fopen 不接受
##       "r" flags 字符串；h5ls() 无 name 形参；h5dim() 未导出。
##
###############################################################################

## PART 0 常量与目录 ----------------------------------------------------------

setwd("F:/RAPA")

H5AD_PATH <- "01_data/ad9f52c0-054b-4dc5-944c-46c376a544cf.h5ad"
RDS_PATH <- "01_data/references/BayesPrism_reference_atherosclerosis.rds"
TAB_DIR <- "03_results/S09a/table"
CKPT_DIR <- "03_results/S09a/checkpoints"
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(RDS_PATH), showWarnings = FALSE, recursive = TRUE)

## 注释列名（CELLxGENE curated schema 标准字段）
CELLTYPE_COL <- "cell_type"     ## /obs/cell_type      —— 实测 12 类
GENE_COL <- "feature_name"  ## /var/feature_name   —— 实测 25,377 唯一 symbol
BIOTYPE_COL <- "feature_type"   ## /var/feature_type   —— 实测 17 类

KEEP_BIOTYPE <- "protein_coding"
MIN_CELLS_PER_TYPE <- 20L       ## 少于该细胞数的类型 GEP 噪声过大，剔除并登记
CHUNK_CELLS <- 10000L           ## 每块细胞数：控制单次 rowsum 的内存峰值

LOG_PATH <- file.path(CKPT_DIR, "S09a_build.log")
f.log <- function(...) cat(..., "\n", sep = "", file = LOG_PATH, append = TRUE)
f.log("=== S09a_build_reference.R started ===")
f.log(paste("R", getRversion(), "| rhdf5", as.character(packageVersion(
        "rhdf5")),
    "| Matrix", as.character(packageVersion("Matrix"))))

## PART 1 打开 h5 并读元信息 ---------------------------------------------------

if (!requireNamespace("rhdf5", quietly = TRUE)) {
  stop("S09a STOP: rhdf5 not installed. Install from Bioconductor first.")
}
if (!file.exists(H5AD_PATH)) {
  stop(paste0("S09a STOP: ", H5AD_PATH, " not found. Get the CELLxGENE Plaque ",
      "Atlas h5ad (collection db70986c-7d91-49fe-a399-a4730be394ac) ",
      "and place it at that path."))
}

h5 <- rhdf5::H5Fopen(H5AD_PATH)
h5.tree <- rhdf5::h5ls(H5AD_PATH, all = TRUE, recursive = TRUE)

## 通用列读取：CELLxGENE 的注释列既可能是纯 dataset，也可能是 categorical
## （组下挂 categories + codes 两个 dataset），两种形态按实测都支持。
ReadH5Column <- function(path) {
  vec.kids <- h5.tree$name[h5.tree$group == path]
  if ("codes" %in% vec.kids) {
    vec.cats <- as.character(rhdf5::h5read(h5, paste0(path, "/categories")))
    vec.codes <- as.integer(rhdf5::h5read(h5, paste0(path, "/codes")))
    vec.out <- vec.cats[vec.codes + 1L]
    vec.out[is.na(vec.codes)] <- NA_character_
    vec.out
  } else {
    as.character(rhdf5::h5read(h5, path))
  }
}

vec.title <- as.character(rhdf5::h5read(h5, "/uns/title"))
vec.citation <- as.character(rhdf5::h5read(h5, "/uns/citation"))
vec.schema <- as.character(rhdf5::h5read(h5, "/uns/schema_version"))
vec.organism <- as.character(rhdf5::h5read(h5, "/uns/organism"))

f.log(paste("file:", H5AD_PATH, "(", format(file.size(H5AD_PATH),
      big.mark = ","), "bytes )"))
f.log(paste("title:", vec.title, "| schema:", vec.schema,
    "| organism:", vec.organism))

## PART 2 读 obs / var 注释，做列名硬断言 --------------------------------------

vec.celltype <- ReadH5Column(paste0("/obs/", CELLTYPE_COL))
vec.gene <- ReadH5Column(paste0("/var/", GENE_COL))
vec.biotype <- ReadH5Column(paste0("/var/", BIOTYPE_COL))

f.log(paste("cell types:", length(unique(vec.celltype)),
    "| features:", length(vec.gene),
    "| biotypes:", length(unique(vec.biotype))))

idx.keep.biotype <- which(vec.biotype %in% KEEP_BIOTYPE)
if (length(idx.keep.biotype) == 0L) {
  stop(sprintf(paste0(
        "S09a STOP: no feature with biotype '%s'. Actual biotypes: %s"),
      KEEP_BIOTYPE, paste(unique(vec.biotype), collapse = ", ")))
}
f.log(paste("biotype filter", KEEP_BIOTYPE, ":",
    length(idx.keep.biotype), "of", length(vec.gene), "features kept"))

## PART 3 读 raw counts CSR 并按细胞类型累加 ------------------------------------

vec.indptr <- as.integer(rhdf5::h5read(h5, "/raw/X/indptr"))
vec.indices <- as.integer(rhdf5::h5read(h5, "/raw/X/indices"))
vec.data <- rhdf5::h5read(h5, "/raw/X/data")
rhdf5::H5Fclose(h5)

n.cell <- length(vec.celltype)
n.gene <- length(vec.gene)
if (length(vec.indptr) != n.cell + 1L) {
  stop(sprintf("S09a STOP: CSR indptr length %d != n.cell + 1 (%d).",
      length(vec.indptr), n.cell + 1L))
}
if (length(vec.data) != length(vec.indices)) {
  stop(sprintf("S09a STOP: CSR data (%d) / indices (%d) length mismatch.",
      length(vec.data), length(vec.indices)))
}
dbl.total.counts <- sum(vec.data)
f.log(paste("raw counts CSR:", n.cell, "cells x", n.gene, "genes,",
    format(length(vec.data), big.mark = ","), "non-zeros"))

vec.celltype.uniq <- sort(unique(vec.celltype))
mat.gep <- matrix(0, nrow = length(vec.celltype.uniq), ncol = n.gene,
  dimnames = list(vec.celltype.uniq, vec.gene))
vec.cells.used <- setNames(
  rep(0L, length(vec.celltype.uniq)), vec.celltype.uniq)

vec.chunk <- split(seq_len(n.cell), ceiling(seq_len(n.cell) / CHUNK_CELLS))
f.log(paste("aggregating in", length(vec.chunk), "cell chunks ..."))

for (.k in seq_along(vec.chunk)) {
  .idx <- vec.chunk[[.k]]
  .rng <- c(vec.indptr[.idx[1L]], vec.indptr[.idx[length(.idx)] + 1L])
  ## CSR indices 为 0 基（实测 /raw/X/indices 范围 0~25378），
  ## 转为 R 的 1 基列下标，否则含 0 的下标会使替换长度错配。
  .gi <- vec.indices[.rng[1L]:.rng[2L]] + 1L
  .vv <- vec.data[.rng[1L]:.rng[2L]]
  ## .vv 是本块所有细胞的非零值拼成的扁平向量，因此每个非零元素还要
  ## 回溯它属于哪个细胞（进而属于哪个类型）；直接用细胞序号去索引 .vv
  ## 会取错元素（本脚本曾因此丢掉 99.9% 的 counts，）。
  ## rep(向量 times) 才按元素展开；rep.int(x, 标量 times) 是"每个元素重复
  ## times 次"，会把本块放大 length(.idx) 倍（曾尝试申请 764 Gb 而失败）。
  ## 逐细胞非零数：长度必须等于 length(.idx)（diff() 会少一个元素，
  ## rep() 要求 times 与 x 等长，故不用 diff）。
  .nnz.per.cell <- vec.indptr[.idx + 1L] - vec.indptr[.idx]
  .el.cell <- rep(seq_along(.idx), times = as.double(.nnz.per.cell))
  .el.type <- match(vec.celltype[.idx], vec.celltype.uniq)[.el.cell]
  for (.u in unique(.el.type)) {
    .sel <- which(.el.type == .u)
    .s <- rowsum(.vv[.sel], group = .gi[.sel], reorder = TRUE)
    .i <- as.integer(rownames(.s))
    .ct <- vec.celltype.uniq[.u]
    mat.gep[.ct, .i] <- mat.gep[.ct, .i] + .s[, 1L]
    vec.cells.used[.ct] <- vec.cells.used[.ct] + length(unique(.el.cell[.sel]))
  }
  if (.k %% 5L == 0L || .k == length(vec.chunk)) {
    f.log(paste("  chunk", .k, "of", length(vec.chunk), "done"))
  }
}

## PART 级断点：此后任何失败都可不经重算聚合而续跑
saveRDS(list(gep = mat.gep, cells.used = vec.cells.used),
  file.path(CKPT_DIR, "checkpoint_gep.rds"))

## QA1 counts 总量守恒（独立由 CSR 重新求和）
dbl.gep.sum <- sum(mat.gep)
if (!isTRUE(all.equal(dbl.gep.sum, dbl.total.counts))) {
  stop(sprintf(paste0("S09a STOP: QA1 total-counts conservation FAILED. ",
        "raw CSR = %.0f, aggregated GEP = %.0f"),
      dbl.total.counts, dbl.gep.sum))
}

## QA2 每类型实际参与累加的细胞数 == obs 中该类型的细胞数
vec.cells.obs <- table(vec.celltype)[vec.celltype.uniq]
vec.cells.bad <- vec.celltype.uniq[vec.cells.used != vec.cells.obs]
if (length(vec.cells.bad) > 0L) {
  stop(sprintf(paste0("S09a STOP: QA2 cell-count consistency FAILED for: %s ",
        "(obs vs. aggregated)"),
      paste(vec.cells.bad, collapse = ", ")))
}
f.log(paste("QA1 total-counts conserved:", format(dbl.total.counts,
      big.mark = ",")))
f.log(paste("QA2 cell counts consistent for all", length(vec.celltype.uniq),
    "types"))

## PART 4 biotype 过滤、同名合并、装配参考对象 ---------------------------------

## mat.gep 是 type-by-gene，biotype 过滤作用于「基因列」而非行
mat.gep <- mat.gep[, idx.keep.biotype, drop = FALSE]
vec.cells.obs <- vec.cells.obs[rownames(mat.gep)]
vec.cells.used <- vec.cells.used[rownames(mat.gep)]

## 同名 symbol（本文件 25,379 特征 -> 25,377 唯一 symbol）合并求和
if (anyDuplicated(colnames(mat.gep)) > 0L) {
  mat.gep <- t(rowsum(t(mat.gep), group = colnames(mat.gep), reorder = TRUE))
}
f.log(paste("after biotype filter + symbol merge:", nrow(mat.gep), "types x",
    ncol(mat.gep), "genes"))

vec.type.low <- rownames(mat.gep)[vec.cells.obs < MIN_CELLS_PER_TYPE]
if (length(vec.type.low) > 0L) {
  stop(sprintf(paste0(
        "S09a STOP: cell type(s) below MIN_CELLS_PER_TYPE=%d: %s ",
        ". Raise biotype coverage or lower the threshold."),
      MIN_CELLS_PER_TYPE, paste(vec.type.low, collapse = ", ")))
}

dat.prov <- data.frame(
  cell.type = rownames(mat.gep),
  n.cells = as.integer(vec.cells.obs),
  total.counts = apply(mat.gep, 1L, sum),
  detected.features = apply(mat.gep, 1L, function(.r) sum(.r > 0)),
  stringsAsFactors = FALSE)

out <- list(
  gep = mat.gep,
  cell.type.labels = rownames(mat.gep),
  cell.state.labels = rownames(mat.gep),
  provenance = list(
    source.file = H5AD_PATH,
    source.bytes = file.size(H5AD_PATH),
    dataset.title = vec.title,
    dataset.citation = vec.citation,
    cellxgene.collection = "db70986c-7d91-49fe-a399-a4730be394ac",
    schema.version = vec.schema,
    organism = vec.organism,
    cells = n.cell,
    features = n.gene,
    nonzeros = length(vec.data),
    biotype.filter = KEEP_BIOTYPE,
    types.kept = nrow(dat.prov),
    symbols.after.merge = ncol(mat.gep),
    r.version = as.character(getRversion()),
    rhdf5.version = as.character(packageVersion("rhdf5")),
    built.at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    cell.types = dat.prov,
    qa.total.counts.raw = dbl.total.counts,
    qa.total.counts.gep = dbl.gep.sum))

saveRDS(out, RDS_PATH)
f.log(paste("reference written ->", RDS_PATH))

dat.summary <- rbind(
  data.frame(section = "source", key = "source.file",
    value = H5AD_PATH, stringsAsFactors = FALSE),
  data.frame(section = "source", key = "dataset.title",
    value = vec.title, stringsAsFactors = FALSE),
  data.frame(section = "source", key = "dataset.citation",
    value = vec.citation, stringsAsFactors = FALSE),
  data.frame(section = "source", key = "cellxgene.collection",
    value = "db70986c-7d91-49fe-a399-a4730be394ac",
    stringsAsFactors = FALSE),
  data.frame(section = "source", key = "schema.version",
    value = vec.schema, stringsAsFactors = FALSE),
  data.frame(section = "source", key = "organism",
    value = vec.organism, stringsAsFactors = FALSE),
  data.frame(section = "dimension", key = "cells",
    value = as.character(n.cell), stringsAsFactors = FALSE),
  data.frame(section = "dimension", key = "features",
    value = as.character(n.gene), stringsAsFactors = FALSE),
  data.frame(section = "dimension", key = "nonzeros",
    value = format(length(vec.data), big.mark = ","),
    stringsAsFactors = FALSE),
  data.frame(section = "dimension", key = "types.kept",
    value = as.character(nrow(dat.prov)), stringsAsFactors = FALSE),
  data.frame(section = "dimension", key = "symbols.after.merge",
    value = as.character(ncol(mat.gep)), stringsAsFactors = FALSE),
  data.frame(section = "filter", key = "biotype.filter",
    value = KEEP_BIOTYPE, stringsAsFactors = FALSE),
  data.frame(section = "filter", key = "min.cells.per.type",
    value = as.character(MIN_CELLS_PER_TYPE),
    stringsAsFactors = FALSE),
  data.frame(section = "qa", key = "qa.total.counts.raw",
    value = format(dbl.total.counts, big.mark = ","),
    stringsAsFactors = FALSE),
  data.frame(section = "qa", key = "qa.total.counts.gep",
    value = format(dbl.gep.sum, big.mark = ","),
    stringsAsFactors = FALSE),
  data.frame(section = "qa", key = "qa.result", value = "PASS",
    stringsAsFactors = FALSE))

write.csv(dat.summary, file.path(TAB_DIR, "reference_summary.csv"),
  row.names = FALSE)
write.csv(dat.prov, file.path(TAB_DIR, "reference_celltypes.csv"),
  row.names = FALSE)
f.log(paste("DONE:", nrow(mat.gep), "types x", ncol(mat.gep),
    "symbols ->", RDS_PATH))
