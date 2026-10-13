#===============================================================
## S12a_prepare_spatial.R —— 单细胞 / 空间数据解压与归位（一次性准备，幂等）
##
## 目的：把 GEO RAW.tar 的原始内容整理成 S12 可直读的标准布局。
##       本脚本**只做数据归位**，不做任何分析；重复运行即跳过已就位项。
##
## 产出目录布局：
##   01_data/GSE237230/<GSMid>_<sample>/barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
##   01_data/GSE247238/<GSMid>_<sample>/barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
##   01_data/GSE212388/<GSMid>_<sample>/barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
##       —— 均为 Seurat::Read10X 可直读的 10X 目录；宽表样本由本脚本转换为 10X 三件套
##   01_data/GSE277170/GSM85157xx_DSP-*.dcc.gz           —— GeoMx WTA dcc，原样保留
#===============================================================

## PART 0 常量 ------------------------------------------------------------------------------

dir.data <- "01_data"
n.expect <- c(GSE237230 = 8L, GSE247238 = 19L, GSE212388 = 8L, GSE277170 = 56L)

## PART 1 解压（目标目录已有内容则跳过）------------------------------------------------------

vec.tar <- c(
  GSE237230 = "GSE237230_RAW.tar", GSE247238 = "GSE247238_RAW.tar",
  GSE212388 = "GSE212388_RAW.tar", GSE277170 = "GSE277170_RAW.tar")
for (acc in names(vec.tar)) {
  dir.acc <- file.path(dir.data, acc)
  dir.create(dir.acc, showWarnings = FALSE, recursive = TRUE)
  n.files <- length(list.files(dir.acc))
  if (n.files >= n.expect[acc]) {
    cat(sprintf("[skip] %s 已有 %d 个文件，无需解压\n", acc, n.files))
    next
  }
  fp.tar <- file.path(dir.data, vec.tar[acc])
  if (!file.exists(fp.tar)) stop(paste0("S12a: 缺 ", vec.tar[acc]))
  ## RAW.tar 是未压缩 tar（首字节即文件名），禁用 -z；gzip 报 "not in gzip format" 是误判信号
  utils::untar(fp.tar, exdir = dir.acc, tar = "internal")
  cat(sprintf("[tar ] %s -> %d 个文件\n", acc,
      length(list.files(dir.acc))))
}

## PART 2 10X 三件套归位 --------------------------------------------------------------------
## Seurat::Read10X 要求目录内文件名严格为 barcodes.tsv.gz / features.tsv.gz /
## matrix.mtx.gz；GEO 文件名带 GSM 前缀，须重命名并移入每样本子目录。

PrepareTenX <- function(acc) {
  dir.acc <- file.path(dir.data, acc)
  vec.files <- list.files(dir.acc,
    pattern = "_barcodes\\.tsv\\.gz$|_features\\.tsv\\.gz$|_matrix\\.mtx\\.gz$")
  if (length(vec.files) == 0L) {
    cat(sprintf("[skip] %s 无待归位三件套（可能已归位）\n", acc))
    return(invisible(NULL))
  }
  vec.gsm <- unique(sub("^(GSM[0-9]+)_[^_]+_.*$", "\\1", vec.files))
  for (gsm in vec.gsm) {
    vec.f <- vec.files[startsWith(vec.files, gsm)]
    if (length(vec.f) == 0L) next
    tag <- sub(paste0("^", gsm, "_"), "", vec.f[1L])
    tag <- sub("_(barcodes|features|matrix)\\..*$", "", tag)
    dir.samp <- file.path(dir.acc, paste0(gsm, "_", tag))
    dir.create(dir.samp, showWarnings = FALSE, recursive = TRUE)
    for (f in vec.f) {
      target <- sub(paste0("^", gsm, "_", tag, "_"), "", f)
      file.rename(file.path(dir.acc, f), file.path(dir.samp, target))
    }
  }
  cat(sprintf("[10x ] %s -> %d 个样本子目录\n", acc,
      length(list.dirs(dir.acc, recursive = FALSE))))
}

PrepareTenX("GSE237230")
PrepareTenX("GSE247238")
PrepareTenX("GSE212388")

## PART 3 GSE212388 宽表 -> 10X 三件套 -------------------------------------------------------
## 4 个样本（GSM6529254/56/58/59）的 `.txt.gz` 是 gene × cell 宽表（首列基因，
## 首行 6000 个细胞 barcode）；另 4 个样本已由 PART 2 归位。

PrepareWideToTenX <- function(acc) {
  dir.acc <- file.path(dir.data, acc)
  ## GEO 文件名如 GSM6529254_101_PhysioFlow_Control_Sample_001.txt.gz（多下划线），
  ## 不能用 [^_]+ 只匹配一段
  vec.wide <- list.files(dir.acc, pattern = "^GSM[0-9]+_.*\\.txt\\.gz$")
  if (length(vec.wide) == 0L) {
    cat(sprintf("[skip] %s 无宽表待转换\n", acc))
    return(invisible(NULL))
  }
  for (f in vec.wide) {
    gsm <- sub("^([GSM0-9]+)_.*$", "\\1", f)
    tag <- sub(paste0("^", gsm, "_"), "", f)
    tag <- sub("\\.txt\\.gz$", "", tag)
    dir.samp <- file.path(dir.acc, paste0(gsm, "_", tag))
    if (file.exists(file.path(dir.samp, "matrix.mtx.gz"))) {
      cat(sprintf("[skip] %s 已转换\n", tag))
      next
    }
    dir.create(dir.samp, showWarnings = FALSE, recursive = TRUE)
    con <- gzfile(file.path(dir.acc, f), "rt")
    dat <- utils::read.delim(con, row.names = 1L, check.names = FALSE)
    close(con)
    cat(sprintf("[wide] %s : %d gene x %d cell\n", tag, nrow(dat), ncol(dat)))
    ## 写 10X：matrix.mtx（dgCMatrix）、features.tsv、barcodes.tsv；gzip 压缩
    mat.sparse <- Matrix::Matrix(as.matrix(dat), sparse = TRUE)
    Matrix::writeMM(mat.sparse, file.path(dir.samp, "matrix.mtx"))
    utils::write.table(data.frame(rownames(dat), rownames(dat),
        "Gene Expression"),
      file.path(dir.samp, "features.tsv"), sep = "\t",
      quote = FALSE, row.names = FALSE, col.names = FALSE)
    ## writeLines 属 base 包，非 utils —— 硬指名会报 "not an exported object"
    writeLines(colnames(dat), file.path(dir.samp, "barcodes.tsv"))
    R.utils::gzip(file.path(dir.samp, "matrix.mtx"), overwrite = TRUE)
    R.utils::gzip(file.path(dir.samp, "features.tsv"), overwrite = TRUE)
    R.utils::gzip(file.path(dir.samp, "barcodes.tsv"), overwrite = TRUE)
    ## 原始宽表已转换为 10X 三件套，删去副本避免重复
    unlink(file.path(dir.acc, f))
  }
}

PrepareWideToTenX("GSE212388")

## PART 4 样本数断言 ------------------------------------------------------------------------

for (acc in names(n.expect)) {
  dir.acc <- file.path(dir.data, acc)
  if (acc == "GSE277170") {
    n.obs <- length(list.files(dir.acc, pattern = "\\.dcc\\.gz$"))
  } else {
    n.obs <- length(list.dirs(dir.acc, recursive = FALSE))
  }
  cat(sprintf("[assert] %-11s 实测 %2d / 期望 %2d  ->  %s\n",
      acc, n.obs, n.expect[acc],
      if (n.obs == n.expect[acc]) "PASS" else "FAIL"))
  if (n.obs != n.expect[acc]) {
    stop(paste0("S12a: ", acc, " 样本数不符（", n.obs, " != ", n.expect[acc],
        "）"))
  }
}

message("S12a OK.")
