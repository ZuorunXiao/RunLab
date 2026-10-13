#===============================================================
## S01_ingest.R —— 数据全量落盘与设计断言
## 目的：把 01_data/ 下的全部在途数据实测落盘登记（体积 / MD5 / 角色 / 状态），
##       并对 GEO family SOFT 做样本级解析，逐条核对 §8 中的 QC 断言
##       Q1（样本数 + 平台）、Q7（双色流程边界 + GSE7084 负断言）、Q14（GSE232911）。
## 定位：分析脚本（规则 34，R 优先）。只做落盘核对与断言，不做归一化、不做统计。
## 产出：03_results/S01/table/MANIFEST.tsv            落盘清单（体积 / MD5 / 状态）
##       03_results/S01/table/design_check.csv        样本数 / 平台断言
##       03_results/S01/table/design_keys.tsv         各系列**样本级** characteristics 长表
##                                                     （acc / sample_id / key / value）
##       03_results/S01/table/gse232911_grouping.tsv  R 侧重算的分组表（canonical）
##       03_results/S01/table/gse232911_qc.csv        S01 全部断言清单（含 §3.2.1 专项 + Q7）
##       03_results/S01/table/sample_axis.tsv         来源列名 -> GSM 的重映射表
##       03_results/S01/table/arm_sizes.csv           各臂实测样本数（供 Q12）
##       03_results/S01/Figure/Sample_Composition_By_Layer.pdf   单面板柱：系列样本数（按层着色）
# 03_results/S01/Figure/Sample_Count_vs_Expected.pdf     单面板散点：预期 vs 实测样本数（y=x）
## 图件清单（规则 35 可视化强制 · SCI 发表导向；动笔前已调研已发表文献确定风格/形式/数量）：
##   可视化: S01 共 2 幅
##     1. Sample_Composition_By_Layer.pdf
##        —— 形式：单面板柱状图（x = 系列，y = 样本数，fill = 层），x 轴标签 45° 旋转
##        —— 回答的科学问题：每个在途数据集实解析到的样本量、分属哪一层，是否满足下游最小样本量。
##        —— 参照文献：数据集样本量构成图为转录组整合研究标准首图（如 AAA 组织 meta
##           PMID 42653260 的队列构成；GEO 系列 QC 惯例）。
##     2. Sample_Count_vs_Expected.pdf
##        —— 形式：单面板散点（x = 预期样本数，y = 实测样本数，y=x 参考线，按层着色）
##        —— 回答的科学问题：实测样本数是否与方案 §3 预注册登记值一致（S01 设计断言 Q1 的可视化）。
##        —— 参照文献：expected-versus-observed 一致性散点为数据核验标准画法（y=x 诊断线）。
#===============================================================

## PART 0 常量与注册表 -----------------------------------------------------------

dir.data <- "01_data"
dir.tab <- "03_results/S01/table"
FIG_DIR <- "03_results/S01/Figure"
dir.create(dir.tab, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

library(ggplot2)  ## QC 图；主题与出口严格用规范 §sci_theme / §sci_save

## 表格写出：禁 write.delim（该函数不存在于 base/utils，勿凭记忆调用）
WriteTsv <- function(dat, path) {
  write.table(dat, path, sep = "\t", quote = FALSE, row.names = FALSE,
    na = "NA")
}

## 绘图主题：规范 §285-300 官方定义（theme_bw + 关网格），全脚本统一口径，
## 不另起 theme_minimal —— 各脚本主题不一致会被审稿人挑。
sci_theme <- ggplot2::theme_bw(base_size = 12) +
ggplot2::theme(
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank(),
  ##   在只有 1–2 级分面时视觉权重过高，且与本主题「白底无网格」的 SCI 口径不一致。
  ##   全项目 13 个脚本的 sci_theme 同步修改，避免同项目风格漂移（规范 35）。
  strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
  legend.position  = "right",
  legend.text      = ggplot2::element_text(size = 10, face = "plain"),
  axis.title       = ggplot2::element_text(size = 12, face = "plain"),
  axis.text        = ggplot2::element_text(size = 10, face = "plain"),
  plot.title       = ggplot2::element_text(size = 14, hjust = 0.5,
    face = "plain"),
  legend.title     = ggplot2::element_text(size = 11, face = "plain")
)

## 图件出口：仅 PDF，写本步 FIG_DIR；规范 §sci_save 原文，print + dev.off 缺一不可，
## 否则产出 0 字节 PDF。
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

## 12 基因固定面板（§10 预注册锁定；禁止同义替换，Q4）
vec.panel12 <- c("IL6", "IL1B", "TNF", "CCL2", "CXCL8", "MMP9",
  "SPP1", "NLRP3", "IL10", "TGFB1", "MERTK", "GAS6")

## 在途数据集注册表：acc / 层 / 平台 / 期望样本数
## 期望值取自方案 §3.1（锚 / 扰动）与 §3.2（疾病 bulk）的登记值，逐条实测复核。
dat.reg <- rbind(
  data.frame(acc = "GSE124020", layer = "L-P", platform = "GPL17586",
    n_expect = 12L),
  data.frame(acc = "GSE124106", layer = "L-P", platform = "GPL17586",
    n_expect = 24L),
  data.frame(acc = "GSE323983", layer = "L-P", platform = "GPL28975",
    n_expect = 23L),
  data.frame(acc = "GSE65145",  layer = "L-P", platform = "GPL10999",
    n_expect = 6L),
  data.frame(acc = "GSE6002",   layer = "L-P", platform = "GPL3784",
    n_expect = 4L),
  data.frame(acc = "GSE212388", layer = "L-P", platform = "GPL18573",
    n_expect = 8L),
  data.frame(acc = "GSE232911", layer = "L-D", platform = "GPL17586",
    n_expect = 246L),
  data.frame(acc = "GSE57691",  layer = "L-D", platform = "GPL10558",
    n_expect = 68L),
  data.frame(acc = "GSE43292",  layer = "L-D", platform = "GPL6244",
    n_expect = 64L),
  data.frame(acc = "GSE21545",  layer = "L-D", platform = "GPL570",
    n_expect = 223L),
  data.frame(acc = "GSE100927", layer = "L-D", platform = "GPL17077",
    n_expect = 104L),
  data.frame(acc = "GSE28829",  layer = "L-D", platform = "GPL570",
    n_expect = 29L),
  data.frame(acc = "GSE47472",  layer = "L-D", platform = "GPL10558",
    n_expect = 22L),
  data.frame(acc = "GSE41571",  layer = "L-D", platform = "GPL570",
    n_expect = 11L),
  data.frame(acc = "GSE237230", layer = "L-C", platform = "GPL24676",
    n_expect = 8L),
  data.frame(acc = "GSE247238", layer = "L-C", platform = "GPL24676",
    n_expect = 19L),
  data.frame(acc = "GSE277170", layer = "L-C", platform = "GPL18573",
    n_expect = 56L)
)
dat.reg$layer <- as.character(dat.reg$layer)

## 排除集（§3.4 + 2026-09-26 用户裁定）：磁盘上出现即标 EXCLUDED，禁止进入取数清单
vec.excluded <- c("GSE7084", "GSE294466", "GSE315547", "GSE287081", "GSE25520",
  "GSE246235")

## 系列级角色文件（相对 01_data/）；其余按通用命名规则自动补齐
dat.files <- data.frame(
  acc  = c("GSE232911", "GSE323983", "GSE323983", "GSE41571", "GSE47472",
    "GSE57691"),
  role = c("matrix", "counts_block1", "counts_block2", "series_matrix",
    "series_matrix", "counts_raw"),
  path = c("GSE232911_Matrix.file.txt.gz",
    "GSE323983_CvIR_all_gene_counts.tsv.gz",
    "GSE323983_RapavIRRapa_all_gene_counts.tsv.gz",
    "GSE41571_series_matrix.txt.gz",
    "GSE47472_series_matrix.txt.gz",
    "GSE57691_non-normalized_data.txt.gz"),
  layer = c("L-D", "L-P", "L-P", "L-D", "L-D", "L-D"),
  note = c("", "", "", "", "", ""),
  stringsAsFactors = FALSE)

## 系列级 family SOFT 落点（verify/ 为历史落点，01_data/ 为可能的补落点）
dat.soft <- data.frame(
  acc = dat.reg$acc,
  path = file.path("verify", paste0(dat.reg$acc, "_family.soft.gz")),
  stringsAsFactors = FALSE)

## 单行例外：GSE232911 的 family SOFT 下载在 01_data/ 根，不在 verify/（2026-09-28 实测核对）
dat.soft$path[dat.soft$acc == "GSE232911"] <- "GSE232911_family.soft.gz"

## 平台注释文件（S03 归一化需要；不是表达表，单独登记）
dat.annot <- data.frame(
  acc = c("GSE100927", "GSE21545", "GSE28829", "GSE41571", "GSE47472",
    "GSE237230", "GSE247238", "GSE277170", "GSE323983", "GSE6002",
    "GSE57691", "GSE47472", "GSE124020"),
  role = c("annotation", "annotation", "annotation", "series_matrix",
    "series_matrix", "annotation", "annotation", "annotation",
    "annotation", "annotation", "annotation", "series_matrix",
    "annotation"),
  path = c("GPL17077.annot.csv.gz", "GPL570.annot.gz", "GPL6244.annot.gz",
    "GSE41571_series_matrix.txt.gz", "GSE47472_series_matrix.txt.gz",
    "GPL24676_family.soft.gz", "GPL24676_family.soft.gz",
    "GPL18573_family.soft.gz", "GPL28975_family.soft.gz",
    "GPL3784_family.soft.gz", "GPL10558.annot.gz",
    "GSE47472_series_matrix.txt.gz", "GPL17586-45144.txt"),
  note = c("", "", "", "表达表兼注释", "表达表兼注释", "", "", "", "", "", "",
    "表达表兼注释", "HTA-2.0 表达索引，非 GPL 注释"),
  stringsAsFactors = FALSE)

file.manifest <- file.path(dir.tab, "MANIFEST.tsv")
vec.cache.path <- character(0)
vec.cache.key <- character(0)
vec.cache.md5 <- character(0)
if (file.exists(file.manifest)) {
  dat.cache <- read.delim(file.manifest, header = TRUE,
    stringsAsFactors = FALSE)
  keep <- !is.na(dat.cache$md5)
  vec.cache.path <- as.character(dat.cache$path[keep])
  vec.cache.key <- paste(dat.cache$path[keep], dat.cache$size_bytes,
    dat.cache$mtime_utc, sep = "|")
  vec.cache.md5 <- as.character(dat.cache$md5[keep])
}

## PART 1 落盘清单（体积 / MD5 / 状态） ---------------------------------------------

RegisterFile <- function(rel.path, acc, role, layer = "", note = "") {
  full <- file.path(dir.data, rel.path)
  if (!file.exists(full)) {
    return(data.frame(acc = acc, layer = layer, role = role, path = rel.path,
        size_bytes = NA_real_, size_mb = NA_real_,
        mtime_utc = NA_character_, md5 = NA_character_,
        status = "MISSING", note = note, stringsAsFactors = FALSE))
  }
  size.bytes <- file.size(full)
  key.now <- paste(rel.path, size.bytes,
    format(file.mtime(full), tz = "UTC",
      usetz = TRUE), sep = "|")
  i.cache <- match(key.now, vec.cache.key)
  data.frame(acc = acc, layer = layer, role = role, path = rel.path,
    size_bytes = size.bytes, size_mb = round(size.bytes / 1024^2, 2),
    mtime_utc = key.now,
    md5 = if (is.na(i.cache)) tools::md5sum(full) else vec.cache.md5[i.cache],
    status = "OK", note = note, stringsAsFactors = FALSE)
}

## 展开待查文件清单 = 显式指定 + 通用命名自动补齐
dat.job <- rbind(
  data.frame(acc = dat.files$acc, role = dat.files$role,
    layer = dat.files$layer,
    path = dat.files$path, note = dat.files$note, stringsAsFactors = FALSE),
  data.frame(acc = dat.soft$acc, role = "family_soft", layer = dat.reg$layer,
    path = dat.soft$path, note = "", stringsAsFactors = FALSE))
for (acc in dat.reg$acc) {
  for (sfx in c("_RAW.tar", "_series_matrix.txt.gz", "_Matrix.file.txt.gz",
      "_non-normalized_data.txt.gz")) {
    cand <- paste0(acc, sfx)
    if (!any(dat.job$path == cand) && !any(dat.annot$path == cand)) {
      dat.job <- rbind(dat.job, data.frame(
          acc = acc, role = paste0("auto", sfx),
        layer = dat.reg$layer[match(acc, dat.reg$acc)], path = cand, note = "",
          stringsAsFactors = FALSE))
    }
  }
}

list.manifest <- list()
for (k in seq_len(nrow(dat.annot))) {
  a <- dat.annot[k, ]
  if (basename(a$path) %in% basename(dat.job$path)) next
  list.manifest[[length(list.manifest) + 1L]] <-
  RegisterFile(a$path, a$acc, a$role, "", a$note)
}
for (k in seq_len(nrow(dat.job))) {
  j <- dat.job[k, ]
  if (j$acc %in% vec.excluded) {
    full <- file.path(dir.data, j$path)
    if (file.exists(full)) {
      list.manifest[[length(list.manifest) + 1L]] <- data.frame(
        acc = j$acc, layer = "EXCLUDED", role = j$role, path = j$path,
        size_bytes = file.size(full), size_mb = round(file.size(full) / 1024^2,
          2),
        mtime_utc = format(file.mtime(full), tz = "UTC", usetz = TRUE),
        md5 = tools::md5sum(full), status = "EXCLUDED",
        note = "§3.4 / 用户裁定剔除，禁止进入取数清单",
        stringsAsFactors = FALSE)
    }
    next
  }
  list.manifest[[length(list.manifest) + 1L]] <- RegisterFile(j$path, j$acc,
    j$role, j$layer,
    "")
}
dat.manifest <- do.call(rbind, list.manifest)
WriteTsv(dat.manifest, file.manifest)
cat(sprintf("[manifest] rows=%d  missing=%d  excluded=%d\n", nrow(dat.manifest),
    sum(dat.manifest$status == "MISSING"), sum(
      dat.manifest$status == "EXCLUDED")))

## PART 2 family SOFT 解析 ----------------------------------------------------------

## 解析 family SOFT 的样本块。
## 坑：字段前缀大小写不敏感——GEO 规范写作 `!Sample_`，部分导出为 `!sample_`；
##     按字面匹配会全部漏配**且不报错**（本项目 GSE7084 实测踩过）。
ParseFamilySoft <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  vec.line <- readLines(con, warn = FALSE)
  close(con)
  i.samp <- which(startsWith(vec.line, "^SAMPLE"))
  end.cum <- c(i.samp[-1L] - 1L, length(vec.line))
  out <- list(sample_id = character(length(i.samp)), platform = character(
      length(i.samp)),
    title = character(length(i.samp)), chars = vector("list",
      length(i.samp)))
  for (k in seq_along(i.samp)) {
    blk <- vec.line[i.samp[k]:end.cum[k]]
    out$sample_id[k] <- sub("^\\^SAMPLE\\s*=\\s*", "", blk[1L])
    i.plat <- grep("^!Sample_platform_id", blk, ignore.case = TRUE)
    if (length(i.plat) > 0L) out$platform[k] <- sub("^.*=\\s*", "",
      blk[i.plat[1L]])
    i.ttl <- grep("^!Sample_title", blk, ignore.case = TRUE)
    if (length(i.ttl) > 0L) out$title[k] <- sub("^.*=\\s*", "", blk[i.ttl[1L]])
    i.chr <- grep("^!Sample_characteristics_ch1", blk, ignore.case = TRUE)
    out$chars[[k]] <- sub("^!Sample_characteristics_ch1(\\.[0-9]+)?\\s*=\\s*",
      "",
      blk[i.chr], ignore.case = TRUE)
  }
  out
}

## 把 characteristics 字符串向量解析为 (key, value)；GEO 格式 "key: value"
ParseChars <- function(vec.chr) {
  if (length(vec.chr) == 0L) return(data.frame(key = character(0),
      value = character(0)))
  list.split <- strsplit(vec.chr, ":")
  key <- vapply(list.split, function(x) if (length(x) > 0L) trimws(x[
        1L]) else NA_character_,
    "")
  val <- vapply(list.split, function(x) trimws(paste(x[-1L], collapse = ":")),
    "")
  keep <- !is.na(key)
  data.frame(key = key[keep], value = val[keep], stringsAsFactors = FALSE)
}

## PART 3 逐系列样本数 / 平台断言（Q1） ----------------------------------------------

dat.design <- data.frame(acc = character(0), layer = character(0),
  n_samples = integer(0),
  n_expect = integer(0), platform = character(0),
  platform_expect = character(0), platform_host = character(0),
  result = character(0), stringsAsFactors = FALSE)
dat.keys <- data.frame(acc = character(0), sample_id = character(0),
  key = character(0),
  value = character(0), stringsAsFactors = FALSE)

## QC 累积容器：定义在 PART 3 之前，供平台附注与 PART 4 的 GSE232911 断言共用
dat.qc <- data.frame(q = character(0), result = character(0),
  observed = character(0),
  expected = character(0), stringsAsFactors = FALSE)
AddQc <- function(q, result, observed, expected) {
  dat.qc <<- rbind(dat.qc, data.frame(q = q, result = result,
      observed = as.character(observed),
      expected = as.character(expected),
      stringsAsFactors = FALSE))
}

vec.soft.avail <- dat.soft$path[file.exists(file.path(dir.data, dat.soft$path))]

for (rel in vec.soft.avail) {
  acc <- sub("_family\\.soft\\.gz$", "", basename(rel))
  par <- ParseFamilySoft(file.path(dir.data, rel))
  list.key <- lapply(seq_along(par$chars), function(i) {
      dat.p <- ParseChars(par$chars[[i]])
      if (nrow(dat.p) == 0L) return(NULL)
      cbind(data.frame(acc = acc, sample_id = par$sample_id[i],
          stringsAsFactors = FALSE), dat.p)
    })
  dat.keys <- rbind(dat.keys, do.call(rbind, list.key))
  i.reg <- match(acc, dat.reg$acc)
  dat.design <- rbind(dat.design, data.frame(
    acc = acc, layer = dat.reg$layer[i.reg], n_samples = length(par$sample_id),
      n_expect = dat.reg$n_expect[i.reg],
      platform = paste(unique(par$platform), collapse = ";"),
      platform_expect = dat.reg$platform[i.reg],
      platform_host = names(sort(table(par$platform), decreasing = TRUE))[1L],
      result = NA_character_, stringsAsFactors = FALSE))
}
## 平台断言按**主体平台**（样本数最多的平台）判定，允许系列内存在少量异平台样本
dat.design$result <- ifelse(dat.design$n_samples != dat.design$n_expect,
  "FAIL_N",
  ifelse(dat.design$platform_host != dat.design$platform_expect,
    "FAIL_PLATFORM", "PASS"))
WriteTsv(dat.design, file.path(dir.tab, "design_check.csv"))
WriteTsv(dat.keys, file.path(dir.tab, "design_keys.tsv"))

## 异平台样本登记：这些样本不是生物样本（如 NTC 无模板对照），须在 S03/S12 剔除
dat.notes <- data.frame(acc = character(0), platform = character(0),
  n_samples = integer(0), sample_id = character(0),
  stringsAsFactors = FALSE)
vec.notes <- character(0)
for (rel in vec.soft.avail) {
  acc <- sub("_family\\.soft\\.gz$", "", basename(rel))
  par <- ParseFamilySoft(file.path(dir.data, rel))
  ## 登记所有**非主体平台**的样本：它们不是生物样本（NTC 之类），下游须剔除。
  ## 判据是与主体平台比对，不是与期望平台比对——后者一旦命中就跳过，会把
  ## 真正该登记的情况漏掉（2026-09-28 首跑 offplatform_samples.tsv 未生成即此因）。
  i.host <- names(sort(table(par$platform), decreasing = TRUE))[1L]
  sel <- par$platform != i.host
  if (sum(sel) == 0L) next
  vec.notes <- c(vec.notes, paste0(acc, " 非主体平台 ", i.host, " 以外样本: ",
      paste(par$sample_id[sel], collapse = ",")))
  dat.notes <- rbind(dat.notes, data.frame(
      acc = acc, platform = paste(unique(par$platform[sel]), collapse = ";"),
    n_samples = sum(sel), sample_id = paste(par$sample_id[sel], collapse = ","),
      stringsAsFactors = FALSE))
  ## 强制登记为 QC 行，避免这条信息只活在 log 里（规则 26）
  AddQc(paste0("NTC_", acc, "_offplatform"), "NOTE",
    paste(sum(sel), "off-host-platform sample(s)"), i.host)
}
if (length(vec.notes) > 0L) {
  cat("[off-platform] ", paste(vec.notes, collapse = " | "), "\n", sep = "")
  WriteTsv(dat.notes, file.path(dir.tab, "offplatform_samples.tsv"))
}

## PART 4 GSE232911 专项断言（§3.2.1 ③ + Q14） ----------------------------------------

if ("GSE232911" %in% dat.design$acc) {
  par <- ParseFamilySoft(file.path(dir.data, "GSE232911_family.soft.gz"))
  dat.chr <- do.call(rbind, lapply(seq_along(par$chars), function(i) {
        d <- ParseChars(par$chars[[i]])
        d$sample_id <- par$sample_id[i]
        d$title <- par$title[i]
        d
      }))

  ## 长表 → 宽表：一个样本一行、每个 characteristics 键成列。
  ## 坑（2026-09-28 实测）：直接 value[match("key", key)] 取回的是**第一个**匹配
  ## 下标，同一值被广播给全部 984 行，交叉表 cells={984}、血栓维度退化成单一取值。
  ## 宽表实测交叉表 = 76/76/34/34/13/13，与下方期望逐字吻合，可反证期望值无误。
  dat.chr.w <- reshape(dat.chr[, c("sample_id", "key", "value")],
    idvar = "sample_id",
    timevar = "key", direction = "wide")
  names(dat.chr.w) <- sub("^value\\.", "", names(dat.chr.w))
  dat.chr.w <- dat.chr.w[match(par$sample_id, dat.chr.w$sample_id), ]
  dat.chr <- data.frame(
    sample_id = dat.chr.w$sample_id,
    title     = par$title,
    state     = dat.chr.w[["disease state"]],
    tunica    = dat.chr.w[["tunica"]],
    thrombus  = dat.chr.w[["thrombus"]],
    stringsAsFactors = FALSE)
  ## 四键缺一即导致交叉表静默错位，此处硬失败（规则 26：禁静默跳过）
  if (sum(complete.cases(dat.chr)) != length(par$sample_id)) {
    stop(paste0("S01: GSE232911 有样本缺 characteristics 键，缺键样本数=",
        length(par$sample_id) - sum(complete.cases(dat.chr))))
  }

  ## (a) 样本总数 = 246
  n.total <- length(unique(dat.chr$sample_id))
  AddQc("Q14a_n_total", if (n.total == 246L) "PASS" else "FAIL",
    n.total, 246L)

  ## (b) 交叉表逐格核对（220 AAA + 26 对照）
  ## 期望八格（2026-09-28 实测反证，非改宽以迁就结果）：
  ##   AAA  220 = 76(Adv x covered) + 76(Media x covered)
  ##              + 34(Adv x free)   + 34(Media x free)
  ##   Donor 26 = 13(Adv x free)    + 13(Media x free)
  ##   余两格（Donor x Thrombus-covered）恒为 0 —— §3.2.1「Donor 全为 Thrombus-free」。
  ## 旧期望只列六个非零数，与八格向量比长度不同，恒不相等（判据本身写错）。
  dat.tab <- table(state = dat.chr$state, tunica = dat.chr$tunica,
    thrombus = dat.chr$thrombus)
  vec.flat <- as.vector(dat.tab)
  vec.flat.expected <- sort(c(76L, 76L, 34L, 34L, 13L, 13L, 0L, 0L))
  AddQc("Q14b_crosstab", if (identical(sort(vec.flat),
        vec.flat.expected)) "PASS" else "FAIL",
    paste0("cells={", paste(vec.flat, collapse = ","), "}"),
    "cells={0,0,13,13,34,34,76,76}")

  ## (c) 对照侧必须全部 Thrombus-free（§3.2.1 ④ 硬约束 1）
  sel.ctrl <- dat.chr$state == "Donor abdominal aorta"
  n.ctrl <- length(unique(dat.chr$sample_id[sel.ctrl]))
  n.ctrl.free <- length(unique(dat.chr$sample_id[sel.ctrl & dat.chr$thrombus ==
        "Thrombus-free"]))
  AddQc("Q14c_control_all_thrombusfree", if (
      n.ctrl == n.ctrl.free) "PASS" else "FAIL",
    paste0(n.ctrl, " control / ", n.ctrl.free, " Thrombus-free"),
    paste0(n.ctrl, " control / ", n.ctrl, " Thrombus-free"))

  ## (d) AAA 侧两个血栓维度都存在
  n.cov <- length(unique(dat.chr$sample_id[
        dat.chr$state == "Abdominal aortic aneurysm" &
        dat.chr$thrombus == "Thrombus-covered"]))
  n.free <- length(unique(dat.chr$sample_id[
        dat.chr$state == "Abdominal aortic aneurysm" &
        dat.chr$thrombus == "Thrombus-free"]))
  AddQc("Q14d_aaa_has_both_thrombus", if (
      n.cov > 0L && n.free > 0L) "PASS" else "FAIL",
    paste0("covered=", n.cov, " free=", n.free), "covered>0 / free>0")

  ## (e) 伪重复：主体 ID 取 title 下划线第 2 段。
  ## 口径实测（2026-09-28）：AAA 侧 `001_3A003_TH_adv` → 3A003；
  ##     Donor 侧 `233_CTRL01_adv` → CTRL01；共 89 主体（76 患者 + 13 供者），
  ##     每个主体 >= 2 份样本，与下方期望逐字吻合。
  ## 旧正则 ^[0-9]+_([A-Za-z0-9]+)_ 会把 `3A003_TH` 整个吃掉而只剩 TH/adv，
  ## 实测主体数塌缩为 6，与 246 样本自相矛盾。
  vec.subj <- vapply(strsplit(dat.chr$title, "_"), function(x) x[2L], "")
  tab.subj <- table(vec.subj)
  n.subj <- length(tab.subj)
  n.dup.ok <- sum(tab.subj >= 2L)
  AddQc("Q14e_subject_pseudorep", if (
      n.subj == 89L && n.dup.ok == n.subj) "PASS" else "FAIL",
    paste0(n.subj, " subjects; ", n.dup.ok, " with >=2 samples"),
    "89 subjects; 89 with >=2 samples")

  ## dat.chr 已是样本级宽表，直接取即可（旧写法在此再套 unique，语义重复）
  dat.group <- dat.chr[, c("sample_id", "title", "state", "tunica", "thrombus")]
  dat.group$subject <- vapply(strsplit(dat.group$title, "_"), function(x) x[2L],
    "")
  WriteTsv(dat.group, file.path(dir.tab, "gse232911_grouping.tsv"))
}

## (f)(g) 矩阵尺度 + 12 基因面板（§3.2.1 ④ 硬约束 2、3；Q4）
mat.rel <- "GSE232911_Matrix.file.txt.gz"
if (file.exists(file.path(dir.data, mat.rel))) {
  dat.mat <- readr::read_tsv(file.path(dir.data, mat.rel),
    show_col_types = FALSE)
  vec.col <- names(dat.mat)
  i.sym <- grep("^Gene.Symbol$", vec.col)
  vec.sym <- as.character(dat.mat[[i.sym[1]]])
  vec.sym.split <- unique(unlist(strsplit(gsub('^"|"$', "", vec.sym), ";")))
  hit <- vec.panel12[vec.panel12 %in% vec.sym.split]
  miss <- setdiff(vec.panel12, hit)
  AddQc("Q14f_panel12_on_matrix", if (length(hit) == 12L) "PASS" else "FAIL",
    paste0(length(hit), "/12 hit; missing=", paste(miss, collapse = ",")),
    "12/12 on Gene.Symbol column")
  
  ## 样本列选择：与 S03 同口径 —— 取「末列向前的最长连续数值列块」，自动跳过前 8 列注释
  ## （4 个数值索引列 ...1/...2/.2/.1 + 4 个文本列 ID/Public.Gene.IDs/Gene.Symbol/Description）。
  ## 根因（2026-09-29 实测）：本文件样本列按**样本标题**命名（001_3A003_TH_adv …），
  ## 不以 GSM 开头；旧写法 grep("^GSM", vec.col) 命中 0 列 → range() 返回 Inf/-Inf
  ## → Q14g 因 "Inf > 0" 误判 PASS（假通过）。此处同时加 is.finite 守卫，杜绝空选择假通过。
  is.num <- vapply(dat.mat, is.numeric, logical(1L))
  k <- length(is.num)
  while (k >= 1L && is.num[k]) k <- k - 1L
  i.num <- (k + 1L):length(is.num)
  mat.val <- as.matrix(dat.mat[, i.num, drop = FALSE])
  vec.range <- range(mat.val, na.rm = TRUE)
  vec.med <- median(mat.val, na.rm = TRUE)
  AddQc("Q14g_matrix_scale",
    if (length(i.num) > 0L && is.finite(vec.range[1L]) && vec.range[1L] > 0) {
      "PASS"
    } else {
      "FAIL"
    },
    paste0("n_sample_col=", length(i.num),
      " min=", format(vec.range[1L], digits = 4),
      " median=", format(vec.med, digits = 4),
      " max=", format(vec.range[2L], digits = 4)),
    "246 sample cols; min > 0  => 已 log2 归一化；S03 禁止二次 log2")
  cat(sprintf(paste0("[info] GSE232911 matrix range %.4f - %.4f (median",
    "%.4f), %d sample cols\n"),
      vec.range[1L], vec.range[2L], vec.med, length(i.num)))
}
## PART 4.5 样本轴重映射表（来源列名 -> GSM）与派生设计字段 ----------------------------------
## 目的：部分系列的表达来源文件**不是** series_matrix（其表头即 GSM），列名是来源文件自带的
##       匿名列名，与 design_keys.tsv 的 sample_id（= GSM）无法按名连接 —— 实测交集 = 0，
##       证据 01_data/verify/check_s05_sample_axis_gap.R。
##       本 PART 从 family SOFT 复原映射，并登记 SOFT 未提供的派生字段（GSE323983 的 donor）。
##
## 注册表：只登记「表达来源列名不是 GSM」的系列（判定依据：S03 的 dat.reg$source 不是
##         series_matrix）。
##   GSE57691  : 来源 = non-normalized_data.txt.gz，列名 "Sample N"；
##               SOFT title = "<group>_Sample N"（实测 68 条全覆盖）
##   GSE323983 : 来源 = 两个 counts tsv，列名 "<arm><donor>_sum"；
##               SOFT title = "HCAEC, <arm>, <age>YO, <sex>"（实测 23 条全 4 段）
vec.axis.remap <- c("GSE57691", "GSE323983")

## GSE323983：SOFT title 里的臂写法 -> counts 列名**实际前缀**。
## 前缀**逐字**取自两个 counts 文件的表头实测（大小写不同，不得统一）：
##   CvIR（file 1）      : `Gene | CTRL1_sum … CTRL6_sum | IR1_sum … IR6_sum`（大写）
##   RapavIRRapa（file 2）: `gene_id | gene_name | rapa1_sum … rapa6_sum |
##                          irrapa1_sum … irrapa6_sum`（小写）
## 2026-09-29 首次实现时统一写成小写，S01 断言 sample_axis_gse323983 立即 FAIL
## （extra_in_counts = CTRL1_sum…IR6_sum 共 11 列），即大小写必须分辨。
dat.arm.title <- data.frame(
  title_arm = c("Ctrl", "Rapa", "IR", "IR with Rapa"),
  prefix = c("CTRL", "rapa", "IR", "irrapa"),
  stringsAsFactors = FALSE)

## 来源文件表头笔误（实测）：counts 文件表头逐字为
##   `... rapa3_sum rara4_sum rapa5_sum ...`，即 rapa4 被写作 rara4。
## 按位置归属 donor 4（= 17YO 的 Rapa 样本），此处显式登记，不静默容错。
vec.axis.typo <- c(rapa4_sum = "rara4_sum")

dat.axis <- data.frame(acc = character(0), source_col = character(0),
  sample_id = character(0), donor = character(0),
  arm = character(0), age_yo = character(0), sex = character(0),
  stringsAsFactors = FALSE)

for (acc in vec.axis.remap) {
  par <- ParseFamilySoft(file.path(dir.data,
      dat.soft$path[dat.soft$acc == acc]))

  if (acc == "GSE57691") {
    ## title 形如 `small AAA_Sample 1` / `Donor_Sample 66`，下划线后即来源列名
    dat.ax <- data.frame(
      acc = acc, source_col = sub("^.*_", "",
        par$title), sample_id = par$sample_id,
      donor = NA_character_, arm = NA_character_, age_yo = NA_character_,
      sex = NA_character_, stringsAsFactors = FALSE)
    ## 断言：映射必须与来源文件表头**同集合**。表头形如
    ##   `ID_REF | Sample 1 | Detection Pval | Sample 2 | …`，剔掉两列非样本列后应余 68。
    con <- gzfile(file.path(dir.data, "GSE57691_non-normalized_data.txt.gz"),
      "rt")
    vec.hd <- strsplit(readLines(con, 1L, warn = FALSE), "\t")[[1L]]
    close(con)
    vec.hd <- setdiff(vec.hd, c("ID_REF", "Detection Pval"))
    AddQc("sample_axis_gse57691",
      if (length(vec.hd) == nrow(dat.ax) &&
        setequal(vec.hd, dat.ax$source_col)) "PASS" else "FAIL",
      sprintf("%d header sample cols vs %d SOFT samples; setequal=%s",
        length(vec.hd), nrow(dat.ax), setequal(vec.hd,
          dat.ax$source_col)),
      "two sets identical; SOFT order == header order")
  } else {
    ## title 形如 `HCAEC, IR with Rapa, 17YO, Female`，按 "," 切成 4 段
    list.prt <- strsplit(par$title, ",")
    vec.arm.txt <- trimws(vapply(list.prt, function(x) x[2L], ""))
    vec.age.txt <- trimws(vapply(list.prt, function(x) x[3L], ""))
    vec.sex.txt <- trimws(vapply(list.prt, function(x) x[4L], ""))
    ## donor 序号 = SOFT 中 donor 的**首次出现序**（54YO 与 54YO(2) 是两位不同供体，
    ## 靠 age 标记区分）。该口径为推断，下一段用「缺 IR 的 donor」做交叉校验。
    vec.donor <- as.character(match(vec.age.txt, unique(vec.age.txt)))
    pref <- dat.arm.title$prefix[match(vec.arm.txt, dat.arm.title$title_arm)]
    if (anyNA(pref)) {
      stop(paste0("S01: GSE323983 title 里出现未登记的臂写法: ",
          paste(unique(vec.arm.txt[is.na(pref)]), collapse = ", ")))
    }
    s.col <- paste0(pref, vec.donor, "_sum")
    i.typo <- match(s.col, names(vec.axis.typo))
    s.col[!is.na(i.typo)] <- unname(vec.axis.typo[i.typo[!is.na(i.typo)]])
    dat.ax <- data.frame(
      acc = acc, source_col = s.col, sample_id = par$sample_id,
      donor = vec.donor,
      arm = dat.arm.title$title_arm[match(pref, dat.arm.title$prefix)],
      age_yo = vec.age.txt, sex = vec.sex.txt, stringsAsFactors = FALSE)

    ## 断言 1：映射集合必须与两个 counts 文件的表头**逐字同集合**
    vec.hd2 <- character(0)
    for (fl in c("GSE323983_CvIR_all_gene_counts.tsv.gz",
        "GSE323983_RapavIRRapa_all_gene_counts.tsv.gz")) {
      con <- gzfile(file.path(dir.data, fl), "rt")
      v <- strsplit(readLines(con, 1L, warn = FALSE), "\t")[[1L]]
      close(con)
      vec.hd2 <- c(vec.hd2, setdiff(v, c("Gene", "gene_id", "gene_name")))
    }
    AddQc("sample_axis_gse323983",
      if (length(vec.hd2) == nrow(dat.ax) &&
        setequal(vec.hd2, dat.ax$source_col)) "PASS" else "FAIL",
      sprintf(
        "%d counts cols vs %d SOFT samples; setequal=%s; extra_in_counts=%s",
        length(vec.hd2), nrow(dat.ax), setequal(vec.hd2,
          dat.ax$source_col),
        paste(setdiff(vec.hd2, dat.ax$source_col), collapse = ",")),
      "23 cols both sides; ctrl/ir/rapa/irrapa x donor")
    ## 断言 2（交叉校验，本映射为推断必须自证）：SOFT 里缺 IR 样本的 donor
    ##   必须与 CvIR 文件里缺失的序号一致 —— 实测两侧都指向 donor 4（17YO）。
    d.miss.soft <- setdiff(sort(unique(as.integer(vec.donor))),
      sort(unique(as.integer(vec.donor[vec.arm.txt == "IR"]))))
    d.miss.file <- setdiff(1:6, as.integer(sub("^IR([0-9]+)_sum$", "\\1",
          grep("^IR[0-9]+_sum$", vec.hd2,
            value = TRUE))))
    AddQc("sample_axis_gse323983_ir_missing_xcheck",
      if (identical(d.miss.soft,
          d.miss.file) && length(d.miss.soft) == 1L) "PASS"
      else "FAIL",
      sprintf("SOFT lacks IR for donor {%s}; CvIR lacks {%s}",
        paste(d.miss.soft, collapse = ","), paste(d.miss.file,
          collapse = ",")),
      "both = donor 4 (17YO), the only donor without an IR sample")
  }
  dat.axis <- rbind(dat.axis, dat.ax)
}
WriteTsv(dat.axis, file.path(dir.tab, "sample_axis.tsv"))
cat(sprintf("[sample_axis] rows=%d  series=%s\n", nrow(dat.axis),
    paste(unique(dat.axis$acc), collapse = ", ")))


## PART 5 Q7 负断言：GSE7084 的样本号不得出现在任何解析结果中 --------------------------

vec.hit7084 <- character(0)
if ("GSE7084" %in% dat.design$acc) {
  par <- ParseFamilySoft(file.path(dir.data, "verify/GSE7084_family.soft.gz"))
  vec.hit7084 <- intersect(par$sample_id, sprintf("GSM%d", 170550:170568))
}
AddQc("Q7_gse7084_negative", if (length(vec.hit7084) == 0L) "PASS" else "FAIL",
  paste0("hits=", paste(vec.hit7084,
      collapse = ",")), "0 hit (GSE7084 已完全剔除)")
if (length(vec.hit7084) > 0L) {
  stop(paste0("Q7 FAILED: GSE7084 sample(s) present: ", paste(vec.hit7084,
        collapse = ","),
      " ; must not enter any extraction list."))
}

## PART 6 臂长登记（供 Q12：n < 10 的臂不进主聚合） ------------------------------------

dat.arm <- data.frame(acc = dat.design$acc, layer = dat.design$layer,
  group = "ALL",
  n_samples = dat.design$n_samples,
  derivation = "family SOFT 样本总数（S01 实测）",
  stringsAsFactors = FALSE)
if (exists("dat.group")) {
  tab.arm <- table(group = paste(dat.group$state, dat.group$tunica,
      dat.group$thrombus, sep = "|"))
  dat.arm.g <- data.frame(acc = "GSE232911", layer = "L-D", group = names(
      tab.arm),
    n_samples = as.integer(tab.arm),
    derivation = "characteristics: disease state + tunica + thrombus",
    stringsAsFactors = FALSE)
  dat.arm <- rbind(dat.arm[dat.arm$acc != "GSE232911", ], dat.arm.g)
}
WriteTsv(dat.arm, file.path(dir.tab, "arm_sizes.csv"))

## PART 6.5 图件 ------------------------------------------------------------------------
## 命名遵守规范（描述性英文，禁 Fig 前缀；图内文字全英文）；不引入拼图包，单面板出图。
## 两幅图分别对应样本量构成与设计断言（见文件头图件清单 / 规则 35）。

dat.design.plot <- dat.design[!is.na(dat.design$result), ]
dat.design.plot$layer <- factor(dat.design.plot$layer,
  levels = c("L-C", "L-P", "L-D"))
## x 轴按「层 + 样本数」排序：同层数据集相邻、柱高单调，便于肉眼核对。
dat.design.plot$acc <- factor(
  dat.design.plot$acc,
  levels = dat.design.plot$acc[order(dat.design.plot$layer,
      dat.design.plot$n_samples)])

## 图 1：各系列实解析样本数（单面板，按层着色）。
## 替代原 facet_wrap(~ layer) 的宽分面（分面后横跨多面板、信息被稀释）。
p.comp <- ggplot(dat.design.plot, aes(x = acc, y = n_samples, fill = layer)) +
geom_col() +
labs(title = "Samples parsed per series, grouped by layer",
  x = "Accession", y = "Samples parsed", fill = "Layer") +
sci_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
sci_save(p.comp, "Sample_Composition_By_Layer", width = 9, height = 5)

## 图 2：实测 vs 预期样本数（单面板散点 + y=x 参考线，按层着色）。
## 替代原 facet_wrap(~ layer)：17 个点分三面板后每面板信息过少，单面板更利于判读一致性。
dat.count.chk <- data.frame(
  layer = dat.design.plot$layer,
  acc = dat.design.plot$acc,
  expected = dat.design.plot$n_expect,
  observed = dat.design.plot$n_samples)
p.chk <- ggplot(dat.count.chk, aes(x = expected, y = observed, color = layer)) +
geom_point(size = 2.5) +
geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
labs(title = "Observed versus expected sample counts",
  x = "Expected samples (protocol registration)",
  y = "Observed samples (family SOFT)", color = "Layer") +
sci_theme
sci_save(p.chk, "Sample_Count_vs_Expected", width = 7.5, height = 5)

## PART 7 收口 -----------------------------------------------------------------------


WriteTsv(dat.qc, file.path(dir.tab, "gse232911_qc.csv"))

n.fail <- sum(dat.design$result != "PASS")
n.qc.fail <- sum(dat.qc$result == "FAIL")
n.missing <- sum(dat.manifest$status == "MISSING")

message(paste0("\n===== S01 ingest summary =====\nfiles registered = ", nrow(
      dat.manifest),
    "  missing = ", n.missing,
    "  excluded-on-disk = ", sum(dat.manifest$status == "EXCLUDED"),
    "\nseries parsed = ", nrow(dat.design),
    "  design FAIL = ", n.fail,
    "  GSE232911 QC FAIL = ", n.qc.fail))
if (n.fail > 0L) print(dat.design[dat.design$result != "PASS", ])
if (n.qc.fail > 0L) print(dat.qc[dat.qc$result == "FAIL", ])
if (n.qc.fail > 0L || n.fail > 0L) {
  stop(paste0("S01 assertions FAILED: design_fail=", n.fail, " qc_fail=",
      n.qc.fail,
      " ; see 03_results/S01/table/ ."))
}
message("S01 OK. tables -> 03_results/S01/table/")
