#===============================================================
## S00_setup.R  —— 环境与接口准入（步骤规格 §8：S00 环境与接口基线）
## 目的：把「本机 R 环境 + 本项目所需 R 包」逐条实测登记，产出 install_log.csv，
##       并作为下游步骤（S01 起）的**准入门禁**：必备包缺失时 stop，下游不得启动。
## 定位：02_scripts/ 下的分析脚本；只做登记与门禁，不做任何生物统计、不安装包。
## 可视化：无。
## 运行：D:/R/R-4.6.1/bin/Rscript.exe --vanilla 02_scripts/S00_setup.R [--check]
##       --check 兼容保留（本脚本已无安装分支，--check 与不传行为一致：仅核对门禁）。
## 产出：03_results/S00/table/install_log.csv（环境基线，含每个包实测版本与状态）
##       退出码 1 表示必备包缺失或门禁失败，下游步骤不得启动。
## 整合说明（2026-10-03）：原 S00* 共 6 个脚本，其中 5 个（S00b_install_github /
##   S00b_install_missing / S00c_restore_immunedeconv / S00d_install_gh_deps /
##   S00e_reinstall_immunedeconv）纯做「装包 / 还原 / 重装」，属一次性环境构建，
##   不进入分析流水线，已移至 02_scripts/_archive_S00_install_2026-10-03/ 备查。
##   本文件保留并整合的是「数据流必须依托」的步骤：包清单登记 + 环境门禁 +
##   基线落盘。安装逻辑如确需复用，从 _archive 取历史脚本（其联网/RTools 坑已实测）。
## 凭据：包清单与来源（CRAN / Bioconductor）为方案 §6 与 §0.2 的登记值；
##       每个包的**实际版本**由 installed.packages() 实测，不采信清单标注。
#===============================================================

## PART 0 环境与常量 ------------------------------------------------------------

args.cmd <- commandArgs(trailingOnly = TRUE)
flag.check <- "--check" %in% args.cmd

dir.log <- "03_results/S00/table"
dir.create(dir.log, showWarnings = FALSE, recursive = TRUE)
file.log <- file.path(dir.log, "install_log.csv")

## 包清单与来源写死，不做运行时推断。
## （教训：曾用 BiocManager::available() 判来源，它返回「带 names 的字符向量」，
##   对它取 rownames() 得到 NULL，Bioc 包被误判为 CRAN 而安装失败。）
vec.bioc.core <- c(
  "limma", "edgeR", "preprocessCore",                  ## 差异表达与归一化
  "org.Hs.eg.db", "clusterProfiler", "enrichplot", "fgsea", "GSVA", "msigdbr",
  "sva", "RUVSeq", "variancePartition",                ## 批次与混杂因子
  "metafor", "RankProd", "RobustRankAggreg", "decoupleR",  ## 元分析与聚合
  "WGCNA",                                             ## 共表达网络
  "affy", "oligo",                                     ## 阵列读取后端
  ## HTA-2.0 / Focus / 133A2 CDF
  "pd.hugene10kplus", "pd.hgfocus", "pd.hgu133a2",
  "hgu133plus2.db", "hgu133a2.db"
)
vec.cran.core <- c(
  "tibble", "readr", "broom", "checks", "optparse", "foreach", "doParallel",
  "ggplot2", "ggpubr", "cowplot", "patchwork", "RColorBrewer", "viridis",
  "circlize", "ComplexHeatmap"
)
vec.bioc.ext <- c("Seurat", "SeuratObject", "harmony", "scDblFinder", "SingleR",
  "scater", "scuttle", "BiocNeighbors",
  "immunedeconv", "quantiseq", "EPIC", "BayesPrism",
  "Scissor", "oncoPredict", "PharmacoGx",
  ## S09a 读斑块 atlas h5ad 的通道：初版走 zellkonverter 的
  ## reader="R"（实测其 .native_reader 忽略 anndata /raw），
  ## 现改为 rhdf5 直读 /raw/X，两包都须登记（均为 Bioconductor）。
  "rhdf5", "zellkonverter")
## timeROC（时依 ROC，生存验证层备料，2026-09-28 装妥 0.4.1）：
## 原 "timerc" 查无此包（六个仓库 + CRAN 均无，2026-09-27 实测），已删；
## timeROC 属 CRAN，须列在 vec.cran.ext 而非 vec.bioc.ext。
vec.cran.ext <- c("xgboost", "mlr3", "mlr3learners", "shapviz", "glmnet",
  "timeROC")

vec.core <- c(vec.bioc.core, vec.cran.core)
vec.all <- c(vec.core, vec.bioc.ext, vec.cran.ext)
vec.source.of <- c(rep("Bioconductor", length(vec.bioc.core)),
  rep("CRAN", length(vec.cran.core)),
  rep("Bioconductor", length(vec.bioc.ext)),
  rep("CRAN", length(vec.cran.ext)))
names(vec.source.of) <- vec.all
vec.target <- vec.all
vec.source.of <- vec.source.of[vec.target]

## 门禁范围：只卡「S01 起真正被脚本载入的包」（base R + tools + readr）。
## 其余 core/ext 包由用户按需补装，缺失只登记为 PENDING，不阻断下游
## （2026-09-26 用户明确：缺失的包后补，不希望 S00 卡住流程）。
vec.gate.core <- c("readr")

## 已装基线只在文件存在**且表结构匹配当前版本**时并入。
## 旧版脚本写过的 install_log.csv 表结构已变（core -> layer），混入会 rbind 失败；
## 此处显式识别并忽略，同时打印提示，不做静默丢弃。
vec.schema <- c("package", "version", "source", "status", "layer", "seconds")
list.log.prev <- list()
if (file.exists(file.log)) {
  dat.prev <- read.csv(file.log, stringsAsFactors = FALSE)
  if (setequal(names(dat.prev), vec.schema)) {
    list.log.prev <- list(dat.prev)
  } else {
    message(paste0("[WARN] existing install_log.csv has schema ",
        paste(names(dat.prev), collapse = "/"),
        " ; superseded script version -> ignored, will be overwritten."))
  }
}

## PART 1 逐包核对 / 登记（只读，不安装） ---------------------------------------

vec.missing.core <- character(0)
vec.missing.ext <- character(0)
list.log.new <- list()

for (p in vec.target) {
  src <- unname(vec.source.of[p])
  ok.now <- p %in% rownames(installed.packages())
  tim.start <- Sys.time()
  status <- if (ok.now) "OK" else "MISSING"
  dat.new <- data.frame(
    package = p,
    version = if (ok.now) as.character(packageVersion(p)) else NA_character_,
    source = src,
    status = status,
    layer = if (p %in% vec.core) "core" else "extend",
    seconds = round(as.numeric(difftime(Sys.time(), tim.start, units = "secs")),
      1)
  )
  list.log.new[[length(list.log.new) + 1L]] <- dat.new
  if (status == "MISSING") {
    if (p %in% vec.gate.core) vec.missing.core <- c(vec.missing.core, p)
    else vec.missing.ext <- c(vec.missing.ext, p)
  }
}
dat.log <- do.call(rbind, c(list.log.prev, list.log.new))
write.csv(dat.log, file.log, row.names = FALSE)

## PART 2 门禁 -------------------------------------------------------------------

message(paste0("\n===== S00 environment gate =====\nmode = ",
    if (flag.check) "check-only" else "gate",
    "   packages checked = ", length(vec.target),
    "  gate-core missing = ", length(vec.missing.core),
    "  other missing = ", length(vec.missing.ext)))
if (length(vec.missing.core) > 0) {
  message(paste0("[GATE CORE MISSING] ", paste(vec.missing.core,
        collapse = ", ")))
}
if (length(vec.missing.ext) > 0) {
  message(paste0("[PENDING, user installs] ", paste(vec.missing.ext,
        collapse = ", ")))
}

msg.start <- paste0("R ", R.version.string, " | ",
  if (file.exists(file.log)) "baseline recorded" else "no baseline")
if (length(vec.missing.core) > 0) {
  stop(paste0("gate-core package(s) missing: ",
      paste(vec.missing.core, collapse = ", "),
      " ; S01 cannot start. Install them first. [", msg.start, "]"))
}
message(paste0("gate passed. log -> ", file.log))
