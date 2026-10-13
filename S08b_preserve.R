# ===============================================================================================
# 功能：跨系列 WGCNA 模块保守性（modulePreservation）
# 定位：S08 补充分析。GSE57691 / GSE43292 未满足 scale-free 拓扑（见
#       wgcna_scale_free_not_met.csv），不硬建网络；改以已建网络的 GSE232911 /
#       GSE323983 模块为参考，评估其在两失败系列中的保守性（Langfelder 2011 官方函数）。
# 可视化：参考模块 x 测试系列 的 Zsummary 保守性热图（SCI 出版导向，规则 35）。
# 方法依据：WGCNA::modulePreservation（PMID 21910903）。Zsummary > 2 保守好、
#       > 0 弱保守、< 0 不保守（该文判定阈值）。
# ===============================================================================================

library(WGCNA)   # 共表达网络与模块保守性官方函数
library(ggplot2)  # 热图可视化
disableWGCNAThreads()   # 本机 SnowParam 14 工并行比串行慢 23-38x，强制单线程（规则 33 性能）

## ---- 常量（顶部集中，规则 28）----
SEED <- 20260925
TAB_DIR <- "03_results/S08/table"
FIG_DIR <- "03_results/S08/Figure"
CKPT_DIR <- file.path(TAB_DIR, "checkpoints")
MAX_WGCNA_GENES <- 5000L
NETWORK_TYPE <- "unsigned"
N_PERM <- 200L   # Zsummary 为置换统计量，须 nPermutations>=200（Langfelder 2011 判定阈值）
vec.ref <- c("GSE232911", "GSE323983")    # 已建网络的参考系列
vec.test <- c("GSE57691", "GSE43292")     # scale-free 失败的测试系列

dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

## 主题/出口：规范 §sci_theme / §sci_save 官方定义（与 S08 逐字一致，禁 source）
sci_theme <- ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
    legend.position = "right",
    legend.text = ggplot2::element_text(size = 10, face = "plain"),
    axis.title = ggplot2::element_text(size = 12, face = "plain"),
    axis.text = ggplot2::element_text(size = 10, face = "plain"),
    plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "plain"),
    legend.title = ggplot2::element_text(size = 11, face = "plain"))
sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  pdf(pdf_file, width = width, height = height)
  print(plot_obj)
  dev.off()
}

## ---- 加载基因级表达矩阵（genes x samples，按 S08 口径取 top SD 基因）----
LoadExpr <- function(acc) {
  obj <- readRDS(file.path("03_results", "S03", "table",
    paste0(acc, "_gene.rds")))
  expr <- obj$expr
  if (nrow(expr) > MAX_WGCNA_GENES) {
    vec.sd <- apply(expr, 1L, sd)
    vec.sd <- vec.sd[!is.na(vec.sd)]
    expr <- expr[names(sort(vec.sd, decreasing = TRUE))[
      seq_len(min(MAX_WGCNA_GENES, length(vec.sd)))], , drop = FALSE]
  }
  expr
}

## ---- 参考模块颜色（基因 -> 模块，含 ME0 grey）----
ref.colors <- list()
for (acc in vec.ref) {
  df <- read.csv(file.path(TAB_DIR, paste0("wgcna_modules_", acc,
    ".csv")), stringsAsFactors = FALSE)
  ref.colors[[acc]] <- setNames(df$module, df$gene)
}

## ---- PART 1 逐参考系列 x 测试系列 计算模块保守性（PART 级 checkpoint，规则 33）----
dat.pres <- data.frame()
for (acc.ref in vec.ref) {
  expr.ref <- LoadExpr(acc.ref)
  col.ref <- ref.colors[[acc.ref]]
  for (acc.test in vec.test) {
    f.out <- file.path(TAB_DIR, paste0("module_preservation_",
      acc.ref, "_in_", acc.test, ".csv"))
    f.ck <- file.path(CKPT_DIR, paste0("S08b_pres_", acc.ref,
      "_", acc.test, "_p", N_PERM, ".rds"))
    if (file.exists(f.ck)) {
      cat(sprintf("[S08b] %s in %s: checkpoint reused\n", acc.ref, acc.test))
      res <- readRDS(f.ck)
    } else {
      ## 共同基因（modulePreservation 按基因名对齐模块成员）
      expr.test <- LoadExpr(acc.test)
      vec.gene <- intersect(intersect(rownames(expr.ref), names(col.ref)),
        rownames(expr.test))
      stopifnot(length(vec.gene) > 0L)
      dat.ref <- t(expr.ref[vec.gene, , drop = FALSE])    # samples x genes
      dat.test <- t(expr.test[vec.gene, , drop = FALSE])
      col.ref.com <- col.ref[vec.gene]
      col.test <- ifelse(vec.gene %in% names(col.ref),
        col.ref[vec.gene], "ME0")   # 测试集基因投影参考色，缺失->grey
      multiExpr <- list(list(data = dat.ref), list(data = dat.test))
      multiColor <- list(col.ref.com, col.test)
      cat(sprintf("[S08b] modulePreservation %s (ref) in %s (test), %d genes\n",
        acc.ref, acc.test, length(vec.gene)))
      res <- modulePreservation(multiExpr, multiColor, reference = 1,
        test = 2, networkType = NETWORK_TYPE, nPermutations = N_PERM,
        randomSeed = SEED)
      saveRDS(res, f.ck)   # 规则 33：PART 级 checkpoint，失败仍 stop
    }
    ## 提取保守性统计量（解析 Zsummary.pres / medianRank.pres）
    ## 官方返回结构（实测本机 WGCNA 版本）：preservation$Z 为按参考集索引的列表
    ## （元素名 ref.Set_k），每个参考集下为按测试集索引的比较列表；参考集自身的
    ## inColumnsAlsoPresentIn.Set_k 为 NA，测试集对应元素为 data.frame，Zsummary
    ## 列名实为 Zsummary.pres、medianRank 列名实为 medianRank.pres（规则 23：以
    ## 实际返回结构为准，不按名硬编码）。按「非 NA 的 data.frame 元素」稳健提取。
    lst.z <- res$preservation$Z[[1L]]
    lst.obs <- res$preservation$observed[[1L]]
    idx.z <- which(sapply(lst.z, is.data.frame))
    idx.obs <- which(sapply(lst.obs, is.data.frame))
    stopifnot(length(idx.z) == 1L, length(idx.obs) == 1L)
    mat.z <- lst.z[[idx.z]]
    mat.obs <- lst.obs[[idx.obs]]
    stopifnot(all(rownames(mat.z) == rownames(mat.obs)))
    df.mod <- data.frame(reference = acc.ref, test = acc.test,
      module = rownames(mat.z),
      Zsummary = mat.z[, "Zsummary.pres"],
      medianRank = mat.obs[, "medianRank.pres"],
      stringsAsFactors = FALSE)
    dat.pres <- rbind(dat.pres, df.mod)
    write.csv(df.mod, f.out, row.names = FALSE)
    cat(sprintf("  %s in %s: %d modules -> %s\n", acc.ref, acc.test,
      nrow(df.mod), f.out))
  }
}
write.csv(dat.pres, file.path(TAB_DIR, "module_preservation_all.csv"),
  row.names = FALSE)

## ---- PART 2 保守性热图（参考模块 x 测试系列，Zsummary 着色）----
dat.pres$label <- paste0(dat.pres$reference, ":", dat.pres$module)
sci_save(ggplot(dat.pres, aes(x = test, y = label, fill = Zsummary)) +
  geom_tile(colour = "grey80", size = 0.2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
    midpoint = 0, name = "Zsummary") +
  labs(title = "Cross-series module preservation (Zsummary)",
    x = "Test series", y = "Reference module") +
  theme(axis.text.y = element_text(size = 6)) +
  sci_theme,
  "ModulePreservation_Zsummary", width = 7, height = 10)

cat(sprintf("\n[S08b] done: %d module-test pairs, figure ModulePreservation_Zsummary.pdf\n",
  nrow(dat.pres)))
