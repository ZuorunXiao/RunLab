# =============================================================================
# 00_build_reference.R
# 功能: 读取单细胞参考 h5ad, 构建 label transfer 参考对象 + ssGSEA 基因集
# 定位: Wu et al. Cancer Discovery 2022 ST 复现流水线 Step 0 (纯 R, 无 Python)
# 可视化: 无 (产出 RDS/CSV 供 02/03 复用)
# 方法: anndataR::read_h5ad(as="InMemoryAnnData") 读 h5ad -> 抽取稀疏表达矩阵与
#       obs(main_cell_type/sub_cell_type) -> CreateSeuratObject -> LogNormalize ->
#       FindVariableFeatures -> ScaleData(仅可变基因) -> RunPCA(RNA) ->
#       SCTransform(RNA) -> RunPCA(SCT, 供 03 SCT label transfer) ->
#       导出 ref_maincell_seurat.rds (含 main+sub 注释 / RNA data / SCT data+PCA) 与 gene_sets.csv
# 依赖(D-B): 本脚本不重新执行论文 scRNA 注释, 仅作 loader。sc.h5ad 必须是论文的
#       CD45+ 注释参考(Harmony 整合 + SingleR + 手动 marker, 论文:1740-1746; 仅保留
#       CD45+ clusters, :1741)。若 h5ad 未经该流程注释, 整条参考链失真, 需先核实出处。
# 幂等: 若 ref_maincell_seurat.rds 已存在, 跳过 h5ad 读取直接加载并构建基因集
# 数据: 输入 raw/reference_signatures/sc.h5ad; 输出 results/00_build_reference/
# 参考: Seurat 官方 FindTransferAnchors/TransferData(SCTransfer); GSVA ssGSEA
# 版本(D-G): 论文 Seurat V3.2.2; 本机 Seurat 为 V5, SCTransfer 行为基本一致但结果
#       可能略有差异, 关键步骤已按 V5 API 校准。
# =============================================================================

# --- 环境(限线程避免 BLAS 段错误, 退出跳过 finalizer) -----------------------
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(stringsAsFactors = FALSE)
suppressMessages({
  library(anndataR)   # 纯 R 读取 h5ad (rhdf5, 无 Python)
  library(Seurat)     # 参考对象与 PCA
  library(Matrix)     # 稀疏矩阵
})

# --- 目录结构 ---------------------------------------------------------------
BASE_DIR <- "G:/CRLM/ST/scCRLM"
REF_DIR  <- file.path(BASE_DIR, "results/00_build_reference")
H5AD     <- file.path(BASE_DIR, "raw/reference_signatures/sc.h5ad")
REF_RDS  <- file.path(REF_DIR, "ref_maincell_seurat.rds")
dir.create(REF_DIR, recursive = TRUE, showWarnings = FALSE)

META_COLS <- c("main_cell_type", "sub_cell_type", "patient", "tissue",
               "doublet", "percent.mt")

# =============================================================================
# PART 1-3: 读取 h5ad 并构建参考对象 (幂等: RDS 已存在则直接加载)
# 目的: 获得含 main/sub 注释与 PCA 的参考 Seurat 对象
# 方法: 见文件头; 仅对可变特征 ScaleData 以省内存(避免 ~12GB 密集化)
# =============================================================================
if (file.exists(REF_RDS)) {
  cat("========== 检测到已有参考 RDS, 跳过 h5ad 读取 (幂等) ==========\n")
  ref <- readRDS(REF_RDS)
  genes    <- rownames(ref)
  meta_keep <- ref@meta.data[, META_COLS]
  cat("已加载参考:", ncol(ref), "cells x", nrow(ref), "genes; PCA:", ncol(ref[["pca"]]), "\n")
} else {
  cat("========== PART 1: 读取参考 h5ad ==========\n")
  ad <- read_h5ad(H5AD, as = "InMemoryAnnData")
  genes <- ad$var_names
  meta  <- ad$obs
  cat("cells:", nrow(meta), " genes:", length(genes), "\n")
  cat("main_cell_type:", paste(sort(unique(meta$main_cell_type)), collapse = "/"), "\n")
  cat("sub_cell_type 数量:", length(unique(meta$sub_cell_type)), "\n")

  # dgRMatrix(cells x genes) 直接转置为 dgCMatrix(genes x cells), 稀疏->稀疏
  mat <- t(ad$X)
  rownames(mat) <- genes
  colnames(mat) <- rownames(meta)
  rm(ad); gc()
  cat("表达矩阵:", nrow(mat), "genes x", ncol(mat), "cells, nnz:", sum(mat@x != 0), "\n")
  meta_keep <- meta[, META_COLS]

  cat("\n========== PART 2: 构建参考 Seurat 对象 ==========\n")
  ref <- CreateSeuratObject(counts = mat, meta.data = meta_keep, assay = "RNA")
  rm(mat); gc()
  # 内存限制: 全量参考 178630 cells 做 SCTransform 会申请 ~12Gb 稠密残差矩阵而 OOM。
  # label transfer 参考只需各主型有代表性抽样(标准做法, 不改变方法学), 故按
  # main_cell_type 下采样到每型 max.ref.cells 个; 残差矩阵随之从 ~12Gb 降到 ~2Gb。
  max.ref.cells <- 3000
  set.seed(42)
  cells.by.type <- split(colnames(ref), ref$main_cell_type)
  keep.cells <- unlist(lapply(cells.by.type, function(x) sample(x, min(length(x), max.ref.cells))))
  ref <- subset(ref, cells = keep.cells)
  cat("下采样后参考:", ncol(ref), "cells (每主型 <=", max.ref.cells, ")\n")
  ref <- NormalizeData(ref, verbose = FALSE)
  ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  # 仅对可变特征做 ScaleData: 密集 scale.data 为 27k x 2000 (~0.4GB), 避免全基因 OOM
  ref <- ScaleData(ref, features = VariableFeatures(ref), verbose = FALSE)
  ref <- RunPCA(ref, features = VariableFeatures(ref), npcs = 50, verbose = FALSE)
  # SCT 表示供 03 的 SCT label transfer(论文 ST 侧 SCTransform, :1841-1843;
  # Seurat SCTransfer 要求参考与查询同为 SCT, 故参考也做 SCTransform)
  ref <- SCTransform(ref, assay = "RNA", verbose = FALSE)
  ref <- RunPCA(ref, assay = "SCT", npcs = 50, verbose = FALSE)  # 覆盖为 SCT PCA, 供 transfer
  # 保留 counts + RNA data(ssGSEA/基因集用) + SCT(data/PCA 供 transfer),
  # 仅丢弃 RNA scale.data (~2.8GB) 释放内存
  LayerData(ref, layer = "scale.data") <- NULL; gc()
  cat("参考 SCT PCA:", ncol(ref[["pca"]]), " PCs; 细胞:", ncol(ref), "\n")

  cat("\n========== PART 3: 导出参考对象 RDS ==========\n")
  saveRDS(ref, REF_RDS)
  cat("Saved ref_maincell_seurat.rds\n")
}
# 依赖校验(D-B): 参考必须是论文 CD45+ 注释结果, 本脚本仅消费 h5ad 已存注释
if (!"main_cell_type" %in% colnames(ref@meta.data))
  stop("ref 缺少 main_cell_type: 请确认 sc.h5ad 是论文 CD45+ 注释参考(:1740-1746)")
cat("参考 main_cell_type:", paste(sort(unique(ref$main_cell_type)), collapse = "/"), "\n")

# =============================================================================
# PART 4: 构建 ssGSEA 基因集 (按 sub_cell_type 取高表达标志基因)
# 目的: 供 03 计算各 spot 的巨噬细胞等亚群富集评分
# 方法: 对每个 sub_cell_type, 取 (群内均值 - 全局均值) 最高的 top 100 基因
# 数据: 基于 ref 的 log-normalized data 层 (稀疏), 向量化避免密集化
# =============================================================================
build_signature_genes <- function(obj, ident, topn = 100) {
  d  <- GetAssayData(obj, layer = "data")       # genes x cells, log-normalized 稀疏
  cl <- obj@meta.data[[ident]]
  ov <- Matrix::rowMeans(d)
  uniq <- sort(unique(cl))
  out <- list()
  for (c in uniq) {
    j <- which(cl == c)
    if (length(j) < 20) next
    m  <- Matrix::rowMeans(d[, j, drop = FALSE])
    sc <- order(m - ov, decreasing = TRUE)
    out[[c]] <- rownames(d)[sc[1:min(topn, length(sc))]]
  }
  out
}
gene_sets <- build_signature_genes(ref, "sub_cell_type", topn = 100)
cat("基因集数量:", length(gene_sets), "\n")
gs_df <- data.frame(
  gene_set = rep(names(gene_sets), lengths(gene_sets)),
  gene     = unlist(gene_sets),
  stringsAsFactors = FALSE
)
write.csv(gs_df, file.path(REF_DIR, "gene_sets.csv"), row.names = FALSE)
cat("Saved gene_sets.csv (", nrow(gs_df), " 行)\n")

# =============================================================================
# PART 5: 导出基因列表与汇总 RDS
# =============================================================================
writeLines(genes, file.path(REF_DIR, "genes.txt"))
ref_data <- list(genes = genes, meta = meta_keep,
                 main_types = sort(unique(meta_keep$main_cell_type)),
                 sub_types  = sort(unique(meta_keep$sub_cell_type)))
saveRDS(ref_data, file.path(REF_DIR, "reference_data.rds"))
cat("Saved genes.txt + reference_data.rds\n")

cat("\n00_build_reference.R 完成 (纯 R, 无 Python)\n")
# 退出时跳过 finalizer, 避免 rhdf5 退出期段错误
quit(save = "no", status = 0, runLast = FALSE)
