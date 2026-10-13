# ================================================================
# 02b_GSE178318_Extract_Malignant_Counts.R
# 功能: 从 Step 1 注释 Seurat 对象中提取恶性细胞 counts 矩阵
#       为 Step 3 (scTenifoldKnk 虚拟敲除) 准备输入数据
# 定位: Step 2 → Step 3 之间的桥梁脚本
# 输入: seurat_obj_annotated.rds (Step 1 输出, 含 malignant/tissue 列)
# 输出: count_tumor_PT.rds  (恶性细胞 PT  counts 矩阵)
#        count_tumor_LM.rds  (恶性细胞 LM  counts 矩阵)
# ================================================================

# --- 包加载 --------------------------------------------------------
library(Seurat)      # Seurat v5 — 数据读取与子集提取
library(Matrix)       # 稀疏矩阵操作

# --- 目录结构定义 ------------------------------------------------
INPUT_DIR  <- "X:/CRLM/XX/01_QC/table"
OUTPUT_DIR <- "X:/CRLM/XX/03_GRN_VirtualKO/table"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ================================================================
# PART 1: 加载 Step 1 注释 Seurat 对象
# ================================================================

seurat_file <- file.path(INPUT_DIR, "seurat_obj_annotated.rds")
if (!file.exists(seurat_file)) {
  stop(sprintf("ERROR: 输入文件不存在: %s\n", seurat_file),
       "请确保 Step 1 (01_GSE178318_Preprocessing_QC.R) 已成功运行。")
}

seurat.obj <- readRDS(seurat_file)
gc()

# Seurat v5: 合并分片 layers, 确保 "counts" layer 可访问
seurat.obj <- JoinLayers(seurat.obj)
gc()

# ================================================================
# PART 2: 验证元数据完整性
# ================================================================


# 检查必要列
required_cols <- c("malignant", "tissue")
for (col in required_cols) {
  if (!col %in% colnames(seurat.obj@meta.data)) {
    stop(sprintf("ERROR: 缺少元数据列 '%s'", col),
         "请确认 Step 1 已运行 InferCNV 恶性鉴定和组织标注。")
  }
}

# 输出基本信息
cat(sprintf("  总细胞数: %d\n", ncol(seurat.obj)))
cat(sprintf("  组织分布: %s\n",
    paste(names(table(seurat.obj$tissue)),
          sprintf("(%d)", table(seurat.obj$tissue)), collapse = ", ")))
cat(sprintf("  恶性标注: %s\n",
    paste(names(table(seurat.obj$malignant)),
          sprintf("(%d)", table(seurat.obj$malignant)), collapse = ", ")))

# ================================================================
# PART 3: 筛选恶性细胞并拆分组织
# ================================================================


# 筛选恶性细胞
malignant_barcodes <- WhichCells(seurat.obj, expression = malignant == "malignant")

# 按组织拆分
pt_barcodes <- intersect(
  malignant_barcodes,
  WhichCells(seurat.obj, expression = tissue == "PT")
)
lm_barcodes <- intersect(
  malignant_barcodes,
  WhichCells(seurat.obj, expression = tissue == "LM")
)

cat(sprintf("  恶性细胞总数: %d\n", length(malignant_barcodes)))
cat(sprintf("    PT (原发灶): %d\n", length(pt_barcodes)))
cat(sprintf("    LM (肝转移): %d\n", length(lm_barcodes)))

# ================================================================
# PART 4: 提取 counts 矩阵并保存
# ================================================================


# Seurat v5: 从 "counts" layer 提取原始 UMI 计数
extract_counts <- function(obj, cells, label) {
  if (length(cells) == 0) {
    stop(sprintf("ERROR: %s 恶性细胞数为 0，无法生成 counts 矩阵。", label))
  }

  # 直接取 counts layer (稀疏矩阵, genes × cells)
  counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
  counts_mat <- counts_mat[, cells, drop = FALSE]

  cat(sprintf("  %s — 基因数: %d, 细胞数: %d\n", label, nrow(counts_mat), ncol(counts_mat)))
  cat(sprintf("    矩阵类型: %s, 总 UMI: %.2e\n",
              class(counts_mat)[1], sum(counts_mat)))
  counts_mat
}

count_tumor_PT <- extract_counts(seurat.obj, pt_barcodes, "count_tumor_PT")
count_tumor_LM <- extract_counts(seurat.obj, lm_barcodes, "count_tumor_LM")

# 保存
pt_file <- file.path(OUTPUT_DIR, "count_tumor_PT.rds")
lm_file <- file.path(OUTPUT_DIR, "count_tumor_LM.rds")

saveRDS(count_tumor_PT, pt_file)
saveRDS(count_tumor_LM, lm_file)

# 验证写入
stopifnot(file.exists(pt_file))
stopifnot(file.exists(lm_file))

cat(sprintf("\n  已保存:\n    %s\n    %s\n", pt_file, lm_file))

# ================================================================
# 完成
# ================================================================

cat("\n========== 02b 完成 — 可运行 Step 3  ==========\n")
cat(sprintf("  输入 Seurat 对象:  %s\n", seurat_file))
cat(sprintf("  输出 counts 矩阵:\n"))
cat(sprintf("    PT: %s (%d cells × %d genes)\n",
            pt_file, ncol(count_tumor_PT), nrow(count_tumor_PT)))
cat(sprintf("    LM: %s (%d cells × %d genes)\n",
            lm_file, ncol(count_tumor_LM), nrow(count_tumor_LM)))
