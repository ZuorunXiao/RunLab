# ============================================================================
# 06B_GSE178318_CellChat_KO_Comparison.R
# 功能: CellChat FBXL12 虚拟 KO vs WT 比较分析
# 定位: GSE178318 scRNA-seq 分析流程 Step 6B
# 依赖: CellChat v2.1.2, Seurat, Matrix, ggplot2, ComplexHeatmap, dplyr
# 官方 workflow（CellChat v2）:
#   subsetData -> identifyOverExpressedGenes(do.fast=FALSE)
#   -> identifyOverExpressedInteractions
#   -> computeCommunProb(type="triMean")
#   -> computeCommunProbPathway -> aggregateNet
# ============================================================================

# --- 包加载 ------------------------------------------------------
library(CellChat)       # v2.1.2
library(Seurat)         # Seurat v5
library(Matrix)         # 稀疏矩阵
library(ggplot2)        # 可视化
library(dplyr)          # 数据操作
library(ComplexHeatmap) # 热图
library(patchwork)      # 多图拼接

# --- 目录结构 ----------------------------------------------------
DATA_DIR <- "X:/CRLM/XX/results/02_Expression/table"
KO_DIR   <- "X:/CRLM/XX/results/03_GRN_VirtualKO/table"
FIG_DIR  <- "X:/CRLM/XX/results/06_CellChat_KO/Figure"
TAB_DIR  <- "X:/CRLM/XX/results/06_CellChat_KO/table"
dir.create(FIG_DIR, recursive = TRUE)
dir.create(TAB_DIR, recursive = TRUE)

options(future.globals.maxSize = Inf)

# --- 分析参数 ----------------------------------------------------
MAX_CELLS  <- 8000
FDR_CUTOFF <- 0.2
# KO 扰动因子：使用 scTenifoldKnk 连续 FC 值，截断到合理范围
# FC > 1: 基因在 KO 后上调; FC < 1: 基因在 KO 后下调
KO_FC_MIN  <- 0.2   # 最低倍数 (log2FC ≈ -2.32)
KO_FC_MAX  <- 5.0   # 最高倍数 (log2FC ≈ +2.32)

# --- SCI 标准绘图 ------------------------------------------------
sci_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.text     = element_text(size = 10, face = "plain"),
    axis.title      = element_text(size = 12, face = "plain"),
    axis.text       = element_text(size = 10, face = "plain"),
    plot.title      = element_text(size = 14, hjust = 0.5, face = "plain"),
    legend.title    = element_text(size = 11, face = "plain")
  )

sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  pdf(pdf_file, width = width, height = height)
  print(plot_obj)
  dev.off()
}

# --- 配色 -------------------------------------------------------
COL_WT <- "#2166AC"
COL_KO  <- "#B2182B"

# ============================================================================
# 辅助函数: CellChat v2 标准 workflow
# data.input: genes × cells 矩阵（无 rownames/col-names）
# meta.sub: data.frame，行名=细胞barcode
# ============================================================================
RunCellChatWorkflow <- function(data.input, meta.sub, label) {

  # Step 1: 对齐 barcode 与基因名
  # 确保 meta.sub 的 rownames 与 data.input 的列名一致
  bcs <- make.names(rownames(meta.sub), unique = TRUE)
  rownames(meta.sub)  <- bcs
  colnames(data.input) <- bcs

  # 设置基因行名（CellChat 需要 rownames 匹配 CellChatDB）
  rownames(data.input) <- gene_names

  # Step 2: 创建 CellChat 对象
  cc <- createCellChat(object = data.input, meta = meta.sub, group.by = "celltype")
  cat(sprintf("  [%s] 细胞类型: %s (%d 细胞)\n",
    label, paste(levels(cc@idents), collapse = ", "), ncol(data.input)))

  # Step 3: 设置 DB（Secreted Signaling）
  CellChatDB.use <- subsetDB(CellChatDB.human,
                             search = "Secreted Signaling", key = "annotation")
  cc@DB <- CellChatDB.use

  # Step 4: 预处理（CellChat v2 标准步骤）
  cc <- subsetData(cc)
  if (nrow(cc@data.signaling) == 0) {
    stop(sprintf("  [%s] data.signaling 为空 — rownames 未正确传递，请检查", label))
  }
  cc <- identifyOverExpressedGenes(cc, do.fast = FALSE)
  cc <- identifyOverExpressedInteractions(cc)
  cat(sprintf("  [%s] LR 对数: %d\n", label, nrow(cc@LR$LRsig)))

  # Step 5: 计算通讯概率
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cat(sprintf("  [%s] 信号通路数: %d\n", label, length(cc@netP$pathways)))

  # Step 6: 保存 CellChat 对象
  saveRDS(cc, file.path(TAB_DIR, paste0("CellChat_", label, ".rds")))

  # Step 7: 可视化 — 通讯网络图
  groupSize <- as.numeric(table(cc@idents))

  tryCatch({
    pdf(file.path(FIG_DIR, paste0("CirclePlot_Count_", label, ".pdf")), 10, 10)
    par(oma = c(0, 0, 2, 0))  # 外边距：bottom, left, top, right（top=2行）
    netVisual_circle(cc@net$count, vertex.weight = groupSize,
      weight.scale = TRUE, label.edge = FALSE,
      title.name = NULL)  # title.name 被 igraph margin 吃掉了，用 mtext
    mtext(paste0(label, " - Number of Interactions"),
      side = 3, line = 0.5, outer = TRUE, cex = 1.1, font = 2)
    dev.off()
  }, error = function(e) {
    cat(sprintf("    [警告] CirclePlot_Count_%s: %s\n", label, e$message))
  })

  tryCatch({
    pdf(file.path(FIG_DIR, paste0("CirclePlot_Weight_", label, ".pdf")), 10, 10)
    par(oma = c(0, 0, 2, 0))
    netVisual_circle(cc@net$weight, vertex.weight = groupSize,
      weight.scale = TRUE, label.edge = FALSE,
      title.name = NULL)
    mtext(paste0(label, " - Interaction Strength"),
      side = 3, line = 0.5, outer = TRUE, cex = 1.1, font = 2)
    dev.off()
  }, error = function(e) {
    cat(sprintf("    [警告] CirclePlot_Weight_%s: %s\n", label, e$message))
  })

  # Step 8: 通讯数据表
  df.comm <- tryCatch(subsetCommunication(cc), error = function(e) NULL)
  if (!is.null(df.comm) && nrow(df.comm) > 0) {
    write.csv(df.comm,
      file.path(TAB_DIR, paste0("Table_Communication_", label, ".csv")),
      row.names = FALSE)

    pathway.summary <- df.comm %>%
      group_by(pathway_name) %>%
      summarise(
        n_interactions = n(),
        mean_prob = mean(prob, na.rm = TRUE),
        max_prob  = max(prob, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(n_interactions))
    write.csv(pathway.summary,
      file.path(TAB_DIR, paste0("Table_Pathway_", label, ".csv")),
      row.names = FALSE)
  }

  gc()
  return(cc)
}

# ============================================================================
# 断点续传：若已有 RDS，跳过 PART1/2/3 的 CellChat 计算，直接进入可视化
# ============================================================================
rds_files <- c(
  counts_wt = file.path(TAB_DIR, "counts_WT.rds"),
  counts_ko = file.path(TAB_DIR, "counts_KO.rds"),
  meta      = file.path(TAB_DIR, "meta.rds"),
  cc_wt     = file.path(TAB_DIR, "CellChat_WT.rds"),
  cc_ko     = file.path(TAB_DIR, "CellChat_KO.rds")
)
all_rds_exist <- all(file.exists(rds_files))

if (all_rds_exist) {
  cat("检测到已有 RDS，跳过 PART1/2/3 CellChat 计算，直接进入 PART4 可视化...\n")
  counts_wt <- readRDS(rds_files["counts_wt"])
  counts_ko <- readRDS(rds_files["counts_ko"])
  meta.sub  <- readRDS(rds_files["meta"])
  cc_wt     <- readRDS(rds_files["cc_wt"])
  cc_ko     <- readRDS(rds_files["cc_ko"])

  # 重新派生 PART4 需要的变量
  pathways_wt  <- cc_wt@netP$pathways
  pathways_ko  <- cc_ko@netP$pathways
  pathways_all <- union(pathways_wt, pathways_ko)
  common_pw    <- intersect(pathways_wt, pathways_ko)

  skip_to_part4 <- TRUE
} else {
  skip_to_part4 <- FALSE
}

# ============================================================================
# PART 1: 数据加载 + 分层下采样 + KO 矩阵构建
# ============================================================================

if (!skip_to_part4) {
  cat("========== PART 1: 数据加载 ==========\n")

seurat_wt    <- readRDS(file.path(DATA_DIR, "seurat_obj_FBXL12_status.rds"))
counts_full  <- seurat_wt@assays$RNA@layers$counts  # genes × cells, 无 rownames/col-names
gene_names   <- rownames(seurat_wt)    # 基因符号向量（用于 lookup）
cell_names   <- colnames(seurat_wt)    # 细胞 barcode 向量
meta_all     <- seurat_wt@meta.data
ncell_all    <- ncol(seurat_wt)
rm(seurat_wt); gc()

cat(sprintf("原始: %d 基因 × %d 细胞\n", nrow(counts_full), ncell_all))

# --- 分层下采样（按 celltype 比例保持） -----------------------
set.seed(42)
keep_list <- lapply(split(seq_len(ncell_all), meta_all$celltype), function(idx) {
  sample(idx, min(length(idx), round(MAX_CELLS * length(idx) / ncell_all)))
})
keep_idx <- unlist(keep_list)
meta.sub <- meta_all[keep_idx, , drop = FALSE]
rm(meta_all); gc()

counts_sub <- counts_full[, keep_idx, drop = FALSE]
rm(counts_full); gc()

cat(sprintf("下采样: %d 细胞\n", nrow(meta.sub)))
cat(sprintf("细胞类型: %s\n", paste(table(meta.sub$celltype), collapse = ", ")))

# --- 识别 FBXL12 KO 调控基因 -----------------------------------
# 优先使用 03 产出的无 FBXL12 版本(Table1/2_KO_*_diffRegulation_nofbxl12.csv);
# 若 03 尚未重跑(无 _nofbxl12), 从 RAW 表剔除被扰动基因 FBXL12 作为兜底
read_ko_nofbxl12 <- function(fname) {
  p_nof <- file.path(KO_DIR, sub("\\.csv$", "_nofbxl12.csv", fname))
  p_raw <- file.path(KO_DIR, fname)
  if (file.exists(p_nof)) {
    cat(sprintf("  载入无 FBXL12 版本: %s\n", basename(p_nof)))
    return(read.csv(p_nof))
  }
  if (file.exists(p_raw)) {
    d <- read.csv(p_raw)
    cat(sprintf("  警告: 无 FBXL12 版本缺失, 已从 RAW 剔除 FBXL12: %s\n", basename(p_raw)))
    return(d[d$gene != "FBXL12", ])
  }
  stop(sprintf("缺少 KO 差异调控表: %s (亦无 %s)", basename(p_raw), basename(p_nof)))
}
diff_pt <- read_ko_nofbxl12("Table1_KO_PT_diffRegulation.csv")
diff_lm <- read_ko_nofbxl12("Table2_KO_LM_diffRegulation.csv")

# 恶性细胞（PT/LM）的细胞索引
is_mal  <- meta.sub$malignant == "malignant" |
            grepl("Epithelial|Tumor", meta.sub$celltype, ignore.case = TRUE)
cells_pt <- which(is_mal & meta.sub$tissue == "PT")
cells_lm <- which(is_mal & meta.sub$tissue == "LM")

cat(sprintf("PT 恶性: %d, LM 恶性: %d\n", length(cells_pt), length(cells_lm)))

# E1：scTenifoldKnk 的 diffRegulation 中 Z=(distance-mean)/sd，非方向性指标；
# 真正方向由 FC 决定：FC>1 上调，FC<1 下调。原代码误用 Z>0/Z<0 判方向，导致
# PT/LM 下调基因数恒为 0（实测 Z 分布：PT 162 上调, 0 下调；LM 27 上调, 0 下调）。
n_pt_up   <- sum(diff_pt$p.adj < FDR_CUTOFF & diff_pt$FC > 1)
n_pt_down <- sum(diff_pt$p.adj < FDR_CUTOFF & diff_pt$FC < 1)
n_lm_up   <- sum(diff_lm$p.adj < FDR_CUTOFF & diff_lm$FC > 1)
n_lm_down <- sum(diff_lm$p.adj < FDR_CUTOFF & diff_lm$FC < 1)
cat(sprintf("PT 调控: 上调 %d, 下调 %d\n", n_pt_up, n_pt_down))
cat(sprintf("LM 调控: 上调 %d, 下调 %d\n", n_lm_up, n_lm_down))

# --- 构建 WT 矩阵 ------------------------------------------------
counts_wt <- counts_sub

# --- 构建 KO 矩阵（向量化，一次处理所有非零元素） -----------
# counts_sub: genes × cells, dgCMatrix
# @i : 0-based gene indices for each non-zero element
# @p : column pointers（@p[c] 到 @p[c+1]-1 是第 c 列的非零元素在 @i/@x 中的索引）
# @x : non-zero values
#
# 构建基因名 → 倍数 lookup（named vector，names=基因符号）
BuildKOFactor <- function(diff_reg, label = "") {
  factor_vec <- rep(1.0, length(gene_names))
  names(factor_vec) <- gene_names
  m <- diff_reg[diff_reg$gene %in% gene_names, ]
  n_sig     <- 0L
  n_na_fc   <- 0L
  n_clipped <- 0L

  for (j in seq_len(nrow(m))) {
    g <- m$gene[j]; fc <- m$FC[j]; p <- m$p.adj[j]
    if (!is.na(p) && p < FDR_CUTOFF) {
      n_sig <- n_sig + 1L
      # 跳过 FC 无效的基因（NA / 零 / 负值），保留默认 1.0
      if (is.na(fc) || fc <= 0) {
        n_na_fc <- n_na_fc + 1L
        next
      }
      # 截断极端 FC 值到合理范围
      if (fc < KO_FC_MIN) {
        n_clipped <- n_clipped + 1L
        fc <- KO_FC_MIN
      } else if (fc > KO_FC_MAX) {
        n_clipped <- n_clipped + 1L
        fc <- KO_FC_MAX
      }
      factor_vec[g] <- fc
    }
  }

  cat(sprintf("  [%s] 显著基因: %d, FC缺失/无效: %d, 截断: %d\n",
      label, n_sig, n_na_fc, n_clipped))
  invisible(factor_vec)
}

cat("构建 KO 矩阵（向量化）...\n")
ko_factor_pt <- BuildKOFactor(diff_pt, "PT")
ko_factor_lm <- BuildKOFactor(diff_lm, "LM")

# Step 1: 构建列 → 组别 lookup（0-based 列索引 → "PT_KO"/"LM_KO"/"WT"）
# 长度为 ncol（细胞数），names() = 0-based 列索引
col_group <- rep("WT", ncol(counts_sub))
names(col_group) <- as.character(seq_len(ncol(counts_sub)) - 1)  # 0-based
col_group[as.character(cells_pt - 1)] <- "PT_KO"
col_group[as.character(cells_lm - 1)] <- "LM_KO"

# Step 2: 对每个非零元素，用 @i（基因索引）+1 查对应的 KO factor
# @i = 0-based gene index；@p[c] 到 @p[c+1]-1 是第 c 列（0-based）的非零元素索引
all_cols_0    <- rep(seq_len(ncol(counts_sub)) - 1, diff(counts_sub@p))
all_genes_0    <- counts_sub@i  # 0-based gene index for each non-zero element
all_gene_names <- gene_names[all_genes_0 + 1]  # 基因符号
all_col_labels <- col_group[as.character(all_cols_0)]

# 一次性查所有 KO factor
ko_factors_all <- rep(1.0, length(all_cols_0))
ko_factors_all[all_col_labels == "PT_KO"] <- ko_factor_pt[all_gene_names[all_col_labels == "PT_KO"]]
ko_factors_all[all_col_labels == "LM_KO"] <- ko_factor_lm[all_gene_names[all_col_labels == "LM_KO"]]

# Step 3: 构建 KO 矩阵
counts_ko <- counts_sub
counts_ko@x <- counts_sub@x * ko_factors_all

# --- 语义补全: FBXL12 在 KO 细胞中置 0 -----------------------------------
# scTenifoldKnk 虚拟 KO 已给出下游靶标倍数, 但被扰动基因 FBXL12 自身
# 在表达矩阵中仍保留原值。真正的"敲除"应将其在 KO 细胞(PT/LM 恶性)中置 0。
fbx_idx <- which(gene_names == "FBXL12")
if (length(fbx_idx) == 1L) {
  ko_cols <- c(cells_pt, cells_lm)            # 1-based 列索引
  counts_ko[fbx_idx, ko_cols] <- 0
  cat(sprintf("  FBXL12 已在 %d 个 KO 细胞中置 0\n", length(ko_cols)))
} else {
  cat("  警告: 未在基因列表中找到 FBXL12, 跳过置 0\n")
}

cat(sprintf("  PT KO: %d 列, LM KO: %d 列\n",
  sum(all_col_labels == "PT_KO"), sum(all_col_labels == "LM_KO")))

cat(sprintf("KO 矩阵完成: %d × %d\n", nrow(counts_ko), ncol(counts_ko)))

# --- 保存矩阵和元数据 -------------------------------------------
# 设置基因行名（CellChat 需要基因名匹配 DB）
rownames(counts_wt) <- gene_names
rownames(counts_ko) <- gene_names
gc()
saveRDS(counts_wt, file.path(TAB_DIR, "counts_WT.rds"))
saveRDS(counts_ko, file.path(TAB_DIR, "counts_KO.rds"))
saveRDS(meta.sub,  file.path(TAB_DIR, "meta.rds"))

ko_info <- data.frame(
  param = c("MAX_CELLS", "KO_FC_MIN", "KO_FC_MAX", "FDR_CUTOFF",
             "PT_cells", "LM_cells",
             "PT_up", "PT_down", "LM_up", "LM_down",
             "FBXL12_zeroed", "nrow", "ncol"),
  value = c(MAX_CELLS, KO_FC_MIN, KO_FC_MAX, FDR_CUTOFF,
             length(cells_pt), length(cells_lm),
             n_pt_up, n_pt_down, n_lm_up, n_lm_down,
             ifelse(exists("fbx_idx") && length(fbx_idx) == 1L, "yes", "no"),
             nrow(counts_ko), ncol(counts_ko))
)
write.csv(ko_info, file.path(TAB_DIR, "ko_params.csv"), row.names = FALSE)
cat("矩阵 + 元数据保存完成\n")

rm(counts_sub); gc()

# ============================================================================
# PART 2/3: WT / KO CellChat
# ============================================================================
if (!skip_to_part4) {
  cat("\n========== PART 2: WT CellChat ==========\n")
  cc_wt <- RunCellChatWorkflow(counts_wt, meta.sub, "WT")

  cat("\n========== PART 3: KO CellChat ==========\n")
  cc_ko <- RunCellChatWorkflow(counts_ko, meta.sub, "KO")
}

# ============================================================================
# PART 4: WT vs KO 比较分析
# ============================================================================

pathways_wt  <- cc_wt@netP$pathways
pathways_ko  <- cc_ko@netP$pathways
pathways_all <- union(pathways_wt, pathways_ko)
cat(sprintf("WT 通路: %d, KO 通路: %d, 合计: %d\n",
  length(pathways_wt), length(pathways_ko), length(pathways_all)))

# --- 4.1 通路权重汇总 ----------------------------------------
pathway_stats <- lapply(pathways_all, function(pw) {
  # M2 修复：cc@netP$prob 是正常 3D 数组，索引不会抛错；移除无意义 tryCatch
  w_wt <- sum(cc_wt@netP$prob[, , pw], na.rm = TRUE)
  w_ko <- sum(cc_ko@netP$prob[, , pw], na.rm = TRUE)
  tibble(
    pathway    = pw,
    weight_WT  = w_wt,
    weight_KO  = w_ko,
    diff       = w_ko - w_wt,
    log2FC     = log2((w_ko + 1e-6) / (w_wt + 1e-6)),
    status     = if (pw %in% pathways_wt && pw %in% pathways_ko) "Common" else
                 if (pw %in% pathways_wt) "WT_only" else "KO_only"
  )
}) %>% bind_rows() %>% arrange(desc(abs(diff)))

write.csv(pathway_stats,
  file.path(TAB_DIR, "Table_Pathway_Comparison.csv"), row.names = FALSE)

# --- 4.2 柱状图（只显示非零差异的通路） ---------------------
# 大量通路 diff=0，画出来会造成"空图"错觉，故过滤后展示
nonzero_pw <- pathway_stats %>% filter(abs(diff) > 0)
top_n      <- min(15, nrow(nonzero_pw))
top15      <- nonzero_pw %>% slice_max(abs(diff), n = top_n)

if (nrow(top15) > 0) {
  p_bar <- ggplot(top15,
    aes(x = reorder(pathway, diff), y = diff)) +
    geom_bar(stat = "identity",
      aes(fill = ifelse(diff > 0, COL_KO, COL_WT)), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_text(aes(label = formatC(diff, format = "e", digits = 1)),
              hjust = -0.15, size = 2.8, color = "grey20") +
    coord_flip() + sci_theme +
    labs(x = "Pathway", y = "KO - WT (Weight Difference)",
         title = "FBXL12 KO vs WT: Pathway Communication Difference",
         subtitle = sprintf("Only %d/%d pathways show non-zero change; differences are very small",
                            nrow(nonzero_pw), nrow(pathway_stats))) +
    scale_fill_identity() +
    theme(axis.text.y = element_text(size = 9)) +
    expand_limits(y = max(abs(top15$diff)) * 1.3)
  sci_save(p_bar, "Pathway_Comparison_WT_vs_KO",
    width = 8, height = max(4, 0.6 + nrow(top15) * 0.35))
}

# --- 4.3 热图（只显示非零 log2FC 的共同通路） ----------------
# 注意：与 pathway_stats 保持一致，使用 netP$prob 计算 log2FC；
# netP$weight 是聚合后的通讯强度，可能因数值修约导致 WT/KO 完全相同。
common_pw <- intersect(pathways_wt, pathways_ko)
if (length(common_pw) > 0) {
  prob_wt <- cc_wt@netP$prob
  prob_ko <- cc_ko@netP$prob

  diff_vec <- sapply(common_pw, function(pw) {
    w_wt <- sum(prob_wt[, , pw], na.rm = TRUE)
    w_ko <- sum(prob_ko[, , pw], na.rm = TRUE)
    log2((w_ko + 1e-10) / (w_wt + 1e-10))
  })
  names(diff_vec) <- common_pw

  # 过滤掉 log2FC 为 0 的通路，避免空白挤压
  diff_vec_nz <- diff_vec[abs(diff_vec) > 0]

  if (length(diff_vec_nz) > 0) {
    diff_vec_nz <- sort(diff_vec_nz, decreasing = TRUE)

    p_heat <- ggplot(tibble(pathway = names(diff_vec_nz), log2FC = diff_vec_nz),
      aes(x = pathway, y = 1, fill = log2FC)) +
      geom_tile() +
      geom_text(aes(label = formatC(log2FC, format = "f", digits = 3)),
                color = "black", size = 3) +
      scale_fill_gradient2(low = COL_WT, high = COL_KO, mid = "white",
        name = "log2FC", midpoint = 0) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
            axis.text.y = element_blank(),
            panel.grid  = element_blank(),
            axis.ticks  = element_blank()) +
      labs(x = "Pathway", y = "", title = "Pathway Communication Change (log2FC)",
           subtitle = sprintf("%d/%d common pathways show non-zero change",
                              length(diff_vec_nz), length(diff_vec)))
    sci_save(p_heat, "Pathway_Diff_Heatmap",
      width = max(6, 2 + length(diff_vec_nz) * 0.8), height = 4.5)
  } else {
    cat("[注意] 所有共同通路的 log2FC 均为 0，跳过热图绘制\n")
  }
}

# --- 4.4 LR 对通讯差异 ---------------------------------------
df_comm_wt <- tryCatch(subsetCommunication(cc_wt), error = function(e) NULL)
df_comm_ko <- tryCatch(subsetCommunication(cc_ko), error = function(e) NULL)

if (!is.null(df_comm_wt) && !is.null(df_comm_ko) &&
    nrow(df_comm_wt) > 0 && nrow(df_comm_ko) > 0) {

  # 同一通路+LR 对可能来自多个细胞类型对，先按 LR 对聚合（取总概率），
  # 再比较 WT vs KO，避免 many-to-many join
  df_wt_agg <- df_comm_wt %>%
    group_by(pathway_name, ligand, receptor) %>%
    summarise(prob_WT = sum(prob, na.rm = TRUE), .groups = "drop")

  df_ko_agg <- df_comm_ko %>%
    group_by(pathway_name, ligand, receptor) %>%
    summarise(prob_KO = sum(prob, na.rm = TRUE), .groups = "drop")

  lr_compare <- inner_join(df_wt_agg, df_ko_agg,
    by = c("pathway_name", "ligand", "receptor")) %>%
    mutate(
      pathway = pathway_name,
      diff    = prob_KO - prob_WT,
      log2FC  = log2((prob_KO + 1e-10) / (prob_WT + 1e-10))
    ) %>%
    select(pathway, ligand, receptor, prob_WT, prob_KO, diff, log2FC) %>%
    # 过滤弱信号对：至少有一方概率 > 1e-8
    filter(pmin(prob_WT, prob_KO) > 1e-8) %>%
    arrange(desc(abs(log2FC)))  # 按 |log2FC| 排序

  write.csv(lr_compare,
    file.path(TAB_DIR, "Table_LR_Comparison.csv"), row.names = FALSE)

  # 只保留非零差异的 LR 对用于可视化
  lr_nz <- lr_compare %>% filter(abs(log2FC) > 0)

  if (nrow(lr_nz) > 0) {
    top_n_lr <- min(20, nrow(lr_nz))
    top_lr <- lr_nz %>% slice_head(n = top_n_lr) %>%
      mutate(label = ifelse(nchar(paste0(ligand, "-", receptor)) > 24,
        paste0(substr(paste0(ligand, "-", receptor), 1, 22), ".."),
        paste0(ligand, "-", receptor)),
        log2FC_label = formatC(log2FC, format = "f", digits = 3))

    # 主图：展示所有非零 LR 对
    p_lr <- ggplot(top_lr,
      aes(x = log2FC, y = reorder(label, log2FC),
          size = abs(log2FC), fill = log2FC)) +
      geom_point(shape = 21, stroke = 0.5, color = "black", alpha = 0.85) +
      geom_text(aes(label = log2FC_label), hjust = -0.35, size = 3) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      scale_fill_gradient2(low = COL_WT, high = COL_KO, name = "log2FC") +
      scale_size_continuous(name = "|log2FC|", range = c(4, 10)) +
      sci_theme +
      labs(x = "log2FC (KO / WT)", y = "Ligand - Receptor",
           title = "LR Pairs with Non-zero |log2FC| (prob > 1e-8)",
           subtitle = sprintf("%d/%d LR pairs show non-zero change; shown top %d",
                              nrow(lr_nz), nrow(lr_compare), top_n_lr)) +
      theme(axis.text.y = element_text(size = 8),
            legend.position = "right") +
      expand_limits(x = max(abs(top_lr$log2FC)) * 1.15)
    sci_save(p_lr, "LR_Pairs_Diff_Bubble",
      width = 8.5, height = max(4, 0.6 + nrow(top_lr) * 0.4))

    # 附图：放大 log2FC 较小（< 0.1）的条目，避免 MIF 点被压缩到看不见
    small_lr <- top_lr %>% filter(abs(log2FC) < 0.1)
    if (nrow(small_lr) > 0) {
      p_lr_zoom <- ggplot(small_lr,
        aes(x = log2FC, y = reorder(label, log2FC),
            size = abs(log2FC), fill = log2FC)) +
        geom_point(shape = 21, stroke = 0.5, color = "black", alpha = 0.85) +
        geom_text(aes(label = log2FC_label), hjust = -0.2, size = 3) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
        scale_fill_gradient2(low = COL_WT, high = COL_KO, name = "log2FC") +
        scale_size_continuous(name = "|log2FC|", range = c(5, 10)) +
        sci_theme +
        labs(x = "log2FC (KO / WT)", y = "Ligand - Receptor",
             title = "Zoom: LR Pairs with small |log2FC|",
             subtitle = sprintf("Focused view for %d LR pair(s) with |log2FC| < 0.1",
                                nrow(small_lr))) +
        theme(axis.text.y = element_text(size = 8),
              legend.position = "right") +
        expand_limits(x = max(abs(small_lr$log2FC)) * 1.3)
      sci_save(p_lr_zoom, "LR_Pairs_Diff_Bubble_Zoom",
        width = 8.5, height = max(3, 0.6 + nrow(small_lr) * 0.5))
    }
  } else {
    cat("[注意] 所有 LR 对的 log2FC 均为 0，跳过气泡图绘制\n")
  }
}

# --- 4.5 共同通路列表 ----------------------------------------
if (length(common_pw) > 0) {
  write.csv(data.frame(pathway = sort(common_pw)),
    file.path(TAB_DIR, "Table_Common_Pathways.csv"), row.names = FALSE)
}

cat("\n========== 完成 ==========\n")
cat(sprintf("图件目录: %s\n", FIG_DIR))
cat(sprintf("表格目录: %s\n", TAB_DIR))
cat(sprintf("WT CellChat: %d 通路\n", length(pathways_wt)))
cat(sprintf("KO CellChat: %d 通路\n", length(pathways_ko)))
cat(sprintf("共同通路: %d\n", length(common_pw)))
