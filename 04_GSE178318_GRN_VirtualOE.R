# ================================================================
# 04_GSE178318_GRN_VirtualOE.R
# 功能: FBXL12虚拟过表达效应推断--基于2×2因子设计的GRN差异网络分析
#       仅在恶性细胞内进行,按组织来源(PT/LM)分层
# 定位: GSE178318 scRNA-seq 分析流程 Step 4
# 方法: scTenifoldNet (X = FBXL12_high, Y = FBXL12_low)
# 参数: Top=1000 (Mean排序), nc_nNet=10, nc_nCells=500, nCores=4
#       PT/LM/ALL 三组独立分析
#
# ---- 统计量语义(官方 dRegulation,务必注意) ----------------------
#   distance = 基因在两条件流形对齐坐标间的欧氏距离
#   FC       = distance^2 / mean(distance^2)
#              官方定义为 "fold change with respect to the expectation",
#              即相对全基因期望的倍数。**该量无方向性**:
#              FC 大 = 该基因调控关系改变幅度大,
#              既不代表上调/激活,也不代表下调/抑制。
#   p.value  = pchisq(FC, df = 1, lower.tail = FALSE)
#   p.adj    = BH 校正
#   推论: p 由 FC 单调决定 => FDR<0.05 等价于 FC > qchisq(0.95,1) = 3.84。
#         因此筛选只需 p.adj < 0.05; 附加 "FC > 1" 恒真且冗余(已实测确认)。
#         同理不可绘制火山图(x=log2FC 与 y=-log10FDR 单调共线),
#         官方标准可视化为 chi-square Q-Q plot。
# ================================================================

# --- 包加载 -------------------------------------------------------
library(Matrix)           # 稀疏矩阵处理
library(Seurat)           # 单细胞数据对象
library(scTenifoldNet)    # GRN构建与比较
library(org.Hs.eg.db)     # 基因注释(ENTREZID映射)
library(clusterProfiler)  # GO富集分析
library(ReactomePA)       # Reactome通路富集(本地数据库)
library(ggplot2)          # 数据可视化
library(ggrepel)          # 标签去重叠
library(ggvenn)           # Venn图(无系统依赖)
library(patchwork)        # 多图排版
library(reshape2)         # melt(热图长格式转换)
library(RhpcBLASctl)      # 限制BLAS线程,避免与nCores争抢

blas_set_num_threads(4)

# --- 辅助函数 ----------------------------------------------------

# Scientific plot theme(SCI出版标准:白色背景、无网格线、统一字号)
SciTheme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.text      = element_text(size = 10, face = "plain"),
    axis.title       = element_text(size = 12, face = "plain"),
    axis.text        = element_text(size = 10, face = "plain"),
    plot.title       = element_text(size = 14, hjust = 0.5, face = "plain"),
    legend.title     = element_text(size = 11, face = "plain")
  )

# Save plot as PDF(仅输出PDF,符合SCI投稿要求)
SciSave <- function(plot_obj, filename, width = 8, height = 6) {
  pdf(file.path(FIG_DIR, paste0(filename, ".pdf")), width = width, height = height)
  print(plot_obj)
  dev.off()
}

# ==============================================================
# GRN 稀疏化: Pearson 相关性 + |rho|>=0.7 硬阈值
# ==============================================================
# 注意: 该 GRN 为"共表达网络"(非因果调控网络),用于 Step 05 中心性分析时
#       hub gene 代表"与 FBXL12 共同响应基因"。
N_SAMPLE_HIGH <- 800
N_SAMPLE_LOW  <- 1200

BuildCorrGRN <- function(counts, group_tag, rho_thresh = 0.7) {
  # counts: genes x cells dgCMatrix
  n_cells  <- ncol(counts)
  n_sample <- if (grepl("_high$", group_tag)) N_SAMPLE_HIGH else N_SAMPLE_LOW
  if (n_cells > n_sample) {
    set.seed(42)
    counts <- counts[, sample(colnames(counts), n_sample), drop = FALSE]
    cat(sprintf("  [%s] 采样 %d -> %d cells\n", group_tag, n_cells, n_sample))
  }
  mat <- log1p(as.matrix(t(counts)))   # cells x genes (logCPM)
  rm(counts); gc()
  rho <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
  rm(mat); gc()
  diag(rho) <- 0
  rho[is.na(rho)] <- 0          # sd=0 基因对 -> 无相关性
  rho[abs(rho) < rho_thresh] <- 0
  rho <- drop0(as(rho, "dgCMatrix"))
  n_edges        <- nnzero(rho) / 2
  total_possible <- nrow(rho) * (nrow(rho) - 1) / 2
  cat(sprintf("  [%s] Pearson GRN: %d edges (density=%.4f, |rho|>=%.1f)\n",
              group_tag, n_edges, n_edges / total_possible, rho_thresh))
  rho
}

# 构建并保存一组(high/low)的 GRN 邻接矩阵,已存在则跳过
SaveGRN <- function(tag, counts_high, counts_low) {
  f_high <- file.path(TAB_DIR, sprintf("GRN_%s_high.rds", tag))
  f_low  <- file.path(TAB_DIR, sprintf("GRN_%s_low.rds",  tag))
  if (!file.exists(f_high)) {
    saveRDS(BuildCorrGRN(counts_high, paste0(tag, "_high")), f_high)
  }
  if (!file.exists(f_low)) {
    saveRDS(BuildCorrGRN(counts_low, paste0(tag, "_low")), f_low)
  }
}

# 官方标准可视化: chi-square Q-Q plot (scTenifoldNet README)
# 横轴为 chi-square(df=1) 理论分位数,纵轴为观测 FC;偏离对角线者即差异调控基因
QQPlotChisq <- function(df, group_label, x_max, y_max) {
  d <- df[order(df$FC), ]
  d$theo <- qchisq(ppoints(nrow(d)), df = 1)
  d$sig  <- ifelse(d$p.adj < 0.05, "FDR < 0.05", "NS")
  lab    <- tail(d[d$sig == "FDR < 0.05", ], 10)

  ggplot(d, aes(x = theo, y = FC)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(aes(colour = sig), size = 1, alpha = 0.75) +
    scale_colour_manual(values = c("FDR < 0.05" = SIG_COL$sig, "NS" = SIG_COL$ns),
                        name = NULL) +
    ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 2.5,
                             fontface = "italic", max.overlaps = Inf,
                             min.segment.length = 0) +
    coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max)) +
    labs(title = group_label,
         x = expression(chi^2 ~ "quantiles (df = 1)"),
         y = expression(FC[GRN])) +
    SciTheme +
    theme(plot.title = element_text(size = 11))
}

# --- 全局配色方案 -------------------------------------------------
# FC 无方向性,故仅区分"显著/不显著",不设 Up/Down
SIG_COL <- list(
  sig = "#D62728",   # FDR < 0.05
  ns  = "grey70"     # 不显著
)

# --- 目录结构 ----------------------------------------------------

DATA_DIR <- "X:/CRLM/XX/results/02_Expression/table"
FIG_DIR  <- "X:/CRLM/XX/results/04_VirtualOE/Figure"
TAB_DIR  <- "X:/CRLM/XX/results/04_VirtualOE/table"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# ================================================================
# PART 1: 数据加载 + 按tissue分层提取counts
# 目的: 提取恶性细胞,按tissue(PT/LM)×FBXL12_status分层
# 数据: seurat_obj_FBXL12_status.rds (来自 Step 02)
# ================================================================

TMR_SCRIPT <- Sys.time()

seurat_obj <- readRDS(file.path(DATA_DIR, "seurat_obj_FBXL12_status.rds"))
seurat_obj <- JoinLayers(seurat_obj)

# 提取恶性细胞(malignant == "malignant" OR celltype含Epithelial/Tumor)
is_malignant <- seurat_obj$malignant == "malignant"
is_epi_tumor <- grepl("Epithelial|Tumor", seurat_obj$celltype, ignore.case = TRUE)
seurat_tumor <- subset(seurat_obj, cells = which(is_malignant | is_epi_tumor))

cat("  总恶性细胞数:", ncol(seurat_tumor), "\n")
cat("  tissue分布:", table(seurat_tumor$tissue), "\n")
cat("  FBXL12_status分布:", table(seurat_tumor$FBXL12_status), "\n\n")

# 按 tissue × FBXL12_status 提取四组细胞
GetGroupCounts <- function(obj, tissue_tag, status_tag) {
  cells <- colnames(obj)[obj$tissue == tissue_tag & obj$FBXL12_status == status_tag]
  GetAssayData(obj, assay = "RNA", layer = "counts")[, cells]
}

counts_PT_high <- GetGroupCounts(seurat_tumor, "PT", "FBXL12_high")
counts_PT_low  <- GetGroupCounts(seurat_tumor, "PT", "FBXL12_low")
counts_LM_high <- GetGroupCounts(seurat_tumor, "LM", "FBXL12_high")
counts_LM_low  <- GetGroupCounts(seurat_tumor, "LM", "FBXL12_low")

cat(sprintf("  PT  FBXL12_high: %d cells\n", ncol(counts_PT_high)))
cat(sprintf("  PT  FBXL12_low:  %d cells\n", ncol(counts_PT_low)))
cat(sprintf("  LM  FBXL12_high: %d cells\n", ncol(counts_LM_high)))
cat(sprintf("  LM  FBXL12_low:  %d cells\n", ncol(counts_LM_low)))

# 组5: ALL(PT+LM合并)
counts_ALL_high <- cbind(counts_PT_high, counts_LM_high)
counts_ALL_low  <- cbind(counts_PT_low,  counts_LM_low)
cat(sprintf("  ALL FBXL12_high: %d cells (PT %d + LM %d)\n",
            ncol(counts_ALL_high), ncol(counts_PT_high), ncol(counts_LM_high)))
cat(sprintf("  ALL FBXL12_low:  %d cells (PT %d + LM %d)\n",
            ncol(counts_ALL_low), ncol(counts_PT_low), ncol(counts_LM_low)))

# 保存五组counts矩阵
saveRDS(counts_PT_high,  file.path(TAB_DIR, "counts_PT_high.rds"))
saveRDS(counts_PT_low,   file.path(TAB_DIR, "counts_PT_low.rds"))
saveRDS(counts_LM_high,  file.path(TAB_DIR, "counts_LM_high.rds"))
saveRDS(counts_LM_low,   file.path(TAB_DIR, "counts_LM_low.rds"))
saveRDS(counts_ALL_high, file.path(TAB_DIR, "counts_ALL_high.rds"))
saveRDS(counts_ALL_low,  file.path(TAB_DIR, "counts_ALL_low.rds"))

# ================================================================
# PART 2: 基因过滤--combined检测率 + Mean排序
# 目的: 筛选高表达代表性基因,确保三组scTenifoldNet比较时
#       high组和low组均有足够的基因表达量
# 方法: combined >=5%检测率 → Top 1000 by Mean
# 依据: scTenifoldNet内部做CPM标准化,基因过滤用Mean而非CV;
#       Mean排序选择"高表达代表性基因",适合GRN拓扑推断
# ================================================================


# 用ALL恶性细胞计算检测率(覆盖所有恶性细胞)
counts_all  <- cbind(counts_PT_high, counts_PT_low, counts_LM_high, counts_LM_low)
n_all_cells <- ncol(counts_all)
cat("  ALL counts:", nrow(counts_all), "genes x", n_all_cells, "cells\n")

# 步骤2A: combined检测率过滤(>=5%细胞表达)
detected_pct <- Matrix::rowSums(counts_all > 0) / n_all_cells
passed_det   <- detected_pct >= 0.05
n_passed_det <- sum(passed_det)
cat(sprintf("  通过检测率过滤 (>=5%% cells): %d / %d genes\n",
            n_passed_det, nrow(counts_all)))

# 步骤2B: 按Mean排序选择Top基因
counts_qc  <- counts_all[passed_det, ]
gene_mean  <- Matrix::rowMeans(counts_qc)  # 不做log2变换,与scTenifoldNet内部CPM一致
N_TOP      <- 1000
mean_order <- order(gene_mean, decreasing = TRUE)
top_genes  <- rownames(counts_qc)[mean_order[seq_len(min(N_TOP, nrow(counts_qc)))]]

fbxl12_mean_rank <- sum(gene_mean > gene_mean["FBXL12"], na.rm = TRUE) + 1
if (!"FBXL12" %in% top_genes) {
  cat(sprintf("  FBXL12不在Top%d中(Mean排名%d/%d),强制纳入\n",
              N_TOP, fbxl12_mean_rank, n_passed_det))
  top_genes <- c("FBXL12", top_genes[seq_len(N_TOP - 1)])
} else {
  cat(sprintf("  FBXL12在Top%d中(Mean排名%d/%d)\n",
              N_TOP, fbxl12_mean_rank, n_passed_det))
}
cat(sprintf("  最终基因集: %d genes (Mean排序)\n", length(top_genes)))

# 步骤2C: 直接按基因子集切片
# 依据: libsize过滤在Seurat QC阶段(Step 01)已完成,此处不额外过滤;
#       scTenifoldNet内部CPM标准化对libsize差异鲁棒
counts_PT_high_filt  <- counts_PT_high[top_genes, , drop = FALSE]
counts_PT_low_filt   <- counts_PT_low[top_genes,  , drop = FALSE]
counts_LM_high_filt  <- counts_LM_high[top_genes, , drop = FALSE]
counts_LM_low_filt   <- counts_LM_low[top_genes,  , drop = FALSE]
counts_ALL_high_filt <- cbind(counts_PT_high_filt, counts_LM_high_filt)
counts_ALL_low_filt  <- cbind(counts_PT_low_filt,  counts_LM_low_filt)

cat(sprintf("  PT组:  high=%d cells, low=%d cells (x %d genes)\n",
            ncol(counts_PT_high_filt), ncol(counts_PT_low_filt), length(top_genes)))
cat(sprintf("  LM组:  high=%d cells, low=%d cells (x %d genes)\n",
            ncol(counts_LM_high_filt), ncol(counts_LM_low_filt), length(top_genes)))
cat(sprintf("  ALL组: high=%d cells (PT %d + LM %d), low=%d cells (PT %d + LM %d)\n",
            ncol(counts_ALL_high_filt), ncol(counts_PT_high_filt), ncol(counts_LM_high_filt),
            ncol(counts_ALL_low_filt),  ncol(counts_PT_low_filt),  ncol(counts_LM_low_filt)))

# 保存过滤后的counts矩阵
saveRDS(counts_PT_high_filt,  file.path(TAB_DIR, "counts_PT_high_filtered.rds"))
saveRDS(counts_PT_low_filt,   file.path(TAB_DIR, "counts_PT_low_filtered.rds"))
saveRDS(counts_LM_high_filt,  file.path(TAB_DIR, "counts_LM_high_filtered.rds"))
saveRDS(counts_LM_low_filt,   file.path(TAB_DIR, "counts_LM_low_filtered.rds"))
saveRDS(counts_ALL_high_filt, file.path(TAB_DIR, "counts_ALL_high_filtered.rds"))
saveRDS(counts_ALL_low_filt,  file.path(TAB_DIR, "counts_ALL_low_filtered.rds"))

# 保存基因列表
write.csv(
  data.frame(gene             = top_genes,
             mean_expr        = gene_mean[top_genes],
             det_rate_overall = detected_pct[top_genes],
             stringsAsFactors = FALSE),
  file.path(TAB_DIR, "gene_list_top1000_byMean.csv"),
  row.names = FALSE
)

# ================================================================
# PART 3: scTenifoldNet--PT组(主分析)
# 目的: PT组织中FBXL12_high vs FBXL12_low的GRN差异
# 方法: scTenifoldNet(X = PT_high, Y = PT_low)
# 参数: nc_nCells=500(官方默认子采样值,与Step 03统一;
#       高亚组不足500时封顶为实际细胞数), nc_nNet=10, nCores=4
# ================================================================


nCells_PT <- min(500L, ncol(counts_PT_high_filt))
tab_PT    <- file.path(TAB_DIR, "Table1_STN_PT_diffRegulation.csv")

if (file.exists(tab_PT)) {
  cat("  PT scTenifoldNet 已有结果,跳过\n")
  diff_PT <- read.csv(tab_PT, stringsAsFactors = FALSE)
} else {
  cat(sprintf("  PT组参数: nc_nCells=%d, nc_nNet=10, Top_genes=%d\n",
              nCells_PT, length(top_genes)))
  cat(sprintf("  [%s] scTenifoldNet PT 开始\n", format(Sys.time(), "%H:%M:%S")))
  tmr_PT <- Sys.time()
  stn_PT <- scTenifoldNet(
    X         = counts_PT_high_filt,
    Y         = counts_PT_low_filt,
    qc        = FALSE,
    nc_nNet   = 10L,
    nc_nCells = nCells_PT,
    nc_nComp  = 3L,
    td_K      = 3L,
    ma_nDim   = 30L,
    nCores    = 4L
  )
  cat(sprintf("  [%s] scTenifoldNet PT 完成 (用时 %.1f min)\n",
              format(Sys.time(), "%H:%M:%S"),
              as.numeric(difftime(Sys.time(), tmr_PT, units = "mins"))))
  diff_PT <- stn_PT$diffRegulation
  rm(stn_PT); gc()
  write.csv(diff_PT, tab_PT, row.names = FALSE)
}

diff_PT$group <- "PT"
cat(sprintf("  PT 显著差异调控靶标 (FDR<0.05): %d / %d\n",
            sum(diff_PT$p.adj < 0.05, na.rm = TRUE), nrow(diff_PT)))

# GRN 邻接矩阵(供 Step 05 中心性分析)
SaveGRN("PT", counts_PT_high_filt, counts_PT_low_filt)

# ================================================================
# PART 4: scTenifoldNet--LM组(探索性)
# ================================================================


nCells_LM <- min(500L, ncol(counts_LM_high_filt))
tab_LM    <- file.path(TAB_DIR, "Table2_STN_LM_diffRegulation.csv")

if (file.exists(tab_LM)) {
  cat("  LM scTenifoldNet 已有结果,跳过\n")
  diff_LM <- read.csv(tab_LM, stringsAsFactors = FALSE)
} else {
  cat(sprintf("  LM组参数: nc_nCells=%d, nc_nNet=10, Top_genes=%d\n",
              nCells_LM, length(top_genes)))
  cat(sprintf("  [%s] scTenifoldNet LM 开始\n", format(Sys.time(), "%H:%M:%S")))
  tmr_LM <- Sys.time()
  stn_LM <- scTenifoldNet(
    X         = counts_LM_high_filt,
    Y         = counts_LM_low_filt,
    qc        = FALSE,
    nc_nNet   = 10L,
    nc_nCells = nCells_LM,
    nc_nComp  = 3L,
    td_K      = 3L,
    ma_nDim   = 30L,
    nCores    = 4L
  )
  cat(sprintf("  [%s] scTenifoldNet LM 完成 (用时 %.1f min)\n",
              format(Sys.time(), "%H:%M:%S"),
              as.numeric(difftime(Sys.time(), tmr_LM, units = "mins"))))
  diff_LM <- stn_LM$diffRegulation
  rm(stn_LM); gc()
  write.csv(diff_LM, tab_LM, row.names = FALSE)
}

diff_LM$group <- "LM"
cat(sprintf("  LM 显著差异调控靶标 (FDR<0.05): %d / %d\n",
            sum(diff_LM$p.adj < 0.05, na.rm = TRUE), nrow(diff_LM)))

SaveGRN("LM", counts_LM_high_filt, counts_LM_low_filt)

# ================================================================
# PART 5: scTenifoldNet--ALL组(合并稳健性验证)
# ================================================================


nCells_ALL <- min(500L, ncol(counts_ALL_high_filt))
tab_ALL    <- file.path(TAB_DIR, "Table3_STN_ALL_diffRegulation.csv")

if (file.exists(tab_ALL)) {
  cat("  ALL scTenifoldNet 已有结果,跳过\n")
  diff_ALL <- read.csv(tab_ALL, stringsAsFactors = FALSE)
} else {
  cat(sprintf("  ALL组参数: nc_nCells=%d, nc_nNet=10, Top_genes=%d\n",
              nCells_ALL, length(top_genes)))
  cat(sprintf("  [%s] scTenifoldNet ALL 开始\n", format(Sys.time(), "%H:%M:%S")))
  tmr_ALL <- Sys.time()
  stn_ALL <- scTenifoldNet(
    X         = counts_ALL_high_filt,
    Y         = counts_ALL_low_filt,
    qc        = FALSE,
    nc_nNet   = 10L,
    nc_nCells = nCells_ALL,
    nc_nComp  = 3L,
    td_K      = 3L,
    ma_nDim   = 30L,
    nCores    = 4L
  )
  cat(sprintf("  [%s] scTenifoldNet ALL 完成 (用时 %.1f min)\n",
              format(Sys.time(), "%H:%M:%S"),
              as.numeric(difftime(Sys.time(), tmr_ALL, units = "mins"))))
  diff_ALL <- stn_ALL$diffRegulation
  rm(stn_ALL); gc()
  write.csv(diff_ALL, tab_ALL, row.names = FALSE)
}

diff_ALL$group <- "ALL"
cat(sprintf("  ALL 显著差异调控靶标 (FDR<0.05): %d / %d\n",
            sum(diff_ALL$p.adj < 0.05, na.rm = TRUE), nrow(diff_ALL)))

SaveGRN("ALL", counts_ALL_high_filt, counts_ALL_low_filt)

# ================================================================
# PART 5b: 剔除被扰动基因 FBXL12 (在所有可视化之前)
# 目的: scTenifoldNet 的被扰动基因 FBXL12 其 distance 因被"移除"而人为最大,
#       是工具假象,须从下游靶标/DEG 列表中剔除(官方 Net FAQ 要求)。
# 方法: 从内存 RAW 对象剔除 FBXL12 -> 写出 _nofbxl12 可用表 ->
#       内存对象替换为剔除版,后续 PART6~10 全部自动使用剔除版。
#       原始 Table1/2/3_STN_*_diffRegulation.csv 保持不变(仍含 FBXL12)。
# 兼容断点续传: PART3/4/5 无论重算或读盘, diff_PT/LM/ALL 均含 FBXL12,
#       本段每次运行幂等(已剔除则无操作)。
# ================================================================


KO_GENE <- "FBXL12"

diff_PT  <- diff_PT[diff_PT$gene  != KO_GENE, , drop = FALSE]
diff_LM  <- diff_LM[diff_LM$gene  != KO_GENE, , drop = FALSE]
diff_ALL <- diff_ALL[diff_ALL$gene != KO_GENE, , drop = FALSE]

write.csv(diff_PT,  file.path(TAB_DIR, "Table1_STN_PT_diffRegulation_nofbxl12.csv"),  row.names = FALSE)
write.csv(diff_LM,  file.path(TAB_DIR, "Table2_STN_LM_diffRegulation_nofbxl12.csv"),  row.names = FALSE)
write.csv(diff_ALL, file.path(TAB_DIR, "Table3_STN_ALL_diffRegulation_nofbxl12.csv"), row.names = FALSE)

cat(sprintf("  PT/LM/ALL 剔除 %s 后基因数: %d / %d / %d\n",
            KO_GENE, nrow(diff_PT), nrow(diff_LM), nrow(diff_ALL)))

# ================================================================
# PART 6: 三组比较可视化--chi-square Q-Q plot
# 目的: 并排展示PT/LM/ALL三组差异调控基因相对 chi-square(df=1) 的偏离
# 方法: 官方README标准做法;三组统一坐标轴范围便于横向比较
# 注意: 不使用火山图。FC 与 p 严格单调共线,火山图会退化为单调曲线,
#       且 FC 下界可达 1e-25 (log2 ≈ -83) 会把坐标轴彻底拉爆。
# ================================================================


x_max <- qchisq(ppoints(nrow(diff_PT)), df = 1)[nrow(diff_PT)]
y_max <- max(c(diff_PT$FC, diff_LM$FC, diff_ALL$FC), na.rm = TRUE)

group_labels <- c(
  PT  = sprintf("PT (n_high=%d, n_low=%d)",  ncol(counts_PT_high_filt),  ncol(counts_PT_low_filt)),
  LM  = sprintf("LM (n_high=%d, n_low=%d)",  ncol(counts_LM_high_filt),  ncol(counts_LM_low_filt)),
  ALL = sprintf("ALL (n_high=%d, n_low=%d)", ncol(counts_ALL_high_filt), ncol(counts_ALL_low_filt))
)

p_qq_combined <-
  (QQPlotChisq(diff_PT,  group_labels[["PT"]],  x_max, y_max) |
   QQPlotChisq(diff_LM,  group_labels[["LM"]],  x_max, y_max) |
   QQPlotChisq(diff_ALL, group_labels[["ALL"]], x_max, y_max)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "FBXL12 Virtual Overexpression -- GRN Differential Regulation (PT vs LM vs ALL)",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 14))
  ) &
  theme(legend.position = "bottom")

SciSave(p_qq_combined, "QQplot_ThreeGroups", width = 14, height = 5.5)

# ================================================================
# PART 7: 靶标分类--Venn图
# 目的: 识别PT特异、LM特异、共同的差异调控靶标
# 方法: 三组 FDR<0.05 的交集/差集(FC>1 恒真,不作为筛选条件)
# ================================================================


targets_PT  <- diff_PT$gene[diff_PT$p.adj   < 0.05]
targets_LM  <- diff_LM$gene[diff_LM$p.adj   < 0.05]
targets_ALL <- diff_ALL$gene[diff_ALL$p.adj < 0.05]

cat(sprintf("  PT  显著差异调控靶标: %d\n", length(targets_PT)))
cat(sprintf("  LM  显著差异调控靶标: %d\n", length(targets_LM)))
cat(sprintf("  ALL 显著差异调控靶标: %d\n", length(targets_ALL)))

targets_common    <- intersect(targets_PT, targets_LM)   # PT∩LM
targets_PT_only   <- setdiff(targets_PT, targets_LM)     # PT\LM
targets_LM_only   <- setdiff(targets_LM, targets_PT)     # LM\PT
targets_PT_in_ALL <- intersect(targets_PT, targets_ALL)  # PT∩ALL(稳健性)
targets_LM_in_ALL <- intersect(targets_LM, targets_ALL)  # LM∩ALL(稳健性)

cat(sprintf("\n  共同靶标 (PT ∩ LM): %d\n", length(targets_common)))
cat(sprintf("  PT特异靶标 (PT \\ LM): %d\n", length(targets_PT_only)))
cat(sprintf("  LM特异靶标 (LM \\ PT): %d\n", length(targets_LM_only)))
cat(sprintf("  PT靶标在ALL中稳健 (PT ∩ ALL): %d / %d\n",
            length(targets_PT_in_ALL), length(targets_PT)))
cat(sprintf("  LM靶标在ALL中稳健 (LM ∩ ALL): %d / %d\n",
            length(targets_LM_in_ALL), length(targets_LM)))

# 分类结果表格
classification <- data.frame(
  gene     = c(targets_common, targets_PT_only, targets_LM_only),
  category = c(rep("Common (PT∩LM)", length(targets_common)),
               rep("PT-specific",    length(targets_PT_only)),
               rep("LM-specific",    length(targets_LM_only))),
  stringsAsFactors = FALSE
)
classification$FC_PT   <- diff_PT$FC[match(classification$gene, diff_PT$gene)]
classification$FC_LM   <- diff_LM$FC[match(classification$gene, diff_LM$gene)]
classification$FC_ALL  <- diff_ALL$FC[match(classification$gene, diff_ALL$gene)]
classification$FDR_PT  <- diff_PT$p.adj[match(classification$gene, diff_PT$gene)]
classification$FDR_LM  <- diff_LM$p.adj[match(classification$gene, diff_LM$gene)]
classification$FDR_ALL <- diff_ALL$p.adj[match(classification$gene, diff_ALL$gene)]

classification <- classification[order(classification$category, -classification$FC_PT), ]
write.csv(classification,
          file.path(TAB_DIR, "Table4_TargetClassification.csv"),
          row.names = FALSE)

# Venn图(ggvenn,无系统依赖)
p_venn <- ggvenn(list(PT = targets_PT, LM = targets_LM),
                 stroke_size   = 0.5,
                 text_size     = 3.5,
                 set_name_size = 4,
                 auto_scale    = FALSE) +
  labs(title = "FBXL12 Differentially Regulated Targets: PT vs LM") +
  theme(plot.title  = element_text(hjust = 0.5, size = 11, margin = margin(b = 10)),
        plot.margin = margin(10, 10, 10, 10))

SciSave(p_venn, "Venn_PT_LM_Targets", width = 7, height = 6)

# ================================================================
# PART 8: 共同靶标热图--三组FC对比
# 目的: 展示共同靶标在三组中的调控改变幅度一致性
# 方法: 共同靶标按ALL组FDR取Top30,热图展示PT/LM/ALL的log2(FC)
# ================================================================


if (length(targets_common) > 0) {
  n_show     <- min(30, length(targets_common))
  top_common <- targets_common[
    order(diff_ALL$p.adj[match(targets_common, diff_ALL$gene)])
  ][seq_len(n_show)]

  # 构建热图矩阵(行=基因,列=组)
  hm_mat <- matrix(NA_real_, nrow = n_show, ncol = 3,
                   dimnames = list(top_common, c("PT", "LM", "ALL")))
  hm_mat[, "PT"]  <- log2(diff_PT$FC[match(top_common, diff_PT$gene)])
  hm_mat[, "LM"]  <- log2(diff_LM$FC[match(top_common, diff_LM$gene)])
  hm_mat[, "ALL"] <- log2(diff_ALL$FC[match(top_common, diff_ALL$gene)])
  hm_mat[!is.finite(hm_mat)] <- NA   # 防 FC=0 时 log2 -> -Inf

  # 转为长格式;行序按 PT 组 log2FC 排序(直接取自矩阵行名,避免依赖melt顺序)
  gene_order <- rownames(hm_mat)[order(hm_mat[, "PT"])]
  hm_long <- reshape2::melt(hm_mat, varnames = c("gene", "group"),
                            value.name = "log2FC")
  hm_long$group <- factor(hm_long$group, levels = c("PT", "LM", "ALL"))
  hm_long$gene  <- factor(hm_long$gene,  levels = rev(gene_order))

  # 动态色标范围(以中位数为锚点,避免全红/全蓝)
  p_heatmap <- ggplot(hm_long, aes(x = group, y = gene, fill = log2FC)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low  = "#1F77B4", mid = "white", high = "#D62728",
                         midpoint  = median(hm_long$log2FC, na.rm = TRUE),
                         limits    = range(hm_long$log2FC, na.rm = TRUE),
                         na.value  = "grey90") +
    labs(title = "FBXL12 GRN Targets: log2(FC) Across PT, LM, and ALL",
         x = "", y = "", fill = expression(log[2](FC[GRN]))) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.y     = element_text(size = 7, face = "italic"),
      axis.text.x     = element_text(size = 10, face = "bold"),
      plot.title      = element_text(hjust = 0.5, size = 11),
      legend.position = "right"
    )

  SciSave(p_heatmap, "Heatmap_CommonTargets",
          width = 5, height = max(6, n_show * 0.2 + 1))
  cat(sprintf("  热图已生成: %d个共同靶标\n", n_show))
} else {
  cat("  无共同靶标,跳过热图\n")
}

# ================================================================
# PART 9: GO/Reactome富集分析--PT组(主分析)
# 目的: 注释PT组FBXL12差异调控靶基因的生物学功能
# 方法: clusterProfiler (GO BP/MF/CC) + ReactomePA
# ================================================================

cat(sprintf("  PT显著差异调控靶标: %d\n", length(targets_PT)))

if (length(targets_PT) >= 3) {
  entrez_PT <- clusterProfiler::bitr(targets_PT,
                                     fromType = "SYMBOL",
                                     toType   = "ENTREZID",
                                     OrgDb    = org.Hs.eg.db)
  cat(sprintf("  成功ENTREZID映射: %d/%d genes\n",
              nrow(entrez_PT), length(targets_PT)))

  if (nrow(entrez_PT) >= 3) {
    for (ont in c("BP", "MF", "CC")) {
      ego <- clusterProfiler::enrichGO(
        gene          = entrez_PT$ENTREZID,
        OrgDb         = org.Hs.eg.db,
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        minGSSize     = 5
      )
      ego_df <- as.data.frame(ego)
      if (nrow(ego_df) > 0) {
        p_go <- clusterProfiler::dotplot(ego, showCategory = 12) +
          labs(title = sprintf("PT FBXL12 Targets -- GO %s Enrichment", ont)) +
          theme(axis.text.y = element_text(size = 8)) +
          SciTheme
        SciSave(p_go, sprintf("GO_%s_PT", ont), width = 10, height = 7)
        write.csv(ego_df, file.path(TAB_DIR, sprintf("Table5_GO_%s_PT.csv", ont)),
                  row.names = FALSE)
      }
      cat(sprintf("  GO %s: %d terms\n", ont, nrow(ego_df)))
    }

    # Reactome通路(本地数据库;KEGG REST API 国内超时,故改用 ReactomePA)
    ereact    <- ReactomePA::enrichPathway(
      gene          = entrez_PT$ENTREZID,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      minGSSize     = 5,
      organism      = "human"
    )
    ereact_df <- as.data.frame(ereact)
    if (nrow(ereact_df) > 0) {
      p_react <- clusterProfiler::dotplot(ereact, showCategory = 12) +
        labs(title = "PT FBXL12 Targets -- Reactome Pathway Enrichment") +
        theme(axis.text.y = element_text(size = 8)) +
        SciTheme
      SciSave(p_react, "Reactome_PT", width = 10, height = 7)
      write.csv(ereact_df, file.path(TAB_DIR, "Table5_Reactome_PT.csv"),
                row.names = FALSE)
    }
    cat(sprintf("  Reactome: %d pathways\n", nrow(ereact_df)))
  } else {
    cat("  ENTREZID映射少于3个,跳过富集分析\n")
  }
} else {
  cat("  PT显著靶标少于3个,跳过富集分析\n")
}

# ================================================================
# PART 10: GO/Reactome富集分析--共同靶标
# 目的: 注释跨组织保守的FBXL12核心调控靶标
# ================================================================


if (length(targets_common) >= 3) {
  entrez_common <- clusterProfiler::bitr(targets_common,
                                         fromType = "SYMBOL",
                                         toType   = "ENTREZID",
                                         OrgDb    = org.Hs.eg.db)
  cat(sprintf("  共同靶标ENTREZID映射: %d/%d genes\n",
              nrow(entrez_common), length(targets_common)))

  if (nrow(entrez_common) >= 3) {
    for (ont in c("BP", "MF", "CC")) {
      ego_df <- as.data.frame(clusterProfiler::enrichGO(
        gene          = entrez_common$ENTREZID,
        OrgDb         = org.Hs.eg.db,
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        minGSSize     = 3
      ))
      if (nrow(ego_df) > 0) {
        write.csv(ego_df, file.path(TAB_DIR, sprintf("Table6_GO_%s_Common.csv", ont)),
                  row.names = FALSE)
      }
      cat(sprintf("  GO %s 共同靶标: %d terms\n", ont, nrow(ego_df)))
    }

    ereact_df <- as.data.frame(ReactomePA::enrichPathway(
      gene          = entrez_common$ENTREZID,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      minGSSize     = 3,
      organism      = "human"
    ))
    if (nrow(ereact_df) > 0) {
      write.csv(ereact_df, file.path(TAB_DIR, "Table6_Reactome_Common.csv"),
                row.names = FALSE)
    }
    cat(sprintf("  Reactome 共同靶标: %d pathways\n", nrow(ereact_df)))
  }
} else {
  cat(sprintf("  共同靶标仅%d个,跳过富集分析\n", length(targets_common)))
}

# ================================================================
# PART 11: 完成汇总
# ================================================================

cat("--- 三组分析结果 ---\n")
cat(sprintf("  PT组:  %d cells (high) / %d cells (low), nc_nCells=%d\n",
            ncol(counts_PT_high_filt), ncol(counts_PT_low_filt), nCells_PT))
cat(sprintf("  LM组:  %d cells (high) / %d cells (low), nc_nCells=%d\n",
            ncol(counts_LM_high_filt), ncol(counts_LM_low_filt), nCells_LM))
cat(sprintf("  ALL组: %d cells (high) / %d cells (low), nc_nCells=%d\n",
            ncol(counts_ALL_high_filt), ncol(counts_ALL_low_filt), nCells_ALL))

cat("\n--- 靶标分类 (FDR < 0.05) ---\n")
cat(sprintf("  PT  显著差异调控靶标: %d\n", length(targets_PT)))
cat(sprintf("  LM  显著差异调控靶标: %d\n", length(targets_LM)))
cat(sprintf("  共同靶标 (PT∩LM):     %d\n", length(targets_common)))
cat(sprintf("  PT特异靶标:           %d\n", length(targets_PT_only)))
cat(sprintf("  LM特异靶标:           %d\n", length(targets_LM_only)))

cat("\n--- 输出文件 ---\n")
cat(sprintf("  表格目录: %s\n", TAB_DIR))
cat(sprintf("  图表目录: %s\n", FIG_DIR))
cat("\n========== 完成 ==========\n")
cat(sprintf("脚本总用时: %.1f min\n",
            as.numeric(difftime(Sys.time(), TMR_SCRIPT, units = "mins"))))
cat(sprintf("参数: Top=%d genes, nc_nNet=10, nc_nCells=500(上限), nCores=4\n", N_TOP))
