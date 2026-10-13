#!/usr/bin/env Rscript
# =============================================================================
# 02_unsupervised_cluster.R
# 功能: 8 个 ST 样本无监督聚类 (FindClusters), 输出空间聚类图 (Fig 1F)
# 定位: Wu et al. Cancer Discovery 2022 ST 复现 Step 3 (Fig 1F)
# 方法: 严格按原文 Methods/Results: SCTransform -> RunPCA -> FindNeighbors ->
#       FindClusters -> RunUMAP (原文 Fig 1F 基础: "Unsupervised clustering
#       analysis classified the samples into distinct regions such as tumor
#       and fibroblast regions")。本步只产出 spot 聚类编号; 逐 cluster 的
#       生物学注释(肿瘤/成纤维/内皮/肝细胞 + 免疫亚群)由 03c 用论文同款
#       "标记基因判区域 + label transfer 填免疫" 逻辑合并为 Fig 1F 全注释图
# 可视化: ClusterMap_STP*.pdf (按 seurat_clusters 着色的空间聚类图, 未注释)
# =============================================================================

# --- 包加载(统计 -> 数据 -> 可视化 -> 领域) --------------------------------
library(Seurat)
library(Matrix)
library(ggplot2)

# --- 全局辅助函数(内联, 规范22 自给自足) ----------------------------------
sci_theme <- theme_bw(base_size = 12) +
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

sci_save <- function(plot_obj, filename, width = 10, height = 8) {
  pdf(file.path(FIG_DIR, paste0(filename, ".pdf")), width = width, height = height)
  print(plot_obj)
  dev.off()
}

# --- 目录结构 ---------------------------------------------------------------
BASE_DIR <- "G:/CRLM/ST/scCRLM"
TAB01_DIR <- file.path(BASE_DIR, "results/01_read_qc/table")
FIG_DIR  <- "G:/CRLM/ST/scCRLM/results/02_unsupervised_cluster/Figure"
TAB_DIR  <- "G:/CRLM/ST/scCRLM/results/02_unsupervised_cluster/table"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

sample.id <- c("ST-P1_liver","ST-P1_colon","ST-P2_liver","ST-P2_colon",
               "ST-P3_liver","ST-P3_colon","ST-P4_liver","ST-P4_colon")
# D-H: 论文 Methods 未给 FindClusters resolution / PC 数, 属透明校准参数。
#   原校准(res=1.0, 称"覆盖全部 8 区")依据旧版 z-score argmax 扫描, 已废弃且证据缺失 -> 注销。
#   本回合用**新正交 ssGSEA 面板**重跑 _scan_resolution.R(RES_LIST=1.0/1.5/2.0/2.5)数据驱动校准,
#   结果(_scan_region_map.csv, 2026-09-07):
#     - 区域集合随分辨率几乎不变, 仅细分; 分辨率越高簇越多但区域数不增 -> 调细不"创造"新区。
#     - 赋值置信 margin 随分辨率升高而下降(res>=2.0 明显碎裂: 如 P4_colon margin 96.5->29.7,
#       P1_colon 1427->221), 故 res>=2.0 过度碎裂。
#     - 最佳折中 = res=1.5: 比 1.0 多恢复 Smooth muscle 等亚区(如 P4_colon 1.0 仅 3 区 -> 1.5 得
#       Normal epithelium + Lamina propria + Smooth muscle + Tumor 共 4 区, margin 仍健康 73.4), 且未碎裂。
#     - 局限: 单样本区域数上限仅 4(P4_colon); P1_colon 全 Tumor / P2-P4_liver 全 Hepatocytes 等
#       单区域样本在所有分辨率下均不恢复其他区 -> 此非分辨率问题, 需 spot 级 ssGSEA(方案 b)或重检数据。
#   - 注: 2026-09-07 用户显式指令将分辨率改为 0.8 (更粗, 偏少簇); 上方 1.5 最佳折中结论保留为历史校准依据,
#       本回合用 res=0.8。
cluster.res <- 0.5   # 2026-09-07 用户指令: 从 0.8 改为 0.5 (更粗聚类, 簇更少)

# =============================================================================
# PART 1: 读取 + SCTransform + PCA + 聚类 + UMAP
# 目的: 每样本无监督聚类, 复现 Fig 1F 的 spot 分区
# 方法: Load10X_Spatial -> SCTransform -> RunPCA -> FindNeighbors/FindClusters -> RunUMAP
# 数据: 输入 raw/ST/ST-Px; 输出带 clusters + umap 的 Seurat 对象列表
# =============================================================================
cat("========== PART 1: 读取 + 聚类 ==========\n")

# D-E: 直接读取 01 已处理对象(SCTransform + RunPCA 的 SCT 表示), 不再重复
#      从 raw 加载与 SCTransform, 保证 02 与 03 的 ST 归一化完全一致的 SCT。
st.list <- list()
for (s in sample.id) {
  rds <- file.path(TAB01_DIR, paste0(s, "_seurat.rds"))
  if (!file.exists(rds)) stop("缺少 ", rds, " —— 请先跑完 01_read_qc.R")
  q <- readRDS(rds)   # 01 已 SCTransform + RunPCA(SCT), PCA 降维在此复用
  q <- FindNeighbors(q, reduction = "pca", dims = 1:50, verbose = FALSE)
  q <- FindClusters(q, resolution = cluster.res, verbose = FALSE)
  q <- RunUMAP(q, reduction = "pca", dims = 1:50, verbose = FALSE)
  st.list[[s]] <- q
}

# =============================================================================
# PART 2: 空间聚类图 (Fig 1F 交付物)
# 目的: 将每样本聚类映射回 HE 切片
# 方法: SpatialDimPlot 按 seurat_clusters 着色, 显示 cluster 编号 label
# 数据: 输出 ClusterMap_STP*.pdf (FIG_DIR), 画布10x8英寸
# =============================================================================
cat("========== PART 2: 空间聚类图 ==========\n")

for (s in sample.id) {
  p <- SpatialDimPlot(st.list[[s]], group.by = "seurat_clusters", label = TRUE) +
    labs(title = paste0(s, " (res=", cluster.res, ")")) + sci_theme
  sci_save(p, paste0("ClusterMap_", s), width = 10, height = 8)
}

# =============================================================================
# PART 3: 导出聚类对象 + 逐 spot 聚类 ID (供 03c 读取并按 02 cluster 注释)
# 目的: 03c 必须严格使用 02 产出的无监督聚类结果, 避免 03c 重算 cluster
#       导致与 02 ClusterMap 不一致。本步保存 Seurat 对象列表与 cluster_id。
# 数据: obj/SeuratList_Clusters.rds; table/cluster_id.csv
# =============================================================================
cat("========== PART 3: 导出聚类对象 ==========\n")
OBJ_DIR <- file.path(BASE_DIR, "results/02_unsupervised_cluster/obj")
dir.create(OBJ_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(st.list, file.path(OBJ_DIR, "SeuratList_Clusters.rds"))

clid <- data.frame()
for (s in sample.id) {
  df <- data.frame(sample = s, spot = colnames(st.list[[s]]),
                   seurat_clusters = st.list[[s]]$seurat_clusters,
                   stringsAsFactors = FALSE)
  clid <- rbind(clid, df)
}
write.csv(clid, file.path(TAB_DIR, "cluster_id.csv"), row.names = FALSE)
cat("Saved SeuratList_Clusters.rds + cluster_id.csv (", nrow(clid), " spots)\n")

cat("ST unsupervised clustering 完成\n")
quit(save = "no", status = 0, runLast = FALSE)
