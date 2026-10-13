# =============================================================================
# ST01_read_qc.R
# 功能: 读取 8 个 Visium ST 样本, SCTransform 归一化, PCA 降维, 核对 spot 数与原文
# 定位: Wu et al. Cancer Discovery 2022 ST 复现流水线 Step 1-2 (Fig 1F 前置)
# 可视化: SpotCount_QC.pdf (实测 vs 原文 spot 数柱状图)
# =============================================================================

# --- 包加载(统计 -> 数据 -> 可视化 -> 领域) --------------------------------
library(Seurat)   # 空间转录组对象与 SCTransform/PCA
library(Matrix)   # 稀疏矩阵读写
library(ggplot2)  # 数据可视化

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

sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf(file.path(FIG_DIR, paste0(filename, ".pdf")), width = width, height = height)
  print(plot_obj)
  dev.off()
}

# --- 目录结构 ---------------------------------------------------------------
BASE_DIR <- "G:/CRLM/ST/scCRLM"
DATA_DIR <- file.path(BASE_DIR, "raw/ST")
FIG_DIR  <- "G:/CRLM/ST/scCRLM/results/01_read_qc/Figure"
TAB_DIR  <- "G:/CRLM/ST/scCRLM/results/01_read_qc/table"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# ST 样本清单(论文:240-242, 8 样本 = P1-P4 的 liver+colon, 2 untreated + 2 NAC-PR)
sample.id <- c("ST-P1_liver","ST-P1_colon","ST-P2_liver","ST-P2_colon",
               "ST-P3_liver","ST-P3_colon","ST-P4_liver","ST-P4_colon")
# 注: 论文 Methods "ST" 段(L278)明确给出 8 样本逐 spot 数(如 ST-P1 colon=3,313,
#   ST-P1 liver=3,826 等), 故下方 PART 2 做"实测 vs 论文报告"逐样本比对(验证数据同一性)。

# =============================================================================
# PART 1: 读取 8 个 Visium ST 样本并 SCTransform 归一化 + PCA 降维
# 目的: 将 Space Ranger 输出读入 Seurat, 归一化后降维, 供后续聚类/注释使用
# 方法: Load10X_Spatial 读 h5 + 空间信息; SCTransform 归一化; RunPCA 取前 50 主成分
# 数据: 输入 ST_matrix2/ST/ST-Px_{liver,colon}; 输出 SCT + PCA 后的 Seurat 对象列表
# =============================================================================
cat("========== PART 1: 读取 + SCTransform + PCA ==========\n")

st.list <- list()
for (s in sample.id) {
  q <- Load10X_Spatial(file.path(DATA_DIR, s),
                       filename = "filtered_feature_bc_matrix.h5", assay = "Spatial")
  q <- SCTransform(q, assay = "Spatial", verbose = FALSE)
  q <- RunPCA(q, assay = "SCT", npcs = 50, verbose = FALSE)
  # 保存SCT对象供后续脚本复用
  saveRDS(q, file.path(TAB_DIR, paste0(s, "_seurat.rds")))
  st.list[[s]] <- q
}

# =============================================================================
# PART 2: spot 数核对(实测 vs 原文)与 QC 图
# 目的: 验证本地数据与论文 Methods ST 段(L278)报告的逐样本 spot 数一致
# 方法: 统计每样本 ncol; 柱状图对比实测与原文
# 数据: 输出 SpotCount_QC.csv (TAB_DIR) + SpotCount_QC.pdf (FIG_DIR)
# =============================================================================
cat("========== PART 2: spot 数核对 ==========\n")

# 论文 Methods ST 段(L278)报告 8 样本逐 spot 数(顺序同 sample.id: liver/colon)
paper.spot <- c(3826, 3313, 4658, 4174, 3695, 4007, 3721, 3902)  # ST-P1 liver, ST-P1 colon, ST-P2 liver, ST-P2 colon, ST-P3 liver, ST-P3 colon, ST-P4 liver, ST-P4 colon
names(paper.spot) <- sample.id
obs.spot <- sapply(st.list, ncol)
qc.df <- data.frame(sample = sample.id,
                    observed = as.integer(obs.spot),
                    paper    = paper.spot[sample.id],
                    match    = as.integer(obs.spot) == paper.spot[sample.id],
                    stringsAsFactors = FALSE)
write.csv(qc.df, file.path(TAB_DIR, "SpotCount_QC.csv"), row.names = FALSE)
cat("  spot 数核对(实测 vs 论文 L278):\n"); print(qc.df)

# 实测 vs 论文分组柱状图
qc.long <- data.frame(sample = rep(sample.id, 2),
                      group  = rep(c("Observed", "Paper (L278)"),
                                   each = length(sample.id)),
                      count  = c(as.integer(obs.spot), paper.spot[sample.id]),
                      stringsAsFactors = FALSE)
qc.plot <- ggplot(qc.long, aes(x = sample, y = count, fill = group)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(aes(label = count), position = position_dodge(width = 0.9),
            vjust = -0.3, size = 3) +
  labs(title = "ST Spot Count: Observed vs Paper-reported (Wu 2022, Methods L278)",
       x = "Sample", y = "Spot number") +
  scale_fill_manual(values = c("Observed" = "#0072B2", "Paper (L278)" = "#999999")) +
  sci_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sci_save(qc.plot, "SpotCount_QC", width = 9, height = 5)

cat("ST read + QC 完成\n")
