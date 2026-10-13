# =============================================================================
# 03_immune_spatial.R
# 功能: Seurat label transfer 注释 ST spot 主细胞型 + ssGSEA 评分全部亚细胞型
# 定位: Wu et al. Cancer Discovery 2022 ST 复现 Step 4-5 (纯 R, 无 Python)
# 方法: 严格按原文 Methods (L725-737):
#       (1) SCTransform->RunPCA 后, FindTransferAnchors/TransferData 由 scRNA
#           参考预测 ST spot 主细胞型 (main cell type);
#       (2) GSVA::ssgsea 对 00 构建的全部亚细胞型基因集打分 (subcell types)
# 可视化: MainCellType_STP*.pdf (Fig 1G, label transfer 主型);
#        MacrophageScore_STP*.pdf (Fig 2C, 巨噬细胞亚群富集);
#        Macrophage_*_ByTreatment.pdf (治疗组间箱线图)
# 数据: 输入 raw/ST + 00 ref_maincell_seurat.rds + gene_sets.csv
# 参考: Seurat 官方 SCTransfer(FindTransferAnchors/TransferData); GSVA 官方 ssGSEA
# 版本(D-G): 论文 Seurat V3.2.2; 本机 Seurat 为 V5。SCT label transfer 按 V5 SCTransfer
#       工作流(参考与查询同为 SCT)实现, 行为基本一致, 结果可能略有差异。
# 注: 本地 scRNA 参考为 CD45+ 免疫-only (main_cell_type 仅 9 免疫主型), 与论文一致——
#     论文 scRNA 处理段(L286)逐字 "We performed the clustering analysis of all cells and
#     only reserved the CD45+ clusters", 故论文 ST label transfer 参考亦为免疫-only, 仅产出
#     免疫主型。Fig 1F 的结构区(肿瘤/肝细胞/成纤维/内皮/平滑肌)来自无监督聚类 + 标记基因
#  [D2] 参考质量依赖 raw/reference_signatures/sc.h5ad 是否已按论文 (Harmony 整合 + SingleR +
#      marker genes, L1740-1746 / L220-222) 注释: 本脚本直接消费 00 的 ref_maincell_seurat.rds,
#      不重做注释也不校验。若 h5ad 出处非论文注释, 整条参考链失真——需用户手动核实 h5ad 来源。
#     (Results L80 + S2), 非 label transfer; 本参考已忠实复现论文, 无需"含癌细胞全参考"。
# =============================================================================

# --- 包加载(统计 -> 数据 -> 可视化 -> 领域) --------------------------------
suppressMessages({
  library(Seurat)    # 空间对象 / label transfer
  library(GSVA)      # ssGSEA 富集评分
  library(ggplot2)   # 可视化
  library(patchwork) # 拼图
})

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
  pdf(file.path(FIG_DIR, paste0(filename, ".pdf")), width = width, height = height,
      family = "Helvetica")
  print(plot_obj)
  dev.off()
}

# --- 目录结构 ---------------------------------------------------------------
BASE_DIR   <- "X:/XX"
REF_DIR    <- file.path(BASE_DIR, "results/00_build_reference")
TAB01_DIR  <- file.path(BASE_DIR, "results/01_read_qc/table")
FIG_DIR    <- file.path(BASE_DIR, "results/03_immune_spatial/Figure")
TAB_DIR    <- file.path(BASE_DIR, "results/03_immune_spatial/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

sample.id <- c("ST-P1_liver","ST-P1_colon","ST-P2_liver","ST-P2_colon",
               "ST-P3_liver","ST-P3_colon","ST-P4_liver","ST-P4_colon")
treatment_map <- c("P1" = "untreated", "P2" = "untreated",
                   "P3" = "NAC-PR",    "P4" = "NAC-PR")
# [D5] 9 主型含独立 "B"; 论文 Fig 1B (L1520) 主型示意 8 类 (CD4 CD8 MAIT Myeloid Neu NK
#      Plasma Treg) 未单列 B (B 标志 CD79A/MZB1 见 L1741-1746, IGHA+/IGHG+ plasma B 见 L1243)。
#      若论文参考将 B 并入 Plasma/其他, 应从 ref 与本节去掉 "B"; 当前保留供参考(多余项无害)。
colors_main <- c("B" = "#1f77b4", "CD4" = "#ff7f0e", "CD8" = "#2ca02c",
                 "MAIT" = "#d62728", "Myeloid" = "#9467bd", "NEU" = "#8c564b",
                 "NK" = "#e377c2", "Plasma" = "#7f7f7f", "Treg" = "#bcbd22")
n.pca <- 30   # D-H: 论文未规定 label transfer 的 PC 数, 本机 V5 本地默认(取 SCT PCA 前 30 维); 属透明校准参数

# =============================================================================
# PART 1: 读取参考对象 + 基因集
# 目的: 从 00 导出的 RDS 重建 label transfer 参考; 读取 ssGSEA 基因集
# 方法: readRDS(ref_maincell_seurat.rds); split(gene_sets.csv)
# =============================================================================
cat("========== PART 1: 读取参考对象 + 基因集 ==========\n")
ref <- readRDS(file.path(REF_DIR, "ref_maincell_seurat.rds"))
cat("参考:", ncol(ref), "cells x", nrow(ref), "genes; PCA:", ncol(ref[["pca"]]), "\n")
cat("参考 main_cell_type:", paste(sort(unique(ref$main_cell_type)), collapse = "/"), "\n")

gs_df <- read.csv(file.path(REF_DIR, "gene_sets.csv"), stringsAsFactors = FALSE)
gene_sets <- split(gs_df$gene, gs_df$gene_set)
paper_gs <- file.path(REF_DIR, "paper_genesets.csv")
if (file.exists(paper_gs)) {
  pg <- read.csv(paper_gs, stringsAsFactors = FALSE)
  gene_sets <- split(pg$gene, pg$gene_set)
  cat("使用论文亚群标志基因集 paper_genesets.csv:", length(gene_sets), "个 (D-D 已对齐原文)\n")
} else {
  cat("基因集(亚细胞型, 本地 top-100 启发式):", length(gene_sets), "\n")
}
mac_names <- grep("Macrophage", names(gene_sets), value = TRUE)
cat("巨噬细胞基因集:", paste(mac_names, collapse = " | "), "\n")

# =============================================================================
# PART 2: Label Transfer (每样本)
# 目的: 将参考主细胞类型标签迁移至 ST spots (原文 main cell type 预测, :1845/:246)
# 方法: 读取 01 的 SCT 对象(SCT 表示) -> FindTransferAnchors(reference=ref[SCT], query=st[SCT], normalization.method="SCT")
#       FindTransferAnchors(reference=ref[SCT], query=st[SCT], normalization.method="SCT")
#       -> TransferData(main_cell_type)
# 说明: 论文 ST 侧 SCTransform(:1841-1843); Seurat SCTransfer 要求参考与查询同为
#       SCT, 故 00 参考与 03 查询均用 SCT(见 00/03 注释), 与论文一致。
# =============================================================================
cat("\n========== PART 2: Label Transfer ==========\n")
all_obj <- list()
for (s in sample.id) {
  cat("  处理", s, "...\n")
  # D-E: 读取 01 已处理对象(SCTransform + RunPCA 的 SCT 表示), 不再重复从 raw 加载,
  #      保证与 02 一致的 SCT 归一化。Seurat V5 下 SCTransform(assay="Spatial")
  #      输出 assay 名为 "SCT"(非 "Spatial_SCT"), 与 00 参考的 "SCT" 对齐, 供 transfer。
  rds <- file.path(TAB01_DIR, paste0(s, "_seurat.rds"))
  if (!file.exists(rds)) stop("缺少 ", rds, " —— 请先跑完 01_read_qc.R")
  st <- readRDS(rds)   # assays: Spatial(原始 counts) + SCT(SCTransform 残差); 
  sp_counts <- LayerData(st, assay = "Spatial", layer = "counts")
  if (!is.null(sp_counts)) LayerData(st, assay = "SCT", layer = "counts") <- sp_counts
  # 论文 ST 经 SCTransform 归一化 (L1842), 全程统一用 SCT 表示: label transfer 用 query.assay="SCT"(见下),
  # ssGSEA 也直接在 SCT assay 的归一化残差上打分 (见 PART3, 与论文 "SCTransform normalize -> ssGSEA" 一致)。
  # 故不再对 Spatial assay 做独立 NormalizeData (旧版 log-normalized 用于 ssGSEA 与论文 SCTransform 不一致, 已移除)。
  # 注: SCT 残差即论文 "normalized data"; 若后续更偏好 log-normalized 表达打分, 改回 GetAssayData(assay="Spatial") 即可。


  # SCT label transfer(论文 L308 / L80 "main immune subpopulations";
  # ST 侧 SCTransform, :1841-1843)。Seurat SCTransfer 工作流: 参考与查询同为 SCT。
  # reduction="pcaproject" 将 query 投影到参考 SCT PCA; weight.reduction 默认
  # "pcaproject"(投影后联合嵌入)。切勿手动传 ref[["pca"]] —— 源码(TransferData
  # L60-73)会把它改名 *_query, 与 anchor set 的 query 细胞名不对齐导致下标越界。
  anchors <- FindTransferAnchors(reference = ref, query = st,
                                 normalization.method = "SCT",
                                 reference.assay = "SCT", query.assay = "SCT",
                                 dims = 1:n.pca, verbose = FALSE)
  # refdata 必须传"参考对象中的注释向量"(命名列表), 不可传字符串(否则
  # 报 "None of the provided refdata elements are valid")。
  # 必须传 query = st: 否则 TransferData 返回 data.frame 而非带注释的 Seurat
  # 对象(源码 L213-218, query=NULL 时返回预测 data.frame)。传 query=st 后,
  # 预测列名为 predicted.main_cell_type(源码 L174: paste0("predicted.",rd.name))。
  st <- TransferData(anchorset = anchors,
                     refdata = list(main_cell_type = ref$main_cell_type),
                     query = st, query.assay = "SCT", dims = 1:n.pca)
  # TransferData 预测列命名: 传 query=st 时为 predicted.main_cell_type
  # (源码 L174)。兼容旧命名 main_cell_type.predicted / main_cell_type。
  pred_col <- c("predicted.main_cell_type", "main_cell_type.predicted", "main_cell_type")[
    c("predicted.main_cell_type", "main_cell_type.predicted", "main_cell_type") %in% colnames(st@meta.data)]
  if (length(pred_col) == 0) stop(paste("样本", s, "TransferData 未生成预测列"))
  colnames(st@meta.data)[colnames(st@meta.data) == pred_col[1]] <- "pred_main"
  # 捕获 label transfer 分类置信度(最大类后验概率), 作参考列写入 CSV;
  # 论文未规定置信阈值, 该列不参与 03c 结构区判类(方案 A 已删免疫叠加规则), 仅供下游核查。
  # 注意: 本机 Seurat v5 的分数列名为 "predicted.main_cell_type.score"(实证: 单样本
  #   TransferData 后 meta.data 实测列名), 非旧版 "prediction.score.main_cell_type";
  #   故以正则 "predicted\\..*\\.score" 匹配, 避免依赖版本特定前缀。
  score_cands <- grep("predicted\\..*\\.score", colnames(st@meta.data), value = TRUE)
  if (length(score_cands)) {
    colnames(st@meta.data)[colnames(st@meta.data) == score_cands[1]] <- "pred_score"
  } else {
    st@meta.data[["pred_score"]] <- rep(NA_real_, ncol(st))
  }
  cat("    ", s, "预测主细胞类型分布:",
      paste(names(sort(table(st$pred_main, useNA = "no"), decreasing = TRUE)),
            collapse = "/"), "\n")
  all_obj[[s]] <- st
  # 主细胞类型空间注释图 (论文: Seurat label transfer 预测主型): 离散着色。
  # 注: 本地 scRNA 参考为 CD45+ 免疫-only, 故该图仅显示免疫主型; 非免疫区
  # (肿瘤/成纤维/内皮/肝细胞) 由 03c 基于 02 聚类 + marker 判区合并呈现。
  # [D7] label transfer 用免疫-only 参考, 每个 spot 均被强制判为某免疫主型(最近邻);
  #      非免疫 spot(肿瘤/肝细胞)的 pred_main 为强制标签, 可能误导。此为免疫-only 参考的
  #      固有限制(论文同构: 参考仅 CD45+), 03c 用 02 聚类补结构区。论文未规定置信阈值,
  #      故不自行加阈值过滤(避免引入自创门控); 如需严格可据 pred_score 过滤低置信 spot。
  p <- SpatialDimPlot(st, group.by = "pred_main", cols = colors_main,
                      label = FALSE, image.alpha = 1) +
    labs(title = paste0(s, " - predicted main immune cell type (label transfer)")) +
    sci_theme + theme(legend.position = "right", legend.title = element_blank())
  sci_save(p, paste0("MainCellType_", s), width = 8, height = 8)
}

# =============================================================================
# PART 3: ssGSEA 全部亚细胞型评分 (GSVA)
# 目的: 计算各 spot 的全部亚细胞型富集评分 (原文 subcell types 评分)
# 方法: gsva(SCT 归一化残差, geneSets, method="ssgsea", kcdf="Gaussian")  # D3: 改 SCT 残差(论文 L1842 SCTransform 归一化); 原 log-normalized 已弃用
# =============================================================================
ssgsea_wide <- data.frame()
for (s in sample.id) {
  cat("  评分", s, "...\n")
  st <- all_obj[[s]]
  expr_mat <- as(GetAssayData(st, assay = "SCT", layer = "data"), "matrix")   # genes x spots, SCTransform 残差(论文 L1842 SCTransform 归一化)
  rownames(expr_mat) <- rownames(st)
  # 仅保留 expr 中存在的基因, 保证 gsva 可匹配
  gs_sub <- lapply(gene_sets, function(g) g[g %in% rownames(expr_mat)])
  gs_sub <- gs_sub[lengths(gs_sub) > 0]
  cat("    ", s, "实际评分基因集:", length(gs_sub), "\n")

  # GSVA 2.x (本机 2.4.9) 新 API: gsva() 仅接收 *Param 对象, 旧式
  # gsva(mat, geneSets, method=, kcdf=) 已不可用于 matrix(报
  # "no inherited method for signature param='matrix'")。ssGSEA 用专用
  # ssgseaParam(exprData, geneSets, alpha=0.25, normalize=TRUE), 返回
  # gene_sets x spots 矩阵(与下游 gsva_res[mn, colnames(st)] 一致)。
  # D-H: ssGSEA alpha=0.25 / normalize=TRUE 为 GSVA 默认, 论文未规定参数; 属透明本地默认。
  gsva_param <- ssgseaParam(expr_mat, geneSets = gs_sub,
                            alpha = 0.25, normalize = TRUE, verbose = FALSE)
  gsva_res <- gsva(gsva_param)                                # gene_sets x spots
  # [FIX-PART3] gs_sub 已按本样本 0 匹配基因剔除空集 -> gsva_res 行数随样本变化,
  #   直接 rbind 跨样本会因列数不一致报 "numbers of columns ... do not match"。
  #   扩展回全部 gene_sets 名(缺失集填 NA), 保证每样本 w 列数恒等, rbind 安全。
  gs_names_full <- names(gene_sets)
  res_full <- matrix(NA_real_, nrow = length(gs_names_full), ncol = ncol(gsva_res),
                     dimnames = list(gs_names_full, colnames(gsva_res)))
  res_full[rownames(gsva_res), ] <- gsva_res
  # 把全部亚细胞型评分写回 meta.data 用于空间出图
  for (mn in rownames(gsva_res)) {
    colnm <- gsub("[^A-Za-z0-9]", "_", mn)
    st@meta.data[[colnm]] <- as.numeric(gsva_res[mn, colnames(st)])
  }
  all_obj[[s]] <- st

  # 巨噬细胞亚群富集空间图 (论文 Fig 2C 等价: ssGSEA 连续打分): 多亚群横向拼图
  mac_grad <- colorRampPalette(c("#F7F7F7", "#2166AC"))(100)
  plist <- list()
  for (mn in mac_names) {
    colnm <- gsub("[^A-Za-z0-9]", "_", mn)
    plist[[mn]] <- SpatialFeaturePlot(st, features = colnm,
                              min.cutoff = "q05", max.cutoff = "q95",
                              image.alpha = 1) +
      scale_color_gradientn(colors = mac_grad) +
      labs(title = mn) + sci_theme + theme(legend.position = "right")
  }
  pcomb <- wrap_plots(plist, nrow = 1) +
    plot_annotation(title = paste0(s, " - macrophage subpopulation enrichment (ssGSEA)"))
  sci_save(pcomb, paste0("MacrophageScore_", s), width = 16, height = 5)

  # [D6] CD8_score = 论文 CD8+ T 亚群 ssGSEA 均值 (CD8A/GZMK/XCL1 等, L420/L1240/L1241), 属衍生可视化,
  #      非论文原图直接对应(论文 Fig 1G 展示 CD8+ T 富集, Fig 5B 展示 GZMK+/XCL1+ CD8+ T)。
  cd8_rows <- grep("CD8", rownames(res_full), value = TRUE)
  cd8_score <- colMeans(gsva_res[cd8_rows, , drop = FALSE])
  st@meta.data[["CD8_score"]] <- as.numeric(cd8_score[colnames(st)])
  cd8_grad <- colorRampPalette(c("#F7F7F7", "#B2182B"))(100)
  p_cd8 <- SpatialFeaturePlot(st, features = "CD8_score",
                     min.cutoff = "q05", max.cutoff = "q95",
                     image.alpha = 1) +
    scale_color_gradientn(colors = cd8_grad) +
    labs(title = paste0(s, " - CD8+ T cell score (ssGSEA)")) + sci_theme +
    theme(legend.position = "right")
  sci_save(p_cd8, paste0("CD8_score_", s), width = 8, height = 8)

  # 逐样本 ssGSEA 宽表 (spots x 亚细胞型); 用 res_full(全 gene_sets 列, 缺失 NA) 保证跨样本列数一致
  w <- data.frame(sample = s, spot = colnames(st), stringsAsFactors = FALSE)
  w <- cbind(w, as.data.frame(t(res_full)))
  ssgsea_wide <- rbind(ssgsea_wide, w)
}

# 逐 spot 主细胞型 + 亚细胞型评分表 (细胞注释交付物)
anno_tab <- data.frame()
for (s in sample.id) {
  st <- all_obj[[s]]
  df <- data.frame(
    sample    = s,
    treatment = treatment_map[substring(s, 4, 5)],
    spot      = colnames(st),
    pred_main = st$pred_main,
    pred_score = st$pred_score,
    # immune_confident: 基于 TransferData 原生最大类后验 pred_score 的透明布尔列(中点 0.5),
    # 供下游/用户按置信筛选免疫 spot; 不改 pred_main(论文免疫-only 参考强制全覆盖, 见 PART2 [D7])。
    immune_confident = st$pred_score >= 0.5,
    stringsAsFactors = FALSE
  )
  anno_tab <- rbind(anno_tab, df)
}
anno_tab <- merge(anno_tab, ssgsea_wide, by = c("sample", "spot"))
write.csv(anno_tab, file.path(TAB_DIR, "spot_celltype_annotation.csv"), row.names = FALSE)
cat("Saved spot_celltype_annotation.csv (", nrow(anno_tab), " spots x ",
    ncol(anno_tab), " cols: pred_main + ", length(gene_sets), " ssGSEA 亚细胞型)\n")
write.csv(ssgsea_wide, file.path(TAB_DIR, "ssGSEA_subcell.csv"), row.names = FALSE)
cat("Saved ssGSEA_subcell.csv (", nrow(ssgsea_wide), " spots x ",
    ncol(ssgsea_wide) - 2, " gene sets)\n")

# =============================================================================
# PART 4: 治疗组间巨噬细胞评分箱线图 (Fig 2C 佐证)
# =============================================================================
for (mn in mac_names) {
  colnm <- gsub("[^A-Za-z0-9]", "_", mn)
  bp_df <- data.frame(
    score     = ssgsea_wide[[mn]],
    Treatment = factor(treatment_map[substring(ssgsea_wide$sample, 4, 5)],
                       levels = c("untreated", "NAC-PR"))
  )
  bp <- ggplot(bp_df, aes(x = Treatment, y = score, fill = Treatment)) +
    geom_boxplot() +
    scale_fill_manual(values = c("untreated" = "#007AC2", "NAC-PR" = "#FF6B00")) +
    labs(title = paste0(mn, " Score by Treatment"),
         x = "Treatment Group", y = "ssGSEA Score") +
    sci_theme
  sci_save(bp, paste0("Macrophage_", colnm, "_ByTreatment"), width = 6, height = 5)
}

cat("\n03_immune_spatial.R 完成 (纯 R, 无 Python)\n")
quit(save = "no", status = 0, runLast = FALSE)
