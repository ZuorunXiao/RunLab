# ===============================================================
# 07c_SPP1_trajectory_SPATA2.R 
# 定位: 补充 07/07b "静态 spot 空间图", 新增"沿空间路径的表达梯度"维度.
# 科学问题(与 CRLM 吻合): 正常肝/邻近肝 -> 转移瘤核心 的路径上,
#       PLAU/PLAUR/SPP1/FBXL12 及免疫亚群签名如何渐变?
# 计划状态: analysis_plan.md (§5/§7.1) 原文空间方法为 BANKSY(Fig 3)+
#       spatstat K-cross + SPARK/Moran's I, **未含 SPATA2 轨迹, 亦无 07c 步骤**.
#       故 07c 为探索增补, 非计划核心流程; 论文须单列"延伸发现".
# 方法学边界: 轨迹分析超出 Wu et al. 2022 原文方法, 属探索性增补,
#       不替代 06b(pseudobulk 区域整体差异).
# 排除项: spot 级 inferCNV (教程作者亦指存疑; 我方 CNV 在 scRNA
#       GSE178318 用 InferCNV 1.23.0 做, ST 不做).
# 轨迹特征(2026-09-10 改): 主特征 = 我们分析 §4 的 4 个核心签名
#       (SPP1_mac / PLAUR_mac / PLAU_tumor / FBXL12_tumor), 经 addFeatures 注册为
#       SPATA2 特征后 spatialTrajectoryScreening 的 variables 可直接接受(自定义程序签名名
#       须先 addFeatures 注册, 否则被 variables 校验判 invalid); 4 轴基因(FBXL12/PLAU/PLAUR/
#       SPP1, 来自 count_mtr 非负化 sct_*) 作基因级验证. 原文本免疫锚点
#       (MRC1_CCL18_Macrophages/FOXP3_Treg) 已按用户指令移除 —— 这是我们"自己的分析".
# ===============================================================

suppressMessages(library(SPATA2))
suppressMessages(library(tidyverse))

# --- 目录 ---------------------------------------------------------------
BASE    <- "X:/XX"
TAB06   <- file.path(BASE, "results/06_PLAUR_SPP1_signature/table")
TAB03   <- file.path(BASE, "results/03_immune_spatial/table")
FIG_DIR <- file.path(BASE, "results/07c_SPP1_trajectory_SPATA2/Figure")
TAB_DIR <- file.path(BASE, "results/07c_SPP1_trajectory_SPATA2/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 常量 ---------------------------------------------------------------
SAMPLE_ID <- c("ST-P1_liver", "ST-P1_colon", "ST-P2_liver", "ST-P2_colon")
SIG_FEATS  <- c("SPP1_mac", "PLAUR_mac", "PLAU_tumor", "FBXL12_tumor")   # §4 hallmark 签名(可视化主特征)
GENE_FEATS <- c("FBXL12", "PLAU", "PLAUR", "SPP1")                       # 4 轴基因(轨迹筛选稳健变量)
TRAJ_NAME <- "normal_to_tumor"           # 轨迹 id (非肿瘤质心 -> 肿瘤质心)

# 高对比可视化配色
TUMOR_GROUP_COLORS <- c(normal = "#1F78B4", tumor = "#E31A1C")  # 蓝 vs 红, 与 07  violin 配色一致且区分度高
FEAT_COLORS <- c(
  SPP1_mac      = "#E41A1C",   # red   (巨噬 SPP1, hallmark)
  PLAUR_mac     = "#377EB8",   # blue  (巨噬 PLAUR, 招募受体)
  PLAU_tumor    = "#FF7F00",   # orange(肿瘤 PLAU, 效应配体)
  FBXL12_tumor  = "#984EA3"    # purple(肿瘤 FBXL12, 癌状态)
)
# 热图填充色: 与 07 空间热力图统一, 黄-橙-红 连续色阶 (低=浅黄, 高=深红; 无蓝/紫)
HEATMAP_COLORS <- colorRampPalette(c("#FFFFCC", "#FED976", "#FEB24C",
                                      "#FD8D3C", "#FC4E2A", "#E31A1C",
                                      "#BD0026", "#800026"))(200)

cat("========== 07c: SPATA2 spatial trajectory (exploratory supplement, v3 API) ==========\n")
master <- readRDS(file.path(TAB06, "signature_master.rds"))
sct_cols <- grep("^sct_", colnames(master), value = TRUE)
gene_names <- sub("^sct_", "", sct_cols)

for (s in SAMPLE_ID) {
  d <- master[master$sample == s & !is.na(master$x), ]
  if (nrow(d) == 0) { cat("WARN: no coords for", s, "\n"); next }
  s_lab <- gsub("_", " ", s)

  # count_mtr: 非必要但 initiateSpataObject 要求; 基因不绘制, 仅作占位.
  M <- Matrix::Matrix(as.matrix(t(d[, sct_cols, drop = FALSE]))); rownames(M) <- gene_names; colnames(M) <- d$spot
  M[M < 0] <- 0
  coords_df <- data.frame(barcodes = d$spot, x = d$x, y = d$y)
  obj <- initiateSpataObject(sample_name = s, modality = "gene",
                              count_mtr = M, coords_df = coords_df)

  feat <- data.frame(
    barcodes        = d$spot,
    # 可视化用: NA = 该 spot 不属于目标子集, 对展示等价于"无签名活性", 故填 0.
    # (统计筛选仍用 GENE_FEATS, 不受此填充影响.)
    SPP1_mac        = ifelse(is.na(d$SPP1_mac), 0, d$SPP1_mac),
    PLAUR_mac       = ifelse(is.na(d$PLAUR_mac), 0, d$PLAUR_mac),
    PLAU_tumor      = ifelse(is.na(d$PLAU_tumor), 0, d$PLAU_tumor),
    FBXL12_tumor    = ifelse(is.na(d$FBXL12_tumor), 0, d$FBXL12_tumor),
    tumor_group     = factor(d$is_tumor, c(FALSE, TRUE), c("normal", "tumor")),
    stringsAsFactors = FALSE
  )
  obj <- addFeatures(obj, feature_df = feat,
                     feature_names = c(SIG_FEATS, "tumor_group"), overwrite = TRUE)

  # ============ Step A: 定义轨迹 (代码式, 非交互, 可审计) ============
  # 以 非肿瘤 spot 质心 -> 肿瘤 spot 质心 为直线轨迹, 对应
  # "正常肝/邻近肝 -> 转移瘤核心" 的科学问题. 免交互 shiny, 符合可审计偏好.
  # 注(addSpatialTrajectory 坐标解析): 官方文档参数说明 start/end 为
  #   "numeric vectors of length two" (x,y); 示例 "mm" 单位串用于已 scale 的
  #   Visium 对象. 本对象 coords 为像素且无 scale_factors, 故传数值像素坐标.
  #   安装完成后首次运行前, 建议 dump addSpatialTrajectory 源码确认坐标解析
  #   (是否强制单位后缀); 若需单位, 改成 c(x="<X>px", y="<Y>px") 形式.
  nrm <- d[!d$is_tumor, ]; tum <- d[d$is_tumor, ]
  if (nrow(nrm) == 0 || nrow(tum) == 0) {
    cat("SKIP", s, ": need both normal and tumor spots for trajectory.\n"); next
  }
  start_xy <- c(x = mean(nrm$x), y = mean(nrm$y))
  end_xy   <- c(x = mean(tum$x),  y = mean(tum$y))
  obj <- addSpatialTrajectory(object = obj, id = TRAJ_NAME,
                              start = start_xy, end = end_xy, overwrite = TRUE)
  # 备选(交互式, 适合弯曲路径): obj <- createSpatialTrajectories(object = obj); 须 <- 存回.

  if (!TRAJ_NAME %in% getSpatialTrajectoryIds(obj)) {
    cat("SKIP", s, ": trajectory '", TRAJ_NAME, "' not defined.\n", sep = ""); next
  }

  # ============ Step B1: 轨迹线 + 肿瘤分组着色(高对比蓝/红) ============
  p <- plotSpatialTrajectories(object = obj, ids = TRAJ_NAME,
                               color_by = "tumor_group",
                               pt_alpha = 0.35, pt_alpha2 = 0.95) +
    ggpLayerSpatialTrajectories(object = obj, ids = TRAJ_NAME) +
    scale_color_manual(values = TUMOR_GROUP_COLORS, name = "Region") +
    theme_bw() + theme(legend.position = "top",
      axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +
    labs(title = paste0("Trajectory | ", s_lab))
  pdf(file.path(FIG_DIR, paste0("Trajectory_", s, ".pdf")), width = 7, height = 7,
      family = "Helvetica"); print(p); dev.off()
  cat("saved Trajectory_", s, ".pdf\n", sep = "")

  # ============ Step B2: 特征沿轨迹趋势线 (我们 §4 核心签名, 高对比分色) ============
  p <- plotStsLineplot(object = obj, variables = SIG_FEATS, id = TRAJ_NAME,
                       nrow = 2) +
    scale_color_manual(values = FEAT_COLORS, name = "Signature") +
    theme_bw() + theme(legend.position = "top")
  pdf(file.path(FIG_DIR, paste0("TrajectoryFeatures_", s, ".pdf")), width = 9, height = 7,
      family = "Helvetica"); print(p); dev.off()
  cat("saved TrajectoryFeatures_", s, ".pdf\n", sep = "")

  # ============ Step B3: 轨迹筛选(4 轴基因, SPATA2 稳健变量) ============
  # 用 GENE_FEATS(4 轴基因)做显著性/最佳模型筛选: spatialTrajectoryScreening 的 loess
  # 梯度拟合对 addFeatures 连续签名会 x 轴退化崩溃, 对词表基因稳健. 签名级梯度见 B2/B4 可视化.
  # rm_zero_infl=FALSE: 跳过内置零膨胀剔除(对全平变量返回 NA 触发 select 崩溃).
  sts <- spatialTrajectoryScreening(object = obj, id = TRAJ_NAME,
                                    variables = GENE_FEATS, rm_zero_infl = FALSE)
  # 显著性表: fdr < 0.05 视为沿轨迹显著变化 (slot 名依官方 trajectory-screening 教程)
  sig_df <- sts@results$significance
  write.csv(sig_df, file.path(TAB_DIR, paste0("TrajectoryScreening_", s, ".csv")),
            row.names = FALSE)
  sig_vars <- sig_df[which(sig_df$fdr < 0.05), "variables"]
  cat(s, " significant (fdr<0.05):", paste(sig_vars, collapse = ","), "\n")
  # 每特征最佳拟合模型 (按 mae 最小); 升/降依 models 列 ^ascending/^descending 判定.
  # model_fits 在无显著/无模型拟合时返回 0x0 空 data.frame(无 variables 列), 须防御.
  fits <- sts@results$model_fits
  if ("variables" %in% colnames(fits)) {
    write.csv(fits, file.path(TAB_DIR, paste0("TrajectoryModelFits_", s, ".csv")), row.names = FALSE)
  }
  if ("variables" %in% colnames(fits) && nrow(fits) > 0 && length(sig_vars) > 0) {
    best <- fits %>% filter(variables %in% sig_vars) %>%
      group_by(variables) %>% slice_min(mae, n = 1)
    asc  <- best %>% filter(grepl("^ascending",  models, ignore.case = TRUE)) %>% pull(variables)
    desc <- best %>% filter(grepl("^descending", models, ignore.case = TRUE)) %>% pull(variables)
    cat(s, " ascending model:", paste(asc, collapse = ","), "\n")
    cat(s, " descending model:", paste(desc, collapse = ","), "\n")
  } else {
    cat(s, " ascending model: (none)\n")
    cat(s, " descending model: (none)\n")
  }

  # ============ Step B4: 轨迹梯度热图(4 轴基因; 签名级 loess 在 heatmap 内崩溃, 故用稳健基因) ============
  p <- plotStsHeatmap(object = obj, variables = GENE_FEATS, id = TRAJ_NAME,
                      arrange_rows = "maxima", smooth_span = 0.5) +
    scale_fill_gradientn(colors = HEATMAP_COLORS, name = "Scaled\nexpression") +  # 覆盖默认 SPATA2 色阶, 统一为黄-橙-红
    theme_bw()
  pdf(file.path(FIG_DIR, paste0("TrajectoryHeatmap_", s, ".pdf")), width = 8, height = 6,
      family = "Helvetica"); print(p); dev.off()
  cat("saved TrajectoryHeatmap_", s, ".pdf\n", sep = "")
}

