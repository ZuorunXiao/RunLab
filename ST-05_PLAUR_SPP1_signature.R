# ===============================================================
# 06_PLAUR_SPP1_signature.R
# 功能: 构建 PLAU-PLAUR-SPP1-FBXL12 分析的逐 spot 主表(4 个非治疗样本)
# 定位: 分析方案 W1 阶段0 门禁 + 全系列脚本的数据基座(后续 07-13 消费本表)
# 可视化: GeneDetectionRate.pdf (4 基因逐样本检出率)
#         CellTypeComposition.pdf (每样本免疫主型占比, 对应原文 Fig 1G)
#         GeneExpressionDistribution.pdf (4 基因按区室箱线图)
#         SignatureDistribution.pdf (4 签名按区室箱线图)
#         MyeloidMaskValidation.pdf (数据驱动髓系掩码纯度验证)
# 数据: 输入 results/01_read_qc/table/ST-*_seurat.rds
#            results/03c_full_annotation/table/full_annotation.csv (region)
#            results/03_immune_spatial/table/spot_celltype_annotation.csv
#       输出 results/06_PLAUR_SPP1_signature/table/signature_master.csv (+ .rds)
#            results/06_PLAUR_SPP1_signature/table/GeneDetectionRate.csv
#            results/06_PLAUR_SPP1_signature/table/MyeloidMaskValidation.csv
# 依据: 方案 analysis_plan.md 第 4 节(4 签名定义) 与第 2 节(仅非治疗样本)
# 官方 API 实测(Seurat 5.5.1 / SeuratObject 5.4.0):
#   - FetchData() 的 slot= 参数在 SeuratObject 5.0.0 起 defunct, 必须用 layer=
#   - GetTissueCoordinates() 返回 data.frame(x, y, cell), 行名=barcode,
#     顺序与 colnames(seu) 完全一致
# ===============================================================

# --- 包加载(统计 -> 数据 -> 可视化 -> 领域) -----------------------------
suppressMessages({
  library(Seurat)      # 空间对象读取与表达提取
  library(Matrix)      # 稀疏矩阵支持
  library(ggplot2)     # 可视化
  library(reshape2)    # 长宽格式转换(melt)
})

# --- 全局辅助函数(内联, 规范22 自给自足, 禁 source()) -------------------
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
  pdf(file.path(FIG_DIR, paste0(filename, ".pdf")), width = width,
      height = height, family = "Helvetica")
  print(plot_obj)
  dev.off()
}

# --- 目录结构 -----------------------------------------------------------
BASE_DIR <- "X:/XX"
TAB01    <- file.path(BASE_DIR, "results/01_read_qc/table")
TAB03    <- file.path(BASE_DIR, "results/03_immune_spatial/table")
TAB03C   <- file.path(BASE_DIR, "results/03c_full_annotation/table")
FIG_DIR  <- file.path(BASE_DIR, "results/06_PLAUR_SPP1_signature/Figure")
TAB_DIR  <- file.path(BASE_DIR, "results/06_PLAUR_SPP1_signature/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 分析常量(方案锁定, 执行前不得改动) ---------------------------------
# 仅非治疗样本: 依据 03_immune_spatial.R L66-67 treatment_map =
#   c(P1="untreated", P2="untreated", P3="NAC-PR", P4="NAC-PR")
SAMPLE_ID <- c("ST-P1_liver", "ST-P1_colon", "ST-P2_liver", "ST-P2_colon")
GENE_SET  <- c("FBXL12", "PLAUR", "PLAU", "SPP1")
# 4 个签名列名(由本脚本构造, 供 07-13 消费): FBXL12/PLAU 限定肿瘤,
#   PLAUR/SPP1 限定巨噬, 子集外置 NA 以保留空间特异性
SIG_SET   <- c("FBXL12_tumor", "PLAU_tumor", "PLAUR_mac", "SPP1_mac")
# --- 数据驱动髓系掩码(方案 v1.1) ----------------------------------------
# 背景: 03 的 label transfer 参考集为 CD45+ 免疫细胞专属(9 类 x 每类 3000
#   强制均衡), 无上皮/肿瘤/肝细胞选项, 故每个 spot 被强制归入某免疫标签。
#   外部验证(Myeloid ssGSEA 的 AUC)显示 pred_main=="Myeloid" 基本不识别
#   髓系 spot: Adjacent colon 0.562(ns) / CRC 0.453(反向) /
#   Adjacent liver 0.514(ns) / LM 0.595(弱)。故弃用 pred_main 定义巨噬。
# 掩码: 泛巨噬标志基因的 SCT 表达 -> 样本内逐基因 z-score -> 行均值得到
#   myeloid_score -> 取样本内 top 25% 为髓系 spot。
#   为何用样本内分位: 肝组织本底髓系信号(Kupffer)远高于结肠, 全局阈值会
#   使肝侧几乎全入选; 样本内分位消除组织间与样本间的整体水平差异。
MAC_GENE  <- c("CD68", "AIF1", "C1QA", "C1QB", "CSF1R", "MRC1", "CD163",
               "TYROBP", "FCER1G", "LAPTM5", "CD14", "ITGAM")
# 阴性对照基因集(掩码纯度验证): 上皮与肝细胞, 髓系 spot 应显著更低
EPI_GENE  <- c("EPCAM", "KRT8", "KRT18", "CDH1")
HEP_GENE  <- c("ALB", "APOA1", "APOA2", "TF", "SERPINA1", "TTR")
MAC_Q     <- 0.25   # 样本内髓系评分 top 25% 定义为髓系 spot
# 一次性提取的全部基因: 目标基因 + 髓系标志 + 阴性对照
PANEL_GENE <- c(GENE_SET, MAC_GENE, EPI_GENE, HEP_GENE)
# 4 区室映射(方案第 4 节): 仅保留上皮/肝实质/肿瘤, 基质(Fibroblast/
#   Smooth muscle/Vascular mesenchymal niche)不计入区室, 置 NA
COMP_MAP  <- c("Tumor" = "Tumor", "Normal epithelium" = "Adjacent",
               "Lamina propria" = "Adjacent", "Hepatocytes" = "Adjacent")
# 4 区室展示顺序(与原文 Fig 2 一致)与配色: 结肠=蓝系, 肝=绿系,
#   肿瘤(PT = primary tumor, colorectal)=深色, 癌旁=浅色
COMP_ORDER <- c("Adjacent colon", "PT", "Adjacent liver", "LM")
COMP_COL <- c("Adjacent colon" = "#A6CEE3", "PT" = "#1F78B4",
              "Adjacent liver" = "#B2DF8A", "LM" = "#33A02C")

# =============================================================================
# PART 1: 读取注释表并按样本裁剪至非治疗样本
# 目的: 取得逐 spot 的 region / pred_main / ssGSEA 亚群评分
# 方法: read.csv 后按 SAMPLE_ID 过滤
# 数据: 输入 full_annotation.csv + spot_celltype_annotation.csv
# =============================================================================

cat("========== PART 1: Load annotations ==========\n")

# region 注释来自 03c(结构区, 已含 8 类正式名)
anno.region <- read.csv(file.path(TAB03C, "full_annotation.csv"),
                        stringsAsFactors = FALSE)
# 细胞型注释来自 03(含 treatment 列与 21 个亚群 ssGSEA 评分)
anno.cell <- read.csv(file.path(TAB03, "spot_celltype_annotation.csv"),
                      stringsAsFactors = FALSE)

anno.region <- anno.region[anno.region$sample %in% SAMPLE_ID, ]
anno.cell   <- anno.cell[anno.cell$sample %in% SAMPLE_ID, ]

cat("region 表行数:", nrow(anno.region), "\n")
cat("细胞型表行数:", nrow(anno.cell), "\n")
cat("treatment 取值:", paste(unique(anno.cell$treatment), collapse = ", "), "\n")

# =============================================================================
# PART 2: 逐样本提取 4 基因表达与空间坐标
# 目的: 组装逐 spot 表达 + 坐标 + 注释的主表
# 方法: readRDS 取 Seurat; FetchData(assay="SCT", layer="data") 取归一化表达;
#       GetTissueCoordinates() 取 x/y; FetchData(Spatial/counts) 算检出率
# 数据: 输入 ST-*_seurat.rds; 输出 signature_master.csv
# =============================================================================

cat("========== PART 2: Extract expression and coordinates ==========\n")

spot.list <- list()
det.list  <- list()

for (s in SAMPLE_ID) {
  st <- readRDS(file.path(TAB01, paste0(s, "_seurat.rds")))

  # SCT 归一化表达(log1p 尺度)用于签名与统计: 目标基因 + 髓系/阴性对照基因
  expr.sct <- FetchData(st, vars = PANEL_GENE, assay = "SCT", layer = "data")
  # Spatial 原始 counts 用于检出率
  expr.cnt <- FetchData(st, vars = GENE_SET, assay = "Spatial",
                        layer = "counts")
  # 坐标: 官方访问器, 返回 x/y/cell, 行名即 barcode
  coord <- GetTissueCoordinates(st)

  df <- data.frame(
    sample = s,
    spot   = colnames(st),
    coord[match(colnames(st), rownames(coord)), c("x", "y")],
    expr.sct[match(colnames(st), rownames(expr.sct)), PANEL_GENE],
    stringsAsFactors = FALSE
  )
  names(df)[3:4] <- c("x", "y")
  names(df)[5:(4 + length(PANEL_GENE))] <- paste0("sct_", PANEL_GENE)

  # 检出率: counts > 0 的 spot 占比与平均 counts
  det.list[[s]] <- data.frame(
    sample    = s,
    gene      = GENE_SET,
    n_detect  = colSums(expr.cnt > 0),
    n_spot    = ncol(st),
    pct       = round(100 * colSums(expr.cnt > 0) / ncol(st), 2),
    mean_cnt  = round(colMeans(expr.cnt), 3),
    max_cnt   = as.integer(apply(expr.cnt, 2, max)),
    stringsAsFactors = FALSE
  )
  rownames(det.list[[s]]) <- NULL

  spot.list[[s]] <- df
}

spot.expr <- do.call(rbind, spot.list)
det.tab   <- do.call(rbind, det.list)

cat("逐 spot 表达表行数:", nrow(spot.expr), "\n")
print(det.tab)

# =============================================================================
# PART 3: 合并注释并映射 4 区室
# 目的: 为每 spot 赋予 region / compartment / pred_main / ssGSEA 评分
# 方法: left_join 按 sample+spot; COMP_MAP 把 region 映射为 Tumor/Adjacent,
#       再按 tissue(colon/liver) 组合成 PT/Adjacent colon/LM/Adjacent liver
# 数据: 输出 signature_master.csv + .rds
# =============================================================================

cat("========== PART 3: Merge annotations and map compartments ==========\n")

master <- anno.region[, c("sample", "spot", "cluster", "region")]
master <- merge(master, spot.expr, by = c("sample", "spot"))

# 细胞型表列众多, 仅取主型 + 论文 Fig2 相关的巨噬亚群评分
keep.cell <- c("sample", "spot", "treatment", "pred_main", "pred_score",
               "Myeloid", "SPP1_Macrophages", "MRC1_CCL18_Macrophages",
               "CXCL10_M1_Macrophages", "S100A9_Monocytes")
master <- merge(master, anno.cell[, keep.cell], by = c("sample", "spot"))

# tissue: liver / colon
master$tissue <- ifelse(grepl("liver", master$sample), "liver", "colon")

# 区室: 先由 region 映射 Tumor/Adjacent, 再组合 tissue
#   原发结直肠肿瘤在文中记为 PT(= primary tumor, colorectal)
master$comp_base <- COMP_MAP[master$region]
master$compartment <- ifelse(is.na(master$comp_base), NA,
  ifelse(master$comp_base == "Tumor",
         ifelse(master$tissue == "colon", "PT", "LM"),
         ifelse(master$tissue == "colon", "Adjacent colon",
                "Adjacent liver")))

# --- 数据驱动髓系掩码(方案 v1.1, 弃用有偏的 pred_main) -----------------
# 步骤: 泛巨噬标志基因 SCT 表达 -> 样本内逐基因 z-score -> 行均值 =
#   myeloid_score(消除组织间整体水平差异) -> 样本内前 MAC_Q 定义为髓系
#   另算上皮 epi_score 与肝实质 hep_score 作阴性对照, 验证掩码纯度
ZrowMean <- function(mat) {
  z <- apply(mat, 2, function(v) {
    s <- sd(v)
    if (s == 0) rep(0, length(v)) else (v - mean(v)) / s
  })
  rowMeans(z)
}
mk.cols <- paste0("sct_", MAC_GENE)
ep.cols <- paste0("sct_", EPI_GENE)
hp.cols <- paste0("sct_", HEP_GENE)
master$myeloid_score <- NA_real_
master$epi_score     <- NA_real_
master$hep_score     <- NA_real_
for (s in SAMPLE_ID) {
  idx <- which(master$sample == s)
  master$myeloid_score[idx] <- ZrowMean(as.matrix(master[idx, mk.cols]))
  master$epi_score[idx]     <- ZrowMean(as.matrix(master[idx, ep.cols]))
  master$hep_score[idx]     <- ZrowMean(as.matrix(master[idx, hp.cols]))
}
master$is_macro_pred <- master$pred_main == "Myeloid"   # 旧标签, 仅作对照
master$is_macro <- FALSE
for (s in SAMPLE_ID) {
  idx <- which(master$sample == s)
  cut <- quantile(master$myeloid_score[idx], 1 - MAC_Q,
                  names = FALSE, na.rm = TRUE)
  master$is_macro[idx] <- master$myeloid_score[idx] >= cut
}
# 肿瘤标记: region 为 Tumor 的 spot
master$is_tumor <- master$region == "Tumor"

# 4 个签名列(方案第 4 节): 子集外置 NA, 使 07-13 的空间与统计检验
# 只在相应细胞型/区域上进行, 避免非目标 spot 稀释信号
master$FBXL12_tumor <- ifelse(master$is_tumor, master$sct_FBXL12, NA)
master$PLAU_tumor   <- ifelse(master$is_tumor, master$sct_PLAU, NA)
master$PLAUR_mac    <- ifelse(master$is_macro, master$sct_PLAUR, NA)
master$SPP1_mac     <- ifelse(master$is_macro, master$sct_SPP1, NA)

cat("区室分布:\n")
print(table(master$compartment, useNA = "ifany"))
cat("巨噬(数据驱动掩码) spot 数:", sum(master$is_macro), "\n")
cat("肿瘤 spot 数:", sum(master$is_tumor), "\n")
cat("签名列缺失比例(子集外 NA):\n")
print(round(100 * colMeans(is.na(master[, SIG_SET])), 2))

# 掩码纯度验证: 髓系 spot 应 myeloid 高、上皮/肝实质显著低(阴性对照)
val.rows <- list()
for (s in SAMPLE_ID) {
  d <- master[master$sample == s, ]
  a <- d$myeloid_score[d$is_macro]; b <- d$myeloid_score[!d$is_macro]
  wt <- wilcox.test(a, b, exact = FALSE)
  val.rows[[length(val.rows) + 1]] <- data.frame(
    sample = s,
    n_masked = sum(d$is_macro), n_total = nrow(d),
    myeloid_auc = round(wt$statistic / (length(a) * length(b)), 4),
    med_myeloid_in = round(median(a), 3),
    med_myeloid_out = round(median(b), 3),
    med_epi_in = round(median(d$epi_score[d$is_macro], na.rm = TRUE), 3),
    med_epi_out = round(median(d$epi_score[!d$is_macro], na.rm = TRUE), 3),
    med_hep_in = round(median(d$hep_score[d$is_macro], na.rm = TRUE), 3),
    med_hep_out = round(median(d$hep_score[!d$is_macro], na.rm = TRUE), 3),
    stringsAsFactors = FALSE)
}
val.tab <- do.call(rbind, val.rows)
write.csv(val.tab, file.path(TAB_DIR, "MyeloidMaskValidation.csv"),
          row.names = FALSE)
cat("髓系掩码验证(阴性对照应显著更低, AUC->1):\n")
print(val.tab)

write.csv(master, file.path(TAB_DIR, "signature_master.csv"),
          row.names = FALSE)
saveRDS(master, file.path(TAB_DIR, "signature_master.rds"))
write.csv(det.tab, file.path(TAB_DIR, "GeneDetectionRate.csv"),
          row.names = FALSE)

# =============================================================================
# PART 4: 4 基因检出率可视化(阶段0 门禁的证据)
# 目的: 直观展示 4 基因在各样本的检出率, 判定是否可做后续分析
# 方法: ggplot2 分组柱状图
# 数据: 输出 GeneDetectionRate.pdf
# =============================================================================

cat("========== PART 4: Plot detection rate ==========\n")

det.tab$gene <- factor(det.tab$gene, levels = GENE_SET)
p.det <- ggplot(det.tab, aes(x = gene, y = pct, fill = sample)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 20, linetype = "dashed", linewidth = 0.4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Detection Rate of Target Genes",
       x = "Gene", y = "Detected spots (%)", fill = "Sample") +
  sci_theme

sci_save(p.det, "GeneDetectionRate", width = 8, height = 6)

# =============================================================================
# PART 5: 数据表征可视化(SCI 发表标准)
# 目的: 展示免疫格局、目标基因与签名的区室分布, 以及数据驱动髓系掩码纯度
# 方法:
#   (1) 细胞型组成: 堆叠比例柱(pred_main, 展示免疫格局, 对应原文 Fig 1G)
#   (2) 4 基因表达分布: 箱线图, 按基因分面、按区室着色
#   (3) 4 签名分布: 箱线图, 验证子集限定签名在各区室间确有差异
#   (4) 髓系掩码纯度: 箱线图对比 髓系/上皮/肝实质 评分在掩码内外, 阴性
#       对照应显著更低, 证明掩码区分真实髓系而非随机
# 规范17: 图件标题与标签全英文; 规范16: 文件名禁 Fig 前缀
# =============================================================================

cat("========== PART 5: Data characterization plots ==========\n")

# (1) 细胞型组成: 03 label transfer 的免疫主型格局
ct.tab <- as.data.frame(table(master$sample, master$pred_main))
names(ct.tab) <- c("sample", "celltype", "n")
p.ct <- ggplot(ct.tab, aes(x = sample, y = n, fill = celltype)) +
  geom_col(position = "fill", width = 0.72) +
  labs(title = "Main Cell Type Composition per Sample", x = "Sample",
       y = "Proportion of spots", fill = "Cell type") +
  sci_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9))
sci_save(p.ct, "CellTypeComposition", width = 8, height = 6)

# (2) 4 基因表达分布: 箱线图(真实表达, 仅检出 spot)
#   未检出的 spot 在 SCT 归一化下记为 0, 会压扁箱体并使真实表达变成离群点;
#   故仅取表达 > 0 的 spot 绘制箱线, 并隐藏离群点, 呈现真实表达分布
gene.long <- melt(master[, c("compartment", paste0("sct_", GENE_SET))],
                  id.vars = "compartment", variable.name = "gene",
                  value.name = "expr")
gene.long$gene <- sub("^sct_", "", gene.long$gene)
gene.long <- gene.long[!is.na(gene.long$compartment) & gene.long$expr > 0, ]
gene.long$compartment <- factor(gene.long$compartment,
                                levels = COMP_ORDER)
# 各基因×区室中位数标注; 分布用 violin 展示(右偏/多 ties 数据下
#   箱线图 IQR 退化为线, violin 直接呈现密度形状);
#   抖动散点揭示真实点密度与 ties 行
#   (SCT data 层为 log1p 归一化, 右偏分布导致部分区室箱体扁至低位)
gene.med <- aggregate(expr ~ gene + compartment, data = gene.long,
                      FUN = median, na.rm = TRUE)
gene.med$label <- sprintf("%.2f", gene.med$expr)

p.gd <- ggplot(gene.long, aes(x = compartment, y = expr,
                              fill = compartment)) +
  geom_violin(width = 0.8, alpha = 0.35, trim = TRUE,
              scale = "width", color = "#666666", linewidth = 0.3) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.12, size = 0.35,
              color = "#333333") +
  geom_text(data = gene.med, aes(x = compartment, y = expr,
                                 label = label, fill = NULL),
            vjust = -0.8, size = 3.2, color = "#1A1A1A",
            show.legend = FALSE) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = COMP_COL) +
  labs(title = "Target Gene Expression by Compartment",
       subtitle = paste0("violin = density; dots = individual Visium spots; ",
                         "label = median (detected spots only)"),
       x = "Compartment", y = "SCT normalized expression") +
  sci_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
        legend.position = "none")
sci_save(p.gd, "GeneExpressionDistribution", width = 11, height = 8)

# (3) 4 签名分布: 箱线图(真实表达, 仅检出 spot), 验证子集签名在区室间差异
sig.long <- melt(master[, c("compartment", SIG_SET)],
                 id.vars = "compartment", variable.name = "signature",
                 value.name = "value")
sig.long <- sig.long[!is.na(sig.long$compartment) & !is.na(sig.long$value) &
                       sig.long$value > 0, ]
sig.long$compartment <- factor(sig.long$compartment,
                               levels = COMP_ORDER)
lab.map <- c(FBXL12_tumor = "FBXL12 tumor", PLAU_tumor = "PLAU tumor",
             PLAUR_mac = "PLAUR macrophage", SPP1_mac = "SPP1 macrophage")
sig.long$signature <- lab.map[as.character(sig.long$signature)]
# 各签名×区室中位数标注; violin 展示密度(签名多 ties 时箱线 IQR 退化为线);
#   抖动散点揭示 ties 行
sig.med <- aggregate(value ~ signature + compartment, data = sig.long,
                     FUN = median, na.rm = TRUE)
sig.med$label <- sprintf("%.2f", sig.med$value)

p.sd <- ggplot(sig.long, aes(x = compartment, y = value,
                             fill = compartment)) +
  geom_violin(width = 0.8, alpha = 0.35, trim = TRUE,
              scale = "width", color = "#666666", linewidth = 0.3) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.12, size = 0.35,
              color = "#333333") +
  geom_text(data = sig.med, aes(x = compartment, y = value,
                                label = label, fill = NULL),
            vjust = -0.8, size = 3.2, color = "#1A1A1A",
            show.legend = FALSE) +
  facet_wrap(~ signature, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = COMP_COL) +
  labs(title = "Signature Score by Compartment",
       subtitle = paste0("violin = density; dots = individual Visium spots; ",
                         "label = median (detected spots only)"),
       x = "Compartment", y = "Signature score (SCT normalized)") +
  sci_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
        legend.position = "none")
sci_save(p.sd, "SignatureDistribution", width = 11, height = 8)

# (4) 髓系掩码纯度: 三套评分(髓系/上皮/肝实质)在掩码内外的箱线对比
#   阴性对照(上皮、肝实质)在掩码内应显著低于掩码外, 证明掩码特异
score.long <- data.frame(
  score = rep(c("myeloid", "epithelial", "hepatocyte"), each = nrow(master)),
  value = c(master$myeloid_score, master$epi_score, master$hep_score),
  mask  = rep(ifelse(master$is_macro, "myeloid mask", "other"), 3),
  compartment = rep(master$compartment, 3),
  stringsAsFactors = FALSE)
score.long <- score.long[!is.na(score.long$value), ]
p.mv <- ggplot(score.long, aes(x = mask, y = value, fill = mask)) +
  geom_boxplot(width = 0.5, outlier.size = 0.2, alpha = 0.85) +
  facet_wrap(~ score, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("myeloid mask" = "#1F78B4",
                               "other" = "#BEBADA")) +
  labs(title = "Myeloid Mask Purity Validation",
       x = "Spot group", y = "Z-scored signature score") +
  sci_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
        legend.position = "none")
sci_save(p.mv, "MyeloidMaskValidation", width = 11, height = 4.5)

cat("06_PLAUR_SPP1_signature.R done: signature_master.csv + ",
    "GeneDetectionRate.pdf + 4 characterization plots\n", sep = "")
