# ===============================================================
# 08_PLAUR_SPP1_heatmap.R
# 功能: 绘制双层区室比较变量(6) + 原文 anchor(2) x 4 区室的缩放热图(Fig 2A)
# 定位: 方案 Fig 2A(复刻原文 Fig 2A 的 scaled proportion 热图), 方案 v1.1
# 可视化: CompartmentSignatureHeatmap.pdf
# 数据: 输入 results/08_PLAUR_SPP1_heatmap/table/signature_master.rds
#       输出 results/08_PLAUR_SPP1_heatmap/table/CompartmentMeans.csv
# 依据: 方案第 7.1 节 Fig 2A; 原文 Fig 2A 行方向缩放后按富集模式聚类
# ===============================================================

# --- 包加载(统计 -> 数据 -> 可视化 -> 领域) -----------------------------
suppressMessages({
  library(pheatmap)      # 热图绘制
  library(RColorBrewer)  # 配色
  library(reshape2)      # 长宽格式转换
  library(ggplot2)       # 可视化
})

# --- 全局辅助函数(内联, 规范22 自给自足) -------------------------------
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
TAB06    <- file.path(BASE_DIR, "results/06_PLAUR_SPP1_signature/table")
FIG_DIR  <- file.path(BASE_DIR, "results/08_PLAUR_SPP1_heatmap/Figure")
TAB_DIR  <- file.path(BASE_DIR, "results/08_PLAUR_SPP1_heatmap/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 分析常量 -----------------------------------------------------------
# 4 区室顺序与原文 Fig 2 一致; 原发结直肠肿瘤记为 PT (primary tumor)
COMP_ORDER <- c("Adjacent colon", "PT", "Adjacent liver", "LM")
# 双层区室比较变量集(方案 v1.1)
# 理由: FBXL12_tumor / PLAU_tumor 为"肿瘤区内"子集限定签名, 在 Adjacent
#       colon / Adjacent liver 的 spot 数 = 0(癌旁按定义无 Tumor 区), 用于
#       区室间比较会产生无意义空列; 故区室比较改用双层设计:
#       全组织层 = 全部 spot 的原始 SCT 表达(4 区室均有值, 含基质/免疫背景)
#       细胞型层 = 仅 Myeloid spot 的每巨噬表达(4 区室均有值, 去细胞组成混杂)
GENE_RAW  <- c("sct_FBXL12", "sct_PLAU", "sct_PLAUR", "sct_SPP1")
MAC_SIG   <- c("PLAUR_mac", "SPP1_mac")
CMP_SET   <- c(GENE_RAW, MAC_SIG)
# 验证性 anchor: SPP1_Macrophages ssGSEA(外部已发表 SPP1+ 巨噬基因集),
# 独立交叉验证本分析签名 SPP1_mac; MRC1_CCL18_Macrophages 已按分析原则移除
ANCHOR_SET <- c("SPP1_Macrophages")
# 英文行标签(规范17)
LAB_MAP <- c(
  "sct_FBXL12" = "FBXL12 (all spots)",
  "sct_PLAU"   = "PLAU (all spots)",
  "sct_PLAUR"  = "PLAUR (all spots)",
  "sct_SPP1"   = "SPP1 (all spots)",
  "PLAUR_mac"  = "PLAUR (per macrophage)",
  "SPP1_mac"   = "SPP1 (per macrophage)"
)

# =============================================================================
# PART 1: 计算各区室内各签名的均值(spot 级, 方案 v1.1)
# 目的: 得到 签名 x 区室 的均值矩阵, 作为热图与后续趋势图的基础
# 方法: 直接对区内全部 spot 取均值(不再做"先样本后患者"两级平均, 以免
#       极小单元被等权放大导致结论方向颠倒); 同步披露每个单元 spot 数,
#       对 n<30 的单元标注 low_n 提示解释需谨慎
# 数据: 输出 CompartmentMeans.csv + CompartmentSpotCounts.csv(n 披露)
# =============================================================================

cat("========== PART 1: Compartment means (spot-level) ==========\n")

master <- readRDS(file.path(TAB06, "signature_master.rds"))
master <- master[!is.na(master$compartment), ]

all.set <- c(CMP_SET, ANCHOR_SET)

# 主估计: spot 级均值
mat <- t(sapply(all.set, function(v)
  sapply(COMP_ORDER, function(cp)
    mean(master[[v]][master$compartment == cp], na.rm = TRUE))))
colnames(mat) <- COMP_ORDER

# 单元 spot 数披露
n.sig <- expand.grid(signature = all.set, compartment = COMP_ORDER,
                     stringsAsFactors = FALSE)
n.sig$n <- sapply(seq_len(nrow(n.sig)), function(i)
  sum(!is.na(master[[n.sig$signature[i]]][
        master$compartment == n.sig$compartment[i]])))
n.sig$low_n <- n.sig$n < 30

cat("签名 x 区室 均值矩阵(spot 级):\n")
print(round(mat, 4))
cat("单元 spot 数(n<30 标注 low_n):\n")
print(n.sig)

write.csv(data.frame(signature = rownames(mat), mat, check.names = FALSE),
          file.path(TAB_DIR, "CompartmentMeans.csv"), row.names = FALSE)
write.csv(n.sig, file.path(TAB_DIR, "CompartmentSpotCounts.csv"),
          row.names = FALSE)

# =============================================================================
# PART 2: 行方向 Z-score 缩放并绘制热图(Fig 2A)
# 目的: 复刻原文 Fig 2A, 直观比较各签名在 4 区室的相对富集
# 方法: 行内 scale() 标准化; pheatmap 绘制, 蓝-白-红配色
# 数据: 输出 CompartmentSignatureHeatmap.pdf
# =============================================================================

cat("========== PART 2: Fig 2A heatmap ==========\n")

mat.z <- t(scale(t(mat)))
colnames(mat.z) <- COMP_ORDER

# 行名改为可读英文标签(规范17: 全英文, 禁省略号)
rownames(mat.z) <- c(unname(LAB_MAP[CMP_SET]), "SPP1 Macrophages ssGSEA")

p.hm <- pheatmap(mat.z,
                 color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
                 breaks = seq(-2, 2, length.out = 101),
                 cluster_rows = TRUE, cluster_cols = FALSE,
                 border_color = NA,
                 main = "Scaled Signature Score by Compartment",
                 fontsize_row = 11, fontsize_col = 12,
                 cellwidth = 55, cellheight = 26)

pdf(file.path(FIG_DIR, "CompartmentSignatureHeatmap.pdf"), width = 11,
    height = 8, family = "Helvetica")
print(p.hm)
dev.off()

# =============================================================================
# PART 3: 区室绝对均值点图(Fig 2A 配套)
# 目的: 热图只展示行内 Z-score 相对模式, 会隐藏绝对表达量; 用点图补充各
#       签名在 4 区室的绝对均值, 使读者同时看到"趋势"与"强度"
# 方法: 8 签名 x 4 区室均值矩阵长宽转换; 每签名一个 facet(scales="free_y")
# 数据: CompartmentSignatureDotplot.pdf
# =============================================================================

cat("========== PART 3: Compartment absolute mean dot plot ==========\n")

# 长格式, 保留热图行序
row.labs <- c(unname(LAB_MAP[CMP_SET]), "SPP1 Macrophages ssGSEA")
df.dot <- reshape2::melt(mat, value.name = "mean_score")
colnames(df.dot)[1:2] <- c("signature", "compartment")
df.dot$signature <- factor(df.dot$signature, levels = rownames(mat),
                           labels = row.labs)
df.dot$compartment <- factor(df.dot$compartment, levels = COMP_ORDER)

p.dot <- ggplot(df.dot, aes(x = compartment, y = mean_score)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = mean_score, size = mean_score), shape = 16) +
  scale_colour_gradientn(
    colours = c("#FFFFCC", "#FED976", "#FEB24C", "#FD8D3C",
                "#FC4E2A", "#E31A1C", "#BD0026", "#800026")
  ) +
  scale_size_continuous(range = c(2, 6)) +
  facet_wrap(~ signature, scales = "free_y", ncol = 4) +
  labs(title = "Mean signature score by compartment (unscaled)",
       x = "Compartment", y = "Mean SCT score") +
  sci_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
        strip.text = element_text(size = 9, face = "plain"))

pdf(file.path(FIG_DIR, "CompartmentSignatureDotplot.pdf"), width = 11,
    height = 7, family = "Helvetica")
print(p.dot)
dev.off()

# =============================================================================
# PART 4: 区室 spot 数条形图(样本量披露)
# 目的: SCI 发表必须展示每个区室的 spot 样本量; 热图/点图均不直接展示 n
# 方法: 按 compartment 统计总 spot 数, 低 n(<30)用颜色警示
# 数据: CompartmentSpotCounts.pdf
# =============================================================================

cat("========== PART 4: Compartment spot counts bar plot ==========\n")

cnt <- as.data.frame(table(master$compartment), stringsAsFactors = FALSE)
colnames(cnt) <- c("compartment", "n")
cnt$compartment <- factor(cnt$compartment, levels = COMP_ORDER)
cnt$low_n <- cnt$n < 30
cnt$fill_col <- ifelse(cnt$low_n, "#E31A1C", "#1F78B4")

p.cnt <- ggplot(cnt, aes(x = compartment, y = n, fill = compartment)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = n), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Adjacent colon" = "#1F78B4",
                               "PT" = "#E31A1C",
                               "Adjacent liver" = "#A6CEE3",
                               "LM" = "#FD8D3C")) +
  labs(title = "Spot counts per compartment",
       x = "Compartment", y = "Number of spots") +
  sci_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
        legend.position = "none")

pdf(file.path(FIG_DIR, "CompartmentSpotCounts.pdf"), width = 6,
    height = 5, family = "Helvetica")
print(p.cnt)
dev.off()

cat("08_PLAUR_SPP1_heatmap.R done:\n")
cat("  - CompartmentSignatureHeatmap.pdf (enlarged)\n")
cat("  - CompartmentSignatureDotplot.pdf\n")
cat("  - CompartmentSpotCounts.pdf\n")
