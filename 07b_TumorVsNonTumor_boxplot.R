# ===============================================================
# 07_TumorVsNonTumor_boxplot.R  (本研究自有分析: 肿瘤 vs 非肿瘤区域表达量箱线比较)
# 功能: 对 4 个轴基因(SPP1/PLAUR/PLAU/FBXL12), 在 4 个样本(P1/P2 colon 与 liver)内,
#       比较 is_tumor=TRUE (肿瘤区域) 与 is_tumor=FALSE (非肿瘤区域) 的 SCT 表达量,
#       画小提琴图并标注 Wilcoxon 显著性(* / ** / *** / ns).
# 定位: 与 07 空间热力图配套, 展示"基因在肿瘤 vs 非肿瘤区域表达差异"; 属本研究自有分析,
#       仅参照原文可视化风格. 不受原文 Fig 2C 仅 LM 限制.
# 可视化选择: 用户反馈低丰度基因的箱线图被压缩在 0 附近; 故改用小提琴图(ggpubr::ggviolin),
#             用密度形状呈现 0 值堆积 + 少数高表达 spot 的分布, 叠加 median+IQR 标记.
# 数据: signature_master.rds (is_tumor, sct_* 基因表达列). 表达量 = SCT 归一化表达(与 07 一致).
# ===============================================================

# --- 包加载 (均为已安装标准包) -------------------------------------
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))     # ggboxplot + stat_compare_means
suppressMessages(library(rstatix))    # wilcox_test 统计表

# --- 目录结构 -----------------------------------------------------------
BASE_DIR <- "X:/XX"
TAB06    <- file.path(BASE_DIR, "results/06_PLAUR_SPP1_signature/table")
FIG_DIR  <- file.path(BASE_DIR, "results/07_PLAUR_SPP1_spatial/Figure")
TAB_DIR  <- file.path(BASE_DIR, "results/07_PLAUR_SPP1_spatial/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 读取 ---------------------------------------------------------------
master <- readRDS(file.path(TAB06, "signature_master.rds"))
GENES  <- c("SPP1", "PLAUR", "PLAU", "FBXL12")
cat("master rows:", nrow(master), " samples:", paste(unique(master$sample), collapse = ","), "\n")

# 肿瘤/非肿瘤分组 (因子, 显式水平保证 NonTumor 在前, 配对比较方向一致)
master$tumor_group <- factor(ifelse(master$is_tumor, "Tumor", "NonTumor"),
                              levels = c("NonTumor", "Tumor"))

# 蓝-红对比: 非肿瘤=蓝, 肿瘤=红 (用户指定非肿瘤蓝色)
PAL <- c(NonTumor = "#4575B4", Tumor = "#E31A1C")

# --- 拼接 long 表 (gene, sample, tumor_group, expr) 供 rstatix 逐组检验 ----
long <- do.call(rbind, lapply(GENES, function(g) {
  col <- paste0("sct_", g)
  data.frame(sample      = master$sample,
             gene        = g,
             tumor_group = master$tumor_group,
             expr        = master[[col]])
}))
long <- long[!is.na(long$expr), ]

# 统计表: 每 (sample x gene) 内 Wilcoxon 检验 Tumor vs NonTumor
# rstatix 风格: 用 dplyr::group_by 分组建模(此版本 wilcox_test 不接受 group.by=)
suppressMessages(library(dplyr))
stat_tbl <- long %>% dplyr::group_by(sample, gene) %>%
  rstatix::wilcox_test(expr ~ tumor_group)
# Holm 校正 p (rstatix 标准函数, 非自造)
stat_tbl <- rstatix::adjust_pvalue(stat_tbl, p.col = "p", method = "holm")
stat_tbl <- stat_tbl[, c("sample", "gene", "group1", "group2",
                         "n1", "n2", "statistic", "p", "p.adj")]
write.csv(stat_tbl, file.path(TAB_DIR, "TumorVsNonTumor_boxplot_stats.csv"),
          row.names = FALSE)
cat("wilcox_test table rows:", nrow(stat_tbl), "\n")

# --- 逐基因画小提琴图 (4 个样本分面, 每面 2 提琴 + 显著性括号) -----
for (g in GENES) {
  col  <- paste0("sct_", g)
  df   <- long[long$gene == g, ]
  p <- ggviolin(df, x = "tumor_group", y = "expr",
                fill = "tumor_group", palette = PAL,
                facet.by = "sample", short.panel.labs = TRUE, nrow = 1,
                scales = "free_y",          # y 轴按各样本自适应, 避免 colon 被 liver 压扁
                add = "median_iqr",         # 显示中位点+IQR 线, 不画抖动黑点
                trim = TRUE, scale = "width",
                ggtheme = theme_bw(base_size = 12),
                ylab = paste0(g, " SCT expression"), xlab = "Region") +
    stat_compare_means(method = "wilcox.test",
                       comparisons = list(c("NonTumor", "Tumor")),
                       label = "p.signif",
                       bracket.size = 0.8, text.size = 4,
                       hide.ns = FALSE) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.18))) +  # y 起点=0 (用户指令); 顶部留 18% 余量, 确保星号/括号不被截断
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "top",
          strip.text = element_text(size = 11, face = "plain"),
          plot.title = element_text(size = 13, hjust = 0.5, face = "plain")) +
    labs(title = paste0(g, ": Tumor vs Non-tumor region (SCT expression)"),
         fill = "Region")
  pdf_path <- file.path(FIG_DIR, paste0("Violin_TumorVsNonTumor_", g, ".pdf"))
  ggsave(pdf_path, p, width = 9, height = 4.5, dpi = 150, family = "Helvetica")
  cat("saved", basename(pdf_path), "\n")
}

cat("07 violin done: 4 PDFs (1 per gene) + TumorVsNonTumor_boxplot_stats.csv\n")
