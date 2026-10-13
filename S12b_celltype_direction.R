# ===============================================================================================
# 功能：细胞层（解卷积）方向代理 —— 每细胞谱系比例 x 每五段段分 的样本级相关
# 定位：S12 补充分析。segment_by_celltype.csv 无 per-sample 方向信息，故以 S09 的
#       逐样本细胞谱系比例（tool_celltype_fractions_<acc>.csv 三工具共识）与 S06 的
#       逐样本段分（segment_score_matrix.rds）做相关，符号定方向，口径同 S09 partial-correlation。
# 可视化：谱系 x 五段 的方向代理热图。
# 方法依据：样本级 Pearson 相关（cor.test），BH 校正；n < 3 判 indeterminate（同 S09 阈值）。
# ===============================================================================================

library(ggplot2)  # 热图可视化

## ---- 常量----
SEED <- 20260925
S09_DIR <- "03_results/S09/table"
S12_DIR <- "03_results/S12/table"
FIG_DIR <- "03_results/S12/Figure"
MIN_N <- 3L   # 相关所需最小样本数（同 S09 partial-correlation 阈值）
vec.segment <- c("S1_initiation", "S2_amplification", "S3_effect",
  "S4_resolution", "S5_chronic")

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

## 主题/出口：规范 §sci_theme / §sci_save 官方定义（与 S08 逐字一致，禁 source）
sci_theme <- ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
    legend.position = "right",
    legend.text = ggplot2::element_text(size = 10, face = "plain"),
    axis.title = ggplot2::element_text(size = 12, face = "plain"),
    axis.text = ggplot2::element_text(size = 10, face = "plain"),
    plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "plain"),
    legend.title = ggplot2::element_text(size = 11, face = "plain"))
sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  pdf(pdf_file, width = width, height = height)
  print(plot_obj)
  dev.off()
}

## ---- 发现含解卷积比例的系列 ----
vec.frac.file <- list.files(S09_DIR, pattern = "^tool_celltype_fractions_.*\\.csv$",
  full.names = TRUE)
stopifnot(length(vec.frac.file) > 0L)
vec.acc <- sub("^tool_celltype_fractions_", "",
  sub("\\.csv$", "", basename(vec.frac.file)))

dat.dir <- data.frame()
for (acc in vec.acc) {
  f.frac <- file.path(S09_DIR, paste0("tool_celltype_fractions_", acc, ".csv"))
  dat.frac <- read.csv(f.frac, stringsAsFactors = FALSE)
  ## 三工具共识：按 (lineage, sample) 取比例均值
  ## （aggregate 返回命名 data.frame；旧 rowsum/除法写法会丢失向量名致 strsplit 失败）
  dat.frac.ag <- aggregate(fraction ~ lineage + sample, data = dat.frac,
    FUN = mean, stringsAsFactors = FALSE)
  ## 段分（segment x sample）
  f.score <- "03_results/S06/table/segment_score_matrix.rds"
  stopifnot(file.exists(f.score))
  mat.score <- as.matrix(readRDS(f.score)[[acc]])
  vec.seg.use <- intersect(vec.segment, rownames(mat.score))
  stopifnot(length(vec.seg.use) > 0L)
  mat.score <- mat.score[vec.seg.use, , drop = FALSE]
  ## 按样本名对齐
  dat.frac.ag <- dat.frac.ag[dat.frac.ag$sample %in% colnames(mat.score), ]
  stopifnot(nrow(dat.frac.ag) > 0L)
  ## 逐 (lineage, segment) 相关
  for (lin in unique(dat.frac.ag$lineage)) {
    vec.sam <- dat.frac.ag$sample[dat.frac.ag$lineage == lin]
    v.frac <- dat.frac.ag$fraction[dat.frac.ag$lineage == lin]
    names(v.frac) <- vec.sam
    for (seg in vec.seg.use) {
      v.score <- as.numeric(mat.score[seg, colnames(mat.score) %in% vec.sam])
      names(v.score) <- colnames(mat.score)[colnames(mat.score) %in% vec.sam]
      vec.com <- intersect(names(v.frac), names(v.score))
      n <- length(vec.com)
      if (n < MIN_N) {
        dat.dir <- rbind(dat.dir, data.frame(acc = acc, lineage = lin,
          segment = seg, n_sample = n, corr = NA_real_, p = NA_real_,
          p_adj = NA_real_, direction = "indeterminate", stringsAsFactors = FALSE))
        next
      }
      r <- cor.test(v.frac[vec.com], v.score[vec.com])
      dat.dir <- rbind(dat.dir, data.frame(acc = acc, lineage = lin,
        segment = seg, n_sample = n, corr = as.numeric(r$estimate),
        p = as.numeric(r$p.value), p_adj = NA_real_,
        direction = "pending", stringsAsFactors = FALSE))
    }
  }
  cat(sprintf("[S12b] %s: %d lineage x %d segment pairs evaluated\n", acc,
    length(unique(dat.frac.ag$lineage)), length(vec.seg.use)))
}
## BH 校正（按系列内检验族）
dat.dir$p_adj <- unlist(by(dat.dir, dat.dir$acc, function(d) {
  padj <- d$p
  padj[!is.na(padj)] <- p.adjust(d$p[!is.na(d$p)], method = "BH")
  padj
}))
dat.dir$direction <- ifelse(!is.na(dat.dir$p_adj) & dat.dir$p_adj < 0.05,
  ifelse(dat.dir$corr > 0, "up", "down"), "indeterminate")
write.csv(dat.dir, file.path(S12_DIR, "celltype_segment_direction.csv"),
  row.names = FALSE)

## ---- 方向代理热图（谱系 x 五段，按符号着色；indeterminate 灰）----
dat.dir$seg_lab <- factor(dat.dir$segment, levels = vec.segment)
dat.dir$dir_sign <- ifelse(dat.dir$direction == "up", 1,
  ifelse(dat.dir$direction == "down", -1, 0))
sci_save(ggplot(dat.dir, aes(x = seg_lab, y = lineage, fill = factor(dir_sign))) +
  geom_tile(colour = "grey80", size = 0.2) +
  scale_fill_manual(name = "direction",
    values = c("-1" = "blue", "0" = "grey90", "1" = "red"),
    labels = c("down", "indeterminate", "up")) +
  facet_wrap(~ acc, scales = "free_y") +
  labs(title = "Cell-layer directional proxy (lineage x segment)",
    x = "Segment", y = "Lineage") +
  sci_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  "Celltype_Segment_Direction_Proxy", width = 10, height = 6)

cat(sprintf("\n[S12b] done: %d lineage-segment pairs -> celltype_segment_direction.csv\n",
  nrow(dat.dir)))
