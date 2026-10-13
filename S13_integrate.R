#===============================================================
## S13_integrate.R —— 整合与对账（步骤规格 §8：S13）
##
## 目的：把层 1（12 基因面板）与层 2（段分）、层 3（跨集方向）、细胞层、空间层
##       按「段 x 数据集 x 细胞类型/区域 x 方向」聚合，产出对账矩阵，
##       并**强制披露 indeterminate 占比**（§4.2 / §10 预注册义务）。
##
## 定位：分析脚本。本步只做聚合与对账，不引入新统计。
##
## 产出：03_results/S13/table/utilization_matrix.csv
##       03_results/S13/table/indeterminate_report.csv
##       03_results/S13/Figure/Segment_Direction_By_Series_Heatmap.pdf 分面：系列 x 层方向
##       03_results/S13/Figure/Indeterminate_Fraction_By_Readout.pdf   分面：读出层 x 占比
##       03_results/S13/Figure/Segment_Direction_Matrix.pdf            段 x 数据集 方向矩阵
##
## 图件清单：
##   可视化: S13 共 3 类
##     1. Segment_Direction_By_Series_Heatmap.pdf
##        —— 形式：分面瓦片（x = 读出层，y = 段，按系列分面，fill = 方向），
##           每格叠加真实效应量文字标注（layer1/2 = signed log2FC 取自 S05/S06；
##           layer3 = pooled SMD 取自 S07；单位不同仅作文字、不混用连续色阶）。
##        —— 科学问题：每一段在三个读出层上的方向判定是否自洽 + 定量幅值语境。
##        —— 参照文献：跨研究方向/竞争性检验口径（Wu & Smyth 2012, Nucleic Acids Res,
##           doi:10.1093/nar/gks461, PMID 22684630）。
##     2. Indeterminate_Fraction_By_Readout.pdf
##        —— 形式：单色柱（x = 读出层，y = indeterminate 占比，柱顶写 coverage + n_total），
##           加 0.5 预注册披露阈值虚线（§4.2）。
##        —— 科学问题：各读出层的 indeterminate 占比（§4.2 硬性披露义务）+ 样本量语境。
##        —— 参照文献：不确定性披露口径见方案 §4.2 / §10；异质性度量参见
##           Higgins 2003, BMJ, doi:10.1136/bmj.327.7414.557, PMID 12958120。
##     3. Segment_Direction_Matrix.pdf
##        —— 形式：单面板瓦片（x = 系列，y = 段，fill = 层 2 方向）
##        —— 科学问题：层 2 段分方向在「段 x 系列」上的对账矩阵。
##        —— 参照文献：同上 Wu & Smyth 2012, PMID 22684630。
##
## 图件口径：sci_theme / sci_save 逐字取自规范 §sci_theme / §sci_save 官方定义；
##       分面一律 facet_wrap / facet_grid，不引入任何拼图包；可读性核验见
##       sci_save 之 RAPA_QC_PNG 钩子。
##
## 前置：S05 / S06 / S07 —— 缺任一层即显式 stop，不挑食、不择优（§4.2）。
##       S02（segment_overlap.csv）缺失仅 WARN，视为段间无重叠（joint_weight = 1）。
##       S12（细胞/空间层）为**可选补充**：仅当其表含 direction 列时才并入
##       indeterminate 披露，否则告警并视为缺席——S13 不得因 S12 而阻塞
##       （panel 级 QC 表无方向语义，原 stop 属契约错位）。
#===============================================================

## PART 0 常量 ------------------------------------------------------------------------------

SEED <- 20260925
set.seed(SEED)
FIG_DIR <- "03_results/S13/Figure"
TAB_DIR <- "03_results/S13/table"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

vec.segment <- c("S1_initiation", "S2_amplification", "S3_effect",
  "S4_resolution", "S5_chronic")

library(ggplot2)  ## 对账图；主题与出口用规范 §sci_theme / §sci_save

## 主题：规范 §sci_theme 官方定义（theme_bw + 关网格），全脚本统一口径
sci_theme <- ggplot2::theme_bw(base_size = 12) +
ggplot2::theme(
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank(),
  strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
  legend.position  = "right",
  legend.text      = ggplot2::element_text(size = 10, face = "plain"),
  axis.title       = ggplot2::element_text(size = 12, face = "plain"),
  axis.text        = ggplot2::element_text(size = 10, face = "plain"),
  plot.title       = ggplot2::element_text(size = 14, hjust = 0.5,
    face = "plain"),
  legend.title     = ggplot2::element_text(size = 11, face = "plain")
)

## 出口：规范 §sci_save 官方定义，仅 PDF，落本步 FIG_DIR
sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  pdf(pdf_file, width = width, height = height)
  print(plot_obj)
  dev.off()
  ## 可读性核验钩子：默认关闭，交付物仍为纯 PDF；
  ## 置 RAPA_QC_PNG=1 时额外导 _qc/ PNG 供逐张目视核验，不改变 PDF 产出。
  if (Sys.getenv("RAPA_QC_PNG") == "1") {
    dir.create(file.path(FIG_DIR, "_qc"), showWarnings = FALSE)
    png(file.path(FIG_DIR, "_qc", paste0(filename, ".png")),
      width = width * 150, height = height * 150, res = 150)
    print(plot_obj)
    dev.off()
  }
}

## PART 1 汇总三层结果 --------------------------------------------------------------------------

vec.need <- c("03_results/S05/table/layer1_gene_dir.csv",
  "03_results/S06/table/layer2_segment_all.csv",
  "03_results/S07/table/layer3_meta_dir.csv")
vec.missing <- vec.need[!file.exists(vec.need)]
if (length(vec.missing) > 0L) {
  stop(paste0("S13 STOP: missing upstream table(s): ", paste(vec.missing,
        collapse = "; ")))
}

## §4.2 段间重叠折扣：Jaccard > 0.5 的参与段，其「联合贡献」x 0.5。
## 判定证据由 S02 在 segment_overlap.csv 给出，本步只读判定值，不重算 Jaccard。
file.ov <- "03_results/S02/table/segment_overlap.csv"
vec.seg.disc <- if (!file.exists(file.ov)) {
  message("[WARN] ", file.ov,
    " absent — treating every segment as non-overlapping ",
    "(joint_weight = 1). Run S02 first for the §4.2 discount.")
  character(0)
} else {
  dat.ov <- read.csv(file.ov, stringsAsFactors = FALSE)
  if (!all(c("row.segment", "col.segment",
        "overlap_discounted") %in% names(dat.ov))) {
    stop(paste0("S13 STOP: segment_overlap.csv lacks one of row.segment",
      "/ col.segment / overlap_discounted"))
  }
  unique(dat.ov$row.segment[dat.ov$overlap_discounted])
}
vec.joint.weight <- ifelse(vec.segment %in% vec.seg.disc, 0.5, 1.0)

dat.l1 <- read.csv(vec.need[1L], stringsAsFactors = FALSE)
dat.l2 <- read.csv(vec.need[2L], stringsAsFactors = FALSE)
dat.l3 <- read.csv(vec.need[3L], stringsAsFactors = FALSE)

## 层 1 自 2026-09-29 起是「系列 x 对比 x 基因」长表（GSE57691 有 3 个对比、
## GSE323983 有 4 个），S13 的方向投票必须限定在**主对比**上，否则同系列多对比会把
## IL6/IL1B/TNF 的票数成倍放大。主对比由 S05 的 dat.contrasts$primary 声明。
if (!all(c("primary", "exploratory", "contrast") %in% colnames(dat.l1))) {
  stop(paste0(
      "S13 STOP: layer1_gene_dir.csv 缺 primary / exploratory / contrast 列；",
      "该表由 S05 产出，请先重跑 S05。"))
}
dat.l1.primary <- dat.l1[dat.l1$primary, , drop = FALSE]

## 对账矩阵：段 x 数据集 的方向（层 2）+ 是否显著（层 1）
dat.util <- data.frame(segment = character(0), acc = character(0),
  layer1_dir = character(0),
  layer1_sig = character(0), layer2_dir = character(0),
  layer3_verdict = character(0), exploratory = logical(0),
  joint_weight = numeric(0), stringsAsFactors = FALSE)

for (acc in unique(dat.l2$acc)) {
  for (seg in vec.segment) {
    a <- dat.l1.primary[dat.l1.primary$acc == acc &
      dat.l1.primary$gene %in% c("IL6", "IL1B", "TNF"), ]
    dir.l1 <- if (nrow(a) == 0L) "n/a" else
    ifelse(all(a$direction == "indeterminate"), "indeterminate",
      ifelse(sum(a$direction == "up") > sum(a$direction == "down"), "up",
        "down"))
    sig.l1 <- if (nrow(a) == 0L) "n/a" else
    ifelse(sum(a$significant == "YES") >= 2L, "YES", "NO")
    b <- dat.l2[dat.l2$acc == acc & dat.l2$segment == seg, ]
    dir.l2 <- if (nrow(b) == 0L) "n/a" else as.character(b$direction[1L])
    c3 <- dat.l3[dat.l3$segment == seg, ]
    verdict <- if (nrow(c3) == 0L) "n/a" else as.character(c3$verdict[1L])
    dat.util <- rbind(dat.util, data.frame(segment = seg, acc = acc,
        layer1_dir = dir.l1,
        layer1_sig = sig.l1, layer2_dir = dir.l2,
        layer3_verdict = verdict,
        ## 探索性标记随行传递：MIN_ARM_N 豁免的系列
        ## （user ruling 2026-09-29）在 C_k 中不得计权重
        exploratory = any(a$exploratory),
        joint_weight = vec.joint.weight[match(seg,
            vec.segment)],
        stringsAsFactors = FALSE))
  }
}
write.csv(dat.util, file.path(TAB_DIR, "utilization_matrix.csv"),
  row.names = FALSE)

## PART 2 indeterminate 披露（§4.2 硬性义务） -----------------------------------------------------

## unit 列显式登记分母语义：三层行数单位不同（基因 x 系列 / 段 x 系列 / 段元分析），
## coverage 只作层内披露用，禁止跨层求平均或排序比较。
dat.report <- data.frame(
  source = c("layer1_12gene", "layer2_segment", "layer3_meta"),
  unit = c("gene x primary-contrast rows", "segment x series rows",
    "segment meta rows"),
  n_total = c(nrow(dat.l1.primary), nrow(dat.l2), nrow(dat.l3)),
  stringsAsFactors = FALSE)
dat.report$n_indeterminate <- c(
  sum(dat.l1.primary$direction == "indeterminate"),
  sum(dat.l2$direction == "indeterminate"),
  sum(dat.l3$verdict == "indeterminate"))

## S12/S14（可选补充层）并入
MergeOptional <- function(file.in, source.name, unit.name) {
  if (!file.exists(file.in)) {
    message("[WARN] ", file.in, " absent; ", source.name, " not merged.")
    return(NULL)
  }
  dat <- read.csv(file.in, stringsAsFactors = FALSE)
  if (!("direction" %in% names(dat))) {
    message("[WARN] ", file.in, " has no 'direction' column (not a ",
      "direction table) — treated as absent; ", source.name, " not merged.")
    return(NULL)
  }
  data.frame(source = source.name, unit = unit.name, n_total = nrow(dat),
    n_indeterminate = sum(dat$direction == "indeterminate"),
    stringsAsFactors = FALSE)
}
## S12 空间 QC 表（panel 级，无方向，预期告警）+ S12 scRNA rapa 方向层（含方向）
dat.opt <- list(
  MergeOptional("03_results/S12/table/spatial_panel_qc.csv", "layer4_spatial",
    "spatial panel QC rows"),
  MergeOptional("03_results/S12/table/gse212388_rapa_direction.csv",
    "layer5_scrna_rapa", "scRNA rapa arm x segment rows"),
  MergeOptional("03_results/S12/table/celltype_segment_direction.csv",
    "layer7_celltype_proxy", "celltype fraction x segment proxy rows"))
## S14 独立免疫臂：fgsea NES + BH 判向（NES>0 up / NES<0 down / 其余 indeterminate）
## 口径：NES 正负是 fgsea 官方定义的方向量，BH<0.05 才判向，否则记 indeterminate；
## NA NES（全同号基因集，fgsea 已知行为）一并计入 indeterminate，不静默丢弃。
file.s14 <- "03_results/S14/table/immune_specific_rapa_enrichment.csv"
if (file.exists(file.s14)) {
  dat.im <- read.csv(file.s14, stringsAsFactors = FALSE)
  sig.im <- !is.na(dat.im$NES) & dat.im$padj < 0.05
  dir.im <- ifelse(sig.im & dat.im$NES > 0, "up",
    ifelse(sig.im & dat.im$NES < 0, "down", "indeterminate"))
  dat.opt[[length(dat.opt) + 1L]] <- data.frame(source = "layer6_immune_rapa",
    unit = "immune gene set x rapa contrast rows", n_total = nrow(dat.im),
    n_indeterminate = sum(dir.im == "indeterminate"),
    stringsAsFactors = FALSE)
} else {
  message("[WARN] ", file.s14, " absent; immune arm not merged.")
}
dat.report <- rbind(dat.report, do.call(rbind, dat.opt))
dat.report$fgsea_BH_significant <- NA_integer_
dat.report$cameraPR_BH_significant <- NA_integer_
file.s14sum <- "03_results/S14/table/immune_arm_summary.csv"
if (file.exists(file.s14sum)) {
  dat.sum <- read.csv(file.s14sum, stringsAsFactors = FALSE)
  if (!all(c("rapa_arm_fgsea_BH", "rapa_arm_cameraPR_BH") %in%
      names(dat.sum))) {
    message("[WARN] ", file.s14sum,
      " lacks rapa_arm_fgsea_BH / rapa_arm_cameraPR_BH; dual-route BH not folded.")
  } else {
    row.imm <- dat.sum[!is.na(dat.sum$rapa_arm_fgsea_BH), ]
    if (nrow(row.imm) != 1L) {
      message("[WARN] ", file.s14sum, " has ", nrow(row.imm),
        " rows with rapa_arm_fgsea_BH (expected 1); dual-route BH not folded.")
    } else {
      idx <- which(dat.report$source == "layer6_immune_rapa")
      if (length(idx) == 1L) {
        dat.report$fgsea_BH_significant[idx] <- row.imm$rapa_arm_fgsea_BH
        dat.report$cameraPR_BH_significant[idx] <-
          row.imm$rapa_arm_cameraPR_BH
      } else {
        message("[WARN] layer6_immune_rapa not found in report; ",
          "dual-route BH not folded.")
      }
    }
  }
} else {
  message("[WARN] ", file.s14sum, " absent; dual-route BH not folded.")
}
dat.report$coverage <- round(dat.report$n_indeterminate / dat.report$n_total, 3)
write.csv(dat.report, file.path(TAB_DIR, "indeterminate_report.csv"),
  row.names = FALSE)

## PART 3 图件 -----------------------------------------------------------------------------
## 命名遵守规范 16（描述性英文，禁 Fig 前缀）；图内文字全英文（规范 17）；
## 分面一律用 facet_wrap，不引入任何拼图包。

## (1) 层方向总览：宽表摊平成 (系列 x 段 x 层 x 方向)，按系列分面
dat.util$segment <- factor(dat.util$segment, levels = vec.segment)
dat.dir.long <- do.call(rbind, lapply(unique(dat.util$acc), function(a) {
      d <- dat.util[dat.util$acc == a, ]
      rbind(
        data.frame(acc = a, segment = d$segment, layer = "layer1_12gene",
          direction = d$layer1_dir, stringsAsFactors = FALSE),
        data.frame(acc = a, segment = d$segment, layer = "layer2_segment",
          direction = d$layer2_dir, stringsAsFactors = FALSE),
        data.frame(acc = a, segment = d$segment, layer = "layer3_meta",
          direction = d$layer3_verdict, stringsAsFactors = FALSE))
    }))
dat.dir.long$direction <- as.character(dat.dir.long$direction)
dat.dir.long$direction[is.na(dat.dir.long$direction)] <- "n/a"
dat.dir.long$direction[!dat.dir.long$direction %in%
  c("up", "down", "indeterminate")] <- "n/a"
dat.dir.long$acc <- factor(dat.dir.long$acc)
dat.eff <- dat.dir.long[, c("acc", "segment", "layer")]
dat.eff$eff <- NA_real_
for (ac in unique(dat.l1.primary$acc)) {
  g <- dat.l1.primary$log2FC[dat.l1.primary$acc == ac &
    dat.l1.primary$gene %in% c("IL6", "IL1B", "TNF")]
  if (length(g) > 0L) {
    dat.eff$eff[dat.eff$acc == ac & dat.eff$layer == "layer1_12gene"] <-
      mean(g, na.rm = TRUE)
  }
}
for (ac in unique(dat.l2$acc)) {
  for (sg in vec.segment) {
    rw <- dat.l2$log2FC[dat.l2$acc == ac & dat.l2$segment == sg]
    if (length(rw) > 0L) {
      dat.eff$eff[dat.eff$acc == ac & dat.eff$segment == sg &
        dat.eff$layer == "layer2_segment"] <- rw[1L]
    }
  }
}
for (sg in vec.segment) {
  rw <- dat.l3$pooled[dat.l3$segment == sg]
  if (length(rw) > 0L) {
    dat.eff$eff[dat.eff$segment == sg & dat.eff$layer == "layer3_meta"] <-
      rw[1L]
  }
}
dat.eff$eff_lab <- ifelse(is.na(dat.eff$eff), "", sprintf("%+.2f", dat.eff$eff))
dat.dir.long <- merge(dat.dir.long,
  dat.eff[, c("acc", "segment", "layer", "eff_lab")],
  by = c("acc", "segment", "layer"), all.x = TRUE)
dat.dir.long$txt_col <- ifelse(dat.dir.long$direction == "indeterminate",
  "black", "white")
p.dir <- ggplot(dat.dir.long, aes(x = layer, y = segment, fill = direction)) +
geom_tile(color = "white") +
geom_text(aes(label = eff_lab), colour = dat.dir.long$txt_col, size = 3) +
facet_wrap(~ acc) +
labs(title = "Direction and effect magnitude per readout layer, segment and series",
  subtitle = paste0("Tile label: signed log2FC (layers 1-2) or pooled SMD ",
    "(layer 3); layer 1 dataset-level, layer 3 cross-dataset pooled"),
  x = "Readout layer", y = "Segment") +
scale_fill_manual(values = c("up" = "#C44E52", "down" = "#4C72B0",
  "indeterminate" = "#DDDDDD", "n/a" = "#FFFFFF"),
  name = "Direction") +
sci_theme
sci_save(p.dir, "Segment_Direction_By_Series_Heatmap", width = 11, height = 5)

## (2) indeterminate 占比：分母语义已在 dat.report$unit 登记，此处只按行出条。
## x / fill / facet 三重编码同一个 source 变量，属完全冗余：改单张
## 单色柱状图 + 柱顶数值，避免同一信息画三遍。
## （§4.2），使「占比高」的事实配齐样本量语境与判定阈值，避免光柱无据。
dat.report$source <- factor(dat.report$source, levels = dat.report$source)
dat.report$lab <- sprintf("%.3f\n(n=%d)", dat.report$coverage, dat.report$n_total)
sci_save(ggplot(dat.report, aes(x = source, y = coverage)) +
  geom_col(fill = "#4C72B0") +
  geom_text(aes(label = lab), vjust = -0.4, size = 3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "#C44E52") +
  expand_limits(y = c(0, 1.08)) +
  labs(title = paste0("Indeterminate fraction per readout layer ",
    "(preregistered disclosure)"),
    subtitle = "Dashed line = 0.5 preregistered disclosure threshold (S4.2)",
    x = "Readout layer", y = "Indeterminate fraction") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  sci_theme,
  "Indeterminate_Fraction_By_Readout", width = 11, height = 4.5)

## PART 4 收口 -------------------------------------------------------------------------------

message(paste0("\n===== S13 integrate summary =====\nutilization rows = ", nrow(
      dat.util),
    "\nsegments discounted x0.5 (Jaccard > 0.5) = ",
    if (length(vec.seg.disc) == 0L) "none" else paste(vec.seg.disc,
      collapse = ", ")))
print(dat.report)
if (any(dat.report$coverage > 0.5)) {
  message(
    "[WARN] a readout layer has >50% indeterminate; report this prominently, ",
    "do not downplay it (§4.2 / §10).")
}
message("S13 OK. tables + figure -> 03_results/S13/")
