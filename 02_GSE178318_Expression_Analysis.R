# ===============================================================
# 02_GSE178318_Expression_Analysis.R
# 功能: FBXL12 表达基础分析
# 定位: GSE178318 scRNA-seq 分析流程 Step 2
# 可视化: FBXL12 表达分布、差异分析火山图、组织间表达箱线图
# ===============================================================

# --- 包加载 ------------------------------------------------------

library(Seurat)            # Seurat 对象操作
library(dplyr)            # 数据框转换
library(ggplot2)         # 绘图
library(ggrepel)         # 散点标签避让
library(scales)          # 坐标轴格式化

# --- 全局辅助函数 ------------------------------------------------

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

SciSave <- function(plot_obj, filename, width = 8, height = 6) {
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  pdf(pdf_file, width = width, height = height)
  on.exit(dev.off(), add = TRUE)   # 出错也关闭设备；失败直接抛出（fail-fast，不静默吞错）
  print(plot_obj)
}

# --- 全局常量 -------------------------------------------------

GENE <- "FBXL12"         # 研究目标基因（集中常量，避免散落硬编码）

# Wilson 比例置信区间，返回百分比（0-100）；PART2/PART6 共用单一实现
wilson_ci <- function(k, n, conf.level = 0.95) {
  z <- qnorm(1 - (1 - conf.level) / 2)
  p_hat <- k / n
  denom <- 1 + z^2 / n
  center <- (p_hat + z^2 / (2 * n)) / denom
  margin <- z * sqrt((p_hat * (1 - p_hat) + z^2 / (4 * n)) / n) / denom
  c(lower = pmax(0, center - margin) * 100,
    upper = pmin(1, center + margin) * 100)
}

# --- 目录结构定义 ------------------------------------------------

DATA_DIR <- "X:/CRLM/X/01_QC/table"
FIG_DIR  <- "X:/CRLM/X/02_Expression/Figure"
TAB_DIR  <- "X:/CRLM/X/02_Expression/table"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# ================================================================
# PART 1: 加载注释数据
# 目的: 读取 Step 1 输出的带注释 Seurat 对象
# 数据: seurat_obj_annotated.rds（含 celltype、tissue、malignant 列）
# ================================================================


seurat.obj <- readRDS(file.path(DATA_DIR, "seurat_obj_annotated.rds"))
gc()
seurat.obj <- JoinLayers(seurat.obj)
gc()

# 仅保留 Singlet（与 01 的 doublet 移除 QC 意图一致；当前数据全为 Singlet，此处为安全闸）
seurat.obj <- subset(seurat.obj, subset = doublet == "Singlet")

# P2 守卫：目标基因必须存在于表达矩阵（缺失 → 友好报错，避免后续 FetchData 全 NA 致 FindMarkers 崩溃）
if (!GENE %in% rownames(seurat.obj)) {
  stop(sprintf("目标基因 %s 不在 seurat.obj 的 rownames；请确认基因名/大小写是否与表达矩阵一致", GENE))
}

# ================================================================
# PART 2: FBXL12 表达可视化
# 目的: FeaturePlot + 堆叠柱状图 + RidgePlot
# 数据: seurat_obj_annotated.rds（celltype 列来自手工注释）
# ================================================================


# FeaturePlot
p.fbxl12.umap <- FeaturePlot(
  seurat.obj,
  features     = GENE,
  reduction    = "umap",
  min.cutoff   = "q10",
  max.cutoff   = "q90"
) + sci_theme
SciSave(p.fbxl12.umap, "FBXL12_UMAP", width = 8, height = 6)

# 堆叠柱状图：各细胞类型中 FBXL12+ 细胞比例
seurat.obj$FBXL12 <- FetchData(seurat.obj, vars = GENE)[, 1]
meta <- seurat.obj@meta.data

fbxl12.prop <- meta %>%
  filter(!is.na(celltype) & celltype != "NA") %>%
  mutate(FBXL12_status = ifelse(FBXL12 > 0, "FBXL12+", "FBXL12-")) %>%
  group_by(celltype, FBXL12_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(celltype) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  filter(total >= 10) %>%
  mutate(
    prop         = n / total * 100,          # P1：统一为百分比，与 PART6 的 prop 单位一致
    wilson_lower = wilson_ci(n, total)["lower"],
    wilson_upper = wilson_ci(n, total)["upper"]
  )

celltype.order.bar <- fbxl12.prop %>%
  filter(FBXL12_status == "FBXL12+") %>%
  arrange(desc(prop)) %>%
  pull(celltype)
fbxl12.prop$celltype <- factor(fbxl12.prop$celltype,
                                levels = rev(celltype.order.bar))

p.fbxl12.bar <- ggplot(fbxl12.prop,
  aes(x = celltype, y = prop, fill = FBXL12_status)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(
    data = filter(fbxl12.prop, FBXL12_status == "FBXL12+"),
    aes(ymin = wilson_lower, ymax = wilson_upper),
    width = 0.2, linewidth = 0.5, color = "black"
  ) +
  geom_text(
    data = filter(fbxl12.prop, FBXL12_status == "FBXL12+"),
    aes(label = sprintf("%.1f%% (n=%d)", prop, n)),
    position = position_stack(vjust = 1),
    hjust = -0.15,
    size = 3, color = "black", fontface = "bold"
  ) +
  scale_fill_manual(values = c("FBXL12-" = "#CCCCCC",
                                "FBXL12+" = "#3C5488")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), expand = expansion(mult = c(0, 0.03))) +
  coord_flip() +
  labs(x = NULL, y = "Cell Proportion", fill = "FBXL12 Status") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "top",
    axis.text.y        = element_text(size = 10),
    axis.text.x        = element_text(size = 9)
  )
SciSave(p.fbxl12.bar, "FBXL12_StackedBar_celltype",
        width = 8, height = 9)

write.csv(fbxl12.prop,
  file.path(TAB_DIR, "Table_FBXL12_proportion_by_celltype.csv"),
  row.names = FALSE)

# RidgePlot（log1p 变换，FBXL12 为低表达基因）
seurat.obj$FBXL12_log <- log1p(seurat.obj$FBXL12)
celltype.order <- seurat.obj@meta.data %>%
  group_by(celltype) %>%
  summarise(mean = mean(FBXL12_log, na.rm = TRUE)) %>%
  arrange(desc(mean))
# 局部绘图用因子（不改全局 celltype 列，避免污染下游 03/04/06 读取的字符列）
seurat.obj$celltype_plot <- factor(seurat.obj$celltype,
                                    levels = celltype.order$celltype)

p.ridge <- RidgePlot(seurat.obj, features = "FBXL12_log",
                     group.by = "celltype_plot") +
  scale_x_continuous(limits = c(0, NA)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.direction = "vertical"
  ) +
  labs(x = "log1p(FBXL12 expression)")
SciSave(p.ridge, "FBXL12_RidgePlot", width = 10, height = 12)

gc()

# ================================================================
# PART 3: FBXL12+ 细胞定义
# 目的: 二元分类 FBXL12_high vs FBXL12_low
# 方法: FBXL12 表达量 > 0 标记为 FBXL12_high
# ================================================================


seurat.obj$FBXL12_status <- ifelse(
  seurat.obj$FBXL12 > 0, "FBXL12_high", "FBXL12_low")

# 恶性上皮标识：celltype 为 Epithelial/Tumor 且 malignant_bin==1
# （"恶性上皮"= 肿瘤实质细胞，FBXL12 肿瘤生物学的主要载体）
seurat.obj$epi_malignant <- (
  seurat.obj$celltype == "Epithelial/Tumor" & seurat.obj$malignant_bin == 1
)

p.status <- DimPlot(
  seurat.obj, reduction = "umap",
  group.by = "FBXL12_status",
  cols = c("FBXL12_low" = "grey80", "FBXL12_high" = "red3")
) + sci_theme
SciSave(p.status, "FBXL12_status_UMAP", width = 8, height = 6)

gc()

# ================================================================
# PART 4: FBXL12+ vs FBXL12- 差异表达分析
# 目的: MAST 差异基因 + 火山图
# 方法: FindMarkers（MAST），阈值 FDR<0.05，|logFC|>0.5
# 数据: seurat_obj_annotated.rds（celltype、tissue 列）
# ================================================================

deg.fbxl12 <- FindMarkers(
  seurat.obj,
  ident.1    = "FBXL12_high",
  ident.2    = "FBXL12_low",
  group.by   = "FBXL12_status",
  test.use   = "MAST",
  logfc.threshold = 0.5,
  min.pct    = 0.1
)
deg.fbxl12$gene <- rownames(deg.fbxl12)
deg.fbxl12$significance <- ifelse(
  deg.fbxl12$p_val_adj < 0.05 & deg.fbxl12$avg_log2FC > 0.5, "Up",
  ifelse(
    deg.fbxl12$p_val_adj < 0.05 & deg.fbxl12$avg_log2FC < -0.5, "Down",
    "NS"
  )
)

write.csv(deg.fbxl12,
  file = file.path(TAB_DIR, "Table1_DEG_FBXL12_high_vs_low.csv"),
  row.names = FALSE)

# 火山图
# sqrt 变换：压缩极端 p 值，避免平头截断
deg.fbxl12$neg_log_p <- sqrt(-log10(deg.fbxl12$p_val_adj))

top.up <- head(deg.fbxl12 %>% filter(significance == "Up") %>%
                 arrange(p_val_adj), 5)
fbxl12.row <- deg.fbxl12 %>% filter(gene == GENE)
label.genes <- bind_rows(top.up, fbxl12.row) %>%
  distinct(gene, .keep_all = TRUE)

p.volcano <- ggplot(deg.fbxl12,
  aes(x = avg_log2FC, y = neg_log_p)) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1) +
  scale_color_manual(
    values = c("Up" = "#D62728", "Down" = "#1F77B4", "NS" = "grey70"),
    drop = FALSE
  ) +
  geom_hline(yintercept = sqrt(-log10(0.05)), linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  ggrepel::geom_text_repel(
    data = label.genes, aes(label = gene),
    size = 3.5, fontface = "italic",
    max.overlaps = 15, box.padding = 0.5,
    point.padding = 0.3, segment.color = "grey60",
    segment.size = 0.3
  ) +
  coord_cartesian(xlim = c(-1.5, 4), ylim = c(0, NA)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "plain", size = 14)) +
  labs(
    title = expression(italic("FBXL12")~high~vs~italic("FBXL12")~low),
    x = bquote(Average~log[2]~Fold~Change),
    y = bquote(sqrt(-log[10](adjusted~italic(p)-value))),
    color = "Significance"
  )
SciSave(p.volcano, "Volcano_FBXL12_DEG", width = 8, height = 6.5)
gc()

# ================================================================
# PART 4b: 恶性上皮内 FBXL12 差异表达（肿瘤内 FBXL12 生物学）
# 目的: 将 FBXL12_high vs FBXL12_low 的 DEG 限定在恶性上皮内
#        （PART 4 为全细胞 DEG；此处提取肿瘤实质特异信号）
# 数据: subset(seurat.obj, subset = epi_malignant)
# ================================================================


# P3 守卫：恶性上皮为空/过少时跳过（避免 subset 空集 → FindMarkers/DimPlot 崩溃）
n_mepi <- sum(seurat.obj$epi_malignant, na.rm = TRUE)
if (n_mepi < 3) {
  warning(sprintf("恶性上皮数 = %d (< 3)，跳过 PART 4b 恶性上皮内 DEG 分析", n_mepi))
} else {
obj_mepi <- subset(seurat.obj, subset = epi_malignant)

deg.mepi <- FindMarkers(
  obj_mepi,
  ident.1    = "FBXL12_high",
  ident.2    = "FBXL12_low",
  group.by   = "FBXL12_status",
  test.use   = "MAST",
  logfc.threshold = 0.5,
  min.pct    = 0.1
)
deg.mepi$gene <- rownames(deg.mepi)
deg.mepi$significance <- ifelse(
  deg.mepi$p_val_adj < 0.05 & deg.mepi$avg_log2FC > 0.5, "Up",
  ifelse(
    deg.mepi$p_val_adj < 0.05 & deg.mepi$avg_log2FC < -0.5, "Down",
    "NS"
  )
)
write.csv(deg.mepi,
  file = file.path(TAB_DIR, "Table1b_DEG_FBXL12_high_vs_low_MalignantEpithelium.csv"),
  row.names = FALSE)

# 火山图（恶性上皮）
deg.mepi$neg_log_p <- sqrt(-log10(deg.mepi$p_val_adj))
top.up.me <- head(deg.mepi %>% filter(significance == "Up") %>%
                    arrange(p_val_adj), 5)
fbxl12.row.me <- deg.mepi %>% filter(gene == GENE)
label.genes.me <- bind_rows(top.up.me, fbxl12.row.me) %>%
  distinct(gene, .keep_all = TRUE)

p.volcano.me <- ggplot(deg.mepi,
  aes(x = avg_log2FC, y = neg_log_p)) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1) +
  scale_color_manual(
    values = c("Up" = "#D62728", "Down" = "#1F77B4", "NS" = "grey70"),
    drop = FALSE
  ) +
  geom_hline(yintercept = sqrt(-log10(0.05)), linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  ggrepel::geom_text_repel(
    data = label.genes.me, aes(label = gene),
    size = 3.5, fontface = "italic",
    max.overlaps = 15, box.padding = 0.5,
    point.padding = 0.3, segment.color = "grey60",
    segment.size = 0.3
  ) +
  coord_cartesian(xlim = c(-1.5, 4), ylim = c(0, NA)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "plain", size = 14)) +
  labs(
    title = expression(italic("FBXL12")~high~vs~low~"(malignant epithelium)"),
    x = bquote(Average~log[2]~Fold~Change),
    y = bquote(sqrt(-log[10](adjusted~italic(p)-value))),
    color = "Significance"
  )
SciSave(p.volcano.me, "Volcano_FBXL12_DEG_MalignantEpithelium",
        width = 8, height = 6.5)

# 恶性上皮 UMAP：按 FBXL12_status 着色，展示肿瘤内 FBXL12 空间模式
p.mepi.status <- DimPlot(
  obj_mepi, reduction = "umap",
  group.by = "FBXL12_status",
  cols = c("FBXL12_low" = "grey80", "FBXL12_high" = "red3")
) + sci_theme + ggtitle("Malignant epithelium: FBXL12 status")
SciSave(p.mepi.status, "FBXL12_status_UMAP_MalignantEpithelium",
        width = 8, height = 6)

}  # end else (n_mepi >= 3)
if (exists("obj_mepi")) rm(obj_mepi)
gc()

# ================================================================
# PART 5: 组织间 FBXL12 表达差异
# 目的: PT vs LM 表达差异比较
# 方法: 小提琴图（全量细胞）+ 箱线图（FBXL12+ 细胞）
# 数据: seurat_obj_annotated.rds（tissue 列）
# ================================================================


fbxl12.stats <- seurat.obj@meta.data %>%
  group_by(tissue) %>%
  summarise(
    mean   = mean(FBXL12),
    median = median(FBXL12),
    n      = n()
  )

write.csv(fbxl12.stats,
  file = file.path(TAB_DIR, "Table2_FBXL12_expression_by_tissue.csv"),
  row.names = FALSE)

# Wilcoxon 秩和检验
tissue.levels <- unique(seurat.obj$tissue)
comparisons.top <- combn(tissue.levels, 2, simplify = FALSE)

.tissue_order <- c("PT", "LM")
.tissue_order <- .tissue_order[.tissue_order %in% tissue.levels]

.fill_colors <- c(PT = "#3B4992", LM = "#008280")
.fill_colors <- .fill_colors[names(.fill_colors) %in% .tissue_order]

# 图B：FBXL12+ 细胞箱线图
fbxl12.pos.meta <- seurat.obj@meta.data %>% filter(FBXL12 > 0)
.max_fbxl12_pos <- max(fbxl12.pos.meta$FBXL12)

.wilcox_list_pos <- lapply(seq_along(comparisons.top), function(i) {
  pair <- comparisons.top[[i]]
  .x <- fbxl12.pos.meta$tissue == pair[1]
  .y <- fbxl12.pos.meta$tissue == pair[2]
  if (sum(.x) < 5 || sum(.y) < 5) {
    return(list(comparisons = list(pair), annotations = "n/a",
                y_position = .max_fbxl12_pos + 0.1 + (i - 1) * 0.2))
  }
  p <- wilcox.test(fbxl12.pos.meta$FBXL12[.x],
                    fbxl12.pos.meta$FBXL12[.y])$p.value
  list(
    comparisons = list(pair),
    annotations = ifelse(p < 0.001, "***",
                    ifelse(p < 0.01, "**",
                    ifelse(p < 0.05, "*", "ns"))),
    y_position  = .max_fbxl12_pos + 0.1 + (i - 1) * 0.2
  )
})

set.seed(20260806)   # P5：固定抽样种子，保证抖动散点图可复现
p.box.pos <- ggplot(fbxl12.pos.meta,
  aes(x = tissue, y = FBXL12)) +
  geom_boxplot(aes(fill = tissue),
               width = 0.6, outlier.shape = NA) +
  geom_jitter(
    data = fbxl12.pos.meta %>% sample_frac(min(1, 500 / nrow(fbxl12.pos.meta))),
    mapping = aes(x = tissue, y = FBXL12),
    alpha = 0.2, size = 0.5, width = 0.15,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = lapply(.wilcox_list_pos, function(.item) {
      pair <- .item$comparisons[[1]]
      data.frame(
        x = which(.tissue_order == pair[1]),
        xend = which(.tissue_order == pair[2]),
        y = .item$y_position, yend = .item$y_position
      )
    }) %>% bind_rows(),
    mapping = aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x = vapply(.wilcox_list_pos, function(.item) {
      pair <- .item$comparisons[[1]]
      (which(.tissue_order == pair[1]) +
         which(.tissue_order == pair[2])) / 2
    }, numeric(1)),
    y = vapply(.wilcox_list_pos, function(.item) .item$y_position + 0.15,
               numeric(1)),
    label = vapply(.wilcox_list_pos, "[[", "", "annotations"),
    size = 4, fontface = "bold"
  ) +
  scale_fill_manual(values = .fill_colors) +
  scale_x_discrete(limits = .tissue_order) +
  theme_classic() +
  labs(
    title    = "FBXL12 expression by tissue type",
    subtitle = "FBXL12+ cells only (expression > 0)",
    x = "", y = "FBXL12 expression (log-normalized counts)",
    caption  = sprintf("n = %d FBXL12+ cells; Wilcoxon rank-sum test",
                       nrow(fbxl12.pos.meta))
  ) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    plot.caption    = element_text(size = 8, color = "grey50", hjust = 1),
    plot.subtitle   = element_text(size = 10, color = "grey40")
  )
SciSave(p.box.pos, "FBXL12_positive_by_tissue_boxplot",
        width = 8, height = 6)

# ================================================================
# PART 6: FBXL12+ 细胞比例统计
# 目的: 按组织统计 FBXL12+ 比例 + Wilson CI + Fisher 精确检验
# 方法: 比例差异 Fisher 检验，Wilson 区间估计
# 数据: seurat_obj_annotated.rds（tissue、FBXL12_status 列）
# ================================================================

fbxl12.prop <- seurat.obj@meta.data %>%
  group_by(tissue, FBXL12_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(total = sum(n), prop = n / total * 100) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    wilson_lower = wilson_ci(n, total)["lower"],
    wilson_upper = wilson_ci(n, total)["upper"]
  ) %>%
  ungroup() %>%
  mutate(tissue = factor(tissue,
                               levels = c("PT", "LM")))

# Fisher 精确检验
prop_tests <- lapply(
  combn(tissue.levels, 2, simplify = FALSE),
  function(pair) {
    m <- matrix(c(
      sum(seurat.obj$tissue == pair[1] &
            seurat.obj$FBXL12_status == "FBXL12_high"),
      sum(seurat.obj$tissue == pair[1] &
            seurat.obj$FBXL12_status == "FBXL12_low"),
      sum(seurat.obj$tissue == pair[2] &
            seurat.obj$FBXL12_status == "FBXL12_high"),
      sum(seurat.obj$tissue == pair[2] &
            seurat.obj$FBXL12_status == "FBXL12_low")
    ), nrow = 2, byrow = TRUE)
    ft <- fisher.test(m)
    data.frame(group1 = pair[1], group2 = pair[2],
               p.value = ft$p.value)
  }
) %>% bind_rows() %>%
  mutate(p.signif = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ "ns"
  ))

write.csv(fbxl12.prop,
  file = file.path(TAB_DIR, "Table3_FBXL12_status_proportion.csv"),
  row.names = FALSE)
write.csv(prop_tests,
  file = file.path(TAB_DIR, "Table3b_FBXL12_proportion_fisher_test.csv"),
  row.names = FALSE)

# 堆叠柱状图
fbxl12.high <- fbxl12.prop %>%
  filter(FBXL12_status == "FBXL12_high") %>%
  mutate(
    y_center = prop / 100,
    ci_lower = wilson_lower / 100,
    ci_upper = wilson_upper / 100,
    y_label  = wilson_upper / 100 + 0.015
  )

p.prop <- ggplot(fbxl12.prop,
  aes(x = tissue, y = prop / 100, fill = FBXL12_status)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_errorbar(
    data = fbxl12.high,
    aes(x = tissue, ymin = ci_lower, ymax = ci_upper),
    width = 0.25, color = "grey30", linewidth = 0.5
  ) +
  geom_text(
    data = fbxl12.high,
    aes(x = tissue, y = y_label,
        label = sprintf("n=%d\n(%.1f%%)", n, prop)),
    size = 2.8, color = "black", vjust = 0
  ) +
  scale_fill_manual(
    values = c("FBXL12_low" = "grey80", "FBXL12_high" = "red3"),
    labels = c("FBXL12_low" = "FBXL12-", "FBXL12_high" = "FBXL12+")
  ) +
  scale_y_continuous(limits = c(0, 1.15),
                     expand = expansion(mult = c(0, 0))) +
  theme_classic() +
  labs(title = "FBXL12+ cell proportion by tissue",
       x = "", y = "Proportion", fill = "Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
SciSave(p.prop, "FBXL12_proportion_barplot", width = 8, height = 6)

# ================================================================
# PART 6b: 恶性上皮内 FBXL12+ 比例（按 tissue 分层）
# 目的: 在恶性上皮亚群内统计 FBXL12_high 比例 + Wilson CI + Fisher
# 数据: seurat.obj@meta.data（filter epi_malignant）
# ================================================================


# P3 守卫：恶性上皮为空/过少时跳过（避免 filter 空集 → wilson_ci/Fisher 无意义或绘图崩溃）
n_mepi6 <- sum(seurat.obj$epi_malignant, na.rm = TRUE)
if (n_mepi6 < 3) {
  warning(sprintf("恶性上皮数 = %d (< 3)，跳过 PART 6b 恶性上皮内比例统计", n_mepi6))
  fbxl12.prop.mepi <- data.frame(
    tissue = factor(), FBXL12_status = character(),
    n = integer(), total = numeric(), prop = numeric(),
    wilson_lower = numeric(), wilson_upper = numeric())
  prop_tests.mepi <- data.frame(
    group1 = character(), group2 = character(),
    p.value = numeric(), p.signif = character())
  # 仍写出空表（仅表头），保证下游 03/04 文件契约不破
  write.csv(fbxl12.prop.mepi,
    file = file.path(TAB_DIR, "Table3c_FBXL12_status_proportion_MalignantEpithelium.csv"),
    row.names = FALSE)
  write.csv(prop_tests.mepi,
    file = file.path(TAB_DIR, "Table3d_FBXL12_proportion_fisher_MalignantEpithelium.csv"),
    row.names = FALSE)
} else {
fbxl12.prop.mepi <- seurat.obj@meta.data %>%
  filter(epi_malignant) %>%
  group_by(tissue, FBXL12_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(total = sum(n), prop = n / total * 100) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    wilson_lower = wilson_ci(n, total)["lower"],
    wilson_upper = wilson_ci(n, total)["upper"]
  ) %>%
  ungroup() %>%
  mutate(tissue = factor(tissue, levels = c("PT", "LM")))

# Fisher 精确检验（恶性上皮内 PT vs LM）
me_pairs <- combn(unique(seurat.obj$tissue[seurat.obj$epi_malignant]), 2,
                  simplify = FALSE)
prop_tests.mepi <- lapply(me_pairs, function(pair) {
  m <- matrix(c(
    sum(seurat.obj$epi_malignant & seurat.obj$tissue == pair[1] &
          seurat.obj$FBXL12_status == "FBXL12_high"),
    sum(seurat.obj$epi_malignant & seurat.obj$tissue == pair[1] &
          seurat.obj$FBXL12_status == "FBXL12_low"),
    sum(seurat.obj$epi_malignant & seurat.obj$tissue == pair[2] &
          seurat.obj$FBXL12_status == "FBXL12_high"),
    sum(seurat.obj$epi_malignant & seurat.obj$tissue == pair[2] &
          seurat.obj$FBXL12_status == "FBXL12_low")
  ), nrow = 2, byrow = TRUE)
  ft <- fisher.test(m)
  data.frame(group1 = pair[1], group2 = pair[2], p.value = ft$p.value)
}) %>% bind_rows() %>%
  mutate(p.signif = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ "ns"
  ))

write.csv(fbxl12.prop.mepi,
  file = file.path(TAB_DIR, "Table3c_FBXL12_status_proportion_MalignantEpithelium.csv"),
  row.names = FALSE)
write.csv(prop_tests.mepi,
  file = file.path(TAB_DIR, "Table3d_FBXL12_proportion_fisher_MalignantEpithelium.csv"),
  row.names = FALSE)

fbxl12.high.me <- fbxl12.prop.mepi %>%
  filter(FBXL12_status == "FBXL12_high") %>%
  mutate(y_label = wilson_upper / 100 + 0.015)

p.prop.me <- ggplot(fbxl12.prop.mepi,
  aes(x = tissue, y = prop / 100, fill = FBXL12_status)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_errorbar(
    data = fbxl12.high.me,
    aes(x = tissue, ymin = wilson_lower / 100, ymax = wilson_upper / 100),
    width = 0.25, color = "grey30", linewidth = 0.5
  ) +
  geom_text(
    data = fbxl12.high.me,
    aes(x = tissue, y = y_label,
        label = sprintf("n=%d\n(%.1f%%)", n, prop)),
    size = 2.8, color = "black", vjust = 0
  ) +
  scale_fill_manual(
    values = c("FBXL12_low" = "grey80", "FBXL12_high" = "red3"),
    labels = c("FBXL12_low" = "FBXL12-", "FBXL12_high" = "FBXL12+")
  ) +
  scale_y_continuous(limits = c(0, 1.15), expand = expansion(mult = c(0, 0))) +
  theme_classic() +
  labs(title = "FBXL12+ cell proportion (malignant epithelium)",
       x = "", y = "Proportion", fill = "Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
SciSave(p.prop.me, "FBXL12_proportion_barplot_MalignantEpithelium",
        width = 8, height = 6)
}  # end else (n_mepi6 >= 3)

# ================================================================
# PART 7: 保存结果
# 目的: 输出分析结果 Seurat 对象和表格
# 输出: seurat_obj_FBXL12_status.rds
# ================================================================


saveRDS(seurat.obj,
  file = file.path(TAB_DIR, "seurat_obj_FBXL12_status.rds"))

