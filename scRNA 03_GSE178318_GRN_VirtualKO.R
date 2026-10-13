# ================================================================
# 03_GSE178318_GRN_VirtualKO.R
# 功能: FBXL12 虚拟敲除分析
# 定位: PART1-5 输出CSV+RDS; PART6 富集分析; PART7 火山图;
#        PART8 散点图; PART9 共同基因可视化(含热图)
# ================================================================

# --- 包加载 --------------------------------------------------------
library(Matrix)           # 稀疏矩阵
library(Seurat)           # PART 9 热图需 LayerData (data 层 log1p 归一化表达)
library(dplyr)           # 数据框处理
library(tidyr)            # PART 9 pivot_longer
library(ggplot2)         # 绘图
library(ggrepel)         # 标签避让
library(grid)             # PART 9 grid.draw（gridExtra 不导出该函数）
library(clusterProfiler)  # 功能富集分析
library(enrichplot)      # 富集可视化
library(org.Hs.eg.db)  # 人类基因注释
library(scTenifoldKnk)   # 虚拟敲除
library(ReactomePA)       # Reactome通路富集（KEGG替代）
library(gridExtra)        # 表格可视化
library(pheatmap)        # 热图绘制

# --- 命名空间冲突消解 ----------------------------------------------
# 实测(R 4.6.1): 按上述加载顺序, org.Hs.eg.db/AnnotationDbi 与 clusterProfiler
# 会用 S4 泛型 select / rename 覆盖 dplyr 同名函数, 而这些泛型对 data.frame
# 没有继承方法, 会使 PART5/PART9 中所有 `%>% select(...)` 抛出
#   "unable to find an inherited method for function 'select' for signature x = data.frame"
# 故在此把动词显式绑回 dplyr。bitr 内部自行调用 AnnotationDbi, 不受影响。
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

# --- 全局常量 ------------------------------------------------------
FDR_CUTOFF     <- 0.2      # 显著性FDR阈值
MIN_SIG_GENES  <- 3L        # 富集最少显著基因数
TOP_LABEL_N    <- 15L        # 火山图标签Top N

# --- 目录结构定义 ------------------------------------------------
# 上游依赖: 01 (seurat_obj_annotated.rds) → 02b (count_tumor_PT/LM.rds) → 本脚本
# 注: 02 (表达分析) 是 FBXL12 描述性旁路, 本脚本不消费其任何产物
SEURAT_RDS <- "E:/CRLM/scRNA/GSE178318/results/01_QC/table/seurat_obj_annotated.rds"
FIG_DIR  <- "E:/CRLM/scRNA/GSE178318/results/03_GRN_VirtualKO/Figure"
TAB_DIR  <- "E:/CRLM/scRNA/GSE178318/results/03_GRN_VirtualKO/table"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

options(warn = 1)
options(timeout = 300)

# --- 全局辅助函数 --------------------------------------------------

# 科学绘图主题 (SCI 标准)
sci_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.text      = element_text(size = 9, face = "plain"),
    axis.title       = element_text(size = 11, face = "plain"),
    axis.text        = element_text(size = 9, face = "plain"),
    plot.title       = element_text(size = 13, hjust = 0.5, face = "plain"),
    legend.title     = element_text(size = 10, face = "plain")
  )

# 富集点图主题 (左侧文字标签清晰)
enrich_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.text      = element_text(size = 9, face = "plain"),
    axis.title       = element_text(size = 11, face = "plain"),
    axis.text        = element_text(size = 9, face = "plain"),
    axis.text.y      = element_text(size = 9, face = "plain", hjust = 0),
    plot.title       = element_text(size = 13, hjust = 0.5, face = "plain"),
    legend.title     = element_text(size = 10, face = "plain"),
    plot.margin      = margin(10, 10, 10, 20)
  )

# 读取 02b 产出的 counts 矩阵（带存在性守卫）
ReadTumorCounts <- function(rds_name, label) {
  f <- file.path(TAB_DIR, rds_name)
  if (!file.exists(f)) {
    stop(sprintf("缺少输入文件: %s\n请先运行 02b_GSE178318_Extract_Malignant_Counts.R 生成 %s。",
                 f, rds_name))
  }
  m <- readRDS(f)
  if (!"FBXL12" %in% rownames(m)) {
    stop(sprintf("%s: 矩阵中不含 FBXL12, 无法执行虚拟敲除。", label))
  }
  cat(sprintf("  %s 恶性细胞: %d, 基因数: %d\n", label, ncol(m), nrow(m)))
  m
}

# 稀疏矩阵按行方差（避免 apply() 把 dgCMatrix 转成 dense, 25122x5274 约 1GB）
RowVarsSparse <- function(m) {
  n  <- ncol(m)
  mu <- Matrix::rowMeans(m)
  (Matrix::rowSums(m^2) - n * mu^2) / (n - 1)
}

# 选 top-N 高变基因并强制保留 gKO（基因总数不足时自动收缩, 避免 NA 索引）
SelectHVG <- function(m, n_hvg = 5000L, gko = "FBXL12") {
  n_hvg <- min(n_hvg, nrow(m))
  ord   <- names(sort(RowVarsSparse(m), decreasing = TRUE))
  hvg   <- ord[seq_len(n_hvg)]
  if (!gko %in% hvg) {
    hvg <- c(setdiff(hvg, gko)[seq_len(n_hvg - 1L)], gko)
  }
  hvg <- hvg[!is.na(hvg)]
  m[hvg, , drop = FALSE]
}

# 保存PDF (仅PDF输出, 动态尺寸)
sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  tryCatch({
    pdf(pdf_file, width = width, height = height)
    print(plot_obj)
  }, error = function(e) {
    cat(sprintf("  保存PDF失败 [%s]: %s\n", filename, e$message))
  }, finally = {
    if (dev.cur() > 1) dev.off()
  })
}

# ================================================================
# PART 2: scTenifoldKnk 虚拟敲除（PT 原发灶）
# ================================================================

count.wt.pt <- ReadTumorCounts("count_tumor_PT.rds", "PT")

# HVG 过滤：选 top5000 高变基因 + 强制保留 FBXL12
count.wt.pt.hvg <- SelectHVG(count.wt.pt, 5000L, "FBXL12")

cat(sprintf("  HVG 过滤后基因数: %d (含 FBXL12: %s)\n",
            nrow(count.wt.pt.hvg), "FBXL12" %in% rownames(count.wt.pt.hvg)))

# 虚拟敲除 FBXL12
# 断点续跑：已有 ko_out_PT.rds 直接载入（网络构建约 15-40 min，避免重复计算）
ko_pt_rds <- file.path(TAB_DIR, "ko_out_PT.rds")
if (file.exists(ko_pt_rds)) {
  cat("  载入已有 KO 结果: ko_out_PT.rds\n")
  ko.out.pt <- readRDS(ko_pt_rds)
} else {
  ko.out.pt <- scTenifoldKnk(
    countMatrix = count.wt.pt.hvg,
    qc          = FALSE,
    gKO         = "FBXL12",
    nc_nNet     = 10L,
    nc_nCells   = 500L, # 官方默认：每网随机子样本 500 细胞，集成更稳健且提速约 10x
    nCores      = 1L
  )
  saveRDS(ko.out.pt, file = ko_pt_rds)
}

# 提取差异调控结果
# NOTE: scTenifoldKnk 的 diffRegulation 已含 gene 列（基因名），无需再赋值 rownames
diff.reg.pt <- ko.out.pt$diffRegulation
diff.reg.pt <- diff.reg.pt %>% arrange(p.adj)

cat(sprintf("  PT KO 显著基因 (FDR<0.2): %d\n",
            sum(diff.reg.pt$p.adj < FDR_CUTOFF)))

write.csv(diff.reg.pt,
          file = file.path(TAB_DIR, "Table1_KO_PT_diffRegulation.csv"),
          row.names = FALSE)

rm(count.wt.pt, ko.out.pt)
gc()

# ================================================================
# PART 3: scTenifoldKnk 虚拟敲除（LM 肝转移灶）
# ================================================================


count.wt.lm <- ReadTumorCounts("count_tumor_LM.rds", "LM")

# HVG 过滤
count.wt.lm.hvg <- SelectHVG(count.wt.lm, 5000L, "FBXL12")

cat(sprintf("  HVG 过滤后基因数: %d (含 FBXL12: %s)\n",
            nrow(count.wt.lm.hvg), "FBXL12" %in% rownames(count.wt.lm.hvg)))

# 虚拟敲除 FBXL12
# 断点续跑：已有 ko_out_LM.rds 直接载入
ko_lm_rds <- file.path(TAB_DIR, "ko_out_LM.rds")
if (file.exists(ko_lm_rds)) {
  cat("  载入已有 KO 结果: ko_out_LM.rds\n")
  ko.out.lm <- readRDS(ko_lm_rds)
} else {
  ko.out.lm <- scTenifoldKnk(
    countMatrix = count.wt.lm.hvg,
    qc          = FALSE,
    gKO         = "FBXL12",
    nc_nNet     = 10L,
    nc_nCells   = 500L, # 与 PT 一致：统一子样本规模，保证 PT-vs-LM 可比
    nCores      = 1L
  )
  saveRDS(ko.out.lm, file = ko_lm_rds)
}

diff.reg.lm <- ko.out.lm$diffRegulation
# NOTE: scTenifoldKnk 的 diffRegulation 已含 gene 列（基因名），无需再赋值 rownames
diff.reg.lm <- diff.reg.lm %>% arrange(p.adj)

cat(sprintf("  LM KO 显著基因 (FDR<0.2): %d\n",
            sum(diff.reg.lm$p.adj < FDR_CUTOFF)))

write.csv(diff.reg.lm,
          file = file.path(TAB_DIR, "Table2_KO_LM_diffRegulation.csv"),
          row.names = FALSE)

rm(count.wt.lm, ko.out.lm)
gc()

# ================================================================
# PART 3b: 剔除被扰动基因 FBXL12（官方 scTenifoldKnk/Net FAQ 要求）
# ================================================================
# 被扰动基因(gerKO)的 distance 因"被移除"而人为最大, 是 scTenifoldKnk 工具假象,
# 必须从下游 DEG 表 / 全部可视化 / 下游脚本(06)中剔除。
# 原始 RAW 表 Table1/Table2/Table3 (含 FBXL12) 由 PART2/3/5 照常写出, 此处仅生成
# 无 FBXL12 的 _nofbxl12 副本供 03 后续可视化与 06 消费, 并保留 RAW 内存对象供 PART5 重建原始 merged。


KO_GENE <- "FBXL12"

# 从原始 RAW 表读取(确保与磁盘一致, 兼容断点续传)
diff.reg.pt.raw <- read.csv(file.path(TAB_DIR, "Table1_KO_PT_diffRegulation.csv"))
diff.reg.lm.raw <- read.csv(file.path(TAB_DIR, "Table2_KO_LM_diffRegulation.csv"))

# 剔除 FBXL12, 作为后续所有可视化与下游脚本的工作数据
diff.reg.pt <- diff.reg.pt.raw %>% filter(gene != KO_GENE)
diff.reg.lm <- diff.reg.lm.raw %>% filter(gene != KO_GENE)

# 写出无 FBXL12 的可用数据(磁盘持久化, 供 03 可视化与 06 消费)
write.csv(diff.reg.pt,
          file = file.path(TAB_DIR, "Table1_KO_PT_diffRegulation_nofbxl12.csv"),
          row.names = FALSE)
write.csv(diff.reg.lm,
          file = file.path(TAB_DIR, "Table2_KO_LM_diffRegulation_nofbxl12.csv"),
          row.names = FALSE)

cat(sprintf("  PART 3b 完成: 已剔除 %s (PT n=%d, LM n=%d), _nofbxl12 已写出\n",
            KO_GENE, nrow(diff.reg.pt), nrow(diff.reg.lm)))

# ================================================================
# PART 4: 可视化（PT + LM 对比）
# ================================================================


diff.reg.pt <- read.csv(file.path(TAB_DIR, "Table1_KO_PT_diffRegulation_nofbxl12.csv"))
diff.reg.lm  <- read.csv(file.path(TAB_DIR, "Table2_KO_LM_diffRegulation_nofbxl12.csv"))

# 添加 log10_FC 列（用于火山图 log10 x轴）
diff.reg.pt$log10_FC <- log10(diff.reg.pt$FC)
diff.reg.lm$log10_FC <- log10(diff.reg.lm$FC)

# --- Manhattan / Volcano 通用绘制 ---------------------------------
# 统一处理：y 轴 sqrt(-log10(FDR))；顶部留 12% 空白避免最高点被截断；
#          标注 FDR 最显著的 Top5 基因（ggrepel 自动避让重叠）
manhattan_plot <- function(df, title) {
  df$xidx  <- seq_len(nrow(df))
  # p.adj=0 会使 -log10(p.adj) -> Inf, 用机器最小正数兜底, 避免 ylim 崩溃
  df$p.adj.plot <- pmax(df$p.adj, .Machine$double.xmin)
  y_max    <- max(sqrt(-log10(df$p.adj.plot)), na.rm = TRUE)
  sig_df   <- df[!is.na(df$p.adj) & df$p.adj < FDR_CUTOFF, ]
  top5     <- if (nrow(sig_df) > 0) head(sig_df[order(sig_df$p.adj), ], 5) else sig_df
  p <- ggplot(df, aes(x = xidx, y = sqrt(-log10(p.adj.plot)))) +
    geom_point(aes(color = ifelse(p.adj < FDR_CUTOFF, "sig", "ns")),
               size = 1.2, alpha = 0.8) +
    scale_color_manual(values = c(sig = "#D62728", ns = "grey60"), name = "FDR<0.2") +
    coord_cartesian(ylim = c(0, y_max * 1.35)) +
    labs(title = title, x = "Gene index", y = expression(sqrt(-log[10](FDR)))) +
    sci_theme
  if (nrow(top5) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = top5, aes(label = gene), size = 2.8, fontface = "italic",
      max.overlaps = Inf, min.segment.length = 0, box.padding = 0.35)
  }
  p
}

volcano_plot <- function(df, title) {
  # p.adj=0 兜底, 避免 ylim 崩溃
  df$p.adj.plot <- pmax(df$p.adj, .Machine$double.xmin)
  y_max  <- max(sqrt(-log10(df$p.adj.plot)), na.rm = TRUE)
  sig_df <- df[!is.na(df$p.adj) & df$p.adj < FDR_CUTOFF, ]
  top5   <- if (nrow(sig_df) > 0) head(sig_df[order(sig_df$p.adj), ], 5) else sig_df
  # FC 极小正数会使 log10(FC) 极度负; 左侧截断保护 x 轴可读性,
  # 右侧保留到最大 log10FC 的 1.15 倍, 避免 GAL 等 outlier 标签被截断
  x_low  <- quantile(df$log10_FC, probs = 0.005, na.rm = TRUE)
  x_high <- max(df$log10_FC, na.rm = TRUE) * 1.15
  p <- ggplot(df, aes(x = log10_FC, y = sqrt(-log10(p.adj.plot)))) +
    geom_point(aes(color = ifelse(p.adj < FDR_CUTOFF, "sig", "ns")),
               size = 1.0, alpha = 0.6) +
    scale_color_manual(values = c(sig = "#D62728", ns = "grey70"), name = "FDR<0.2") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    coord_cartesian(xlim = c(x_low, x_high), ylim = c(0, y_max * 1.5)) +
    labs(title = title,
         x = expression(log[10]~"(Fold Change (KO/control))"),
         y = expression(sqrt(-log[10](FDR)))) +
    sci_theme
  if (nrow(top5) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = top5, aes(label = gene), size = 2.8, fontface = "italic",
      max.overlaps = Inf, min.segment.length = 0, box.padding = 0.35)
  }
  p
}

# 曼哈顿图：PT / LM
p.manh.pt <- manhattan_plot(diff.reg.pt, "FBXL12 KO: PT tumor cells")
p.manh.lm <- manhattan_plot(diff.reg.lm, "FBXL12 KO: LM tumor cells")
sci_save(p.manh.pt, "Manhattan_KO_PT", width = 8, height = 7)
sci_save(p.manh.lm, "Manhattan_KO_LM", width = 8, height = 7)

# 火山图：PT / LM
p.vol.pt <- volcano_plot(diff.reg.pt, "FBXL12 KO: PT tumor cells")
p.vol.lm <- volcano_plot(diff.reg.lm, "FBXL12 KO: LM tumor cells")
sci_save(p.vol.pt, "Volcano_KO_PT", width = 8, height = 7)
sci_save(p.vol.lm, "Volcano_KO_LM", width = 8, height = 7)

cat("  PART 4 完成: 曼哈顿图 + 火山图已保存（Top5 标注 + y轴防截断）\n")

# ================================================================
# PART 5: PT vs LM 差异对比
# ================================================================


# 读取差异调控结果
diff.reg.pt <- read.csv(file.path(TAB_DIR, "Table1_KO_PT_diffRegulation_nofbxl12.csv"))
diff.reg.lm  <- read.csv(file.path(TAB_DIR, "Table2_KO_LM_diffRegulation_nofbxl12.csv"))

# 合并 PT 和 LM 结果
sig.pt <- diff.reg.pt$gene[diff.reg.pt$p.adj < FDR_CUTOFF]
sig.lm <- diff.reg.lm$gene[diff.reg.lm$p.adj < FDR_CUTOFF]
common.genes <- intersect(sig.pt, sig.lm)

cat(sprintf("  PT 显著基因: %d\n", length(sig.pt)))
cat(sprintf("  LM 显著基因: %d\n", length(sig.lm)))
cat(sprintf("  共同基因: %d\n", length(common.genes)))

# 合并统计表
diff.reg.pt.sub <- diff.reg.pt %>%
  filter(gene %in% sig.pt) %>%
  select(gene, FC, p.adj)
colnames(diff.reg.pt.sub) <- c("gene", "FC_PT", "p.adj_PT")

diff.reg.lm.sub <- diff.reg.lm %>%
  filter(gene %in% sig.lm) %>%
  select(gene, FC, p.adj)
colnames(diff.reg.lm.sub) <- c("gene", "FC_LM", "p.adj_LM")

merged <- full_join(diff.reg.pt.sub, diff.reg.lm.sub, by = "gene") %>%
  mutate(category = case_when(
    !is.na(FC_PT) & !is.na(FC_LM) ~ "Common",
    !is.na(FC_PT)                 ~ "PT_only",
    TRUE                           ~ "LM_only"
  ))

write.csv(merged,
          file = file.path(TAB_DIR, "Table3_KO_PT_vs_LM_merged_nofbxl12.csv"),
          row.names = FALSE)

# 保留原始 merged(含 FBXL12)以兼容历史产物: 从 PART 3b 的 RAW 内存对象重建
sig.pt.raw <- diff.reg.pt.raw$gene[diff.reg.pt.raw$p.adj < FDR_CUTOFF]
sig.lm.raw <- diff.reg.lm.raw$gene[diff.reg.lm.raw$p.adj < FDR_CUTOFF]
diff.reg.pt.sub.raw <- diff.reg.pt.raw %>%
  filter(gene %in% sig.pt.raw) %>% select(gene, FC, p.adj)
colnames(diff.reg.pt.sub.raw) <- c("gene", "FC_PT", "p.adj_PT")
diff.reg.lm.sub.raw <- diff.reg.lm.raw %>%
  filter(gene %in% sig.lm.raw) %>% select(gene, FC, p.adj)
colnames(diff.reg.lm.sub.raw) <- c("gene", "FC_LM", "p.adj_LM")
merged.raw <- full_join(diff.reg.pt.sub.raw, diff.reg.lm.sub.raw, by = "gene") %>%
  mutate(category = case_when(
    !is.na(FC_PT) & !is.na(FC_LM) ~ "Common",
    !is.na(FC_PT)                 ~ "PT_only",
    TRUE                           ~ "LM_only"
  ))
write.csv(merged.raw,
          file = file.path(TAB_DIR, "Table3_KO_PT_vs_LM_merged.csv"),
          row.names = FALSE)

# 散点图：PT vs LM FC
p.scatter <- ggplot(merged %>% filter(!is.na(FC_PT) & !is.na(FC_LM)),
                    aes(x = FC_PT, y = FC_LM, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey70") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey70") +
  ggrepel::geom_text_repel(aes(label = gene), size = 2.8, fontface = "italic",
                           max.overlaps = Inf, min.segment.length = 0,
                           box.padding = 0.35, show.legend = FALSE) +
  scale_color_manual(values = c(Common = "#D62728",
                                 PT_only = "#E69F00",
                                 LM_only = "#56B4E9")) +
  labs(title = "FBXL12 KO: PT vs LM tumor cells",
       x = "Fold Change (PT)", y = "Fold Change (LM)",
       color = "") +
  sci_theme
sci_save(p.scatter, "Scatter_KO_PT_vs_LM", width = 8, height = 7)

cat("  PART 5 完成: PT vs LM 对比散点图已保存\n")

# ================================================================
# PART 6: 功能富集分析（GO BP/MF/CC + Reactome通路）
# ================================================================


# 富集分析函数
RunEnrichment <- function(gene.list, label) {
  cat(sprintf("\n=== 富集分析: %s (%d genes) ===\n", label, length(gene.list)))
  
  entrez.ids <- bitr(gene.list, fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)
  if (nrow(entrez.ids) < MIN_SIG_GENES) {
    cat(sprintf("  %s: 映射后基因数不足 (%d)，跳过\n", label, nrow(entrez.ids)))
    return(NULL)
  }
  cat(sprintf("  %s: %d/%d 基因成功映射\n", label, nrow(entrez.ids), length(gene.list)))
  
  # --- GO BP ---
  ego.bp <- enrichGO(gene         = entrez.ids$ENTREZID,
                     OrgDb        = org.Hs.eg.db,
                     ont          = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     minGSSize    = 5,
                     maxGSSize    = 500)
  
  ego.bp.df <- as.data.frame(ego.bp)
  if (nrow(ego.bp.df) > 0) {
    n_terms_bp <- min(15, nrow(ego.bp.df))
    p.bp <- dotplot(ego.bp, showCategory = n_terms_bp) + enrich_theme +
      labs(title = paste0(label, " - GO BP"))
    sci_save(p.bp, paste0("Enrich_GO_BP_", label),
             width = 8, height = max(8, n_terms_bp * 0.5))
    write.csv(ego.bp.df,
              file = file.path(TAB_DIR, paste0("Table_Enrich_GO_BP_", label, ".csv")),
              row.names = FALSE)
    cat(sprintf("  %s: GO BP 富集完成（%d terms）\n", label, nrow(ego.bp.df)))
  } else {
    cat(sprintf("  %s: GO BP 无显著富集（跳过PDF）\n", label))
    write.csv(data.frame(),
              file = file.path(TAB_DIR, paste0("Table_Enrich_GO_BP_", label, ".csv")),
              row.names = FALSE)
  }
  
  # --- GO MF ---
  ego.mf <- enrichGO(gene         = entrez.ids$ENTREZID,
                     OrgDb        = org.Hs.eg.db,
                     ont          = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     minGSSize    = 5,
                     maxGSSize    = 500)
  
  ego.mf.df <- as.data.frame(ego.mf)
  if (nrow(ego.mf.df) > 0) {
    n_terms_mf <- min(15, nrow(ego.mf.df))
    p.mf <- dotplot(ego.mf, showCategory = n_terms_mf) + enrich_theme +
      labs(title = paste0(label, " - GO MF"))
    sci_save(p.mf, paste0("Enrich_GO_MF_", label),
             width = 8, height = max(8, n_terms_mf * 0.5))
    write.csv(ego.mf.df,
              file = file.path(TAB_DIR, paste0("Table_Enrich_GO_MF_", label, ".csv")),
              row.names = FALSE)
    cat(sprintf("  %s: GO MF 富集完成（%d terms）\n", label, nrow(ego.mf.df)))
  } else {
    cat(sprintf("  %s: GO MF 无显著富集（跳过PDF）\n", label))
    write.csv(data.frame(),
              file = file.path(TAB_DIR, paste0("Table_Enrich_GO_MF_", label, ".csv")),
              row.names = FALSE)
  }
  
  # --- GO CC ---
  ego.cc <- enrichGO(gene         = entrez.ids$ENTREZID,
                     OrgDb        = org.Hs.eg.db,
                     ont          = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     minGSSize    = 5,
                     maxGSSize    = 500)
  
  ego.cc.df <- as.data.frame(ego.cc)
  if (nrow(ego.cc.df) > 0) {
    n_terms_cc <- min(15, nrow(ego.cc.df))
    p.cc <- dotplot(ego.cc, showCategory = n_terms_cc) + enrich_theme +
      labs(title = paste0(label, " - GO CC"))
    sci_save(p.cc, paste0("Enrich_GO_CC_", label),
             width = 8, height = max(8, n_terms_cc * 0.5))
    write.csv(ego.cc.df,
              file = file.path(TAB_DIR, paste0("Table_Enrich_GO_CC_", label, ".csv")),
              row.names = FALSE)
    cat(sprintf("  %s: GO CC 富集完成（%d terms）\n", label, nrow(ego.cc.df)))
  } else {
    cat(sprintf("  %s: GO CC 无显著富集（跳过PDF）\n", label))
    write.csv(data.frame(),
              file = file.path(TAB_DIR, paste0("Table_Enrich_GO_CC_", label, ".csv")),
              row.names = FALSE)
  }
  
  # --- Reactome通路 (本地数据库,KEGG替代) ---
  er <- enrichPathway(
    gene          = entrez.ids$ENTREZID,
    organism      = "human",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    minGSSize     = 5,
    readable      = TRUE
  )
  
  er.df <- as.data.frame(er)
  if (nrow(er.df) > 0) {
    n_terms_react <- min(15, nrow(er.df))
    p.react <- dotplot(er, showCategory = n_terms_react) + enrich_theme +
      labs(title = paste0(label, " - Reactome Pathway"))
    sci_save(p.react, paste0("Reactome_dotplot_", label),
             width = 10, height = max(7, n_terms_react * 0.5 + 1))
    write.csv(er.df,
              file = file.path(TAB_DIR, paste0("Table_Reactome_", label, ".csv")),
              row.names = FALSE)
    cat(sprintf("  %s: Reactome通路富集完成（%d terms）\n", label, nrow(er.df)))
  } else {
    cat(sprintf("  %s: Reactome通路无显著富集\n", label))
    write.csv(data.frame(),
              file = file.path(TAB_DIR, paste0("Table_Reactome_", label, ".csv")),
              row.names = FALSE)
  }
}

# 运行富集分析
sig.pt <- diff.reg.pt$gene[diff.reg.pt$p.adj < FDR_CUTOFF]
sig.lm <- diff.reg.lm$gene[diff.reg.lm$p.adj < FDR_CUTOFF]

RunEnrichment(sig.pt, "PT_all")
RunEnrichment(sig.lm, "LM_all")
RunEnrichment(common.genes, "Common_PT_LM")

cat("\n========== PART 6 完成 ==========\n")

# ================================================================
# PART 9: 共同基因可视化
# ================================================================


# 读取共同基因列表（排除 FBXL12）
common.genes.clean <- setdiff(common.genes, "FBXL12")
cat(sprintf("  共同基因数(排除FBXL12): %d\n", length(common.genes.clean)))

# 空集守卫: 无共同基因时后续条形图/表格/热图均无意义, tableGrob 亦会报错
if (length(common.genes.clean) == 0) {
  warning("PART 9: PT 与 LM 无共同显著基因(FDR<", FDR_CUTOFF,
          "), 跳过全部共同基因可视化。")
  write.csv(data.frame(gene = character(0), FC_PT = numeric(0), FC_LM = numeric(0),
                       p.adj_PT = numeric(0), p.adj_LM = numeric(0)),
            file = file.path(TAB_DIR, "Table_Common_Genes_Detailed.csv"),
            row.names = FALSE)
}

if (length(common.genes.clean) > 0) {

common.df <- data.frame(
  gene = common.genes.clean,
  FC_PT = NA,
  FC_LM = NA,
  p.adj_PT = NA,
  p.adj_LM = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(common.genes.clean)) {
  g <- common.genes.clean[i]
  if (g %in% diff.reg.pt$gene) {
    r <- diff.reg.pt %>% filter(gene == g)
    common.df$FC_PT[i] <- r$FC
    common.df$p.adj_PT[i] <- r$p.adj
  }
  if (g %in% diff.reg.lm$gene) {
    r <- diff.reg.lm %>% filter(gene == g)
    common.df$FC_LM[i] <- r$FC
    common.df$p.adj_LM[i] <- r$p.adj
  }
}

# 保存详细表格
write.csv(common.df,
          file = file.path(TAB_DIR, "Table_Common_Genes_Detailed.csv"),
          row.names = FALSE)
cat("  详细表格已保存: Table_Common_Genes_Detailed.csv\n")

# --- 图1: Top 20 共同基因 FC 对比条形图 ---
TOP_N_GENES <- min(20L, nrow(common.df))
common.top <- common.df %>%
  arrange(desc(FC_PT + FC_LM)) %>%
  head(TOP_N_GENES) %>%
  mutate(gene = factor(gene, levels = gene[order(FC_PT + FC_LM, decreasing = TRUE)]))

common.long <- common.top %>%
  select(gene, FC_PT, FC_LM) %>%
  tidyr::pivot_longer(
    cols = c(FC_PT, FC_LM),
    names_to = "Tissue",
    values_to = "FC"
  ) %>%
  mutate(
    Tissue = gsub("FC_", "", Tissue),
    Tissue = factor(Tissue, levels = c("PT", "LM")),
    # FC 截断（避免极端值如 MALAT1/IGKC 的 FC>100 导致图表不可读）
    FC = ifelse(FC > 100, 100, FC)
  )

p.bar <- ggplot(common.long, aes(x = gene, y = FC, fill = Tissue)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black", size = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c(PT = "#1F77B4", LM = "#D62728")) +
  coord_flip() +
  labs(
    title = paste0("FBXL12 KO: Top ", TOP_N_GENES, " common DEGs (excl. FBXL12)"),
    x = "Gene symbol",
    y = "Fold Change (KO / control)",
    fill = "Tissue"
  ) +
  sci_theme +
  theme(
    legend.position = "top"
  )
sci_save(p.bar, "Common_Genes_FC_Barplot",
         width = 8, height = max(6, TOP_N_GENES * 0.35))
cat(sprintf("  已保存: Common_Genes_FC_Barplot.pdf (Top %d genes)\n", TOP_N_GENES))

# --- 图2: PT vs LM 散点图（高亮共同基因，带坐标轴截断）---
vol.df <- full_join(
  diff.reg.pt %>% select(gene, FC, p.adj) %>% rename(FC_PT = FC, p.adj_PT = p.adj),
  diff.reg.lm  %>% select(gene, FC, p.adj) %>% rename(FC_LM = FC, p.adj_LM = p.adj),
  by = "gene"
) %>%
  mutate(is.common = ifelse(gene %in% common.genes.clean, "Common", "Specific"))

# 计算坐标轴截断值
FC_CAP_SCATTER <- min(
  quantile(vol.df$FC_PT[vol.df$FC_PT > 0 & !is.na(vol.df$FC_PT)], 0.99, na.rm = TRUE),
  quantile(vol.df$FC_LM[vol.df$FC_LM > 0 & !is.na(vol.df$FC_LM)], 0.99, na.rm = TRUE),
  100
)

p.scatter2 <- ggplot(vol.df %>% filter(!is.na(FC_PT) & !is.na(FC_LM)),
                    aes(x = FC_PT, y = FC_LM, color = is.common)) +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  coord_cartesian(xlim = c(0, FC_CAP_SCATTER), ylim = c(0, FC_CAP_SCATTER)) +
  scale_color_manual(values = c(Common = "#D62728", Specific = "grey70"),
                    name = "Gene type") +
  labs(
    title = "FBXL12 KO: PT vs LM tumor cells",
    x = "Fold Change (PT)",
    y = "Fold Change (LM)"
  ) +
  sci_theme +
  theme(legend.position = "top")
sci_save(p.scatter2, "Common_Genes_Scatter", width = 8, height = 7)
cat("  已保存: Common_Genes_Scatter.pdf\n")

# --- 图3: 详细统计表格 ---
TABLE_TOP_N <- min(30L, nrow(common.df))
table.data <- common.df %>%
  arrange(desc(FC_PT + FC_LM)) %>%
  head(TABLE_TOP_N) %>%
  mutate(
    FC_PT = sprintf("%.2f", FC_PT),
    FC_LM = sprintf("%.2f", FC_LM),
    p.adj_PT = ifelse(p.adj_PT < 0.001, "<0.001", sprintf("%.3f", p.adj_PT)),
    p.adj_LM = ifelse(p.adj_LM < 0.001, "<0.001", sprintf("%.3f", p.adj_LM))
  ) %>%
  select(Gene = gene, `FC (PT)` = FC_PT, `p.adj (PT)` = p.adj_PT,
         `FC (LM)` = FC_LM, `p.adj (LM)` = p.adj_LM)

p.table <- tableGrob(table.data,
                     theme = ttheme_minimal(
                       base_size = 10,
                       padding = unit(c(4, 4), "mm"),
                       colhead = list(fg_params = list(fontface = "bold"))
                     ))
pdf(file.path(FIG_DIR, "Common_Genes_Table.pdf"),
    width = 10, height = max(3, TABLE_TOP_N * 0.28) + 1)
grid.draw(p.table)
dev.off()
cat(sprintf("  已保存: Common_Genes_Table.pdf (Top %d genes)\n", TABLE_TOP_N))

# --- 图4: 单细胞中表达热图（PT vs LM 恶性细胞）---


# 加载 Seurat 对象（路径须与 01 的 tab_dir 一致: results/01_QC/table/）
if (file.exists(SEURAT_RDS)) {
  seurat_obj <- readRDS(SEURAT_RDS)
  cat("  Seurat 对象已加载\n")

  # 组织归属必须用 tissue 列, 不能用 patient:
  # 本数据集每位患者同时贡献 PT 与 LM 样本(COL07/COL12/COL16 均如此),
  # 按 patient 划分会把同一患者的 PT 与 LM 细胞混入同一组, 造成串组。
  # 此处筛选标准与 02b(malignant=="malignant" & tissue=="PT"/"LM")严格一致。
  meta_df <- seurat_obj@meta.data
  need_cols <- c("tissue", "malignant")
  miss <- setdiff(need_cols, colnames(meta_df))
  if (length(miss) > 0) {
    stop(sprintf("Seurat 对象缺少元数据列: %s (应由 01 脚本写入)",
                 paste(miss, collapse = ", ")))
  }

  is_mal   <- meta_df$malignant == "malignant"
  pt_cells <- rownames(meta_df)[is_mal & meta_df$tissue == "PT"]
  lm_cells <- rownames(meta_df)[is_mal & meta_df$tissue == "LM"]

  cat(sprintf("  PT 恶性细胞: %d, LM 恶性细胞: %d\n", length(pt_cells), length(lm_cells)))

  if (length(pt_cells) > 0 && length(lm_cells) > 0) {
    # Seurat v5 对象在 01 中已 JoinLayers (layers: counts, data, scale.data)
    # 热图输入取 data 层 = log1p(CP10K)：既做测序深度归一，又用 log1p 压缩离群值。
    # 若用 counts 原始 UMI，稀疏矩阵的 Z-score 会被极端值拉到 0~20 量级，
    # 绝大多数细胞挤在色标最左端 → 整张图呈深蓝，是视觉失真而非表达真低。
    # NOTE: JoinLayers 返回新对象，不修改原对象（验证通过）。
    expr_mat <- LayerData(seurat_obj, assay = "RNA", layer = "data")

    # 筛选共同基因（在表达矩阵中存在的）
    common_genes_in_mat <- intersect(common.genes.clean, rownames(expr_mat))
    cat(sprintf("  共同基因在 data 层: %d/%d\n", length(common_genes_in_mat), length(common.genes.clean)))

    if (length(common_genes_in_mat) == 0) {
      cat("  警告: 共同基因不在 data 层中，跳过热图\n")
    } else {
      # 细胞降采样: PT+LM 恶性细胞合计数千, 对全部列做层次聚类(cluster_cols)
      # 需要 O(n^2) 距离矩阵且 hclust 极慢, 同时热图列过密不可读。
      # 每组随机抽至多 HEAT_MAX_CELLS 个细胞, 固定随机种子保证可复现。
      HEAT_MAX_CELLS <- 500L
      set.seed(20260807)
      pt_use <- if (length(pt_cells) > HEAT_MAX_CELLS)
                  sample(pt_cells, HEAT_MAX_CELLS) else pt_cells
      lm_use <- if (length(lm_cells) > HEAT_MAX_CELLS)
                  sample(lm_cells, HEAT_MAX_CELLS) else lm_cells
      cat(sprintf("  热图取样: PT %d/%d, LM %d/%d 个细胞\n",
                  length(pt_use), length(pt_cells),
                  length(lm_use), length(lm_cells)))

      # 提取表达量（log1p 归一化表达）
      expr_sub <- as.matrix(expr_mat[common_genes_in_mat, c(pt_use, lm_use)])

      # 添加细胞类型注释
      cell_anno <- data.frame(
        Tissue = c(rep("PT", length(pt_use)), rep("LM", length(lm_use))),
        row.names = c(pt_use, lm_use)
      )

      # 行缩放（Z-score）
      expr_scaled <- t(apply(expr_sub, 1, function(x) {
        if (sd(x) > 0) (x - mean(x)) / sd(x) else rep(0, length(x))
      }))

      # 色标对称截断：以 |Z| 的 99% 分位数封顶，把白色锚定在 Z=0。
      # 否则 pheatmap 将 [min, max] 线性映射，少数离群高值会把白色推离 0，
      # 使多数细胞落进蓝端 → 出现"整片深蓝"的假象。
      zlim <- max(as.numeric(quantile(abs(expr_scaled), 0.99, na.rm = TRUE)), 1e-6)
      expr_scaled[expr_scaled >  zlim] <-  zlim
      expr_scaled[expr_scaled < -zlim] <- -zlim
      heat_breaks <- seq(-zlim, zlim, length.out = 101)
      cat(sprintf("  色标范围: Z ∈ [%.2f, %.2f] (99%% 分位截断)\n", -zlim, zlim))

      # 保存热图为 PDF
      pdf_file <- file.path(FIG_DIR, "Heatmap_Common_Genes_PT_LM.pdf")
      pdf(pdf_file, width = 10, height = max(6, length(common_genes_in_mat) * 0.15))
      pheatmap::pheatmap(
        expr_scaled,
        annotation_col = cell_anno,
        show_rownames = TRUE,
        show_colnames = FALSE,
        fontsize_row = 8,
        fontsize_col = 8,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        color = colorRampPalette(c("#3C5488", "white", "#B93C50"))(100),
        breaks = heat_breaks,
        main = sprintf("Common FBXL12 KO genes: PT (n=%d) vs LM (n=%d) tumor cells\n(row Z-score of log-normalized expression)",
                       length(pt_use), length(lm_use))
      )
      dev.off()
      cat(sprintf("  已保存: Heatmap_Common_Genes_PT_LM.pdf (%d genes x %d cells)\n",
                  length(common_genes_in_mat), ncol(expr_scaled)))
    }
  } else {
    cat("  警告: PT 或 LM 恶性细胞数为0，跳过热图\n")
  }
} else {
  warning("PART 9: Seurat 对象不存在, 跳过热图: ", SEURAT_RDS)
}

cat("\n========== PART 9 完成 ==========\n")
cat("生成的图件:\n")
cat(sprintf("  1. Common_Genes_FC_Barplot.pdf (Top %d 条形图)\n", TOP_N_GENES))
cat("  2. Common_Genes_Scatter.pdf (散点图)\n")
cat(sprintf("  3. Common_Genes_Table.pdf (Top %d 统计表格)\n", TABLE_TOP_N))
cat("  4. Heatmap_Common_Genes_PT_LM.pdf (单细胞中表达热图)\n")

}  # end if (length(common.genes.clean) > 0)


