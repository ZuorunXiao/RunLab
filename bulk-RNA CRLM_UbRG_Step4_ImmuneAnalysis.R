##############################################################
##  CRLM 泛素相关基因研究 - Step 4: 免疫微环境分析
##
##  【分析逻辑】
##   基于训练集(GSE39582) Tumor样本的风险评分高低分组，分析免疫微环境特征：
##   1. ssGSEA: 28种免疫细胞浸润评分（高低风险组比较）
##   2. GSEA: Hallmark通路（高低风险组比较）
##   3. 免疫检查点基因差异表达
##   4. ESTIMATE算法：免疫评分/基质评分
##   5. 免疫细胞与风险评分相关性热图
##############################################################

setwd("X:/CRLM")
library(msigdbr)
library(GSVA)
library(dplyr)
library(ggplot2)
library(tibble)
library(pheatmap)
library(tidyr)
library(RColorBrewer)
library(ggpubr)
library(fgsea)

select <- dplyr::select  # 防止被AnnotationDbi遮蔽

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)

# ================================================================
# Part 0: 加载数据
# ================================================================
expr39582_norm  <- readRDS("data/processed/expr39582_norm.rds")
surv_df_train   <- readRDS("data/processed/surv_df_train.rds")
group_39582     <- readRDS("data/processed/group_39582.rds")
final_genes     <- readRDS("data/processed/final_genes.rds")
# MSigDB Hallmark gene sets (auto-download GMT if not exists)
hallmark_file <- "data/processed/hallmark_h_gene_sets.rds"
if (!file.exists(hallmark_file)) {
  cat("  Downloading MSigDB Hallmark gene sets...")
  gmt_url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt"
  gmt_path <- "data/processed/h.all.v2024.1.Hs.symbols.gmt"
  download.file(gmt_url, gmt_path, method = "wininet", quiet = TRUE)
  lines <- readLines(gmt_path)
  hallmark_h <- list()
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    set_name <- parts[1]
    genes <- parts[3:length(parts)]
    hallmark_h[[set_name]] <- genes[nchar(genes) > 0]
  }
  saveRDS(hallmark_h, hallmark_file)
  hallmark_gsets <- hallmark_h
  cat(sprintf(" %d Hallmark gene sets loaded\n", length(hallmark_gsets)))
} else {
  hallmark_gsets <- readRDS(hallmark_file)
  cat(sprintf(" %d Hallmark gene sets loaded\n", length(hallmark_gsets)))
}

# 仅使用Tumor样本
tumor_samps <- surv_df_train$sample_id
expr_tumor  <- expr39582_norm[, tumor_samps, drop = FALSE]
risk_group  <- setNames(surv_df_train$risk_group, surv_df_train$sample_id)
risk_score  <- setNames(surv_df_train$risk_score,  surv_df_train$sample_id)

cat("Tumor样本数:", ncol(expr_tumor), "\n")
cat("High组:", sum(risk_group == "High"), "| Low组:", sum(risk_group == "Low"), "\n")

# ================================================================
# Part 1: ssGSEA - 28种免疫细胞浸润
# ================================================================

# 28种免疫细胞特征基因集（来自Charoentong等2017 Cell Reports）
immune_sigs <- list(
  "B_cell"            = c("CD19","CD79A","MS4A1","CD38","CD22"),
  "CD4_T_naive"       = c("CCR7","LEF1","TCF7","SELL"),
  "CD4_T_memory"      = c("CD4","CXCR3","GATA3","RORC"),
  "CD8_T"             = c("CD8A","CD8B","GZMB","PRF1"),
  "T_reg"             = c("FOXP3","IKZF2","IL2RA","ENTPD1"),
  "NK_cell"           = c("NKG7","GNLY","NCAM1","CD244"),
  "M1_macrophage"     = c("CD86","HLA-DRA","NOS2","TNF"),
  "M2_macrophage"     = c("CD163","MRC1","MSR1","TGFB1"),
  "Neutrophil"        = c("FCGR3B","CSF3R","SIGLEC5","FPR1"),
  "Mast_cell"         = c("KIT","TPSAB1","CPA3","GATA1"),
  "DC_activated"      = c("LAMP3","CXCL9","CXCL10","CCL19"),
  "DC_immature"       = c("ITGAX","ITGAM","HLA-DRB1","CD1C"),
  "pDC"               = c("LILRA4","CLEC4C","IGKC","TCF4"),
  "Eosinophil"        = c("CCR3","SIGLEC8","IL5RA","PRSS33"),
  "Monocyte"          = c("CD14","CD68","FCGR3A","ADGRE1"),
  "TFH"               = c("BCL6","CXCR5","ICOS","PDCD1"),
  "Th1"               = c("TBX21","IFNG","IL12RB2","CXCR3"),
  "Th2"               = c("GATA3","IL4","IL13","CCR4"),
  "Th17"              = c("RORC","IL17A","IL23R","CCR6"),
  "CTL"               = c("CD8A","GZMB","PRF1","IFNG","FASLG"),
  "Exhausted_T"       = c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","CD244"),
  "Activated_B"       = c("CD27","CD38","XBP1","IRF4"),
  "Memory_B"          = c("CD27","IL4R","CD80","BACH2"),
  "Plasma"            = c("IGHG1","IGHG4","IGKC","IGLC2"),
  "ILC"               = c("RORC","IL22","NCR1","GATA3"),
  "Stromal"           = c("FAP","PDPN","THY1","CDH11","COL1A1","ACTA2"),
  "EndothelialCell"   = c("PECAM1","VWF","CLDN5","EMCN","CDH5"),
  "CAF"               = c("ACTA2","FAP","PDPN","S100A4","COL1A1","COL3A1")
)

# 过滤基因集中在表达矩阵中不存在的基因
immune_sigs_filt <- lapply(immune_sigs, function(g) intersect(g, rownames(expr_tumor)))
immune_sigs_filt <- immune_sigs_filt[sapply(immune_sigs_filt, length) >= 2]
cat("有效免疫细胞特征数:", length(immune_sigs_filt), "\n")

# ssGSEA评分（兼容新旧版GSVA）
ssgsea_res <- tryCatch({
  # 新版GSVA (>=1.50) API
  param <- ssgseaParam(as.matrix(expr_tumor), immune_sigs_filt, normalize = TRUE)
  gsva(param, verbose = FALSE)
}, error = function(e1) {
  tryCatch({
    # 旧版API
    gsva(as.matrix(expr_tumor), immune_sigs_filt, method = "ssgsea",
         ssgsea.norm = TRUE, verbose = FALSE)
  }, error = function(e2) {
    cat("Both GSVA APIs failed:", e2$message, "\n")
    cat("Using manual ssGSEA approximation...\n")
    # 手动计算均值作为代理
    sapply(immune_sigs_filt, function(g) {
      colMeans(expr_tumor[g, , drop = FALSE])
    }) %>% t()
  })
})

# 将结果矩阵转置：行=样本，列=免疫细胞
ssgsea_df <- as.data.frame(t(ssgsea_res)) %>%
  rownames_to_column("sample_id") %>%
  left_join(data.frame(sample_id = names(risk_group), risk_group = risk_group),
            by = "sample_id")

# 高低风险组免疫细胞差异（Wilcoxon）
cell_types <- names(immune_sigs_filt)
diff_res <- lapply(cell_types, function(ct) {
  high_v <- ssgsea_df[[ct]][ssgsea_df$risk_group == "High"]
  low_v  <- ssgsea_df[[ct]][ssgsea_df$risk_group == "Low"]
  if (length(high_v) < 5 | length(low_v) < 5) return(NULL)
  wt <- wilcox.test(high_v, low_v, exact = FALSE)
  rho <- cor(as.numeric(risk_score[ssgsea_df$sample_id]), ssgsea_df[[ct]],
             method = "spearman", use = "complete.obs")
  data.frame(
    cell_type  = ct,
    mean_high  = mean(high_v, na.rm = TRUE),
    mean_low   = mean(low_v, na.rm = TRUE),
    p_value    = wt$p.value,
    rho        = rho,
    stringsAsFactors = FALSE
  )
})
diff_df <- do.call(rbind, diff_res) %>%
  mutate(adj_p = p.adjust(p_value, method = "BH"),
         significance = ifelse(adj_p < 0.05, "Sig.", "NS"),
         direction  = ifelse(mean_high > mean_low, "Higher in High", "Higher in Low")) %>%
  arrange(p_value)

cat("显著差异免疫细胞 (adj_p<0.05):", sum(diff_df$adj_p < 0.05), "\n")
write.csv(diff_df, "results/tables/Table17_Immune_Diff.csv", row.names = FALSE)
saveRDS(ssgsea_df, "data/processed/ssgsea_df.rds")

# ================================================================
# Part 2: GSEA - Hallmark通路
# ================================================================

# 计算差异t统计量（High vs Low）
high_samps <- names(risk_group)[risk_group == "High"]
low_samps  <- names(risk_group)[risk_group == "Low"]

gene_stats <- apply(expr_tumor, 1, function(x) {
  h <- x[high_samps]
  l <- x[low_samps]
  tryCatch(t.test(h, l)$statistic, error = function(e) 0)
})
gene_rank <- sort(gene_stats, decreasing = TRUE)

# fgsea
set.seed(42)
gsea_res <- fgsea(pathways = hallmark_gsets, stats = gene_rank,
                  minSize = 15, maxSize = 500, nPermSimple = 1000, eps = 0)
gsea_df  <- gsea_res %>%
  as.data.frame() %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = "; ")) %>%
  arrange(pval)

cat("显著Hallmark通路 (FDR<0.05):", sum(gsea_df$padj < 0.05, na.rm = TRUE), "\n")
write.csv(gsea_df %>% select(-leadingEdge), "results/tables/Table18_GSEA_Hallmark.csv", row.names = FALSE)

# ================================================================
# Part 3: 免疫检查点差异
# ================================================================

checkpoint_genes <- c("PDCD1","CD274","CTLA4","LAG3","HAVCR2","TIGIT","SIGLEC15",
                       "CD276","VTCN1","IDO1","IDO2","ENTPD1","NT5E","ADORA2A",
                       "BTLA","CD244","CD160","LAIR1","KLRG1","LGALS9","FOXP3",
                       "TNFSF4","TNFRSF4","CD27","CD70","ICOS","ICOSLG")
checkpoint_avail <- intersect(checkpoint_genes, rownames(expr_tumor))
cat("可用免疫检查点基因:", length(checkpoint_avail), "\n")

cp_res <- lapply(checkpoint_avail, function(g) {
  high_v <- expr_tumor[g, high_samps]
  low_v  <- expr_tumor[g, low_samps]
  wt <- wilcox.test(high_v, low_v, exact = FALSE)
  data.frame(
    gene      = g,
    mean_high = mean(high_v),
    mean_low  = mean(low_v),
    log2FC    = mean(high_v) - mean(low_v),  # 已log2转换的表达，差即log2FC
    p_value   = wt$p.value,
    stringsAsFactors = FALSE
  )
})
cp_df <- do.call(rbind, cp_res) %>%
  mutate(adj_p = p.adjust(p_value, method = "BH"),
         sig   = adj_p < 0.05) %>%
  arrange(p_value)
write.csv(cp_df, "results/tables/Table19_Checkpoint_Diff.csv", row.names = FALSE)

# ================================================================
# Part 4: ESTIMATE 免疫/基质评分（使用 ESTIMATE 包）
# ================================================================
library(estimate)

# --- 4a: 写 GCT 格式临时文件 ---

filt_tmp <- tempfile(pattern = "filt_", fileext = ".gct")
est_tmp  <- tempfile(pattern = "est_",  fileext = ".gct")

# filterCommonGenes 接受 TSV 格式（非GCT），输出 GCT
tsv_tmp <- tempfile(pattern = "expr_", fileext = ".txt")
expr_df <- data.frame(NAME = rownames(expr_tumor),
                       Description = rep("na", nrow(expr_tumor)),
                       expr_tumor, check.names = FALSE, stringsAsFactors = FALSE)
write.table(expr_df, file = tsv_tmp, sep = "\t", quote = FALSE, row.names = FALSE)

# --- 4b: filterCommonGenes ---
filterCommonGenes(input.f = tsv_tmp, output.f = filt_tmp, id = "GeneSymbol")

# --- 4c: estimateScore ---
estimateScore(input.ds = filt_tmp, output.ds = est_tmp, platform = "affymetrix")

# --- 4d: 解析 ESTIMATE 输出 ---
# CORRECT parsing of GCT output:
# Row 1: #1.1 (version) - skip
# Row 2: 4\t566 (dimensions) - skip
# Row 3: NAME, Description, Description.1, GSM971957... (header row) - becomes column names
# Rows 4-7: 4 score rows
# Note: col 1=NAME, col 2=Description, col 3=Description.1 (duplicate), cols 4:569=566 real GSM samples
est_raw <- read.table(est_tmp, skip = 2, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "")
score_cols <- c("ImmuneScore", "StromalScore", "ESTIMATEScore", "TumorPurity")
est_scores <- est_raw[est_raw$NAME %in% score_cols, ]

# Sample IDs from column names (cols 4:569 = 566 GSM IDs)
sample_ids <- colnames(est_scores)[4:ncol(est_scores)]
cat("  ESTIMATE sample count:", length(sample_ids), "\n")

# Build score matrix: rows 4-7, cols 4:569 (all 4 score types x 566 samples)
sm <- as.matrix(est_scores[, 4:ncol(est_scores), drop = FALSE])
rownames(sm) <- est_scores$NAME
cat("  Score matrix:", nrow(sm), "x", ncol(sm), "\n")

# VERIFY: score matrix colnames vs expr_tumor colnames (should match by position)
sm_names <- colnames(sm)
et_names <- colnames(expr_tumor)
cat(sprintf("  ESTIMATE score columns: %d, Expression Tumor samples: %d\n", ncol(sm), ncol(expr_tumor)))
if (ncol(sm) != ncol(expr_tumor)) {
  warning(sprintf("WARNING: ESTIMATE has %d cols but expr_tumor has %d cols!\n", ncol(sm), ncol(expr_tumor)))
  cat("  Only using first", ncol(expr_tumor), "columns for analysis.\n")
}
cat("  Position match (col last):", sm_names[length(sm_names)], "==", et_names[length(et_names)], "=", sm_names[length(sm_names)] == et_names[length(et_names)], "\n")
cat("  All columns match:", sum(sm_names[seq_len(ncol(expr_tumor))] == et_names), "/", ncol(expr_tumor), "\n")

# Build est_df using position-based matching (ESTIMATE preserves column order)
est_df <- data.frame(
  sample_id     = et_names,
  ImmuneScore   = as.numeric(sm["ImmuneScore", ]),
  StromalScore  = as.numeric(sm["StromalScore", ]),
  ESTIMATEScore = as.numeric(sm["ESTIMATEScore", ]),
  TumorPurity   = as.numeric(sm["TumorPurity", ]),
  risk_group    = risk_group[et_names],
  risk_score    = as.numeric(risk_score[et_names]),
  stringsAsFactors = FALSE
)
cat(sprintf("  ESTIMATE完成: ImmuneScore range [%.1f, %.1f], StromalScore range [%.1f, %.1f]\n",
            min(est_df$ImmuneScore, na.rm=TRUE), max(est_df$ImmuneScore, na.rm=TRUE),
            min(est_df$StromalScore, na.rm=TRUE), max(est_df$StromalScore, na.rm=TRUE)))

write.csv(est_df, "results/tables/Table20_ESTIMATE_Scores.csv", row.names = FALSE)
saveRDS(est_df, "data/processed/estimate_df.rds")

# 清理临时文件
unlink(c(tsv_tmp, filt_tmp, est_tmp))

# Part 5: 可视化
# ================================================================

## --- Fig3A: 免疫细胞浸润差异箱线图 ---
cat("绘制 IA1: 免疫细胞差异箱线图\n")
# 选显著差异的细胞类型
sig_cells <- diff_df %>% filter(adj_p < 0.05) %>% pull(cell_type)
if (length(sig_cells) == 0) sig_cells <- diff_df %>% head(8) %>% pull(cell_type)

box_data <- ssgsea_df %>%
  select(sample_id, risk_group, all_of(sig_cells)) %>%
  pivot_longer(cols = all_of(sig_cells), names_to = "cell_type", values_to = "score")

fig3a <- ggplot(box_data, aes(x = risk_group, y = score, fill = risk_group)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, alpha = 0.85) +
  scale_fill_manual(values = c(High = "#E64B35", Low = "#4DBBD5")) +
  facet_wrap(~cell_type, scales = "free_y", ncol = min(4, length(sig_cells))) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     comparisons = list(c("High", "Low")), size = 3) +
  labs(title = "Immune Cell Infiltration: High vs Low Risk",
       x = "", y = "ssGSEA Score") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 8),
        axis.text.x = element_text(angle = 0))

n_row <- ceiling(length(sig_cells) / 4)
pdf("results/figures/IA1_Immune_Boxplot.pdf",
    width = 4 * min(4, length(sig_cells)),
    height = 4 * n_row)
print(fig3a)
dev.off()
cat("  -> IA1 完成\n")

## --- Fig3B: 免疫细胞与风险评分相关性条形图 ---
cat("绘制 IA2: 相关性条形图\n")
# 所有细胞类型的相关性
cor_mat <- cor(as.numeric(risk_score[ssgsea_df$sample_id]),
               ssgsea_df[, cell_types[cell_types %in% colnames(ssgsea_df)], drop=FALSE],
               method = "spearman")
cor_df  <- data.frame(
  cell_type = colnames(cor_mat),
  rho       = as.numeric(cor_mat[1,]),
  p_value   = sapply(colnames(cor_mat), function(ct) {
    if (!(ct %in% colnames(ssgsea_df))) return(NA)
    cor.test(as.numeric(risk_score[ssgsea_df$sample_id]),
             ssgsea_df[[ct]], method = "spearman")$p.value
  })) %>%
  mutate(adj_p = p.adjust(p_value, method = "BH"),
         sig   = cut(adj_p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                     labels = c("***", "**", "*", ""))) %>%
  arrange(rho)

fig3b <- ggplot(cor_df, aes(x = reorder(cell_type, rho), y = rho, fill = rho)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sig, y = ifelse(rho >= 0, rho + 0.01, rho - 0.02)),
            size = 4, color = "black") +
  coord_flip() +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                       midpoint = 0, name = "Spearman rho") +
  labs(title = "Immune Cell - Risk Score Correlation",
       x = "Immune Cell Type", y = "Spearman rho") +
  theme_classic(base_size = 10) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")

pdf("results/figures/IA2_Immune_Correlation.pdf", width = 7, height = 8)
print(fig3b)
dev.off()
cat("  -> IA2 完成\n")

## --- Fig3C: GSEA气泡图 ---
cat("绘制 IA3: GSEA气泡图\n")
gsea_plot_df <- gsea_df %>%
  arrange(NES) %>%
  mutate(
    direction  = ifelse(NES > 0, "Activated", "Suppressed"),
    pathway_s  = gsub("HALLMARK_", "", pathway)
  )

# Take top 30
if (nrow(gsea_plot_df) > 30) gsea_plot_df <- slice(gsea_plot_df, 1:30)

fig3c <- ggplot(gsea_plot_df, aes(x = NES, y = reorder(pathway_s, NES),
                                   color = direction, size = -log10(pval + 1e-10))) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = c("Activated" = "#E64B35", "Suppressed" = "#4DBBD5")) +
  scale_size(range = c(2, 8), name = "-log10(p)") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(title = "GSEA: Hallmark Pathways (High vs Low Risk)",
       x = "Normalized Enrichment Score (NES)",
       y = "", color = "Direction") +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(size = 8))

pdf("results/figures/IA3_GSEA_Hallmark.pdf", width = 9, height = 8)
print(fig3c)
dev.off()
cat("  -> IA3 完成\n")

## --- Fig3D: 免疫检查点差异 ---
cat("绘制 IA4: 检查点差异\n")
cp_sig <- cp_df %>% filter(adj_p < 0.05)
if (nrow(cp_sig) == 0) {
  cat("  无显著差异检查点基因 (adj_p<0.05)，跳过Fig3D\n")
  cp_sig <- cp_df  # 显示全部用于参考
}

cp_plot_df <- expr_tumor[cp_sig$gene, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expr") %>%
  mutate(risk_group = risk_group[sample_id])

fig3d <- ggplot(cp_plot_df, aes(x = risk_group, y = expr, fill = risk_group)) +
  geom_boxplot(outlier.size = 0.4, width = 0.6, alpha = 0.85) +
  scale_fill_manual(values = c(High = "#E64B35", Low = "#4DBBD5")) +
  facet_wrap(~gene, scales = "free_y", ncol = min(5, nrow(cp_sig))) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     comparisons = list(c("High", "Low")), size = 3) +
  labs(title = "Immune Checkpoint Gene Expression",
       subtitle = "High vs Low Risk Group",
       x = "", y = "Expression") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 9))

pdf("results/figures/IA4_Checkpoint_Diff.pdf",
    width = 3 * min(5, nrow(cp_sig)) + 1,
    height = 3 * ceiling(nrow(cp_sig)/5) + 1)
print(fig3d)
dev.off()
cat("  -> IA4 完成\n")

## --- Fig3E: ESTIMATE评分箱线图 ---
cat("绘制 IA5: ESTIMATE评分\n")
est_long <- est_df %>%
  pivot_longer(cols = c("ImmuneScore", "StromalScore", "ESTIMATEScore", "TumorPurity"),
               names_to = "Score_type", values_to = "Score")

p_immune  <- wilcox.test(
  est_df$ImmuneScore[est_df$risk_group=="High"],
  est_df$ImmuneScore[est_df$risk_group=="Low"], exact=FALSE)$p.value
p_stromal <- wilcox.test(
  est_df$StromalScore[est_df$risk_group=="High"],
  est_df$StromalScore[est_df$risk_group=="Low"], exact=FALSE)$p.value
p_estimate <- wilcox.test(
  est_df$ESTIMATEScore[est_df$risk_group=="High"],
  est_df$ESTIMATEScore[est_df$risk_group=="Low"], exact=FALSE)$p.value
p_purity  <- wilcox.test(
  est_df$TumorPurity[est_df$risk_group=="High"],
  est_df$TumorPurity[est_df$risk_group=="Low"], exact=FALSE)$p.value

# 格式化p值避免0.000
format_p <- function(p) {
  if (p < 0.001) return("p<0.001")
  else return(sprintf("p=%.3f", p))
}

fig3e <- ggplot(est_long, aes(x = risk_group, y = Score, fill = risk_group)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, alpha = 0.85) +
  scale_fill_manual(values = c(High = "#E64B35", Low = "#4DBBD5")) +
  facet_wrap(~Score_type, scales = "free_y", ncol = 4) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     comparisons = list(c("High", "Low")), size = 4) +
  labs(title = "ESTIMATE Scores (Yoshihara et al., 2012)",
       subtitle = sprintf("Immune: %s | Stromal: %s | ESTIMATE: %s | Purity: %s",
                          format_p(p_immune), format_p(p_stromal), format_p(p_estimate), format_p(p_purity)),
       x = "", y = "Score") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

pdf("results/figures/IA5_ESTIMATE_Scores.pdf", width = 12, height = 5)
print(fig3e)
dev.off()

cat("Step4 完成！\n")
cat("图表: IA1-IA5 (5张)\n")
cat("表格: Table17-20\n")

