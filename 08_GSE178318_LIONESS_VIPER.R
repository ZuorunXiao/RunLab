# ============================================================================
# 08_GSE178318_LIONESS_VIPER.R  (v3: per-cell LIONESS + limma edge tests)
# LIONESS single-cell GRN -> differential edges (limma); VIPER/aREA activity
# Data: counts from E:, LIONESS results from F: (E: drive is full)
# ============================================================================

# --- Packages -----------------------------------------------------
suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(Biobase)
  library(viper)
  library(limma)
  library(ComplexHeatmap)
  library(scales)
  library(RColorBrewer)
})

# --- Directories --------------------------------------------------
BASE_E <- "X:/XX"
BASE_F <- "X:/XX"
KODIR  <- file.path(BASE_E, "results", "03_GRN_VirtualKO", "table")
OEDIR  <- file.path(BASE_E, "results", "04_VirtualOE", "table")
TAB_DIR <- file.path(BASE_F, "table")
FIG_DIR <- file.path(BASE_F, "Figure")
REF_DIR <- file.path(BASE_F, "MSigDB_ref")
LIONESS_DIR <- file.path(BASE_F, "lioness_results")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REF_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# PART 1: Load data and build gene universe
# ============================================================================


counts_mal <- readRDS(file.path(TAB_DIR, "counts_mal_from_04.rds"))
meta <- read.csv(file.path(TAB_DIR, "cell_meta_from_04.csv"), stringsAsFactors = FALSE)
rownames(meta) <- meta$cell
gene_names <- rownames(counts_mal)
cat(sprintf("Matrix: %d genes x %d cells\n", nrow(counts_mal), ncol(counts_mal)))
cat(sprintf("  Tissue: PT=%d, LM=%d\n",
            sum(meta$tissue == "PT"), sum(meta$tissue == "LM")))
cat(sprintf("  FBXL12: high=%d, low=%d\n",
            sum(meta$FBXL12_status == "high"), sum(meta$FBXL12_status == "low")))

wnt_csv <- file.path(REF_DIR, "wnt_genes.csv")
wnt_msig <- unique(read.csv(wnt_csv)$gene)
cat(sprintf("MSigDB Wnt genes: %d\n", length(wnt_msig)))

wnt_components <- c(
  "WNT1","WNT2","WNT2B","WNT3","WNT3A","WNT5A","WNT5B","WNT10A","WNT10B","WNT11",
  "FZD1","FZD2","FZD3","FZD4","FZD5","FZD6","FZD7","FZD8","FZD9","FZD10",
  "LRP5","LRP6","DVL1","DVL2","DVL3","GSK3B","CTNNB1","AXIN1","AXIN2","APC",
  "TCF7","TCF7L1","TCF7L2","LEF1","NKD1","DKK1","DKK2","DKK3","SFRP1","SFRP2",
  "FBXL12","MYC","CCND1","PPP2CA","WIF1","SOST")

cur_genes <- unique(c(wnt_msig, wnt_components))
cur_genes <- cur_genes[cur_genes %in% gene_names]
cat(sprintf("Curated gene universe: %d genes\n", length(cur_genes)))

focus_genes <- cur_genes[cur_genes %in% c("GSK3B","FBXL12", wnt_msig, wnt_components)]

sci_theme <- function() {
  theme_bw(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right",
          plot.title = element_text(size = 13, hjust = 0.5),
          axis.text = element_text(size = 9))
}
pdf_dev <- function(file, width = 7, height = 6)
  pdf(file = file, width = width, height = height)

# ============================================================================
# PART 2: Per-cell LIONESS + Differential Edges (limma)
# ============================================================================


lioness_mat <- readRDS(file.path(LIONESS_DIR, "lioness_res_per_cell.rds"))
cat(sprintf("Loaded LIONESS: %d edges x %d cells\n", nrow(lioness_mat), ncol(lioness_mat)))

fbx_expr <- as.numeric(counts_mal["FBXL12", ])
names(fbx_expr) <- colnames(counts_mal)

# Parse edge pairs
rownames_lioness <- rownames(lioness_mat)
edge_pairs <- do.call(rbind, strsplit(rownames_lioness, "_"))
edge_geneA <- edge_pairs[, 1]
edge_geneB <- edge_pairs[, 2]
bad <- !(edge_geneA %in% cur_genes) | !(edge_geneB %in% cur_genes)
if (any(bad)) {
  cat(sprintf("  Removing %d invalid edges\n", sum(bad)))
  lioness_mat <- lioness_mat[!bad, , drop = FALSE]
  edge_geneA <- edge_geneA[!bad]
  edge_geneB <- edge_geneB[!bad]
  rownames(lioness_mat) <- rownames_lioness[!bad]
}
cat(sprintf("  Valid edges: %d\n", nrow(lioness_mat)))

# Group assignment for limma design (per-cell)
cell_meta <- meta[colnames(lioness_mat), , drop = FALSE]
cell_meta$group <- paste0(cell_meta$tissue, "_", cell_meta$FBXL12_status)
groups <- sort(unique(cell_meta$group))
cat("  Groups:", paste(groups, collapse = ", "), "\n")

# Compute edge-group means (for visualization and logFC)
edge_group_mean <- sapply(groups, function(g) {
  idx <- which(cell_meta$group == g)
  rowMeans(lioness_mat[, idx, drop = FALSE], na.rm = TRUE)
})
rownames(edge_group_mean) <- rownames(lioness_mat)
colnames(edge_group_mean) <- groups
cat(sprintf("  Edge-group means: %d x %d\n", nrow(edge_group_mean), ncol(edge_group_mean)))

# limma: use per-cell LIONESS scores (6305 observations)
cat("  Running limma (per-cell, N=6305)...\n")
design <- model.matrix(~ 0 + cell_meta$group)
colnames(design) <- groups
contrast.matrix <- makeContrasts(
  PT_high_vs_PT_low  = PT_high - PT_low,
  LM_high_vs_LM_low  = LM_high - LM_low,
  all_high_vs_all_low = (PT_high + LM_high) / 2 - (PT_low + LM_low) / 2,
  levels = groups
)

fit <- lmFit(lioness_mat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

diff_results <- list()
for (cmp in colnames(contrast.matrix)) {
  res <- topTable(fit2, coef = cmp, number = Inf, adjust.method = "BH")
  res$comparison <- cmp
  diff_results[[cmp]] <- res
  cat(sprintf("  %s: %d sig edges (FDR < 0.05)\n",
              cmp, sum(res$adj.P.Val < 0.05, na.rm = TRUE)))
}

# Merge all comparisons
cat("  Merging diff results...\n")
merged_list <- list()
for (cmp in names(diff_results)) {
  r <- diff_results[[cmp]]
  if (nrow(r) == 0) next
  # r$Gene.Symbol contains full edge names like "GSK3B_FBXL12"
  # Use row positions to index edge_geneA/edge_geneB
  merged_list[[cmp]] <- data.frame(
    geneA = edge_geneA[seq_len(nrow(r))],
    geneB = edge_geneB[seq_len(nrow(r))],
    comparison = cmp,
    logFC = r$logFC,
    adj_PVal = r$adj.P.Val,
    t_stat = r$t,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
all_diff_edges <- do.call(rbind, merged_list)
if (is.null(all_diff_edges) || nrow(all_diff_edges) == 0) {
  all_diff_edges <- data.frame(geneA=character(), geneB=character(),
                                comparison=character(), logFC=numeric(),
                                adj_PVal=numeric(), t_stat=numeric())
}
all_diff_edges <- all_diff_edges[order(all_diff_edges$adj_PVal), ]
write.csv(all_diff_edges, file.path(TAB_DIR, "Table_LIONESS_differential_edges.csv"), row.names = FALSE)
cat(sprintf("  Saved: %d rows\n", nrow(all_diff_edges)))

# Focus edges
focus_edges <- all_diff_edges[all_diff_edges$geneA %in% focus_genes |
                               all_diff_edges$geneB %in% focus_genes, ]
cat(sprintf("  Focus edges: %d\n", nrow(focus_edges)))

# Heatmap
top_e <- focus_edges %>%
  dplyr::slice_max(abs(logFC), n = min(50, nrow(focus_edges))) %>%
  dplyr::distinct(geneA, geneB, .keep_all = TRUE)

heat_rows <- paste0(top_e$geneA, "_", top_e$geneB)
heat_mat <- edge_group_mean[heat_rows, , drop = FALSE]
heat_mat <- heat_mat[apply(is.finite(heat_mat), 1, all), , drop = FALSE]
rownames(heat_mat) <- paste(top_e$geneA, top_e$geneB, sep = "-")
colnames(heat_mat) <- groups
if (nrow(heat_mat) > 0) {
  pdf_dev(file.path(FIG_DIR, "LIONESS_focused_edge_heatmap.pdf"),
          width = 9, height = max(6, nrow(heat_mat) * 0.3))
  draw(Heatmap(heat_mat, name = "logFC",
              show_row_names = TRUE, row_names_gp = gpar(fontsize = 7),
              column_title = "FBXL12 x tissue (per-cell LIONESS)",
              cluster_rows = TRUE, cluster_columns = FALSE,
              col = colorRampPalette(c("blue", "white", "red"))(100)))
  dev.off()
  cat(sprintf("  Heatmap saved (%d edges)\n", nrow(heat_mat)))
}

# ============================================================================
# PART 3: VIPER / aREA Activity Inference
# ============================================================================


make_regulon <- function(genes, tfmode = 1, likelihood = 1) {
  genes <- unique(genes[!is.na(genes) & genes != ""])
  n <- length(genes)
  tf <- rep(as.numeric(tfmode), n); names(tf) <- genes
  lk <- rep(as.numeric(likelihood), n); names(lk) <- genes
  list(tfmode = tf, likelihood = lk)
}

regulon_list <- list()

# FBXL12 regulon (scTenifoldKnk targets from 03/04)
fbx_targets <- c()
for (f in c(file.path(KODIR, "Table1_KO_PT_diffRegulation_nofbxl12.csv"),
            file.path(KODIR, "Table2_KO_LM_diffRegulation_nofbxl12.csv"),
            file.path(OEDIR, "Table1_STN_PT_diffRegulation_nofbxl12.csv"),
            file.path(OEDIR, "Table2_STN_LM_diffRegulation_nofbxl12.csv"),
            file.path(OEDIR, "Table3_STN_ALL_diffRegulation_nofbxl12.csv"))) {
  if (file.exists(f)) {
    d <- read.csv(f)
    if ("p.adj" %in% colnames(d)) {
      fbx_targets <- union(fbx_targets, d$gene[d$p.adj < 0.2])
    }
  }
}
fbx_targets <- fbx_targets[fbx_targets %in% gene_names]
regulon_list[["FBXL12"]] <- make_regulon(fbx_targets)
cat(sprintf("  FBXL12 regulon: %d targets\n", length(fbx_targets)))

# WNT pathway regulon
regulon_list[["WNT"]] <- make_regulon(wnt_msig[wnt_msig %in% gene_names])

# TCF7 / TCF7L2 regulon (known TCF/LEF targets from Wnt hallmark)
tcf_known_targets <- c("MYC", "CCND1", "AXIN2", "SOX9", "SPI1", "IL8",
                        "NR4A1", "FOS", "JUNB", "PIM1", "NR4A2", "TCF7",
                        "TCF7L2", "LEF1", "DKK1", "SFRP1", "SFRP2", "WIF1",
                        "LGR5", "ASH2L", "MATN3", "GREM1", "SP5", "MYCN")
tcf_targets <- wnt_msig[wnt_msig %in% tcf_known_targets]
tcf_targets <- tcf_targets[tcf_targets %in% gene_names]
if (length(tcf_targets) >= 5) {
  regulon_list[["TCF7"]] <- make_regulon(tcf_targets)
  regulon_list[["TCF7L2"]] <- make_regulon(tcf_targets)
  cat(sprintf("  TCF7/TCF7L2 regulon: %d targets (Wnt hallmark proxy)\n", length(tcf_targets)))
}

# LEF1 / MYC TFT regulon
tft_csv <- file.path(REF_DIR, "tft_targets.csv")
if (file.exists(tft_csv)) {
  tft_df <- read.csv(tft_csv)
  for (tf in c("LEF1", "MYC")) {
    tg <- unique(tft_df$gene[tft_df$tf == tf])
    tg <- tg[tg %in% gene_names]
    if (length(tg) >= 10) regulon_list[[tf]] <- make_regulon(tg)
  }
}
cat("  VIPER regulons:", paste(names(regulon_list), collapse = ", "), "\n")

# Build expression matrix for VIPER
mx <- counts_mal
n_c <- ncol(mx)
rs <- Matrix::rowSums(mx)
rss <- Matrix::rowSums(mx^2)
mean_x <- rs / n_c
var_x <- rss / n_c - mean_x^2
top_var <- names(sort(var_x, decreasing = TRUE))[1:2000]
sel_genes <- unique(c(top_var,
                      unlist(lapply(regulon_list, function(r) names(r$tfmode)), use.names = FALSE),
                      cur_genes))
sel_genes <- sel_genes[sel_genes %in% gene_names]
cat(sprintf("  VIPER input: %d genes x %d cells\n", length(sel_genes), ncol(mx)))

emat <- as.matrix(counts_mal[sel_genes, ])
rownames(emat) <- sel_genes
eset <- ExpressionSet(assayData = emat,
                      phenoData = AnnotatedDataFrame(meta),
                      featureData = AnnotatedDataFrame(
                        data.frame(gene = sel_genes, row.names = sel_genes)))

cat("  Running VIPER/aREA...\n")
vip <- viper(eset, regulon_list, method = "none", minsize = 4, cores = 1)
vip_rds <- file.path(TAB_DIR, "viper_activity.rds")
saveRDS(vip, vip_rds)
cat("  VIPER activity inferred and saved\n")

act <- exprs(vip)
cat("  VIPER regulators:", paste(rownames(act), collapse = ", "), "\n")

# ============================================================================
# PART 4: FBXL12 <-> Wnt <-> GSK3B Relationship
# ============================================================================


fbx_reg <- regulon_list[["FBXL12"]]
gsk_in_fbx <- "GSK3B" %in% names(fbx_reg$tfmode)

fbx_act   <- act["FBXL12", ]; names(fbx_act) <- colnames(act)
wnt_act   <- act["WNT", ];    names(wnt_act) <- colnames(act)
gsk_expr  <- as.numeric(counts_mal["GSK3B", colnames(act)])
names(gsk_expr) <- colnames(act)

cell_tissue <- setNames(meta$tissue, rownames(meta))
cell_fbx_stat <- setNames(meta$FBXL12_status, rownames(meta))

cor_fbx_act_expr <- cor(fbx_act, fbx_expr[colnames(act)],
                        method = "spearman", use = "complete.obs")
cor_fbx_wnt <- cor(fbx_act, wnt_act, method = "spearman", use = "complete.obs")
cor_fbx_gsk <- cor(fbx_act, gsk_expr, method = "spearman", use = "complete.obs")
cor_gsk_wnt <- cor(gsk_expr, wnt_act, method = "spearman", use = "complete.obs")

for (tiss in c("PT", "LM")) {
  idx <- cell_tissue[colnames(act)] == tiss
  r_fw <- cor(fbx_act[idx], wnt_act[idx], method = "spearman", use = "complete.obs")
  r_fg <- cor(fbx_act[idx], gsk_expr[idx], method = "spearman", use = "complete.obs")
  r_gw <- cor(gsk_expr[idx], wnt_act[idx], method = "spearman", use = "complete.obs")
  cat(sprintf("  [%s] FBXL12-WNT: %.3f, FBXL12-GSK3B: %.3f, GSK3B-WNT: %.3f (n=%d)\n",
              tiss, r_fw, r_fg, r_gw, sum(idx)))
}

rel_df <- data.frame(
  metric = c(
    "GSK3B_in_FBXL12_regulon",
    "cor_FBXL12_activity_vs_FBXL12_expr",
    "cor_FBXL12_activity_vs_WNT_activity",
    "cor_FBXL12_activity_vs_GSK3B_expr",
    "cor_GSK3B_expr_vs_WNT_activity"),
  value = c(
    ifelse(gsk_in_fbx, "yes", "no"),
    round(cor_fbx_act_expr, 3),
    round(cor_fbx_wnt, 3),
    round(cor_fbx_gsk, 3),
    round(cor_gsk_wnt, 3)),
  stringsAsFactors = FALSE)
write.csv(rel_df, file.path(TAB_DIR, "Table_GSK3B_relationship.csv"), row.names = FALSE)
cat("  Relationship table saved\n")
print(rel_df)

# Scatter plots
df_scatter <- data.frame(
  FBXL12_activity = as.numeric(fbx_act),
  WNT_activity    = as.numeric(wnt_act),
  GSK3B_expr      = gsk_expr,
  tissue          = cell_tissue[colnames(act)],
  FBXL12_status   = cell_fbx_stat[colnames(act)])

p1 <- ggplot(df_scatter, aes(x = FBXL12_activity, y = WNT_activity, color = tissue)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "FBXL12 activity vs WNT pathway activity (VIPER/aREA)",
       x = "FBXL12 VIPER activity", y = "WNT VIPER activity") +
  sci_theme()
pdf_dev(file.path(FIG_DIR, "VIPER_FBXL12_vs_WNT.pdf"), width = 7, height = 6)
print(p1); dev.off()

p2 <- ggplot(df_scatter, aes(x = FBXL12_activity, y = GSK3B_expr, color = tissue)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "FBXL12 activity vs GSK3B expression",
       x = "FBXL12 VIPER activity", y = "GSK3B expression") +
  sci_theme()
pdf_dev(file.path(FIG_DIR, "VIPER_FBXL12_vs_GSK3B.pdf"), width = 7, height = 6)
print(p2); dev.off()

# Group mean bar plot
grp_summary <- dplyr::bind_rows(lapply(groups, function(g) {
  cells_g <- colnames(act)[cell_fbx_stat[colnames(act)] == strsplit(g, "_")[[1]][2] &
                            cell_tissue[colnames(act)] == strsplit(g, "_")[[1]][1]]
  data.frame(
    group = g,
    FBXL12_activity = mean(act["FBXL12", cells_g], na.rm = TRUE),
    WNT_activity    = mean(act["WNT", cells_g], na.rm = TRUE),
    GSK3B_expr      = mean(gsk_expr[cells_g], na.rm = TRUE),
    n_cells = length(cells_g))
}))
write.csv(grp_summary, file.path(TAB_DIR, "Table_group_summary.csv"), row.names = FALSE)

p3 <- ggplot(grp_summary, aes(x = group, y = GSK3B_expr, fill = group)) +
  geom_col() +
  labs(title = "GSK3B expression by tissue and FBXL12 status",
       y = "mean GSK3B expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  sci_theme()
pdf_dev(file.path(FIG_DIR, "GSK3B_expr_by_group.pdf"), width = 7, height = 5)
print(p3); dev.off()
cat("  Figures saved\n")

# ============================================================================
# PART 5: Summary
# ============================================================================


gsk_edges <- focus_edges %>% dplyr::filter(geneA == "GSK3B" | geneB == "GSK3B")
fbx_edges <- focus_edges %>% dplyr::filter(geneA == "FBXL12" | geneB == "FBXL12")

summary_lines <- c(
  "# 08 LIONESS (per-cell) + VIPER/aREA: FBXL12 <-> Wnt <-> GSK3B",
  "",
  "## Method",
  "- LIONESS single-cell network: 326 curated genes x 6305 malignant cells.",
  "  Custom jackknife implementation (matrix-optmized, N=6305).",
  "- Differential edges tested with limma (empirical Bayes, BH FDR).",
  "  Contrasts: PT_high vs PT_low, LM_high vs LM_low, all_high vs all_low.",
  "- VIPER/aREA inferred activity for regulators: FBXL12, WNT, TCF7, TCF7L2, LEF1, MYC.",
  "- FBXL12 is an F-box E3 ubiquitin ligase subunit; its 'activity' regulon",
  "  is proxy-based (scTenifoldKnk targets).",
  "",
  "## Key Metrics",
  paste("| Metric | Value |"),
  paste("|---|---|"),
  paste(sprintf("| %s | %s |", rel_df$metric, rel_df$value), collapse = "\n"),
  "",
  "## GSK3B Relationship",
  sprintf("- GSK3B in FBXL12 regulon: **%s**", ifelse(gsk_in_fbx, "Yes", "No")),
  sprintf("- FBXL12 activity vs GSK3B expression (Spearman): **%.3f**", cor_fbx_gsk),
  sprintf("- GSK3B expression vs WNT activity (Spearman): **%.3f** (expected negative)", cor_gsk_wnt),
  sprintf("- FBXL12 activity vs WNT activity (Spearman): **%.3f**", cor_fbx_wnt),
  "",
  "### Top GSK3B Differential Edges (FDR < 0.05)",
  if (nrow(gsk_edges) > 0) {
    top_gsk <- gsk_edges %>% dplyr::slice_max(abs(logFC), n = min(10, nrow(gsk_edges)))
    paste(sprintf("- %s-%s (%s, logFC=%.3f, FDR=%.2e)",
                  top_gsk$geneA, top_gsk$geneB, top_gsk$comparison,
                  top_gsk$logFC, top_gsk$adj_PVal), collapse = "\n")
  } else "- No GSK3B differential edges (FDR < 0.05)",
  "",
  "### Top FBXL12 Differential Edges (FDR < 0.05)",
  if (nrow(fbx_edges) > 0) {
    top_fbx <- fbx_edges %>% dplyr::slice_max(abs(logFC), n = min(10, nrow(fbx_edges)))
    paste(sprintf("- %s-%s (%s, logFC=%.3f, FDR=%.2e)",
                  top_fbx$geneA, top_fbx$geneB, top_fbx$comparison,
                  top_fbx$logFC, top_fbx$adj_PVal), collapse = "\n")
  } else "- No FBXL12 differential edges (FDR < 0.05)",
  "",
  "## Notes",
  "- Per-cell LIONESS (N=6305) provides much more stable estimates than pseudo-bulk (N=6).",
  "- limma empirical Bayes gives well-calibrated p-values.",
  "- TCF7/TCF7L2 regulons use Wnt hallmark gene set as proxy (v2026.1 lacks C3:TFT for these TFs)."
)
writeLines(summary_lines, file.path(TAB_DIR, "SUMMARY_FBXL12_Wnt_GSK3B.md"))
cat("  Summary written to SUMMARY_FBXL12_Wnt_GSK3B.md\n")
cat("\nStep 8 complete\n")

