# =============================================================================
# 05_GSE178318_VirtualScreening.R
# 功能: 基于 GRN 邻接矩阵的基因中心性分析 + 富集分析（GO BP/MF/CC + Reactome）
# 定位: GSE178318 scRNA-seq 分析流程 Step 5
# =============================================================================

# --- 包加载 ------------------------------------------------------
library(Matrix)          # 稀疏矩阵（GRN 邻接矩阵操作）
library(igraph)          # 网络构建与中心性计算
library(reshape2)        # 长宽格式转换
library(org.Hs.eg.db)   # 基因 ID 映射
library(clusterProfiler) # GO 富集分析（本地数据库）
library(ReactomePA)      # Reactome 通路富集（本地数据库，无需网络）
library(ggplot2)        # 可视化
library(ggrepel)        # 标签去重叠
library(dplyr)          # 数据操作

# --- 配色 -------------------------------------------------------
.COL_STABLE   <- "#2166AC"   # 蓝：Stable / high 组
.COL_UNSTABLE <- "#B2182B"   # 红：Unstable / low 组
.COL_PT  <- "#1B9E77"       # 绿：PT 组
.COL_LM  <- "#E7298A"        # 洋红：LM 组
.COL_ALL <- "#7570B3"        # 紫：ALL 组

# --- 辅助函数 ----------------------------------------------------
sci_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position  = "right",
        legend.text      = element_text(size = 10, face = "plain"),
        axis.title       = element_text(size = 12, face = "plain"),
        axis.text        = element_text(size = 10, face = "plain"),
        plot.title       = element_text(size = 14, hjust = 0.5, face = "plain"),
        legend.title     = element_text(size = 11, face = "plain"))

sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf(file.path(FIG_DIR, paste0(filename, ".pdf")), width, height)
  print(plot_obj)
  dev.off()
}

CalcCent <- function(g, grp) {
  deg <- igraph::degree(g, mode = "all", normalized = TRUE)
  bet <- igraph::betweenness(g, normalized = TRUE)
  clo <- igraph::closeness(g, mode = "all", normalized = TRUE)
  eig <- igraph::eigen_centrality(g, directed = FALSE)$vector
  df  <- data.frame(gene = names(deg), degree = deg, betweenness = bet,
                    closeness = clo, eigenvector = eig, stringsAsFactors = FALSE)
  df$group <- grp
  rownames(df) <- NULL
  df
}

CalcDiffCent <- function(cent.hi, cent.lo, grp) {
  df <- merge(
    cent.hi %>% dplyr::select(gene, degree, betweenness, closeness, eigenvector),
    cent.lo  %>% dplyr::select(gene, degree, betweenness, closeness, eigenvector),
    by = "gene", suffix = c(".hi", ".lo")
  )
  df$degree.diff      <- df$degree.hi - df$degree.lo
  df$betweenness.diff <- df$betweenness.hi - df$betweenness.lo
  df$closeness.diff  <- df$closeness.hi - df$closeness.lo
  df$eigenvector.diff <- df$eigenvector.hi - df$eigenvector.lo
  # Z-score 标准化（scale 返回 n×1 矩阵，取第一列）
  df$degree.z      <- as.numeric(scale(df$degree.diff))
  df$betweenness.z <- as.numeric(scale(df$betweenness.diff))
  df$closeness.z   <- as.numeric(scale(df$closeness.diff))
  df$eigenvector.z  <- as.numeric(scale(df$eigenvector.diff))
  df$composite.z <- rowMeans(
    df[, c("degree.z", "betweenness.z", "closeness.z", "eigenvector.z")], na.rm = TRUE)
  df$rank  <- rank(-df$composite.z, ties.method = "min")
  df$group <- grp
  df %>% arrange(rank)
}

# GO 富集（自动放宽阈值，结果写入 CSV；返回 NULL 时打印诊断）
RunGOEnrich <- function(entrez.vec, ont, grp, ont_label) {
  if (length(entrez.vec) < 3) {
    cat(sprintf("    [跳过] %s %s: 仅 %d 个基因，跳过\n", grp, ont, length(entrez.vec)))
    return(NULL)
  }
  res <- enrichGO(gene          = entrez.vec,
                  OrgDb         = org.Hs.eg.db,
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  minGSSize     = 3)
  if (nrow(as.data.frame(res)) == 0) {
    cat(sprintf("    [诊断] %s %s: FDR<0.05 无显著富集\n", grp, ont))
    return(NULL)
  }
  fname <- paste0("Table_", grp, "_", ont, "_Top20Central.csv")
  write.csv(as.data.frame(res), file.path(TAB_DIR, fname), row.names = FALSE)
  if (nrow(as.data.frame(res)) > 0) {
    n_show <- min(15, nrow(as.data.frame(res)))
    p <- dotplot(res, showCategory = n_show) +
      labs(title = paste0(grp, " ", ont_label, " (n=", n_show, ")")) +
      theme(axis.text.y = element_text(size = 10, hjust = 0))
    height_dynamic <- max(6, n_show * 0.7 + 1.5)
    sci_save(p, paste0(grp, "_", ont, "_Top20Central"), width = 9, height = height_dynamic)
  }
  res
}

# ReactomePA 富集（本地数据库，无需网络）
RunPathwayEnrich <- function(entrez.vec, grp) {
  if (length(entrez.vec) < 3) {
    cat(sprintf("    [跳过] %s ReactomePA: 仅 %d 个基因，跳过\n", grp, length(entrez.vec)))
    return(NULL)
  }
  res <- enrichPathway(gene          = entrez.vec,
                       organism      = "human",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       minGSSize     = 3)
  if (nrow(as.data.frame(res)) == 0) {
    cat(sprintf("    [诊断] %s ReactomePA: FDR<0.05 无显著富集\n", grp))
    return(NULL)
  }
  fname <- paste0("Table_", grp, "_ReactomePA_Top20Central.csv")
  write.csv(as.data.frame(res), file.path(TAB_DIR, fname), row.names = FALSE)
  if (nrow(as.data.frame(res)) > 0) {
    n_show <- min(20, nrow(as.data.frame(res)))
    p <- dotplot(res, showCategory = n_show) +
      labs(title = paste0(grp, " ReactomePA Pathway (n=", n_show, ")")) +
      theme(axis.text.y = element_text(size = 10, hjust = 0))
    height_dynamic <- max(6, n_show * 0.7 + 1.5)
    sci_save(p, paste0(grp, "_ReactomePA_Top20Central"), width = 9, height = height_dynamic)
  }
  res
}

# --- 目录 --------------------------------------------------------
TAB_DIR_04 <- "X:/CRLM/XX/results/04_VirtualOE/table"
TAB_DIR    <- "X:/CRLM/XX/results/05_VirtualScreening/table"
FIG_DIR    <- "X:/CRLM/XX/results/05_VirtualScreening/Figure"
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,  recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# PART 1: 加载 GRN 邻接矩阵
# =============================================================================
grn.PT.high  <- readRDS(file.path(TAB_DIR_04, "GRN_PT_high.rds"))
grn.PT.low   <- readRDS(file.path(TAB_DIR_04, "GRN_PT_low.rds"))
grn.LM.high  <- readRDS(file.path(TAB_DIR_04, "GRN_LM_high.rds"))
grn.LM.low   <- readRDS(file.path(TAB_DIR_04, "GRN_LM_low.rds"))
grn.ALL.high <- readRDS(file.path(TAB_DIR_04, "GRN_ALL_high.rds"))
grn.ALL.low  <- readRDS(file.path(TAB_DIR_04, "GRN_ALL_low.rds"))

# =============================================================================
# PART 2: 构建 igraph 网络对象
# =============================================================================
g.PT.high  <- igraph::graph_from_adjacency_matrix(grn.PT.high,  weighted = TRUE, mode = "undirected")
g.PT.low   <- igraph::graph_from_adjacency_matrix(grn.PT.low,   weighted = TRUE, mode = "undirected")
g.LM.high  <- igraph::graph_from_adjacency_matrix(grn.LM.high,  weighted = TRUE, mode = "undirected")
g.LM.low   <- igraph::graph_from_adjacency_matrix(grn.LM.low,   weighted = TRUE, mode = "undirected")
g.ALL.high <- igraph::graph_from_adjacency_matrix(grn.ALL.high, weighted = TRUE, mode = "undirected")
g.ALL.low  <- igraph::graph_from_adjacency_matrix(grn.ALL.low,  weighted = TRUE, mode = "undirected")

# =============================================================================
# PART 3: 计算中心性指标
# =============================================================================
cent.PT.high  <- CalcCent(g.PT.high,  "PT_high")
cent.PT.low   <- CalcCent(g.PT.low,   "PT_low")
cent.LM.high  <- CalcCent(g.LM.high,  "LM_high")
cent.LM.low   <- CalcCent(g.LM.low,   "LM_low")
cent.ALL.high <- CalcCent(g.ALL.high, "ALL_high")
cent.ALL.low  <- CalcCent(g.ALL.low,  "ALL_low")
cent.all <- rbind(cent.PT.high, cent.PT.low, cent.LM.high, cent.LM.low,
                  cent.ALL.high, cent.ALL.low)
write.csv(cent.all, file.path(TAB_DIR, "Table1_centrality_all.csv"), row.names = FALSE)

# =============================================================================
# PART 4: GRN 网络拓扑统计
# P2 ：先用 reshape2::melt 将 topo.df 转换为长格式再作图
# =============================================================================
topo.list <- list(
  PT_high  = g.PT.high,  PT_low  = g.PT.low,
  LM_high  = g.LM.high,  LM_low  = g.LM.low,
  ALL_high = g.ALL.high, ALL_low = g.ALL.low
)
topo.stats <- lapply(topo.list, function(g) {
  data.frame(
    network     = NA,
    n_nodes     = igraph::vcount(g),
    n_edges     = igraph::ecount(g),
    density     = igraph::edge_density(g, loops = FALSE),
    transitivity = igraph::transitivity(g, type = "global"),
    avg_degree  = mean(igraph::degree(g, mode = "all"), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
topo.df <- do.call(rbind, topo.stats)
topo.df$network <- rownames(topo.df)
topo.df$group   <- sub("_(high|low)$", "", topo.df$network)
topo.df$side    <- ifelse(grepl("high$", topo.df$network), "high", "low")
write.csv(topo.df, file.path(TAB_DIR, "Table2_Network_Topology_Summary.csv"), row.names = FALSE)

# melt 转换为长格式后再作图（修复 df.topo.plot 未定义 bug）
df.topo.plot <- reshape2::melt(topo.df[, c("network", "n_edges", "density", "transitivity", "side")],
                               id.vars = c("network", "side"),
                               measure.vars = c("n_edges", "density", "transitivity"),
                               variable.name = "metric",
                               value.name = "value")
df.topo.plot$metric <- factor(df.topo.plot$metric,
                               levels = c("n_edges", "density", "transitivity"),
                               labels = c("Edge Count", "Edge Density", "Transitivity"))

col.side <- c("high" = .COL_STABLE, "low" = .COL_UNSTABLE)
p.topo <- ggplot(df.topo.plot, aes(x = network, y = value, fill = side)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = col.side, name = "Group") +
  labs(title = "GRN Network Topology", x = "", y = "") +
  sci_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text = element_text(size = 11, face = "plain"))
sci_save(p.topo, "Network_Topology_Barplot", width = 9, height = 10)

# =============================================================================
# PART 5: 综合中心性得分（PT / LM / ALL 三组）
# =============================================================================
cent.PT  <- CalcDiffCent(cent.PT.high,  cent.PT.low,  "PT")
cent.LM  <- CalcDiffCent(cent.LM.high,  cent.LM.low,  "LM")
cent.ALL <- CalcDiffCent(cent.ALL.high, cent.ALL.low, "ALL")
write.csv(cent.PT,  file.path(TAB_DIR, "Table3_centrality_PT_ranking.csv"),  row.names = FALSE)
write.csv(cent.LM,  file.path(TAB_DIR, "Table3_centrality_LM_ranking.csv"),  row.names = FALSE)
write.csv(cent.ALL, file.path(TAB_DIR, "Table3_centrality_ALL_ranking.csv"), row.names = FALSE)

# 三组 Top 基因热图（排除已知 scRNA 噪声基因）
# 核糖体亚基为 RPS*/RPL* 家族，精确 %in% 写 "RPS"/"RPL" 会漏掉（如 RPS2、RPL3）；
# 改用前缀 grepl 匹配。同时抽一次过滤复用，避免 PART 5/6 多处重复。
noise_genes_exact <- c("MALAT1", "NEAT1", "ZFAS1", "MT-CO1", "MT-CO2", "MT-CO3",
                       "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5",
                       "MT-CYB", "MT-ATP6", "EEF1A1")
is_noise_gene <- function(g) {
  (g %in% noise_genes_exact) | grepl("^RPS|^RPL", g)
}
pt_filt <- cent.PT[!is_noise_gene(cent.PT$gene), ]
lm_filt <- cent.LM[!is_noise_gene(cent.LM$gene), ]
all_filt <- cent.ALL[!is_noise_gene(cent.ALL$gene), ]
top10.PT  <- head(pt_filt,  10)
top10.LM  <- head(lm_filt,  10)
top10.ALL <- head(all_filt, 10)
union.genes <- unique(c(top10.PT$gene, top10.LM$gene, top10.ALL$gene))
union.genes <- head(union.genes[!is.na(union.genes)], 30)
df.fig1b <- data.frame(
  gene  = union.genes,
  PT    = cent.PT$composite.z[match(union.genes, cent.PT$gene)],
  LM    = cent.LM$composite.z[match(union.genes, cent.LM$gene)],
  ALL   = cent.ALL$composite.z[match(union.genes, cent.ALL$gene)],
  stringsAsFactors = FALSE
)
df.fig1b$gene <- factor(df.fig1b$gene, levels = rev(union.genes))
df.long <- reshape2::melt(df.fig1b, id.vars = "gene",
                           variable.name = "group", value.name = "composite.z")
df.long$group <- factor(df.long$group, levels = c("PT", "LM", "ALL"))
p.top10.group <- ggplot(df.long, aes(x = group, y = gene, fill = composite.z)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = .COL_UNSTABLE, mid = "white",
                       high = .COL_STABLE, midpoint = 0, name = "Composite Z") +
  labs(title = "Top Central Genes: PT vs LM vs ALL\n(excluding scRNA noise genes)",
       x = "", y = "") +
  sci_theme +
  theme(axis.text.y = element_text(size = 8, face = "italic"),
        axis.text.x = element_text(size = 10),
        plot.title  = element_text(hjust = 0.5, size = 11))
sci_save(p.top10.group, "Top10_Centrality_GroupComp",
         width = 6, height = max(5, length(union.genes) * 0.22))

# =============================================================================
# PART 6: 关键节点可视化（PT 组，排除噪声基因）
# =============================================================================
# Top 10 柱状图（复用 PART 5 已算的 top10.PT，避免重复过滤）
p.top10 <- ggplot(top10.PT, aes(x = reorder(gene, -composite.z), y = composite.z)) +
  geom_bar(stat = "identity", fill = .COL_STABLE, width = 0.7) +
  coord_flip() +
  labs(title = "Top 10 Key Nodes (PT high vs low, composite centrality Z-score)",
       x = "", y = "Composite Z-score") +
  sci_theme +
  theme(axis.text.y = element_text(size = 10, face = "italic"))
sci_save(p.top10, "Top10_Centrality", width = 8, height = 5)

# 四种中心性热图（Top 30 PT 组，复用 pt_filt）
top30.pt <- head(pt_filt, 30)
heat.mat <- top30.pt[, c("degree.z", "betweenness.z", "closeness.z", "eigenvector.z")]
rownames(heat.mat) <- top30.pt$gene
heat.long <- reshape2::melt(as.matrix(heat.mat))
colnames(heat.long) <- c("gene", "centrality", "z.score")
heat.long$gene <- factor(heat.long$gene, levels = rev(top30.pt$gene))
p.heat <- ggplot(heat.long, aes(x = centrality, y = gene, fill = z.score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = .COL_UNSTABLE, mid = "white",
                      high = .COL_STABLE, midpoint = 0, name = "Z-score") +
  labs(title = "Top 30 Genes: Centrality Z-scores (PT high vs low)", x = "", y = "") +
  sci_theme +
  theme(axis.text.y = element_text(size = 8, face = "italic"),
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
        plot.title  = element_text(hjust = 0.5, size = 12))
sci_save(p.heat, "Centrality_Heatmap_Top30", width = 7, height = 8)

# =============================================================================
# PART 7: 中心性与差异调控整合（PT / LM / ALL 三组）
# P0 修复：
#   1. x 轴改名为"网络距离偏离度"（非 fold-change）
#   2. 标注逻辑改为按 |composite.z| 排序 Top 15，排除噪声基因
#   3. FC 中性点为 1（distance=均值），非 0
#   4. MALAT1 等噪声基因加注 * 提醒
# =============================================================================

# 优先使用 04 PART5b 产出的剔除被扰动基因 FBXL12 的可用数据;
# 避免将 FBXL12(被扰动基因, distance 人为最大)当作差异调控靶标纳入整合。
read_stn_nofbxl12 <- function(base_csv) {
  p_nf <- sub("\\.csv$", "_nofbxl12.csv", base_csv)
  if (file.exists(p_nf)) return(read.csv(p_nf, stringsAsFactors = FALSE))
  d <- read.csv(base_csv, stringsAsFactors = FALSE)
  d[d$gene != "FBXL12", , drop = FALSE]
}
diff.PT  <- read_stn_nofbxl12(file.path(TAB_DIR_04, "Table1_STN_PT_diffRegulation.csv"))
diff.LM  <- read_stn_nofbxl12(file.path(TAB_DIR_04, "Table2_STN_LM_diffRegulation.csv"))
diff.ALL <- read_stn_nofbxl12(file.path(TAB_DIR_04, "Table3_STN_ALL_diffRegulation.csv"))

IntegrateOne <- function(cent.df, diff.df, grp) {
  out <- merge(cent.df %>% dplyr::select(gene, rank, composite.z, degree.diff),
               diff.df %>% dplyr::select(gene, FC, p.adj),
               by = "gene", all.x = TRUE)
  out$group <- grp
  out
}

cent.PT.full  <- IntegrateOne(cent.PT,  diff.PT,  "PT")
cent.LM.full  <- IntegrateOne(cent.LM,  diff.LM,  "LM")
cent.ALL.full <- IntegrateOne(cent.ALL, diff.ALL, "ALL")
cent.full <- rbind(cent.PT.full, cent.LM.full, cent.ALL.full)

# 按 is_noise 标记；标注基因由下面 lapply 各取 Top 5 by |composite.z| 生成，避免重复过滤
cent.full <- cent.full %>%
  mutate(is_noise = is_noise_gene(gene))

# 每组（PT/LM/ALL）各取 Top 5 by |composite.z|，过滤噪声基因，避免跨组重复
top_label_each <- lapply(list(
  PT  = cent.PT.full,
  LM  = cent.LM.full,
  ALL = cent.ALL.full
), function(df) {
  df %>% filter(!is_noise_gene(gene)) %>%
    arrange(desc(abs(composite.z))) %>%
    slice_head(n = 5) %>%
    pull(gene)
})
top_label <- unique(unlist(top_label_each))

# 噪声基因中 FC > 5 的标注为参考（提醒存在噪声）
noise_highFC <- cent.full %>%
  filter(is_noise & FC > 5) %>%
  pull(gene)

# 标签文字：噪声基因加注 *
cent.full <- cent.full %>%
  mutate(
    label = ifelse(gene %in% top_label, gene, NA_character_),
    label = ifelse(gene %in% noise_highFC,
                   paste0(gene, "*"),  # * = 可能的 scRNA 噪声
                   label)
  )

write.csv(cent.full, file.path(TAB_DIR, "Table4_Centrality_vs_DiffReg.csv"), row.names = FALSE)

# P0 修复：x 轴标签改为"网络距离偏离度"
p.scatter.all <- ggplot(cent.full,
  aes(x = FC, y = -log10(pmax(p.adj, 1e-30)),
      color = group, alpha = !is.na(label))) +
  geom_point(size = 1.2) +
  scale_color_manual(values = c("PT" = .COL_PT, "LM" = .COL_LM, "ALL" = .COL_ALL),
                     name = "Group") +
  scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.25), guide = "none") +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 50,
                  box.padding = 0.3, point.padding = 0.2,
                  fontface = "italic", show.legend = FALSE) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +  # FC=1 为中性（distance=均值），非 0
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  labs(title = "GRN Network Distance Deviation vs Centrality\n(Top 15 by |composite Z| + high-FC noise genes marked with *)",
       x = "Network Distance Deviation (FC = scaled network distance, higher = larger deviation)",
       y = expression(-log[10](FDR))) +
  sci_theme +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 10, face = "plain"))
sci_save(p.scatter.all, "DiffReg_Volcano_AllGroups", width = 8, height = 6)

# 输出标注基因清单供检查
cat(sprintf("  标注基因（每组 Top 5 去重 = 15 unique）: %s\n",
  paste(top_label, collapse = ", ")))
cat(sprintf("  噪声基因（FC>5 供参考标注）: %s\n",
  paste(unique(noise_highFC), collapse = ", ")))

# =============================================================================
# PART 8: 富集分析（PT / LM / ALL 三组分别做 GO BP/MF/CC + ReactomePA）
# P1 修复：分别对三组做富集（不限于 PT 组）
# P1 修复：报告 nominal p vs adjusted p，明确说明哪些结果显著
# =============================================================================

RunGroupEnrich <- function(cent.df, grp) {
  # 取 Top 20 基因，排除噪声基因
  top20 <- cent.df %>%
    filter(!is_noise_gene(gene)) %>%
    slice_head(n = 20) %>%
    pull(gene)
  entrez.map <- bitr(top20, fromType = "SYMBOL", toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db, drop = FALSE)
  entrez.vec <- entrez.map$ENTREZID[!is.na(entrez.map$ENTREZID)]
  write.csv(entrez.map, file.path(TAB_DIR,
    paste0("Table_", grp, "_bitr_mapping.csv")), row.names = FALSE)
  cat(sprintf("  %s: bitr 映射 %d/%d 基因\n", grp, length(entrez.vec), length(top20)))

  # GO 富集
  cat(sprintf("  [%s] GO BP...\n", grp))
  ego.bp <- RunGOEnrich(entrez.vec, "BP", grp, "GO Biological Process")
  cat(sprintf("  [%s] GO MF...\n", grp))
  ego.mf <- RunGOEnrich(entrez.vec, "MF", grp, "GO Molecular Function")
  cat(sprintf("  [%s] GO CC...\n", grp))
  ego.cc <- RunGOEnrich(entrez.vec, "CC", grp, "GO Cellular Component")

  # ReactomePA
  cat(sprintf("  [%s] ReactomePA...\n", grp))
  pathway.res <- RunPathwayEnrich(entrez.vec, grp)

  # 汇总报告（写入单个汇总表）
  all_res <- list()
  if (!is.null(ego.bp) && nrow(as.data.frame(ego.bp)) > 0) {
    d <- as.data.frame(ego.bp); d$ont <- "BP"; all_res[["BP"]] <- d
  }
  if (!is.null(ego.mf) && nrow(as.data.frame(ego.mf)) > 0) {
    d <- as.data.frame(ego.mf); d$ont <- "MF"; all_res[["MF"]] <- d
  }
  if (!is.null(ego.cc) && nrow(as.data.frame(ego.cc)) > 0) {
    d <- as.data.frame(ego.cc); d$ont <- "CC"; all_res[["CC"]] <- d
  }
  if (!is.null(pathway.res) && nrow(as.data.frame(pathway.res)) > 0) {
    d <- as.data.frame(pathway.res); d$ont <- "ReactomePA"; all_res[["ReactomePA"]] <- d
  }

  if (length(all_res) > 0) {
    summary_df <- bind_rows(all_res) %>%
      dplyr::select(group = ont, ID, Description, pvalue, p.adjust, GeneRatio, Count) %>%
      mutate(sig_note = ifelse(p.adjust < 0.05, "*** FDR<0.05 ***",
                               ifelse(pvalue < 0.05, "nominal p<0.05",
                                      "not significant")))
    write.csv(summary_df,
              file.path(TAB_DIR, paste0("Table_", grp, "_Enrichment_Summary.csv")),
              row.names = FALSE)
    cat(sprintf("  %s 富集汇总: FDR<0.05: %d, nominal p<0.05: %d\n",
      grp,
      sum(summary_df$p.adjust < 0.05, na.rm = TRUE),
      sum(summary_df$pvalue < 0.05 & summary_df$p.adjust >= 0.05, na.rm = TRUE)))
  }
}

RunGroupEnrich(cent.PT,  "PT")
RunGroupEnrich(cent.LM,  "LM")
RunGroupEnrich(cent.ALL, "ALL")

cat("\n========== 05 完成 ==========\n")
