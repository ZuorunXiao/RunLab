# ===============================================================
# 07_GSE178318_Trajectory_Monocle3.R
# 功能: Monocle3 轨迹分析
# 定位: GSE178318 scRNA-seq 分析流程 Step 7
# 可视化: 轨迹图、伪时间图、差异表达基因热图、趋势图
# ===============================================================

# --- 包加载 ------------------------------------------------

library(dplyr)           # 数据操作
library(ggplot2)         # 数据可视化
library(patchwork)        # 多面板合并图
library(ComplexHeatmap)   # 复杂热图
library(circlize)         # 热图颜色映射（colormap）
library(Seurat)           # 单细胞数据对象
library(monocle3)         # 轨迹分析与伪时间排序
library(igraph)           # 网络分析（principal graph 节点/边统计）

# --- 辅助函数 ----------------------------------------------

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
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  pdf(pdf_file, width = width, height = height)
  print(plot_obj)
  dev.off()
}

# --- 目录结构 ----------------------------------------------

DATA_DIR <- "X:/CRLM/XX/results/02_Expression/table"
FIG_DIR  <- "X:/CRLM/XX/results/07_Trajectory/Figure"
TAB_DIR  <- "X:/CRLM/XX/results/07_Trajectory/table"
dir.create(FIG_DIR, recursive = TRUE)
dir.create(TAB_DIR, recursive = TRUE)

# ================================================================
# PART 1: 加载数据
# 目的：读取预处理后的 Seurat 对象
# 方法：readRDS
# 数据：data/02_Expression/table/seurat_obj_FBXL12_status.rds
# ================================================================

seurat.obj <- readRDS(file.path(DATA_DIR, "seurat_obj_FBXL12_status.rds"))

# ================================================================
# PART 2: Seurat → CellDataSet 转换
# 目的：将 Seurat 对象转为 Monocle3 的 cell_data_set 格式
# 方法：new_cell_data_set() + 赋值 reducedDims / clusters / colData
# 数据：seurat.obj → cds
# ================================================================


cds <- new_cell_data_set(
  expression_data = GetAssayData(seurat.obj, layer = "counts"),
  cell_metadata   = seurat.obj@meta.data,
  gene_metadata   = data.frame(
    gene_short_name = rownames(seurat.obj),
    row.names       = rownames(seurat.obj)
  )
)

# 添加 UMAP 坐标
rd <- reducedDims(cds)
rd$UMAP <- seurat.obj@reductions$umap@cell.embeddings
reducedDims(cds) <- rd

# 添加聚类信息
cds@clusters$UMAP$clusters <- seurat.obj@meta.data$seurat_clusters

# 添加元数据
colData(cds)$tissue         <- seurat.obj$tissue
colData(cds)$patient        <- seurat.obj$patient
colData(cds)$celltype       <- seurat.obj@meta.data$celltype
colData(cds)$FBXL12_status  <- seurat.obj$FBXL12_status

# ================================================================
# PART 3: 轨迹学习
# 目的：基于 UMAP 坐标拟合主曲线
# 方法：cluster_cells() → learn_graph()
# 数据：cds → 输出 6 张轨迹图 PDF
# ================================================================


# ---- 断点 1：轨迹学习结果 ----
# cluster_cells + learn_graph 是最耗时步骤，中断后从此续跑
ckpt3 <- file.path(TAB_DIR, "ckpt3_graph_learned.rds")
if (file.exists(ckpt3)) {
  cat("  [Checkpoint 1] 发现已有轨迹图，载入:", basename(ckpt3), "\n")
  cds <- readRDS(ckpt3)
} else {
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  saveRDS(cds, ckpt3)
}

# ---- 3A: 基础轨迹图（无标注）-------------------------------

p.trajectory <- plot_cells(
  cds,
  color_cells_by = "celltype",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 3
) + labs(title = "Trajectory - Cell Type") + sci_theme

sci_save(p.trajectory, "Trajectory_celltype", width = 10, height = 8)

p.trajectory.tissue <- plot_cells(
  cds,
  color_cells_by = "tissue",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 3
) + labs(title = "Trajectory - Tissue") + sci_theme

sci_save(p.trajectory.tissue, "Trajectory_tissue", width = 10, height = 8)

p.trajectory.fbxl12 <- plot_cells(
  cds,
  color_cells_by = "FBXL12_status",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 3
) + labs(title = "Trajectory - FBXL12 Status") + sci_theme

sci_save(p.trajectory.fbxl12, "Trajectory_FBXL12_status", width = 10, height = 8)

# ---- 3B: Branch point 统计分析（SCI 标准：数字标注降级为附表）----
# 主图不显示节点编号，拓扑信息输出为 CSV 供补充材料
# 提取 principal graph 节点信息
pr_graph <- principal_graph(cds)[["UMAP"]]

# 计算每个节点的细胞类型组成（closest_vertex 映射）
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
node_celltype <- data.frame(
  node_id = as.numeric(closest_vertex[, 1]),
  celltype = colData(cds)$celltype,
  tissue = colData(cds)$tissue
) %>%
  group_by(node_id, celltype) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(node_id) %>%
  mutate(prop = n_cells / sum(n_cells))

# 识别分支点（degree > 2 的节点）
node_degree <- igraph::degree(pr_graph)
branch_nodes <- names(node_degree[node_degree > 2])

branch_summary <- node_celltype %>%
  filter(node_id %in% as.numeric(branch_nodes)) %>%
  arrange(node_id, desc(prop))

write.csv(branch_summary,
  file = file.path(TAB_DIR, "Table_S3_BranchPoint_Celltype.csv"),
  row.names = FALSE)

cat("  Branch point 统计表已保存: Table_S3_BranchPoint_Celltype.csv\n")

# 可选：生成简化版 branch point 图（仅显示前 5 个主要分支点）
if (length(branch_nodes) > 0) {
  top5_branches <- head(branch_summary %>% filter(prop > 0.1) %>% pull(node_id) %>% unique(), 5)
  # 使用 Monocle3 内置标注但限制数量
  p.trajectory.annotated <- plot_cells(
    cds,
    color_cells_by = "celltype",
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = length(top5_branches) <= 8,  # 仅当分支点少时标注
    graph_label_size = 3
  ) + labs(title = "Trajectory (Top Branch Points in Table S3)") +
    theme(legend.position = "right", plot.title = element_text(size = 10, hjust = 0.5))

  sci_save(p.trajectory.annotated, "Trajectory_TopBranches", width = 12, height = 9)
}

# ---- 3C: 3 面板轨迹图（celltype / tissue / FBXL12）-------------
# 注意：pseudotime 需 order_cells 后才存在，完整 4 面板图见 PART 4C

# D2/D3 修复：重复 theme 调用提取为局部变量，避免 7+2 处重复代码
.theme.small <- theme(legend.position = "right",
                      plot.title = element_text(size = 10, hjust = 0.5))
.theme.tag   <- theme(plot.tag = element_text(size = 14, face = "bold"))

p.trajectory.celltype <- plot_cells(
  cds, color_cells_by = "celltype",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  graph_label_size = 2
) + .theme.small

p.trajectory.tissue.p <- plot_cells(
  cds, color_cells_by = "tissue",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  graph_label_size = 2
) + .theme.small

p.trajectory.fbxl12.p <- plot_cells(
  cds, color_cells_by = "FBXL12_status",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  graph_label_size = 2
) + .theme.small

p.multi.3c <- (p.trajectory.celltype | p.trajectory.tissue.p) /
              (p.trajectory.fbxl12.p | plot_spacer()) +
              plot_annotation(tag_levels = "A") &
              .theme.tag

sci_save(p.multi.3c, "Trajectory_MultiPanel", width = 14, height = 12)

# ================================================================
# PART 4: 伪时间排序
# 目的：以 EPCAM+ 上皮细胞为根计算沿轨迹的伪时间
# 方法：order_cells(root_pr_nodes) + closest_vertex 图节点方案
# 数据：cds → 输出 Pseudotime PDF + 轨迹统计表
# ================================================================


umap.coords <- reducedDims(cds)$UMAP
closest_vertex <- as.matrix(
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex)
closest_vertex <- closest_vertex[colnames(cds), , drop = FALSE]

# ---- 根节点选择：最大 partition ----
# 选细胞数最多的 principal graph partition 作为轨迹主干，
# 以其最接近图顶点的众数为根，最大化伪时间覆盖（避免仅覆盖小分支）
partitions <- cds@clusters$UMAP$partitions
partitions <- partitions[!is.na(partitions)]
part.tab <- sort(table(partitions), decreasing = TRUE)
largest.part <- as.numeric(names(part.tab)[1])
cat("  最大 partition:", largest.part,
    " (", part.tab[1], "细胞 / ", length(partitions), "总细胞 )\n")
cells.in.largest <- which(partitions == largest.part)
root.pr.vertices <- as.vector(closest_vertex[cells.in.largest, 1])
chosen.node <- as.numeric(names(which.max(table(root.pr.vertices))))
root.pr.nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[chosen.node]

# ---- 断点 2：伪时间排序结果 ----
# order_cells 是降维后非线性映射，cds 序列化后近邻索引不保留，
# 因此 ckpt4 必须在 order_cells 计算完成后才保存，且 resume 时不重跑
ckpt4 <- file.path(TAB_DIR, "ckpt4_ordered.rds")
if (file.exists(ckpt4)) {
  cat("  [Checkpoint 2] 发现已有伪时间排序，载入:", basename(ckpt4), "\n")
  cds <- readRDS(ckpt4)
} else {
  cds <- order_cells(cds, root_pr_nodes = root.pr.nodes)
  saveRDS(cds, ckpt4)
}

p.pseudotime <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  graph_label_size = 3
) + labs(title = "Pseudotime") +
  scale_color_viridis_c() + sci_theme

sci_save(p.pseudotime, "Pseudotime", width = 10, height = 8)

# ---- 4C: 完整 2x2 多面板轨迹图（含 pseudotime）-----------------
# order_cells 完成后 colData(cds)$pseudotime 已赋值，安全调用

p.trajectory.celltype2 <- plot_cells(
  cds, color_cells_by = "celltype",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  graph_label_size = 2
) + .theme.small

p.trajectory.tissue2 <- plot_cells(
  cds, color_cells_by = "tissue",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  graph_label_size = 2
) + .theme.small

p.trajectory.fbxl12.2 <- plot_cells(
  cds, color_cells_by = "FBXL12_status",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  graph_label_size = 2
) + .theme.small

p.trajectory.pseudotime2 <- plot_cells(
  cds, color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE,
  graph_label_size = 2
) + scale_color_viridis_c() + .theme.small

p.multi.4c <- (p.trajectory.celltype2 | p.trajectory.tissue2) /
               (p.trajectory.fbxl12.2 | p.trajectory.pseudotime2) +
               plot_annotation(tag_levels = "A") &
               .theme.tag

sci_save(p.multi.4c, "Trajectory_MultiPanel_Full", width = 14, height = 12)

# ---- 4B: 轨迹统计表（P0-5）-----------------------------------
# 提取 principal graph 的节点数、边数、各分支细胞数


pg <- principal_graph(cds)[["UMAP"]]

n.nodes <- igraph::vcount(pg)
n.edges <- igraph::ecount(pg)

# 提取各分支（partition）
partitions <- cds@clusters$UMAP$partitions
n.partitions <- length(unique(partitions[!is.na(partitions)]))

# 各分支细胞数
partition_counts <- table(partitions[!is.na(partitions)])
names(partition_counts) <- paste0("Partition_", names(partition_counts))

# 主图节点名称
node.names <- igraph::V(pg)$name
n.pg.nodes <- length(node.names)

# 分支节点（度数 > 2 的节点 = 分叉点）
node.degrees <- igraph::degree(pg)
branch.nodes <- names(node.degrees[node.degrees > 2])
n.branch.points <- length(branch.nodes)

# 叶节点（度数 = 1 的节点）
leaf.nodes <- names(node.degrees[node.degrees == 1])
n.leaves <- length(leaf.nodes)

stats.list <- list(
  Metric = c(
    "Number of Principal Graph Nodes",
    "Number of Principal Graph Edges",
    "Number of Partitions (Branches)",
    "Number of Branch Points (Degree>2)",
    "Number of Leaf Nodes (Degree=1)",
    "Root Graph Node (ID)",
    "Total Cells"
  ),
  Value = c(
    n.pg.nodes,
    n.edges,
    n.partitions,
    n.branch.points,
    n.leaves,
    root.pr.nodes,
    ncol(cds)
  )
)

traj.stats <- as.data.frame(stats.list)
write.csv(traj.stats,
  file = file.path(TAB_DIR, "Table3_Trajectory_Stats.csv"),
  row.names = FALSE)

# ================================================================
# PART 5: 伪时间依赖的差异表达基因
# 目的：鉴定沿伪时间显著变化的基因
# 方法：graph_test(q_value < 0.05)
# 数据：cds → table/Table1_pseudotime_DEG.csv
# ================================================================

# ---- 断点 3：graph_test 结果 ----
# graph_test 是第二耗时的步骤（~2h），中断后从此续跑
ckpt5 <- file.path(TAB_DIR, "ckpt5_graph_test.rds")
if (file.exists(ckpt5)) {
  cat("  [Checkpoint 3] 发现已有 graph_test 结果，载入:", basename(ckpt5), "\n")
  res <- readRDS(ckpt5)
  cds.graph.test <- res$cds.graph.test
  cds.branch.test <- res$cds.branch.test
  deg.pseudo      <- res$deg.pseudo
  deg.branch      <- res$deg.branch
} else {
  cat("========== PART 5: 伪时间DEG（子抽样3000细胞）==========\n")

  # 子抽样 3000 细胞（平衡覆盖与速度）
  # 仅对有限伪时间细胞抽样；全部 32995 细胞仍用于伪时间排序和热图
  pt <- pseudotime(cds)
  finite.cells <- names(pt)[is.finite(pt)]
  set.seed(42)
  if (length(finite.cells) > 3000) {
    sampled.cells <- sample(finite.cells, 3000)
  } else {
    sampled.cells <- finite.cells
  }
  cat("  有限伪时间细胞:", length(finite.cells),
      " → 抽样:", length(sampled.cells), "\n")
  cds.sub <- cds[, sampled.cells]

  cat("  [Graph Test 1/2] 伪时间依赖基因（", length(sampled.cells), "细胞）...\n")
  cds.graph.test <- graph_test(cds.sub,
                               neighbor_graph = "principal_graph",
                               cores = 4)
  deg.pseudo <- cds.graph.test %>%
    filter(q_value < 0.05) %>%
    arrange(q_value)
  write.csv(deg.pseudo,
    file = file.path(TAB_DIR, "Table1_pseudotime_DEG.csv"),
    row.names = FALSE)

  # ---- 5B: 分支特异性差异基因 -------------------------
  cat("========== PART 5B: 分支特异性DEG ==========\n")
  cat("  [Graph Test 2/2] 分支特异性基因（", length(sampled.cells), "细胞）...\n")
  cds.branch.test <- graph_test(cds.sub,
                                 neighbor_graph = "principal_graph",
                                 alternative = "branching",
                                 cores = 4)
  deg.branch <- cds.branch.test %>%
    filter(q_value < 0.05) %>%
    arrange(q_value)
  write.csv(deg.branch,
    file = file.path(TAB_DIR, "Table2_Branch_DEG.csv"),
    row.names = FALSE)

  saveRDS(list(
    cds.graph.test = cds.graph.test,
    cds.branch.test = cds.branch.test,
    deg.pseudo      = deg.pseudo,
    deg.branch      = deg.branch
  ), ckpt5)
}

# ================================================================
# PART 6: 伪时间 DEG 热图（含列注释，P0-2）
# 目的：展示 Top 30 沿伪时间变化的基因 + 细胞类型/组织/伪时间注释
# 方法：ComplexHeatmap，行 Z-score，top_annotation 标注元数据
# 数据：cds + deg.pseudo → Figure/Pseudotime_DEG_Heatmap.pdf
# ================================================================


top.genes <- head(deg.pseudo$gene, 30)

# 沿伪时间排序细胞
pt.full <- pseudotime(cds)
finite.cells <- is.finite(pt.full)
cat("  有限伪时间细胞:", sum(finite.cells), "/", length(pt.full), "\n")

pt.order <- order(pt.full)
cells.ordered <- colnames(cds)[pt.order]

# 仅保留有限伪时间细胞（Inf 细胞在另一 partition，无伪时间距离）
cells.ordered <- intersect(cells.ordered, colnames(cds)[finite.cells])
# 位置索引：不依赖 colData 列是否带 names（factor/字符按名索引会返回 NA）
pos.ordered <- match(cells.ordered, colnames(cds))
pt.ordered <- pt.full[pos.ordered]

# 提取表达矩阵并缩放
# 过滤零方差的基因（scale 遇到 sd=0 会返回 NaN，导致 hclust 报错）
expr.mat <- exprs(cds)[top.genes, cells.ordered, drop = FALSE]
expr.mat <- as.matrix(expr.mat)
gene.sd <- apply(expr.mat, 1, sd, na.rm = TRUE)
if (any(gene.sd == 0 | is.na(gene.sd))) {
  keep <- which(gene.sd > 0 & !is.na(gene.sd))
  cat("  排除零方差基因:", length(top.genes) - length(keep), "\n")
  top.genes <- top.genes[keep]
  expr.mat <- expr.mat[keep, , drop = FALSE]
}
expr.scaled <- t(scale(t(expr.mat)))

# ---- 构建列注释元数据 ----
# 细胞类型
celltype.ordered <- colData(cds)$celltype[pos.ordered]
# 组织来源
tissue.ordered <- colData(cds)$tissue[pos.ordered]
# 细胞类型颜色映射：使用高饱和、高区分度的鲜艳配色
# 注意：celltype.ordered 是 factor，必须用 as.character 后再按名称索引，
# 否则 factor 会被当作整数位置索引，导致颜色完全错位；同时必须覆盖所有 level。
CELLTYPE_COLS <- c(
  "Epithelial/Tumor" = "#E41A1C",
  "Mono/Macrophage"  = "#377EB8",
  "B cells"          = "#4DAF4A",
  "T cells"          = "#984EA3",
  "Plasma cells"     = "#FF7F00",
  "Mast cells"       = "#F781BF",
  "pDC"              = "#A65628",
  "NK cells"         = "#999999",
  "CAFs"             = "#66C2A5",
  "Endothelial"      = "#FFD92F",
  "DC"               = "#8DA0CB"
)
celltype.cols <- CELLTYPE_COLS[unique(as.character(celltype.ordered))]

# 组织颜色映射
TISSUE_COLS <- c(PT = "#E69F00", LM = "#56B4E9")

# 列注释对象（显式传入颜色映射，避免 ComplexHeatmap 自动分配高相似色）
ha.top <- HeatmapAnnotation(
  Celltype = celltype.ordered,
  Tissue   = tissue.ordered,
  Pseudotime = anno_points(pt.ordered,
                           gp = gpar(col = "#4DAF4A"),
                           size = unit(1, "pt"),
                           axis = TRUE),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 8),
  col = list(Celltype = celltype.cols, Tissue = TISSUE_COLS)
)

# 热图（pseudotime 为列级信息，不在行注释中重复展示）
hm <- Heatmap(expr.scaled,
  name = "Z-score",
  top_annotation = ha.top,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  column_title = "Top 30 Pseudotime DEG (ordered by pseudotime)",
  use_raster = TRUE,
  raster_quality = 3
)
sci_save(hm, "Pseudotime_DEG_Heatmap", width = 14, height = 10)

# ---- 6B: 分支 DEG 热图（如果存在显著分支基因）----------------

if (nrow(deg.branch) > 0) {
  top.branch.genes <- head(deg.branch$gene, 30)
  expr.branch <- exprs(cds)[top.branch.genes, cells.ordered, drop = FALSE]
  expr.branch <- as.matrix(expr.branch)
  bg.sd <- apply(expr.branch, 1, sd, na.rm = TRUE)
  if (any(bg.sd == 0 | is.na(bg.sd))) {
    keep.bg <- which(bg.sd > 0 & !is.na(bg.sd))
    cat("  排除零方差分支基因:", length(top.branch.genes) - length(keep.bg), "\n")
    top.branch.genes <- top.branch.genes[keep.bg]
    expr.branch <- expr.branch[keep.bg, , drop = FALSE]
  }
  expr.branch.scaled <- t(scale(t(expr.branch)))

  hm.branch <- Heatmap(expr.branch.scaled,
    name = "Z-score",
    top_annotation = ha.top,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 7),
    column_title = "Top 30 Branch-dependent DEG",
    use_raster = TRUE,
    raster_quality = 3
  )
  sci_save(hm.branch, "Branch_DEG_Heatmap", width = 14, height = 10)
}

# ================================================================
# PART 7: FBXL12 伪时间趋势线图（P0-4）
# 目的：展示 FBXL12 沿伪时间的表达量变化趋势（含误差棒）
# 方法：伪时间分箱 → dplyr 计算均值 ± SE → ggplot2 折线图
# 数据：cds → Figure/FBXL12_Expression_Trend.pdf
# ================================================================


# 提取 FBXL12 在各细胞的表达量（log-transformed counts）
fbxl12.expr <- log1p(exprs(cds)["FBXL12", ])
pt.vals <- pseudotime(cds)

# 构造数据框（排除 NA）
df.fbxl12 <- data.frame(
  expression = as.numeric(fbxl12.expr),
  pseudotime = as.numeric(pt.vals),
  celltype   = colData(cds)$celltype,
  tissue     = colData(cds)$tissue,
  fbxl12.st  = colData(cds)$FBXL12_status,
  stringsAsFactors = FALSE
) %>% filter(is.finite(pseudotime))

# 分箱（20 个 bin）
df.fbxl12 <- df.fbxl12 %>%
  mutate(pt_bin = cut(pseudotime, breaks = 20, labels = FALSE, include.lowest = TRUE)) %>%
  group_by(pt_bin) %>%
  mutate(
    pt_mid   = mean(pseudotime),
    pt_mean  = mean(expression, na.rm = TRUE),
    pt_se    = sd(expression, na.rm = TRUE) / sqrt(n()),
    n_cells  = n()
  ) %>%
  ungroup()

# 全局趋势线
p.trend <- ggplot(df.fbxl12, aes(x = pseudotime, y = expression)) +
  geom_point(aes(color = fbxl12.st), alpha = 0.15, size = 0.8) +
  geom_smooth(method = "loess", se = TRUE,
              color = "#E63946", fill = "#E63946", alpha = 0.2) +
  labs(
    title = "FBXL12 Expression along Pseudotime",
    x = "Pseudotime", y = "log(FBXL12 + 1)",
    color = "FBXL12 Status"
  ) +
  sci_theme

sci_save(p.trend, "FBXL12_Expression_Trend", width = 9, height = 6)

# 分组织趋势（叠加折线 ± SE）
p.trend.tissue <- df.fbxl12 %>%
  group_by(pt_bin, tissue) %>%
  summarise(
    pt_mid = dplyr::first(pt_mid),
    pt_mean = mean(expression, na.rm = TRUE),
    pt_se = sd(expression, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = pt_mid, y = pt_mean, color = tissue)) +
  geom_line(size = 1.2) +
  geom_errorbar(aes(ymin = pt_mean - pt_se, ymax = pt_mean + pt_se),
                 width = 0.05, alpha = 0.7) +
  geom_point(size = 1.5) +
  labs(
    title = "FBXL12 Expression Trend by Tissue",
    x = "Pseudotime (binned, mean +/- SE)",
    y = "log(FBXL12 + 1)",
    color = "Tissue"
  ) +
  scale_color_manual(values = TISSUE_COLS) +
  sci_theme

sci_save(p.trend.tissue, "FBXL12_Expression_Trend_Tissue", width = 9, height = 6)

# 空间分布图（备选）
p.fbxl12.space <- plot_cells(
  cds, genes = "FBXL12",
  label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE
) + labs(title = "FBXL12 Spatial Distribution") +
  scale_color_viridis_c() + sci_theme

sci_save(p.fbxl12.space, "FBXL12_Spatial", width = 10, height = 8)

# ================================================================
# PART 8: 伪时间构成比变化图（P1-6）
# 目的：展示各细胞类型 / 组织在伪时间各分箱中的比例变化
# 方法：伪时间分箱 → dplyr 计数 → ggplot2 stacked barplot
# 数据：cds → Figure/Pseudotime_Composition.pdf
# ================================================================


# 细胞类型构成比
df.comp.ct <- df.fbxl12 %>%
  group_by(pt_bin, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(pt_bin) %>%
  mutate(freq = n / sum(n) * 100)

p.comp.celltype <- df.comp.ct %>%
  ggplot(aes(x = factor(pt_bin), y = freq, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.85) +
  labs(
    title = "Cell Type Composition across Pseudotime Bins",
    x = "Pseudotime Bin (1=early, 20=late)",
    y = "Cell Proportion (%)",
    fill = "Cell Type"
  ) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

sci_save(p.comp.celltype, "Pseudotime_Composition_Celltype", width = 12, height = 6)

# 组织来源构成比
df.comp.tissue <- df.fbxl12 %>%
  group_by(pt_bin, tissue) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(pt_bin) %>%
  mutate(freq = n / sum(n) * 100)

p.comp.tissue <- df.comp.tissue %>%
  ggplot(aes(x = factor(pt_bin), y = freq, fill = tissue)) +
  geom_bar(stat = "identity", width = 0.85) +
  labs(
    title = "Tissue Composition across Pseudotime Bins",
    x = "Pseudotime Bin (1=early, 20=late)",
    y = "Cell Proportion (%)",
    fill = "Tissue"
  ) +
  scale_fill_manual(values = TISSUE_COLS) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

sci_save(p.comp.tissue, "Pseudotime_Composition_Tissue", width = 12, height = 6)

# ================================================================
# PART 9: Top DEG Feature Plot（P1-7）
# 目的：单独展示 Top 6 伪时间 DEG 的 UMAP 空间分布
# 方法：Seurat::FeaturePlot
# 数据：seurat.obj → Figure/TopDEG_FeaturePlot.pdf（2x3 多面板）
# ================================================================


top6.genes <- head(deg.pseudo$gene, 6)

# 过滤存在于 seurat.obj 中的基因
top6.avail <- top6.genes[top6.genes %in% rownames(seurat.obj)]

if (length(top6.avail) >= 4) {
  n.rows <- ceiling(length(top6.avail) / 3)

  # Feature plot 列表
  fp.list <- lapply(top6.avail, function(g) {
    FeaturePlot(seurat.obj, features = g, reduction = "umap",
                cols = c("grey90", "#E63946")) &
      theme(plot.title = element_text(size = 11, hjust = 0.5, face = "plain"),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank())
  })

  if (length(fp.list) >= 4) {
    p.fp <- wrap_plots(fp.list, ncol = 3) +
      plot_annotation(title = "Top Pseudotime DEG - Spatial Distribution",
                     tag_levels = "A") &
      theme(plot.tag = element_text(size = 12, face = "bold"))
    sci_save(p.fp, "TopDEG_FeaturePlot", width = 12, height = 4 * n.rows)
  }
}

# ================================================================
# PART 10: 保存结果
# 目的：保存带伪时间信息的 cell_data_set 对象
# 方法：saveRDS
# 数据：cds → table/cds_trajectory.rds
# ================================================================

saveRDS(cds, file = file.path(TAB_DIR, "cds_trajectory.rds"))


