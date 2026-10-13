library(monocle3)
library(ComplexHeatmap)
library(circlize)

# 常量
DIR_07 <- "X:/CRLM/XX/results/07_Trajectory"
FIG_DIR <- file.path(DIR_07, "Figure")
TAB_DIR <- file.path(DIR_07, "table")

sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  pdf(pdf_file, width = width, height = height)
  draw(plot_obj)
  dev.off()
  cat("  已保存:", pdf_file, "\n")
}

cds <- readRDS(file.path(TAB_DIR, "cds_trajectory.rds"))
deg.pseudo <- read.csv(file.path(TAB_DIR, "Table1_pseudotime_DEG.csv"), stringsAsFactors = FALSE)

top.genes <- head(deg.pseudo$gene, 30)

pt.full <- pseudotime(cds)
finite.cells <- is.finite(pt.full)
pt.order <- order(pt.full)
cells.ordered <- colnames(cds)[pt.order]
cells.ordered <- intersect(cells.ordered, colnames(cds)[finite.cells])
pos.ordered <- match(cells.ordered, colnames(cds))
pt.ordered <- pt.full[pos.ordered]

expr.mat <- exprs(cds)[top.genes, cells.ordered, drop = FALSE]
expr.mat <- as.matrix(expr.mat)
gene.sd <- apply(expr.mat, 1, sd, na.rm = TRUE)
keep <- which(gene.sd > 0 & !is.na(gene.sd))
if (length(keep) < length(top.genes)) {
  cat("  排除零方差基因:", length(top.genes) - length(keep), "\n")
  top.genes <- top.genes[keep]
  expr.mat <- expr.mat[keep, , drop = FALSE]
}
expr.scaled <- t(scale(t(expr.mat)))

celltype.ordered <- colData(cds)$celltype[pos.ordered]
tissue.ordered <- colData(cds)$tissue[pos.ordered]

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
TISSUE_COLS <- c(PT = "#E69F00", LM = "#56B4E9")

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
