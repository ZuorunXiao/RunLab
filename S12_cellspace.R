#===============================================================
## S12_cellspace.R —— 细胞层与空间层（步骤规格 §8：S12）
##
## 目的：
##   (a) scRNA：Seurat + harmony 跨供者整合（GSE237230 AAA / GSE247238 斑块 /
##       GSE212388 人主动脉 EC 血流 ± 抗增殖药）；canonical marker 簇级注释 +
##       五段模块打分 x 细胞类型（方案 L156「免疫亚群的五段活性」，喂 S13 细胞层）；
##       GSE212388 伪批量 + limma 出抗增殖药五段方向（方案 L133 Rapa 证据臂，
##       与 S11 同口径）。
##   (b) 空间：GeoMx DSP 的 dcc.gz 做区域级段分。
##       §0.3 实测：GSE277170 分 Panel1 / Panel2，须分别 QC 并通过段基因数断言。
##
## 定位：分析脚本。本步只出「亚群 / 区域」层表，不合并到层 1-3 的结论（那是 S13）。
##
## 产出：03_results/S12/table/scRNA_clusters.csv       每簇一行（簇 × 细胞数）
##       03_results/S12/table/scRNA_markers.csv       每簇正标记基因（聚类 QC 依据）
##       03_results/S12/table/scRNA_celltype_annotation.csv  簇级细胞类型注释 [74]
##       03_results/S12/table/segment_by_celltype.csv 五段活性 × 细胞类型 [74]
##       03_results/S12/table/gse212388_rapa_direction.csv   抗增殖药五段方向 [74]
##       03_results/S12/table/spatial_panel_qc.csv
##       03_results/S12/table/spatial_roi_qc.csv
##       03_results/S12/Figure/ScRNA_UMAP_{acc}_CellType.pdf    注释 UMAP（图内标题带 acc）
##       03_results/S12/Figure/Segment_By_CellType_Dot_{acc}.pdf 段活性点图 [74]
##       03_results/S12/Figure/EC_Segment_Direction_Forest.pdf   Rapa 方向点图 [74/74d]
##       03_results/S12/Figure/EC_SegmentScore_Forest.pdf        样本级段分 log2FC 森林图 [75]
##       03_results/S12/Figure/EC_Cycling_Fraction_By_Treatment.pdf 增殖比例 [75]
##       03_results/S12/Figure/ScRNA_Cluster_Size_By_Series.pdf   分面：系列 x 簇细胞数
##       03_results/S12/table/gse212388_segment_score_stats.csv   样本级段分 limma [75]
##       03_results/S12/table/gse212388_segment_score_by_sample.csv 每样本段分 [75]
##       03_results/S12/table/gse212388_cycling_fraction.csv     增殖亚群比例 [75]
##       03_results/S12/table/segment_by_celltype_stats.csv      段分×细胞类型 KW [75]
##   （空间分支无图：GSE277170 无 ROI 区域元数据，见 与 Limitations）
##
## 图件清单：
##   可视化: S12 共 7 类（空间 2 张 QC 图已按）
##     1. ScRNA_Cluster_Size_By_Series.pdf
##        —— 形式：分面柱（x = 簇编号，y = 细胞数，按系列分面）
##        —— 科学问题：整合后每个系列的簇大小分布（仅保留 >= MIN_CELLS_PER_CLUSTER 的簇）。
##        —— 参照文献：Seurat 单细胞整合与聚类（Hao 2021, Cell,
##           doi:10.1016/j.cell.2021.04.048, PMID 34062119）。
##     2. ScRNA_UMAP_{acc}_CellType.pdf  (每系列各一张)
##        —— 形式：UMAP 按注释细胞类型着色 + label=TRUE 标类型名。
##        —— 科学问题：各系列的亚群构成（方案 L546「亚群表」）；GSE212388 应以
##           内皮为主（单一细胞类型扰动数据，簇间连续谱是数据本质非聚类缺陷）。
##        —— 参照文献：AddModuleScore（Tirosh 2016, Science,
##           doi:10.1126/science.aad0501, PMID 26849431）+ Seurat 可视化标准。
##     3. Segment_By_CellType_Dot_{acc}.pdf  (每系列各一张)
##        —— 形式：点图（x = 五段，y = 细胞类型，size/color = 平均模块分）。
##        —— 科学问题：五段炎症活性落在哪类细胞（方案 L144 核心问题之「哪类细胞」）；
##           显著性见 segment_by_celltype_stats.csv（Kruskal-Wallis + BH）。
##        —— 参照文献：模块打分点图惯例（Tirosh 2016, PMID 26849431）。
##     4. EC_Segment_Direction_Forest.pdf
##        —— 形式：方向点图（y = 段，x = 段基因平均 log2FC，色 = 药臂）。
##        —— 科学问题：Rapa/Paclitaxel 在 EC 中对五段的位移方向（方案 L133
##           Rapa 证据臂）；BH 后无显著基因时为描述性点估计（indeterminate
##           是合法输出，方案 L15），显著度读数在 CSV 的 n.sig/min.BH/Wilson CI。
##        —— 参照文献：方案 L15 三层递进设计；Wilson CI 记 Newcombe 1998
##           (PMID 9699278)（CSV 口径）。
##     5. EC_SegmentScore_Forest.pdf
##        —— 形式：森林图（y = 段，x = 样本级段分 log2FC + 95% CI，色 = 药臂）。
##        —— 科学问题：Rapa/Paclitaxel 是否位移 EC 的五段活性（样本为重复单位，
##           规避基因级 BH 在 n=2 vs 4 下零检验力的缺陷）。
##        —— 参照文献：limma  moderated t（Smyth 2004, Stat Appl Genet Mol Biol,
##           doi:10.2202/1544-6115.1026, PMID 16646809）。
##     6. EC_Segment_Effect_By_Flow.pdf，替代旧
##        EC_SegmentScore_By_FlowCondition 绝对水平点图，旧图归档删除）
##        —— 形式：patchwork 双联效应图。A = 各 flow 下药效 vs 匹配 Control
##           （limma 95%CI，段名短码 S1-S5，PhysioFlow 在前）；B = 交互
##           （药效在 Altered - Physio 之差）；Rapa 红 / Pac 蓝，0 线灰虚线。
##        —— 科学问题：药物效应是否依赖炎症背景（方案 L144「哪种炎症背景下」）——
##           GSE212388 的设计即生理流 / 紊乱流 x Control/Rapa/Paclitaxel。
##        —— 参照文献：紊乱流致炎模型（Davies 2013 综述框架）；交互项仅探索性
##           （残差 df = 2），显著性读数落 gse212388_segment_flow_interaction.csv。
##     7. EC_Cycling_Fraction_By_Treatment.pdf
##        —— 形式：抖动点图（x = 处理，y = 增殖亚群占比，每点一样本）。
##        —— 科学问题：抗增殖药（Rapa/Paclitaxel）是否压低 EC 增殖亚群比例——
##           Rapa 抗增殖效应的单细胞直接读数（方案 L133）。
##        —— 参照文献：细胞周期标记（Tirosh 2016, PMID 26849431）与增殖亚群
##           比例比较惯例（同项目 S11 扰动臂方向仲裁口径）。
##
## 图件口径：sci_theme / sci_save 逐字取自规范 §sci_theme / §sci_save 官方定义；
##       分面一律 facet_wrap / facet_grid，不引入任何拼图包；可读性核验见
##       sci_save 之 RAPA_QC_PNG 钩子。
##
## PART 1（scRNA）说明：harmony 整合 + 无监督聚类（QC 层），再对齐主线——
##       canonical marker 簇级注释（含 Cycling 两级命名）、五段模块打分 x 细胞类型
##       （+ Kruskal-Wallis 统计）、GSE212388 伪批量 limma 基因级方向、
##       样本级段分 limma 与增殖亚群比例。分样本图经用户裁定删除（不必要保留）。
## PART 2（空间）说明：**只做 Panel 覆盖度 QC**（Q11：panel 是否覆盖五段基因集），
##       判定落 spatial_panel_qc.csv 与收口 message，**不产图**。
##       本步刻意**不做「区域级段分打分」**：实测 GSE277170 无 ROI 级区域元数据
##       （series_matrix 仅 1 个 NTC 样本，tissue/roi/subset/grade/localisation
##       全 NA），无区域标签可解释；且跨 ROI 归一需有出处的空间评分方法。
##       宁可少一个表，也不自造公式——该缺口写入 Limitations。
##
## 前置：S02（基因集）。空间分支解析 pkc 需 jsonlite；scRNA 分支需 Seurat；
##       两者缺失即显式 stop，不降级。无需 GeoMxTools（pkc 用 jsonlite 直接解析）。
##       签名依据见方案 §0.2-measured，取不到预期槽位即 stop，不猜槽位继续。
#===============================================================

## PART 0 常量 ------------------------------------------------------------------------------

SEED <- 20260925
set.seed(SEED)
FIG_DIR <- "03_results/S12/Figure"
TAB_DIR <- "03_results/S12/table"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
## 断点落 checkpoints/（命名 S12_<acc>_seurat.rds）；TAB_DIR 须先落定，
## 否则 file.path 在其赋值前求值而立即报错
CKPT_DIR <- file.path(TAB_DIR, "checkpoints")
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)

library(ggplot2)   ## 簇与 panel 图；主题出口用规范 §sci_theme / §sci_save
library(patchwork) ## 效应图双联拼接（A 效应 + B 交互）
## 空间分支解析 NanoString pkc（RTS_ID -> 基因符号），替代缺失的 GeoMxTools
library(jsonlite)

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

MIN_CELLS_PER_CLUSTER <- 10L          ## 少于此细胞数的簇不参与下游
MIN_REGION_GENES <- 20L               ## §9 Q11：Panel 命中的段基因数下限
MIN_ROI_RTS <- 10000L
vec.segment <- c("S1_initiation", "S2_amplification", "S3_effect",
  "S4_resolution", "S5_chronic")

dat.final <- read.csv("03_results/S02/table/inflammation_geneset_final.csv",
  stringsAsFactors = FALSE)
vec.col.req <- c("segment", "gs_name", "members")
vec.col.miss <- setdiff(vec.col.req, names(dat.final))
if (length(vec.col.miss) > 0L) {
  stop("S12 FAILED: inflammation_geneset_final.csv lacks column(s): ",
    paste(vec.col.miss, collapse = ", "))
}
## 每段基因向量：S02 已用 msigdbr 实测成员，空成员串表示未命中当前 MSigDB 版本
list.seg.gene <- lapply(vec.segment, function(s) {
    dat.s <- dat.final[dat.final$segment == s & dat.final$members != "", ]
    unique(unlist(lapply(dat.s$members, function(x) strsplit(x, ";")[[1L]])))
  })
names(list.seg.gene) <- vec.segment
## Q3 一致性：五段都必须在锁定表中出现，否则下游有段会静默缺席
vec.seg.empty <- vec.segment[lengths(list.seg.gene) == 0L]
if (length(vec.seg.empty) > 0L) {
  stop(paste0("S12 FAILED: segment(s) have no genes after MSigDB matching: ",
      paste(vec.seg.empty, collapse = ", "),
      " (see 03_results/S02/table/segment_member_count.csv)."))
}

ACC_RAPA <- "GSE212388"
##   S11 同口径：PVAL_CUTOFF / LOGFC_MIN 与 summarise_direction 三分支判定
##   （S11_perturbation.R L102-103 / L148-161），禁自造阈值。
PVAL_CUTOFF <- 0.05
LOGFC_MIN <- 0.585
##   细胞类型 canonical marker 注释参考（零外部依赖、离线可审计；
##   血管/免疫标准标记，来自 CellMarker 与血管生物学通识；簇级 argmax）。
MIN_MARKER_CLUSTER_MEAN <- 0.25
lst.ct.mk <- list(
  Endothelial = c("PECAM1", "VWF", "KDR", "CDH5", "CLDN5", "NPR3", "EMCN"),
  SMC_Pericyte = c("ACTA2", "TAGLN", "MYH11", "CNN1", "RGS5", "LMOD1"),
  Fibroblast = c("COL1A1", "DCN", "LUM", "COL3A1"),
  Macrophage_Mono = c("LYZ", "CD68", "C1QA", "C1QB", "APOC1", "FCN1",
    "SPP1"),
  Neutrophil = c("CXCR2", "FCGR3B", "CSF3R", "S100A8", "S100A9"),
  T_cell = c("CD3D", "CD3E", "CD3G", "TRBC2", "IL32"),
  NK_cytotoxic = c("NKG7", "GNLY", "GZMB", "PRF1"),
  pDC = c("LILRA4", "PLD4", "SERPINF1", "LILRB2"),
  B_cell = c("MS4A1", "CD79A", "BANK1", "CD79B"),
  Plasma = c("IGLL5", "MZB1", "JCHAIN", "IGHG1"),
  Mast = c("HDC", "CPA3", "KIT"),
  Cycling = c("MKI67", "TOP2A", "CCNB1", "PCNA", "ASPM", "CEP55", "TYMS",
    "UBE2C"))

## PART 1 scRNA：Seurat + harmony ----------------------------------------------------------------

vec.sc <- c("GSE237230", "GSE247238", "GSE212388")
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop(
    "S12 FAILED: 'Seurat' not installed ; scRNA branch aborted (no fallback).")
}
library(Seurat)                           ## 空间 / 单细胞对象构建、聚类与整合
library(harmony)  ## Harmony 整合（可选；包缺失即 stop，不静默降级）

dat.cluster <- data.frame(acc = character(0), cluster = integer(0),
  n_cells = integer(0),
  stringsAsFactors = FALSE)

lst.markers <- list()   ## 每 acc 的簇标记表，循环后合并写盘
dat.ann <- data.frame(acc = character(0), cluster = integer(0),
  celltype = character(0), n_cells = integer(0), stringsAsFactors = FALSE)
dat.seg.ct <- data.frame(acc = character(0), celltype = character(0),
  segment = character(0), mean_score = numeric(0),
  median_score = numeric(0), sd_score = numeric(0), n_cells = integer(0),
  stringsAsFactors = FALSE)
dat.seg.sample <- data.frame(acc = character(0), sample = character(0),
  segment = character(0), mean_score = numeric(0), n_cells = integer(0),
  stringsAsFactors = FALSE)
dat.cycl <- data.frame(acc = character(0), sample = character(0),
  n_cells = integer(0), n_cycling = integer(0), frac_cycling = numeric(0),
  stringsAsFactors = FALSE)
dat.ct.stat <- data.frame(acc = character(0), segment = character(0),
  kw_p = numeric(0), top_celltype = character(0), max_mean = numeric(0),
  min_mean = numeric(0), n_celltypes = integer(0), stringsAsFactors = FALSE)
dat.dir <- data.frame(acc = character(0), arm = character(0),
  segment = character(0), n_genes_tested = integer(0), n.sig = integer(0),
  n.up = integer(0), n.down = integer(0), dir.index = numeric(0),
  direction = character(0), mean.logFC.seg = numeric(0),
  median.logFC.seg = numeric(0), min.BH = numeric(0),
  stringsAsFactors = FALSE)
for (acc in vec.sc) {
  cat(sprintf("\n[S12] scRNA %s (Seurat + harmony)\n", acc))
  dir.acc <- file.path("01_data", acc)
  vec.sdir <- list.dirs(dir.acc, recursive = FALSE)
  if (length(vec.sdir) == 0L) {
    stop(paste0("S12 STOP: no per-sample 10X dirs for ", acc,
        "; run S12a first."))
  }
  ## 整条 Seurat 处理链是重时步骤，按 acc 粒度续传；
  ## 恢复失败即报错，不静默丢弃损坏对象。
  f.ck <- file.path(CKPT_DIR, paste0("S12_", acc, "_seurat.rds"))
  if (file.exists(f.ck)) {
    obj <- readRDS(f.ck)
    obj <- SeuratObject::JoinLayers(obj)
    red.use <- if ("harmony" %in% names(obj@reductions)) "harmony" else "pca"
  } else {
    ## 每样本读 10X 目录（barcodes/features/matrix），合并前各自过最低门槛
    lst.obj <- lapply(vec.sdir, function(d) {
        counts <- Seurat::Read10X(d)
        Seurat::CreateSeuratObject(counts = counts, project = basename(d),
          min.cells = 3L, min.features = 200L)
      })
    obj <- Reduce(function(a, b) merge(a, y = b), lst.obj)
    obj <- Seurat::NormalizeData(obj)
    obj <- Seurat::FindVariableFeatures(obj, n.top.genes = 2000)
    obj <- Seurat::ScaleData(obj)
    obj <- Seurat::RunPCA(obj, npcs = 30)
    ## 跨供者整合：FindNeighbors 必须在 PCA（及 harmony）之后，否则 'pca' not found；
    ## harmony 可用则用 harmony 归约建邻居图（真正落地整合），否则退回 pca 归约。
    if (requireNamespace("harmony", quietly = TRUE)) {
      obj <- harmony::RunHarmony(obj, group.by.vars = "orig.ident")
      red.use <- "harmony"
    } else {
      red.use <- "pca"
    }
    obj <- Seurat::FindNeighbors(obj, reduction = red.use, dims = 1:20)
    obj <- Seurat::FindClusters(obj, resolution = 0.5)
    obj <- Seurat::RunUMAP(obj, reduction = red.use, dims = 1:20)
    saveRDS(obj, f.ck)
  }
  ## UMAP 底图（harmony/pca 归约；Stuart 2019 PMID 31178118 / Hao 2021 PMID 34062119）。
  if (is.null(obj@reductions$umap)) {
    obj <- Seurat::RunUMAP(obj, reduction = red.use, dims = 1:20)
  }
  vec.ct <- names(lst.ct.mk)
  vec.cl <- as.character(obj$seurat_clusters)
  mat.mk.cl <- matrix(0, nrow = length(vec.ct), ncol = length(unique(vec.cl)),
    dimnames = list(vec.ct, sort(unique(vec.cl))))
  for (ct in vec.ct) {
    g.ct <- intersect(lst.ct.mk[[ct]], rownames(obj))
    if (length(g.ct) >= 2L) {
      ## 只取标记基因行再转密集（整矩阵密集会爆内存：23754 x 23614）
      mat.sub <- as.matrix(SeuratObject::LayerData(obj[["RNA"]],
        layer = "data")[g.ct, , drop = FALSE])
      mat.mk.cl[ct, ] <- colMeans(sapply(split(seq_len(ncol(mat.sub)),
        vec.cl), function(idx) rowMeans(mat.sub[, idx, drop = FALSE])))
    }
    max.cl <- max(mat.mk.cl[ct, ])
    if (length(g.ct) >= 2L && max.cl >= MIN_MARKER_CLUSTER_MEAN) {
      set.seed(SEED)
      obj <- Seurat::AddModuleScore(obj, features = list(g.ct),
        name = paste0("ctscore_", ct, "_"))
    } else {
      obj@meta.data[[paste0("ctscore_", ct, "_1")]] <-
        rep(NA_real_, ncol(obj))
    }
  }
  mat.ct.sc <- as.matrix(obj@meta.data[, paste0("ctscore_", vec.ct, "_1"),
    drop = FALSE])
  colnames(mat.ct.sc) <- vec.ct
  means.ct <- sapply(split(seq_along(vec.cl), vec.cl), function(idx)
    colMeans(mat.ct.sc[idx, , drop = FALSE], na.rm = TRUE))
  ## 两级命名：谱系 = 非 Cycling 类型的 argmax；若 Cycling 同时为全体最高分，
  ## 则记为 <谱系>_Cycling（增殖亚群保留谱系归属）
  vec.ln <- setdiff(vec.ct, "Cycling")
  vec.ct.cl <- vapply(seq_len(ncol(means.ct)), function(k) {
      v <- means.ct[, k]
      ok.ct <- mat.mk.cl[, k] >= MIN_MARKER_CLUSTER_MEAN
      ok.ln <- vec.ln[ok.ct[vec.ln]]
      if (length(ok.ln) == 0L) return("Unknown")
      w.ln <- v[ok.ln]
      if (all(is.nan(w.ln))) return("Unknown")
      ln <- names(w.ln)[which.max(w.ln)]
      if (ok.ct["Cycling"] && !is.na(v["Cycling"]) &&
          v["Cycling"] > max(w.ln, na.rm = TRUE)) {
        return(paste0(ln, "_Cycling"))
      }
      ln
    }, character(1L))
  names(vec.ct.cl) <- colnames(means.ct)
  obj$celltype <- factor(unname(vec.ct.cl[vec.cl]))
  tbl.cl <- table(obj$seurat_clusters)
  dat.ann <- rbind(dat.ann, data.frame(acc = acc,
    cluster = as.integer(names(vec.ct.cl)), celltype = unname(vec.ct.cl),
    n_cells = as.integer(tbl.cl[names(vec.ct.cl)]),
    stringsAsFactors = FALSE))
  p.ct <- Seurat::DimPlot(obj, reduction = "umap", group.by = "celltype",
    label = TRUE, label.size = 4) + sci_theme +
    ggplot2::labs(title = paste(acc, "CellType"))
  sci_save(p.ct, paste0("ScRNA_UMAP_", acc, "_CellType"), width = 9,
    height = 7)


  n.seg.hit <- vapply(vec.segment, function(sg)
    length(intersect(list.seg.gene[[sg]], rownames(obj))), integer(1L))
  if (any(n.seg.hit < MIN_REGION_GENES)) {
    stop(paste0("S12 STOP: ", acc, " segment gene match < ",
        MIN_REGION_GENES, ": ",
        paste(names(n.seg.hit)[n.seg.hit < MIN_REGION_GENES],
          collapse = ", ")))
  }
  cat(sprintf("[S12] %s segment gene hits: %s\n", acc,
    paste(paste0(sub("_.*", "", names(n.seg.hit)), "=", n.seg.hit),
      collapse = " ")))
  for (sg in vec.segment) {
    set.seed(SEED)
    obj <- Seurat::AddModuleScore(obj,
      features = list(intersect(list.seg.gene[[sg]], rownames(obj))),
      name = paste0("segscore_", sg, "_"))
  }
  ## 五段活性 × 细胞类型（每细胞段分按注释类型聚合；喂 S13 细胞层）
  dat.seg.ct.acc <- do.call(rbind, lapply(vec.segment, function(sg) {
      v <- obj@meta.data[[paste0("segscore_", sg, "_1")]]
      lst.byct <- split(v, as.character(obj$celltype))
      do.call(rbind, lapply(seq_along(lst.byct), function(i)
        data.frame(acc = acc, celltype = names(lst.byct)[i], segment = sg,
          mean_score = mean(lst.byct[[i]]),
          median_score = stats::median(lst.byct[[i]]),
          sd_score = stats::sd(lst.byct[[i]]),
          n_cells = length(lst.byct[[i]]), stringsAsFactors = FALSE)))
    }))
  dat.seg.ct <- rbind(dat.seg.ct, dat.seg.ct.acc)
  dat.ct.stat <- rbind(dat.ct.stat,
    do.call(rbind, lapply(vec.segment, function(sg) {
        v <- obj@meta.data[[paste0("segscore_", sg, "_1")]]
        ct <- as.character(obj$celltype)
        kw <- stats::kruskal.test(v ~ ct)
        mm <- tapply(v, ct, mean)
        data.frame(acc = acc, segment = sg, kw_p = kw$p.value,
          top_celltype = names(mm)[which.max(mm)], max_mean = max(mm),
          min_mean = min(mm), n_celltypes = length(mm),
          stringsAsFactors = FALSE)
      })))
  dat.dot <- dat.seg.ct.acc
  dat.dot$segment <- factor(dat.dot$segment, levels = vec.segment)
  p.dot <- ggplot2::ggplot(dat.dot,
      ggplot2::aes(x = segment, y = celltype, size = mean_score,
        color = mean_score)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "white",
      high = "#B2182B", midpoint = 0) +
    ggplot2::scale_size_continuous(range = c(1, 7)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
      hjust = 1)) + sci_theme +
    ggplot2::labs(x = "Inflammation segment", y = "Cell type",
      size = "Mean module score", color = "Mean module score",
      title = paste0(acc, ": segment activity by cell type"))
  sci_save(p.dot, paste0("Segment_By_CellType_Dot_", acc), width = 9,
    height = 7)

  if (acc == ACC_RAPA) {
    dat.seg.sample <- rbind(dat.seg.sample,
      do.call(rbind, lapply(vec.segment, function(sg) {
          v <- obj@meta.data[[paste0("segscore_", sg, "_1")]]
          lst.s <- split(v, as.character(obj$orig.ident))
          do.call(rbind, lapply(seq_along(lst.s), function(i)
            data.frame(acc = acc, sample = names(lst.s)[i], segment = sg,
              mean_score = mean(lst.s[[i]]), n_cells = length(lst.s[[i]]),
              stringsAsFactors = FALSE)))
        })))
    is.cyc <- grepl("_Cycling$", as.character(obj$celltype))
    lst.cy <- split(is.cyc, as.character(obj$orig.ident))
    dat.cycl <- rbind(dat.cycl,
      do.call(rbind, lapply(seq_along(lst.cy), function(i)
        data.frame(acc = acc, sample = names(lst.cy)[i],
          n_cells = length(lst.cy[[i]]), n_cycling = sum(lst.cy[[i]]),
          frac_cycling = mean(lst.cy[[i]]), stringsAsFactors = FALSE))))
  }
  ## 簇标记基因（聚类 QC 与注释审计依据；只取正标记，min.pct/logfc 过滤）。
  f.mk <- file.path(CKPT_DIR, paste0("S12_", acc, "_markers.rds"))
  if (file.exists(f.mk)) {
    dat.markers <- readRDS(f.mk)
  } else {
    dat.markers <- Seurat::FindAllMarkers(obj, assay = "RNA", only.pos = TRUE,
      min.pct = 0.25, logfc.threshold = 0.25)
    saveRDS(dat.markers, f.mk)
  }
  dat.markers$acc <- acc
  lst.markers[[acc]] <- dat.markers
  if (acc == ACC_RAPA) {
    mat.cnt <- SeuratObject::LayerData(obj[["RNA"]], layer = "counts")
    lst.p <- strsplit(as.character(obj$orig.ident), "_", fixed = TRUE)
    vec.flow <- vapply(lst.p, `[`, character(1L), 3L)
    vec.cond <- vapply(lst.p, `[`, character(1L), 4L)
    vec.sid <- vapply(lst.p, function(x) paste(x[seq_len(4L)],
      collapse = "_"), character(1L))
    vec.sid.f <- factor(vec.sid)
    mat.pb <- sapply(levels(vec.sid.f), function(s)
      rowSums(mat.cnt[, vec.sid.f == s, drop = FALSE]))
    ## 样本属性按名取回（禁位置对齐：rowsum/sapply 结果序 ≠ 细胞序）
    dat.attr <- data.frame(sid = vec.sid, flow = vec.flow, cond = vec.cond,
      stringsAsFactors = FALSE)
    dat.attr <- dat.attr[!duplicated(dat.attr$sid), , drop = FALSE]
    dat.attr <- dat.attr[match(colnames(mat.pb), dat.attr$sid), , drop = FALSE]
    stopifnot(!any(is.na(dat.attr$sid)))
    dat.attr$cond <- factor(dat.attr$cond,
      levels = c("Control", "Rapamycin", "Paclitaxel"))
    dat.attr$flow <- factor(dat.attr$flow,
      levels = c("PhysioFlow", "AlteredFlow"))
    stopifnot(!any(is.na(dat.attr$cond)), !any(is.na(dat.attr$flow)))
    lib.sz <- colSums(mat.pb)
    mat.lcpm <- log2(t(t(mat.pb) / lib.sz) * 1e6 + 1)
    des <- stats::model.matrix(~ flow + cond, data = dat.attr)
    fit <- limma::eBayes(limma::lmFit(mat.lcpm, des))
    cm.gene <- limma::makeContrasts(
      Rapamycin_vs_Control = condRapamycin,
      Paclitaxel_vs_Control = condPaclitaxel,
      Rapamycin_vs_Paclitaxel = condRapamycin - condPaclitaxel,
      levels = des)
    fit.g <- limma::eBayes(limma::contrasts.fit(fit, cm.gene))
    for (cmp in colnames(cm.gene)) {
      tt <- limma::topTable(fit.g, coef = cmp, number = Inf,
        adjust.method = "BH")
      for (sg in vec.segment) {
        g.sg <- intersect(list.seg.gene[[sg]], rownames(tt))
        tt.sg <- tt[g.sg, , drop = FALSE]
        sig <- tt.sg$adj.P.Val < PVAL_CUTOFF
        n.up <- sum(sig & tt.sg$logFC >= LOGFC_MIN, na.rm = TRUE)
        n.down <- sum(sig & tt.sg$logFC <= -LOGFC_MIN, na.rm = TRUE)
        n.sig <- n.up + n.down
        dir <- if (n.sig == 0L) "indeterminate"
          else if (n.up > n.down) "up"
          else if (n.down > n.up) "down"
          else "mixed"
        dat.dir <- rbind(dat.dir, data.frame(acc = acc,
          arm = cmp, segment = sg,
          n_genes_tested = length(g.sg), n.sig = n.sig, n.up = n.up,
          n.down = n.down,
          dir.index = if (n.sig > 0L) (n.up - n.down) / n.sig else NA_real_,
          direction = dir,
          mean.logFC.seg = mean(tt.sg$logFC, na.rm = TRUE),
          median.logFC.seg = stats::median(tt.sg$logFC, na.rm = TRUE),
          ## 臂级最小 BH（全空时 disclose 检验力上限，随 CSV 走）
          min.BH = signif(min(tt$adj.P.Val, na.rm = TRUE), 3),
          stringsAsFactors = FALSE))
      }
    }
  }

 
  dat.sub <- data.frame(acc = acc, cluster = as.integer(names(tbl.cl)),
    n_cells = as.integer(tbl.cl), stringsAsFactors = FALSE)
  dat.cluster <- rbind(dat.cluster, dat.sub)
}
## 合并写簇标记表（聚类 QC 依据；命名注释需另指定参考，不在本步）
dat.markers.all <- do.call(rbind, lst.markers)
write.csv(dat.markers.all, file.path(TAB_DIR, "scRNA_markers.csv"),
  row.names = FALSE)
## MIN_CELLS_PER_CLUSTER 生效：小簇标记后剔除，并在收口处披露被剔除的簇
dat.cluster$pass_min_cells <- dat.cluster$n_cells >= MIN_CELLS_PER_CLUSTER
dat.cluster <- dat.cluster[dat.cluster$pass_min_cells, , drop = FALSE]
write.csv(dat.cluster, file.path(TAB_DIR, "scRNA_clusters.csv"),
  row.names = FALSE)
write.csv(dat.ann, file.path(TAB_DIR, "scRNA_celltype_annotation.csv"),
  row.names = FALSE)
write.csv(dat.seg.ct, file.path(TAB_DIR, "segment_by_celltype.csv"),
  row.names = FALSE)
dat.dir$p.down <- ifelse(dat.dir$n.sig > 0L, dat.dir$n.down / dat.dir$n.sig,
  NA_real_)
ok.ci <- !is.na(dat.dir$p.down)
z.ci <- stats::qnorm(0.975)
center <- (dat.dir$p.down[ok.ci] + z.ci^2 / (2 * dat.dir$n.sig[ok.ci])) /
  (1 + z.ci^2 / dat.dir$n.sig[ok.ci])
half <- z.ci * sqrt(dat.dir$p.down[ok.ci] * (1 - dat.dir$p.down[ok.ci]) /
    dat.dir$n.sig[ok.ci] + z.ci^2 / (4 * dat.dir$n.sig[ok.ci]^2)) /
  (1 + z.ci^2 / dat.dir$n.sig[ok.ci])
dat.dir$ci.lo <- NA_real_
dat.dir$ci.hi <- NA_real_
dat.dir$ci.lo[ok.ci] <- 1 - 2 * pmin(1, center + half)
dat.dir$ci.hi[ok.ci] <- 1 - 2 * pmax(0, center - half)
write.csv(dat.dir, file.path(TAB_DIR, "gse212388_rapa_direction.csv"),
  row.names = FALSE)
dat.dir$lab <- paste0(dat.dir$arm, " (n.sig=", dat.dir$n.sig,
  ", min BH=", dat.dir$min.BH, ")")
p.forest <- ggplot2::ggplot(dat.dir,
    ggplot2::aes(x = mean.logFC.seg, y = segment, color = arm)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
    color = "grey55") +
  ggplot2::geom_point(size = 2.8,
    position = ggplot2::position_dodge(width = 0.45)) +
  ggplot2::scale_color_manual(values = c(Rapamycin_vs_Control = "#B2182B",
    Paclitaxel_vs_Control = "#2166AC", Rapamycin_vs_Paclitaxel = "#4C4C4C"),
    name = "Contrast") +
  sci_theme +
  ggplot2::labs(x = "Mean log2FC of segment genes (descriptive)", y = NULL,
    title = "GSE212388 EC: segment-level shift under anti-proliferative drugs",
    subtitle = paste0("Pseudobulk per sample; limma ~ flow + condition; ",
      "no BH-significant genes (all arms indeterminate)"),
    caption = paste0("Significance readout: n.sig / dir.index / Wilson CI ",
      "in gse212388_rapa_direction.csv; do not cite point estimates alone."))
sci_save(p.forest, "EC_Segment_Direction_Forest", width = 9.5, height = 5)

dat.ss <- dat.seg.sample
lst.sp <- strsplit(dat.ss$sample, "_", fixed = TRUE)
dat.ss$flow <- factor(vapply(lst.sp, `[`, character(1L), 3L),
  levels = c("PhysioFlow", "AlteredFlow"))
dat.ss$cond <- factor(vapply(lst.sp, `[`, character(1L), 4L),
  levels = c("Control", "Rapamycin", "Paclitaxel"))
stopifnot(!any(is.na(dat.ss$flow)), !any(is.na(dat.ss$cond)))
vec.samp <- unique(dat.ss$sample)
## 长转宽：样本 x 段（按名对齐，禁位置对齐）
mat.sc <- sapply(vec.segment, function(sg) {
    d <- dat.ss[dat.ss$segment == sg, , drop = FALSE]
    d$mean_score[match(vec.samp, d$sample)]
  })
rownames(mat.sc) <- vec.samp
dat.design <- dat.ss[match(vec.samp, dat.ss$sample), c("sample", "flow",
  "cond")]
stopifnot(identical(rownames(mat.sc), dat.design$sample))
des.m <- stats::model.matrix(~ flow + cond, data = dat.design)
cm.m <- limma::makeContrasts(
  Rapamycin_vs_Control = condRapamycin,
  Paclitaxel_vs_Control = condPaclitaxel,
  Rapamycin_vs_Paclitaxel = condRapamycin - condPaclitaxel, levels = des.m)
fit.m <- limma::eBayes(limma::contrasts.fit(limma::lmFit(t(mat.sc), des.m),
  cm.m))
dat.seg.stat <- do.call(rbind, lapply(colnames(cm.m), function(cmp) {
    tt <- limma::topTable(fit.m, coef = cmp, number = Inf,
      adjust.method = "BH", confint = TRUE)
    data.frame(arm = cmp, segment = rownames(tt), logFC = tt$logFC,
      CI.L = tt$CI.L, CI.R = tt$CI.R, P.Value = tt$P.Value,
      adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
  }))
map.ref <- c(Rapamycin_vs_Control = "Control",
  Paclitaxel_vs_Control = "Control", Rapamycin_vs_Paclitaxel = "Paclitaxel")
map.tst <- c(Rapamycin_vs_Control = "Rapamycin",
  Paclitaxel_vs_Control = "Paclitaxel",
  Rapamycin_vs_Paclitaxel = "Rapamycin")
dat.seg.stat$mean_ref <- vapply(seq_len(nrow(dat.seg.stat)), function(i)
  mean(dat.ss$mean_score[dat.ss$segment == dat.seg.stat$segment[i] &
    dat.ss$cond == map.ref[[dat.seg.stat$arm[i]]]]), numeric(1L))
dat.seg.stat$mean_test <- vapply(seq_len(nrow(dat.seg.stat)), function(i)
  mean(dat.ss$mean_score[dat.ss$segment == dat.seg.stat$segment[i] &
    dat.ss$cond == map.tst[[dat.seg.stat$arm[i]]]]), numeric(1L))
write.csv(dat.seg.stat, file.path(TAB_DIR,
  "gse212388_segment_score_stats.csv"), row.names = FALSE)
write.csv(dat.seg.sample, file.path(TAB_DIR,
  "gse212388_segment_score_by_sample.csv"), row.names = FALSE)


des.i <- stats::model.matrix(~ flow * cond, data = dat.design)

vec.raw <- colnames(des.i)
vec.ic <- grep("^flow.*:cond", vec.raw, value = TRUE)
colnames(des.i) <- gsub(":", "_", vec.raw, fixed = TRUE)
stopifnot(length(vec.ic) == nlevels(dat.design$cond) - 1L)
fit.raw <- limma::lmFit(t(mat.sc), des.i)
cm.inter <- limma::makeContrasts(
  Rapamycin = flowAlteredFlow_condRapamycin,
  Paclitaxel = flowAlteredFlow_condPaclitaxel, levels = des.i)
fit.i <- limma::eBayes(limma::contrasts.fit(fit.raw, cm.inter))
dat.inter <- do.call(rbind, lapply(seq_along(vec.ic), function(k) {
    tt <- limma::topTable(fit.i, coef = colnames(cm.inter)[k],
      number = Inf, adjust.method = "BH", confint = TRUE)
    data.frame(term = vec.ic[k], segment = rownames(tt), logFC = tt$logFC,
      CI.L = tt$CI.L, CI.R = tt$CI.R, P.Value = tt$P.Value,
      adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
  }))
dat.inter$drug <- sub("^.*:cond", "", dat.inter$term)
stopifnot(nrow(dat.inter) == length(vec.ic) * length(vec.segment))
write.csv(dat.inter, file.path(TAB_DIR,
  "gse212388_segment_flow_interaction.csv"), row.names = FALSE)


cm.eff <- limma::makeContrasts(
  Rapa_Physio = condRapamycin,
  Rapa_Altered = condRapamycin + flowAlteredFlow_condRapamycin,
  Pac_Physio = condPaclitaxel,
  Pac_Altered = condPaclitaxel + flowAlteredFlow_condPaclitaxel,
  Rapa_interact = flowAlteredFlow_condRapamycin,
  Pac_interact = flowAlteredFlow_condPaclitaxel, levels = des.i)
fit.e <- limma::eBayes(limma::contrasts.fit(fit.raw, cm.eff))
dat.eff <- do.call(rbind, lapply(colnames(cm.eff), function(cc) {
    tt <- limma::topTable(fit.e, coef = cc, number = Inf,
      adjust.method = "BH", confint = TRUE)
    data.frame(contrast = cc, segment = rownames(tt), logFC = tt$logFC,
      CI.L = tt$CI.L, CI.R = tt$CI.R, P.Value = tt$P.Value,
      adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
  }))
stopifnot(nrow(dat.eff) == ncol(cm.eff) * length(vec.segment))
write.csv(dat.eff, file.path(TAB_DIR,
  "gse212388_segment_effect_by_flow.csv"), row.names = FALSE)

## 分组均值（flow x cond 共 6 组）：分层描述的直接读数
dat.grp <- do.call(rbind, lapply(vec.segment, function(sg) {
    do.call(rbind, lapply(levels(dat.design$flow), function(fl) {
        do.call(rbind, lapply(levels(dat.design$cond), function(cd) {
            v <- dat.ss$mean_score[dat.ss$segment == sg &
              dat.ss$flow == fl & dat.ss$cond == cd]
            data.frame(segment = sg, flow = fl, cond = cd,
              n_samples = length(v),
              mean_score = if (length(v) > 0L) mean(v) else NA_real_,
              stringsAsFactors = FALSE)
          }))
      }))
  }))
write.csv(dat.grp, file.path(TAB_DIR,
  "gse212388_segment_score_by_group.csv"), row.names = FALSE)
p.seg <- ggplot2::ggplot(dat.seg.stat,
    ggplot2::aes(x = logFC, y = segment, color = arm)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
    color = "grey55") +
  ggplot2::geom_linerange(ggplot2::aes(xmin = CI.L, xmax = CI.R),
    linewidth = 0.7,
    position = ggplot2::position_dodge(width = 0.45)) +
  ggplot2::geom_point(size = 2.8,
    position = ggplot2::position_dodge(width = 0.45)) +
  ggplot2::scale_color_manual(values = c(Rapamycin_vs_Control = "#B2182B",
    Paclitaxel_vs_Control = "#2166AC", Rapamycin_vs_Paclitaxel = "#4C4C4C"),
    name = "Contrast") +
  sci_theme +
  ggplot2::labs(x = "Module score change (limma estimate)", y = NULL,
    title = "GSE212388 EC: sample-level segment score shift",
    subtitle = "limma ~ flow + condition; per-sample segment module scores",
    caption = paste0("Exploratory: n = 2 drug vs 4 control samples; ",
      "drug-vs-drug contrast is the specificity check (Rapa vs Paclitaxel)"))
sci_save(p.seg, "EC_SegmentScore_Forest", width = 9.5, height = 5)


seg.short <- c(S1_initiation = "S1", S2_amplification = "S2",
  S3_effect = "S3", S4_resolution = "S4", S5_chronic = "S5")
dat.a1 <- dat.eff[grepl("_(Physio|Altered)$", dat.eff$contrast), ]
dat.a1$drug <- ifelse(grepl("^Rapa_", dat.a1$contrast), "Rapamycin",
  "Paclitaxel")
dat.a1$flow <- factor(ifelse(grepl("_Physio$", dat.a1$contrast),
  "PhysioFlow", "AlteredFlow"), levels = c("PhysioFlow", "AlteredFlow"))
dat.a1$seg <- unname(seg.short[dat.a1$segment])
stopifnot(!any(is.na(dat.a1$seg)),
  nrow(dat.a1) == 4L * length(vec.segment))
dat.a2 <- dat.eff[grepl("_interact$", dat.eff$contrast), ]
dat.a2$drug <- ifelse(grepl("^Rapa_", dat.a2$contrast), "Rapamycin",
  "Paclitaxel")
dat.a2$seg <- unname(seg.short[dat.a2$segment])
stopifnot(!any(is.na(dat.a2$seg)),
  nrow(dat.a2) == 2L * length(vec.segment))
lst.pl <- lapply(list(dat.a1, dat.a2), function(d) {
    ggplot2::ggplot(d, ggplot2::aes(x = seg, y = logFC, color = drug,
        group = drug)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
        color = "grey40", linewidth = 0.5) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = CI.L, ymax = CI.R),
        width = 0.12, linewidth = 0.6,
        position = ggplot2::position_dodge(width = 0.5)) +
      ggplot2::geom_point(size = 2.6,
        position = ggplot2::position_dodge(width = 0.5)) +
      ggplot2::scale_color_manual(values = c(Rapamycin = "#B2182B",
        Paclitaxel = "#2166AC"), name = "Treatment")
  })
p.eff <- lst.pl[[1]] + ggplot2::facet_grid(~ flow) +
  ggplot2::labs(x = NULL, y = "Module score change vs Control")
p.int <- lst.pl[[2]] +
  ggplot2::labs(x = "Inflammation segment",
    y = "Effect difference (module score)\nAltered - Physio flow")
p.flow <- (p.eff / p.int + patchwork::plot_layout(guides = "collect")) &
  sci_theme &
  ggplot2::theme(legend.position = "right")
p.flow <- p.flow + patchwork::plot_annotation(
  tag_levels = "A",
  title = paste0("Drug effects on inflammatory segment scores under ",
    "physiological vs disturbed flow"),
  subtitle = paste0("GSE212388 EC; linear model on sample-level module ",
    "scores; n = 2 vs 4 (Control)"),
  caption = paste0(
    "A, drug effect vs matched Control within each flow; ",
    "B, interaction (effect change between flow backgrounds).\n",
    "Segments: S1 initiation / S2 amplification / S3 effect / ",
    "S4 resolution / S5 chronic.  AlteredFlow = pro-inflammatory flow.\n",
    "Error bars = 95% CI.  No contrast reached FDR < 0.05; ",
    "interaction model residual df = 2 (exploratory)."))
sci_save(p.flow, "EC_Segment_Effect_By_Flow", width = 9.5, height = 7)


dat.cy <- dat.cycl
lst.cp <- strsplit(dat.cy$sample, "_", fixed = TRUE)
dat.cy$flow <- factor(vapply(lst.cp, `[`, character(1L), 3L),
  levels = c("PhysioFlow", "AlteredFlow"))
dat.cy$cond <- factor(vapply(lst.cp, `[`, character(1L), 4L),
  levels = c("Control", "Rapamycin", "Paclitaxel"))
kw.cy <- stats::kruskal.test(frac_cycling ~ cond, data = dat.cy)
dat.cy.sum <- do.call(rbind, lapply(levels(dat.cy$cond), function(cd)
  data.frame(cond = cd, n_samples = sum(dat.cy$cond == cd),
    mean_frac = mean(dat.cy$frac_cycling[dat.cy$cond == cd]),
    stringsAsFactors = FALSE)))
dat.cy.sum$kw_p <- signif(kw.cy$p.value, 3)
write.csv(dat.cy.sum, file.path(TAB_DIR,
  "gse212388_cycling_fraction.csv"), row.names = FALSE)
p.cy <- ggplot2::ggplot(dat.cy,
    ggplot2::aes(x = cond, y = frac_cycling, color = cond)) +
  ggplot2::geom_point(size = 3,
    position = ggplot2::position_jitter(width = 0.12)) +
  ggplot2::scale_color_manual(values = c(Control = "#4C72B0",
    Rapamycin = "#B2182B", Paclitaxel = "#2166AC")) +
  ggplot2::guides(color = "none") + sci_theme +
  ggplot2::labs(x = "Condition", y = "Fraction of cycling cells",
    title = "GSE212388 EC: cycling fraction by treatment",
    subtitle = "Cycling = clusters annotated *_Cycling",
    caption = paste0("Kruskal-Wallis p = ", signif(kw.cy$p.value, 3),
      "; n = ", nrow(dat.cy), " samples (2 drugs x 2 flow x 1-2 replicates)"))
sci_save(p.cy, "EC_Cycling_Fraction_By_Treatment", width = 7, height = 5)


dat.ct.stat$bh <- ave(dat.ct.stat$kw_p, dat.ct.stat$acc,
  FUN = function(x) stats::p.adjust(x, method = "BH"))
write.csv(dat.ct.stat, file.path(TAB_DIR, "segment_by_celltype_stats.csv"),
  row.names = FALSE)

## PART 2 空间：GeoMx dcc.gz ---------------------------------------------------------------------

## q11 的判定单位是「Panel 命中的段基因个数」（MIN_REGION_GENES），
## 不是「有命中的段数」——旧写法把两者混用，量纲不同，此处已更正。
dat.spatial.qc <- data.frame(acc = character(0), panel = character(0),
  n_gene_panel = integer(0),
  n_segment_genes_hit = integer(0),
  n_segment_covered = integer(0),
  q11 = character(0), stringsAsFactors = FALSE)


f.pkc <- "01_data/GSE277170/GSE326968_Hs_R_NGS_WTA_v1.0.pkc.gz"
if (!file.exists(f.pkc)) {
  stop(paste0("S12 STOP: GeoMx WTA pkc 缺失: ", f.pkc,
    " ; 请从 GEO GSE326968 suppl 下载后放入该路径。"))
}
con.pkc <- gzfile(f.pkc, "rt")
txt.pkc <- readLines(con.pkc)
close(con.pkc)
pkg.pkc <- jsonlite::fromJSON(paste(txt.pkc, collapse = "\n"))
if (is.null(pkg.pkc$Targets)) {
  stop(paste0("S12 STOP: pkc 解析失败，Targets 键缺失（pkc schema",
    "可能已变更，须复核）。"))
}
probes.list <- pkg.pkc$Targets$Probes
vec.rts.id <- character(0L)
vec.rts.sym <- character(0L)
for (i in seq_along(probes.list)) {
  rid <- probes.list[[i]]$RTS_ID
  sym <- pkg.pkc$Targets$DisplayName[i]
  if (length(rid) == 1L && length(sym) == 1L) {
    vec.rts.id <- c(vec.rts.id, rid)
    vec.rts.sym <- c(vec.rts.sym, sym)
  }
}
map.sym <- setNames(vec.rts.sym,
  vec.rts.id)   ## RTS_ID -> 干净基因符号（命名向量查表）
cat(sprintf("[S12] pkc map built: %d RTS -> %d unique symbols\n",
  length(map.sym), length(unique(vec.rts.sym))))

## 返回单个 dcc 的 Code_Summary 探针行数（RTS 号条数）。
CountRts <- function(f.dcc) {
  con <- gzfile(f.dcc, "rt")
  vec.ln <- readLines(con)
  close(con)
  i.beg <- which(startsWith(vec.ln, "<Code_Summary>"))
  i.end <- which(startsWith(vec.ln, "</Code_Summary>"))
  if (!length(i.beg) || !length(i.end)) return(0L)
  as.integer(i.end[1L] - i.beg[1L] - 1L)
}

for (acc in "GSE277170") {
  vec.dcc <- list.files(file.path("01_data", acc), pattern = "\\.dcc\\.gz$")
  if (length(vec.dcc) == 0L) {
    stop(paste0("S12 STOP: no dcc.gz under 01_data/", acc, "; run S12a first."))
  }
  
  vec.n.rts <- vapply(file.path("01_data", acc, vec.dcc), CountRts,
    integer(1L), USE.NAMES = FALSE)
  dat.roi.qc <- data.frame(roi = sub("\\.dcc\\.gz$", "", vec.dcc),
    rts_id_count = vec.n.rts,
    pass_min_rts = vec.n.rts >= MIN_ROI_RTS, stringsAsFactors = FALSE)
  write.csv(dat.roi.qc, file.path(TAB_DIR, "spatial_roi_qc.csv"),
    row.names = FALSE)
  cat(sprintf(
    "[S12] ROI QC: %d dcc, rts min=%d median=%d max=%d; %d below %d\n",
    length(vec.dcc), min(vec.n.rts),
    as.integer(stats::median(vec.n.rts)), max(vec.n.rts),
    sum(!dat.roi.qc$pass_min_rts), MIN_ROI_RTS))
  if (!any(dat.roi.qc$pass_min_rts)) {
    stop(paste0("S12 STOP: no dcc with >= ", MIN_ROI_RTS,
        " RTS probes under 01_data/", acc, "; all ROI failed sequencing."))
  }
  f.dcc <- file.path("01_data", acc, vec.dcc[which.max(vec.n.rts)])
  vec.rts <- character(0)
  in.code <- FALSE
  con <- gzfile(f.dcc, "rt")
  repeat {
    ln <- readLines(con, n = 1L, warn = FALSE)
    if (length(ln) == 0L) break
    if (startsWith(ln, "<Code_Summary>")) { in.code <- TRUE; next }
    if (startsWith(ln, "</Code_Summary>")) break
    if (in.code) vec.rts <- c(vec.rts, strsplit(ln, ",")[[1L]][1L])
  }
  close(con)
  ## pkc 映射（map.sym 已在上方构建：RTS_ID -> 干净基因符号，命名向量查表）
  vec.gene <- unique(stats::na.omit(map.sym[match(vec.rts, names(map.sym))]))
  n.seg.hit <- sum(vapply(list.seg.gene, function(g) any(g %in% vec.gene),
      TRUE))
  n.seg.gene.hit <- length(intersect(unique(unlist(list.seg.gene)), vec.gene))
  dat.spatial.qc <- rbind(dat.spatial.qc, data.frame(
      acc = acc, panel = "WTA", n_gene_panel = length(vec.gene),
      n_segment_genes_hit = n.seg.gene.hit,
      n_segment_covered = as.integer(n.seg.hit),
      q11 = if (n.seg.gene.hit >= MIN_REGION_GENES) "PASS" else "FAIL",
      stringsAsFactors = FALSE))
}
write.csv(dat.spatial.qc, file.path(TAB_DIR, "spatial_panel_qc.csv"),
  row.names = FALSE)

## PART 2.5 图件 -------------------------------------------------------------------------------
## 命名遵守规范 16（描述性英文，禁 Fig 前缀）；图内文字全英文（规范 17）；
## 分面一律用 facet_wrap，不引入任何拼图包。

dat.cluster$acc <- factor(dat.cluster$acc)
sci_save(ggplot(dat.cluster, aes(x = cluster, y = n_cells)) +
  geom_col(fill = "#4C72B0") +
  facet_wrap(~ acc, scales = "free_x") +
  labs(title = "scRNA clusters kept per series",
    x = "Cluster", y = "Cells") +
  sci_theme,
  "ScRNA_Cluster_Size_By_Series", width = 10, height = 4.5)



## PART 3 收口 -------------------------------------------------------------------------------

vec.q11.fail <- dat.spatial.qc$panel[dat.spatial.qc$q11 == "FAIL"]
message(paste0("\n===== S12 cell/space summary =====\nscRNA clusters kept = ",
    nrow(dat.cluster), " (>= ", MIN_CELLS_PER_CLUSTER, " cells)",
    "  spatial panels checked = ", nrow(dat.spatial.qc),
    "\nQ11 failures = ", if (length(vec.q11.fail) == 0L) "none"
    else paste(vec.q11.fail, collapse = ", ")))
if (length(vec.q11.fail) > 0L) {
  stop(paste0("S12 FAILED: spatial panel(s) failing Q11: ", paste(vec.q11.fail,
        collapse = ", "),
      " ; see 03_results/S12/table/spatial_panel_qc.csv ."))
}
message("S12 OK. tables -> 03_results/S12/table/")
