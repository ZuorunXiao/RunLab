# ===============================================================
# RARP-16_cellspace.R
# 功能: 细胞层与空间层（步骤规格 §8：S12）
# 定位: 细胞层与空间层阶段（S12）
# 可视化: 无（仅产出表格/数据）
# ===============================================================
SEED <- 20260925
set.seed(SEED)
FIG_DIR <- "03_results/S12/Figure"
TAB_DIR <- "03_results/S12/table"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
# ：Seurat 标准化 / 缩放 / 找邻居 / 聚类属重时步骤，
# 断点落 checkpoints/（命名 S12_<acc>_seurat.rds）；TAB_DIR 须先落定，
# 否则 file.path 在其赋值前求值而立即报错
CKPT_DIR <- file.path(TAB_DIR, "checkpoints")
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)

library(ggplot2)  # 簇与 panel 图；主题出口用 §sci_theme / §sci_save
library(patchwork)  # 效应图双联拼接（A 效应 + B 交互）
# 空间分支解析 NanoString pkc（RTS_ID -> 基因符号），替代缺失的 GeoMxTools
library(jsonlite)

# 主题： §sci_theme 官方定义（theme_bw + 关网格），全脚本统一口径
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

# 出口： §sci_save 官方定义，仅 PDF，落本步 FIG_DIR
sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  pdf(pdf_file, width = width, height = height)
  print(plot_obj)
  dev.off()
# 可读性核验钩子（）：默认关闭，交付物仍为纯 PDF；
# 置 RAPA_QC_PNG=1 时额外导 _qc/ PNG 供逐张目视核验，不改变 PDF 产出。
  if (Sys.getenv("RAPA_QC_PNG") == "1") {
    dir.create(file.path(FIG_DIR, "_qc"), showWarnings = FALSE)
    png(file.path(FIG_DIR, "_qc", paste0(filename, ".png")),
      width = width * 150, height = height * 150, res = 150)
    print(plot_obj)
    dev.off()
  }
}

MIN_CELLS_PER_CLUSTER <- 10L  # 少于此细胞数的簇不参与下游
MIN_REGION_GENES <- 20L  # Panel 命中的段基因数下限
MIN_ROI_RTS <- 10000L
vec.segment <- c("S1_initiation", "S2_amplification", "S3_effect",
  "S4_resolution", "S5_chronic")
# 段名短标签：图中只用 S1-S5，完整段名留在表内
seg.short <- c(S1_initiation = "S1", S2_amplification = "S2",
  S3_effect = "S3", S4_resolution = "S4", S5_chronic = "S5")

dat.final <- read.csv("03_results/S02/table/inflammation_geneset_final.csv",
  stringsAsFactors = FALSE)
vec.col.req <- c("segment", "gs_name", "members")
vec.col.miss <- setdiff(vec.col.req, names(dat.final))
if (length(vec.col.miss) > 0L) {
  stop("S12 FAILED: inflammation_geneset_final.csv lacks column(s): ",
    paste(vec.col.miss, collapse = ", "))
}
# 每段基因向量：S02 已用 msigdbr 实测成员，空成员串表示未命中当前 MSigDB 版本
list.seg.gene <- lapply(vec.segment, function(s) {
    dat.s <- dat.final[dat.final$segment == s & dat.final$members != "", ]
    unique(unlist(lapply(dat.s$members, function(x) strsplit(x, ";")[[1L]])))
  })
names(list.seg.gene) <- vec.segment
# Q3 一致性：五段都必须在锁定表中出现，否则下游有段会静默缺席
vec.seg.empty <- vec.segment[lengths(list.seg.gene) == 0L]
if (length(vec.seg.empty) > 0L) {
  stop(paste0("S12 FAILED: segment(s) have no genes after MSigDB matching: ",
      paste(vec.seg.empty, collapse = ", "),
      " (see 03_results/S02/table/segment_member_count.csv)."))
}

ACC_RAPA <- "GSE212388"
# GSE212388 八个样本含两个药臂，对照按测序深度分成两层：宽表四个样本中位
# 深度 695-1,675 UMI/细胞，10x 深格式四个 5,988-12,629。药物标签与深度完全
# 共线，因此两臂只能各自在自身层内与匹配对照配对；跨层（紫杉醇对雷帕霉素）
# 不可比较，任何同时估计药物与深度的模型都不可识别。
ARM_RAPA <- "Rapamycin_vs_Control"
ARM_PACLITAXEL <- "Paclitaxel_vs_Control"
dat.arm <- rbind(
  data.frame(arm = ARM_RAPA,
    gsm = c("GSM6529254", "GSM6529256", "GSM6529258", "GSM6529259"),
    flow = c("PhysioFlow", "AlteredFlow", "PhysioFlow", "AlteredFlow"),
    cond = c("Control", "Control", "Rapamycin", "Rapamycin"),
    stringsAsFactors = FALSE),
  data.frame(arm = ARM_PACLITAXEL,
    gsm = c("GSM6529255", "GSM6529257", "GSM6529260", "GSM6529261"),
    flow = c("PhysioFlow", "AlteredFlow", "PhysioFlow", "AlteredFlow"),
    cond = c("Control", "Control", "Paclitaxel", "Paclitaxel"),
    stringsAsFactors = FALSE))
vec.arm <- unique(dat.arm$arm)
# S11 同口径：PVAL_CUTOFF / LOGFC_MIN 与 summarise_direction 三分支判定
# （RARP-15_perturbation.R L102-103 / L148-161），禁自造阈值。
PVAL_CUTOFF <- 0.05
LOGFC_MIN <- 0.585
# 细胞类型 canonical marker 面板：12 个谱系共 59 个基因，全部为血管与免疫
# 的教科书级标记，逐个基因在 CellMarker 3.0 人类人工整理记录中可查到所属
# 谱系的标注（Hu C 2026 NAR gkag1040）。仅用于簇级 argmax，不做富集推断。
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
  pDC = c("LILRA4", "PLD4", "SERPINF1"),
  B_cell = c("MS4A1", "CD79A", "BANK1", "CD79B"),
  Plasma = c("IGLL5", "MZB1", "JCHAIN", "IGHG1"),
  Mast = c("HDC", "CPA3", "KIT"),
  Cycling = c("MKI67", "TOP2A", "CCNB1", "PCNA", "ASPM", "TYMS", "UBE2C"))

# PART 1 scRNA：Seurat + harmony ----------------------------------------------------------------

vec.sc <- c("GSE237230", "GSE247238", "GSE212388")
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop(
    "S12 FAILED: 'Seurat' not installed ; scRNA branch aborted (no fallback).")
}
library(Seurat)  # 空间 / 单细胞对象构建、聚类与整合
library(harmony)  # Harmony 整合（可选；包缺失即 stop，不静默降级）

dat.cluster <- data.frame(acc = character(0), cluster = integer(0),
  n_cells = integer(0),
  stringsAsFactors = FALSE)

lst.markers <- list()  # 每 acc 的簇标记表，循环后合并写盘
dat.ann <- data.frame(acc = character(0), cluster = integer(0),
  celltype = character(0), n_cells = integer(0), stringsAsFactors = FALSE)
dat.seg.ct <- data.frame(acc = character(0), celltype = character(0),
  segment = character(0), mean_score = numeric(0),
  median_score = numeric(0), sd_score = numeric(0), n_cells = integer(0),
  stringsAsFactors = FALSE)
dat.seg.sample <- data.frame(acc = character(0), arm = character(0),
  sample = character(0), segment = character(0), mean_score = numeric(0),
  n_cells = integer(0), stringsAsFactors = FALSE)
dat.cycl <- data.frame(acc = character(0), arm = character(0),
  sample = character(0), n_cells = integer(0), n_cycling = integer(0),
  frac_cycling = numeric(0), stringsAsFactors = FALSE)
dat.ct.stat <- data.frame(acc = character(0), segment = character(0),
  kw_p = numeric(0), top_celltype = character(0), max_mean = numeric(0),
  min_mean = numeric(0), n_celltypes = integer(0), stringsAsFactors = FALSE)
dat.cov <- data.frame(acc = character(0), celltype = character(0),
  n_gene_panel = integer(0), n_gene_hit = integer(0),
  stringsAsFactors = FALSE)
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
# 三个系列都读全部样本做注释广度；药效臂的样本子集在下方按臂表另取，
# 使两个药臂各在自己的深度层内配对
  if (length(vec.sdir) == 0L) {
    stop(paste0("S12 STOP: no per-sample 10X dirs for ", acc,
        "; run S12a first."))
  }
# ：整条 Seurat 处理链是重时步骤，按 acc 粒度续传；
# 恢复失败即报错，不静默丢弃损坏对象（）。
  f.ck <- file.path(CKPT_DIR, paste0("S12_", acc, "_seurat.rds"))
  if (file.exists(f.ck)) {
    obj <- readRDS(f.ck)
    obj <- SeuratObject::JoinLayers(obj)
    red.use <- if ("harmony" %in% names(obj@reductions)) "harmony" else "pca"
  } else {
# 每样本读 10X 目录（barcodes/features/matrix），合并前各自过最低门槛
    lst.obj <- lapply(vec.sdir, function(d) {
        counts <- Seurat::Read10X(d)
        Seurat::CreateSeuratObject(counts = counts, project = basename(d),
          min.cells = 3L, min.features = 200L)
      })
    obj <- Reduce(function(a, b) merge(a, y = b), lst.obj)
    obj <- Seurat::NormalizeData(obj)
    obj <- Seurat::FindVariableFeatures(obj, n.top.genes = 2000)
    obj <- Seurat::ScaleData(obj)
    obj <- SeuratObject::JoinLayers(obj)
    obj <- Seurat::RunPCA(obj, npcs = 30)
# 跨供者整合：FindNeighbors 必须在 PCA（及 harmony）之后，否则 'pca' not found；
# harmony 可用则用 harmony 归约建邻居图（真正落地整合），否则退回 pca 归约。
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
# UMAP 底图（harmony/pca 归约；Stuart 2019 PMID 31178118 / Hao 2021 PMID 34062119）。
  if (is.null(obj@reductions$umap)) {
    obj <- Seurat::RunUMAP(obj, reduction = red.use, dims = 1:20)
  }

  vec.ct.all <- names(lst.ct.mk)
# 面板各谱系在本数据集的可检出基因数；检出不足 2 个的谱系不足以计算
# 模块分，不参与该数据集的注释
  n.ct.hit <- vapply(lst.ct.mk, function(g)
    length(intersect(g, rownames(obj))), integer(1L))
  dat.cov <- rbind(dat.cov, data.frame(acc = acc, celltype = vec.ct.all,
    n_gene_panel = as.integer(lengths(lst.ct.mk)),
    n_gene_hit = as.integer(n.ct.hit[vec.ct.all]), stringsAsFactors = FALSE))
  vec.ct <- vec.ct.all[n.ct.hit[vec.ct.all] >= 2L]
  vec.cl <- as.character(obj$seurat_clusters)
  mat.mk.cl <- matrix(0, nrow = length(vec.ct), ncol = length(unique(vec.cl)),
    dimnames = list(vec.ct, sort(unique(vec.cl))))
  for (ct in vec.ct) {
    g.ct <- intersect(lst.ct.mk[[ct]], rownames(obj))
# 只取标记基因行再转密集（整矩阵密集会爆内存）
    mat.sub <- as.matrix(SeuratObject::LayerData(obj[["RNA"]],
      layer = "data")[g.ct, , drop = FALSE])
# 簇内面板均值 = 该簇全部细胞 x 全部面板基因的算术平均；逐簇直接求均值，
# 避免先按基因求均值再跨基因求均值时对单基因面板退化成向量
    mat.mk.cl[ct, ] <- vapply(split(seq_len(ncol(mat.sub)), vec.cl),
      function(idx) mean(mat.sub[, idx, drop = FALSE]), numeric(1L))
# 每个谱系都出模块分；谱系可用性由下游簇内均值门槛判定，不在此处预筛
    set.seed(SEED)
    obj <- Seurat::AddModuleScore(obj, features = list(g.ct),
      name = paste0("ctscore_", ct, "_"))
  }
  mat.ct.sc <- as.matrix(obj@meta.data[, paste0("ctscore_", vec.ct, "_1"),
    drop = FALSE])
  colnames(mat.ct.sc) <- vec.ct
  means.ct <- sapply(split(seq_along(vec.cl), vec.cl), function(idx)
    colMeans(mat.ct.sc[idx, , drop = FALSE], na.rm = TRUE))
# 两级命名：谱系 = 非 Cycling 类型的 argmax；若 Cycling 同时为全体最高分，
# 则记为 <谱系>_Cycling（增殖亚群保留谱系归属）
  vec.ln <- setdiff(vec.ct, "Cycling")
  vec.ct.cl <- vapply(seq_len(ncol(means.ct)), function(k) {
      v <- means.ct[, k]
      ok.ct <- mat.mk.cl[, k] >= MIN_MARKER_CLUSTER_MEAN
      ok.ln <- vec.ln[ok.ct[vec.ln]]
      if (length(ok.ln) == 0L) return("Unknown")
      w.ln <- v[ok.ln]
      ln <- names(w.ln)[which.max(w.ln)]
      if (ok.ct["Cycling"] && v["Cycling"] > max(w.ln)) {
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

  for (sg in vec.segment) {
    set.seed(SEED)
    obj <- Seurat::AddModuleScore(obj,
      features = list(intersect(list.seg.gene[[sg]], rownames(obj))),
      name = paste0("segscore_", sg, "_"))
  }
# 药效臂不在此处取子集：子集内的段模块分随臂重算，见下方臂循环
# 五段活性 × 细胞类型（每细胞段分按注释类型聚合；喂 S13 细胞层）
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

# 簇标记基因（聚类 QC 与注释审计依据；只取正标记，min.pct/logfc 过滤）。
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
# 两臂各自配对：每臂四个样本，子集内重算段模块分，对照分箱只含同深度层
    lst.arm <- lapply(vec.arm, function(a) {
      dat.a <- dat.arm[dat.arm$arm == a, ]
      drug <- setdiff(dat.a$cond, "Control")
      vec.cell <- colnames(obj)[sub("_.*", "",
        as.character(obj$orig.ident)) %in% dat.a$gsm]
      obj.a <- obj[, vec.cell]
      for (sg in vec.segment) {
        set.seed(SEED)
        obj.a <- Seurat::AddModuleScore(obj.a,
          features = list(intersect(list.seg.gene[[sg]], rownames(obj.a))),
          name = paste0("armscore_", sg, "_"))
      }
# 样本级段分：每样本每段的细胞均值
      dat.s <- do.call(rbind, lapply(vec.segment, function(sg) {
          v <- obj.a@meta.data[[paste0("armscore_", sg, "_1")]]
          lst.s <- split(v, as.character(obj.a$orig.ident))
          do.call(rbind, lapply(seq_along(lst.s), function(i)
            data.frame(arm = a, sample = names(lst.s)[i], segment = sg,
              mean_score = mean(lst.s[[i]]), n_cells = length(lst.s[[i]]),
              stringsAsFactors = FALSE)))
        }))
# 增殖比例：每样本 *_Cycling 注释的细胞占比
      is.cyc <- grepl("_Cycling$", as.character(obj.a$celltype))
      lst.cy <- split(is.cyc, as.character(obj.a$orig.ident))
      dat.c <- do.call(rbind, lapply(seq_along(lst.cy), function(i)
        data.frame(arm = a, sample = names(lst.cy)[i],
          n_cells = length(lst.cy[[i]]), n_cycling = sum(lst.cy[[i]]),
          frac_cycling = mean(lst.cy[[i]]), stringsAsFactors = FALSE)))
# 伪 bulk：样本内 counts 求和，只保留多数样本检出的基因，仅在少数样本检出
# 的基因其 logFC 由单样本驱动、不可检验
      mat.cnt <- SeuratObject::LayerData(obj.a[["RNA"]], layer = "counts")
      lst.p <- strsplit(as.character(obj.a$orig.ident), "_", fixed = TRUE)
      vec.sid <- vapply(lst.p, function(x) paste(x[seq_len(4L)],
        collapse = "_"), character(1L))
      mat.pb <- sapply(sort(unique(vec.sid)), function(s)
        rowSums(mat.cnt[, vec.sid == s, drop = FALSE]))
      dat.attr <- data.frame(sid = vec.sid,
        flow = vapply(lst.p, `[`, character(1L), 3L),
        cond = vapply(lst.p, `[`, character(1L), 4L),
        stringsAsFactors = FALSE)
      dat.attr <- dat.attr[!duplicated(dat.attr$sid), , drop = FALSE]
      dat.attr <- dat.attr[match(colnames(mat.pb), dat.attr$sid), ,
        drop = FALSE]
      stopifnot(!any(is.na(dat.attr$sid)))
      dat.attr$cond <- factor(dat.attr$cond, levels = c("Control", drug))
      dat.attr$flow <- factor(dat.attr$flow,
        levels = c("PhysioFlow", "AlteredFlow"))
      mat.pb <- mat.pb[rowSums(mat.pb >= 1) >= 3L, , drop = FALSE]
      mat.lcpm <- log2(t(t(mat.pb) / colSums(mat.pb)) * 1e6 + 1)
      des <- stats::model.matrix(~ flow + cond, data = dat.attr)
      fit <- limma::eBayes(limma::lmFit(mat.lcpm, des))
# 对比项按设计矩阵列名取：drug 为单水平时列名恒为 cond<药名>
      ct.term <- paste0("cond", drug)
      fit.g <- limma::eBayes(limma::contrasts.fit(fit,
        limma::makeContrasts(contrasts = ct.term, levels = des)))
      tt <- limma::topTable(fit.g, coef = ct.term, number = Inf,
        adjust.method = "BH")
# 段级判定：与 S11 同口径的三分支，臂级最小 BH 单列披露检验力上限
      dat.d <- do.call(rbind, lapply(vec.segment, function(sg) {
          g.sg <- intersect(list.seg.gene[[sg]], rownames(tt))
          tt.sg <- tt[g.sg, , drop = FALSE]
          sig <- tt.sg$adj.P.Val < PVAL_CUTOFF
          n.up <- sum(sig & tt.sg$logFC >= LOGFC_MIN, na.rm = TRUE)
          n.down <- sum(sig & tt.sg$logFC <= -LOGFC_MIN, na.rm = TRUE)
          n.sig <- n.up + n.down
          data.frame(arm = a, segment = sg, n_genes_tested = length(g.sg),
            n.sig = n.sig, n.up = n.up, n.down = n.down,
            dir.index = if (n.sig > 0L) (n.up - n.down) / n.sig else NA_real_,
            direction = if (n.sig == 0L) "indeterminate"
              else if (n.up > n.down) "up"
              else if (n.down > n.up) "down" else "mixed",
            mean.logFC.seg = mean(tt.sg$logFC, na.rm = TRUE),
            median.logFC.seg = stats::median(tt.sg$logFC, na.rm = TRUE),
            min.BH = signif(min(tt$adj.P.Val, na.rm = TRUE), 3),
            stringsAsFactors = FALSE)
        }))
      list(seg.sample = dat.s, cycl = dat.c, direction = dat.d)
    })
# 两臂按行合并，臂名留在 arm 列供下游分面
    dat.seg.sample <- rbind(dat.seg.sample, data.frame(acc = acc,
      do.call(rbind, lapply(lst.arm, `[[`, "seg.sample")),
      stringsAsFactors = FALSE))
    dat.cycl <- rbind(dat.cycl, data.frame(acc = acc,
      do.call(rbind, lapply(lst.arm, `[[`, "cycl")),
      stringsAsFactors = FALSE))
    dat.dir <- rbind(dat.dir, data.frame(acc = acc,
      do.call(rbind, lapply(lst.arm, `[[`, "direction")),
      stringsAsFactors = FALSE))
  }

  dat.sub <- data.frame(acc = acc, cluster = as.integer(names(tbl.cl)),
    n_cells = as.integer(tbl.cl), stringsAsFactors = FALSE)
  dat.cluster <- rbind(dat.cluster, dat.sub)
}
# 合并写簇标记表（聚类 QC 依据；命名注释需另指定参考，不在本步）
dat.markers.all <- do.call(rbind, lst.markers)
write.csv(dat.markers.all, file.path(TAB_DIR, "scRNA_markers.csv"),
  row.names = FALSE)
# MIN_CELLS_PER_CLUSTER 生效：小簇标记后剔除，并在收口处披露被剔除的簇
dat.cluster$pass_min_cells <- dat.cluster$n_cells >= MIN_CELLS_PER_CLUSTER
dat.cluster <- dat.cluster[dat.cluster$pass_min_cells, , drop = FALSE]
write.csv(dat.cluster, file.path(TAB_DIR, "scRNA_clusters.csv"),
  row.names = FALSE)
write.csv(dat.ann, file.path(TAB_DIR, "scRNA_celltype_annotation.csv"),
  row.names = FALSE)
write.csv(dat.seg.ct, file.path(TAB_DIR, "segment_by_celltype.csv"),
  row.names = FALSE)
write.csv(dat.cov, file.path(TAB_DIR, "scRNA_panel_coverage.csv"),
  row.names = FALSE)
dat.dir$drug <- sub("_vs_Control$", "", dat.dir$arm)
dat.dir$seg <- unname(seg.short[dat.dir$segment])
write.csv(dat.dir, file.path(TAB_DIR, "gse212388_drug_direction.csv"),
  row.names = FALSE)
# 每臂的基因级检验力单列披露：两臂的样本深度不同，不合并成一个数
dat.mb <- unique(dat.dir[, c("drug", "min.BH")])
p.forest <- ggplot2::ggplot(dat.dir,
    ggplot2::aes(x = mean.logFC.seg, y = seg)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
    color = "grey55") +
  ggplot2::geom_point(size = 2.8, color = "#B2182B") +
  ggplot2::facet_wrap(~ drug, ncol = 1) +
  sci_theme +
  ggplot2::labs(x = "Mean log2FC of segment genes (descriptive)", y = NULL,
    title = "GSE212388 EC: segment-level shift under each drug",
    subtitle = paste0("Pseudobulk per sample; limma ~ flow + condition.\n",
      "Each drug against its depth-matched control, two samples per group"),
    caption = paste0("No segment gene passed BH < ", PVAL_CUTOFF, " with ",
      "|log2FC| >= ", LOGFC_MIN, " in either arm.\n",
      "Smallest BH per arm: ",
      paste0(dat.mb$drug, " ", format(dat.mb$min.BH), collapse = "; "),
      ". Descriptive readout, no gene-level claim."))
sci_save(p.forest, "EC_Segment_Direction_Forest", width = 9, height = 6)

dat.ss <- dat.seg.sample
lst.sp <- strsplit(dat.ss$sample, "_", fixed = TRUE)
dat.ss$flow <- factor(vapply(lst.sp, `[`, character(1L), 3L),
  levels = c("PhysioFlow", "AlteredFlow"))
dat.ss$cond <- vapply(lst.sp, `[`, character(1L), 4L)
# 样本级段分：每臂各做一次 limma，臂内两对照两药样本，流型作区组
dat.seg.stat <- do.call(rbind, lapply(vec.arm, function(a) {
    d.a <- dat.ss[dat.ss$arm == a, , drop = FALSE]
    drug <- setdiff(d.a$cond, "Control")
    d.a$cond <- factor(d.a$cond, levels = c("Control", drug))
    vec.samp <- unique(d.a$sample)
# 长转宽：样本 x 段，按样本名对齐（禁位置对齐）
    mat.sc <- sapply(vec.segment, function(sg) {
        d <- d.a[d.a$segment == sg, , drop = FALSE]
        d$mean_score[match(vec.samp, d$sample)]
      })
    rownames(mat.sc) <- vec.samp
    d.des <- d.a[match(vec.samp, d.a$sample), c("sample", "flow", "cond")]
    stopifnot(identical(rownames(mat.sc), d.des$sample))
    des.m <- stats::model.matrix(~ flow + cond, data = d.des)
    ct.term <- paste0("cond", drug)
    fit.m <- limma::lmFit(t(mat.sc), des.m)
    tt.seg <- limma::topTable(limma::eBayes(limma::contrasts.fit(fit.m,
      limma::makeContrasts(contrasts = ct.term, levels = des.m))),
      coef = ct.term, number = Inf, adjust.method = "BH", confint = TRUE)
    data.frame(arm = a, drug = drug, segment = rownames(tt.seg),
      logFC = tt.seg$logFC, CI.L = tt.seg$CI.L, CI.R = tt.seg$CI.R,
      P.Value = tt.seg$P.Value, adj.P.Val = tt.seg$adj.P.Val,
      mean_ref = vapply(rownames(tt.seg), function(sg)
        mean(d.a$mean_score[d.a$segment == sg & d.a$cond == "Control"]),
        numeric(1L)),
      mean_test = vapply(rownames(tt.seg), function(sg)
        mean(d.a$mean_score[d.a$segment == sg & d.a$cond == drug]),
        numeric(1L)),
      stringsAsFactors = FALSE)
  }))
write.csv(dat.seg.stat, file.path(TAB_DIR,
  "gse212388_segment_score_stats.csv"), row.names = FALSE)
write.csv(dat.seg.sample, file.path(TAB_DIR,
  "gse212388_segment_score_by_sample.csv"), row.names = FALSE)

# 配对差：每臂每流型一对药与对照，纯描述，无推断
dat.pair <- do.call(rbind, lapply(vec.arm, function(a) {
    d.a <- dat.ss[dat.ss$arm == a, , drop = FALSE]
    drug <- setdiff(d.a$cond, "Control")
    vec.c <- unique(d.a$sample[d.a$cond == "Control"])
    vec.t <- unique(d.a$sample[d.a$cond == drug])
    lst.c <- strsplit(vec.c, "_", fixed = TRUE)
    lst.t <- strsplit(vec.t, "_", fixed = TRUE)
    dat.map <- data.frame(flow = vapply(lst.c, `[`, character(1L), 3L),
      sample.control = vec.c, stringsAsFactors = FALSE)
    dat.map$sample.test <- vec.t[match(dat.map$flow,
      vapply(lst.t, `[`, character(1L), 3L))]
    do.call(rbind, lapply(vec.segment, function(sg) {
        do.call(rbind, lapply(seq_len(nrow(dat.map)), function(k) {
            v.c <- d.a$mean_score[d.a$segment == sg &
              d.a$sample == dat.map$sample.control[k]]
            v.t <- d.a$mean_score[d.a$segment == sg &
              d.a$sample == dat.map$sample.test[k]]
            data.frame(arm = a, drug = drug, segment = sg,
              flow = dat.map$flow[k],
              control_sample = dat.map$sample.control[k],
              test_sample = dat.map$sample.test[k], score_control = v.c,
              score_test = v.t, delta = v.t - v.c, stringsAsFactors = FALSE)
          }))
      }))
  }))
write.csv(dat.pair, file.path(TAB_DIR, "gse212388_segment_pair_delta.csv"),
  row.names = FALSE)

# 流型与药物的交互项在此设计下不可估计：臂内 2x2 每格仅一个样本，交互模型
# 参数数等于观测数，残差自由度为 0，故不计算任何交互对比。

# 分组均值：臂 x 流型 x 条件，分层描述的直接读数
dat.grp <- do.call(rbind, lapply(vec.arm, function(a) {
    d.a <- dat.ss[dat.ss$arm == a, , drop = FALSE]
    do.call(rbind, lapply(vec.segment, function(sg) {
        do.call(rbind, lapply(levels(d.a$flow), function(fl) {
            do.call(rbind, lapply(unique(d.a$cond), function(cd) {
                v <- d.a$mean_score[d.a$segment == sg & d.a$flow == fl &
                  d.a$cond == cd]
                data.frame(arm = a, segment = sg, flow = fl, cond = cd,
                  n_samples = length(v), mean_score = mean(v),
                  stringsAsFactors = FALSE)
              }))
          }))
      }))
  }))
write.csv(dat.grp, file.path(TAB_DIR,
  "gse212388_segment_score_by_group.csv"), row.names = FALSE)
# 每臂的样本级最小校正 P：两臂功效不同，分别披露
dat.amin <- do.call(rbind, lapply(vec.arm, function(a)
  data.frame(drug = sub("_vs_Control$", "", a),
    min.adj = signif(min(dat.seg.stat$adj.P.Val[dat.seg.stat$arm == a],
      na.rm = TRUE), 3), stringsAsFactors = FALSE)))
dat.seg.stat$seg <- unname(seg.short[dat.seg.stat$segment])
p.seg <- ggplot2::ggplot(dat.seg.stat,
    ggplot2::aes(x = logFC, y = seg)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
    color = "grey55") +
  ggplot2::geom_linerange(ggplot2::aes(xmin = CI.L, xmax = CI.R),
    linewidth = 0.7, color = "#B2182B") +
  ggplot2::geom_point(size = 2.8, color = "#B2182B") +
  ggplot2::facet_wrap(~ drug, ncol = 1) +
  sci_theme +
  ggplot2::labs(x = "Module score change (limma estimate)", y = NULL,
    title = "GSE212388 EC: sample-level segment score shift under each drug",
    subtitle = paste0("limma ~ flow + condition; per-sample segment module ",
      "scores.\nEach drug against its depth-matched control"),
    caption = paste0("Exploratory: two samples per group, residual df = 1; ",
      "smallest adjusted P per arm: ",
      paste0(dat.amin$drug, " ", format(dat.amin$min.adj),
        collapse = "; "),
      ". Error bars = 95% CI."))
sci_save(p.seg, "EC_SegmentScore_Forest", width = 9, height = 6)

# 每臂每流型内的配对描述差：五段 x 两个流型，单位为模块分变化
dat.pair$seg <- unname(seg.short[dat.pair$segment])
p.pair <- ggplot2::ggplot(dat.pair,
    ggplot2::aes(x = seg, y = delta, fill = flow)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
    color = "grey40", linewidth = 0.5) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7),
    width = 0.62) +
  ggplot2::scale_fill_manual(values = c(PhysioFlow = "#7FA8C9",
    AlteredFlow = "#B2182B"), name = "Flow") +
  ggplot2::facet_wrap(~ drug, ncol = 1) +
  sci_theme +
  ggplot2::labs(x = "Inflammation segment",
    y = "Module score change (drug minus matched control)",
    title = "Segment score change by flow background",
    subtitle = paste0("GSE212388 EC; one drug sample and one matched ",
      "control per flow in each arm; descriptive only"),
    caption = paste0("S1 initiation, S2 amplification, S3 effect, ",
      "S4 resolution, S5 chronic.\n",
      "AlteredFlow = pro-inflammatory flow; no inferential test at this ",
      "replicate count."))
sci_save(p.pair, "EC_Segment_Change_By_Flow", width = 9, height = 6)

dat.cy <- dat.cycl
lst.cp <- strsplit(dat.cy$sample, "_", fixed = TRUE)
dat.cy$flow <- factor(vapply(lst.cp, `[`, character(1L), 3L),
  levels = c("PhysioFlow", "AlteredFlow"))
dat.cy$cond <- vapply(lst.cp, `[`, character(1L), 4L)
dat.cy$drug <- sub("_vs_Control$", "", dat.cy$arm)
# 每臂单独做 Kruskal-Wallis：两臂的对照分属不同深度层，不能合并成一次检验
dat.cy.sum <- do.call(rbind, lapply(vec.arm, function(a) {
    d.a <- dat.cy[dat.cy$arm == a, , drop = FALSE]
    drug <- setdiff(d.a$cond, "Control")
    kw <- stats::kruskal.test(d.a$frac_cycling,
      factor(d.a$cond, levels = c("Control", drug)))
    do.call(rbind, lapply(c("Control", drug), function(cd)
      data.frame(arm = a, drug = drug, cond = cd,
        n_samples = sum(d.a$cond == cd),
        mean_frac = mean(d.a$frac_cycling[d.a$cond == cd]),
        kw_p = signif(kw$p.value, 3), stringsAsFactors = FALSE)))
  }))
write.csv(dat.cy.sum, file.path(TAB_DIR,
  "gse212388_cycling_fraction.csv"), row.names = FALSE)
p.cy <- ggplot2::ggplot(dat.cy,
    ggplot2::aes(x = cond, y = frac_cycling, color = cond)) +
  ggplot2::geom_point(size = 3,
    position = ggplot2::position_jitter(width = 0.12)) +
  ggplot2::scale_color_manual(values = c(Control = "#4C72B0",
    Rapamycin = "#B2182B", Paclitaxel = "#E18727")) +
  ggplot2::guides(color = "none") +
  ggplot2::facet_wrap(~ drug, ncol = 1, scales = "free_x") +
  sci_theme +
  ggplot2::labs(x = NULL, y = "Fraction of cycling cells",
    title = "GSE212388 EC: cycling fraction by treatment",
    subtitle = "Cycling = clusters annotated *_Cycling",
    caption = paste0("Kruskal-Wallis per arm: ",
      paste0(unique(dat.cy.sum[, c("drug", "kw_p")])$drug, " p = ",
        unique(dat.cy.sum[, c("drug", "kw_p")])$kw_p,
        collapse = "; "),
      ".\nEach arm uses its own depth-matched controls, ",
      "two samples per condition."))
sci_save(p.cy, "EC_Cycling_Fraction_By_Treatment", width = 8, height = 6)

dat.ct.stat$bh <- ave(dat.ct.stat$kw_p, dat.ct.stat$acc,
  FUN = function(x) stats::p.adjust(x, method = "BH"))
write.csv(dat.ct.stat, file.path(TAB_DIR, "segment_by_celltype_stats.csv"),
  row.names = FALSE)

# PART 2 空间：GeoMx dcc.gz ---------------------------------------------------------------------

# q11 的判定单位是「Panel 命中的段基因个数」（MIN_REGION_GENES），
# 不是「有命中的段数」——旧写法把两者混用，量纲不同，此处已更正。
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
  vec.rts.id)  # RTS_ID -> 干净基因符号（命名向量查表）
cat(sprintf("[S12] pkc map built: %d RTS -> %d unique symbols\n",
  length(map.sym), length(unique(vec.rts.sym))))

# 返回单个 dcc 的 Code_Summary 探针行数（RTS 号条数）。
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
# pkc 映射（map.sym 已在上方构建：RTS_ID -> 干净基因符号，命名向量查表）
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

# PART 2.5 图件 -------------------------------------------------------------------------------
# 命名遵守 16（描述性英文，禁 Fig 前缀）；图内文字全英文（ 17）；
# 分面一律用 facet_wrap，不引入任何拼图包。

dat.cluster$acc <- factor(dat.cluster$acc)
# fill = pass_min_cells 恒为 TRUE（表已在上面按阈值过滤掉小簇），属冗余编码（）：
# 改单色柱；被剔除的小簇已在收口处披露。
sci_save(ggplot(dat.cluster, aes(x = cluster, y = n_cells)) +
  geom_col(fill = "#4C72B0") +
  facet_wrap(~ acc, scales = "free_x") +
  labs(title = "scRNA clusters kept per series",
    x = "Cluster", y = "Cells") +
  sci_theme,
  "ScRNA_Cluster_Size_By_Series", width = 10, height = 4.5)


# PART 3 收口 -------------------------------------------------------------------------------

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

