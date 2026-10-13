# ===============================================================
# 06b_PLAUR_SPP1_pseudobulk.R
# 功能: 4 目标基因的区室 pseudobulk 差异分析(sum 聚合 raw counts + edgeR)
# 定位: 替代 06 中 spot 水平 SCT log1p 箱线/检验; 规避检测下限离散化
#       (SCT data 在 PLAU/FBXL12 上退化为 0/0.693/1.099 三档), 用负二项
#       模型在聚合计数上给出可解释的 log2FC
# 依据: Nat Commun 2023 基准推荐 Pseudobulk_edgeR; 本数据 n=2 患者,
#       残差 df=1, 故以 log2FC 效应量为主证据, p 值仅作参考
# 可视化: PseudobulkLog2FC.pdf (各比较的 log2FC 效应量)
# 数据: 输入 results/01_read_qc/table/ST-*_seurat.rds
#            results/06_PLAUR_SPP1_signature/table/signature_master.rds
#       输出 results/06b_PLAUR_SPP1_pseudobulk/table/PseudobulkCounts.csv
#            results/06b_PLAUR_SPP1_pseudobulk/table/PseudobulkLog2FC.csv
# ===============================================================

suppressMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
})

# --- 全局辅助函数(内联, 规范22 自给自足, 禁 source()) -------------------
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
  pdf(file.path(FIG_DIR, paste0(filename, ".pdf")), width = width,
      height = height, family = "Helvetica")
  print(plot_obj)
  dev.off()
}

# --- 目录结构 -----------------------------------------------------------
BASE_DIR <- "X:/XX"
TAB01    <- file.path(BASE_DIR, "results/01_read_qc/table")
TAB06    <- file.path(BASE_DIR, "results/06_PLAUR_SPP1_signature/table")
FIG_DIR  <- file.path(BASE_DIR, "results/06b_PLAUR_SPP1_pseudobulk/Figure")
TAB_DIR  <- file.path(BASE_DIR, "results/06b_PLAUR_SPP1_pseudobulk/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 分析常量(与 06 一致, 不得改动) ------------------------------------
SAMPLE_ID <- c("ST-P1_liver", "ST-P1_colon", "ST-P2_liver", "ST-P2_colon")
GENE_SET  <- c("FBXL12", "PLAUR", "PLAU", "SPP1")
COMP_ORDER <- c("Adjacent colon", "PT", "Adjacent liver", "LM")
# 每个比较: 名称 -> 参与的两区室(前者为基线, 后者为对照方向末位)
#   log2FC 定义为 log2(区室2 / 区室1), 即 comp_b 相对 comp_a
CMP <- list(
  list(name = "PT_vs_AdjColon",  tissue = "colon", a = "Adjacent colon", b = "PT"),
  list(name = "LM_vs_AdjLiver",  tissue = "liver", a = "Adjacent liver", b = "LM"),
  list(name = "LM_vs_PT",        tissue = "cross", a = "PT",             b = "LM")
)
# 掩码与基因的对应(沿用 06 的签名定义): FBXL12/PLAU 限定肿瘤,
#   PLAUR/SPP1 限定巨噬; "all" = 区室内全部 spot
GENE_MASK <- c(FBXL12 = "tumor", PLAU = "tumor",
               PLAUR  = "macro", SPP1 = "macro")

# =============================================================================
# PART 1: 构建 pseudobulk 计数矩阵(sample × 区室 × 掩码)
# 目的: 把 spot 级 raw counts 聚合到区室级, 脱离检测下限
# 方法: Spatial assay counts 层按 spot 子集 sum 聚合; 三种掩码
#       all(全部 spot) / tumor(is_tumor) / macro(is_macro)
# 数据: 输出 PseudobulkCounts.csv(含 n_spots 便于审计)
# =============================================================================

cat("========== PART 1: Build pseudobulk counts ==========\n")

master <- readRDS(file.path(TAB06, "signature_master.rds"))
master <- master[master$sample %in% SAMPLE_ID, ]

pb.list   <- list()   # 计数向量
meta.list <- list()   # 分组元信息

for (s in SAMPLE_ID) {
  st  <- readRDS(file.path(TAB01, paste0(s, "_seurat.rds")))
  cnt <- GetAssayData(st, assay = "Spatial", layer = "counts")
  md  <- master[master$sample == s, ]
  idx <- match(md$spot, colnames(st))
  stopifnot(!anyNA(idx))

  for (cp in COMP_ORDER) {
    in.cp <- which(md$compartment == cp)
    if (length(in.cp) == 0) next
    for (mk in c("all", "tumor", "macro")) {
      keep <- switch(mk,
        all   = in.cp,
        tumor = in.cp[md$is_tumor[in.cp]],
        macro = in.cp[md$is_macro[in.cp]])
      if (length(keep) == 0) next
      grp <- paste(s, cp, mk, sep = "|")
      pb.list[[grp]] <- Matrix::rowSums(cnt[, idx[keep]])
      meta.list[[grp]] <- data.frame(
        group = grp, sample = s, compartment = cp, mask = mk,
        n_spots = length(keep),
        lib_size = sum(Matrix::rowSums(cnt[, idx[keep]])),
        stringsAsFactors = FALSE)
    }
  }
}

pb.meta <- do.call(rbind, meta.list)
rownames(pb.meta) <- NULL
pb.meta$patient <- sub("^ST-(P[0-9]).*$", "\\1", pb.meta$sample)
pb.meta$tissue  <- ifelse(grepl("liver", pb.meta$sample), "liver", "colon")

pb.mat <- do.call(cbind, pb.list)
pb.mat <- pb.mat[, pb.meta$group]

cat("pseudobulk 分组数:", nrow(pb.meta), "\n")
print(pb.meta[, c("group", "compartment", "mask", "n_spots", "lib_size")])

# 目标基因的聚合计数与检出 spot 数(审计用)
for (g in intersect(GENE_SET, rownames(pb.mat))) {
  pb.meta[[paste0("cnt_", g)]] <- as.integer(pb.mat[g, pb.meta$group])
}
write.csv(pb.meta, file.path(TAB_DIR, "PseudobulkCounts.csv"), row.names = FALSE)

# =============================================================================
# PART 2: edgeR 差异分析(log2FC 效应量为主, p 值仅参考)
# 目的: 在区室级聚合计数上估计 log2FC
# 方法: TMM 归一化 -> estimateDisp -> glmFit -> glmLRT; predFC 提供
#       低计数的收缩 log2FC(prior.count), 与 DESeq2 lfcShrink 同旨
# 数据: 输出 PseudobulkLog2FC.csv
# =============================================================================

cat("========== PART 2: edgeR differential expression ==========\n")

gene.idx <- match(GENE_SET, rownames(pb.mat))
stopifnot(!anyNA(gene.idx))

de.list <- list()

for (cm in CMP) {
  for (mk in c("all", "tumor", "macro")) {
    # 取出该比较两侧、该掩码下的分组
    sub <- pb.meta[pb.meta$mask == mk & pb.meta$compartment %in% c(cm$a, cm$b), ]
    if (nrow(sub) < 2) next
    if (length(unique(sub$compartment)) < 2) next
    if (length(unique(sub$patient)) < 2) next

    y <- DGEList(counts = pb.mat[, sub$group],
                 group  = factor(sub$compartment, levels = c(cm$a, cm$b)))
    # 过滤全零/近零基因以保证离散度可估; 阈值取总计数 >= 5,
    #   远低于 4 目标基因的聚合量级(数百 UMI), 确保目标基因不被剔除
    #   (filterByExpr 默认 min.count=10 会误删贴近检测下限的目标基因)
    keepg <- Matrix::rowSums(y$counts) >= 5
    y <- y[keepg, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y, method = "TMM")

    patient <- factor(sub$patient)
    comp    <- factor(sub$compartment, levels = c(cm$a, cm$b))
    design  <- model.matrix(~ patient + comp)

    y <- estimateDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef = ncol(design))

    # 收缩 log2FC: 低计数基因向 0 收缩, 避免 1 UMI 差异被放大
    pfc <- predFC(y, design, prior.count = 1, dispersion = fit$dispersion)

    tt <- topTags(lrt, n = Inf, sort.by = "none")$table
    for (g in GENE_SET) {
      if (!g %in% rownames(tt)) next
      k <- which(rownames(tt) == g)
      de.list[[length(de.list) + 1]] <- data.frame(
        comparison = cm$name, mask = mk, gene = g,
        base = cm$a, target = cm$b,
        n_group = nrow(sub),
        n_spots_min = min(sub$n_spots), n_spots_max = max(sub$n_spots),
        log2FC = as.numeric(tt$logFC[k]),
        log2FC_shrunk = as.numeric(pfc[k, ncol(design)]),
        logCPM = as.numeric(tt$logCPM[k]),
        p_value = as.numeric(tt$PValue[k]),
        fdr = as.numeric(tt$FDR[k]),
        stringsAsFactors = FALSE)
    }
  }
}

de.tab <- do.call(rbind, de.list)
# 残差 df: 分组数 - (截距 + patient + comp)
de.tab$resid_df <- de.tab$n_group - 3
write.csv(de.tab, file.path(TAB_DIR, "PseudobulkLog2FC.csv"), row.names = FALSE)
print(de.tab)

# =============================================================================
# PART 3: log2FC 效应量可视化
# 目的: 以效应量为主证据呈现区室间差异
# 方法: 分组柱状图, 标注 p 值仅供参考(df=1)
# 数据: 输出 PseudobulkLog2FC.pdf
# =============================================================================

cat("========== PART 3: Plot log2FC effect sizes ==========\n")

d <- de.tab[de.tab$mask == "all", ]
d$gene <- factor(d$gene, levels = GENE_SET)
d$comparison <- factor(d$comparison,
                       levels = c("PT_vs_AdjColon", "LM_vs_AdjLiver", "LM_vs_PT"))

p.fc <- ggplot(d, aes(x = gene, y = log2FC_shrunk, fill = comparison)) +
  geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.4,
             color = "grey40") +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  labs(title = "Pseudobulk log2FC by Compartment (all spots)",
       subtitle = paste0("edgeR shrunken log2FC (prior.count=1); ",
                         "n=2 patients, residual df=1, p-values indicative only"),
       x = "Gene", y = "log2 fold change (target vs base)",
       fill = "Comparison") +
  sci_theme
sci_save(p.fc, "PseudobulkLog2FC", width = 9, height = 6)

cat("06b_PLAUR_SPP1_pseudobulk.R done: PseudobulkCounts.csv + ",
    "PseudobulkLog2FC.csv + PseudobulkLog2FC.pdf\n")
