# ================================================================
# 01_GSE178318_Preprocessing_QC.R
# 功能: GSE178318 单细胞数据预处理与质量控制
# 定位: GSE178318 scRNA-seq 分析流程 Step 1
# 可视化: QC 小提琴图、UMAP 聚类/组织/患者、DotPlot 双参考注释
# ================================================================
# 数据来源: Che et al., Cell Discovery 2021, PMID 33510120
        
        
        
        
        
        
# 患者: COL07/12/16 (treatment-naïve only)
# 组织: PT (原发灶) / LM (肝转移灶)
# ================================================================
# 方法论: 遵循原文 Seurat FindClusters + 标记基因注释
# 优化: DotPlot 双参考验证 + InferCNV 免疫参考细胞
# 断点续跑: 各 PART 检查对应 checkpoint，存在则跳过
# ================================================================

# --- 包加载 ------------------------------------------------------

library(dplyr)        # 数据操作
library(Matrix)       # 稀疏矩阵读取
library(ggplot2)      # 数据可视化
library(patchwork)    # 图形拼接
library(viridis)      # 色盲友好连续配色
library(Seurat)      # 单细胞分析核心包 (v5)
library(harmony)      # 批次效应校正
library(infercnv)    # 拷贝数变异分析
library(SingleR) 
options(warn = 1)
set.seed(12345)   # 全局随机种子：UMAP / 聚类 / 降采样可复现（P2-⑬）

# --- 全局常量 ----------------------------------------------------

nfeature_min     <- 500L    # 原文阈值：基因数 < 500 剔除
percent_mt_max  <- 15L     # 原文：线粒体 > 15% 剔除
n_hvg           <- 1511L   # 原文：HVG 由 cutoff(mean∈[0.0125,3] 且 dispersion>0.5) 选出，结果 1511 个；此处作 nfeatures 上限（P2-⑧）
npcs_total      <- 50L     # PCA 计算 50 PCs（供肘部图查看）
npcs_use        <- 50L     # 原文：FindClusters on 50 PCs（与 npcs_total 一致）
cluster_res     <- 1.0     # 原文：resolution = 1.0

# --- 流程参数（集中管理）--------------------
ncores_nonwin    <- 4L      # 非 Windows 平台并行核数（mclapply / InferCNV 线程）
df_nfeatures     <- 3000L   # DoubletFinder 内部 HVG 数（官方推荐 2000-3000）
df_pN            <- 0.25    # DoubletFinder 标准扰动倍数（官方推荐）
df_chemistry     <- "10xV2" # 10x 化学版本：原文 Methods Page 16 明确用 Chromium 3' Reagent Kits v2，故默认 V2（inv=1.74）；V3 inv=1.61（见 ExpectedDoubletNExp）
infercnv_threads <- 4L      # InferCNV 计算线程数
subcluster_res   <- 1.0     # PART 8 亚型细分 resolution（复现原文 same protocol = res 1.0，主聚类用 cluster_res）
# 最小样本量门控
min_cells_doublet <- 50L    # 患者细胞数低于此则 DoubletFinder 不可靠，默认 Singlet
min_epcam         <- 10L    # 患者上皮观测细胞下限
min_cells_heatmap <- 10L    # 热图可用细胞下限
min_cells_subtype <- 50L    # 亚型细分细胞下限
# 降采样种子基准（PART 7b，固定可复现）
seed_mal_base    <- 2021L   # 恶性细胞降采样种子基数
seed_ref_base    <- 3072L   # 参考细胞降采样种子基数
# PART 7b 热图降采样参数（集中管理，避免散落魔法数字）
heatmap_max_mal     <- 1000L  # 热图恶性细胞最多取 1000
heatmap_max_ref     <- 500L   # 热图参考细胞最多取 500
heatmap_n_top_genes <- 2000L  # 热图选取 top 2000 高变基因（局部别名，避免遮蔽全局 n_hvg）

naive_patients  <- c("COL07", "COL12", "COL16")

# 上皮观测/参照配置（PART 7 恶性鉴定仅在上皮区室内进行）
epi_types <- c("Epithelial/Tumor")   # 上皮 celltype 标签（PART6 paper 注释 + Table3 确认）
# InferCNV 参照 = 本患者所有非上皮细胞（免疫/基质，确定二倍体）。分层下采样上限：
#   参照过多会让 leiden 子聚类与 parallelDist 距离矩阵(O(n^2)) 显著变慢；
#   1500 已足以稳定估计 hspike（HMM 状态分布经验基准，由参照细胞构建，
#   见 infercnv:::.build_and_add_hspike）。COL12 实测 1502 参照 → 7 分钟跑完。
infercnv_ref_cap <- 1500L

# --- 目录结构定义 ------------------------------------------------

data_dir <- "X:/CRLM/raw"

# GENCODE 人类注释 GTF（GRCh38, release 44）
# 用于从真实染色体坐标构建 gene_order.txt（InferCNV 必需，禁止伪造坐标）
# 下载: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gencode_gtf <- file.path(data_dir, "gencode.v44.annotation.gtf.gz")

fig_dir  <- "E:/CRLM/scRNA/GSE178318/results/01_QC/Figure"
tab_dir  <- "E:/CRLM/scRNA/GSE178318/results/01_QC/table"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

checkpoint_dir <- tab_dir

# 组织配色（PT = 原发灶，LM = 肝转移）
tissue_colors <- c("PT" = "#E69F00", "LM" = "#56B4E9")

# --- 全局辅助函数 ------------------------------------------------

FlushLog <- function(msg) {
  cat(msg)
  flush.console()
}

MemReport <- function(label = "") {
  mem <- gc()                 # gc() 返回当前内存占用表（同时触发回收）
  used.mb <- sum(mem[, 2])
  FlushLog(sprintf("[内存%s] 当前占用: %.0f MB\n", label, used.mb))
}

# 保证 seurat.obj 含下游绘图所需的全部降维层（pca / harmony / umap）及其全部
# 前置派生层（data / scale.data / VariableFeatures）。这是 resume / 各 PART 跳过
# 分支的【单一汇聚点】保证。
# 关键根因（实读 checkpoint_part035_doublet.rds 坐实）：降维嵌入依赖 PART 4/5 产出的
# 派生层——pca 依赖 scale.data（ScaleData 产物），scale.data 又依赖 data
# （NormalizeData 产物）。若 resume 落到 PART 4 之前的对象（如 part035：data /
# scale.data 均为空壳、无 VariableFeatures），仅补 PCA 会因缺 scale.data 崩
# （Seurat v5 PrepDR5 报 “No layer matching pattern 'scale.data'”）。故此处从始终
# 存在的原始 counts 出发，完整重建整条派生链（每层幂等，参数与 PART 4/5 完全一致、
# 不引入偏差）：
#   NormalizeData(LogNormalize, sf=1e4) → FindVariableFeatures(dispersion, n_hvg)
#   → ScaleData(regress nCount_RNA/percent.mt) → RunPCA(npcs=npcs_total)
#   → RunHarmony(group.by.vars=patient) → RunUMAP(reduction=harmony, dims=1:npcs_use)
# 下游依赖：PART 6/7.5/7b 用 umap，PART 8 用 harmony。完整对象此处为 no-op。
# 注：LayerData 对“存在但空”的层返回非空对象，故用 dim()>0 判别真正可用层。
ensure_reductions <- function(obj) {
  if (is.null(obj)) return(obj)
  da <- DefaultAssay(obj)
  layer_ok <- function(nm) {      # 层存在且非“空壳”
    m <- tryCatch(LayerData(obj, assay = da, layer = nm), error = function(e) NULL)
    !is.null(m) && all(dim(m) > 0L)
  }
  # 1) 归一化层（RunPCA/ScaleData 的前置，PART 4 产物）
  if (!layer_ok("data")) {
    FlushLog("    data 层缺失/空，补算 NormalizeData（LogNormalize, sf=10000）\n")
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  # 2) 高变基因（RunPCA 默认用 VariableFeatures）
  if (length(VariableFeatures(obj)) == 0L) {
    FlushLog(sprintf("    VariableFeatures 缺失，补算 FindVariableFeatures（dispersion, n=%d）\n", n_hvg))
    obj <- FindVariableFeatures(obj, selection.method = "dispersion", nfeatures = n_hvg, verbose = FALSE)
  }
  # 3) 缩放层（RunPCA 依赖 scale.data，PART 4 产物）
  #    按 PART 4 原参（回归 nCount_RNA/percent.mt）补算；失败即报错（对象应已含这两个变量）。
  if (!layer_ok("scale.data")) {
    FlushLog("    scale.data 层缺失/空，补算 ScaleData（regress nCount_RNA/percent.mt）\n")
    obj <- ScaleData(obj, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
  }
  # 4) 降维嵌入：pca → harmony → umap（参数与 PART 4/5 完全一致）
  if (!"pca" %in% Reductions(obj)) {
    FlushLog(sprintf("    pca 缺失，补算 PCA（npcs=%d）\n", npcs_total))
    obj <- RunPCA(obj, npcs = npcs_total, verbose = FALSE)
  }
  if (!"harmony" %in% Reductions(obj)) {
    FlushLog("    harmony 缺失，从 pca 补算 Harmony（group.by.vars=patient）\n")
    obj <- RunHarmony(obj, "patient", reduction.save = "harmony")
  }
  if (!"umap" %in% Reductions(obj)) {
    FlushLog(sprintf("    从 harmony 重建 UMAP（dims=1:%d）\n", npcs_use))
    obj <- RunUMAP(obj, reduction = "harmony", dims = 1:npcs_use)
  }
  obj
}

# 科学绘图主题
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

# 保存图件为 PDF
sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf_file <- file.path(fig_dir, paste(filename, "pdf", sep = "."))
  ok <- tryCatch({
    pdf(pdf_file, width = width, height = height)
    print(plot_obj)
    TRUE
  }, error = function(e) {
    FlushLog(sprintf("  保存PDF失败 [%s]: %s\n", filename, e$message))
    FALSE
  }, finally = {
    if (dev.cur() > 1) dev.off()
  })
  # 出错时删除残缺/空白 PDF，让"缺失文件"显式暴露失败，而非伪装成成功空白图
  if (!ok && file.exists(pdf_file)) file.remove(pdf_file)
  invisible(NULL)
}

# 带聚类中心标签的 UMAP 图（DimPlot + ggrepel 避让；ggrepel 缺失时退回 geom_text）。
# 抽自 PART 6 细胞类型 UMAP 与 SingleR UMAP 两段同构逻辑，消除复制粘贴。
umap_labeled_plot <- function(obj, group_col, cols, fname, width = 10, height = 7) {
  p.base <- DimPlot(obj, reduction = "umap",
                    group.by = group_col, cols = cols, label = FALSE) +
    theme(aspect.ratio = 1) + sci_theme
  coord <- Embeddings(obj, reduction = "umap")
  label_vec <- obj[[group_col]][[1L]]   # Seurat v5: obj[[col]] 返回单列 data.frame，[[1L]] 取向量；等价于原代码 obj$col，且不丢 NA 行
  label_df <- data.frame(
    UMAP_1 = coord[, 1],
    UMAP_2 = coord[, 2],
    label   = label_vec,
    stringsAsFactors = FALSE
  )
  centroids <- stats::aggregate(
    . ~ label, data = label_df[, c("UMAP_1", "UMAP_2", "label")],
    FUN = median)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p.out <- p.base +
      ggrepel::geom_text_repel(
        data  = centroids,
        aes(x = UMAP_1, y = UMAP_2, label = label),
        size        = 3.5,
        box.padding = 0.4,
        max.overlaps = Inf,
        segment.color = "grey30",
        show.legend  = FALSE
      )
  } else {
    p.out <- p.base +
      geom_text(data = centroids,
        aes(x = UMAP_1, y = UMAP_2, label = label),
        size = 3.5, show.legend = FALSE)
  }
  sci_save(p.out, fname, width = width, height = height)
  invisible(p.out)
}

# 双细胞决策：DoubletFinder 原生率驱动（nExp = 期望双细胞数，见下 ExpectedDoubletNExp）
#   pANN 降序前 nExp 个判为 Doublet（源码 doubletFinder() 纯率驱动，与 pANN 分布形态无关）。
# 诊断仅保留期望率参考线（PART 3.5 的 cut_exp = pANN 降序第 nExp 位），不再做 pANN 谷检测
#   （旧方案已废弃：谷检测对双峰不明显的捕获会误判，且决策本就不依赖 pANN 形态）。
# ---- 期望双细胞数（DoubletFinder 原生率驱动分类所需的 nExp）----
# 依据（非自创，可复核）：
#   · 10x Genomics User Guide CG000592（Next GEM 3' v3.1）多重率表：~1000恢复→0.8%、~10000→7.6%
#   · BimberLab/cellhashR::EstimateMultipletRate（生产级实现）：
#       numCellsLoaded = inv * numCellsRecovered;  multipletRate = 4.597701e-6 * numCellsLoaded
#       10xV3: inv = 1.61 ; 10xV2: inv = 1.74
#   · DoubletFinder 设计要求 nExp = 期望双细胞数（McGinnis 2019, Cell Systems 8(4):329; README FAQ "anticipated doublet rate"）；
#     源码 doubletFinder() L130 = pANN 降序前 nExp 个判为 Doublet —— 纯率驱动，与 pANN 分布形态无关。
#   本数据按患者拆分跑 DF，期望率随各患者恢复细胞数自适应（README: rate depends on number of input cells），
#   故 COL07(26k) 期望率远高于 COL12(12k)，杜绝"大样本静默判 0%"。
ExpectedDoubletNExp <- function(n_recovered, chemistry = df_chemistry) {
  inv <- if (chemistry == "10xV2") 1.74 else 1.61   # 仅 V2/V3 两档；未知默认 V3
  m   <- 4.597701e-06
  rate <- m * inv * n_recovered
  rate <- min(rate, 0.25)          # 安全上限：单泳道极深加载也罕见超过 25%
  max(1L, round(rate * n_recovered))
}

# 从 GENCODE GTF 构建 InferCNV 所需的 gene_order.txt（真实染色体坐标）
# 仅运行一次：若输出文件已存在则跳过（符合本流程 checkpoint 哲学）
# 输出格式（官方 wiki）：gene_symbol \t chr \t start \t end，无表头
BuildGeneOrderFromGTF <- function(gtf_path, out_path, chunk = 500000L) {
  if (file.exists(out_path)) {
    FlushLog(sprintf("  gene_order.txt 已存在，跳过构建 [%s]\n",
                     basename(out_path)))
    return(invisible(NULL))
  }
  if (!file.exists(gtf_path)) {
    stop("GENCODE GTF 不存在: ", gtf_path,
         "\n  请先从以下地址下载并放到 data_dir：",
         "\n  https://ftp.ebi.ac.uk/pub/databases/gencode/",
         "Gencode_human/release_44/gencode.v44.annotation.gtf.gz")
  }
  FlushLog(sprintf("  从 GTF 构建 gene_order.txt: %s\n", basename(gtf_path)))

  con <- gzfile(gtf_path, "rt")
  g_chr <- c(); g_start <- c(); g_end <- c(); g_name <- c()
  repeat {
    lines <- readLines(con, n = chunk)
    if (length(lines) == 0L) break
    is_gene <- grepl("\tgene\t", lines)
    if (!any(is_gene)) next
    gl <- lines[is_gene]
    f  <- strsplit(gl, "\t")
    g_chr   <- c(g_chr,   vapply(f, `[[`, "", 1L))
    g_start <- c(g_start, vapply(f, function(x) as.integer(x[4L]), 0L))
    g_end   <- c(g_end,   vapply(f, function(x) as.integer(x[5L]), 0L))
    attr    <- vapply(f, `[[`, "", 9L)
    m       <- regexec('gene_name "([^"]+)"', attr)
    hit     <- regmatches(attr, m)
    g_name  <- c(g_name, vapply(hit, function(x)
      if (length(x) >= 2L) x[2L] else "", ""))
  }
  close(con)

  # 仅保留标准染色体（排除 alt / patch / UN 等 scaffold）
  std_chr <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
  keep <- g_chr %in% std_chr
  df <- data.frame(
    gene  = g_name[keep],
    chr   = g_chr[keep],
    start = g_start[keep],
    end   = g_end[keep],
    stringsAsFactors = FALSE
  )
  df <- df[df$gene != "", ]
  df$gene <- toupper(df$gene)        # 先统一大小写，再按 symbol 去重
  df <- df[!duplicated(df$gene), ]  # 去重（修正：原去重在 toupper 前，会漏掉大小写导致的重复 symbol）
  df <- df[order(df$chr, df$start), ]

  write.table(df, out_path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  rm(g_chr, g_start, g_end, g_name, df); gc()
  FlushLog(sprintf("  gene_order.txt 已构建: %d 基因 (真实坐标)\n",
                   nrow(df)))
}

# 从 InferCNV 输出目录提取 HMM 离散态矩阵（基因×细胞，值∈{1..6}）
# 兼容 infercnv 版本差异（关键！）：
#   - 旧版(<~1.14)：HMM 在 step17，文件名 "17_HMM_*.infercnv_obj"，且 run() 会把离散态回写主对象 @expr.data
#   - 1.23.0（canonical）：HMM 在 step18，文件名 "18_HMM_pred*.infercnv_obj"；
#   - 但本机用户实际安装的 infercnv 补丁版少一步，HMM 落在 step17，落盘 "17_HMM_pred*.infercnv_obj"
#     （磁盘实证：infercnv_out_COL07/17_HMM_predHMMi6...infercnv_obj）。helper 优先 18_、回退 17_，故两者皆可用。
#             run() 返回值 @expr.data 始终是 pre-HMM 平滑连续值（HMM 态只活在局部 hmm.infercnv_obj，
#             不回写主对象，且 step21 还会把态转成代表表达量），故必须从此 RDS 取态，不能信 run() 返回值
# 统一策略：无论版本，都从该 HMM-pred RDS 读取 @expr.data（=离散态），避免版本语义差异
get_infercnv_hmm_states <- function(out_dir) {
  if (!dir.exists(out_dir)) stop("InferCNV 输出目录不存在: ", out_dir)
  # 1.23.0: 18_HMM_pred*.infercnv_obj（含 raw/filtered HMM 态）；优先 step 号最大者
  cand <- list.files(out_dir, pattern = "HMM_pred.*\\.infercnv_obj$", full.names = TRUE)
  if (length(cand) == 0L) {
    # 旧版回退: 17_HMM_*.infercnv_obj
    cand <- list.files(out_dir, pattern = "^17_HMM.*\\.infercnv_obj$", full.names = TRUE)
  }
  if (length(cand) == 0L) {
    stop("未找到 InferCNV HMM 结果 RDS（", out_dir, "），请确认 up_to_step 含 HMM 步（1.23.0 需 >=18）")
  }
  # 取 _HMM_pred 紧邻前缀的数字 step 号（兼容「18_HMM_pred...」或「infercnv.18_HMM_pred...」两种命名），取最大者
  cand <- cand[which.max(as.integer(sub(".*?(\\d+)_HMM_pred.*", "\\1", basename(cand))))]  # step 号最大=最完整（18_HMM_pred 优先于 17_HMM_pred，避免旧版残留干扰）
  obj  <- readRDS(cand)
  ed   <- as.matrix(obj@expr.data)
  if (!all(ed %in% 1:6)) {
    stop("HMM 态矩阵非预期离散值(1-6)，疑似取到非 HMM 中间对象: ", basename(cand))
  }
  ed
}

# 基于 InferCNV HMM i6 六态逐细胞判定恶性细胞
# 输入：hmm_states（基因×细胞 的离散 HMM 状态矩阵，值∈{1..6}，由 get_infercnv_hmm_states 取得）、
#       该患者的上皮观测细胞(barcode)、参考细胞(barcode)
#   注：InferCNV 的 *.pred_cnv_genes.dat 是按亚群(cell_group)聚合的【长格式】，并非单细胞矩阵，
#       旧代码误将其当宽矩阵读取导致全部判 0 恶性；逐细胞判定须直接取 HMM 离散态矩阵
#       （已含单细胞分辨率的 HMM 状态），列名即细胞 barcode。
# HMM i6 状态：1=0x缺失, 2=0.5x缺失, 3=正常(1x), 4=1.5x扩增, 5=2x扩增, 6>=2.5x扩增
#   非中性 = 状态 != 3；逐细胞统计非中性基因组比例，作为 aneuploidy 分数。
#   阈值 = 文献式“相对分离”：合并 二倍体参照 + 上皮观测细胞的分数分布，用核密度找
#   双峰之间的谷(valley)——谷左=二倍体正常，谷右=非整倍体肿瘤。这是对肿瘤/正常
#   两个群体的相对分离判定，而非对参照分布的固定偏移(ref+3SD)。
#   注：infercnv 内部过滤可能丢弃少量细胞，仅对 @expr.data 中实际存在的细胞判分；
#       未进入 infercnv 的细胞由 PART 7.4 的 ifelse 默认归为 non-malignant。
# 返回：list(malignant=恶性细胞barcode, threshold=阈值, prop=各细胞非中性比例)
CallMalignantFromHMM <- function(hmm_states, obs_cells, ref_cells,
                                 abs_fallback = 0.15) {
  ed <- hmm_states
  if (!all(ed %in% 1:6)) {
    stop("hmm_states 非预期的 HMM 离散状态(1-6)，请检查 InferCNV 版本/输出")
  }

  all_cells <- colnames(ed)
  # 逐细胞非中性比例（状态 != 3 的基因占比）
  nonneut <- (ed != 3)
  n_genes <- nrow(ed)                       # 每细胞基因数（HMM 矩阵为完整基因×细胞，无 NA）
  prop    <- setNames(colSums(nonneut) / n_genes, all_cells)

  ref_present <- intersect(all_cells, ref_cells)
  obs_present <- intersect(all_cells, obs_cells)

  # 相对分离阈值：合并 二倍体参照 + EPCAM+ 细胞 的分数分布，
  # 用核密度找双峰之间的谷(valley)——谷左=二倍体正常，谷右=非整倍体肿瘤。
  # 这是对肿瘤/正常两个群体的相对分离判定（文献主流），而非对参照的固定偏移。
  if (length(ref_present) == 0L) {
    thr <- abs_fallback
    FlushLog(sprintf("    [HMM] 无参考细胞，使用绝对阈值 %.3f\n", thr))
  } else {
    combo <- c(prop[ref_present], prop[obs_present])   # 参照定义二倍体峰，上皮提供肿瘤峰
    d     <- density(combo, from = 0, to = max(combo), n = 2048)
    peak  <- which.min(abs(d$x - mean(prop[ref_present])))   # 锚定二倍体正常峰(参照均值处)
    ry    <- d$y[(peak + 1L):length(d$y)]
    if (length(ry) > 1L) {
      d1      <- diff(ry)
      locmin  <- which(d1[-length(d1)] < 0 & d1[-1] > 0) + 1L   # 峰右侧第一个局部极小 = 谷
    } else {
      locmin <- integer(0)
    }
    if (length(locmin) > 0L) {
      thr <- d$x[locmin[1L] + peak]
      FlushLog(sprintf("    [HMM] 相对分离(谷检测): 二倍体峰@%.3f, 谷阈值=%.4f (n_ref=%d, n_obs=%d)\n",
                       d$x[peak], thr, length(ref_present), length(obs_present)))
    } else {
      # 无清晰双峰结构（肿瘤/正常未分离）：退回绝对 aneuploidy 底限，明确告警
      thr <- abs_fallback
      FlushLog(sprintf("    [HMM] WARNING: 分数分布无清晰谷（肿瘤/正常未分离），退回绝对阈值 %.3f\n",
                       thr))
    }
  }
  FlushLog(sprintf("    [HMM] 上皮观测 非中性比例: median=%.3f range=[%.3f,%.3f] (n=%d)\n",
                   median(prop[obs_present]), min(prop[obs_present]),
                   max(prop[obs_present]), length(obs_present)))
  mal <- obs_present[prop[obs_present] > thr]
  list(malignant = mal, threshold = thr,
       prop = prop[c(obs_present, ref_present)])
}

# 选取本患者的 InferCNV 参照细胞：所有非上皮细胞（免疫/基质，确定二倍体），
# 按 celltype 分层下采样到 cap 以内（保各类型比例，避免大类型独占参照）。
# 每类至少保留 3 个，避免小类型被采没。
PickInferCNVRef <- function(cells, celltypes, cap) {
  if (length(cells) <= cap) return(cells)
  frac <- cap / length(cells)
  unlist(lapply(split(cells, celltypes), function(v) {
    if (length(v) <= 3L) v else sample(v, max(3L, round(length(v) * frac)))
  }), use.names = FALSE)
}

# QC 小提琴图
PlotQCViolin <- function(obj, feature, title) {
  data <- FetchData(obj, vars = c(feature, "tissue"))
  colnames(data)[1] <- "value"
  ggplot(data, aes(x = tissue, y = value, fill = tissue)) +
    geom_violin(color = "black", linewidth = 0.5,
                 scale = "width", width = 0.8) +
    geom_boxplot(width = 0.15, fill = "white", color = "black",
                 linewidth = 0.5, outlier.shape = NA, alpha = 0.8) +
    scale_fill_manual(values = tissue_colors) +
    labs(x = "", y = title, title = title) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")
}

# ================================================================
# 断点续跑管理
# ================================================================

CheckSkip <- function(label) {
  cp_file <- file.path(checkpoint_dir,
    paste("checkpoint_", label, ".rds", sep = ""))
  if (file.exists(cp_file)) {
    FlushLog(sprintf("检查点已存在 [%s]，跳过\n", label))
    return(TRUE)
  }
  return(FALSE)
}

MarkDone <- function(seurat_obj, label) {
  cp_file <- file.path(checkpoint_dir,
    paste("checkpoint_", label, ".rds", sep = ""))
  saveRDS(seurat_obj, cp_file)
  FlushLog(sprintf("  检查点已保存 [%s]\n", label))
}

# 查找首个缺失的断点：从它开始续跑。
# 关键：必须用"第一个缺失"而非"最后一个存在"决定起点——
# 否则当中间 PART(如 part035) 的 checkpoint 缺失、但其后 PART 的 checkpoint 存在时，
# 起点会被跳到末尾，中间 PART 的修复永远套不上。改为首个缺失后，part_start=缺失项索引，
# 仅用于各 PART 的 else 重加载守卫(part_start<=N)：确在此续跑时才用各自 checkpoint 覆盖，
# 避免覆盖更新的 checkpoint（运行与否由各 PART 的 CheckSkip 单独决定）。
ResumeFromLast <- function() {
  cp_order <- c("part03_qc", "part035_doublet", "part04_pca", "part05_harmony",
                "part06_annotation", "part08_subtype")
  cp_paths <- file.path(checkpoint_dir,
                        paste("checkpoint_", cp_order, ".rds", sep = ""))
  exists_vec    <- file.exists(cp_paths)
  first_missing <- which(!exists_vec)
  if (length(first_missing) == 0L) {
    # 所有 checkpoint 都在：全流程已完成，恢复最后一个并令 part_start 超过所有守卫，
    # 仅由 PART 9 收尾刷新最终成品(覆盖旧的 seurat_obj_annotated.rds)。
    FlushLog("所有检查点均存在，全流程已完成，仅刷新最终成品\n")
    obj <- readRDS(cp_paths[length(cp_paths)])
    return(list(obj = obj, part = length(cp_order) + 1L))
  }
  idx <- first_missing[1L]
  if (idx == 1L) {
    FlushLog("无任何检查点，从头开始\n")
    return(list(obj = NULL, part = 1L))
  }
  prev_cp <- cp_order[idx - 1L]
  obj     <- readRDS(cp_paths[idx - 1L])
  FlushLog(sprintf("从检查点恢复 [%s]: %d 细胞（首个待重算: %s）\n",
                   prev_cp, ncol(obj), cp_order[idx]))
  return(list(obj = obj, part = idx))
}

# ================================================================
# 运行
# ================================================================

# 注：不再以最终成品(seurat_obj_annotated.rds)存在与否决定是否跳过——
# 否则成品一旦在盘上，任何"重跑以套用修复"的尝试都会被秒退(quit)。
# 是否重跑由各 PART 的 CheckSkip 决定；末尾 PART 9 仍会刷新该最终成品。
resume <- ResumeFromLast()
seurat.obj <- resume$obj
part_start <- resume$part

# ================================================================
# PART 1-2: 数据加载 + Barcode 解析（无检查点，必须从头）
# ================================================================

if (is.null(seurat.obj)) {

  FlushLog("========== PART 1-2: 数据加载 + Barcode 解析 ==========\n")

  matrix.path   <- file.path(data_dir, "GSE178318_matrix.mtx.gz")
  barcodes.path <- file.path(data_dir, "GSE178318_barcodes.tsv.gz")
  genes.path    <- file.path(data_dir, "GSE178318_genes.tsv.gz")

  counts <- ReadMtx(
    mtx     = matrix.path,
    cells   = barcodes.path,
    features = genes.path,
    feature.column = 2  # 基因符号（匹配 MT- 线粒体基因）
  )

  FlushLog(sprintf("  原始矩阵: %d 基因 x %d 细胞\n",
                   nrow(counts), ncol(counts)))

  seurat.obj <- CreateSeuratObject(
    counts   = counts,
    project  = "GSE178318",
    min.cells    = 3,
    min.features = 1
  )
  rm(counts)
  gc()

  # Barcode 解析
  bc.split <- strsplit(colnames(seurat.obj), "_")
  seurat.obj$patient <- sapply(bc.split,
    function(x) x[length(x) - 1])
  seurat.obj$tissue  <- sapply(bc.split,
    function(x) x[length(x)])

  # CRC→PT 统一命名（barcode 原值为 CRC 表示原发灶；LM/其他保持不变）
  seurat.obj$tissue <- ifelse(seurat.obj$tissue == "CRC", "PT", seurat.obj$tissue)

  # 仅保留 treatment-naïve + PT/LM
  seurat.obj <- subset(seurat.obj,
    cells = colnames(seurat.obj)[
      seurat.obj$patient %in% naive_patients &
      seurat.obj$tissue %in% c("PT", "LM")])

  FlushLog(sprintf("  naive + PT/LM: %d 细胞\n", ncol(seurat.obj)))
  patient_tab <- table(seurat.obj$patient)
  patient_str <- paste(names(patient_tab),
    sprintf("(%d)", patient_tab), collapse = ", ")
  FlushLog(sprintf("  患者: %s\n", patient_str))
  tissue_tab <- table(seurat.obj$tissue)
  tissue_str <- paste(names(tissue_tab),
    sprintf("(%d)", tissue_tab), collapse = ", ")
  FlushLog(sprintf("  组织: %s\n", tissue_str))

  # 基因名大写
  rownames(seurat.obj) <- toupper(rownames(seurat.obj))

}  # end PART 1-2

# ================================================================
# PART 3: QC 过滤（原文四条规则）
# ================================================================

if (!CheckSkip("part03_qc")) {

  FlushLog("========== PART 3: QC 过滤 ==========\n")

  seurat.obj[["percent.mt"]]   <- PercentageFeatureSet(seurat.obj,
    pattern = "^MT-")
  seurat.obj[["log10_UMI"]]   <- log10(seurat.obj$nCount_RNA + 1)

  # 规则 (2)/(3): 每批 ±3SD —— 原文"for each batch" = 每捕获(patient×tissue)
  # 逐捕获计算 log10(UMI) 与基因数 均值±3SD（PBMC 已剔，剩 6 个肿瘤捕获）
  cap_key <- interaction(seurat.obj$patient, seurat.obj$tissue,
                         sep = "_", drop = TRUE)
  seurat.obj$outlier_3sd  <- FALSE
  seurat.obj$outlier_gene <- FALSE
  for (cap in unique(cap_key)) {
    idx     <- cap_key == cap
    mean_l  <- mean(seurat.obj$log10_UMI[idx])
    sd_l    <- sd(seurat.obj$log10_UMI[idx])
    mean_g  <- mean(seurat.obj$nFeature_RNA[idx])
    sd_g    <- sd(seurat.obj$nFeature_RNA[idx])
    seurat.obj$outlier_3sd[idx]  <- seurat.obj$log10_UMI[idx]  < (mean_l - 3 * sd_l) |
                                     seurat.obj$log10_UMI[idx]  > (mean_l + 3 * sd_l)
    seurat.obj$outlier_gene[idx] <- seurat.obj$nFeature_RNA[idx] < (mean_g - 3 * sd_g) |
                                     seurat.obj$nFeature_RNA[idx] > (mean_g + 3 * sd_g)
  }

  # QC 图（过滤前）
  p1 <- PlotQCViolin(seurat.obj, "nFeature_RNA",
                       "nFeature (before)")
  p2 <- PlotQCViolin(seurat.obj, "log10_UMI",
                       "log10(UMI) (before)")
  p3 <- PlotQCViolin(seurat.obj, "percent.mt",
                       "MT percent (before)")
  p.qc.pre <- (p1 | p2 | p3) + plot_layout(ncol = 3)
  sci_save(p.qc.pre, "QC_violin_before_filter", width = 12, height = 5)

  # 执行四条规则
  # 诊断各过滤条件的通过数
  n_passed_ngene   <- sum(seurat.obj$nFeature_RNA >= nfeature_min, na.rm = TRUE)
  n_passed_mt      <- sum(seurat.obj$percent.mt < percent_mt_max, na.rm = TRUE)
  n_passed_3sd     <- sum(!seurat.obj$outlier_3sd, na.rm = TRUE)
  n_passed_gene    <- sum(!seurat.obj$outlier_gene, na.rm = TRUE)
  n_all            <- ncol(seurat.obj)
  FlushLog(sprintf("  QC条件通过数: nFeature>=%d=%d/%d, MT<%d=%d/%d, !outlier_3sd=%d/%d, !outlier_gene=%d/%d\n",
    nfeature_min, n_passed_ngene, n_all,
    percent_mt_max, n_passed_mt, n_all,
    n_passed_3sd, n_all,
    n_passed_gene, n_all))
  
  # 执行四条规则 — 使用 Seurat 原生 subset 表达式
  seurat.obj <- subset(seurat.obj,
    subset = nFeature_RNA >= nfeature_min &
             percent.mt < percent_mt_max &
             !outlier_3sd &
             !outlier_gene)
  FlushLog(sprintf("  QC过滤后: %d 细胞\n", ncol(seurat.obj)))

  # QC 图（过滤后）
  p4 <- PlotQCViolin(seurat.obj, "nFeature_RNA",
                       "nFeature (after)")
  p5 <- PlotQCViolin(seurat.obj, "log10_UMI",
                       "log10(UMI) (after)")
  p6 <- PlotQCViolin(seurat.obj, "percent.mt",
                       "MT percent (after)")
  p.qc.post <- (p4 | p5 | p6) + plot_layout(ncol = 3)
  sci_save(p.qc.post, "QC_violin_after_filter", width = 12, height = 5)

  MarkDone(seurat.obj, "part03_qc")
  MemReport("QC完成")

} else {
  if (part_start <= 3) {
    seurat.obj <- readRDS(file.path(checkpoint_dir,
      "checkpoint_part03_qc.rds"))
    FlushLog(sprintf("  QC检查点跳过，加载: %d 细胞\n",
                      ncol(seurat.obj)))
  }
}

# ================================================================
# PART 3.5: 强制双细胞检测与过滤（DoubletFinder，必选）
# ================================================================
# 双细胞会污染 InferCNV 恶性判定：
#   - 混入恶性集 → 污染肿瘤 GRN（03/04 的 VirtualKO/VirtualOE）；
#   - 混入参考集(免疫+CAF) → 改变 InferCNV 居中尺度，整列 CNV 漂移。
# 故设为必选 QC 步骤。DoubletFinder 未安装时直接终止，不允许静默跳过。
# 续跑接入统一 CheckSkip / MarkDone 链（checkpoint_dir 下 checkpoint_part035_doublet.rds）；
# 与其它 PART 一致：checkpoint 存在即跳过，缺失才重算。未来若改了 DoubletFinder 算法/
# 参数需重算，直接删 checkpoint_part035_doublet.rds 即可（标准续跑重算方式）。
cp035 <- file.path(checkpoint_dir, "checkpoint_part035_doublet.rds")
if (!CheckSkip("part035_doublet")) {
  FlushLog("========== PART 3.5: DoubletFinder 双细胞检测与过滤（逐患者×组织=逐10x捕获）==========\n")
  library(DoubletFinder)
  # 双细胞是 10x 通道(捕获)内产物，必须按「患者×组织」(=每条 10x 捕获)分别检测，
  # 不可按患者合并 PT+LM 跑：合并后 nExp=m·inv·N² 会把 10x 多重率放大 ~2×，且跨组织
  # 杂交体虚高 pANN（把另一组织的真实 singlet 误判为 doublet）。逐捕获后距离矩阵仅
  # n_cap×nExp_cap（最大 ~1.4 万×1 千 ≈ 0.1GB），内存可控且方法学正确。
  # Windows 不支持 mclapply 多核（mc.cores>1 报错），故 Windows 仅单核。
  # DoubletFinder 2.0.6 已合并版本命名：paramSweep()/doubletFinder()（小写 d）/
  # find.pK() 沿用不变，无 reuse.pca 参数。
  ncores_df    <- if (.Platform$OS.type == "windows") 1L else ncores_nonwin
  # 捕获单元 = 患者×组织（每条 10x 通道独立检测），循环见下
  dbl_call     <- structure(rep(NA_character_, ncol(seurat.obj)),
                            names = colnames(seurat.obj))
  pk_list      <- list()   # 各患者 pK 扫描结果（出图用）
  umap_list    <- list()   # 各患者双细胞 UMAP（出图用）

  captures <- unique(data.frame(patient = seurat.obj$patient,
                                tissue  = seurat.obj$tissue))
  FlushLog(sprintf("  捕获单元数（患者×组织）: %d\n", nrow(captures)))
  for (i in seq_len(nrow(captures))) {
    pat <- as.character(captures$patient[i])
    tis <- as.character(captures$tissue[i])
    cap <- paste0(pat, "_", tis)            # 捕获单元标签（=单条 10x 通道）
    idx <- seurat.obj$patient == pat & seurat.obj$tissue == tis
    n_cap <- sum(idx)
    if (n_cap < min_cells_doublet) {   # 捕获过小，DoubletFinder 不可靠，默认 Singlet
      dbl_call[idx] <- "Singlet"
      FlushLog(sprintf("  [%s] 细胞数 %d 过小，跳过 DoubletFinder，默认 Singlet\n", cap, n_cap))
      next
    }
    FlushLog(sprintf("  [%s] 双细胞检测（%d 细胞）...\n", cap, n_cap))
    sub <- subset(seurat.obj, subset = patient == pat & tissue == tis)
    sub <- NormalizeData(sub, verbose = FALSE)
    sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = df_nfeatures, verbose = FALSE)
    sub <- ScaleData(sub, verbose = FALSE)
    sub <- RunPCA(sub, npcs = npcs_use, verbose = FALSE)
    # 自动扫描择优 pK（逐捕获；fail-fast：pK 扫描失败即报错，不再静默回退）
    sweep.res   <- paramSweep(sub, PCs = 1:npcs_use, sct = FALSE, num.cores = ncores_df)
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn       <- find.pK(sweep.stats)
    rm(sweep.res, sweep.stats); gc()
    pk_opt <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    bcmvn$capture <- cap; pk_list[[cap]] <- bcmvn
    # 期望双细胞数：10x 官方多重率（逐捕获恢复细胞数自适应）。
    # 双细胞是 10x 通道内产物，必须逐捕获算 nExp 再求和；合并患者总数算会因
    # nExp=m·inv·N² 把多重率放大 ~2×，且跨组织杂交体虚高 pANN（误判 singlet）。
    nExp_exp <- ExpectedDoubletNExp(n_cap, chemistry = df_chemistry)
    FlushLog(sprintf("  [%s] 期望双细胞数 nExp=%.0f (10x %s, 期望率 %.2f%%)\n",
                     cap, nExp_exp, df_chemistry, 100 * nExp_exp / n_cap))
    sub <- doubletFinder(sub, PCs = 1:npcs_use, pN = df_pN, pK = pk_opt,
                         nExp = nExp_exp, sct = FALSE)
    dcol <- colnames(sub@meta.data)[grepl("^DF.classifications", colnames(sub@meta.data))]
    dcol <- dcol[length(dcol)]
    pcol <- colnames(sub@meta.data)[grepl("^pANN", colnames(sub@meta.data))]
    pcol <- pcol[length(pcol)]
    if (length(pcol) == 0L || length(dcol) == 0L) {
      stop(sprintf("[%s] DoubletFinder 未生成 pANN/DF.classifications 列（预期 pANN_* / DF.classifications_*），请检查 DoubletFinder 版本",
                   cap))
    }
    pann <- sub@meta.data[[pcol]]
    # 分类由 DoubletFinder 原生率驱动（nExp = ExpectedDoubletNExp 期望双细胞数），
    # pANN 降序前 nExp 个判为 Doublet，与 pANN 分布形态无关。诊断仅保留期望率参考线。
    cls <- sub@meta.data[[dcol]]
    FlushLog(sprintf("  [%s] DF原生分类 nExp=%.0f → Doublet %d/%d (%.1f%%)\n",
                     cap, nExp_exp, sum(cls == "Doublet"), n_cap, 100 * sum(cls == "Doublet") / n_cap))

    # pANN 分布诊断图：期望率 cutoff 线（蓝点线，DF 原生分类位置）
    cut_exp <- sort(pann, decreasing = TRUE)[min(nExp_exp, length(pann))]
    pann_diag <- ggplot(data.frame(pANN = pann), aes(x = pANN)) +
      geom_histogram(aes(y = after_stat(density)), bins = 100,
                     fill = "steelblue", color = NA, alpha = 0.6) +
      geom_density(color = "darkred", linewidth = 0.6) +
      geom_vline(xintercept = cut_exp,  linetype = "dotted", color = "blue", linewidth = 0.8) +
      annotate("text", x = cut_exp - 0.02, y = Inf,
               label = sprintf("expected nExp=%d: %.1f%%\n[DF-native]",
                               nExp_exp, 100*mean(cls=="Doublet")),
               hjust = 1, vjust = 1.5, size = 3, color = "blue") +
      sci_theme + labs(title = sprintf("%s pANN distribution (N=%d)", cap, n_cap),
                        x = "pANN", y = "Density") +
      theme(aspect.ratio = 0.65)
    sci_save(pann_diag, sprintf("03.5_pANN_%s", cap), width = 8, height = 5)
    # 保存 pANN 诊断数据，便于以后单独重绘/调样式而无需重跑 DoubletFinder
    saveRDS(list(pann = pann, cls = cls, nExp_exp = nExp_exp,
                 cut_exp = cut_exp, n_cap = n_cap),
            file.path(tab_dir, sprintf("pann_data_%s.rds", cap)))
    rm(pann_diag); gc()
    sub  <- RunUMAP(sub, dims = 1:npcs_use, verbose = FALSE)
    p_umap <- DimPlot(sub, group.by = dcol, reduction = "umap",
                      label = FALSE, pt.size = 0.4) + sci_theme +
      theme(aspect.ratio = 1,
            axis.title = element_blank(),
            axis.text  = element_blank(),
            axis.ticks = element_blank()) +
      labs(title = sprintf("%s (Doublet %.1f%%)", cap,
                           100 * sum(cls == "Doublet", na.rm = TRUE) / n_cap))
    res_cap <- list(cls = cls, pK = pk_opt, umap = p_umap,
             n_d = sum(cls == "Doublet", na.rm = TRUE))
    dbl_call[idx] <- res_cap$cls
    umap_list[[cap]] <- res_cap$umap   # p_umap 在上方已生成，必为非空 ggplot
    FlushLog(sprintf("  [%s] pK=%s | Doublet=%d/%d (%.1f%%)\n",
      cap, if (is.na(res_cap$pK)) "-" else sprintf("%.3f", res_cap$pK),
      res_cap$n_d, n_cap, 100 * res_cap$n_d / n_cap))
    gc()
  }

  seurat.obj$doublet <- dbl_call
  n_dbl <- sum(seurat.obj$doublet == "Doublet", na.rm = TRUE)
  FlushLog(sprintf("  汇总: Doublet=%d / %d (%.1f%%)\n",
                   n_dbl, ncol(seurat.obj), 100 * n_dbl / ncol(seurat.obj)))

  # 可视化 1：各捕获 pK 扫描曲线（BCmetric 随 pK 变化，峰值即所选 pK）
  if (length(pk_list) > 0L) {
    pk_all <- do.call(rbind, pk_list)
    pk_plot <- ggplot(pk_all, aes(x = as.numeric(as.character(pK)), y = BCmetric,
                                  colour = capture)) +
      geom_point(size = 1.5) + geom_line() + sci_theme +
      labs(title = "DoubletFinder pK sweep (BCmetric, by capture)",
           x = "pK", y = "BCmetric")
    sci_save(pk_plot, "03.5_pK_sweep")
    FlushLog("  已保存 pK 扫描图: 03.5_pK_sweep.pdf\n")
  }
  # 可视化 2：各捕获双细胞分类 UMAP（过滤前，肉眼确认 doublet 分布合理）
  # 优化：去掉 panel 内冗余文字(label/轴标题/轴刻度)，6 图共享一个 legend，避免重复占位
  if (length(umap_list) > 0L) {
    dbl_umap <- patchwork::wrap_plots(umap_list, ncol = 3, guides = "collect")
    sci_save(dbl_umap, "03.5_doublet_UMAP", width = 12, height = 6)
    FlushLog("  已保存双细胞分类 UMAP（逐捕获）: 03.5_doublet_UMAP.pdf\n")
    # 同时保存 umap_list，便于以后单独重绘/调样式而无需重跑 DoubletFinder
    saveRDS(umap_list, file.path(tab_dir, "umap_list_doublet.rds"))
  }

  # 强制过滤：移除双细胞（仅保留 Singlet）
  seurat.obj <- subset(seurat.obj, subset = doublet == "Singlet")
  FlushLog(sprintf("  已移除 %d 个双细胞，剩余 %d 细胞\n", n_dbl, ncol(seurat.obj)))
  MemReport("DoubletFinder完成")
  MarkDone(seurat.obj, "part035_doublet")   # 统一续跑链（替代旧 saveRDS(dbl_rds)+writeLines(dbl_flag)）
} else if (part_start <= 4) {
  # 续跑：checkpoint 命中且本就在此续跑时恢复已过滤对象（避免覆盖更新的 checkpoint）
  seurat.obj <- readRDS(cp035)
  FlushLog("  [DoubletFinder] 已完成，恢复已过滤对象\n")
}

# ================================================================
# PART 4: 标准化 + HVG + 回归 + PCA
# ================================================================

if (!CheckSkip("part04_pca")) {

  FlushLog("========== PART 4: 标准化 / HVG / 回归 / PCA ==========\n")

  seurat.obj <- NormalizeData(seurat.obj,
    normalization.method = "LogNormalize", scale.factor = 10000)

  seurat.obj <- FindVariableFeatures(seurat.obj,
      selection.method = "dispersion",
      mean.cutoff = c(0.0125, 3),
      dispersion.cutoff = 0.5,
      nfeatures = n_hvg)          # 原文：mean∈[0.0125,3] 且 dispersion>0.5 → 1511 个 HVG（P2-⑧）

  seurat.obj <- ScaleData(seurat.obj,
    vars.to.regress = c("nCount_RNA", "percent.mt"))

  seurat.obj <- RunPCA(seurat.obj, npcs = npcs_total)

  # 图
  p.elbow <- ElbowPlot(seurat.obj, ndims = npcs_total) + sci_theme
  sci_save(p.elbow, "ElbowPlot", width = 6, height = 4)

  top10.hvg <- head(VariableFeatures(seurat.obj), 10)
  p.hvg <- LabelPoints(VariableFeaturePlot(seurat.obj),
                        points = top10.hvg, repel = TRUE) + sci_theme
  sci_save(p.hvg, "VariableFeatures", width = 8, height = 6)

  FlushLog(sprintf("  HVG: %d\n", length(VariableFeatures(seurat.obj))))

  MarkDone(seurat.obj, "part04_pca")
  MemReport("PCA完成")

} else {
  if (part_start <= 4) {
    seurat.obj <- readRDS(file.path(checkpoint_dir,
      "checkpoint_part04_pca.rds"))
    FlushLog(sprintf("  PCA检查点跳过，加载: %d 细胞\n",
                      ncol(seurat.obj)))
  }
}

# ================================================================
# PART 5: Harmony 批次校正 + UMAP + 聚类
# ================================================================

if (!CheckSkip("part05_harmony")) {

  FlushLog("========== PART 5: Harmony 批次校正 / UMAP / 聚类 ==========\n")

  seurat.obj <- RunHarmony(seurat.obj, "patient",
                           reduction.save = "harmony")
  seurat.obj <- RunUMAP(seurat.obj, reduction = "harmony",
                         dims = 1:npcs_use)
  seurat.obj <- FindNeighbors(seurat.obj, reduction = "harmony",
                              dims = 1:npcs_use)
  seurat.obj <- FindClusters(seurat.obj, resolution = cluster_res)

  ncluster <- length(levels(seurat.obj@meta.data$seurat_clusters))
  FlushLog(sprintf("  %d clusters (res=%.1f)\n",
                   ncluster, cluster_res))

  # UMAP 图
  p.umap.cluster <- DimPlot(seurat.obj, reduction = "umap",
    group.by = "seurat_clusters", label = TRUE) +
    theme(aspect.ratio = 1) + sci_theme
  sci_save(p.umap.cluster, "UMAP_clusters", width = 8, height = 6)

  p.umap.tissue <- DimPlot(seurat.obj, reduction = "umap",
    group.by = "tissue", cols = tissue_colors) +
    theme(aspect.ratio = 1) + sci_theme
  sci_save(p.umap.tissue, "UMAP_tissue", width = 8, height = 6)

  p.umap.patient <- DimPlot(seurat.obj, reduction = "umap",
    group.by = "patient") +
    theme(aspect.ratio = 1) + sci_theme
  sci_save(p.umap.patient, "UMAP_patient", width = 8, height = 6)

  MarkDone(seurat.obj, "part05_harmony")
  MemReport("聚类完成")

} else {
  if (part_start <= 5) {
    seurat.obj <- readRDS(file.path(checkpoint_dir,
      "checkpoint_part05_harmony.rds"))
    FlushLog(sprintf("  聚类检查点跳过，加载: %d 细胞\n",
                      ncol(seurat.obj)))
  }
}

# ================================================================
# PART 6: 细胞类型注释 — DotPlot 双参考验证
# ================================================================

# ---- 降维层可用性保证（resume / 各 PART 跳过分支的【单一汇聚点】）----
# 见顶部 helper ensure_reductions()：保证 pca/harmony/umap 三者齐备——下游
# PART 6/7.5/7b 依赖 umap，PART 8 依赖 harmony。从磁盘恢复的 checkpoint 若因
# 旧版脚本残留而缺任一降维层，此处就地从最深的缺失层逐级补算（pca→harmony→umap），
# 参数与 PART 4/5 完全一致，不重跑聚类/注释、不引入偏差。完整对象此处为 no-op。
seurat.obj <- ensure_reductions(seurat.obj)

if (!CheckSkip("part06_annotation")) {

  FlushLog("========== PART 6: DotPlot 双参考注释 ==========\n")

  # 参考1: 原文标记基因
  paper_markers <- list(
    "Epithelial/Tumor" = c("EPCAM", "KRT8", "KRT18"),
    "T cells"          = c("CD3D", "CD3E", "TRAC"),
    "B cells"          = c("CD19", "CD79A", "MS4A1"),
    "Plasma cells"     = c("IGHG1", "IGHA1", "MZB1"),
    "Mono/Macrophage"  = c("CD68", "CD163", "CD14", "LYZ"),
    "NK cells"         = c("KLRF1", "KLRD1", "FGFBP2", "PRF1"),
    "CAFs"             = c("FAP", "ACTA2", "TAGLN"),
    "Endothelial"      = c("CLDN5", "CDH5", "VWF"),
    "pDC"              = c("IL3RA", "CLEC4C"),
    "Mast cells"       = c("TPSAB1", "TPSB2", "MS4A2"),
    "DC"               = c("CD1C", "CLEC9A", "LILRA4", "LAMP3")
  )

  # 参考2: HCA/Azimuth
  hca_markers <- list(
    "Epithelial"   = c("EPCAM", "KRT5", "KRT14"),
    "T cells"      = c("CD3D", "CD3G", "IL7R"),
    "B cells"      = c("CD79A", "MS4A1", "TCL1A"),
    "Plasma cells" = c("MZB1", "IGKC", "JCHAIN"),
    "Mono/Mac"     = c("CD14", "FCGR3A", "CD68"),
    "NK cells"     = c("KLRF1", "KLRF2", "NKG7"),
    "CAFs"         = c("COL1A1", "DCN", "FAP", "THY1"),
    "Endothelial"  = c("CLDN5", "PECAM1", "VWF"),
    "pDC"          = c("IL3RA", "CLEC4C"),
    "Mast cells"   = c("TPSB2", "TPSAB1", "HDC"),
    "DC"           = c("HLA-DPA1", "HLA-DPB1", "CD1C")
  )

  # 仅保留存在基因
  paper_markers <- lapply(paper_markers,
    function(g) g[g %in% rownames(seurat.obj)])
  hca_markers   <- lapply(hca_markers,
    function(g) g[g %in% rownames(seurat.obj)])

  # DotPlot 原文参考
  p.dot1 <- DotPlot(seurat.obj,
    features = paper_markers, group.by = "seurat_clusters",
    cols = c("lightgrey", "blue")) +
    RotatedAxis() +
    scale_y_discrete(limits = function(x) rev(x)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 9),
      panel.grid.major = element_line(color = "grey95")
    ) +
    labs(title = "DotPlot — Ref1: Paper Markers",
         x = "", y = "Cluster")
  sci_save(p.dot1, "DotPlot_PaperMarkers", width = 20, height = 8)

  # DotPlot HCA 参考
  p.dot2 <- DotPlot(seurat.obj,
    features = hca_markers, group.by = "seurat_clusters",
    cols = c("lightgrey", "red")) +
    RotatedAxis() +
    scale_y_discrete(limits = function(x) rev(x)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 9),
      panel.grid.major = element_line(color = "grey95")
    ) +
    labs(title = "DotPlot — Ref2: HCA/Azimuth Markers",
         x = "", y = "Cluster")
  sci_save(p.dot2, "DotPlot_HCA_Markers", width = 20, height = 8)

  # FeaturePlot 核心标记
  core.markers <- c("EPCAM", "CD3D", "CD19", "MZB1", "CD68",
                     "KLRF1", "FAP", "CLDN5", "TPSB2")
  core.present <- core.markers[core.markers %in% rownames(seurat.obj)]
  p.feat.core <- FeaturePlot(seurat.obj, features = core.present,
    ncol = 3, min.cutoff = "q10", max.cutoff = "q90") +
    scale_color_viridis_c(option = "plasma") +
    theme(aspect.ratio = 1) + sci_theme
  sci_save(p.feat.core, "FeaturePlot_CoreMarkers",
            width = 12, height = 12)

  # 自动注释（基于两参考各类型的平均表达）
  AssignCelltype <- function(marker_list, obj) {
    # 获取平均表达量（cluster x genes）
    mat <- AverageExpression(obj,
      features = unlist(marker_list),
      group.by = "seurat_clusters",
      layer = "data")[[1]]
    if (nrow(mat) == 0 || ncol(mat) == 0) {
      clusters <- sort(unique(as.character(obj$seurat_clusters)))
      return(setNames(rep("Unknown", length(clusters)), clusters))
    }
    # 对每个细胞类型，计算标记基因平均表达量（按 cluster）
    score_list <- lapply(marker_list, function(g) {
      g2 <- g[g %in% rownames(mat)]
      if (length(g2) == 0) {
        rep(0, ncol(mat))
      } else {
        colMeans(mat[g2, , drop = FALSE])
      }
    })
    # 转为矩阵（cluster x cell_type）
    score_mat <- do.call(cbind, score_list)
    colnames(score_mat) <- names(marker_list)
    rownames(score_mat) <- colnames(mat)
    # 每个 cluster 分配得分最高的细胞类型
    apply(score_mat, 1, function(x) {
      if (all(is.na(x))) return(NA_character_)
      names(which.max(x))
    })
  }

  cluster_labels_paper <- AssignCelltype(paper_markers, seurat.obj)
  cluster_labels_hca   <- AssignCelltype(hca_markers,   seurat.obj)

  # 写回每个 cell (使用 AddMetaData 规避 Seurat v5 [[<- 问题)
  seurat.obj <- AddMetaData(seurat.obj,
    unname(cluster_labels_paper[paste0("g", 
      as.character(seurat.obj$seurat_clusters))]),
    col.name = "celltype_ref1")
  seurat.obj <- AddMetaData(seurat.obj,
    unname(cluster_labels_hca[paste0("g",
      as.character(seurat.obj$seurat_clusters))]),
    col.name = "celltype_ref2")
  FlushLog(sprintf("  celltype_ref1 NAs=%d, ref2 NAs=%d\n",
    sum(is.na(seurat.obj$celltype_ref1)),
    sum(is.na(seurat.obj$celltype_ref2))))

  # 注释表（人工填写用）：先读入已有 final_annotation / notes，避免重跑覆盖人工校正
  tbl3_path <- file.path(tab_dir, "Table3_ClusterAnnotation_template.csv")
  prev_final <- NULL
  prev_notes <- NULL
  if (file.exists(tbl3_path)) {
    prev <- read.csv(tbl3_path, stringsAsFactors = FALSE)
    if (!is.null(prev$final_annotation))
      prev_final <- setNames(as.character(prev$final_annotation),
                             sub("^g", "", as.character(prev$cluster)))
    if (!is.null(prev$notes))
      prev_notes <- setNames(as.character(prev$notes),
                             sub("^g", "", as.character(prev$cluster)))
  }
  cl_names <- sub("^g", "", names(cluster_labels_paper))
  annotation_csv <- data.frame(
    cluster       = cl_names,
    ref1_paper    = cluster_labels_paper,
    ref2_hca      = cluster_labels_hca,
    ref1_eq_ref2  = cluster_labels_paper == cluster_labels_hca,
    final_annotation = if (!is.null(prev_final)) {
                          v <- prev_final[cl_names]
                          ifelse(is.na(v), "", v)
                        } else "",
    notes = if (!is.null(prev_notes)) {
              v <- prev_notes[cl_names]
              ifelse(is.na(v), "", v)
            } else "",
    stringsAsFactors = FALSE
  )
  write.csv(annotation_csv, file = tbl3_path, row.names = FALSE)
  FlushLog("  Table3_ClusterAnnotation_template.csv 已保存（保留既有 final_annotation/notes）\n")

  # 最终 celltype：优先采用人工校正 final_annotation（贴合原文 "canonical marker ... confirmed"），
  # 未填写则回退自动标注 celltype_ref1
  final_ct <- seurat.obj@meta.data[["celltype_ref1"]]
  # tbl3_path 已于上方 write.csv 写出，此处必然存在，直接读取应用人工校正
  ann_tbl <- read.csv(tbl3_path, stringsAsFactors = FALSE)
  filled <- ann_tbl[!is.na(ann_tbl$final_annotation) &
                    trimws(as.character(ann_tbl$final_annotation)) != "", ,
                  drop = FALSE]
  if (nrow(filled) > 0L) {
    lut <- setNames(as.character(filled$final_annotation),
                    sub("^g", "", as.character(filled$cluster)))
    cl <- as.character(seurat.obj$seurat_clusters)
    mapped <- lut[cl]
    applied <- !is.na(mapped)
    final_ct[applied] <- mapped[applied]
    FlushLog(sprintf("  应用人工校正 final_annotation：%d/%d 细胞\n",
                     sum(applied), length(cl)))
  } else {
    FlushLog("  Table3 final_annotation 未填写，沿用自动 celltype_ref1\n")
  }
  seurat.obj <- AddMetaData(seurat.obj, metadata = final_ct, col.name = "celltype")
  FlushLog(sprintf("  celltype unique: %s (NAs=%d)\n",
    paste(unique(seurat.obj$celltype), collapse=", "),
    sum(is.na(seurat.obj$celltype))))

  # 标注不一致
  inconsistent <- which(cluster_labels_paper != cluster_labels_hca)
  inconsistent_clusters <- names(cluster_labels_paper)[inconsistent]
  seurat.obj <- AddMetaData(seurat.obj,
    ifelse(seurat.obj$celltype_ref1 != seurat.obj$celltype_ref2,
      "CONFLICT", "OK"),
    col.name = "annotation_conflict")
  if (length(inconsistent_clusters) > 0) {
    conflict_str <- paste(inconsistent_clusters, collapse = ", ")
  } else {
    conflict_str <- "无"
  }
  FlushLog(sprintf("  标注冲突 clusters: %s\n", conflict_str))

  # UMAP 细胞类型（ggrepel 避让标签）
  lv <- levels(factor(seurat.obj$celltype))
  p.cols <- setNames(scales::hue_pal()(length(lv)), lv)
  umap_labeled_plot(seurat.obj, "celltype", p.cols, "UMAP_celltype",
                    width = 10, height = 7)

  p.umap.ct.split <- DimPlot(seurat.obj, reduction = "umap",
    group.by = "celltype", split.by = "tissue",
    cols = p.cols, ncol = 2) +
    theme(aspect.ratio = 1) + sci_theme
  sci_save(p.umap.ct.split, "UMAP_celltype_by_tissue",
            width = 14, height = 6)

  # ================================================================
  # 补充：SingleR 独立注释（与原文注释互相对照）
  # 参考数据集优先用本地预缓存 ref_HPCA_cached.rds（Human Primary Cell Atlas，
  # 已落为纯内存矩阵 SummarizedExperiment，避免 celldex 在线拉的悬空依赖问题）；
  # 缓存缺失时才回退到 celldex::HumanPrimaryCellAtlasData() 在线获取。
  # 对 scrna-seq 数据进行单细胞水平的无偏标签转移。
  # 输出：UMAP_SingleR_Labeled.pdf + celltype_singler 列
  # ================================================================
  # SingleR 为必出产物（celltype_singler 列 + UMAP_SingleR_Labeled.pdf +
  # Table_S2_Paper_vs_SingleR.csv），不允许静默失败——故移除 tryCatch 与
  # "未安装则跳过"分支：任何运行期错误直接抛出终止；依赖缺失同样报错
  # （脚本顶部 library(SingleR) 已硬依赖，celldex 缺失时 celldex:: 调用会
  # 抛出清晰错误）。把 NA 列伪装成"跑过"最危险，必须 fail-fast。
  FlushLog("  [SingleR] 开始独立注释（必出产物，失败即终止）...\n")
  # 参考数据集：优先读本地预缓存 ref_HPCA_cached.rds。
  # 原因：celldex::HumanPrimaryCellAtlasData() 在线拉取在部分环境下会因
  # 悬空依赖 'scrapper'（该包不在 CRAN）导致 loadNamespace 失败；预缓存为
  # 纯内存矩阵 SummarizedExperiment，readRDS 即可，不再触发 celldex 加载。
  hpca_cache <- file.path(tab_dir, "ref_HPCA_cached.rds")
  if (file.exists(hpca_cache)) {
    FlushLog(sprintf("  [SingleR] 载入本地缓存参考: %s\n", basename(hpca_cache)))
    ref <- readRDS(hpca_cache)
  } else {
    FlushLog("  [SingleR] 本地无 HPCA 缓存，回退到 celldex:: 在线获取...\n")
    ref <- celldex::HumanPrimaryCellAtlasData()
  }
  singler.res <- SingleR::SingleR(
    test = GetAssayData(seurat.obj, layer = "data"),
    ref = ref, labels = ref$label.main,
    de.method = "wilcox")
  # 用 AddMetaData 替代 $<-，避免 Seurat v5 [[<- 内部 names 检查失败
  seurat.obj <- AddMetaData(seurat.obj,
    metadata = as.character(singler.res$pruned.labels),
    col.name = "celltype_singler")
  nct_s <- length(unique(na.omit(singler.res$pruned.labels)))
  # SingleR UMAP：手动计算 cluster 中心 + ggrepel 避让标签
  umap_labeled_plot(seurat.obj, "celltype_singler",
                    scales::hue_pal()(nct_s), "UMAP_SingleR_Labeled",
                    width = 12, height = 7)
  ct_compare <- table(paper = seurat.obj$celltype,
    singler = seurat.obj$celltype_singler)
  write.csv(as.data.frame.matrix(ct_compare),
    file = file.path(tab_dir,
      "Table_S2_Paper_vs_SingleR.csv"))
  FlushLog(sprintf(
    "  [SingleR] 完成：%d 种类型，UMAP + 对照表已保存\n",
    nct_s))

  MarkDone(seurat.obj, "part06_annotation")   # 规范 checkpoint：PART 6 完成态的唯一真相源
  # 清理已废弃的旁路副本 seurat_obj_celltyped.rds（下游 02/03/04/06 只读取 PART 9 的
  # seurat_obj_annotated.rds，该副本已无任何读取方；双对象并存只会制造 staleness 隐患）
  celltyped_legacy <- file.path(tab_dir, "seurat_obj_celltyped.rds")
  if (file.exists(celltyped_legacy)) file.remove(celltyped_legacy)
  MemReport("注释完成")

} else {
  if (part_start <= 6) {
    # 跳过分支加载【规范 checkpoint】而非旁路副本：消除双对象 staleness 隐患
    seurat.obj <- readRDS(file.path(tab_dir, "checkpoint_part06_annotation.rds"))
    FlushLog(sprintf("  PART 6 跳过，加载 checkpoint_part06_annotation.rds: %d 细胞\n",
                      ncol(seurat.obj)))
    FlushLog(sprintf("  celltype 列: %s\n",
                      paste(unique(seurat.obj$celltype), collapse = ", ")))
  }
}

# ================================================================
# PART 7: InferCNV — 恶性细胞鉴定（修复版）
# ================================================================
# 依据：Che et al. 2021 原文 Page 18（InferCNV on EPCAM+ epithelial cells,
#   denoise=T, HMM=T, cutoff=0.1, 窗口101基因；参照为两个非恶性上皮簇）
#
# analysis_mode 必须为 "subclusters"（2026-08-04 源码实证的关键根因）：
#   曾用 "samples" 复现原文，导致三样本 prop_mal 全 = 1.000。真因见
#   infercnv:::predict_CNV_via_HMM_on_whole_tumor_samples —— 该函数对整个观测组
#   取 rowMeans 得一条均值谱、只跑一次 Viterbi，再把同一状态向量广播回组内每个细胞
#   （源码 L23/L30/L31），故组内所有细胞 CNV 谱逐基因完全相同、aneuploidy sd=0，
#   逐细胞谷检测在数学上不可能分离。且 .get_state_emission_params 的状态 sd 按
#   num_cells 拟合（∝1/√n），n=5050 时 sd 极窄，令均值谱上微弱偏移即判非中性。
#   "subclusters"（包默认值）改为在每个 leiden 子簇上各跑一次 HMM，上皮内部得到
#   多条 CNV 谱，恶性/正常方可分离。COL12 实测：唯一 aneuploidy 值 1→6，
#   分布由单点 0.5502 变为双峰（74 个 <0.05 正常 + 153 个 >0.40 恶性，中间空白）。
#
# 参照 = 本患者非上皮细胞（免疫/基质，确定二倍体），分层下采样至 infercnv_ref_cap。
#   注：曾改用"上皮表达共识合成参照"，实测无效（COL07 aneuploidy 仅 0.3058→0.3063）
#   且有害——hspike（HMM 状态分布基准）直接由参照细胞构建
#   （infercnv:::.build_and_add_hspike L9），全同副本的合成细胞会令状态分布失真。
#   注释文件：col.names=FALSE（官方格式，无表头）
# 按患者拆分运行，控制内存

if (!CheckSkip("part07_infercnv")) {

  FlushLog("========== PART 7: InferCNV 恶性鉴定（修复版） ==========\n")

  # ----------------------------------------------------------------
  # 7.1 上皮观测集（仅上皮区室；杜绝免疫细胞经 EPCAM>1 泄漏进观测 — 根因已实证）
  # ----------------------------------------------------------------
  # 直接用 metadata 列选型（与 L1460/L1563 同惯用法）；
  # 不用 WhichCells(expression=) —— 其第 4 形参 expression 易被位置匹配误绑到 cells 导致报错
  epi_obs <- colnames(seurat.obj)[seurat.obj$celltype %in% epi_types]
  FlushLog(sprintf("  上皮观测细胞: %d\n", length(epi_obs)))

  # ----------------------------------------------------------------
  # 7.2 参照池：非上皮细胞（免疫/基质，确定二倍体），逐患者取用
  # ----------------------------------------------------------------
  ref_pool_all <- colnames(seurat.obj)[!(seurat.obj$celltype %in% epi_types)]
  FlushLog(sprintf("  非上皮参照池(全体): %d\n", length(ref_pool_all)))

  # gene_order.txt：从 GENCODE GTF 构建真实染色体坐标（见 BuildGeneOrderFromGTF）
  gene_order_path <- file.path(tab_dir, "gene_order.txt")
  BuildGeneOrderFromGTF(gencode_gtf, gene_order_path)

  # ----------------------------------------------------------------
  # 7.3 分患者跑 InferCNV
  # ----------------------------------------------------------------
  malignant_cells <- c()
  patients <- naive_patients

  for (p in patients) {
    FlushLog(sprintf("\n  --- 患者 %s ---\n", p))

    epi_p <- epi_obs[grepl(p, epi_obs)]
    if (length(epi_p) < min_epcam) {
      FlushLog(sprintf("    跳过: 上皮观测细胞不足 (%d)\n", length(epi_p)))
      next
    }
    FlushLog(sprintf("    上皮观测 (observation): %d\n", length(epi_p)))

    # 参照：本患者非上皮细胞（免疫/基质，确定二倍体），按 celltype 分层下采样
    ref_pool_p <- ref_pool_all[grepl(p, ref_pool_all)]
    set.seed(seed_ref_base + which(patients == p))
    ref_p <- PickInferCNVRef(ref_pool_p,
                             seurat.obj$celltype[match(ref_pool_p, colnames(seurat.obj))],
                             infercnv_ref_cap)
    FlushLog(sprintf("    非上皮参照 (reference): %d / 池 %d\n",
                     length(ref_p), length(ref_pool_p)))

    counts_p <- GetAssayData(seurat.obj, layer = "counts")[, c(epi_p, ref_p), drop = FALSE]

    # 注释文件（无表头，符合官方 wiki 格式）
    anno_df <- data.frame(
      cell_id   = colnames(counts_p),
      cell_type = ifelse(colnames(counts_p) %in% ref_p, "ref", "epithelial"),
      stringsAsFactors = FALSE
    )
    anno_file <- file.path(tab_dir,
      paste("infercnv_anno_", gsub(" ", "_", p), ".txt", sep = ""))
    write.table(anno_df, anno_file,
                sep = "\t", quote = FALSE, row.names = FALSE,
                col.names = FALSE)   # ← 官方格式无表头
    FlushLog(sprintf("    注释文件: %s\n", basename(anno_file)))

    out_dir_p <- file.path(tab_dir,
      paste("infercnv_out_", gsub(" ", "_", p), sep = ""))
    dir.create(out_dir_p, recursive = TRUE, showWarnings = FALSE)

    # 不再内嵌 run_done 闸：续跑统一由 CheckSkip("part07_infercnv") 控制。
    # 进入此分支即 PART 7 主动重算，故每次均清空旧输出并重跑 InferCNV，
    # 避免「删主 checkpoint 却漏删 infercnv_out_* 时复用陈旧结果」（P2-⑯）。
    unlink(out_dir_p, recursive = TRUE)
    dir.create(out_dir_p, recursive = TRUE, showWarnings = FALSE)

    infercnv_obj <- CreateInfercnvObject(
      raw_counts_matrix = counts_p,
      annotations_file = anno_file,
      gene_order_file   = gene_order_path,
      ref_group_names   = "ref"
    )

    FlushLog(sprintf("    开始 InferCNV (%d 上皮 + %d 非上皮参照)...\n",
                     length(epi_p), length(ref_p)))
    # up_to_step=17：本机 infercnv 补丁版 HMM 落盘在 STEP17（非 1.23.0 文档的18），
    #   跑完 HMM 即停，跳过 BayesNet/Filter（恶性判定只用原始 HMM 六态，见 get_infercnv_hmm_states）。
    #   step_count 的递增点都在 analysis_mode 判断之外，故 samples/subclusters 下 HMM 同为第 17 步。
    infercnv_obj <- infercnv::run(infercnv_obj,
      num_threads        = infercnv_threads,
      cutoff            = 0.1,
      out_dir           = out_dir_p,
      cluster_by_groups = TRUE,
      denoise           = TRUE,
      HMM               = TRUE,
      analysis_mode     = "subclusters",   # 包默认值；"samples" 会把组均值谱广播给每个细胞致判定退化（见 PART7 抬头注释）
      window_length     = 101,
      up_to_step        = 17,
      save_rds          = TRUE,
      no_plot           = TRUE
    )

    # 基于 HMM i6 六态逐细胞判定恶性（相对分离/谷检测，见 CallMalignantFromHMM）
    hmm_res <- CallMalignantFromHMM(get_infercnv_hmm_states(out_dir_p), epi_p, ref_p)
    mal_p <- hmm_res$malignant
    malignant_cells <- c(malignant_cells, mal_p)
    FlushLog(sprintf("    恶性 (HMM 非中性>%.4f): %d / %d\n",
      hmm_res$threshold, length(mal_p), length(epi_p)))
  }

  # ----------------------------------------------------------------
  # 7.4 汇总，写入 seurat metadata
  # ----------------------------------------------------------------
  n_mal <- length(malignant_cells)
  FlushLog(sprintf("\n  总恶性细胞: %d / %d 上皮观测\n", n_mal, length(epi_obs)))

  is_mal <- is.element(colnames(seurat.obj), malignant_cells)
  seurat.obj[["malignant"]]     <- ifelse(is_mal, "malignant", "non-malignant")
  seurat.obj[["malignant_bin"]] <- ifelse(is_mal, 1L, 0L)  # 跨脚本契约(P0-①)

  # 各患者恶性比例汇总表（生物学验证，P2-⑪）
  mal_tbl <- data.frame(
    patient   = unique(seurat.obj$patient),
    epi_total = NA_integer_,
    malignant = NA_integer_,
    prop_mal  = NA_real_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(mal_tbl))) {
    pat     <- mal_tbl$patient[i]
    cells_p <- colnames(seurat.obj)[seurat.obj$patient == pat]
    ep_p    <- intersect(cells_p, epi_obs)
    mal_p   <- intersect(cells_p, malignant_cells)
    mal_tbl$epi_total[i] <- length(ep_p)
    mal_tbl$malignant[i]   <- length(mal_p)
    mal_tbl$prop_mal[i]    <- ifelse(length(ep_p) > 0,
                                     round(length(mal_p) / length(ep_p), 3), 0)
  }
  write.csv(mal_tbl, file.path(tab_dir, "Table_malignant_by_patient.csv"),
            row.names = FALSE)
  FlushLog("  各患者恶性比例表已保存: Table_malignant_by_patient.csv\n")
  print(mal_tbl)

  # UMAP_malignant 出图统一在 PART 7 末尾（7.5 节）绘制：无论重算或跳过，
  # malignant 列均已写入 seurat.obj，保证图与当前数据一致。

  # 保存 PART 7 checkpoint（供跳过时恢复 malignant 列）；统一走 CheckSkip/MarkDone 链，
  # 不再单独写 infercnv_flag 文本旗标——否则删 rds 留 txt 会导致 malignant 列静默未定义。
  MarkDone(list(
    malignant_cells = malignant_cells,
    epi_obs = epi_obs
  ), "part07_infercnv")

  MemReport("InferCNV 完成（修复版）")

} else {
  FlushLog("  InferCNV 检查点已存在，跳过\n")
  # 从 PART 7 checkpoint 恢复 malignant 列
  # CheckSkip("part07_infercnv") 为真已保证 checkpoint 存在，直接读取
  part07_ckpt <- file.path(tab_dir, "checkpoint_part07_infercnv.rds")
  infercnv_res <- readRDS(part07_ckpt)
  is_mal  <- is.element(colnames(seurat.obj), infercnv_res$malignant_cells)
  seurat.obj[["malignant"]]     <- ifelse(is_mal, "malignant", "non-malignant")
  seurat.obj[["malignant_bin"]] <- ifelse(is_mal, 1L, 0L)  # P0-①
  FlushLog(sprintf("  恢复 malignant 列: %d malignant / %d total\n",
                    sum(is_mal), ncol(seurat.obj)))
}

# ----------------------------------------------------------------
# 7.5 UMAP_malignant 可视化
#   重算与跳过两个分支均已把 malignant 列写入 seurat.obj，故统一在此绘制，
#   保证图始终反映当前数据。
#   注意："空白 PDF" 曾被误判为参数大小写问题；实机验证显示 pdf(width=,height=)
#   小写正确，根因是点尺寸/透明度过低，在 PDF 整页缩放下几乎不可见。
#   红盖灰修复（2026-08-04）：恶性用空心红环(shape=21, fill=NA)，重叠处灰底始终透出；
#   非恶性灰底实心(shape=16)。灰底先画、红环后画置顶，密集肿瘤簇不再饱和成实心红块。
#   Seurat v5 中 obj[["malignant"]] 返回 data.frame，须用 as.vector(obj$malignant)。
# ----------------------------------------------------------------
umap_coords <- as.data.frame(Embeddings(seurat.obj, reduction = "umap"))
umap_coords$malignant <- as.vector(seurat.obj$malignant)   # 强制为向量

mal_vals <- unique(umap_coords$malignant)
if (!all(mal_vals %in% c("malignant", "non-malignant"))) {
  stop("UMAP_malignant: malignant 列出现意外值: ", paste(setdiff(mal_vals, c("malignant", "non-malignant")), collapse = ", "))
}
n_mal <- sum(umap_coords$malignant == "malignant")
n_ref <- sum(umap_coords$malignant == "non-malignant")
if (n_mal == 0L || n_ref == 0L) {
  stop(sprintf("UMAP_malignant: 子集为空 (malignant=%d, non-malignant=%d)", n_mal, n_ref))
}

p.mal <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = malignant)) +
  # 非恶性（灰）：先画、作为灰底；size/alpha 适中，保证整页 PDF 可见
  geom_point(data = umap_coords[umap_coords$malignant == "non-malignant", ],
             shape = 16, alpha = 0.55, size = 0.55) +
  # 恶性（红）：后画置顶；用空心圆(shape=21, fill=NA)——重叠处灰底始终可见，
  # 避免实心红块完全覆盖灰色；stroke 较细、size 较小，防止在密集肿瘤簇饱和成 blob
  geom_point(data = umap_coords[umap_coords$malignant == "malignant", ],
             shape = 21, fill = NA, alpha = 0.85,
             size = 0.60, stroke = 0.35) +
  scale_color_manual(name = NULL,
                     values = c("non-malignant" = "grey70", "malignant" = "#D62728")) +
  theme_bw() +
  theme(aspect.ratio = 1) + sci_theme +
  theme(legend.position = "right") +
  labs(title = sprintf("UMAP_malignant (n=%d, %.1f%% malignant)",
                       nrow(umap_coords),
                       100 * n_mal / nrow(umap_coords)))
sci_save(p.mal, "UMAP_malignant", width = 8, height = 6)
FlushLog(sprintf("  [UMAP_malignant] 已出图: %d malignant / %d total\n",
                 n_mal, nrow(umap_coords)))

# ================================================================
# PART 7b: InferCNV 热图可视化（降采样方案B）
# 目的: 生成每患者 InferCNV CNV 热图，验证恶性 vs 参考差异
# 方法: 降采样（最多 1000 恶性 + 500 参考），pheatmap 绘图
# 数据: infercnv_out_*/（本机实际安装版落盘 17_HMM_pred*.infercnv_obj；canonical 1.23.0 为 18_HMM_pred*），由 get_infercnv_hmm_states() 自适应提取离散态
# ================================================================

if (!CheckSkip("part07b_heatmap")) {
  FlushLog("========== PART 7b: InferCNV 热图 ==========\n")

  library(pheatmap)  # 热图绘制

  # 降采样参数
  max_mal     <- heatmap_max_mal     # 恶性细胞最多取 1000（顶部常量）
  max_ref     <- heatmap_max_ref     # 参考细胞最多取 500
  n_top_genes <- heatmap_n_top_genes # 选取 top 2000 高变基因（局部别名，避免遮蔽全局 n_hvg，P2-⑨）

  patients <- naive_patients

  n_heatmap_done <- 0L
  for (p in patients) {
    FlushLog(sprintf("\n  --- 患者 %s 热图 ---\n", p))

    # 动态提取 HMM 离散态矩阵（基因×细胞，值∈{1..6}）；文件名随 infercnv 版本变化
    # （1.23.0: 18_HMM_pred*.infercnv_obj；旧版: 17_HMM_*.infercnv_obj），由 helper 自适应定位
    out_dir_p <- file.path(tab_dir,
      paste0("infercnv_out_", gsub(" ", "_", p)))
    # PART 7 已为本患者产出 HMM 对象；缺失/损坏即硬报错（不静默跳过，错误在运行时显式暴露）
    cnv_mat <- get_infercnv_hmm_states(out_dir_p)

    # 恶性/参考定义统一取自 seurat.obj$malignant（与 PART 7 metadata 一致，P1-④）
    # 不再依赖 infercnv cell_groupings 的分组标签（P2-⑭）
    mal_barcodes <- colnames(seurat.obj)[seurat.obj$malignant == "malignant"]
    # 参照 = 非恶性上皮（与 PART 7 同谱系，不再用免疫类型作参照 — 根因修复）
    ref_barcodes <- colnames(seurat.obj)[seurat.obj$celltype %in% epi_types &
                                         seurat.obj$malignant == "non-malignant"]
    # 取交集：仅保留该患者、且存在于本热图 CNV 矩阵中的细胞
    pat_cells <- colnames(cnv_mat)
    mal_barcodes <- intersect(mal_barcodes, pat_cells)
    ref_barcodes <- intersect(ref_barcodes, pat_cells)
    mal_barcodes <- unique(mal_barcodes)
    ref_barcodes <- unique(ref_barcodes)

    FlushLog(sprintf("    malignant: %d, reference: %d\n",
      length(mal_barcodes), length(ref_barcodes)))

    # 降采样（固定种子可复现，P1-⑤）
    set.seed(seed_mal_base + which(patients == p))
    if (length(mal_barcodes) > max_mal) {
      mal_barcodes <- sample(mal_barcodes, max_mal)
    }
    set.seed(seed_ref_base + which(patients == p))
    if (length(ref_barcodes) > max_ref) {
      ref_barcodes <- sample(ref_barcodes, max_ref)
    }

    # 提取 CNV 矩阵（已为离散态矩阵，无需再 as.matrix）
    # 仅保留目标细胞
    cells_use <- intersect(c(mal_barcodes, ref_barcodes), colnames(cnv_mat))
    if (length(cells_use) < min_cells_heatmap) {
      FlushLog(sprintf("    跳过: 可用细胞不足 (%d)\n", length(cells_use)))
      next
    }
    mat_use <- cnv_mat[, cells_use, drop = FALSE]

    # 选 top 高变基因
    gene_vars <- apply(mat_use, 1, var, na.rm = TRUE)
    top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(n_top_genes, nrow(mat_use))]
    mat_use <- mat_use[top_genes, , drop = FALSE]

    # 以非恶性上皮的均值居中（减参照向量，非行均值：行均值会令 CNV 谱同质的矩阵
    # 退化为全 0，致 pheatmap generate_breaks 产生重复断点而 cut.default 报错）
    ref_in_mat <- intersect(ref_barcodes, colnames(mat_use))
    if (length(ref_in_mat) > 0) {
      ref_mean <- rowMeans(mat_use[, ref_in_mat, drop = FALSE], na.rm = TRUE)
      mat_use  <- mat_use - ref_mean
    } else {
      # 全部上皮均判恶性（无正常上皮可作对照）：不居中，直接画原始 HMM 离散态
      FlushLog(sprintf("    无非恶性上皮参照，直接绘制原始 HMM 离散态（不做居中）\n"))
    }

    # 退化保护：若矩阵所有值相同（无变异），pheatmap generate_breaks 会生成重复断点致
    # cut.default 报错；此时热图无信息量，显式跳过并告警（非静默吞错）
    if (min(mat_use) == max(mat_use)) {
      FlushLog(sprintf("    跳过: 矩阵无变异 (min==max=%s)，无法绘制热图\n", min(mat_use)))
      next
    }

    # 行/列注释
    cell_type <- ifelse(colnames(mat_use) %in% mal_barcodes,
                        "malignant", "reference")
    ann_col <- data.frame(type = cell_type)
    rownames(ann_col) <- colnames(mat_use)
    ann_colors <- list(type = c(malignant = "#E64B35",
                                reference = "#4DBBD5"))

    # 绘制热图
    pdf_file <- file.path(fig_dir,
      paste0("InferCNV_heatmap_", gsub(" ", "_", p), ".pdf"))
    pheatmap(mat_use,
      cluster_rows    = TRUE,
      cluster_cols    = TRUE,
      annotation_col  = ann_col,
      annotation_colors = ann_colors,
      show_rownames   = FALSE,
      show_colnames   = FALSE,
      fontsize        = 8,
      main            = paste0("InferCNV: ", p),
      filename        = pdf_file,
      width           = 12,
      height          = 8
    )
    FlushLog(sprintf("    热图已保存: %s\n", basename(pdf_file)))
    n_heatmap_done <- n_heatmap_done + 1L
  }

  if (n_heatmap_done == 0L) {
    FlushLog("  WARNING: 未生成任何热图（infercnv_obj 均缺失？），不写完成标记，下次将重试\n")
  } else {
    MarkDone(TRUE, "part07b_heatmap")   # 统一续跑旗标（替代 part07b_heatmap_done.txt）
    FlushLog(sprintf("  共生成 %d 张 InferCNV 热图\n", n_heatmap_done))
  }
  MemReport("InferCNV 热图完成")
}

# ================================================================

if (!CheckSkip("part08_subtype")) {

  FlushLog("========== PART 8: 分层亚型细分（逐区室重聚类） ==========\n")
  FlushLog("  贴合原文 Che 2021：for each compartment, re-centered, scaled, normalized, re-clustered\n")
  FlushLog("  原文 compartment = 5 大主要细胞类型：T/NK, B, Myeloid, CAFs(基质), Epithelial；各自独立重聚类得亚型\n")

  # 原文 5 大 compartment：每个 compartment 用一组主 celltype 值选取其全部细胞
  # （Myeloid 为整体区室含 TAM/DC/mono/MDSC/mast；B 区室含 plasma；T/NK 合并）
  compartment_map <- list(
    "T_NK"            = c("T cells", "NK cells"),
    "B"               = c("B cells", "Plasma cells"),
    "Myeloid"         = c("Mono/Macrophage", "DC", "pDC", "Mast cells"),
    "CAFs_mesenchyme" = "CAFs",
    "Epithelial"      = "Epithelial/Tumor"
  )

  for (ctype in names(compartment_map)) {
    matching_types <- compartment_map[[ctype]]
    cells_sub <- colnames(seurat.obj)[seurat.obj$celltype %in% matching_types]
    if (length(cells_sub) < min_cells_subtype) {
      FlushLog(sprintf("  [%s] < %d cells, 跳过子分\n", ctype, min_cells_subtype))
      next
    }
    obj_sub <- subset(seurat.obj, cells = cells_sub)
    # 原文：re-centered, scaled, normalized, re-clustered
    obj_sub <- NormalizeData(obj_sub, verbose = FALSE)
    obj_sub <- FindVariableFeatures(obj_sub, selection.method = "dispersion",
                                    mean.cutoff = c(0.0125, 3),
                                    dispersion.cutoff = 0.5,
                                    nfeatures = n_hvg, verbose = FALSE)
    obj_sub <- ScaleData(obj_sub, vars.to.regress = c("nCount_RNA", "percent.mt"),
                         verbose = FALSE)
    obj_sub <- RunPCA(obj_sub, npcs = npcs_total, verbose = FALSE)
    # 跨患者则 Harmony（与全局整合一致），否则直接用子集 PCA
    if (length(unique(obj_sub$patient)) > 1L) {
      obj_sub <- RunHarmony(obj_sub, group.by.vars = "patient",
                            reduction.save = "harmony")
      red_use <- "harmony"
    } else {
      red_use <- "pca"
    }
    dim_use <- 1:min(npcs_use, npcs_total)
    obj_sub <- RunUMAP(obj_sub, reduction = red_use, dims = dim_use, verbose = FALSE)
    obj_sub <- FindNeighbors(obj_sub, reduction = red_use, dims = dim_use,
                             verbose = FALSE)
    obj_sub <- FindClusters(obj_sub, resolution = subcluster_res)

    n_sub <- length(levels(obj_sub@meta.data$seurat_clusters))
    FlushLog(sprintf("  [%s] → %d sub-clusters\n", ctype, n_sub))

    p_sub <- DimPlot(obj_sub, reduction = "umap",
      group.by = "seurat_clusters", label = TRUE) +
      theme(aspect.ratio = 1) + sci_theme +
      labs(title = paste("Subtypes: ", ctype, sep = ""))
    subtype_name <- gsub("/", "_", gsub(" ", "_", ctype))
    sci_save(p_sub, paste("UMAP_subtype_", subtype_name, sep = ""),
              width = 8, height = 6)

    # 写回 subtype 列（不改变主 celltype，符合原文：重聚类产出亚型）
    new_col <- paste0("subtype_", subtype_name)
    meta_sub <- data.frame(row.names = colnames(obj_sub),
                           tmp_clusters = as.character(obj_sub$seurat_clusters))
    colnames(meta_sub) <- new_col
    seurat.obj <- AddMetaData(seurat.obj, meta_sub)

    # 基质 compartment（CAFs）：用原文 marker 区分内皮/肌成纤维（均为亚型，不升级主类型）
    # 原文：10 mesenchymal clusters (1 endothelial, 6 CAFs, 3 myofibroblasts)
    if (ctype == "CAFs_mesenchyme") {
      dat_sub <- GetAssayData(obj_sub, layer = "data")
      sc <- as.character(obj_sub$seurat_clusters)
      # 重聚类子簇 id：保留基质区室的亚型细分结构（贴合原文"重聚类得 10 mesenchymal"）
      meta_stromal_id <- data.frame(row.names = colnames(obj_sub),
                                    subtype_stromal_id = sc)
      seurat.obj <- AddMetaData(seurat.obj, meta_stromal_id)
      # 簇级 marker 判定内皮/肌成纤维（复现原文：10 mesenchymal clusters 中 1 为 endothelial，按簇识别而非逐细胞）
      # 原文 Endothelial 标记 CLDN5/CDH5/VWF；CAFs 标记 FAP/ACTA2/TAGLN
      endo_genes  <- c("CLDN5", "CDH5", "VWF")
      fibro_genes <- c("FAP", "ACTA2", "TAGLN")
      endo_genes  <- endo_genes[endo_genes  %in% rownames(obj_sub)]
      fibro_genes <- fibro_genes[fibro_genes %in% rownames(obj_sub)]
      sc <- as.character(obj_sub$seurat_clusters)
      cl_endo  <- tapply(colMeans(dat_sub[endo_genes, , drop = FALSE]), sc, mean)
      cl_fibro <- tapply(colMeans(dat_sub[fibro_genes, , drop = FALSE]), sc, mean)
      cl_label <- ifelse(cl_endo > cl_fibro, "Endothelial", "CAF_myofibroblast")
      cell_label <- cl_label[sc]
      meta_stromal <- data.frame(row.names = colnames(obj_sub),
                                 subtype_stromal = cell_label)
      seurat.obj <- AddMetaData(seurat.obj, meta_stromal)
      n_endo <- sum(cell_label == "Endothelial")
      FlushLog(sprintf("  基质亚型(簇级 marker, res=1.0): Endothelial=%d, CAF/肌成纤维=%d (主类型仍 CAFs)\n",
                       n_endo, length(cell_label) - n_endo))
    }
  }

  MarkDone(seurat.obj, "part08_subtype")
  MemReport("亚型细分完成")

}

# ================================================================
# PART 9: 保存最终结果
# ================================================================

FlushLog("========== PART 9: 保存最终结果 ==========\n")

saveRDS(seurat.obj, file.path(tab_dir, "seurat_obj_annotated.rds"))

# 计数汇总
tissue_stat <- seurat.obj@meta.data %>%
  group_by(tissue) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(proportion = round(n / sum(n) * 100, 1))
write.csv(tissue_stat,
  file = file.path(tab_dir, "Table1_tissue_cell_counts.csv"),
  row.names = FALSE)

celltype_stat <- seurat.obj@meta.data %>%
  group_by(tissue, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = round(n / sum(n) * 100, 1))
write.csv(celltype_stat,
  file = file.path(tab_dir, "Table2_celltype_by_tissue.csv"),
  row.names = FALSE)

cluster_stat <- seurat.obj@meta.data %>%
  group_by(seurat_clusters, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(as.numeric(seurat_clusters))
write.csv(cluster_stat,
  file = file.path(tab_dir, "Table4_cluster_cell_counts.csv"),
  row.names = FALSE)

FlushLog(sprintf("========== Step 1 全部完成 ==========\n"))
FlushLog(sprintf("  最终: %d 细胞, %d clusters\n",
  ncol(seurat.obj),
  length(levels(seurat.obj$seurat_clusters))))
MemReport("全部完成")
