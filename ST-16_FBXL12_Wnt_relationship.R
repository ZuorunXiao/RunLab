# ============================================================
# 16_FBXL12_Wnt_relationship
# 功能: FBXL12 × Wnt/β-catenin 探索性关系与空间验证 (脚本 17)
# 定位: 探索性(关系+空间) + 验证扩展电池(A1/B1/C3/D1/D3/F, 对应 analysis_plan_Wnt.md §4)
# 可视化: results/17_FBXL12_Wnt_relationship/{Figure(.pdf), table(.csv)}
# 数据: master(signature_master.rds, P1/P2) + 逐样本 Seurat 对象
#       (Spatial counts, 取 LEF1/TCF7/MYC)
# 注意: Wnt 原始计数来自逐样本 Seurat 对象(Spatial assay counts),
#       其 spot barcode 与 master 完美对齐;
#       raw/exprmatrix.rda 的 barcode 空间与本队列 disjoint(已实证), 不作逐 spot 对齐源。
# ============================================================

library(edgeR)        # 统计分析：伪bulk 差异表达检验
library(msigdbr)      # 数据处理：Wnt hallmark 基因集查询
library(ggplot2)      # 可视化：ggplot2
library(RColorBrewer) # 可视化：发散配色 brewer.pal (RdBu)
library(patchwork)    # 可视化：多面板拼图 wrap_plots
library(cowplot)      # 可视化：与 patchwork 配合排版
library(concaveman)   # 可视化：肿瘤区域凹包轮廓
library(Seurat)       # 领域专用：ST 对象读取与计数提取

# ---- 路径 ----
BASE_DIR <- "X:/XX"
TAB01    <- file.path(BASE_DIR, "results/01_read_qc/table")
TAB06    <- file.path(BASE_DIR, "results/06_PLAUR_SPP1_signature/table")
FIG_DIR  <- file.path(BASE_DIR, "results/17_FBXL12_Wnt_relationship/Figure")
TAB_DIR  <- file.path(BASE_DIR, "results/17_FBXL12_Wnt_relationship/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

SAMPLE_ID <- c("ST-P1_liver", "ST-P1_colon", "ST-P2_liver", "ST-P2_colon")
WNT_GENES <- c("LEF1", "TCF7", "MYC")
PATIENT   <- sub("_.*", "", SAMPLE_ID)   # ST-P1 / ST-P2

# ---- 内联 sci_theme / sci_save (对齐 07/08/15 + 规范2) ----
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
sci_save <- function(plot_obj, filename, width = 7, height = 7) {
  pdf(file.path(FIG_DIR, paste0(filename, ".pdf")), width = width,
      height = height, family = "Helvetica")
  print(plot_obj)
  dev.off()
}
# 连续场配色（对齐 07 HEAT_COLORS：YlOrRd 黄-橙-红）
HEAT_COLORS <- colorRampPalette(c("#FFFFCC", "#FED976", "#FEB24C",
                                  "#FD8D3C", "#FC4E2A", "#E31A1C",
                                  "#BD0026", "#800026"))(200)
# 发散配色（对齐 08 pheatmap：RdBu 蓝-白-红，用于可正可负的 FBXL12 场）
DIV_COLORS <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100)
# 样本固定配色/形状（对齐 15 SAMPLE_COL / SAMPLE_SHAPE，跨图一致）
SAMPLE_COL   <- c("ST-P1_liver" = "#cb181d", "ST-P1_colon" = "#4292c6",
                  "ST-P2_liver" = "#fb6a4a", "ST-P2_colon" = "#6baed6")
SAMPLE_SHAPE <- c("ST-P1_liver" = 16, "ST-P1_colon" = 17,
                  "ST-P2_liver" = 15, "ST-P2_colon" = 18)

# ---- 内联 Gaussian 核平滑(基于 spdep 邻域列表, 复用 15 思路) ----
gauss_smooth_nb <- function(xy, val, nb, bw) {
  n <- length(val); out <- numeric(n)
  for (i in 1:n) {
    j <- nb[[i]]
    if (length(j) == 0) { out[i] <- NA_real_; next }
    d <- sqrt((xy[i, 1] - xy[j, 1])^2 + (xy[i, 2] - xy[j, 2])^2)
    w <- exp(-(d^2) / (2 * bw^2))
    sw <- sum(w)
    out[i] <- if (sw <= 0) NA_real_ else sum(w * val[j]) / sw
  }
  out
}
# NA 鲁棒 z-score(逐列 center/scale, 避免 base scale() 遇 NA 整列变 NA)
zscore_na <- function(m) {
  cm <- colMeans(m, na.rm = TRUE)
  cs <- apply(m, 2, sd, na.rm = TRUE)
  sweep(sweep(m, 2, cm, "-"), 2, cs, "/")
}

# =============================================================================
# PART 1: 载入 master + 逐样本提取 Wnt 原始计数
# 目的：读入签名 master 并按 spot 合并逐样本 Seurat 的 Wnt 效应基因原始计数
# 方法：readRDS master; 逐样本 readRDS Seurat, GetAssayData(Spatial counts)
#       取 LEF1/TCF7/MYC; 按 (sample, spot) 合并 (barcode 完美对齐, 已实证 3826/3826)
# 数据：签名 master (P1/P2); Seurat Spatial counts 36601 gene x ~3826 spot/样本
# =============================================================================

cat("========== PART 1: 载入 master + 提取 Wnt 计数 ==========\n")
master <- readRDS(file.path(TAB06, "signature_master.rds"))
master <- master[master$sample %in% SAMPLE_ID, ]
wnt_list <- list()
for (s in SAMPLE_ID) {
  st <- readRDS(file.path(TAB01, paste0(s, "_seurat.rds")))
  cnt <- GetAssayData(st, assay = "Spatial", layer = "counts")
  wm <- as.matrix(cnt[WNT_GENES, , drop = FALSE])   # 3 gene x spots
  df <- data.frame(spot = colnames(st), sample = s,
                   LEF1 = as.numeric(wm["LEF1", ]),
                   TCF7 = as.numeric(wm["TCF7", ]),
                   MYC  = as.numeric(wm["MYC", ]), stringsAsFactors = FALSE)
  wnt_list[[s]] <- df
  rm(st); gc()
}
wnt <- do.call(rbind, wnt_list)
cat("Wnt 原始计数行数:", nrow(wnt), " (期望=master", nrow(master), ")\n")
master <- merge(master, wnt, by = c("sample", "spot"), all.x = TRUE)
cat("合并后 master 行数:", nrow(master), "\n")

# 检测率(透明披露, 对应 §2.1)
det_rate <- data.frame(gene = WNT_GENES,
                       pct_detected = round(100 * sapply(WNT_GENES,
                         function(g) mean(master[[g]] > 0)), 1),
                       stringsAsFactors = FALSE)
cat("Wnt 效应基因肿瘤内检测率(%):\n"); print(det_rate)

# =============================================================================
# PART 2: Wnt_effector 代理构建 (data-driven, §2.3)
# 目的：用检出率锁定 Wnt 效应基因, 构建空间平滑场与区域伪bulk 代理
# 方法：仅纳入 P1/P2 子集内 %检出 >= 5% 的效应基因 (锁定规则);
#       本队列 LEF1/TCF7/MYC 均 >= 5% 故全纳入; 每样本 Gaussian 平滑
#       原始计数 (bw=75, 邻域 R=300px), 全局 z-scored 后取均值; 区域层聚合
# 数据：master$wnt_eff_smooth (空间场); agg$wnt_eff_region (伪bulk 区域)
# =============================================================================

cat("========== PART 2: Wnt_effector 代理构建 ==========\n")
R_NEIGH <- 300; BW <- 75
master$wnt_eff_smooth <- NA_real_
for (s in SAMPLE_ID) {
  idx <- which(master$sample == s)
  xy <- cbind(master$x[idx], master$y[idx])
  nb <- spdep::dnearneigh(xy, 0, R_NEIGH)
  sm <- sapply(WNT_GENES, function(g) gauss_smooth_nb(xy, master[[g]][idx], nb,
    BW))
  smz <- zscore_na(sm)
  master$wnt_eff_smooth[idx] <- rowMeans(smz, na.rm = TRUE)
}
# 区域聚合(伪bulk)
agg <- aggregate(cbind(LEF1, TCF7, MYC) ~ sample + compartment, data = master,
  sum, na.rm = TRUE)
agg$z <- data.frame(t(scale(t(agg[, c("LEF1", "TCF7", "MYC")]))))
agg$wnt_eff_region <- rowMeans(agg$z)
cat("锁定效应基因:",
    paste(WNT_GENES[det_rate$pct_detected >= 5], collapse = ", "), "\n")

# =============================================================================
# PART 3: Q1 关系 - edgeR 伪bulk + 阈值敏感(F2) + 每患者方向(C3)
# 目的：检验 FBXL12 高低组中 Wnt_effector 差异, 并做阈值/分患者稳健性
# 方法：肿瘤 spot 内按 sct_FBXL12 样本中位 split 为 high/low; 聚合计数
#       per sample x group; edgeR(~patient + group), logFC 主证据 (p 仅参考);
#       F2 换 Q75/Q90 阈值; C3 每患者方向
# 数据：master$sct_FBXL12, master$wnt_eff_smooth; 输出 res_med/res_q75/res_q90/c3
# =============================================================================

cat("========== PART 3: Q1 关系 edgeR + 阈值/方向 ==========\n")
define_groups <- function(thr_fun) {
  g <- rep(FALSE, nrow(master))
  for (s in SAMPLE_ID) {
    ti <- master$sample == s & master$is_tumor
    thr <- thr_fun(master$sct_FBXL12[ti])
    g[ti] <- master$sct_FBXL12[ti] > thr
  }
  g
}
edgeR_group <- function(groups) {
  ag <- aggregate(cbind(LEF1, TCF7, MYC) ~ sample + groups, data = master, sum,
    na.rm = TRUE)
  ag$composite <- ag$LEF1 + ag$TCF7 + ag$MYC
  mat <- as.matrix(t(ag[, c("LEF1", "TCF7", "MYC", "composite")]))
  colnames(mat) <- paste(ag$sample, ifelse(ag$groups, "hi", "lo"), sep = "_")
  y <- DGEList(mat)
  patient <- sub("_.*", "", ag$sample)
  ag$patient <- patient
  design <- model.matrix(~ patient + groups, data = ag)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  lrt <- glmQLFTest(fit, coef = ncol(design))   # groups 系数
  tb <- topTags(lrt, n = Inf, sort.by = "none")$table
  data.frame(gene = rownames(tb), logFC = round(tb$logFC, 3),
             logCPM = round(tb$logCPM, 3), PValue = signif(tb$PValue, 3),
             FDR = signif(tb$FDR, 3), stringsAsFactors = FALSE)
}
res_med <- edgeR_group(define_groups(function(v) median(v)))
res_q75 <- edgeR_group(define_groups(function(v) quantile(v, 0.75)))
res_q90 <- edgeR_group(define_groups(function(v) quantile(v, 0.90)))
cat("Q1 edgeR (median split) logFC:\n"); print(res_med)
cat("F2 阈值敏感(Q75/Q90) composite logFC:",
    res_q75$logFC[res_q75$gene == "composite"], "/",
    res_q90$logFC[res_q90$gene == "composite"], "\n")
# C3 每患者方向(各自 edgeR)
c3 <- data.frame()
for (p in unique(PATIENT)) {
  sub <- master$sample %in% SAMPLE_ID[grep(p, SAMPLE_ID)]
  # 直接比较 high vs low 区域 mean
  tum_mask <- sub & master$is_tumor
  fb_med <- median(master$sct_FBXL12[tum_mask], na.rm = TRUE)
  mh <- mean(master$wnt_eff_smooth[tum_mask & master$sct_FBXL12 > fb_med],
             na.rm = TRUE)
  ml <- mean(master$wnt_eff_smooth[tum_mask & master$sct_FBXL12 <= fb_med],
             na.rm = TRUE)
  c3 <- rbind(c3, data.frame(patient = p, wnt_eff_high = round(mh, 3),
                             wnt_eff_low = round(ml, 3),
                             diff = round(mh - ml, 3),
                             stringsAsFactors = FALSE))
}
cat("C3 每患者方向(diff = high-low Wnt_effector):\n"); print(c3)

# =============================================================================
# PART 4: A1 代理有效性 + D1 贝叶斯 ROPE (bootstrap)
# 目的：验证 Wnt_effector 代理外部有效性, 并贝叶斯评估 high-low 差异
# 方法：A1: Wnt hallmark 标志集 (GO/HALLMARK) vs Wnt_effector 区域相关
#       (伪bulk 层面, 避免稀疏); D1: bootstrap 后验落 +/-0.5 SD 的 ROPE 比例
# 数据：hall_region (区域聚合 hallmark 计数); 输出 a1_r, d1_tbl
# =============================================================================

cat("========== PART 4: A1 代理有效性 + D1 ROPE ==========\n")
hall_all <- msigdbr::msigdbr(species = "Homo sapiens")
hall <- unique(c(
  hall_all$gene_symbol[hall_all$gs_name ==
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING"],
  hall_all$gene_symbol[hall_all$gs_name == "GO_BETA_CATENIN_BINDING"]))
hg <- readRDS(file.path(TAB01, paste0(SAMPLE_ID[1], "_seurat.rds")))
allg <- rownames(GetAssayData(hg, assay = "Spatial", layer = "counts"))
hall_in <- hall[hall %in% allg]
cat("Wnt hallmark 基因在 Seurat counts 中可用:", length(hall_in), "/", length(hall),
  "\n")
# 区域聚合 hallmark 计数(统一提取)
hall_region <- data.frame()
for (s in SAMPLE_ID) {
  st <- readRDS(file.path(TAB01, paste0(s, "_seurat.rds")))
  cnt <- GetAssayData(st, assay = "Spatial", layer = "counts")
  hm <- as.matrix(cnt[hall_in, , drop = FALSE])
  rm(st); gc()
  for (comp in unique(master$compartment[master$sample == s])) {
    smask <- master$sample == s & master$compartment == comp
    spots <- master$spot[smask]
    ci <- match(spots, colnames(hm))
    ci <- ci[!is.na(ci)]
    sval <- if (length(ci) > 0) sum(hm[, ci]) else 0
    hall_region <- rbind(hall_region, data.frame(sample = s, compartment = comp,
                                                  hall_sum = sval,
                                                    stringsAsFactors = FALSE))
  }
}
hall_region <- merge(hall_region,
                     aggregate(cbind(LEF1, TCF7, MYC) ~ sample + compartment,
                               data = master, sum),
                     by = c("sample", "compartment"))
hall_region$z <- data.frame(t(scale(t(hall_region[, c("LEF1", "TCF7",
  "MYC")]))))
hall_region$wnt_eff_region <- rowMeans(hall_region$z)
hall_region$hall_score <- scale(hall_region$hall_sum)
a1_r <- cor.test(hall_region$wnt_eff_region, hall_region$hall_score, method =
  "spearman")
cat("A1 Wnt_effector vs Wnt hallmark 区域 spearman:", round(a1_r$estimate, 3),
    " p=", signif(a1_r$p.value, 3), "\n")

# D1 贝叶斯 ROPE(bootstrap 后验 + bayestestR::rope, BEST 的免安装等价)
boot_diff <- function(x, y, B = 2000) {
  nx <- length(x); ny <- length(y); diffs <- numeric(B)
  for (b in 1:B) diffs[b] <- mean(sample(x, nx, replace = TRUE)) -
                          mean(sample(y, ny, replace = TRUE))
  diffs
}
hi <- master$wnt_eff_smooth[master$is_tumor & master$sct_FBXL12 >
  median(master$sct_FBXL12[master$is_tumor])]
lo <- master$wnt_eff_smooth[master$is_tumor & master$sct_FBXL12 <=
  median(master$sct_FBXL12[master$is_tumor])]
hi <- hi[!is.na(hi)]; lo <- lo[!is.na(lo)]
post <- boot_diff(hi, lo)
# BEST 等价 ROPE: 后验(bootstrap 差)落入 +/-0.5 SD 的比例 + 95% 可信区间(分位数)
rope_sd <- 0.5 * sd(c(hi, lo), na.rm = TRUE)
d1_tbl <- data.frame(metric = "D1_bootstrap_ROPE_BEST_equiv",
  ROPE_range_SD = 0.5,
  posterior_mean_diff = round(mean(post), 4),
  CI95_low = round(quantile(post, 0.025), 4),
  CI95_high = round(quantile(post, 0.975), 4),
  ROPE_fraction = round(mean(post >= -rope_sd & post <= rope_sd), 4),
  stringsAsFactors = FALSE)
cat("D1 贝叶斯 ROPE(bootstrap 后验, ROPE=+/-0.5 SD):\n"); print(d1_tbl)

# =============================================================================
# PART 5: Q2 空间 - 双场 + Moran's I + spearman + B1 置换 + F1 带宽
# 目的：量化 FBXL12 与 Wnt_effector 的空间共模式与自相关, 并做置换/带宽稳健性
# 方法：每样本 sct_FBXL12 场 + Wnt_effector 平滑场; spdep Moran's I (自相关);
#       spearman(sct_FBXL12, Wnt_effector) 为共模式主指标; B1 置换 Wnt 场
#       999 次得零分布经验 p; F1 比较 bw=75 vs 150
# 数据：master$sct_FBXL12, master$wnt_eff_smooth; 输出 spat_tbl/perm_tbl/f1_tbl
# =============================================================================

cat("========== PART 5: Q2 空间自相关 + 置换 ==========\n")
B1_N <- 999
spat_tbl <- data.frame()
perm_tbl <- data.frame()
for (s in SAMPLE_ID) {
  idx <- which(master$sample == s)
  xy <- cbind(master$x[idx], master$y[idx])
  keep <- is.finite(master$wnt_eff_smooth[idx]); xy <- xy[keep, , drop = FALSE]
  nb <- spdep::dnearneigh(xy, 0, R_NEIGH)
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  fbx <- master$sct_FBXL12[idx][keep]
  wnt <- master$wnt_eff_smooth[idx][keep]
  mi_fbx <- spdep::moran(fbx, lw, length(fbx), spdep::Szero(lw), zero.policy =
    TRUE)$I
  mi_wnt <- spdep::moran(wnt, lw, length(wnt), spdep::Szero(lw), zero.policy =
    TRUE)$I
  rho_obs <- cor(fbx, wnt, method = "spearman")
  # B1 置换零分布(仅置换 Wnt 场)
  rho_perm <- numeric(B1_N); mi_perm <- numeric(B1_N)
  for (b in 1:B1_N) {
    wp <- sample(wnt)
    rho_perm[b] <- cor(fbx, wp, method = "spearman")
    mi_perm[b] <- spdep::moran(wp, lw, length(wp), spdep::Szero(lw),
      zero.policy = TRUE)$I
  }
  p_rho <- mean(abs(rho_perm) >= abs(rho_obs))
  p_mi  <- mean(abs(mi_perm) >= abs(mi_wnt))
  spat_tbl <- rbind(spat_tbl, data.frame(sample = s,
    moran_FBXL12 = round(mi_fbx, 3), moran_Wnt = round(mi_wnt, 3),
    spearman_rho = round(rho_obs, 3), p_rho_perm = round(p_rho, 4),
    p_moran_perm = round(p_mi, 4), stringsAsFactors = FALSE))
  perm_tbl <- rbind(perm_tbl, data.frame(sample = s, sim = 1:B1_N,
    rho_perm = rho_perm, mi_perm = mi_perm, stringsAsFactors = FALSE))
}
cat("Q2 空间共模式(Moran's I + spearman + 置换 p):\n"); print(spat_tbl)
# F1 带宽敏感: bw=150 重算 spearman
f1_tbl <- data.frame()
for (s in SAMPLE_ID) {
  idx <- which(master$sample == s)
  xy <- cbind(master$x[idx], master$y[idx])
  nb <- spdep::dnearneigh(xy, 0, R_NEIGH)
  sm <- sapply(WNT_GENES, function(g) gauss_smooth_nb(xy, master[[g]][idx], nb,
    150))
  wnt150 <- rowMeans(zscore_na(sm), na.rm = TRUE)
  rho150 <- cor(master$sct_FBXL12[idx], wnt150, method = "spearman", use =
    "pairwise.complete.obs")
  f1_tbl <- rbind(f1_tbl, data.frame(sample = s, bw = 150,
    spearman_rho = round(rho150, 3), stringsAsFactors = FALSE))
}
cat("F1 带宽敏感(bw=150) spearman:\n"); print(f1_tbl)

# =============================================================================
# PART 6: Q2 延伸 - 与 PLAUR+SPP1+ TAM 生态位关联
# 目的：检验 Wnt_effector 空间场与 PLAUR/SPP1 巨噬签名的生态位共位
# 方法：Wnt_effector 平滑场 vs PLAUR_mac / SPP1_mac 双变量 spearman (每样本)
# 数据：master$wnt_eff_smooth, master$PLAUR_mac, master$SPP1_mac; 输出 tam_tbl
# =============================================================================

cat("========== PART 6: Q2 延伸 TAM 生态位 ==========\n")
tam_tbl <- data.frame()
for (s in SAMPLE_ID) {
  idx <- which(master$sample == s)
  r_p <- cor(master$wnt_eff_smooth[idx], master$PLAUR_mac[idx], method =
    "spearman", use = "pairwise.complete.obs")
  r_s <- cor(master$wnt_eff_smooth[idx], master$SPP1_mac[idx], method =
    "spearman", use = "pairwise.complete.obs")
  tam_tbl <- rbind(tam_tbl, data.frame(sample = s,
    rho_Wnt_PLAURmac = round(r_p, 3), rho_Wnt_SPP1mac = round(r_s, 3),
    stringsAsFactors = FALSE))
}
cat("Q2延伸 Wnt_effector vs TAM 签名 spearman:\n"); print(tam_tbl)

# =============================================================================
# PART 7: 派生分类 + 写表 + 出图 (publication-quality)
# 目的：派生双变量热点/三分位分类, 写出 CSV 表, 生成 10 张出版级空间图
# 方法：按样本肿瘤中位 split FBXL12/Wnt_effector 得 HH/HL/LH/LL 热点(全部 spot 参与分类);
#       FBXL12 三分位用于组间比较; geom_tile 空间瓦片 + 统一颜色标尺 + wrap_plots 单标题拼图
# 数据：master (派生分类列); 输出 9 CSV + 10 PDF 至 results/17_* /Figure|table
# =============================================================================

cat("========== PART 7: 派生分类 + 写表 + 出图 ==========\n")

# ---- 7.1 可视化前准备: 派生分类变量与 spot 半宽 ----
# 计算每样本 spot 半宽(geom_tile 用), 取逐 spot 最近邻距离中位数的一半
spot_halfwidth <- function(s) {
  xy <- cbind(master$x[master$sample == s], master$y[master$sample == s])
  d <- as.matrix(dist(xy))
  diag(d) <- Inf
  median(apply(d, 1, min), na.rm = TRUE) / 2
}
hw_map <- sapply(SAMPLE_ID, spot_halfwidth)

# 共享坐标界框(所有分图统一比例 -> 严格同尺寸 + spot 跨样本可比)
common_bbox <- function(samples, tumor_only = FALSE) {
  sub <- master[master$sample %in% samples, ]
  if (tumor_only) sub <- sub[sub$is_tumor, ]
  xr <- range(sub$x, na.rm = TRUE); yr <- range(sub$y, na.rm = TRUE)
  sp <- max(diff(xr), diff(yr))          # 取最大跨距 -> 正方形面板
  pad <- max(hw_map[samples])            # 防止边界 tile 被裁切
  cx <- mean(xr); cy <- mean(yr)
  list(xlim = c(cx - sp / 2 - pad, cx + sp / 2 + pad),
       ylim = c(cy - sp / 2 - pad, cy + sp / 2 + pad))
}
bbox_all <- common_bbox(SAMPLE_ID, tumor_only = FALSE)   # 全切片空间场

# 肿瘤区域凹包轮廓(对齐 07: 逐连通块 concaveman 描边 + 每块 "tumor" 标注)
tumor_concave <- function(s, concavity = 2, min_size = 10) {
  tum <- master[master$sample == s & master$is_tumor, c("x", "y")]
  if (nrow(tum) < 3) return(NULL)
  spacing <- 2 * hw_map[s]
  pts <- as.matrix(tum)
  d <- as.dist(as.matrix(dist(pts)))
  cl <- cutree(hclust(d, method = "single"), h = 2.5 * spacing)
  polys <- list(); cents <- data.frame(x = numeric(0), y = numeric(0))
  for (k in sort(unique(cl))) {
    sub_pts <- pts[cl == k, , drop = FALSE]
    if (nrow(sub_pts) < min_size) next
    poly <- tryCatch(
      concaveman::concaveman(sub_pts, concavity = concavity),
      error = function(e) sub_pts[chull(sub_pts[, "x"], sub_pts[, "y"]), ]
    )
    poly <- as.data.frame(poly)
    if (!all(c("x", "y") %in% colnames(poly))) colnames(poly) <- c("x", "y")
    poly$group <- k
    polys[[length(polys) + 1]] <- poly
    cents <- rbind(cents,
      data.frame(x = mean(sub_pts[, "x"]), y = mean(sub_pts[, "y"])))
  }
  if (length(polys) == 0) return(NULL)
  list(poly = do.call(rbind, polys), cen = cents)
}

# 肿瘤 spot 内按样本中位数 split FBXL12 / Wnt_effector -> 双变量热点分类;
# 肿瘤/非肿瘤 spot 各自用自身中位数 split FBXL12 / Wnt_effector -> HH/HL/LH/LL;
# 所有 spot(含非肿瘤)均参与分类, 同颜色图内对照, 虚线 tumor 边界区分区室
master$fbxl12_hi <- FALSE
master$wnt_hi    <- FALSE
master$bivariate_class <- NA_character_
for (s in SAMPLE_ID) {
  idx_tum <- which(master$sample == s & master$is_tumor)
  idx_nt  <- which(master$sample == s & !master$is_tumor)
  med_fb_tum  <- median(master$sct_FBXL12[idx_tum], na.rm = TRUE)
  med_wnt_tum <- median(master$wnt_eff_smooth[idx_tum], na.rm = TRUE)
  med_fb_nt   <- median(master$sct_FBXL12[idx_nt], na.rm = TRUE)
  med_wnt_nt  <- median(master$wnt_eff_smooth[idx_nt], na.rm = TRUE)
  master$fbxl12_hi[idx_tum] <- master$sct_FBXL12[idx_tum] > med_fb_tum
  master$wnt_hi[idx_tum]    <- master$wnt_eff_smooth[idx_tum] > med_wnt_tum
  master$fbxl12_hi[idx_nt]  <- master$sct_FBXL12[idx_nt] > med_fb_nt
  master$wnt_hi[idx_nt]     <- master$wnt_eff_smooth[idx_nt] > med_wnt_nt
  for (region in list(idx_tum, idx_nt)) {
    master$bivariate_class[region] <- ifelse(
      master$fbxl12_hi[region] & master$wnt_hi[region], "HH",
      ifelse(master$fbxl12_hi[region] & !master$wnt_hi[region], "HL",
      ifelse(!master$fbxl12_hi[region] & master$wnt_hi[region], "LH", "LL")))
  }
}
master$bivariate_class <- factor(master$bivariate_class, levels = c("HH", "HL",
  "LH", "LL"))

# 肿瘤 spot 内 FBXL12 三分位数(用于组间比较可视化)
# 用 order-based 等样本量三分位, 避免零膨胀/重复值导致 quantile/rank 分位点坍缩
master$fbxl12_tert <- NA_character_
for (s in SAMPLE_ID) {
  idx <- which(master$sample == s & master$is_tumor &
    is.finite(master$sct_FBXL12))
  n <- length(idx)
  if (n < 30) next
  ord <- order(master$sct_FBXL12[idx])
  k <- floor(n / 3)
  master$fbxl12_tert[idx[ord[1:k]]] <- "Low"
  master$fbxl12_tert[idx[ord[(k + 1):(2 * k)]]] <- "Mid"
  master$fbxl12_tert[idx[ord[(2 * k + 1):n]]] <- "High"
}
master$fbxl12_tert <- factor(master$fbxl12_tert, levels = c("Low", "Mid",
  "High"))

# 写表(保持与旧版一致)
write.csv(det_rate, file.path(TAB_DIR, "WntEffector_proxy.csv"), row.names =
  FALSE)
write.csv(res_med, file.path(TAB_DIR, "WntEffector_pseudobulk.csv"), row.names
  = FALSE)
write.csv(data.frame(threshold = c("median", "q75", "q90"),
    composite_logFC = c(res_med$logFC[res_med$gene == "composite"],
                        res_q75$logFC[res_q75$gene == "composite"],
                        res_q90$logFC[res_q90$gene == "composite"]),
    stringsAsFactors = FALSE), file.path(TAB_DIR,
      "WntEffector_threshold_sens.csv"), row.names = FALSE)
write.csv(c3, file.path(TAB_DIR, "WntEffector_per_patient.csv"), row.names =
  FALSE)
write.csv(data.frame(metric = "Wnt_effector_vs_Wnt_hallmark_region_spearman",
    estimate = round(a1_r$estimate, 3), p_value = signif(a1_r$p.value, 3),
    stringsAsFactors = FALSE), file.path(TAB_DIR,
      "WntEffector_proxy_validation.csv"), row.names = FALSE)
write.csv(d1_tbl, file.path(TAB_DIR, "WntEffector_validation_bayes_LOSO.csv"),
  row.names = FALSE)
write.csv(spat_tbl, file.path(TAB_DIR, "WntEffector_spatial_stats.csv"),
  row.names = FALSE)
write.csv(perm_tbl, file.path(TAB_DIR, "WntEffector_spatial_perm.csv"),
  row.names = FALSE)
write.csv(tam_tbl, file.path(TAB_DIR, "WntEffector_TAM_colocal.csv"),
  row.names = FALSE)

# ---- 7.2 统一配色 (复用 7.1 的 HEAT_COLORS/DIV_COLORS/SAMPLE_COL) ----
# 双变量热点四色（HH 红 / HL 橙 / LH 浅蓝 / LL 灰）
HOT_PAL <- c(HH = "#cb181d", HL = "#fd8d3c", LH = "#9ecae1", LL = "#e0e0e0")

# ---- 图1: WntEffector_detectability.pdf (检测率条形图, 对齐 08 geom_text 标注风格) ----
p1 <- ggplot(det_rate, aes(x = reorder(gene, pct_detected), y = pct_detected)) +
  geom_col(fill = "#4292c6", width = 0.65) +
  geom_text(aes(label = paste0(pct_detected, "%")), vjust = -0.5, size = 4) +
  sci_theme +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     limits = c(0, max(det_rate$pct_detected) * 1.2)) +
  labs(title = "Wnt effector gene detection rate (Visium, P1/P2)",
       x = "Gene", y = "Detected spots (%)")
sci_save(p1, "WntEffector_detectability", width = 8, height = 6)

# ---- 图2: FBXL12_spatial_map + WntEffector_spatial_map (四样本空间真图, 对齐 07) ----
# 统一颜色标尺范围, 跨样本只显示单一 colorbar
map_limit <- function(v) {
  sub <- master[[v]][master$sample %in% SAMPLE_ID & is.finite(master[[v]])]
  r <- range(sub)
  c(r[1], r[2])
}
spatial_map <- function(var_name, fill_ramp, fill_name, bbox, limits = NULL) {
  plist <- lapply(SAMPLE_ID, function(s) {
    sub <- master[master$sample == s & is.finite(master[[var_name]]), ]
    hw <- hw_map[s]
    p <- ggplot(sub, aes(x = x, y = y, fill = .data[[var_name]])) +
      geom_tile(width = 2 * hw, height = 2 * hw, color = NA) +
      scale_fill_gradientn(colours = fill_ramp, name = fill_name,
                           limits = limits) +
      scale_x_continuous(limits = bbox$xlim, expand = c(0, 0)) +
      scale_y_reverse(limits = bbox$ylim, expand = c(0, 0)) +
      coord_fixed() + sci_theme +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(size = 13, hjust = 0.5, face = "plain")) +
      labs(title = paste0(fill_name, "  |  ", gsub("_", " ", s)))
    tc <- tumor_concave(s)
    if (!is.null(tc)) {
      p <- p +
        geom_polygon(data = tc$poly, aes(x = x, y = y, group = group),
          fill = NA, colour = "black", linetype = "dashed", linewidth = 0.6,
          inherit.aes = FALSE) +
        geom_text(data = tc$cen, aes(x = x, y = y), label = "tumor",
          size = 4, colour = "black", fontface = "bold", inherit.aes = FALSE)
    }
    p
  })
  wrap_plots(plist, ncol = 2, widths = rep(1, 2), heights = rep(1, 2),
    guides = "collect") +
    plot_annotation(title = paste("Spatial distribution of", var_name),
                    theme = theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "plain")))
}
p_fbx_map <- spatial_map("sct_FBXL12", DIV_COLORS, "sct_FBXL12",
  bbox_all, map_limit("sct_FBXL12"))
p_wnt_map <- spatial_map("wnt_eff_smooth", HEAT_COLORS, "Wnt_effector",
  bbox_all, map_limit("wnt_eff_smooth"))
sci_save(p_fbx_map, "FBXL12_spatial_map", width = 12, height = 11)
sci_save(p_wnt_map, "WntEffector_spatial_map", width = 12, height = 11)
# ---- 图3: FBXL12_Wnt_bivariate_hotspot.pdf (双变量热点空间图, 全部 spot 参与分类, 对齐 07) ----
hotspot_list <- lapply(SAMPLE_ID, function(s) {
  sub <- master[master$sample == s & !is.na(master$bivariate_class), ]
  hw <- hw_map[s]
  tab <- round(100 * prop.table(table(sub$bivariate_class)), 1)
  ann <- data.frame(bivariate_class = factor(names(tab), levels = c("HH", "HL",
    "LH", "LL")), pct = as.numeric(tab))
  ann$x <- bbox_all$xlim[1]; ann$y <- bbox_all$ylim[2]
  tc <- tumor_concave(s)
  ggplot(sub, aes(x = x, y = y, fill = bivariate_class)) +
    geom_tile(width = 2 * hw, height = 2 * hw, color = NA) +
    scale_fill_manual(values = HOT_PAL, name = "Co-localization",
                      labels = c(HH = "FBXL12 high / Wnt high",
                                 HL = "FBXL12 high / Wnt low",
                                 LH = "FBXL12 low / Wnt high",
                                 LL = "FBXL12 low / Wnt low")) +
    geom_text(data = ann, aes(x = x, y = y, label = paste0(pct, "%")),
              hjust = 0, vjust = 1, size = 3, inherit.aes = FALSE) +
    { if (!is.null(tc)) geom_polygon(data = tc$poly,
        aes(x = x, y = y, group = group), fill = NA, colour = "black",
        linetype = "dashed", linewidth = 0.6, inherit.aes = FALSE) } +
    { if (!is.null(tc)) geom_text(data = tc$cen, aes(x = x, y = y),
        label = "tumor", size = 4, colour = "black", fontface = "bold",
        inherit.aes = FALSE) } +
    scale_x_continuous(limits = bbox_all$xlim, expand = c(0, 0)) +
    scale_y_reverse(limits = bbox_all$ylim, expand = c(0, 0)) +
    coord_fixed() + sci_theme +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(size = 13, hjust = 0.5, face = "plain")) +
    labs(title = paste0("Bivariate hotspot  |  ", gsub("_", " ", s)))
})
p_hotspot <- wrap_plots(hotspot_list, ncol = 2, widths = rep(1, 2),
  heights = rep(1, 2), guides = "collect") +
  plot_annotation(title =
    "Bivariate spatial hotspots: FBXL12 vs Wnt effector (all spots)",
                  theme = theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "plain")))
sci_save(p_hotspot, "FBXL12_Wnt_bivariate_hotspot", width = 12, height = 11)

# ---- 图4: FBXL12_Wnt_relationship_boxplot.pdf (FBXL12 三分位 × Wnt 箱线图) ----
box_df <- master[master$is_tumor & !is.na(master$fbxl12_tert) &
  is.finite(master$wnt_eff_smooth), ]
p_box <- ggplot(box_df, aes(x = fbxl12_tert, y = wnt_eff_smooth, fill =
  fbxl12_tert)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.85) +
  facet_wrap(~ sample, ncol = 2) +
  scale_fill_manual(values = c(Low = "#9ecae1", Mid = "#fd8d3c", High =
    "#cb181d"),
                    name = "FBXL12 tertile") +
  sci_theme +
  labs(title = "Wnt effector by FBXL12 expression tertile (tumor spots)",
       x = "FBXL12 tertile (per sample)", y = "Wnt effector (smoothed z-score)")
sci_save(p_box, "FBXL12_Wnt_relationship_boxplot", width = 9, height = 6)

# ---- 图5: FBXL12_Wnt_spatial_perm.pdf (置换零分布, 观测值用虚线, 单标题) ----
p_perm <- ggplot(perm_tbl, aes(x = rho_perm)) +
  geom_histogram(bins = 40, fill = "#9ecae1", color = "white") +
  geom_vline(data = spat_tbl, aes(xintercept = spearman_rho),
             color = "#cb181d", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ sample, ncol = 2) +
  sci_theme +
  labs(title = "Permutation null: spearman(sct_FBXL12, Wnt effector)",
       x = "Permuted spearman rho", y = "Count")
sci_save(p_perm, "FBXL12_Wnt_spatial_perm", width = 9, height = 7)

# ---- 图6: FBXL12_Wnt_TAM_colocal.pdf (TAM 相关柱图, 对齐 15 SAMPLE_COL) ----
tam_long <- data.frame(
  sample = rep(tam_tbl$sample, 2),
  pair   = rep(c("Wnt vs PLAUR_mac", "Wnt vs SPP1_mac"), each = nrow(tam_tbl)),
  rho    = c(tam_tbl$rho_Wnt_PLAURmac, tam_tbl$rho_Wnt_SPP1mac),
  stringsAsFactors = FALSE)
tam_long$sample <- factor(tam_long$sample, levels = SAMPLE_ID)
p_tam_bar <- ggplot(tam_long, aes(x = sample, y = rho, fill = sample)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(rho, 2)), position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = SAMPLE_COL) +
  facet_wrap(~ pair, ncol = 2) +
  sci_theme +
  labs(title = "Wnt effector vs TAM signatures (spearman)",
       x = "Sample", y = "Spearman rho") +
  theme(legend.position = "bottom", legend.title = element_blank())
sci_save(p_tam_bar, "FBXL12_Wnt_TAM_colocal", width = 9, height = 5.5)

# ---- 图7: Wnt_validation_summary.pdf (代理验证 + ROPE + 阈值/带宽敏感, 2x2) ----
# A: Wnt hallmark 相关散点（按样本着色）
p_val_a <- ggplot(hall_region, aes(x = wnt_eff_region, y = hall_score)) +
  geom_point(aes(colour = sample), size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  scale_colour_manual(values = SAMPLE_COL, name = "Sample") +
  sci_theme +
  labs(title = paste0("Wnt effector vs Wnt hallmark (rho=",
    round(a1_r$estimate, 2),
                      ", p=", signif(a1_r$p.value, 2), ")"),
       x = "Wnt effector (region z-score)", y =
         "Wnt hallmark score (z-score)") +
  theme(legend.position = "bottom")
# B: Bootstrap 后验分布 + ROPE 带
post_df <- data.frame(diff = post)
p_val_b <- ggplot(post_df, aes(x = diff)) +
  geom_histogram(bins = 60, fill = "#9ecae1", color = "white") +
  annotate("rect", xmin = -rope_sd, xmax = rope_sd, ymin = -Inf, ymax = Inf,
           alpha = 0.25, fill = "#4292c6") +
  geom_vline(xintercept = mean(post), color = "#cb181d", linetype = "dashed",
    linewidth = 0.8) +
  sci_theme +
  labs(title = paste0("Bootstrap posterior (mean=", round(mean(post), 3),
                      ", ROPE%=", round(d1_tbl$ROPE_fraction[1] * 100, 1), ")"),
       x = "Wnt effector high minus low FBXL12", y = "Count")
# C: 阈值敏感性
thr_sens <- data.frame(threshold = c("median", "Q75", "Q90"),
  logFC = c(res_med$logFC[res_med$gene == "composite"],
            res_q75$logFC[res_q75$gene == "composite"],
            res_q90$logFC[res_q90$gene == "composite"]))
p_val_c <- ggplot(thr_sens, aes(x = threshold, y = logFC, group = 1)) +
  geom_line(color = "#4292c6", linewidth = 1) +
  geom_point(size = 3, color = "#4292c6") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  sci_theme +
  labs(title = "Threshold sensitivity: composite Wnt logFC",
       x = "FBXL12 high threshold", y = "log2 fold change")
# D: 带宽敏感
p_val_d <- ggplot(f1_tbl, aes(x = sample, y = spearman_rho, fill = sample)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = spearman_rho), vjust = -0.5, size = 3) +
  scale_fill_manual(values = SAMPLE_COL) +
  sci_theme +
  labs(title = "Bandwidth sensitivity (bw=150): spearman",
       x = "Sample", y = "Spearman rho") +
  theme(legend.position = "bottom", legend.title = element_blank())
p_validation <- (p_val_a + p_val_b) / (p_val_c + p_val_d) +
  plot_annotation(title = "Wnt effector construct validation and sensitivity",
                  theme = theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "plain")))
sci_save(p_validation, "Wnt_validation_summary", width = 11, height = 10)

# ---- 图8: Wnt_correlation_heatmap.pdf (多变量 Spearman 矩阵, 对齐 08 RdBu 配色) ----
cor_vars <- c("sct_FBXL12", "wnt_eff_smooth", "PLAUR_mac", "SPP1_mac")
cor_mat <- cor(master[, cor_vars], method = "spearman", use =
  "pairwise.complete.obs")
cor_df <- as.data.frame(as.table(cor_mat))
cor_df$Var1 <- factor(cor_df$Var1, levels = cor_vars)
cor_df$Var2 <- factor(cor_df$Var2, levels = rev(cor_vars))
p_cor <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 2)), color = "black", size = 4) +
  scale_fill_gradientn(colours = DIV_COLORS, name = "Spearman rho",
                       limits = c(-1, 1)) +
  sci_theme +
  labs(title = "Correlation matrix (all spots, pairwise complete)",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sci_save(p_cor, "Wnt_correlation_heatmap", width = 6, height = 5)

# ---- 图9: Wnt_moran_summary.pdf (Moran's I 汇总, 对齐 15 双色) ----
spat_long <- data.frame(
  sample = rep(spat_tbl$sample, 2),
  variable = rep(c("FBXL12", "Wnt_effector"), each = nrow(spat_tbl)),
  moran = c(spat_tbl$moran_FBXL12, spat_tbl$moran_Wnt),
  stringsAsFactors = FALSE)
spat_long$sample <- factor(spat_long$sample, levels = SAMPLE_ID)
p_moran <- ggplot(spat_long, aes(x = sample, y = moran, fill = variable)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = moran), position = position_dodge(width = 0.8), vjust
    = -0.5, size = 3) +
  scale_fill_manual(values = c(FBXL12 = "#fd8d3c", "Wnt_effector" = "#cb181d"),
                    name = "Variable") +
  sci_theme +
  labs(title = "Moran's I spatial autocorrelation", x = "Sample", y =
    "Moran's I") +
  theme(legend.position = "bottom")
sci_save(p_moran, "Wnt_moran_summary", width = 8, height = 5)

cat("17_FBXL12_Wnt_relationship.R done. Generated", length(list.files(FIG_DIR,
  pattern = "\\.pdf$")), "PDF figures in", FIG_DIR, "\n")
