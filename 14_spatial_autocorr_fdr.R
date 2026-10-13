# =============================================================================
# 14_spatial_autocorr_fdr.R
# 功能: 统计校正层——(②) 空间自相关量化 + 空间块 bootstrap 校正 Cliff's delta
#       的置信区间; (④) 对 09/10/11/12 既有 p 值做跨家族 BH(FDR) 校正
# 定位: n=2(P1/P2 四切片)下提升结论严谨性, 不扩大样本、不改写既有脚本
# 可视化: MoranI_Signature.pdf (逐样本逐签名空间自相关)
#         SpatialBlockBootstrap_CI.pdf (naive vs 空间块 bootstrap 的 CI 对比)
# 数据: 输入 results/06.../table/signature_master.csv
#            results/09.../table/CompartmentComparisons.csv
#            results/10.../table/{MyeloidCoexpression,DoubleHighBoundary,
#                                  NearestNeighbourPermutation}.csv
#            results/11.../table/CellChatPLAUR.csv
#            results/12.../table/{MediationResults,MediationSensitivity,
#                                  PerSampleACME}.csv
#       输出 results/14.../Figure/ + results/14.../table/
# 依据: spot 级检验把 ~3000-3500 spot 当 iid, 但 Visium spot 空间自相关使
#       有效样本量远低于名义值, p 值反保守; 本脚本量化自相关并给出空间块
#       校正的 CI, 同时统一各脚本 p 值的 FDR 阈值
# 说明: 本脚本只读既有产物 + 重算, 不修改 00-13 任何文件
# =============================================================================

# --- 包加载(统计 -> 可视化) ----------------------------------------------
suppressMessages({
  library(ggplot2)  # 可视化
})

# --- 全局辅助函数(内联, 规范22 自给自足) -------------------------------
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
  invisible(dev.off())
}

# CliffDelta: 由 Mann-Whitney U 换算 Cliff's delta = 2*AUC - 1
# AUC = P(a > b) = U_a / (n1*n2), U_a = 组内 a 的秩和 - n1*(n1+1)/2
# 显式剔除非有限值(如 SPP1_mac 仅在巨噬 spot 有值, 其余为 NA),
# 分母用有限样本量; 避免 wilcox.test 静默去 NA 但分母仍含 NA 计数导致的
# 符号翻转(旧公式用全向量长度作分母 -> 比值趋零 -> delta 误为 -1)
CliffDelta <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  n1 <- length(a); n2 <- length(b)
  if (n1 < 1 || n2 < 1) return(NA_real_)
  ra <- rank(c(a, b))[seq_len(n1)]
  Ua <- sum(ra) - n1 * (n1 + 1) / 2
  2 * Ua / (n1 * n2) - 1
}

# BuildEdges: 距离阈值无向边表(上三角, j>i 去重), 仅依赖坐标
# 返回 list(i, j) 两向量; 邻接对称 -> 后续分子用 2*sum(z_i*z_j)
# 预分配整数向量 + 下标填充, 避免循环内向量追加的 O(n^2) 复制
BuildEdges <- function(xy, d.thresh) {
  n <- nrow(xy)
  cap <- n * 50
  ei <- integer(cap); ej <- integer(cap)
  cnt <- 0L
  for (i in 1:n) {
    idx.other <- (1:n)[-i]
    di <- sqrt((xy[idx.other, 1] - xy[i, 1])^2 +
               (xy[idx.other, 2] - xy[i, 2])^2)
    js <- idx.other[di <= d.thresh & idx.other > i]
    e <- length(js)
    if (e > 0L) {
      ei[(cnt + 1L):(cnt + e)] <- i
      ej[(cnt + 1L):(cnt + e)] <- js
      cnt <- cnt + e
    }
  }
  list(i = ei[1:cnt], j = ej[1:cnt])
}

# MoranI: 用边表矢量化计算 Moran's I 与置换伪 p(双尾), 避免逐 spot 循环
# 输入 x 单签名向量, edges 对应边表, S0 = 2*n_edges; 内部对 NA 子集自包含建边
# 返回 c(I, p)
MoranI <- function(x, xy, d.thresh, B) {
  keep <- !is.na(x)
  x <- x[keep]; xy <- xy[keep, , drop = FALSE]
  n <- length(x)
  if (n < 3) return(c(I = NA_real_, p = NA_real_))
  edges <- BuildEdges(xy, d.thresh)
  S0 <- 2 * length(edges$i)
  z <- x - mean(x)
  den <- sum(z^2)
  if (den == 0 || S0 == 0) return(c(I = NA_real_, p = NA_real_))
  num <- 2 * sum(z[edges$i] * z[edges$j])
  I.obs <- (n / S0) * num / den
  p <- 0
  for (b in 1:B) {
    zp <- sample(z)
    np <- 2 * sum(zp[edges$i] * zp[edges$j])
    Ip <- (n / S0) * np / den
    if (abs(Ip) >= abs(I.obs)) p <- p + 1
  }
  c(I = I.obs, p = (p + 1) / (B + 1))
}

# BlockBootCliff: 同一比较的 naive(spot 级) 与空间块 bootstrap 的 Cliff's
# delta 分布; 返回 list(naive, spatial) 两向量, 用于对比 CI 宽度
BlockBootCliff <- function(sig, g1, g2, master, block.size, B) {
  naive <- replicate(B, {
    a <- sample(master[[sig]][master$compartment == g1], replace = TRUE)
    b <- sample(master[[sig]][master$compartment == g2], replace = TRUE)
    CliffDelta(a, b)
  })
  spatial <- replicate(B, {
    idx <- c()
    for (s in unique(master$sample)) {
      d <- master[master$sample == s, c("x", "y")]
      bx <- floor(d$x / block.size); by <- floor(d$y / block.size)
      bid <- factor(paste(bx, by, sep = "_"))
      spots.by.block <- split(seq_len(nrow(d)), bid)
      blk.ids <- levels(bid)
      pick <- sample(blk.ids, replace = TRUE)
      sel <- which(master$sample == s)
      idx <- c(idx, sel[unlist(spots.by.block[pick])])
    }
    sub <- master[idx, ]
    CliffDelta(sub[[sig]][sub$compartment == g1],
               sub[[sig]][sub$compartment == g2])
  })
  list(naive = naive, spatial = spatial)
}

# --- 目录结构 -----------------------------------------------------------
BASE_DIR <- "X:/XX"
TAB06    <- file.path(BASE_DIR, "results/06_PLAUR_SPP1_signature/table")
TAB09    <- file.path(BASE_DIR, "results/09_PLAUR_SPP1_trends/table")
TAB10    <- file.path(BASE_DIR, "results/10_PLAUR_SPP1_colocalization/table")
TAB11    <- file.path(BASE_DIR, "results/11_PLAUR_SPP1_communication/table")
TAB12    <- file.path(BASE_DIR, "results/12_PLAUR_SPP1_mediation/table")
FIG_DIR  <- file.path(BASE_DIR, "results/14_spatial_autocorr_fdr/Figure")
TAB_DIR  <- file.path(BASE_DIR, "results/14_spatial_autocorr_fdr/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 分析常量 -----------------------------------------------------------
SIG_MI   <- c("sct_PLAU", "sct_PLAUR", "sct_SPP1", "SPP1_mac", "PLAUR_mac")
D_THRESH <- 300      # Moran's I 邻接距离阈值(数组单位, 非 um)
                    # 实测 spot 最近邻间距≈200 单位(1 间距≈200 单位),
                    # 200 恰好卡在边界使 colon 近邻数=0 -> I 误为 0;
                    # 取 300(≈1.5 间距)稳健捕获两样本第一邻环
N_PERM   <- 999      # Moran's I 置换次数
N_BOOT   <- 1000     # 块 bootstrap 次数
BLOCK    <- 500      # 空间块边长(um)
SEED     <- 42

# =============================================================================
# PART 1: 逐样本逐签名 Moran's I(空间自相关量化)
# 目的: 证明 spot 表达存在强空间自相关 -> 名义样本量虚高, 需校正
# 方法: 距离阈值邻接 + 置换伪 p(近邻表仅依赖坐标, 同样本多签名复用)
# 数据: 输出 MoranI_by_signature_sample.csv
# =============================================================================

cat("========== PART 1: Moran's I spatial autocorrelation ==========\n")

set.seed(SEED)
master <- read.csv(file.path(TAB06, "signature_master.csv"))

mi.rows <- list()
for (s in unique(master$sample)) {
  xy <- as.matrix(master[master$sample == s, c("x", "y")])
  for (sig in SIG_MI) {
    x <- master[[sig]][master$sample == s]
    res <- MoranI(x, xy, D_THRESH, N_PERM)
    mi.rows[[length(mi.rows) + 1]] <- data.frame(
      sample = s, signature = sig, morans_I = round(res[1], 4),
      p_perm = signif(res[2], 3), stringsAsFactors = FALSE)
  }
}
mi.tab <- do.call(rbind, mi.rows)
write.csv(mi.tab, file.path(TAB_DIR, "MoranI_by_signature_sample.csv"),
          row.names = FALSE)
print(mi.tab)

p.mi <- ggplot(mi.tab, aes(x = signature, y = morans_I,
                           fill = p_perm < 0.05)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ sample, ncol = 2) +
  scale_fill_manual(values = c("FALSE" = "#9ecae1", "TRUE" = "#cb181d"),
                    labels = c("FALSE" = "p>=0.05", "TRUE" = "p<0.05"),
                    name = "Permutation") +
  labs(title = "Moran's I of axis signatures by sample",
       subtitle = "Positive I = spatial clustering of expression (effective n < nominal)",
       x = "Signature", y = "Moran's I") +
  sci_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1,
                                               size = 8))
sci_save(p.mi, "MoranI_Signature", width = 9, height = 7)

# =============================================================================
# PART 2: 空间块 bootstrap vs naive bootstrap 的 Cliff's delta CI
# 目的: 量化"忽略空间相关"对 CI 的压缩程度(反保守)
# 方法: 同比较(LM vs PT)下, naive 重抽样 spot 独立 vs 空间块重抽样
# 数据: 输出 SpatialBlockBootstrap_CI.csv
# =============================================================================

cat("========== PART 2: Spatial block vs naive bootstrap CI ==========\n")

set.seed(SEED)
boot.rows <- list()
for (sig in c("SPP1_mac", "sct_SPP1", "sct_PLAUR")) {
  bb <- BlockBootCliff(sig, "LM", "PT", master, BLOCK, N_BOOT)
  boot.rows[[length(boot.rows) + 1]] <- data.frame(
    signature = sig, method = "naive",
    estimate = median(bb$naive),
    ci_low = quantile(bb$naive, 0.025),
    ci_high = quantile(bb$naive, 0.975), stringsAsFactors = FALSE)
  boot.rows[[length(boot.rows) + 1]] <- data.frame(
    signature = sig, method = "spatial.block",
    estimate = median(bb$spatial),
    ci_low = quantile(bb$spatial, 0.025),
    ci_high = quantile(bb$spatial, 0.975), stringsAsFactors = FALSE)
}
boot.tab <- do.call(rbind, boot.rows)
write.csv(boot.tab, file.path(TAB_DIR, "SpatialBlockBootstrap_CI.csv"),
          row.names = FALSE)
print(boot.tab)

p.boot <- ggplot(boot.tab, aes(x = signature, y = estimate, colour = method,
                               ymin = ci_low, ymax = ci_high)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_pointrange(position = position_dodge(width = 0.4),
                  linewidth = 0.7, size = 1.2) +
  scale_colour_manual(values = c("naive" = "#999999",
                                "spatial.block" = "#007AC2"),
                      labels = c("naive" = "naive (iid spots)",
                                "spatial.block" = "spatial block"),
                      name = "Bootstrap") +
  labs(title = "Cliff's delta CI: naive vs spatial-block bootstrap",
       subtitle = "Wider CI under spatial blocking = honest effective sample size",
       x = "Signature", y = "Cliff's delta (LM vs PT)") +
  sci_theme
sci_save(p.boot, "SpatialBlockBootstrap_CI", width = 8, height = 6)

# =============================================================================
# PART 3: 既有 p 值家族 BH(FDR) 校正
# 目的: 统一 09/10/11/12 多比较的误报控制, 不修改原脚本产物
# 方法: 逐家族读 p 值, p.adjust(method="BH"), 输出合并表
# 数据: 输出 FDR_adjusted_summary.csv
# =============================================================================

cat("========== PART 3: Family-wise FDR (BH) correction ==========\n")

read.p <- function(path, pcol, label.col, script, family) {
  df <- read.csv(path)
  data.frame(script = script, family = family,
             label = as.character(df[[label.col]]),
             raw_p = df[[pcol]], stringsAsFactors = FALSE)
}
fdr.list <- list(
  read.p(file.path(TAB09, "CompartmentComparisons.csv"), "p_value",
         "signature", "09", "compartment"),
  read.p(file.path(TAB10, "MyeloidCoexpression.csv"), "p", "sample",
         "10", "myeloid_coexpr"),
  read.p(file.path(TAB10, "DoubleHighBoundary.csv"), "W_p", "sample",
         "10", "double_high_boundary"),
  read.p(file.path(TAB10, "NearestNeighbourPermutation.csv"), "p_two",
         "sample", "10", "nn_permutation"),
  read.p(file.path(TAB11, "CellChatPLAUR.csv"), "pval", "interaction",
         "11", "cellchat"),
  read.p(file.path(TAB12, "MediationResults.csv"), "p_value", "model",
         "12", "mediation"),
  read.p(file.path(TAB12, "MediationSensitivity.csv"), "p_value", "model",
         "12", "mediation_sensitivity"),
  read.p(file.path(TAB12, "PerSampleACME.csv"), "p_value", "sample",
         "12", "persample_acme")
)
fdr.tab <- do.call(rbind, fdr.list)
fdr.tab$fdr_p <- NA_real_
for (fam in unique(fdr.tab$family)) {
  m <- fdr.tab$family == fam
  fdr.tab$fdr_p[m] <- p.adjust(fdr.tab$raw_p[m], method = "BH")
}
fdr.tab$sig <- ifelse(fdr.tab$fdr_p < 0.05, "sig", "ns")
write.csv(fdr.tab, file.path(TAB_DIR, "FDR_adjusted_summary.csv"),
          row.names = FALSE)
cat("FDR 校正后各家族显著数 / 总数:\n")
print(aggregate(sig ~ family, data = fdr.tab, FUN = function(v) sum(v == "sig")))
cat("原始 p 与 FDR p 前 6 行:\n")
print(head(fdr.tab[, c("script", "family", "label", "raw_p", "fdr_p", "sig")],
           6))

cat("14_spatial_autocorr_fdr.R done: MoranI_Signature.pdf + ",
    "SpatialBlockBootstrap_CI.pdf + MoranI_by_signature_sample.csv + ",
    "SpatialBlockBootstrap_CI.csv + FDR_adjusted_summary.csv\n", sep = "")
