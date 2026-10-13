# =============================================================================
# 15_spatial_structure_deepen.R
# 功能: 空间结构深化层——(⑪) SPP1+ TAM 生态位的径向衰减; (⑯) PLAU/PLAUR/SPP1
#       三节点共定位; (⑰) 肿瘤侵袭前沿的连续梯度(loess)
# 定位: n=2(P1/P2 四切片)下补充空间维度的严谨证据, 不扩大样本、不改写既有脚本
# 可视化: RadialDecay_sctSPP1.pdf / RadialDecay_PLAURmac.pdf (生态位径向衰减, 双面板)
#         TripleColocalization_compartment.pdf / _autocorr.pdf / _triplets.pdf
#         (三节点共定位: 区室占比 + 自相关 + 三元组)
#         TumorFrontGradient.pdf (前沿连续梯度: Gaussian 核平滑 + 分箱均值)
# 数据: 输入 results/06.../table/signature_master.rds
#            results/12.../table/Triplets_tumor.csv / Triplets_myeloid.csv
#       输出 results/15.../Figure/ + results/15.../table/
# 依据: ⑪ 量化 SPP1+ TAM 聚集体空间尺度(随距离衰减); ⑯ 把 10 的双高共位升级为
#       三基因共定位并量化其空间聚集; ⑰ 用"到肿瘤前沿距离"作连续变量, loess
#       刻画 SPP1+ TAM 沿侵袭前沿的富集梯度, 避免用离散区室掩盖连续结构
# 说明: 本脚本只读既有产物 + 重算, 不修改 00-13 任何文件
# =============================================================================

# --- 包加载(可视化) ------------------------------------------------------
suppressMessages({
  library(ggplot2)  # 可视化
})

# Gaussian 核平滑(替代 loess: 本机 R 构建 loess/predict.loess 返回 NA, 不可靠)
# x,y: 观测向量; xs: 网格; bw: 带宽(um)
gauss_smooth <- function(x, y, xs, bw) {
  out <- numeric(length(xs))
  for (j in seq_along(xs)) {
    w <- exp(-((x - xs[j])^2) / (2 * bw^2))
    sw <- sum(w)
    out[j] <- if (sw <= 0) NA_real_ else sum(w * y) / sw
  }
  out
}

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
    plot.subtitle    = element_text(size = 11, hjust = 0.5, face = "plain"),
    legend.title     = element_text(size = 11, face = "plain")
  )

sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf(file.path(FIG_DIR, paste0(filename, ".pdf")), width = width,
      height = height, family = "Helvetica")
  print(plot_obj)
  invisible(dev.off())
}

# UmPerPixel: 以标称 spot 间距 PITCH_UM=100 绑定最近邻中值, 把像素距离标定成 um
UmPerPixel <- function(xy) {
  d <- as.matrix(dist(xy)); diag(d) <- Inf
  PITCH_UM / median(apply(d, 1, min))
}

# moranI: k 近邻 Moran's I(内联, 与 12 同口径)
moranI <- function(xy, val, k = 6) {
  d <- as.matrix(dist(xy)); diag(d) <- Inf
  n <- nrow(xy); W <- matrix(0, n, n)
  for (i in 1:n) W[i, order(d[i, ])[1:k]] <- 1
  W <- W / rowSums(W)
  xbar <- mean(val, na.rm = TRUE)
  num <- sum(W * outer(val - xbar, val - xbar), na.rm = TRUE)
  den <- sum((val - xbar)^2, na.rm = TRUE)
  (n / sum(W)) * num / den
}

# --- 目录结构 -----------------------------------------------------------
BASE_DIR <- "X:/XX"
TAB06    <- file.path(BASE_DIR, "results/06_PLAUR_SPP1_signature/table")
TAB12    <- file.path(BASE_DIR, "results/12_PLAUR_SPP1_mediation/table")
FIG_DIR  <- file.path(BASE_DIR, "results/15_spatial_structure_deepen/Figure")
TAB_DIR  <- file.path(BASE_DIR, "results/15_spatial_structure_deepen/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 分析常量 -----------------------------------------------------------
SAMPLE_ID  <- c("ST-P1_liver", "ST-P1_colon", "ST-P2_liver", "ST-P2_colon")
PITCH_UM   <- 100    # Visium 相邻 spot 中心距标称值(um)
R_MAX      <- 600    # ⑪ 径向衰减最大半径(um)
BIN_W      <- 100    # ⑪ 径向分箱宽(um)
RAD_BINS   <- c(seq(0, R_MAX, by = BIN_W), Inf)
K_FRONT    <- 6      # ⑰ 肿瘤前沿判定的近邻数
WIN_FRONT  <- 1000   # ⑰ 到前沿距离上限(um)
SEED       <- 42

CMP_ORDER <- c("Adjacent colon", "Adjacent liver", "PT", "LM")
CMP_FILL  <- c("Adjacent colon" = "#9ecae1", "Adjacent liver" = "#4292c6",
               "PT" = "#fb6a4a", "LM" = "#cb181d")
SAMPLE_COL <- c("ST-P1_liver" = "#cb181d", "ST-P1_colon" = "#4292c6",
                "ST-P2_liver" = "#fb6a4a", "ST-P2_colon" = "#6baed6")
SAMPLE_SHAPE <- c("ST-P1_liver" = 16, "ST-P1_colon" = 17,
                  "ST-P2_liver" = 15, "ST-P2_colon" = 18)

# =============================================================================
# PART 1: (⑪) SPP1+ TAM 生态位径向衰减
# 目的: 量化 SPP1+ TAM 聚集体的空间尺度——以焦点为中心, 向外壳层表达是否衰减
# 方法: 逐样本取 SPP1+ TAM 焦点(巨噬且 sct_SPP1 > 样本中位), 计算其余巨噬 spot
#       到焦点的像素距离并标定 um, 分箱聚合 sct_SPP1 / SPP1_mac 均值
# 数据: 输出 RadialDecay.csv
# =============================================================================

cat("========== PART 1: (11) Radial decay of SPP1+ TAM niche ==========\n")

set.seed(SEED)
master <- readRDS(file.path(TAB06, "signature_master.rds"))
master$patient <- substr(master$sample, 4, 5)

decay.rows <- list()
for (s in SAMPLE_ID) {
  d <- master[master$sample == s & !is.na(master$x) & !is.na(master$y), ]
  dm <- d[d$is_macro & !is.na(d$sct_SPP1),
          c("x", "y", "sct_SPP1", "SPP1_mac", "PLAUR_mac", "compartment")]
  if (nrow(dm) < 20) { cat(s, ": skipped (n_macro=", nrow(dm), ")\n"); next }
  xy.m <- as.matrix(dm[, c("x", "y")])
  DD.m <- as.matrix(dist(xy.m)) * UmPerPixel(xy.m)   # um
  n.m <- nrow(dm)
  foci <- which(dm$sct_SPP1 > median(dm$sct_SPP1, na.rm = TRUE))
  nb <- length(RAD_BINS) - 1
  spp1_sum <- numeric(nb); spp1_cnt <- integer(nb)
  mac_sum  <- numeric(nb); mac_cnt  <- integer(nb)
  for (fi in foci) {
    oi <- (1:n.m)[-fi]
    di <- DD.m[fi, oi]
    bi <- cut(di, RAD_BINS, include.lowest = TRUE, labels = FALSE)
    spp1_sum[bi] <- spp1_sum[bi] + dm$sct_SPP1[oi]
    spp1_cnt[bi] <- spp1_cnt[bi] + 1
    mm <- !is.na(dm$PLAUR_mac[oi])
    mac_sum[bi[mm]] <- mac_sum[bi[mm]] + dm$PLAUR_mac[oi][mm]
    mac_cnt[bi[mm]] <- mac_cnt[bi[mm]] + 1
  }
  for (b in 1:nb) {
    decay.rows[[length(decay.rows) + 1]] <- data.frame(
      sample = s,
      radius_bin_um = ifelse(b < nb, paste0((b-1)*BIN_W, "-", b*BIN_W),
                             paste0((b-1)*BIN_W, "+")),
      radius_center_um = (b - 0.5) * BIN_W,
      mean_sct_SPP1 = ifelse(spp1_cnt[b] > 0, spp1_sum[b] / spp1_cnt[b], NA_real_),
      mean_PLAUR_mac = ifelse(mac_cnt[b] > 0, mac_sum[b] / mac_cnt[b], NA_real_),
      n_pairs = spp1_cnt[b], stringsAsFactors = FALSE)
  }
  cat(sprintf("%s: foci=%d  macrophage spots=%d\n", s, length(foci), n.m))
}
decay.tab <- do.call(rbind, decay.rows)
write.csv(decay.tab, file.path(TAB_DIR, "RadialDecay.csv"), row.names = FALSE)
print(decay.tab)

decay.tab$sample <- factor(decay.tab$sample, levels = SAMPLE_ID)
p.decay <- ggplot(decay.tab, aes(x = radius_center_um, y = mean_sct_SPP1,
                                 colour = sample, group = sample)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3) +
  scale_colour_manual(values = SAMPLE_COL, name = "Sample") +
  labs(title = "Radial decay of SPP1 transcript around SPP1+ TAM foci",
       subtitle = "Mean sct_SPP1 of macrophages in distance shells from each focus (um-scaled)",
       x = "Distance from focus (um)", y = "Mean sct_SPP1 (macrophages)") +
  sci_theme + theme(legend.position = "bottom")
sci_save(p.decay, "RadialDecay_sctSPP1", width = 8, height = 5)

p.decay2 <- ggplot(decay.tab, aes(x = radius_center_um, y = mean_PLAUR_mac,
                                  colour = sample, group = sample)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3) +
  scale_colour_manual(values = SAMPLE_COL, name = "Sample") +
  labs(title = "Radial decay of PLAUR macrophage signature around SPP1+ TAM foci",
       subtitle = "Mean PLAUR_mac of macrophages in distance shells from each focus (um-scaled)",
       x = "Distance from focus (um)", y = "Mean PLAUR_mac (macrophages)") +
  sci_theme + theme(legend.position = "bottom")
sci_save(p.decay2, "RadialDecay_PLAURmac", width = 8, height = 5)

# =============================================================================
# PART 2: (⑯) PLAU/PLAUR/SPP1 三节点共定位
# 目的: 把 10 的双高共位升级为三基因共定位, 并量化其空间聚集程度
# 方法: (a) spot 级三阳性比例(逐样本中位阈值)按区室; (b) 三阳性指示变量的 k-NN
#       Moran's I; (c) 三元组表三阳性共表达(plau_src & plaur_mac & spp1_mac >0)
# 数据: 输出 TriplePositive_by_compartment.csv / TriplePositive_autocorr.csv /
#              TripleNode_triplets.csv
# =============================================================================

cat("========== PART 2: (16) Three-node colocalization ==========\n")

# (a) 区室三阳性比例
tp.comp <- list()
for (s in SAMPLE_ID) {
  d <- master[master$sample == s & !is.na(master$x) & !is.na(master$y), ]
  tp <- (d$sct_PLAU > median(d$sct_PLAU, na.rm = TRUE)) &
        (d$sct_PLAUR > median(d$sct_PLAUR, na.rm = TRUE)) &
        (d$sct_SPP1 > median(d$sct_SPP1, na.rm = TRUE))
  for (c in CMP_ORDER) {
    m <- d$compartment == c & !is.na(d$compartment)
    if (sum(m) == 0) next
    tp.comp[[length(tp.comp) + 1]] <- data.frame(
      sample = s, compartment = c,
      pct_triple_positive = round(100 * mean(tp[m]), 1),
      n_spots = sum(m), stringsAsFactors = FALSE)
  }
}
tp.comp.tab <- do.call(rbind, tp.comp)
tp.comp.tab$compartment <- factor(tp.comp.tab$compartment, levels = CMP_ORDER)
tp.comp.tab$sample <- factor(tp.comp.tab$sample, levels = SAMPLE_ID)
write.csv(tp.comp.tab, file.path(TAB_DIR, "TriplePositive_by_compartment.csv"),
          row.names = FALSE)
print(tp.comp.tab)

# (b) 三阳性指示变量空间自相关
tp.aut <- list()
for (s in SAMPLE_ID) {
  d <- master[master$sample == s & !is.na(master$x) & !is.na(master$y), ]
  tp <- as.numeric((d$sct_PLAU > median(d$sct_PLAU, na.rm = TRUE)) &
                   (d$sct_PLAUR > median(d$sct_PLAUR, na.rm = TRUE)) &
                   (d$sct_SPP1 > median(d$sct_SPP1, na.rm = TRUE)))
  if (sum(tp) < 5) { cat(s, ": skipped (triple-pos n<5)\n"); next }
  tp.aut[[s]] <- data.frame(sample = s,
    morans_I = round(moranI(as.matrix(d[, c("x", "y")]), tp, k = 6), 3),
    n_triple = sum(tp), stringsAsFactors = FALSE)
}
tp.aut.tab <- do.call(rbind, tp.aut)
tp.aut.tab$sample <- factor(tp.aut.tab$sample, levels = SAMPLE_ID)
write.csv(tp.aut.tab, file.path(TAB_DIR, "TriplePositive_autocorr.csv"),
          row.names = FALSE)
print(tp.aut.tab)

# (c) 三元组三阳性共表达
triple_node <- function(csv_path, src_label) {
  tr <- read.csv(csv_path)
  tr$triple <- (tr$plau_src > 0) & (tr$plaur_mac > 0) & (tr$spp1_mac > 0)
  sp <- sort(unique(tr$sample))
  res <- lapply(sp, function(s) {
    x <- tr$triple[tr$sample == s]
    data.frame(sample = s,
               pct_triple = round(100 * mean(x), 1),
               n_triplets = length(x),
               source = src_label, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, res)
  res$sample <- factor(res$sample, levels = SAMPLE_ID)
  res
}
tn.tumor <- triple_node(file.path(TAB12, "Triplets_tumor.csv"), "Tumor source")
tn.mye   <- triple_node(file.path(TAB12, "Triplets_myeloid.csv"), "Myeloid source")
tn.tab   <- rbind(tn.tumor, tn.mye)
write.csv(tn.tab, file.path(TAB_DIR, "TripleNode_triplets.csv"), row.names = FALSE)
print(tn.tab)

# 三节点共定位图: 三面板(区室占比 / 自相关 / 三元组)
p.tp.comp <- ggplot(tp.comp.tab, aes(x = compartment, y = pct_triple_positive,
                                     fill = compartment)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", pct_triple_positive)), vjust = -0.3,
            size = 3, colour = "black") +
  scale_fill_manual(values = CMP_FILL, name = "Compartment") +
  facet_wrap(~ sample, ncol = 4) +
  scale_y_continuous(limits = c(0, max(tp.comp.tab$pct_triple_positive) * 1.2)) +
  labs(title = "Triple-positive (PLAU/PLAUR/SPP1) spot fraction by compartment",
       x = "Compartment", y = "Percent triple-positive (%)") +
  sci_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 8),
                    legend.position = "none", panel.spacing = unit(0.6, "lines"))

p.tp.aut <- ggplot(tp.aut.tab, aes(x = sample, y = morans_I, fill = sample)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.3f\n(n=%d)", morans_I, n_triple)), vjust = -0.3,
            size = 3, colour = "black") +
  scale_fill_manual(values = SAMPLE_COL, name = "Sample") +
  scale_y_continuous(limits = c(0, max(tp.aut.tab$morans_I) * 1.25)) +
  labs(title = "Spatial autocorrelation of triple-positive spots (k=6)",
       x = "Sample", y = "Moran's I") +
  sci_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
                    legend.position = "none")

p.tn <- ggplot(tn.tab, aes(x = sample, y = pct_triple, fill = source)) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_text(aes(label = sprintf("%.1f%%", pct_triple)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3,
            colour = "black") +
  scale_fill_manual(values = c("Tumor source" = "#999999",
                               "Myeloid source" = "#007AC2"), name = "Source") +
  scale_y_continuous(limits = c(0, max(tn.tab$pct_triple) * 1.2)) +
  labs(title = "Triple co-expression within spatial triplets",
       subtitle = "Triplets where plau_src & plaur_mac & spp1_mac all > 0",
       x = "Sample", y = "Percent triple co-expressed (%)") +
  sci_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
                    legend.position = "bottom")

# 三面板分别存为独立 PDF(与 12 的逐图惯例一致, 避免额外依赖)
sci_save(p.tp.comp, "TripleColocalization_compartment", width = 14, height = 4)
sci_save(p.tp.aut, "TripleColocalization_autocorr", width = 8, height = 5)
sci_save(p.tn, "TripleColocalization_triplets", width = 8, height = 5)

# =============================================================================
# PART 3: (⑰) 肿瘤侵袭前沿连续梯度
# 目的: 用"到肿瘤前沿距离"作连续变量, loess 刻画 SPP1+ TAM 沿前沿的富集梯度
# 方法: 逐样本定肿瘤前沿(肿瘤 spot 且 k 近邻含非肿瘤), 算各髓系 spot 到前沿 um
#       距离, 分箱均值 + Gaussian 核平滑 SPP1_mac ~ 距离
# 数据: 输出 TumorFront_gradient.csv
# =============================================================================

cat("========== PART 3: (17) Tumor-front continuous gradient ==========\n")

front.rows <- list()
loess.rows <- list()
for (s in SAMPLE_ID) {
  d <- master[master$sample == s & !is.na(master$x) & !is.na(master$y), ]
  xy <- as.matrix(d[, c("x", "y")])
  DD <- as.matrix(dist(xy))
  nn <- t(apply(DD, 1, function(v) order(v)[2:(K_FRONT + 1)]))  # 排除自身
  is_front <- logical(nrow(d))
  for (i in 1:nrow(d)) {
    if (d$is_tumor[i] && any(!d$is_tumor[nn[i, ]])) is_front[i] <- TRUE
  }
  front_idx <- which(is_front)
  if (length(front_idx) == 0) { cat(s, ": no front detected\n"); next }
  um.px <- UmPerPixel(xy)
  dist_front <- apply(DD[, front_idx, drop = FALSE], 1, min) * um.px
  dm <- d[d$is_macro & !is.na(d$sct_SPP1), ]
  # dist_front 按 d 的行位命名, 用 dm 在 d 中的行位对齐(非 spot ID)
  dm_pos <- which(d$is_macro & !is.na(d$sct_SPP1))
  dfront_m <- dist_front[dm_pos]
  dfront_m[dfront_m > WIN_FRONT] <- WIN_FRONT
  # 分箱均值
  bins <- cut(dfront_m, seq(0, WIN_FRONT, by = 100), include.lowest = TRUE,
              labels = FALSE)
  for (b in 1:(WIN_FRONT / 100)) {
    m <- bins == b & !is.na(dm$SPP1_mac)
    if (sum(m) == 0) next
    front.rows[[length(front.rows) + 1]] <- data.frame(
      sample = s,
      dist_bin_um = paste0((b-1)*100, "-", b*100),
      dist_center_um = (b - 0.5) * 100,
      mean_SPP1_mac = mean(dm$SPP1_mac[m], na.rm = TRUE),
      n_macro = sum(m), stringsAsFactors = FALSE)
  }
  # Gaussian 核平滑(仅用有值点; 替代 loess, 本机构建失效)
  m2 <- !is.na(dm$SPP1_mac)
  if (sum(m2) > 10) {
    xs <- seq(0, WIN_FRONT, by = 25)
    bw <- 75  # Gaussian 核带宽(um)
    pr <- gauss_smooth(dfront_m[m2], dm$SPP1_mac[m2], xs, bw)
    for (k in seq_along(xs)) {
      loess.rows[[length(loess.rows) + 1]] <- data.frame(
        sample = s, dist_center_um = xs[k], loess_SPP1_mac = pr[k],
        stringsAsFactors = FALSE)
    }
  }
  cat(sprintf("%s: front spots=%d  myeloid spots=%d\n", s, length(front_idx),
              nrow(dm)))
}
front.tab <- do.call(rbind, front.rows)
loess.tab <- do.call(rbind, loess.rows)
write.csv(front.tab, file.path(TAB_DIR, "TumorFront_gradient.csv"), row.names = FALSE)
print(front.tab)

front.tab$sample <- factor(front.tab$sample, levels = SAMPLE_ID)
loess.tab$sample <- factor(loess.tab$sample, levels = SAMPLE_ID)

p.front <- ggplot() +
  geom_point(data = front.tab, aes(x = dist_center_um, y = mean_SPP1_mac,
                                   colour = sample, shape = sample), size = 2) +
  geom_line(data = loess.tab, aes(x = dist_center_um, y = loess_SPP1_mac,
                                  colour = sample, group = sample), linewidth = 0.9) +
  scale_colour_manual(values = SAMPLE_COL, name = "Sample") +
  scale_shape_manual(values = SAMPLE_SHAPE, name = "Sample") +
  labs(title = "SPP1+ TAM enrichment gradient from tumor invasive front",
       subtitle = "Mean SPP1_mac of myeloid spots vs distance to tumor front (points: 100 um bin means; line: Gaussian kernel, bw=75 um)",
       x = "Distance to tumor front (um)", y = "Mean SPP1_mac") +
  sci_theme + theme(legend.position = "bottom")
sci_save(p.front, "TumorFrontGradient", width = 9, height = 6)

cat("15_spatial_structure_deepen.R done: RadialDecay_sctSPP1.pdf + ",
    "RadialDecay_PLAURmac.pdf + TripleColocalization_compartment.pdf + ",
    "TripleColocalization_autocorr.pdf + TripleColocalization_triplets.pdf + ",
    "TumorFrontGradient.pdf + RadialDecay.csv + TriplePositive_by_compartment.csv + ",
    "TriplePositive_autocorr.csv + TripleNode_triplets.csv + TumorFront_gradient.csv\n",
    sep = "")
