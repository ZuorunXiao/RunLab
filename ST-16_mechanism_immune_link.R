# =============================================================================
# 16_mechanism_immune_link.R
# 功能: 机制与免疫串联层——肿瘤源+髓系源双发送者联合中介; 距离分层 ACME;
#       SPP1+ TAM 与 CD8/Treg 的空间免疫排斥检验
# 定位: n=2(P1/P2 四切片)下把 12 的中介升级为双源、并串联免疫微环境, 不扩大样本
# 可视化: JointMediation.pdf (双源调整 ACME 对比)
#         DistanceStratifiedACME.pdf (按空间距离分层的 ACME)
#         ImmuneExclusion_radial.pdf (SPP1+ TAM 周边 CD8/Treg 径向轮廓)
#         ImmuneExclusion_test.pdf (SPP1+ TAM 邻域 vs 对照的免疫排斥检验)
# 数据: 输入 results/06.../table/signature_master.rds
#            results/12.../table/Triplets_tumor.csv / Triplets_myeloid.csv
#            results/03_immune_spatial/table/spot_celltype_annotation.csv
#       输出 results/16.../Figure/ + results/16.../table/
# 依据: ⑧ 在同一模型内同时纳入肿瘤 PLAU 与髓系 PLAU, 估计各自独立 ACME(控制对方),
#       直接回答"哪一源驱动 PLAUR->SPP1"; ⑲ 检验效应是否随配对距离衰减; ⑬ 把
#       SPP1+ TAM 与细胞毒性 CD8 / 调节性 Treg 的空间关系量化, 提示免疫抑制生态位
# 说明: 本脚本只读既有产物 + 重算, 不修改 00-13 任何文件
# =============================================================================

# --- 包加载(统计 -> 可视化) ----------------------------------------------
suppressMessages({
  library(mediation)   # 中介分析
  library(ggplot2)     # 可视化
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
    plot.subtitle    = element_text(size = 11, hjust = 0.5, face = "plain"),
    legend.title     = element_text(size = 11, face = "plain")
  )

sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf(file.path(FIG_DIR, paste0(filename, ".pdf")), width = width,
      height = height, family = "Helvetica")
  print(plot_obj)
  invisible(dev.off())
}

SigStar <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

# CliffDelta: 由 Mann-Whitney U 换算(与 14 同口径, 显式剔除非有限值)
CliffDelta <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  n1 <- length(a); n2 <- length(b)
  if (n1 < 1 || n2 < 1) return(NA_real_)
  ra <- rank(c(a, b))[seq_len(n1)]
  Ua <- sum(ra) - n1 * (n1 + 1) / 2
  2 * Ua / (n1 * n2) - 1
}

# --- 目录结构 -----------------------------------------------------------
BASE_DIR <- "X:/XX"
TAB06    <- file.path(BASE_DIR, "results/06_PLAUR_SPP1_signature/table")
TAB12    <- file.path(BASE_DIR, "results/12_PLAUR_SPP1_mediation/table")
TAB03    <- file.path(BASE_DIR, "results/03_immune_spatial/table")
FIG_DIR  <- file.path(BASE_DIR, "results/16_mechanism_immune_link/Figure")
TAB_DIR  <- file.path(BASE_DIR, "results/16_mechanism_immune_link/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 分析常量 -----------------------------------------------------------
SAMPLE_ID <- c("ST-P1_liver", "ST-P1_colon", "ST-P2_liver", "ST-P2_colon")
N_SIM     <- 1000   # mediation bootstrap 次数
N_PERM    <- 999    # ⑬ 置换次数
SEED      <- 42
SHELL_UM  <- 200    # ⑬ SPP1+ TAM 邻域半径(um)

# =============================================================================
# PART 1: (⑧) 双源联合中介(肿瘤 PLAU + 髓系 PLAU 同模型, 互相校正)
# 目的: 在同一模型内估计肿瘤源 / 髓系源各自的独立 ACME, 直接比较谁驱动 PLAUR->SPP1
# 方法: 按 macro_spot 合并肿瘤源与髓系源三元组, 得 plau_tumor + plau_mye;
#       lm(plaur_mac ~ 两源 + patient); lm(spp1_mac ~ 两源 + plaur_mac + patient);
#       mediate 分别以 plau_tumor / plau_mye 为 treat, 控制另一源
# 数据: 输出 JointMediation.csv
# =============================================================================

cat("========== PART 1: (8) Dual-source joint mediation ==========\n")

set.seed(SEED)
trip_tumor <- read.csv(file.path(TAB12, "Triplets_tumor.csv"))
trip_mye   <- read.csv(file.path(TAB12, "Triplets_myeloid.csv"))

joint <- merge(
  trip_tumor[, c("macro_spot", "plau_src", "plaur_mac", "spp1_mac",
                 "sample", "patient")],
  trip_mye[,  c("macro_spot", "plau_src")],
  by = "macro_spot", suffixes = c("_tumor", "_mye"))
cat("Joint macrophage spots (both sources within 300 um):", nrow(joint), "\n")

# 两源 PLAU 在同一巨噬 spot 上高度相关; 若强制放入单方程共线性会使 a 路径被
# 回归抑制翻负(见 12 前因: PLAU 与 PLAUR/SPP1 正相关, 偏出一方另一方翻负) ->
# 联立方程给出的是伪影而非生物学。故在"同一批共享 spot"上分别估计各源独立 ACME,
# 做公平比较, 而非单方程同放两源。记录两源相关以佐证该选择。
cor_both <- cor(joint$plau_src_tumor, joint$plau_src_mye,
                use = "pairwise.complete.obs")
cat("Correlation plau_tumor vs plau_mye (shared spots):", round(cor_both, 3), "\n")
write.csv(data.frame(metric = "cor_plau_tumor_plau_mye", value = round(cor_both, 4)),
          file.path(TAB_DIR, "JointSourceCorrelation.csv"), row.names = FALSE)

run_single_joint <- function(trip, label, contrast = "P75", with_patient = TRUE) {
  fits <- list(
    m = lm(plaur_mac ~ plau_src, data = trip),
    y = lm(spp1_mac ~ plau_src + plaur_mac, data = trip))
  nz <- trip$plau_src[trip$plau_src > 0]
  ctrl <- 0
  trt <- if (length(nz) >= 10) quantile(nz, 0.75) else max(trip$plau_src, na.rm = TRUE)
  if (is.na(trt) || trt <= ctrl) trt <- max(trip$plau_src, na.rm = TRUE)
  med <- tryCatch(mediate(fits$m, fits$y, treat = "plau_src",
                          mediator = "plaur_mac", boot = TRUE, sims = N_SIM,
                          control.value = ctrl, treat.value = trt),
                  error = function(e) stop("mediate() failed: ", conditionMessage(e)))
  stopifnot(!is.null(med), "mediate() returned NULL — 中介模型拟合失败")
  data.frame(model = label, effect = c("ACME", "ADE", "Total", "Prop.mediated"),
             estimate = c(med$d1, med$z1, med$tau.coef, med$n1),
             ci_low = c(med$d1.ci[1], med$z1.ci[1], med$tau.ci[1], med$n1.ci[1]),
             ci_high = c(med$d1.ci[2], med$z1.ci[2], med$tau.ci[2], med$n1.ci[2]),
             p_value = c(med$d1.p, med$z1.p, med$tau.p, med$n1.p),
             stringsAsFactors = FALSE)
}

# 同一批共享 spot, 各源独立 ACME(单源模型, 避免共线性)
jt <- joint[, c("plau_src_tumor", "plaur_mac", "spp1_mac", "sample", "patient")]
names(jt)[1] <- "plau_src"
jm <- joint[, c("plau_src_mye", "plaur_mac", "spp1_mac", "sample", "patient")]
names(jm)[1] <- "plau_src"
joint.tab <- rbind(
  run_single_joint(jt, "Tumor source (shared spots)", with_patient = TRUE),
  run_single_joint(jm, "Myeloid source (shared spots)", with_patient = TRUE))
joint.tab$star <- SigStar(joint.tab$p_value)
cat("--- Dual-source ACME compared on shared macrophage spots ---\n")
print(joint.tab, digits = 4)
write.csv(joint.tab, file.path(TAB_DIR, "JointMediation.csv"), row.names = FALSE)

joint.tab$effect <- factor(joint.tab$effect,
                           levels = c("ACME", "ADE", "Total", "Prop.mediated"))
joint.tab$model  <- factor(joint.tab$model,
   levels = c("Tumor source (shared spots)", "Myeloid source (shared spots)"))
p.joint <- ggplot(joint.tab[joint.tab$effect %in% c("ACME", "ADE"), ],
                  aes(x = effect, y = estimate, colour = model, group = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15,
                linewidth = 0.7, position = position_dodge(width = 0.4)) +
  geom_point(aes(shape = model), size = 4,
             position = position_dodge(width = 0.4)) +
  geom_text(aes(label = star), colour = "black", vjust = -0.9, size = 4,
            position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("Tumor source (shared spots)" = "#999999",
                                 "Myeloid source (shared spots)" = "#007AC2"),
                      name = "PLAU source") +
  scale_shape_manual(values = c("Tumor source (shared spots)" = 16,
                                "Myeloid source (shared spots)" = 17),
                     name = "PLAU source") +
  labs(title = "Dual-source ACME compared on shared macrophage spots",
       subtitle = "Single-source models on same spots (joint equation avoided: r_source correlation shown in CSV)",
       x = "Effect", y = "Estimate with 95% CI") +
  sci_theme + theme(legend.position = "bottom")
sci_save(p.joint, "JointMediation", width = 9, height = 5.5)

# =============================================================================
# PART 2: (⑲) 距离分层 ACME(髓系源, 检验效应是否随空间距离衰减)
# 目的: 把 12 的主结论按空间配对距离分层, 看中介效应是否在近距离更强
# 方法: 逐距离层(<=100, 100-200, 200-300 um)在 trip_mye 子集上跑 mediate
# 数据: 输出 DistanceStratifiedACME.csv
# =============================================================================

cat("========== PART 2: (19) Distance-stratified ACME ==========\n")

run_single <- function(trip, label, setting, contrast = "P75",
                       with_patient = FALSE) {
  fits <- list(
    m = lm(plaur_mac ~ plau_src, data = trip),
    y = lm(spp1_mac ~ plau_src + plaur_mac, data = trip))
  nz <- trip$plau_src[trip$plau_src > 0]
  ctrl <- 0
  trt <- if (length(nz) >= 10) quantile(nz, 0.75) else max(trip$plau_src, na.rm = TRUE)
  if (is.na(trt) || trt <= ctrl) trt <- max(trip$plau_src, na.rm = TRUE)
  med <- tryCatch(mediate(fits$m, fits$y, treat = "plau_src",
                          mediator = "plaur_mac", boot = TRUE, sims = N_SIM,
                          control.value = ctrl, treat.value = trt),
                  error = function(e) stop("mediate() failed: ", conditionMessage(e)))
  stopifnot(!is.null(med), "mediate() returned NULL — 中介模型拟合失败")
  data.frame(model = label, setting = setting, effect = "ACME",
             estimate = med$d1, ci_low = med$d1.ci[1], ci_high = med$d1.ci[2],
             p_value = med$d1.p, stringsAsFactors = FALSE)
}

layers <- list(c(0, 100), c(100, 200), c(200, 300))
dist.tab <- list()
for (lr in layers) {
  sub <- subset(trip_mye, dist_um > lr[1] & dist_um <= lr[2])
  if (nrow(sub) < 30) { cat("skip layer", lr[1], "-", lr[2], "n=", nrow(sub), "\n"); next }
  res <- run_single(sub, "Myeloid source",
                    setting = paste0("Dist ", lr[1], "-", lr[2], " um"))
  if (!is.null(res)) { res$n <- nrow(sub); dist.tab[[length(dist.tab) + 1]] <- res }
}
dist.tab <- do.call(rbind, dist.tab)
dist.tab$star <- SigStar(dist.tab$p_value)
dist.tab$setting <- factor(dist.tab$setting,
  levels = c("Dist 0-100 um", "Dist 100-200 um", "Dist 200-300 um"))
cat("--- Distance-stratified ACME (myeloid source) ---\n"); print(dist.tab, digits = 4)
write.csv(dist.tab, file.path(TAB_DIR, "DistanceStratifiedACME.csv"), row.names = FALSE)

p.dist <- ggplot(dist.tab, aes(x = setting, y = estimate, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15,
                linewidth = 0.7) +
  geom_point(size = 4, colour = "#007AC2") +
  geom_text(aes(label = star), colour = "black", vjust = -0.9, size = 4) +
  geom_text(aes(label = sprintf("n=%d", n), y = ci_low), vjust = 1.6, size = 3,
            colour = "black") +
  labs(title = "ACME of PLAU -> PLAUR -> SPP1 by spatial pair distance",
       subtitle = "Myeloid source; tests whether mediation decays with distance",
       x = "Pair distance layer", y = "ACME estimate with 95% CI") +
  sci_theme + theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 10))
sci_save(p.dist, "DistanceStratifiedACME", width = 8, height = 5)

# =============================================================================
# PART 3: (⑬) SPP1+ TAM 与 CD8/Treg 的空间免疫排斥
# 目的: 检验 SPP1+ TAM 生态位是否排斥细胞毒性 CD8(或富集 Treg)——免疫抑制生态位
# 方法: (a) 径向轮廓: SPP1+ TAM 焦点周边 CD8_T_cells / FOXP3_Treg 均值随距离变化;
#       (b) 邻域检验: SPP1+ TAM 焦点 r<SHELL_UM 邻域 vs 等数对照巨噬焦点邻域的
#       免疫评分差异, Wilcoxon + 置换伪 p
# 数据: 输出 ImmuneExclusion_radial.csv / ImmuneExclusion_test.csv
# =============================================================================

cat("========== PART 3: (13) Immune exclusion around SPP1+ TAM ==========\n")

master <- readRDS(file.path(TAB06, "signature_master.rds"))
master$patient <- substr(master$sample, 4, 5)
imm <- read.csv(file.path(TAB03, "spot_celltype_annotation.csv"))
imm <- imm[, c("sample", "spot", "CD8_T_cells", "FOXP3_Treg")]
mg <- merge(master, imm, by = c("sample", "spot"))
mg <- mg[!is.na(mg$x) & !is.na(mg$y), ]

rad.rows <- list(); test.rows <- list()
test.cells <- c("CD8_T_cells", "FOXP3_Treg")
for (s in SAMPLE_ID) {
  d <- mg[mg$sample == s, ]
  xy <- as.matrix(d[, c("x", "y")])
  DD <- as.matrix(dist(xy))   # 像素距离, 仅需相对远近, 不标定 um
  # SPP1+ TAM 焦点(巨噬且 sct_SPP1 > 样本中位)
  is_focus <- d$is_macro & !is.na(d$sct_SPP1) & (d$sct_SPP1 > median(d$sct_SPP1, na.rm = TRUE))
  foci <- which(is_focus)
  # 对照巨噬焦点: 非 SPP1+ 的巨噬 spot(等数随机抽)
  ctrl_pool <- which(d$is_macro & !is.na(d$sct_SPP1) & !is_focus)
  if (length(foci) < 5 || length(ctrl_pool) < 5) {
    cat(s, ": skipped (foci=", length(foci), " ctrl=", length(ctrl_pool), ")\n"); next
  }
  set.seed(SEED)
  ctrl <- sample(ctrl_pool, size = length(foci), replace = FALSE)
  # (a) 径向轮廓: 焦点周边 shell 内免疫评分均值
  for (fi in foci) {
    oi <- (1:nrow(d))[-fi]
    di <- DD[fi, oi]
    for (cell in test.cells) {
      # 用 shell 100 um 像素尺度: 最近邻中值近似 spot 间距, 分 5 壳
      pix_pitch <- median(di[di > 0])
      sh <- ceiling(di / (2 * pix_pitch))   # 壳序号(1=最近)
      sh[sh > 5] <- 5
      for (k in 1:5) {
        m <- !is.na(d[[cell]][oi]) & sh == k
        if (sum(m) == 0) next
        rad.rows[[length(rad.rows) + 1]] <- data.frame(
          sample = s, immune_cell = cell, shell = k,
          mean_score = mean(d[[cell]][oi][m]), n_spots = sum(m),
          stringsAsFactors = FALSE)
      }
    }
  }
  # (b) 邻域检验: r<SHELL_UM 邻域免疫评分(SPP1+ TAM vs 对照)
  # um.px = 像素/um(最近邻中值像素距 / 100um 标称); SHELL_UM 换算成像素阈值
  um.px <- median(apply(DD, 1, function(v) min(v[v > 0]))) / 100
  thr_pix <- SHELL_UM * um.px
  # 预计算每 spot 邻域索引(一次), 置换内仅重排邻域得分列表, 避免重复建距离子集
  nb_idx <- vector("list", nrow(d))
  for (i in 1:nrow(d)) {
    oi <- (1:nrow(d))[-i]
    nb_idx[[i]] <- oi[DD[i, oi] <= thr_pix]
  }
  for (cell in test.cells) {
    nb_score <- lapply(nb_idx, function(nb) {
      v <- d[[cell]][nb]; v[!is.na(v)]
    })
    near_focus <- unlist(nb_score[foci])
    near_ctrl  <- unlist(nb_score[ctrl])
    cd <- CliffDelta(near_ctrl, near_focus)  # >0: 对照邻域免疫更高 = SPP1+ TAM 排斥
    # 置换伪 p: 打乱 focus/ctrl 标签(仅重排邻域得分列表)
    all_idx <- c(foci, ctrl); grp <- c(rep(1, length(foci)), rep(0, length(ctrl)))
    obs <- mean(near_ctrl, na.rm = TRUE) - mean(near_focus, na.rm = TRUE)
    p <- 0
    for (b in 1:N_PERM) {
      gp <- sample(grp)
      a  <- unlist(nb_score[all_idx[gp == 1]])
      c2 <- unlist(nb_score[all_idx[gp == 0]])
      if (abs(mean(c2, na.rm = TRUE) - mean(a, na.rm = TRUE)) >= abs(obs)) p <- p + 1
    }
    test.rows[[length(test.rows) + 1]] <- data.frame(
      sample = s, immune_cell = cell,
      cliff_delta = round(cd, 4),
      p_perm = (p + 1) / (N_PERM + 1),
      mean_near_SPP1TAM = round(mean(near_focus, na.rm = TRUE), 4),
      mean_near_control = round(mean(near_ctrl, na.rm = TRUE), 4),
      stringsAsFactors = FALSE)
  }
  cat(sprintf("%s: foci=%d  ctrl=%d\n", s, length(foci), length(ctrl)))
}
rad.tab <- do.call(rbind, rad.rows)
test.tab <- do.call(rbind, test.rows)
write.csv(rad.tab, file.path(TAB_DIR, "ImmuneExclusion_radial.csv"), row.names = FALSE)
write.csv(test.tab, file.path(TAB_DIR, "ImmuneExclusion_test.csv"), row.names = FALSE)
print(test.tab)

rad.tab$sample <- factor(rad.tab$sample, levels = SAMPLE_ID)
rad.tab$immune_cell <- factor(rad.tab$immune_cell,
  levels = c("CD8_T_cells", "FOXP3_Treg"),
  labels = c("CD8 T cells", "FOXP3 Treg"))
p.rad <- ggplot(rad.tab, aes(x = shell, y = mean_score, colour = sample, group = sample)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_colour_manual(values = c("ST-P1_liver" = "#cb181d", "ST-P1_colon" = "#4292c6",
    "ST-P2_liver" = "#fb6a4a", "ST-P2_colon" = "#6baed6"), name = "Sample") +
  facet_wrap(~ immune_cell, ncol = 2) +
  labs(title = "Immune cell score around SPP1+ TAM foci (radial shells)",
       subtitle = "Shell 1 = nearest neighbours; depletion near focus suggests exclusion",
       x = "Neighbour shell (1=nearest)", y = "Mean immune-cell score") +
  sci_theme + theme(legend.position = "bottom")
sci_save(p.rad, "ImmuneExclusion_radial", width = 9, height = 5)

test.tab$sample <- factor(test.tab$sample, levels = SAMPLE_ID)
test.tab$immune_cell <- factor(test.tab$immune_cell,
  levels = c("CD8_T_cells", "FOXP3_Treg"),
  labels = c("CD8 T cells", "FOXP3 Treg"))
p.test <- ggplot(test.tab, aes(x = sample, y = cliff_delta, fill = immune_cell)) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f\np=%.3g", cliff_delta, p_perm)),
            position = position_dodge(width = 0.7), vjust = ifelse(test.tab$cliff_delta > 0, -0.3, 1.3),
            size = 3, colour = "black") +
  scale_fill_manual(values = c("CD8 T cells" = "#cb181d", "FOXP3 Treg" = "#7a0177"),
                    name = "Immune cell") +
  labs(title = "Immune exclusion test: SPP1+ TAM neighbourhood vs control",
       subtitle = "Cliff's delta > 0 = control neighbourhood has higher immune score (SPP1+ TAM depletes/excludes)",
       x = "Sample", y = "Cliff's delta (control - SPP1+ TAM)") +
  sci_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
                    legend.position = "bottom")
sci_save(p.test, "ImmuneExclusion_test", width = 9, height = 5)

cat("16_mechanism_immune_link.R done: JointMediation.pdf + DistanceStratifiedACME.pdf + ",
    "ImmuneExclusion_radial.pdf + ImmuneExclusion_test.pdf + JointMediation.csv + ",
    "DistanceStratifiedACME.csv + ImmuneExclusion_radial.csv + ImmuneExclusion_test.csv\n",
    sep = "")
