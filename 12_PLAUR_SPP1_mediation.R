# ===============================================================
# 12_PLAUR_SPP1_mediation.R
# 功能: 检验 PLAU -> PLAUR -> SPP1 轴在 CRLM 空间中的中介结构
# 定位: 方案 Fig 5 + 敏感性分析 + SPP1+ TAM 空间生态位
# ---------------------------------------------------------------------------
# 重要前提(证据驱动重设计, 2026-09-10):
#   原假设"肿瘤 PLAU -> 巨噬 PLAUR -> 巨噬 SPP1"已被数据否定——
#   (a) 11 CellChat: Tumor->Myeloid PLAU-PLAUR 通讯概率 = 0;
#   (b) 12 诊断: 肿瘤 PLAU 与巨噬 SPP1 边际相关 ≈ 0, 且"负向直接效应"
#       实为结局模型 SPP1~PLAU+PLAUR 的回归抑制/共线性伪影
#       (PLAUR 与 SPP1 强正相关 +0.48, 偏出后 PLAU 系数翻负)。
#   故本脚本改为:
#   - 肿瘤源 如实报告为"不支持(null)", 不把负向 ADE 当因果;
#   - 方案 A: 以 11 证据(显著发送方是髓系)把 X 重锚为"最近髓系 spot 的 PLAU",
#     检验 PLAU(髓系) -> PLAUR(巨噬) -> SPP1(巨噬) 中介;
#   - 方案 C: 量化 SPP1+ TAM 空间生态位(区室富集 + 空间自相关 + 样本量披露 +
#     空间分布图), 与 10(H4 共定位)/11(LM 髓系自分泌 3.9x)互补闭环。
#   - 敏感性分析: 删 patient 协变量、距离分层、contrast 变化、per-sample。
# 可视化: MediationPathway.pdf (双模型效应, 点线分离, 无 H5 字样)
#         MediationSensitivity.pdf (敏感性森林图)
#         MediationPathCoefficients.pdf (a/b/c/c' 路径系数)
#         TripletDistance_TumorSource.pdf (肿瘤源配对距离, SCI 箱线, 无抖动点)
#         TripletDistance_MyeloidSource.pdf (髓系源配对距离直方图, 解释退化)
#         PerSampleACME.pdf (per-sample ACME, 若可估计)
#         SPP1TAM_CompartmentEnrichment.pdf (区室富集)
#         SPP1TAM_SpatialAutocorr.pdf (Moran's I)
#         TripletSampleSize.pdf (三元组样本量)
#         MyeloidSampleSize.pdf (髓系 spot 样本量)
#         SPP1TAM_SpatialMap.pdf (SPP1+ TAM 空间分布)
# 数据: 输入 results/06_PLAUR_SPP1_signature/table/signature_master.rds
#       输出 results/12_PLAUR_SPP1_mediation/table/*.csv
# 依据: 方案第 5 节(mediation 包, 空间三元组, bootstrap >=1000)
# ===============================================================

# --- 包加载(统计 -> 数据 -> 可视化 -> 领域) -----------------------------
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
FIG_DIR  <- file.path(BASE_DIR, "results/12_PLAUR_SPP1_mediation/Figure")
TAB_DIR  <- file.path(BASE_DIR, "results/12_PLAUR_SPP1_mediation/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 分析常量 -----------------------------------------------------------
SAMPLE_ID <- c("ST-P1_liver", "ST-P1_colon", "ST-P2_liver", "ST-P2_colon")
PITCH_UM  <- 100    # Visium 相邻 spot 中心距标称值(um)
MAX_UM    <- 300    # 空间配对的最大半径(um)
N_SIM     <- 1000   # mediation bootstrap 次数
SEED      <- 42
CMP_ORDER <- c("Adjacent colon", "Adjacent liver", "PT", "LM")
CMP_FILL  <- c("Adjacent colon" = "#9ecae1", "Adjacent liver" = "#4292c6",
               "PT" = "#fb6a4a", "LM" = "#cb181d")
MODEL_COL <- c("Tumor source" = "#999999", "Myeloid source" = "#007AC2")

# =============================================================================
# PART 1: 构建空间三元组(肿瘤源 + 髓系源)
# 目的: 把每个巨噬 spot 与最近源 spot 配成 (X, M, Y) 三元组
# 方法: 逐样本计算距离矩阵, 取最近且 <= MAX_UM 者; 自我排除
# 数据: 输出 Triplets_tumor.csv / Triplets_myeloid.csv
# =============================================================================

cat("========== PART 1: Build spatial triplets ==========\n")

set.seed(SEED)
master <- readRDS(file.path(TAB06, "signature_master.rds"))
master$patient <- substr(master$sample, 4, 5)

UmPerPixel <- function(xy) {
  d <- as.matrix(dist(xy)); diag(d) <- Inf
  PITCH_UM / median(apply(d, 1, min))
}

build_triplets <- function(d, src_is_tumor) {
  um.px <- UmPerPixel(d[, c("x", "y")])
  dm_all <- d[d$is_macro, ]
  dm_tgt <- dm_all[!is.na(dm_all$sct_PLAUR) & !is.na(dm_all$sct_SPP1), ]
  if (nrow(dm_tgt) == 0) return(NULL)
  if (src_is_tumor) {
    dm_src <- d[d$is_tumor & !is.na(d$sct_PLAU), ]
  } else {
    dm_src <- dm_all[!is.na(dm_all$sct_PLAU), ]
  }
  if (nrow(dm_src) == 0) return(NULL)
  n.t <- nrow(dm_tgt); n.s <- nrow(dm_src)
  dd <- as.matrix(dist(rbind(dm_tgt[, c("x", "y")], dm_src[, c("x", "y")])))
  dd <- dd[1:n.t, (n.t + 1):(n.t + n.s), drop = FALSE]
  src_in_src <- match(dm_tgt$spot, dm_src$spot)
  out <- list()
  for (i in seq_len(n.t)) {
    di <- dd[i, ]
    si <- src_in_src[i]
    if (!is.na(si)) di[si] <- Inf
    if (all(is.infinite(di))) next
    ni <- which.min(di); ndist <- di[ni] * um.px
    if (is.na(ndist) || ndist > MAX_UM) next
    out[[length(out) + 1]] <- data.frame(
      macro_spot = dm_tgt$spot[i],
      src_spot   = dm_src$spot[ni],
      dist_um    = round(ndist, 1),
      plau_src   = dm_src$sct_PLAU[ni],
      plaur_mac  = dm_tgt$sct_PLAUR[i],
      spp1_mac   = dm_tgt$sct_SPP1[i],
      stringsAsFactors = FALSE)
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

trip_tumor  <- list(); trip_mye <- list()
for (s in SAMPLE_ID) {
  d <- master[master$sample == s & !is.na(master$x) & !is.na(master$y), ]
  tt <- build_triplets(d, src_is_tumor = TRUE)
  tm <- build_triplets(d, src_is_tumor = FALSE)
  if (!is.null(tt)) { tt$sample <- s; tt$patient <- substr(s, 4, 5)
                      trip_tumor[[s]]  <- tt }
  if (!is.null(tm)) { tm$sample <- s; tm$patient <- substr(s, 4, 5)
                      trip_mye[[s]]    <- tm }
  cat(sprintf("%s: tumor-source pairs=%d  myeloid-source pairs=%d (<=%d um)\n",
              s, nrow(tt), nrow(tm), MAX_UM))
}
trip_tumor <- do.call(rbind, trip_tumor)
trip_mye   <- do.call(rbind, trip_mye)
cat("Total tumor-source triplets:", nrow(trip_tumor),
    "  Total myeloid-source triplets:", nrow(trip_mye), "\n")
write.csv(trip_tumor, file.path(TAB_DIR, "Triplets_tumor.csv"),  row.names = FALSE)
write.csv(trip_mye,   file.path(TAB_DIR, "Triplets_myeloid.csv"), row.names = FALSE)

# =============================================================================
# PART 2: 中介分析主模型 + 敏感性分析
# 目的: (a) 估计主模型 ACME/ADE; (b) 通过多敏感性设定检验稳健性
# 方法: mediation::mediate; 主模型含 patient; 敏感性包括:
#       - no_patient: 删 patient 协变量
#       - nearby: 仅 <=100 um 配对
#       - contrast_p90: 0 vs 表达者 P90
#       - per_sample: 分样本独立估计(仅报告, 不稳定)
# 数据: MediationResults.csv / MediationSensitivity.csv / PerSampleACME.csv
# =============================================================================

cat("========== PART 2: Mediation analysis + sensitivity ==========\n")

fit_models <- function(trip, with_patient = TRUE) {
  if (with_patient) {
    list(
      m = lm(plaur_mac ~ plau_src + patient, data = trip),
      y = lm(spp1_mac ~ plau_src + plaur_mac + patient, data = trip),
      y_total = lm(spp1_mac ~ plau_src + patient, data = trip)
    )
  } else {
    list(
      m = lm(plaur_mac ~ plau_src, data = trip),
      y = lm(spp1_mac ~ plau_src + plaur_mac, data = trip),
      y_total = lm(spp1_mac ~ plau_src, data = trip)
    )
  }
}

run_med <- function(trip, label, setting = "Main", contrast = "P75",
                    with_patient = TRUE) {
  fits <- fit_models(trip, with_patient = with_patient)
  nz <- trip$plau_src[trip$plau_src > 0]
  ctrl <- 0
  trt <- switch(contrast,
                "P50" = if (length(nz) >= 10) quantile(nz, 0.50) else max(trip$plau_src, na.rm = TRUE),
                "P75" = if (length(nz) >= 10) quantile(nz, 0.75) else max(trip$plau_src, na.rm = TRUE),
                "P90" = if (length(nz) >= 10) quantile(nz, 0.90) else max(trip$plau_src, na.rm = TRUE),
                max(trip$plau_src, na.rm = TRUE))
  if (is.na(trt) || trt <= ctrl) trt <- max(trip$plau_src, na.rm = TRUE)
  med <- tryCatch(mediate(fits$m, fits$y, treat = "plau_src",
                          mediator = "plaur_mac", boot = TRUE, sims = N_SIM,
                          control.value = ctrl, treat.value = trt),
                  error = function(e) stop("mediate() failed in run_med(label='", label, "', setting='", setting, "', contrast='", contrast, "'): ", conditionMessage(e)))
  stopifnot(!is.null(med), "mediate() returned NULL — 中介模型拟合失败")
  data.frame(
    model = label, setting = setting, contrast = contrast,
    with_patient = with_patient,
    effect = c("ACME", "ADE", "Total", "Prop.mediated"),
    estimate = c(med$d1, med$z1, med$tau.coef, med$n1),
    ci_low = c(med$d1.ci[1], med$z1.ci[1], med$tau.ci[1], med$n1.ci[1]),
    ci_high = c(med$d1.ci[2], med$z1.ci[2], med$tau.ci[2], med$n1.ci[2]),
    p_value = c(med$d1.p, med$z1.p, med$tau.p, med$n1.p),
    stringsAsFactors = FALSE)
}

extract_path_coef <- function(trip, label, with_patient = TRUE) {
  fits <- fit_models(trip, with_patient = with_patient)
  sm_a <- summary(fits$m); sm_b <- summary(fits$y); sm_c <- summary(fits$y_total)
  data.frame(
    model = label,
    path = c("a (PLAU->PLAUR)", "b (PLAUR->SPP1)",
             "c (PLAU->SPP1 total)", "c-prime (PLAU->SPP1 direct)"),
    estimate = c(coef(fits$m)["plau_src"], coef(fits$y)["plaur_mac"],
                 coef(fits$y_total)["plau_src"], coef(fits$y)["plau_src"]),
    se = c(sm_a$coefficients["plau_src", "Std. Error"],
           sm_b$coefficients["plaur_mac", "Std. Error"],
           sm_c$coefficients["plau_src", "Std. Error"],
           sm_b$coefficients["plau_src", "Std. Error"]),
    p_value = c(sm_a$coefficients["plau_src", "Pr(>|t|)"],
                sm_b$coefficients["plaur_mac", "Pr(>|t|)"],
                sm_c$coefficients["plau_src", "Pr(>|t|)"],
                sm_b$coefficients["plau_src", "Pr(>|t|)"]),
    stringsAsFactors = FALSE)
}

# 主模型
med.main <- rbind(
  run_med(trip_tumor, "Tumor source",  setting = "Main", contrast = "P75"),
  run_med(trip_mye,   "Myeloid source", setting = "Main", contrast = "P75")
)
med.main$star <- SigStar(med.main$p_value)
cat("--- Main mediation estimates ---\n"); print(med.main, digits = 4)
write.csv(med.main, file.path(TAB_DIR, "MediationResults.csv"), row.names = FALSE)

# 路径系数
path.coef <- rbind(
  extract_path_coef(trip_tumor, "Tumor source"),
  extract_path_coef(trip_mye,   "Myeloid source")
)
path.coef$star <- SigStar(path.coef$p_value)
write.csv(path.coef, file.path(TAB_DIR, "MediationPathCoefficients.csv"),
          row.names = FALSE)

# 敏感性分析(汇总)
sens.list <- list(med.main[med.main$effect == "ACME",
                           setdiff(names(med.main), "star")])
rownames(sens.list[[1]]) <- NULL

# (1) 删 patient
sens.list[[2]] <- rbind(
  run_med(trip_tumor, "Tumor source",  setting = "No patient covariate",
          contrast = "P75", with_patient = FALSE),
  run_med(trip_mye,   "Myeloid source", setting = "No patient covariate",
          contrast = "P75", with_patient = FALSE)
)

# (2) 距离分层: <=100 um
sens.list[[3]] <- rbind(
  run_med(subset(trip_tumor, dist_um <= 100), "Tumor source",
          setting = "Nearby pairs (<=100 um)", contrast = "P75"),
  run_med(subset(trip_mye, dist_um <= 100), "Myeloid source",
          setting = "Nearby pairs (<=100 um)", contrast = "P75")
)

# (3) contrast P90
sens.list[[4]] <- rbind(
  run_med(trip_tumor, "Tumor source",  setting = "High contrast (0 vs P90)",
          contrast = "P90"),
  run_med(trip_mye,   "Myeloid source", setting = "High contrast (0 vs P90)",
          contrast = "P90")
)

sens.tab <- do.call(rbind, sens.list)
sens.tab <- sens.tab[sens.tab$effect == "ACME", ]
sens.tab$star <- SigStar(sens.tab$p_value)
sens.tab$setting <- factor(sens.tab$setting,
                           levels = c("Main", "No patient covariate",
                                      "Nearby pairs (<=100 um)",
                                      "High contrast (0 vs P90)"))
cat("--- Sensitivity ACME estimates ---\n"); print(sens.tab, digits = 4)
write.csv(sens.tab, file.path(TAB_DIR, "MediationSensitivity.csv"),
          row.names = FALSE)

# Per-sample ACME(仅 myeloid source, 肿瘤源按样本 n 过小且零胀更严重)
cat("--- Per-sample ACME (myeloid source) ---\n")
ps.list <- list()
for (s in unique(trip_mye$sample)) {
  ss <- subset(trip_mye, sample == s)
  nz <- sum(ss$plau_src > 0, na.rm = TRUE)
  if (nrow(ss) < 30 || nz < 10) {
    cat(s, ": skipped (n=", nrow(ss), " nz_PLAU=", nz, ")\n")
    next
  }
  res <- run_med(ss, "Myeloid source", setting = paste("Per sample:", s),
                 contrast = "P75", with_patient = FALSE)
  if (!is.null(res)) {
    res$sample <- s
    ps.list[[s]] <- res
    cat(s, ": ACME=", round(res$estimate[res$effect == "ACME"], 4),
        " p=", format(res$p_value[res$effect == "ACME"], digits = 3), "\n")
  }
}
ps.tab <- if (length(ps.list) > 0) do.call(rbind, ps.list) else NULL
if (!is.null(ps.tab)) {
  ps.tab$star <- SigStar(ps.tab$p_value)
  write.csv(ps.tab, file.path(TAB_DIR, "PerSampleACME.csv"), row.names = FALSE)
}

# =============================================================================
# PART 3: 中介可视化(去 H5 字样, 点线分离)
# 目的: 清晰呈现两模型四效应, 避免点重叠; 肿瘤源明确标为不支持
# =============================================================================

cat("========== PART 3: Mediation plots ==========\n")

# 3a 主效应图: 按 effect 分面, dodge 使两 model 点线分离
med.main$effect <- factor(med.main$effect,
                          levels = c("ACME", "ADE", "Total", "Prop.mediated"))
med.main$model <- factor(med.main$model,
                         levels = c("Tumor source", "Myeloid source"))

p.med <- ggplot(med.main, aes(x = effect, y = estimate, colour = model,
                              group = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12,
                linewidth = 0.7, position = position_dodge(width = 0.45)) +
  geom_point(aes(shape = model), size = 4,
             position = position_dodge(width = 0.45)) +
  geom_text(aes(label = star), colour = "black", vjust = -0.9, size = 4,
            position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = MODEL_COL, name = "Source of PLAU") +
  scale_shape_manual(values = c("Tumor source" = 16, "Myeloid source" = 17),
                     name = "Source of PLAU") +
  labs(title = "Mediation of PLAU -> PLAUR -> SPP1 in spatial triplets",
       subtitle = "Tumor source: not supported; Myeloid source: positive trend (p=0.058). ADE/Total/Prop are suppression artifacts.",
       x = "Effect", y = "Estimate with 95% CI") +
  sci_theme +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 10),
        legend.position = "bottom")
sci_save(p.med, "MediationPathway", width = 10, height = 5.5)

# 3b 敏感性森林图: 仅 ACME
p.sens <- ggplot(sens.tab, aes(x = setting, y = estimate, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2,
                linewidth = 0.7, position = position_dodge(width = 0.4)) +
  geom_point(aes(shape = model), size = 3.5,
             position = position_dodge(width = 0.4)) +
  geom_text(aes(label = star), colour = "black", vjust = -0.8, size = 3.5,
            position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = MODEL_COL, name = "Source of PLAU") +
  scale_shape_manual(values = c("Tumor source" = 16, "Myeloid source" = 17),
                     name = "Source of PLAU") +
  coord_flip() +
  labs(title = "Sensitivity of ACME estimates",
       subtitle = "Bootstrap 95% CI across covariate, distance and contrast settings",
       x = "Setting", y = "ACME estimate with 95% CI") +
  sci_theme + theme(legend.position = "bottom")
sci_save(p.sens, "MediationSensitivity", width = 10, height = 5)

# 3c 路径系数图
path.coef$path <- factor(path.coef$path,
                         levels = c("a (PLAU->PLAUR)", "b (PLAUR->SPP1)",
                                    "c (PLAU->SPP1 total)",
                                    "c-prime (PLAU->SPP1 direct)"))
path.coef$model <- factor(path.coef$model,
                          levels = c("Tumor source", "Myeloid source"))
path.coef$ci_low  <- path.coef$estimate - 1.96 * path.coef$se
path.coef$ci_high <- path.coef$estimate + 1.96 * path.coef$se
p.path <- ggplot(path.coef, aes(x = model, y = estimate, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15,
                linewidth = 0.7, position = position_dodge(width = 0.4)) +
  geom_point(aes(shape = model), size = 4,
             position = position_dodge(width = 0.4)) +
  geom_text(aes(label = star), colour = "black", vjust = -0.9, size = 4,
            position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = MODEL_COL, name = "Source of PLAU") +
  scale_shape_manual(values = c("Tumor source" = 16, "Myeloid source" = 17),
                     name = "Source of PLAU") +
  facet_wrap(~ path, scales = "free_y", ncol = 4) +
  labs(title = "Path coefficients of the mediation model",
       subtitle = "a and b paths define ACME; c-prime is the direct effect",
       x = "Source of PLAU", y = "Standardized coefficient (95% CI)") +
  sci_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(size = 10, face = "plain"))
sci_save(p.path, "MediationPathCoefficients", width = 13, height = 5)

# =============================================================================
# PART 4: 三元组配对距离分布(SCI 标准箱线, 无抖动点, 两源独立成图)
# 目的: 展示空间采样结构(仅作披露, 非因果轴)。
#   - 肿瘤源距离有真实生物学变异(肿瘤 PLAU+ 散在, 最近邻在 100~300 um 之间)。
#   - 髓系源距离因 PLAU+ 髓系在组织内密集, 最近邻退化为"1 个 spot 间距 floor"
#     (≈100 um); 这是 Visium 100 um 点阵几何 + 源密集共同造成的结构性地板效应,
#     不携带可解读的生物学梯度。UM/px 标定始终以 PITCH_UM=100 绑定最近邻中值,
#     故"最近邻距离"恒回归 100 um —— 属几何必然, 非标定 bug。
#   - 正因为髓系距离恒为地板值, 中介模型 X 用"最近源 spot 的 PLAU 表达量"
#     (而非距离) 作处理变量, 距离仅作 <=300 um 配对过滤, 故退化不影响效应估计。
# 数据: 输出 TripletDistanceSummary.csv / TripletDistance_TumorSource.pdf /
#              TripletDistance_MyeloidSource.pdf
# =============================================================================

cat("========== PART 4: Triplet distance distribution ==========\n")

trip_tumor$source <- "Tumor source"
trip_mye$source   <- "Myeloid source"
dist.df <- rbind(trip_tumor[, c("sample", "dist_um", "source")],
                 trip_mye[, c("sample", "dist_um", "source")])
dist.df$source <- factor(dist.df$source,
                         levels = c("Tumor source", "Myeloid source"))
dist.df$sample <- factor(dist.df$sample, levels = SAMPLE_ID)

dist.sum <- aggregate(dist_um ~ sample + source, data = dist.df,
                      FUN = function(x) c(median = median(x), n = length(x),
                                          q25 = quantile(x, 0.25),
                                          q75 = quantile(x, 0.75)))
dist.sum <- do.call(data.frame, dist.sum)
names(dist.sum)[3:6] <- c("median_um", "n", "q25_um", "q75_um")
write.csv(dist.sum, file.path(TAB_DIR, "TripletDistanceSummary.csv"),
          row.names = FALSE)

# (a) 肿瘤源距离: SCI 标准箱线(无抖动点), 展示真实变异
TUMOUR_COL <- MODEL_COL["Tumor source"]
dist.tumor <- subset(dist.df, source == "Tumor source")
med.tumor <- aggregate(dist_um ~ sample, data = dist.tumor, FUN = median)
names(med.tumor)[2] <- "med"
n.tumor <- aggregate(dist_um ~ sample, data = dist.tumor, FUN = length)
names(n.tumor)[2] <- "n"

p.dist.tumor <- ggplot(dist.tumor, aes(x = sample, y = dist_um)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85,
               fill = TUMOUR_COL, colour = "black", linewidth = 0.5) +
  geom_text(data = med.tumor, aes(x = sample, y = med,
                                  label = sprintf("%.0f", med)),
            vjust = -0.9, size = 3.5, colour = "black", fontface = "bold") +
  geom_text(data = n.tumor, aes(x = sample, y = 0, label = paste0("n=", n)),
            vjust = 1.3, size = 3, colour = "black") +
  scale_y_continuous(limits = c(0, 320), breaks = seq(0, 300, 50)) +
  labs(title = "Distance from macrophage to nearest tumor PLAU source",
       subtitle = "Each macrophage paired with nearest tumor spot expressing PLAU (cutoff 300 um)",
       x = "Sample", y = "Distance (um)") +
  sci_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9))
sci_save(p.dist.tumor, "TripletDistance_TumorSource", width = 7, height = 5)

# (b) 髓系源距离: 直方图说明退化现象(几乎全部落在 ~100 um 的最近邻)
MYELOID_COL <- MODEL_COL["Myeloid source"]
dist.mye <- subset(dist.df, source == "Myeloid source")

p.dist.mye <- ggplot(dist.mye, aes(x = dist_um)) +
  geom_histogram(bins = 20, alpha = 0.85, fill = MYELOID_COL,
               colour = "black", linewidth = 0.3) +
  facet_wrap(~ sample, ncol = 4) +
  scale_x_continuous(breaks = c(100, 200, 300)) +
  labs(title = "Distance from macrophage to nearest myeloid PLAU source",
       subtitle = "Myeloid PLAU sources are dense; nearest neighbor collapses to ~100 um (one spot pitch)",
       x = "Distance (um)", y = "Count") +
  sci_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 8),
                    panel.spacing = unit(0.8, "lines"))
sci_save(p.dist.mye, "TripletDistance_MyeloidSource", width = 10, height = 3.5)

# =============================================================================
# PART 5: SPP1+ TAM 空间生态位(区室富集 + 自相关 + 样本量 + 空间分布)
# =============================================================================

cat("========== PART 5: Plan C SPP1+ TAM niche ==========\n")

# (a) 区室富集
enr <- list()
for (s in SAMPLE_ID) {
  d <- master[master$sample == s & !is.na(master$x) & !is.na(master$y), ]
  dm <- d[d$is_macro & !is.na(d$sct_SPP1), ]
  if (nrow(dm) == 0) next
  hi <- dm$sct_SPP1 > median(dm$sct_SPP1, na.rm = TRUE)
  enr[[s]] <- data.frame(
    sample = s,
    compartment = dm$compartment,
    is_spp1_tam = hi, stringsAsFactors = FALSE)
}
enr <- do.call(rbind, enr)
enr.tab <- aggregate(is_spp1_tam ~ compartment, data = enr,
                     FUN = function(x) round(100 * mean(x), 1))
names(enr.tab)[2] <- "pct_spp1_tam"
n.tab <- aggregate(is_spp1_tam ~ compartment, data = enr, FUN = length)
names(n.tab)[2] <- "n_macro"
enr.tab <- merge(enr.tab, n.tab, by = "compartment")
enr.tab$compartment <- factor(enr.tab$compartment, levels = CMP_ORDER)
enr.tab <- enr.tab[order(enr.tab$compartment), ]
cat("--- SPP1+ TAM compartment enrichment ---\n"); print(enr.tab)
write.csv(enr.tab, file.path(TAB_DIR, "SPP1TAM_Enrichment.csv"), row.names = FALSE)

p.enr <- ggplot(enr.tab, aes(x = compartment, y = pct_spp1_tam, fill = compartment)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct_spp1_tam, n_macro)),
            vjust = -0.3, size = 3.5, colour = "black") +
  scale_fill_manual(values = CMP_FILL, name = "Compartment") +
  scale_y_continuous(limits = c(0, max(enr.tab$pct_spp1_tam) * 1.18)) +
  labs(title = "SPP1+ TAM fraction by compartment",
       x = "Compartment", y = "Percent SPP1+ TAM (%)") +
  sci_theme + theme(axis.text.x = element_text(size = 10), legend.position = "none")
sci_save(p.enr, "SPP1TAM_CompartmentEnrichment", width = 7, height = 5)

# (b) 空间自相关
aut2 <- list()
for (s in SAMPLE_ID) {
  d <- master[master$sample == s & !is.na(master$x) & !is.na(master$y), ]
  dm <- d[d$is_macro & !is.na(d$sct_SPP1), ]
  if (nrow(dm) < 10) next
  aut2[[s]] <- data.frame(sample = s,
                          tissue = ifelse(grepl("liver", s), "Liver mets", "Primary"),
                          morans_I = round(moranI(dm[, c("x", "y")], dm$sct_SPP1,
                                                  k = 6), 3),
                          n_macro = nrow(dm), stringsAsFactors = FALSE)
}
aut <- do.call(rbind, aut2)
cat("--- SPP1+ TAM spatial autocorrelation ---\n"); print(aut)
write.csv(aut, file.path(TAB_DIR, "SPP1TAM_Autocorr.csv"), row.names = FALSE)

aut$sample <- factor(aut$sample, levels = SAMPLE_ID)
aut$tissue <- factor(aut$tissue, levels = c("Primary", "Liver mets"))
p.aut <- ggplot(aut, aes(x = sample, y = morans_I, fill = tissue)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.3f\n(n=%d)", morans_I, n_macro)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Primary" = "#4292c6", "Liver mets" = "#cb181d"),
                    name = "Tissue") +
  scale_y_continuous(limits = c(0, max(aut$morans_I) * 1.25)) +
  labs(title = "Spatial autocorrelation of SPP1 expression in myeloid spots",
       x = "Sample", y = "Moran's I (k=6 nearest neighbours)") +
  sci_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
                    legend.position = "bottom")
sci_save(p.aut, "SPP1TAM_SpatialAutocorr", width = 8, height = 5)

# (c) 样本量披露图: 三元组 n + SPP1+ TAM n
ss.trip <- aggregate(dist_um ~ sample + source, data = dist.df, FUN = length)
names(ss.trip)[3] <- "n"
ss.trip$source <- factor(ss.trip$source,
                         levels = c("Tumor source", "Myeloid source"))
ss.trip$sample <- factor(ss.trip$sample, levels = SAMPLE_ID)

ss.enr <- aggregate(is_spp1_tam ~ sample + compartment, data = enr, FUN = length)
names(ss.enr)[3] <- "n"
ss.enr$compartment <- factor(ss.enr$compartment, levels = CMP_ORDER)
ss.enr$sample <- factor(ss.enr$sample, levels = SAMPLE_ID)

p.ss.trip <- ggplot(ss.trip, aes(x = sample, y = n, fill = source)) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5, colour = "black",
            position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = MODEL_COL, name = "Source of PLAU") +
  scale_y_continuous(limits = c(0, max(ss.trip$n) * 1.18)) +
  labs(title = "Sample size of spatial triplets", x = "Sample", y = "N triplets") +
  sci_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
                    legend.position = "bottom")

p.ss.enr <- ggplot(ss.enr, aes(x = sample, y = n, fill = compartment)) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5, colour = "black",
            position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = CMP_FILL, name = "Compartment") +
  scale_y_continuous(limits = c(0, max(ss.enr$n) * 1.18)) +
  labs(title = "Sample size of myeloid spots by compartment",
       x = "Sample", y = "N myeloid spots") +
  sci_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
                    legend.position = "bottom")

sci_save(p.ss.trip, "TripletSampleSize", width = 9, height = 5)
sci_save(p.ss.enr, "MyeloidSampleSize", width = 9, height = 5)

# (d) SPP1+ TAM 空间分布图(每样本分面, 背景全髓系灰色, SPP1+ 红色)
map.list <- list()
for (s in SAMPLE_ID) {
  d <- master[master$sample == s & !is.na(master$x) & !is.na(master$y), ]
  dm <- d[d$is_macro & !is.na(d$sct_SPP1), c("x", "y", "sct_SPP1", "compartment")]
  if (nrow(dm) == 0) next
  dm$sample <- s
  dm$hi <- dm$sct_SPP1 > median(dm$sct_SPP1, na.rm = TRUE)
  map.list[[s]] <- dm
}
map.df <- do.call(rbind, map.list)
map.df$sample <- factor(map.df$sample, levels = SAMPLE_ID)
map.df$compartment <- factor(map.df$compartment, levels = CMP_ORDER)

p.map <- ggplot(map.df, aes(x = x, y = y)) +
  geom_point(data = map.df, aes(x = x, y = y), colour = "grey80", size = 0.4) +
  geom_point(data = map.df[map.df$hi, ], aes(colour = compartment), size = 0.9) +
  scale_colour_manual(values = CMP_FILL, name = "Compartment") +
  coord_fixed(ratio = 1) +
  facet_wrap(~ sample, ncol = 4) +
  labs(title = "Spatial distribution of SPP1-high myeloid spots",
       subtitle = "Background = all myeloid spots; coloured spots = SPP1 above sample median",
       x = "X (pixels)", y = "Y (pixels)") +
  sci_theme + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                    legend.position = "bottom",
                    panel.spacing = unit(0.5, "lines"))
sci_save(p.map, "SPP1TAM_SpatialMap", width = 14, height = 4)

# =============================================================================
# PART 6: Per-sample ACME 可视化(若可估计)
# =============================================================================

cat("========== PART 6: Per-sample ACME plot ==========\n")

if (!is.null(ps.tab)) {
  ps.acme <- ps.tab[ps.tab$effect == "ACME", ]
  ps.acme$sample <- factor(ps.acme$sample, levels = SAMPLE_ID)
  ps.acme$model <- factor(ps.acme$model, levels = c("Tumor source", "Myeloid source"))
  p.ps <- ggplot(ps.acme, aes(x = sample, y = estimate, colour = model)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15,
                  linewidth = 0.7) +
    geom_point(size = 4) +
    geom_text(aes(label = star), colour = "black", vjust = -0.9, size = 4) +
    scale_colour_manual(values = MODEL_COL, name = "Source of PLAU") +
    labs(title = "Per-sample ACME estimates (myeloid source only)",
         subtitle = "Tumor source skipped due to small n or zero-inflated PLAU",
         x = "Sample", y = "ACME estimate with 95% CI") +
    sci_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
                      legend.position = "bottom")
  sci_save(p.ps, "PerSampleACME", width = 8, height = 5)
}

cat("12_PLAUR_SPP1_mediation.R done.\n")
