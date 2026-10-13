# ===============================================================
# 13_PLAUR_SPP1_robustness.R
# 功能: 统一稳健性层——覆盖 09(区室趋势)/10(共定位)/11(CellChat 方向)/12(中介)
#       的阳性阈值敏感性、留一样本(LOSO)、留一患者(patient-LOSO)、
#       效应量 bootstrap 95% CI、关键比较敏感性森林
# 定位: 方案第 8 节(稳健性)与 S-fig 1-5
# 可视化: RobustnessLOSO.pdf (09 区室 LOSO 稳定性)
#         RobustnessEffectCI.pdf (09 区室 bootstrap 95% CI)
#         RobustnessThreshold.pdf (10 共定位阈值敏感性)
#         RobustnessCellChat.pdf (11 PLAUR 通讯方向稳健性)
#         RobustnessMediation.pdf (12 中介 ACME 稳健性森林)
# 数据: 输入 results/06_.../table/signature_master.rds
#            results/11_.../table/CellChatPLAUR.csv
#            results/12_.../table/Triplets_tumor.csv / Triplets_myeloid.csv
#            results/12_.../table/MediationSensitivity.csv
#       输出 results/13_.../table/Robustness{Threshold,LOSO,EffectCI,CellChat,Mediation}.csv
# 依据: 方案第 8 节: 三阈值敏感性(q=0.50/0.75/0.90)、LOSO、patient-LOSO、
#       样本级聚类 bootstrap 1000 次、效应量 + 95% CI 优先于 p 值
# 说明: 本脚本是 09/10/11/12 的统一稳健性层, 以"改变分析条件后结论是否稳定"
#       为唯一目的, 不产生新的生物学结论; 11 的方向稳健性为结论确认 + LM/PT
#       比值(逐样本 CellChat 重跑属 11 管线, 不在本脚本重复)
# ===============================================================

# --- 包加载(统计 -> 数据 -> 可视化 -> 领域) -----------------------------
suppressMessages({
  library(ggplot2)     # 可视化
  library(matrixStats) # 行最小值(最近邻距离)高效计算
  library(mediation)   # 中介分析 bootstrap
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

# CliffDelta: 由 Wilcoxon 秩和统计量换算 Cliff's delta = 2*AUC - 1
CliffDelta <- function(a, b) {
  wt <- wilcox.test(a, b, exact = FALSE)
  2 * (wt$statistic / (length(a) * length(b))) - 1
}

# SigStar: p 值 -> 显著性标记
SigStar <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

# --- 目录结构 -----------------------------------------------------------
BASE_DIR <- "X:/XX"
TAB06    <- file.path(BASE_DIR, "results/06_PLAUR_SPP1_signature/table")
TAB11    <- file.path(BASE_DIR, "results/11_PLAUR_SPP1_communication/table")
TAB12    <- file.path(BASE_DIR, "results/12_PLAUR_SPP1_mediation/table")
FIG_DIR  <- file.path(BASE_DIR, "results/13_PLAUR_SPP1_robustness/Figure")
TAB_DIR  <- file.path(BASE_DIR, "results/13_PLAUR_SPP1_robustness/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 分析常量 -----------------------------------------------------------
# 原发结直肠肿瘤记为 PT (primary tumor)
SAMPLE_ID <- c("ST-P1_liver", "ST-P1_colon", "ST-P2_liver", "ST-P2_colon")
# 双层区室比较变量集(方案 v1.1): 与 08/09 保持一致
# 理由: FBXL12_tumor / PLAU_tumor 在 Adjacent colon / Adjacent liver 的
#       spot 数 = 0, 无法用于区室间比较(LOSO / bootstrap CI 同理)
GENE_RAW  <- c("sct_FBXL12", "sct_PLAU", "sct_PLAUR", "sct_SPP1")
MAC_SIG   <- c("PLAUR_mac", "SPP1_mac")
CMP_SET   <- c(GENE_RAW, MAC_SIG)
Q_GRID    <- c(0.50, 0.75, 0.90)   # 阳性阈值敏感性网格
N_BOOT    <- 1000
N_PERM    <- 999
SEED      <- 42
N_SIM     <- 1000   # mediation bootstrap 次数(与 12 对齐)
MODEL_COL <- c("Tumor source" = "#999999", "Myeloid source" = "#007AC2")

# =============================================================================
# PART 1: 阳性阈值敏感性(共定位的方向是否随阈值改变)
# 目的: 检验 H2/H4 的共定位结论对"高表达"定义是否稳健
# 方法: 在 q = 0.50/0.75/0.90 三档下重算平均最近邻跨类距离与置换 p 值
# 数据: 输出 RobustnessThreshold.csv
# =============================================================================

cat("========== PART 1: Threshold sensitivity ==========\n")

set.seed(SEED)
master <- readRDS(file.path(TAB06, "signature_master.rds"))

# MeanNnDist 已内联为 precompute-D-once 版本(见 PART 1 spec 循环):
# 完整成对距离矩阵仅计算一次, 置换重抽样只做索引, 否则 O(n^2) 被 999 次放大致超时

thr.rows <- list()
for (s in SAMPLE_ID) {
  d <- master[master$sample == s, ]
  d <- d[!is.na(d$x) & !is.na(d$y), ]
  for (q in Q_GRID) {
    dt <- d[d$is_tumor & !is.na(d$sct_FBXL12), ]
    dm <- d[d$is_macro & !is.na(d$sct_PLAUR) & !is.na(d$sct_SPP1), ]
    pa <- dt[dt$sct_FBXL12 >= quantile(dt$sct_FBXL12, q, names = FALSE),
             c("x", "y")]
    pb <- dm[dm$sct_PLAUR >= quantile(dm$sct_PLAUR, q, names = FALSE),
             c("x", "y")]
    pc <- dm[dm$sct_SPP1 >= quantile(dm$sct_SPP1, q, names = FALSE),
             c("x", "y")]

    spec <- list(list(id = "H2_FBXL12_PLAUR", p1 = pa, p2 = pb),
                 list(id = "H4_PLAUR_SPP1", p1 = pb, p2 = pc))
    for (sp in spec) {
      pool <- rbind(sp$p1, sp$p2)
      n1 <- nrow(sp$p1); n2 <- nrow(pool) - n1
      if (n1 < 2 || n2 < 2) {
        obs <- NA_real_; null <- NA_real_; p.lo <- NA_real_
      } else {
        # 预先计算完整成对距离矩阵一次, 置换仅做索引(避免每次重算 dist)
        D <- as.matrix(dist(pool))
        obs <- mean(rowMins(D[1:n1, (n1 + 1):nrow(D), drop = FALSE]))
        null <- replicate(N_PERM, {
          idx <- sample(nrow(pool))
          a <- idx[1:n1]; b <- idx[(n1 + 1):nrow(pool)]
          mean(rowMins(D[a, b, drop = FALSE]))
        })
        p.lo <- (sum(null <= obs, na.rm = TRUE) + 1) / (N_PERM + 1)
      }
      thr.rows[[length(thr.rows) + 1]] <- data.frame(
        sample = s, q_cut = q, test = sp$id, n1 = n1,
        n2 = nrow(pool) - n1, nn_obs = round(obs, 1),
        nn_null = round(mean(null), 1),
        ratio = round(obs / mean(null), 3),
        p_closer = signif(p.lo, 4), stringsAsFactors = FALSE
      )
    }
  }
}
thr.tab <- do.call(rbind, thr.rows)
write.csv(thr.tab, file.path(TAB_DIR, "RobustnessThreshold.csv"),
          row.names = FALSE)
print(thr.tab)

# =============================================================================
# PART 2: LOSO(逐样本留一)后重算关键区室比较
# 目的: 检验 H1/H6 的结论是否被单个样本驱动
# 方法: 每次剔除一个样本, 重算 PT vs Adjacent colon、LM vs Adjacent
#       liver、LM vs PT 的 Cliff's delta(spot 级, 方案 v1.1)
# 数据: 输出 RobustnessLOSO.csv
# =============================================================================

cat("========== PART 2: Leave-one-sample-out ==========\n")

PAIRS <- list(c("PT", "Adjacent colon"), c("LM", "Adjacent liver"),
              c("LM", "PT"))

loso.rows <- list()
for (drop.s in SAMPLE_ID) {
  d <- master[master$sample != drop.s, ]
  for (sig in CMP_SET) {
    for (pr in PAIRS) {
      a <- d[[sig]][d$compartment == pr[1]]
      b <- d[[sig]][d$compartment == pr[2]]
      a <- a[!is.na(a)]
      b <- b[!is.na(b)]
      loso.rows[[length(loso.rows) + 1]] <- data.frame(
        dropped = drop.s, signature = sig,
        comparison = paste(pr[1], "vs", pr[2]),
        n1 = length(a), n2 = length(b),
        cliff_delta = round(CliffDelta(a, b), 4),
        stringsAsFactors = FALSE
      )
    }
  }
}
loso.tab <- do.call(rbind, loso.rows)
write.csv(loso.tab, file.path(TAB_DIR, "RobustnessLOSO.csv"),
          row.names = FALSE)
cat("LOSO Cliff's delta 汇总(按比较分组的范围):\n")
print(aggregate(cliff_delta ~ signature + comparison, data = loso.tab,
                FUN = function(v) round(range(v), 3)))

# =============================================================================
# PART 3: 关键比较的效应量 bootstrap 95% CI(样本级聚类重抽样, 方案 v1.1)
# 目的: n=4 下以效应量区间替代单一 p 值作为结论依据
# 方法: 样本级聚类 bootstrap(每次抽 4 样本)重算 Cliff's delta, 取分位
# 数据: 输出 RobustnessEffectCI.csv
# =============================================================================

cat("========== PART 3: Bootstrap CI for effect size ==========\n")

masterc <- master[!is.na(master$compartment), ]

# BootCliffClust: 一次样本级聚类重抽样后计算指定比较的 Cliff's delta
BootCliffClust <- function(sig, g1, g2) {
  bs <- sample(SAMPLE_ID, length(SAMPLE_ID), replace = TRUE)
  av <- unlist(lapply(bs, function(s) {
    masterc[[sig]][masterc$sample == s & masterc$compartment == g1]
  }))
  bv <- unlist(lapply(bs, function(s) {
    masterc[[sig]][masterc$sample == s & masterc$compartment == g2]
  }))
  av <- av[!is.na(av)]; bv <- bv[!is.na(bv)]
  if (length(av) < 2 || length(bv) < 2) NA_real_
  else CliffDelta(av, bv)
}

ci.rows <- list()
for (sig in CMP_SET) {
  for (pr in PAIRS) {
    reps <- replicate(N_BOOT, BootCliffClust(sig, pr[1], pr[2]))
    ci.rows[[length(ci.rows) + 1]] <- data.frame(
      signature = sig, comparison = paste(pr[1], "vs", pr[2]),
      cliff_delta = round(mean(reps, na.rm = TRUE), 4),
      ci_low = round(quantile(reps, 0.025, names = FALSE, na.rm = TRUE), 4),
      ci_high = round(quantile(reps, 0.975, names = FALSE, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )
  }
}
ci.tab <- do.call(rbind, ci.rows)
write.csv(ci.tab, file.path(TAB_DIR, "RobustnessEffectCI.csv"),
          row.names = FALSE)
print(ci.tab)

# =============================================================================
# PART 4: CellChat PLAUR 通讯方向稳健性(脚本 11)
# 目的: 确认"髓系内 uPAR 自分泌为 PLAUR 主方向、肝转移富集(LM/PT 比值)"结论稳定
# 方法: 读 11 的 CellChatPLAUR.csv, 检验每区室显著方向排名 + LM/PT 概率比值
# 数据: 输出 RobustnessCellChat.csv / RobustnessCellChat.pdf
# 说明: 11 按区室汇总(pooling 样本), 逐样本 LOSO 需重跑 11 的 CellChat 管线,
#       不在本脚本重复; 此处为结论确认 + 比值稳健性
# =============================================================================

cat("========== PART 4: CellChat direction robustness (script 11) ==========\n")

cc <- read.csv(file.path(TAB11, "CellChatPLAUR.csv"))
cc$sig <- cc$pval < 0.05
cc$dir <- paste(cc$source, "->", cc$target)

# 主方向 Myeloid->Myeloid 的 LM/PT 概率比值(肝转移富集倍数)
mm <- cc[cc$dir == "Myeloid -> Myeloid" & cc$sig, ]
if (nrow(mm) == 2) {
  ratio_lm_pt <- mm$prob[mm$compartment == "LM"] / mm$prob[mm$compartment == "PT"]
  cat(sprintf("Myeloid->Myeloid PLAUR 通讯 LM/PT 比值 = %.2f (LM=%.5f, PT=%.5f)\n",
              ratio_lm_pt, mm$prob[mm$compartment == "LM"],
              mm$prob[mm$compartment == "PT"]))
} else {
  ratio_lm_pt <- NA_real_
  cat("Myeloid->Myeloid 在 PT/LM 未同时显著, 比值不适用\n")
}

# 方向稳健性: 各方向是否在两区室均显著(排名稳定性)
cc.rob <- aggregate(sig ~ dir, data = cc, FUN = function(v) all(v))
names(cc.rob)[2] <- "sig_both_compartments"
write.csv(cc, file.path(TAB_DIR, "RobustnessCellChat.csv"), row.names = FALSE)
cat("CellChat 方向稳健性(各方向在 PT/LM 均显著):\n"); print(cc.rob)

p.cc <- ggplot(cc, aes(x = dir, y = prob, fill = sig)) +
  geom_col() +
  geom_text(aes(label = ifelse(sig, sprintf("%.4f", prob), "n.s.")),
            vjust = -0.3, size = 3, colour = "black") +
  facet_wrap(~ compartment, ncol = 2) +
  scale_fill_manual(values = c("FALSE" = "#9ecae1", "TRUE" = "#cb181d"),
                    labels = c("FALSE" = "n.s.", "TRUE" = "p<0.05"),
                    name = "Significance") +
  labs(title = "Robustness of PLAUR communication direction",
       subtitle = ifelse(is.na(ratio_lm_pt),
                         "Myeloid->Myeloid dominant axis",
                         paste0("Myeloid->Myeloid dominant; LM/PT ratio = ",
                                round(ratio_lm_pt, 2))),
       x = "Direction (source -> target)", y = "Communication probability") +
  sci_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 8),
                    legend.position = "bottom")
sci_save(p.cc, "RobustnessCellChat", width = 9, height = 6)

# =============================================================================
# PART 5: 中介 ACME 稳健性(脚本 12)
# 目的: 检验 H5(肿瘤源)被拒、方案 A(髓系源)正向结论对 patient 混杂是否稳健
# 方法: 读 12 的三元组 CSV, patient-LOSO(剔除整个 P1 或 P2 患者)重跑 mediate;
#       与 12 既有敏感性(Main/无 patient/≤100um/P90)合并为 ACME 森林
# 数据: 输出 RobustnessMediation.csv / RobustnessMediation.pdf
# =============================================================================

cat("========== PART 5: Mediation ACME robustness (script 12) ==========\n")

trip_tumor <- read.csv(file.path(TAB12, "Triplets_tumor.csv"))
trip_mye  <- read.csv(file.path(TAB12, "Triplets_myeloid.csv"))
sens12 <- read.csv(file.path(TAB12, "MediationSensitivity.csv"))

# 内联中介估计器(与 12 同源, 规范22 自给自足)
fit_models <- function(trip, with_patient = TRUE) {
  if (with_patient) {
    list(m = lm(plaur_mac ~ plau_src + patient, data = trip),
         y = lm(spp1_mac ~ plau_src + plaur_mac + patient, data = trip),
         y_total = lm(spp1_mac ~ plau_src + patient, data = trip))
  } else {
    list(m = lm(plaur_mac ~ plau_src, data = trip),
         y = lm(spp1_mac ~ plau_src + plaur_mac, data = trip),
         y_total = lm(spp1_mac ~ plau_src, data = trip))
  }
}
run_med <- function(trip, label, setting = "Main", contrast = "P75",
                    with_patient = TRUE) {
  fits <- fit_models(trip, with_patient = with_patient)
  nz <- trip$plau_src[trip$plau_src > 0]
  ctrl <- 0
  trt <- switch(contrast,
                "P75" = if (length(nz) >= 10) quantile(nz, 0.75) else max(trip$plau_src, na.rm = TRUE),
                "P90" = if (length(nz) >= 10) quantile(nz, 0.90) else max(trip$plau_src, na.rm = TRUE),
                max(trip$plau_src, na.rm = TRUE))
  if (is.na(trt) || trt <= ctrl) trt <- max(trip$plau_src, na.rm = TRUE)
  med <- tryCatch(mediate(fits$m, fits$y, treat = "plau_src", mediator = "plaur_mac",
                          boot = TRUE, sims = N_SIM, control.value = ctrl,
                          treat.value = trt),
                  error = function(e) {
                    warning("mediate failed (", label, "/", setting, "): ",
                            conditionMessage(e)); NULL
                  })
  if (is.null(med)) {
    warning("run_med returned NULL (", label, "/", setting, "); filling NA row")
    return(data.frame(model = label, setting = setting, contrast = contrast,
                      estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                      p_value = NA_real_, stringsAsFactors = FALSE))
  }
  data.frame(model = label, setting = setting, contrast = contrast,
             estimate = med$d1, ci_low = med$d1.ci[1], ci_high = med$d1.ci[2],
             p_value = med$d1.p, stringsAsFactors = FALSE)
}

# patient-LOSO: 剔除单一患者(该患者的 liver+colon 两样本一并剔除)
# 单患者子集内 patient 仅 1 层 -> 强制 with_patient = FALSE
set.seed(SEED)   # 固定 bootstrap 随机性, 保证 patient-LOSO 结果可复现
med.rob <- list()
for (pat in c("DROP_P1", "DROP_P2")) {
  pid <- sub("DROP_", "", pat)
  tt <- trip_tumor[trip_tumor$patient != pid, ]
  tm <- trip_mye[trip_mye$patient != pid, ]
  med.rob[[pat]] <- rbind(
    run_med(tt, "Tumor source",  setting = pat, with_patient = FALSE),
    run_med(tm, "Myeloid source", setting = pat, with_patient = FALSE))
}
med.rob <- do.call(rbind, med.rob)
med.rob$star <- SigStar(med.rob$p_value)
write.csv(med.rob, file.path(TAB_DIR, "RobustnessMediation.csv"), row.names = FALSE)
cat("--- Mediation ACME under patient-LOSO ---\n"); print(med.rob, digits = 4)

# 合并 12 既有敏感性 + 本 PART patient-LOSO -> 统一 ACME 森林
s12 <- sens12[sens12$effect == "ACME",
              c("model", "setting", "estimate", "ci_low", "ci_high", "p_value")]
s12$star <- SigStar(s12$p_value)
s12$setting <- as.character(s12$setting)
med.rob$setting <- as.character(med.rob$setting)
forest <- rbind(
  s12,
  med.rob[, c("model", "setting", "estimate", "ci_low", "ci_high", "p_value", "star")]
)
forest$setting <- factor(forest$setting,
  levels = c("Main", "No patient covariate", "DROP_P1", "DROP_P2",
             "Nearby pairs (<=100 um)", "High contrast (0 vs P90)"))
forest$model <- factor(forest$model, levels = c("Tumor source", "Myeloid source"))

p.med <- ggplot(forest, aes(x = setting, y = estimate, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15,
                linewidth = 0.6, position = position_dodge(width = 0.4)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_text(aes(label = star), colour = "black", vjust = -0.9, size = 3.5,
            position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = MODEL_COL, name = "Source of PLAU") +
  coord_flip() +
  labs(title = "Robustness of mediation ACME",
       subtitle = "Main / no-patient / patient-LOSO / distance / contrast settings",
       x = "Setting", y = "ACME estimate with 95% CI") +
  sci_theme + theme(legend.position = "bottom")
sci_save(p.med, "RobustnessMediation", width = 10, height = 5)

# =============================================================================
# PART 6: 稳健性图件汇总
# 目的: 输出全部 5 张独立 PDF(无 merge)
# 数据: RobustnessLOSO / RobustnessEffectCI / RobustnessThreshold /
#       RobustnessCellChat / RobustnessMediation
# =============================================================================

cat("========== PART 6: Robustness plots ==========\n")

lab.map <- c(
  "sct_FBXL12" = "FBXL12 (all spots)",
  "sct_PLAU"   = "PLAU (all spots)",
  "sct_PLAUR"  = "PLAUR (all spots)",
  "sct_SPP1"   = "SPP1 (all spots)",
  "PLAUR_mac"  = "PLAUR (per macrophage)",
  "SPP1_mac"   = "SPP1 (per macrophage)"
)

# 参考: 全样本(不剔除)Cliff's delta, 用于 LOSO 图对照
full.rows <- list()
for (sig in CMP_SET) {
  for (pr in PAIRS) {
    a <- master[[sig]][master$compartment == pr[1]]; a <- a[!is.na(a)]
    b <- master[[sig]][master$compartment == pr[2]]; b <- b[!is.na(b)]
    full.rows[[length(full.rows) + 1]] <- data.frame(
      dropped = "all (full)", signature = sig,
      comparison = paste(pr[1], "vs", pr[2]),
      n1 = length(a), n2 = length(b),
      cliff_delta = round(CliffDelta(a, b), 4), stringsAsFactors = FALSE)
  }
}
full.tab <- do.call(rbind, full.rows)
loso.all <- rbind(loso.tab, full.tab)
loso.all$sig_label <- lab.map[loso.all$signature]

# (a) LOSO 稳定性图: 每比较分面, x=签名, 色=剔除样本, 含全样本对照
p.loso <- ggplot(loso.all, aes(x = sig_label, y = cliff_delta, colour = dropped)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 3, position = position_dodge(width = 0.55)) +
  facet_wrap(~ comparison, ncol = 1) +
  labs(title = "Cliff delta under leave-one-sample-out (LOSO)",
       subtitle = "Points = estimate when one sample dropped; 'all' = full-sample estimate",
       x = "Signature", y = "Cliff delta", colour = "Dropped sample") +
  sci_theme + theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 9),
                    legend.position = "bottom")
sci_save(p.loso, "RobustnessLOSO", width = 9, height = 9)

# (b) bootstrap 95% CI 效应量图
ci.tab$sig_label <- lab.map[ci.tab$signature]
p.ci <- ggplot(ci.tab, aes(x = sig_label, y = cliff_delta)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15,
                linewidth = 0.6) +
  geom_point(size = 3.2, shape = 21, fill = "#007AC2") +
  facet_wrap(~ comparison, ncol = 1) +
  labs(title = "Cliff delta with bootstrap 95 percent CI",
       subtitle = "Sample-clustered resampling (n=4); effect-size interval over single p-value",
       x = "Signature", y = "Cliff delta") +
  sci_theme + theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 9))
sci_save(p.ci, "RobustnessEffectCI", width = 8, height = 9)

# (c) 阈值敏感性图
thr.tab$q_cut <- factor(thr.tab$q_cut)
thr.tab$test_label <- c("H2_FBXL12_PLAUR" = "FBXL12 tumor vs PLAUR mac",
                        "H4_PLAUR_SPP1" = "PLAUR mac vs SPP1 mac"
                        )[thr.tab$test]
p.thr <- ggplot(thr.tab, aes(x = q_cut, y = ratio, colour = sample,
                             group = sample)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  facet_wrap(~ test_label, ncol = 2) +
  labs(title = "Threshold sensitivity of spatial colocalization",
       subtitle = "Observed to expected nearest distance ratio across high-expression cutoffs",
       x = "Quantile cutoff for high expression",
       y = "Observed to expected nearest distance ratio", colour = "Sample") +
  sci_theme
sci_save(p.thr, "RobustnessThreshold", width = 9, height = 6)

cat("13_PLAUR_SPP1_robustness.R done: RobustnessLOSO.pdf + RobustnessEffectCI.pdf + ",
    "RobustnessThreshold.pdf + RobustnessCellChat.pdf + RobustnessMediation.pdf\n",
    sep = "")
