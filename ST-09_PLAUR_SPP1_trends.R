# ===============================================================
# 09_PLAUR_SPP1_trends.R
# 功能: 核心假设链(FBXL12-high 癌 -> 分泌 PLAU -> 招募 PLAUR+ 巨噬 ->
#       上调 SPP1 -> PLAUR+SPP1+ 共定位巨噬)在 4 区室的协调上调证据
# 定位: 方案 Fig 2B(区室趋势) + 节点协调/效应端验证; 非"复刻原文"
# 可视化:
#   SignatureCompartmentTrends.pdf   - 6 签名区室 bootstrap 趋势(Fig 2B)
#   NodeCorrelationHeatmap.pdf       - 4 轴节点 Spearman 相关(协调证据)
#   ChainActivationComposite.pdf      - 链综合激活分提琴图(全程序同步)
#   MacrophageEffectorScatter.pdf     - PLAUR_mac vs SPP1_mac 散点(效应端)
#   NodeCoordinationDotplot.pdf       - 各对比 Cliff's delta 同向性(同步上调)
# 数据: 输入 results/06_PLAUR_SPP1_signature/table/signature_master.rds
#       输出 table/BootstrapTrends.csv / CompartmentComparisons.csv
#            table/NodeCorrelation.csv / ChainCoordination.csv
# 设计: bootstrap 采用"样本级聚类重抽样"(方案 v1.1), 每次有放回抽 4 样本
#       取各样本区内 spot 均值后再平均, 保留样本层相关性; 假设重述为链术语
#       (H1'=肿瘤区链节点协调上调 vs 癌旁; H6'=转移灶效应端较原发进一步升高
#        且全节点呈转移方向同步上调), 不再以"复刻原文 Fig 2B"为框架
# ===============================================================

# --- 包加载(统计 -> 数据 -> 可视化 -> 领域) -----------------------------
suppressMessages({
  library(ggplot2)     # 可视化
  library(reshape2)    # 长宽格式转换
  library(RColorBrewer)# 配色
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
  dev.off()
}

# --- 目录结构 -----------------------------------------------------------
BASE_DIR <- "X:/XX"
TAB06    <- file.path(BASE_DIR, "results/06_PLAUR_SPP1_signature/table")
FIG_DIR  <- file.path(BASE_DIR, "results/09_PLAUR_SPP1_trends/Figure")
TAB_DIR  <- file.path(BASE_DIR, "results/09_PLAUR_SPP1_trends/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 分析常量 -----------------------------------------------------------
COMP_ORDER <- c("Adjacent colon", "PT", "Adjacent liver", "LM")
# 核心链 4 轴节点(机制骨干) + 2 巨噬效应签名(效应端)
GENE_RAW  <- c("sct_FBXL12", "sct_PLAU", "sct_PLAUR", "sct_SPP1")
MAC_SIG   <- c("PLAUR_mac", "SPP1_mac")
CMP_SET   <- c(GENE_RAW, MAC_SIG)
# 短标签(规范17: 全英文)
SHORT_LAB <- c(
  "sct_FBXL12" = "FBXL12", "sct_PLAU" = "PLAU", "sct_PLAUR" = "PLAUR",
  "sct_SPP1"   = "SPP1", "PLAUR_mac" = "PLAUR_mac", "SPP1_mac" = "SPP1_mac"
)
N_BOOT     <- 1000
SEED       <- 42
COMP_FILL  <- c("Adjacent colon" = "#1F78B4", "PT" = "#E31A1C",
               "Adjacent liver" = "#A6CEE3", "LM" = "#FD8D3C")

# =============================================================================
# PART 1: 样本级聚类 bootstrap 估计各区室签名均值与 95% 分位区间
# 目的: 为每个签名在每个区室给出 bootstrap 均值与上下分位数(Fig 2B 主数据)
# 方法: 每次有放回抽 4 个样本, 取各样本区内 spot 均值后再平均
# 数据: BootstrapTrends.csv
# =============================================================================

cat("========== PART 1: Sample-clustered bootstrap ==========\n")

set.seed(SEED)
master <- readRDS(file.path(TAB06, "signature_master.rds"))
master <- master[!is.na(master$compartment), ]
master$patient <- substr(master$sample, 4, 5)
SAMPLE_ID <- unique(master$sample)

BootMeanClust <- function(sig, cp) {
  bs <- sample(SAMPLE_ID, length(SAMPLE_ID), replace = TRUE)
  vals <- sapply(bs, function(s)
    mean(master[[sig]][master$sample == s & master$compartment == cp],
         na.rm = TRUE))
  mean(vals, na.rm = TRUE)
}

boot.out <- list()
for (sig in CMP_SET) {
  for (cp in COMP_ORDER) {
    reps <- replicate(N_BOOT, BootMeanClust(sig, cp))
    boot.out[[length(boot.out) + 1]] <- data.frame(
      signature = sig, compartment = cp,
      mean_boot = mean(reps, na.rm = TRUE),
      lo = quantile(reps, 0.025, names = FALSE, na.rm = TRUE),
      hi = quantile(reps, 0.975, names = FALSE, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}
boot.tab <- do.call(rbind, boot.out)
boot.tab$compartment <- factor(boot.tab$compartment, levels = COMP_ORDER)
boot.tab$sig_short <- unname(SHORT_LAB[as.character(boot.tab$signature)])

obs.pt <- aggregate(master[, CMP_SET],
                    by = list(patient = master$patient,
                              compartment = master$compartment),
                    FUN = mean, na.rm = TRUE)

cat("Bootstrap 结果:\n")
print(boot.tab, digits = 4)
write.csv(boot.tab, file.path(TAB_DIR, "BootstrapTrends.csv"),
          row.names = FALSE)

# =============================================================================
# PART 2: 区室间比较, 假设重述为链术语(H1'/H6')
# 目的: 检验肿瘤区链节点协调上调(H1') 与转移灶效应端较原发进一步升高(H6')
# 方法: spot 级 Wilcoxon + Cliff's delta + BH 校正; 标注对比对应的链假设
# 数据: CompartmentComparisons.csv
# =============================================================================

cat("========== PART 2: Compartment comparisons (chain hypotheses) ==========\n")

# 对比对 -> 链假设标签(显式映射, 清晰优先)
PAIR_LAB <- list(
  list(c("PT", "Adjacent colon"), "H1: tumor vs adjacent"),
  list(c("LM", "Adjacent liver"), "H1: metastasis vs adjacent"),
  list(c("LM", "PT"), "H6: metastasis vs primary"),
  list(c("Adjacent liver", "Adjacent colon"), "control: liver vs colon adjacent")
)
comp.out <- list()
for (it in PAIR_LAB) {
  pr <- it[[1]]; hyp <- it[[2]]
  for (sig in CMP_SET) {
    a <- master[[sig]][master$compartment == pr[1]]; a <- a[!is.na(a)]
    b <- master[[sig]][master$compartment == pr[2]]; b <- b[!is.na(b)]
    wt <- wilcox.test(a, b, exact = FALSE)
    auc <- wt$statistic / (length(a) * length(b))
    comp.out[[length(comp.out) + 1]] <- data.frame(
      signature = sig, group1 = pr[1], group2 = pr[2], hypothesis = hyp,
      n1 = length(a), n2 = length(b),
      median1 = median(a), median2 = median(b),
      auc = round(auc, 4), cliff_delta = round(2 * auc - 1, 4),
      p_value = signif(wt$p.value, 4), stringsAsFactors = FALSE
    )
  }
}
comp.tab <- do.call(rbind, comp.out)
comp.tab$p_adj <- p.adjust(comp.tab$p_value, method = "BH")
comp.tab$sig_short <- unname(SHORT_LAB[as.character(comp.tab$signature)])

cat("区室比较(Wilcoxon + Cliff's delta + BH):\n")
print(comp.tab, digits = 4)
write.csv(comp.tab, file.path(TAB_DIR, "CompartmentComparisons.csv"),
          row.names = FALSE)

# =============================================================================
# PART 3: Fig 2B 趋势折线图(保留并放大)
# 目的: 带 bootstrap 95% 分位阴影的区室趋势, 叠加患者水平折线
# 数据: SignatureCompartmentTrends.pdf
# =============================================================================

cat("========== PART 3: Fig 2B trend plot ==========\n")

obs.long <- melt(obs.pt, id.vars = c("patient", "compartment"),
                 variable.name = "signature", value.name = "value")
obs.long$compartment <- factor(obs.long$compartment, levels = COMP_ORDER)
boot.tab$signature <- factor(boot.tab$signature, levels = CMP_SET)
obs.long$signature <- factor(obs.long$signature, levels = CMP_SET)

p.tr <- ggplot(boot.tab, aes(x = compartment, y = mean_boot, group = 1)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70", alpha = 0.5) +
  geom_line(linewidth = 1, colour = "black") +
  geom_point(size = 2.2, colour = "black") +
  geom_line(data = obs.long,
            aes(x = compartment, y = value, group = patient, colour = patient),
            linewidth = 0.5, linetype = "dashed", inherit.aes = FALSE) +
  facet_wrap(~ signature, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c(P1 = "#007AC2", P2 = "#FF6B00")) +
  labs(title = "Chain-node signature trends across compartments",
       x = "Compartment", y = "Signature score (SCT normalized)",
       colour = "Patient") +
  sci_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
        strip.text = element_text(size = 10))
sci_save(p.tr, "SignatureCompartmentTrends", width = 11, height = 9)

# =============================================================================
# PART 4: 4 轴节点 Spearman 相关热图(协调证据)
# 目的: 证明 FBXL12/PLAU/PLAUR/SPP1 在 spot 水平正相关 = 协调程序(支持链)
# 方法: 全 spot Spearman 相关矩阵; geom_tile 绘制(RdBu 发散色)
# 数据: NodeCorrelationHeatmap.pdf + NodeCorrelation.csv
# =============================================================================

cat("========== PART 4: Node correlation heatmap ==========\n")

corr <- cor(master[, GENE_RAW], use = "pairwise.complete.obs",
            method = "spearman")
corr.m <- melt(corr)
colnames(corr.m) <- c("node1", "node2", "rho")
corr.m$node1 <- factor(unname(SHORT_LAB[as.character(corr.m$node1)]),
                       levels = unname(SHORT_LAB[GENE_RAW]))
corr.m$node2 <- factor(unname(SHORT_LAB[as.character(corr.m$node2)]),
                       levels = unname(SHORT_LAB[GENE_RAW]))

p.corr <- ggplot(corr.m, aes(x = node1, y = node2, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 6,
            colour = ifelse(abs(corr.m$rho) > 0.5, "white", "black")) +
  scale_fill_gradientn(
    colours = rev(colorRampPalette(brewer.pal(9, "RdBu"))(100)),
    limits = c(-1, 1)) +
  labs(title = "Spearman correlation among chain axis nodes (all spots)",
       x = "Chain node", y = "Chain node", fill = "rho") +
  sci_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1))

sci_save(p.corr, "NodeCorrelationHeatmap", width = 7, height = 6)
write.csv(corr.m, file.path(TAB_DIR, "NodeCorrelation.csv"), row.names = FALSE)
cat("4 轴节点相关矩阵:\n"); print(round(corr, 3))

# =============================================================================
# PART 4b: 仅肿瘤 spot 内 4 轴节点 Spearman 相关(对 PART4 的补充, 防混淆)
# 目的: PART4 全 spot 相关被零胀(FBXL12/PLAU 多数 spot=0)与跨细胞稀释(上游
#       在肿瘤细胞/下游在髓系)削弱; 这里仅取 is_tumor spot 重算, 在同质肿瘤
#       微环境内评估节点协调。注意: 此图是 PART4 的补充, 并非替代; 两图含义
#       不同(全 spot vs 肿瘤内), 文中引用时需明确区分, 切勿混为一谈
# 方法: is_tumor & 4 节点无缺的 spot 子集; Spearman 相关矩阵 + geom_tile
# 数据: NodeCorrelationTumorHeatmap.pdf + NodeCorrelationTumor.csv
# =============================================================================

cat("========== PART 4b: Tumor-spot-only node correlation (supplement) ==========\n")

tum <- master[master$is_tumor & complete.cases(master[, GENE_RAW]), GENE_RAW]
corr_t <- cor(tum, use = "pairwise.complete.obs", method = "spearman")
corr_t.m <- melt(corr_t)
colnames(corr_t.m) <- c("node1", "node2", "rho")
corr_t.m$node1 <- factor(unname(SHORT_LAB[as.character(corr_t.m$node1)]),
                         levels = unname(SHORT_LAB[GENE_RAW]))
corr_t.m$node2 <- factor(unname(SHORT_LAB[as.character(corr_t.m$node2)]),
                         levels = unname(SHORT_LAB[GENE_RAW]))

p.corr_t <- ggplot(corr_t.m, aes(x = node1, y = node2, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)),
            colour = ifelse(abs(corr_t.m$rho) > 0.5, "white", "black"),
            size = 6) +
  scale_fill_gradientn(
    colours = rev(colorRampPalette(brewer.pal(9, "RdBu"))(100)),
    limits = c(-1, 1)) +
  labs(title = "SUPPLEMENT: chain-node correlation within TUMOR spots only",
       subtitle = "Complement to all-spot Fig; excludes normal-adjacent dilution",
       x = "Chain node", y = "Chain node", fill = "rho") +
  sci_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1))

sci_save(p.corr_t, "NodeCorrelationTumorHeatmap", width = 7, height = 6)
write.csv(corr_t.m, file.path(TAB_DIR, "NodeCorrelationTumor.csv"),
          row.names = FALSE)
cat(sprintf("肿瘤 spot 内 4 轴节点相关(n=%d):\n", nrow(tum)))
print(round(corr_t, 3))

# =============================================================================
# PART 5: 链综合激活分提琴图(全程序同步上调证据)
# 目的: 将 4 轴节点行内 Z 化后求和得综合激活分, 显示转移灶全程序最高
# 方法: 每 spot 对 4 节点 Z 化求和; 提琴 + 箱线按区室
# 数据: ChainActivationComposite.pdf + ChainCoordination.csv
# =============================================================================

cat("========== PART 5: Composite chain activation violin ==========\n")

sub <- master[complete.cases(master[, GENE_RAW]), c(GENE_RAW, "compartment")]
zc <- as.data.frame(scale(sub[, GENE_RAW]))
sub$chain_score <- rowSums(zc)
sub$compartment <- factor(sub$compartment, levels = COMP_ORDER)

p.comp <- ggplot(sub, aes(x = compartment, y = chain_score, fill = compartment)) +
  geom_violin(alpha = 0.6, colour = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.3) +
  scale_fill_manual(values = COMP_FILL) +
  labs(title = "Composite chain activation score by compartment",
       x = "Compartment", y = "Sum of 4 node Z-scores") +
  sci_theme + theme(legend.position = "none",
                    axis.text.x = element_text(angle = 20, hjust = 1, size = 9))

sci_save(p.comp, "ChainActivationComposite", width = 8, height = 6)

# 综合分统计 + 节点同向性汇总(协调证据)
comp_score_sum <- aggregate(chain_score ~ compartment, sub, function(x)
  c(median = median(x), mean = mean(x), n = length(x)))
comp_score_sum <- as.data.frame(comp_score_sum)
comp_score_sum$median <- round(comp_score_sum$chain_score[, "median"], 3)
comp_score_sum$mean <- round(comp_score_sum$chain_score[, "mean"], 3)
comp_score_sum$n <- comp_score_sum$chain_score[, "n"]
comp_score_sum$chain_score <- NULL

# 节点方向同向性: 在 3 个方向性对比中, 各节点 cliff_delta 为正的比例
dir.tab <- comp.tab[comp.tab$hypothesis != "control: liver vs colon adjacent", ]
node_coord <- aggregate(cliff_delta ~ sig_short, dir.tab, function(x)
  c(frac_pos = mean(x > 0), mean_delta = mean(x)))
node_coord <- as.data.frame(node_coord)
node_coord$frac_pos <- round(node_coord$cliff_delta[, "frac_pos"], 3)
node_coord$mean_delta <- round(node_coord$cliff_delta[, "mean_delta"], 3)
node_coord$cliff_delta <- NULL

write.csv(comp_score_sum, file.path(TAB_DIR, "ChainCoordination.csv"),
          row.names = FALSE)
write.csv(node_coord, file.path(TAB_DIR, "NodeDirectionSummary.csv"),
          row.names = FALSE)
cat("综合激活分(区室中位数):\n"); print(comp_score_sum)
cat("节点方向同向性(3 方向对比中正 delta 比例):\n"); print(node_coord)

# =============================================================================
# PART 6: 效应端散点 PLAUR_mac vs SPP1_mac(巨噬细胞水平)
# 目的: 直接显示效应端双签名在髓系 spot 内的共表达, 并富集于转移灶
# 方法: is_macro spot 散点, 颜色=区室, 叠加整体 lm 趋势
# 数据: MacrophageEffectorScatter.pdf
# =============================================================================

cat("========== PART 6: Macrophage effector-end scatter ==========\n")

dm <- master[master$is_macro & !is.na(master$PLAUR_mac) &
             !is.na(master$SPP1_mac), ]
dm$compartment <- factor(dm$compartment, levels = COMP_ORDER)

p.sc <- ggplot(dm, aes(x = PLAUR_mac, y = SPP1_mac, colour = compartment)) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = COMP_FILL) +
  facet_wrap(~ sample, ncol = 2) +
  labs(title = "Effector-end: PLAUR_mac vs SPP1_mac in myeloid spots",
       x = "PLAUR_mac (per macrophage)", y = "SPP1_mac (per macrophage)",
       colour = "Compartment") +
  sci_theme + theme(strip.text = element_text(size = 9))

sci_save(p.sc, "MacrophageEffectorScatter", width = 10, height = 7)

# =============================================================================
# PART 7: 节点同向性点图(同步上调证据)
# 目的: 各签名在 3 个方向性对比的 Cliff's delta 同向为正 = 链同步上调
# 方法: x=签名(链序), y=cliff_delta, 颜色=对比; 0 参考线
# 数据: NodeCoordinationDotplot.pdf
# =============================================================================

cat("========== PART 7: Node coordination dot plot ==========\n")

dir.tab$sig_short <- factor(dir.tab$sig_short,
                           levels = unname(SHORT_LAB[CMP_SET]))
dir.tab$contrast <- factor(paste(dir.tab$group1, "vs", dir.tab$group2),
  levels = c("PT vs Adjacent colon", "LM vs Adjacent liver",
             "LM vs PT"))

p.dot <- ggplot(dir.tab, aes(x = sig_short, y = cliff_delta,
                             colour = contrast, group = contrast)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_line(position = position_dodge(width = 0.5), linewidth = 0.6) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Directional coordination: Cliff's delta across contrasts",
       x = "Chain node", y = "Cliff's delta (group1 > group2)") +
  sci_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9))

sci_save(p.dot, "NodeCoordinationDotplot", width = 9, height = 6)

cat("09_PLAUR_SPP1_trends.R done:\n")
cat("  - SignatureCompartmentTrends.pdf (Fig 2B, enlarged)\n")
cat("  - NodeCorrelationHeatmap.pdf (ALL spots)\n")
cat("  - NodeCorrelationTumorHeatmap.pdf (SUPPLEMENT: tumor spots only)\n")
cat("  - ChainActivationComposite.pdf\n")
cat("  - MacrophageEffectorScatter.pdf\n")
cat("  - NodeCoordinationDotplot.pdf\n")
