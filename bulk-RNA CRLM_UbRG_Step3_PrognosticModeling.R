##############################################################
##  CRLM 泛素相关基因研究 - Step 3: 预后模型构建与验证
##
##  【分析逻辑】
##   基于Step2筛出的独立预后Hub Genes(数量由数据动态决定),
##   在GSE39582建立风险评分模型:
##   1. 风险评分公式(基于LASSO系数)
##   2. 训练集KM曲线 + 时间依赖ROC
##   3. 风险分布图(风险评分 + 生存状态)
##   4. 独立预后因素分析(风险评分 + 临床变量多因素Cox)
##   5. Nomogram构建
##   6. 校准曲线
##   7. 外部验证集(GSE159216)KM + ROC
##   8. 分层分析(TNM分期亚组)
##   9. 单基因KM + 箱线图
##   10. Hub基因功能验证
##############################################################

setwd("X:/CRLM")
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(scales)
library(timeROC)
library(patchwork)
library(rms)
library(pheatmap)
library(RColorBrewer)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)

# ================================================================
# Part 0: 加载数据
# ================================================================
cat("=== 加载数据 ===\n")
expr39582_norm <- readRDS("data/processed/expr39582_norm.rds")
clin_39582     <- readRDS("data/processed/clin_39582.rds")
group_39582    <- readRDS("data/processed/group_39582.rds")
final_genes    <- readRDS("data/processed/final_genes.rds")
final_coef_df  <- readRDS("data/processed/final_coef_df.rds")
multi_res      <- readRDS("data/processed/multi_cox_res.rds")

# 验证集数据
expr159216_norm <- readRDS("data/processed/expr159216_norm.rds")
pheno159216     <- readRDS("data/processed/pheno159216.rds")

cat("最终预后基因:", paste(final_genes, collapse = ", "), "\n")
cat("最终预后基因数:", length(final_genes), "\n")

# ================================================================
# Part 1: 计算风险评分 (GSE39582 训练集)
# ================================================================
cat("\n=== Part1: 计算风险评分 ===\n")

# [设计说明]
# 风险评分使用Hub Genes的LASSO系数计算(与Step2一致)
# Step2多因素Cox筛出的独立预后基因 = Hub Genes = 风险评分基因
multicox_coef_df <- readRDS("data/processed/multicox_coef.rds")
gene_coef_all    <- setNames(multicox_coef_df$coef, multicox_coef_df$gene)

# 确保final_genes与multicox_coef一致
final_genes <- intersect(final_genes, names(gene_coef_all))
cat("风险评分基因数:", length(gene_coef_all), "\n")
cat("独立预后展示基因:", paste(final_genes, collapse = ", "), "\n")

# 风险评分公式: RiskScore = sum(beta_i * expr_i), beta来自LASSO
lasso_genes_use <- names(gene_coef_all)
lasso_genes_in_expr  <- intersect(lasso_genes_use, rownames(expr39582_norm))
gene_coef_used       <- gene_coef_all[lasso_genes_in_expr]
final_genes_in_expr  <- intersect(final_genes, rownames(expr39582_norm))

cat("在表达矩阵中匹配的风险评分基因:", length(lasso_genes_in_expr), "\n")

risk_mat   <- expr39582_norm[lasso_genes_in_expr, , drop = FALSE]
risk_score <- as.numeric(colSums(sweep(risk_mat, 1, gene_coef_used, FUN = "*")))
names(risk_score) <- colnames(expr39582_norm)

# 构建生存数据框(仅Tumor样本)
os_event  <- suppressWarnings(as.numeric(clin_39582$os_event))
os_delay  <- suppressWarnings(as.numeric(clin_39582$os_delay))
rfs_event <- suppressWarnings(as.numeric(clin_39582$rfs_event))
rfs_delay <- suppressWarnings(as.numeric(clin_39582$rfs_delay))

# 过滤有效生存数据的Tumor样本
valid_os  <- group_39582 == "Tumor" & !is.na(os_event)  & !is.na(os_delay)  & os_delay > 0
valid_rfs <- group_39582 == "Tumor" & !is.na(rfs_event) & !is.na(rfs_delay) & rfs_delay > 0

surv_df <- data.frame(
  sample_id  = colnames(expr39582_norm)[valid_os],
  risk_score = risk_score[valid_os],
  OS_event   = os_event[valid_os],
  OS_time    = os_delay[valid_os],
  stringsAsFactors = FALSE
)

# 添加临床信息
tnm_stage <- as.character(clin_39582$tnm_stage)[valid_os]
mol_sub   <- as.character(clin_39582$molecular_subtype)[valid_os]
surv_df$TNM_stage <- tnm_stage
surv_df$Subtype   <- mol_sub

# 中位数分组
tnm_stage_for_filter <- as.character(clin_39582$tnm_stage)[valid_os]
stage0_exclude <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", tnm_stage_for_filter))) >= 1
median_rs <- median(surv_df$risk_score[stage0_exclude], na.rm = TRUE)
surv_df$risk_group <- ifelse(surv_df$risk_score >= median_rs, "High", "Low")

cat("训练集总有效样本:", nrow(surv_df), "\n")
cat("中位风险评分:", round(median_rs, 4), "\n")
cat("High风险:", sum(surv_df$risk_group == "High"), "\n")
cat("Low风险:", sum(surv_df$risk_group == "Low"), "\n")

# 保存
write.csv(surv_df, "results/tables/Table8_RiskScore_TrainingSet.csv", row.names = FALSE)
saveRDS(surv_df,   "data/processed/surv_df_train.rds")

# ================================================================
# Part 2: KM曲线 (训练集)
# ================================================================
cat("\n=== Part2: KM曲线 ===\n")

fit_km <- survfit(Surv(OS_time, OS_event) ~ risk_group, data = surv_df)
pval_km <- surv_pvalue(fit_km, data = surv_df)$pval

fig2a <- ggsurvplot(
  fit_km,
  data      = surv_df,
  pval      = TRUE,
  pval.size = 4,
  risk.table = TRUE,
  risk.table.height = 0.22,
  risk.table.col = "strata",
  palette   = c("#E64B35", "#4DBBD5"),
  title     = "Kaplan-Meier: Risk Score Groups (Training Set)",
  subtitle  = paste0("GSE39582 (n=", nrow(surv_df), ")"),
  xlab      = "Time (months)",
  ylab      = "Overall Survival Probability",
  legend.labs = c("High Risk", "Low Risk"),
  conf.int  = FALSE,
  surv.median.line = "hv",
  ggtheme   = theme_classic(base_size = 12)
)

# ggsurvplot 返回列表对象,直接用 ggsave 会丢失 risk table
# 正确方式:pdf() + print(newpage=FALSE) 或用 survminer:::.save_ggsurvplot
pdf("results/figures/PM1_KM_TrainingSet.pdf", width = 7, height = 7)
print(fig2a, newpage = FALSE)
dev.off()
cat("  -> PM1 完成\n")

# ================================================================
# Part 3: 时间依赖ROC曲线
# ================================================================
cat("=== Part3: 时间依赖ROC ===\n")
roc_res <- tryCatch({
  timeROC::timeROC(
    T        = surv_df$OS_time,
    delta    = surv_df$OS_event,
    marker   = surv_df$risk_score,
    cause    = 1,
    times    = c(12, 36, 60),
    iid      = FALSE
  )
}, error = function(e) {
  cat("timeROC error:", e$message, "\n")
  NULL
})

if (!is.null(roc_res)) {
  auc12 <- roc_res$AUC[1]; auc36 <- roc_res$AUC[2]; auc60 <- roc_res$AUC[3]
  cat("AUC 1年:", round(auc12, 3), "\n")
  cat("AUC 3年:", round(auc36, 3), "\n")
  cat("AUC 5年:", round(auc60, 3), "\n")

  # 构建ROC plot数据
  roc_plot_df <- rbind(
    data.frame(sens=roc_res$TP[,1], spec=1-roc_res$FP[,1], time="1 year (AUC=%.3f)" %>% sprintf(auc12)),
    data.frame(sens=roc_res$TP[,2], spec=1-roc_res$FP[,2], time=sprintf("3 year (AUC=%.3f)", auc36)),
    data.frame(sens=roc_res$TP[,3], spec=1-roc_res$FP[,3], time=sprintf("5 year (AUC=%.3f)", auc60))
  )

  fig2b <- ggplot(roc_plot_df, aes(x = 1 - spec, y = sens, color = time)) +
    geom_line(linewidth = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    scale_color_manual(values = c("#E64B35", "#0072B2", "#009E73")) +
    labs(title = "Time-dependent ROC (Training Set)",
         x = "1 - Specificity", y = "Sensitivity",
         color = "") +
    theme_classic(base_size = 12) +
    theme(legend.position = c(0.65, 0.25),
          legend.background = element_rect(fill = "white", color = "grey80"))

  auc_df <- data.frame(time = c("1y", "3y", "5y"), AUC = c(auc12, auc36, auc60))
  write.csv(auc_df, "results/tables/Table9_AUC_Training.csv", row.names = FALSE)
} else {
  # 备用方案: 使用survivalROC
  cat("使用备用ROC方法...\n")
  fig2b <- ggplot(data.frame(x=0:1, y=0:1), aes(x,y)) +
    geom_line(linetype="dashed", color="grey") +
    labs(title="Time-dependent ROC (unavailable)", x="1-Specificity", y="Sensitivity") +
    theme_classic()
}

pdf("results/figures/Supp_S2O_ROC_TrainingSet.pdf", width = 6, height = 6)
print(fig2b)
dev.off()
cat("  -> Supp_S2O 完成\n")

# ================================================================
# Part 4: 风险分布图(三合一)
# ================================================================
cat("=== Part4: 风险分布图 ===\n")
surv_sorted <- surv_df %>% arrange(risk_score) %>% mutate(idx = seq_len(nrow(.)))

# 最终基因表达热图(用于Fig2C的第三部分,但当前仅保留前两部分)

# 计算风险分组统计
n_high <- sum(surv_sorted$risk_group == "High")
n_low <- sum(surv_sorted$risk_group == "Low")
n_dead <- sum(surv_sorted$OS_event == 1)
n_alive <- sum(surv_sorted$OS_event == 0)

# 只保留Risk Score Distribution和Survival Status两部分
# 使用ggplot2绘制,上下排列填充整个画布

# 1. 风险评分分布图
p1 <- ggplot(surv_sorted, aes(x = idx, y = risk_score, color = risk_group)) +
  geom_point(size = 1.2, alpha = 0.8) +
  geom_hline(yintercept = median_rs, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  scale_color_manual(values = c(High = "#E64B35", Low = "#4DBBD5"),
                     labels = c(High = "High Risk", Low = "Low Risk")) +
  labs(title = "A  Risk Score Distribution",
       subtitle = sprintf("High Risk n=%d  |  Low Risk n=%d", n_high, n_low),
       x = "Patient Index (sorted by risk score)",
       y = "Risk Score") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dotted")
  )

# 2. 生存状态图
p2 <- ggplot(surv_sorted, aes(x = idx, y = OS_time, color = factor(OS_event))) +
  geom_point(size = 1.2, alpha = 0.8) +
  scale_color_manual(values = c("1" = "#E64B35", "0" = "#4DBBD5"),
                     labels = c("0" = "Alive", "1" = "Dead")) +
  labs(title = "B  Survival Status",
       subtitle = sprintf("Dead n=%d  |  Alive n=%d", n_dead, n_alive),
       x = "Patient Index (sorted by risk score)",
       y = "Overall Survival Time (months)",
       color = "Status") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dotted")
  )

# 3. Hub基因表达热图(按风险评分排序,与前两图样本顺序一致)
library(pheatmap)
library(scales)

# 准备热图数据:hub基因 + 按风险评分排序的样本
hm_expr <- expr39582_norm[final_genes_in_expr, surv_sorted$sample_id]
hm_expr_z <- t(scale(t(hm_expr)))  # 行标准化(基因内Z-score)
rownames(hm_expr_z) <- final_genes_in_expr

# 热图annotation bar - 添加行名以匹配列顺序
ha_samples <- data.frame(
  `Risk Group` = factor(surv_sorted$risk_group, levels = c("Low", "High")),
  `OS Status`  = factor(ifelse(surv_sorted$OS_event == 1, "Dead", "Alive"),
                        levels = c("Alive", "Dead")),
  row.names = surv_sorted$sample_id  # 关键修复:添加行名匹配热图列顺序
)
ann_colors <- list(
  `Risk Group` = c(High = "#E64B35", Low = "#4DBBD5"),
  `OS Status`  = c(Alive = "#4DBBD5", Dead = "#E64B35")
)

# 列标题(样本太多时设为空)
n_samp <- ncol(hm_expr_z)
show_colnames <- ifelse(n_samp <= 100, TRUE, FALSE)

# 热图对象(不直接打印,由print触发)
fig2c_heatmap <- pheatmap(
  hm_expr_z,
  annotation = ha_samples,  # 直接使用已构造好的annotation(含行名)
  annotation_colors = ann_colors,
  color             = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
  scale             = "none",
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  show_colnames     = show_colnames,
  fontsize_col      = ifelse(n_samp <= 100, 8, 3),
  fontsize_row      = 10,
  annotation_legend = TRUE,
  legend            = TRUE,
  main              = NA,
  silent            = TRUE
)

# NOTE: pheatmap cannot be combined with ggplot via patchwork.
# pheatmap is saved separately, then risk score + survival are combined.

# 1. Save heatmap separately
pdf("results/figures/PM3_Heatmap_only.pdf", width = 14, height = 10)
print(fig2c_heatmap)
dev.off()
cat("  -> PM3_Heatmap_only.pdf saved\n")


# 2. Save risk score + survival (without heatmap)
library(cowplot)
fig2c_combined <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1, 1),
                             align = "v", axis = "l")

# 保存
pdf("results/figures/PM4_RiskDistribution.pdf", width = 8, height = 10)
print(fig2c_combined)
dev.off()
cat("  -> PM4 done (risk+survival, heatmap saved separately)\n")

# Part 5: 独立预后因素分析
# ================================================================
cat("=== Part5: 独立预后因素分析 ===\n")

# 准备包含临床变量的数据框
# TNM_stage编码: 0=Stage0, 1=StageI, 2=StageII, 3=StageIII, 4=StageIV
# 排除Stage 0(4例)和不明确的值,仅保留I-IV期用于分析
indep_df <- surv_df %>%
  mutate(
    TNM_numeric = suppressWarnings(as.numeric(gsub("[^0-9.]", "", TNM_stage))),
    # 【修复P0】用命名subetting（valid_os是命名逻辑向量，clin_39582字段是命名向量）
    # 避免which()返回数字索引与样本ID名称不匹配的问题
    Age         = suppressWarnings(as.numeric(clin_39582$age[valid_os])),
    Sex         = clin_39582$sex[valid_os],
    TNM_binary  = case_when(
      TNM_numeric %in% 1:2 ~ "I/II",
      TNM_numeric %in% 3:4 ~ "III/IV",
      TRUE ~ NA_character_   # Stage 0和NA排除
    )
  )

# 单因素(风险评分)
uni_indep <- coxph(Surv(OS_time, OS_event) ~ risk_score, data = indep_df)
uni_sum   <- summary(uni_indep)
cat("单因素 - 风险评分: HR=", round(uni_sum$coefficients["risk_score", "exp(coef)"], 3),
    "p=", round(uni_sum$coefficients["risk_score", "Pr(>|z|)"], 4), "\n")

# 多因素 (风险评分 + TNM分期 + Age)
# 排除TNM_stage=0(Stage 0, 4例)和NA, 仅保留I-IV期
indep_df_valid <- indep_df %>% filter(!is.na(TNM_numeric), TNM_numeric >= 1, TNM_numeric <= 4)
cat("独立预后分析有效样本:", nrow(indep_df_valid), "\n")

if (nrow(indep_df_valid) > 50) {
  multi_indep  <- coxph(Surv(OS_time, OS_event) ~ risk_score + TNM_numeric + Age,
                        data = indep_df_valid)
  multi_indep_sum <- summary(multi_indep)

  indep_res <- data.frame(
    Variable = c("Risk Score", "TNM Stage", "Age"),
    HR       = multi_indep_sum$coefficients[, "exp(coef)"],     # conf.int中已是HR尺度
    CI_low   = multi_indep_sum$conf.int[, "lower .95"],         # HR尺度, 勿再exp()
    CI_high  = multi_indep_sum$conf.int[, "upper .95"],         # HR尺度, 勿再exp()
    p_value  = multi_indep_sum$coefficients[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  print(indep_res)
  write.csv(indep_res, "results/tables/Table10_Independent_Prognosis.csv", row.names = FALSE)

  # 森林图 - 与PM19完全一致的双面板风格
  indep_res <- indep_res %>%
    mutate(
      Variable = factor(Variable, levels = rev(Variable)),
      risk = ifelse(HR > 1, "High risk", "Low risk"),
      hr_label = sprintf("%.2f", HR),
      ci_label = sprintf("(%.2f-%.2f)", CI_low, CI_high),
      p_label = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
    )
  log_lo <- floor(log10(min(indep_res$CI_low, na.rm = TRUE) * 0.8))
  log_hi <- ceiling(log10(max(indep_res$CI_high, na.rm = TRUE) * 1.3))
  log_breaks <- 10^seq(log_lo, log_hi, by = 1)
  log_breaks <- log_breaks[log_breaks >= min(indep_res$CI_low) * 0.8 &
                             log_breaks <= max(indep_res$CI_high) * 1.4]
  if (length(log_breaks) < 3) log_breaks <- scales::log_breaks()(range(c(indep_res$CI_low, indep_res$CI_high)))
  # 森林图主体
  pm5_forest <- ggplot(indep_res, aes(x = HR, y = Variable)) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high, color = risk),
                  width = 0.3, linewidth = 1, orientation = "y") +
    geom_point(aes(color = risk), size = 4) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.8) +
    scale_color_manual(values = c("High risk" = "#E64B35", "Low risk" = "#4DBBD5")) +
    scale_x_log10(breaks = log_breaks, labels = sprintf("%.2g", log_breaks)) +
    labs(x = "Hazard Ratio (95% CI)", y = "") +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 11),
      panel.grid.major.x = element_line(color = "grey90", linetype = "dotted")
    )
  # 右侧表格 - 与PM19完全一致
  table_pm5 <- indep_res %>%
    select(Variable, hr_label, ci_label, p_label) %>%
    pivot_longer(cols = c(hr_label, ci_label, p_label), names_to = "col", values_to = "val") %>%
    mutate(x_pos = case_when(col == "hr_label" ~ 1, col == "ci_label" ~ 2, col == "p_label" ~ 3))
  table_ann_pm5 <- ggplot(table_pm5, aes(x = x_pos, y = Variable, label = val)) +
    geom_text(size = 3.5, vjust = 0.5) +
    scale_x_continuous(limits = c(0.5, 3.5), breaks = 1:3,
                       labels = c("HR", "95% CI", "P"), expand = c(0, 0)) +
    theme_void(base_size = 11) +
    theme(axis.text.x = element_text(face = "bold", size = 10))
  pm5_combined <- pm5_forest + table_ann_pm5 + plot_layout(widths = c(2, 1)) +
    plot_annotation(
      title = "Independent Prognostic Factors (OS)",
      subtitle = "Multivariable Cox regression for overall survival",
      theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
    )
  pdf("results/figures/PM5_Forest_Independent.pdf", width = 9, height = 5)
  print(pm5_combined)
  dev.off()
  cat("  -> PM5 完成\n")
} else {
  cat("TNM数据不足,跳过多因素独立预后分析\n")
}

# ================================================================
# Part 6: Nomogram
# ================================================================
cat("=== Part6: Nomogram ===\n")
tryCatch({
  library(rms)

  # 准备临床变量
  surv_df$Age       <- suppressWarnings(as.numeric(clin_39582$age[valid_os]))
  surv_df$TNM_num   <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", surv_df$TNM_stage)))

  # Nomogram只用risk_score + TNM_num + Age,datadist只含这些列
  nomo_data <- surv_df[!is.na(surv_df$TNM_num) & surv_df$TNM_num >= 1,
                       c("OS_time", "OS_event", "risk_score", "TNM_num", "Age")]
  nomo_data <- nomo_data[complete.cases(nomo_data), ]
  dd <- datadist(nomo_data)
  options(datadist = "dd")

  # Weibull参数生存模型构建Nomogram(含RiskScore + TNM + Age)
  # 注:GSE39582未收录Sex,不纳入Nomogram
  nomo_fit <- psm(Surv(OS_time, OS_event) ~ risk_score + TNM_num + Age,
                  data = nomo_data,
                  dist = "weibull", x = TRUE, y = TRUE)
  surv_fun <- Survival(nomo_fit)

  nom <- nomogram(nomo_fit,
                  fun       = list(function(x) surv_fun(36, x), function(x) surv_fun(60, x)),
                  funlabel  = c("3-year Survival", "5-year Survival"),
                  lp        = FALSE)

  pdf("results/figures/PM6_Nomogram.pdf", width = 10, height = 6)
  plot(nom, xfrac = 0.35, cex.main = 1.2, main = "Nomogram for OS Prediction")
  dev.off()
  cat("  -> PM6 完成\n")
}, error = function(e) {
  cat("Nomogram error:", e$message, "\n")
  # 绘制风险评分的生存校准图
  pdf("results/figures/PM6_Nomogram.pdf", width = 6, height = 4)
  plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", main="Nomogram (unavailable)")
  text(0.5, 0.5, "Nomogram generation failed\nSee Risk Score KM instead", cex=1.2)
  dev.off()
})

# ================================================================
# Part 7: 校准曲线
# ================================================================
cat("=== Part7: 校准曲线 ===\n")
tryCatch({
  library(rms)
  cal_df <- surv_df[, c("OS_time", "OS_event", "risk_score")]
  cal_df <- cal_df[complete.cases(cal_df), ]
  dd_cal <- datadist(cal_df)
  options(datadist = "dd_cal")

  cal_fit <- psm(Surv(OS_time, OS_event) ~ risk_score, data = cal_df, dist = "weibull", x = TRUE, y = TRUE)
  cal3    <- calibrate(cal_fit, u = 36, B = 100)
  cal5    <- calibrate(cal_fit, u = 60, B = 100)

  pdf("results/figures/PM7_Calibration.pdf", width = 6, height = 6)
  plot(cal3, subtitles = FALSE, main = "Calibration Curve (3-year OS)",
       lty = 1, col = "#E64B35", lwd = 2)
  plot(cal5, subtitles = FALSE, add = TRUE, lty = 2, col = "#0072B2", lwd = 2)
  legend("bottomright", legend = c("3-year OS", "5-year OS"),
         col = c("#E64B35", "#0072B2"), lty = c(1, 2), lwd = 2, bty = "n")
  dev.off()
  cat("  -> PM7 完成\n")
}, error = function(e) {
  cat("Calibration error:", e$message, "\n")
  pdf("results/figures/PM7_Calibration.pdf", width = 6, height = 5)
  plot(0, type="n", main="Calibration (error)", xlab="", ylab="")
  text(0, 0, e$message, cex=0.8)
  dev.off()
})

# ================================================================
# Part 8: 分层分析 (TNM分期亚组KM)
# ================================================================
cat("=== Part8: 分层KM ===\n")

if (sum(!is.na(suppressWarnings(as.numeric(surv_df$TNM_stage)))) > 50) {
  strat_df <- surv_df %>%
    mutate(TNM_num = suppressWarnings(as.numeric(gsub("[^0-9.]", "", TNM_stage))),
           TNM_grp = case_when(
             TNM_num %in% 1:2 ~ "Stage I/II",
             TNM_num %in% 3:4 ~ "Stage III/IV",
             TRUE ~ NA_character_   # 排除Stage 0
           ))

  groups_strat <- unique(strat_df$TNM_grp[!is.na(strat_df$TNM_grp)])

  # 为每个分层生成ggsurvplot对象,然后用patchwork拼接
  km_plots <- list()

  for (g in groups_strat) {
    sub_df <- strat_df %>% filter(TNM_grp == g, !is.na(TNM_num))
    if (sum(sub_df$risk_group == "High") < 5 || sum(sub_df$risk_group == "Low") < 5) next

    fit_s <- survfit(Surv(OS_time, OS_event) ~ risk_group, data = sub_df)
    pv    <- survdiff(Surv(OS_time, OS_event) ~ risk_group, data = sub_df)
    pv_num <- 1 - pchisq(pv$chisq, df = 1)

    fig_strat <- ggsurvplot(
      fit_s, data = sub_df,
      pval = TRUE, pval.size = 4,
      risk.table = TRUE, risk.table.height = 0.22, risk.table.col = "strata",
      palette = c("#E64B35", "#4DBBD5"),
      title = paste0("Kaplan-Meier: Risk Score Groups"),
      subtitle = paste0(g, " (n=", nrow(sub_df), ", p=", signif(pv_num, 3), ")"),
      xlab = "Time (months)",
      ylab = "Overall Survival Probability",
      legend.labs = c("High Risk", "Low Risk"),
      conf.int = FALSE,
      surv.median.line = "hv",
      ggtheme = theme_classic(base_size = 12)
    )

    # 保存单独的分层数据
    km_plots[[g]] <- fig_strat
  }

  # 使用 patchwork 拼接多个 ggsurvplot
  if (length(km_plots) > 0) {
    # 每个ggsurvplot需要用arrange_ggsurvplots或单独保存
    # 简化方案:保存为单独文件或使用gridExtra
    library(gridExtra)

    # 提取plot组件并用patchwork拼接
    plot_list <- lapply(km_plots, function(x) x$plot)

    # 使用patchwork拼接
    combined_plot <- wrap_plots(plot_list, ncol = length(plot_list))

    pdf("results/figures/PM8_Stratified_KM.pdf", width = 7 * length(plot_list), height = 7)
    print(combined_plot)
    dev.off()
    cat("  -> PM8 完成\n")
  }
}

# ================================================================
# Part 9: 外部验证 (GSE159216)
# ================================================================
cat("\n=== Part9: 外部验证集 (GSE159216) ===\n")

# 提取生存信息 (GSE159216两个独立列)
status_col <- "60 months overall survival, status:ch1"
time_col   <- "60 months overall survival, time:ch1"

if (status_col %in% colnames(pheno159216) && time_col %in% colnames(pheno159216)) {
  os_event_all <- suppressWarnings(as.numeric(pheno159216[[status_col]]))
  os_time_all  <- suppressWarnings(as.numeric(pheno159216[[time_col]]))

  valid_sv   <- !is.na(os_event_all) & !is.na(os_time_all) & os_time_all > 0
  os_event_v <- os_event_all[valid_sv]
  os_time_v  <- os_time_all[valid_sv]
  cat("有效生存数据:", sum(valid_sv), "/", length(valid_sv), "\n")

  samp_v <- colnames(expr159216_norm)
  valid_samp <- rownames(pheno159216)[valid_sv]
  matched    <- intersect(samp_v, valid_samp)
  cat("匹配样本:", length(matched), "\n")

  # 有效生存数据(按matched顺序)
  os_event_v2 <- os_event_v[match(matched, valid_samp)]
  os_time_v2  <- os_time_v[match(matched, valid_samp)]
  valid_v <- !is.na(os_event_v2) & !is.na(os_time_v2) & os_time_v2 > 0

  # 计算风险评分 (验证集) - 用全部LASSO基因(训练集与验证集共有的)
  # 注意: 验证集可能缺少部分LASSO基因,仅用共有基因计算
  genes_v       <- intersect(lasso_genes_in_expr, rownames(expr159216_norm))
  gene_coef_v   <- gene_coef_used[genes_v]
  cat("验证集中可用LASSO基因:", length(genes_v), "/", length(lasso_genes_in_expr), "\n")
  if (length(genes_v) < length(lasso_genes_in_expr)) {
    cat("  缺失基因:", setdiff(lasso_genes_in_expr, rownames(expr159216_norm)), "\n")
  }

  rs_mat_v <- expr159216_norm[genes_v, matched, drop = FALSE]
  rs_v     <- as.numeric(colSums(sweep(rs_mat_v, 1, gene_coef_v, FUN = "*")))

  # 用验证集自身中位数分组(训练集与验证集分布不同,不宜用训练集中位数切分)
  median_rs_v <- median(rs_v, na.rm = TRUE)
  surv_df_v  <- data.frame(
    sample_id  = matched,
    risk_score = rs_v,
    OS_event   = os_event_v2,
    OS_time    = os_time_v2,
    stringsAsFactors = FALSE
  ) %>% filter(valid_v)
  surv_df_v$risk_group <- ifelse(surv_df_v$risk_score >= median_rs_v, "High", "Low")

  cat("验证集有效样本:", nrow(surv_df_v), "\n")
  cat("High:", sum(surv_df_v$risk_group == "High"), "\n")
  cat("Low:",  sum(surv_df_v$risk_group == "Low"), "\n")

  write.csv(surv_df_v, "results/tables/Table11_RiskScore_ValidationSet.csv", row.names = FALSE)
  saveRDS(surv_df_v, "data/processed/surv_df_valid.rds")

  if (sum(surv_df_v$risk_group == "High") >= 5 && sum(surv_df_v$risk_group == "Low") >= 5) {
    # KM (验证集)
    fit_v <- survfit(Surv(OS_time, OS_event) ~ risk_group, data = surv_df_v)
    pval_v <- surv_pvalue(fit_v, data = surv_df_v)$pval

    fig2h <- ggsurvplot(
      fit_v,
      data       = surv_df_v,
      pval       = TRUE,
      risk.table = TRUE,
      risk.table.height = 0.22,
      palette    = c("#E64B35", "#4DBBD5"),
      title      = "Kaplan-Meier: Risk Score Groups (Validation Set)",
      subtitle   = paste0("GSE159216 (n=", nrow(surv_df_v), ")"),
      xlab       = "Time (months)",
      ylab       = "Overall Survival Probability",
      legend.labs = c("High Risk", "Low Risk"),
      conf.int   = FALSE,
      surv.median.line = "hv",
      ggtheme    = theme_classic(base_size = 12)
    )

    # ggsurvplot 返回列表对象,用 newpage=FALSE 避免空白页
    pdf("results/figures/PM9_KM_ValidationSet.pdf", width = 7, height = 7)
    print(fig2h, newpage = FALSE)
    dev.off()
    cat("  -> PM9 完成\n")

    # ROC (验证集)
    tryCatch({
      roc_v <- timeROC::timeROC(
        T = surv_df_v$OS_time, delta = surv_df_v$OS_event,
        marker = surv_df_v$risk_score, cause = 1,
        times = c(12, 36, 60), iid = FALSE
      )
      auc_v_df <- data.frame(time = c("1y", "3y", "5y"),
                              AUC = c(roc_v$AUC[1], roc_v$AUC[2], roc_v$AUC[3]))
      cat("验证集AUC:", paste(round(roc_v$AUC, 3), collapse=" / "), "\n")
      write.csv(auc_v_df, "results/tables/Table12_AUC_Validation.csv", row.names = FALSE)

      roc_plot_v <- rbind(
        data.frame(sens=roc_v$TP[,1], spec=1-roc_v$FP[,1], time=sprintf("1y AUC=%.3f", roc_v$AUC[1])),
        data.frame(sens=roc_v$TP[,2], spec=1-roc_v$FP[,2], time=sprintf("3y AUC=%.3f", roc_v$AUC[2])),
        data.frame(sens=roc_v$TP[,3], spec=1-roc_v$FP[,3], time=sprintf("5y AUC=%.3f", roc_v$AUC[3]))
      )

      fig2i <- ggplot(roc_plot_v, aes(x = 1-spec, y = sens, color = time)) +
        geom_line(linewidth = 1.2) +
        geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
        scale_color_manual(values = c("#E64B35", "#0072B2", "#009E73")) +
        labs(title = "Time-dependent ROC (Validation Set)", x="1-Specificity", y="Sensitivity", color="") +
        theme_classic(base_size = 12) +
        theme(legend.position = c(0.65, 0.25))

      pdf("results/figures/PM10_ROC_ValidationSet.pdf", width = 6, height = 6)
      print(fig2i)
      dev.off()
      cat("  -> PM10 完成\n")
    }, error = function(e) {
      cat("验证集ROC error:", e$message, "\n")
    })
  }
} else {
  cat("未找到GSE159216生存列,跳过外部验证\n")
  cat("已检测到的列包括:", paste(colnames(pheno159216)[grepl("overall survival", colnames(pheno159216), ignore.case=TRUE)], collapse="; "), "\n")
}

# ================================================================
# Part 10: 单个独立预后基因的KM + 箱线图 (经典展示)
# ================================================================
cat("\n=== Part10: 单基因KM + 箱线图 ===\n")

# 单基因KM (全部Hub基因)
genes_show <- final_genes_in_expr
km_dfs <- list()
km_plot_list <- list()

for (g in genes_show) {
  expr_g  <- as.numeric(expr39582_norm[g, valid_os])
  med_g   <- median(expr_g, na.rm = TRUE)
  grp_g   <- ifelse(expr_g >= med_g, "High", "Low")
  km_dfs[[g]] <- data.frame(
    OS_time  = os_delay[valid_os],
    OS_event = os_event[valid_os],
    group    = grp_g,
    gene     = g
  )

  # 生成ggsurvplot对象
  df_g <- km_dfs[[g]]
  fit_g <- survfit(Surv(OS_time, OS_event) ~ group, data = df_g)
  pv_g  <- survdiff(Surv(OS_time, OS_event) ~ group, data = df_g)
  pv_n  <- 1 - pchisq(pv_g$chisq, df = 1)

  fig_g <- ggsurvplot(
    fit_g, data = df_g,
    pval = TRUE, pval.size = 4,
    risk.table = TRUE, risk.table.height = 0.22, risk.table.col = "strata",
    palette = c("#E64B35", "#4DBBD5"),
    title = paste0("Kaplan-Meier: ", g),
    subtitle = paste0("GSE39582 (n=", nrow(df_g), ", p=", signif(pv_n, 3), ")"),
    xlab = "Time (months)",
    ylab = "Overall Survival Probability",
    legend.labs = c("High Risk", "Low Risk"),
    conf.int = FALSE,
    surv.median.line = "hv",
    ggtheme = theme_classic(base_size = 12)
  )

  km_plot_list[[g]] <- fig_g$plot
}

# 用patchwork拼接单基因KM图
if (length(km_plot_list) > 0) {
  combined_km <- wrap_plots(km_plot_list, ncol = 2)
  pdf("results/figures/PM11_SingleGene_KM.pdf", width = 14, height = 7 * ceiling(length(km_plot_list)/2))
  print(combined_km)
  dev.off()
  cat("  -> PM11 完成\n")
}

# 单基因箱线图 (Tumor vs Normal 表达) + 添加P值
box_df_list <- lapply(final_genes_in_expr, function(g) {
  data.frame(
    expr  = as.numeric(expr39582_norm[g, ]),
    group = group_39582,
    gene  = g
  )
})
box_df <- do.call(rbind, box_df_list) %>%
  filter(group %in% c("Tumor", "Normal"))

# 计算每个基因的wilcoxon p值
pval_df <- box_df %>%
  group_by(gene) %>%
  summarise(
    p_val = wilcox.test(expr ~ group, exact = FALSE)$p.value,
    p_label = ifelse(p_val < 0.001, "p<0.001",
                     ifelse(p_val < 0.01, sprintf("p=%.3f", p_val),
                            sprintf("p=%.3f", p_val))),
    y_pos = max(expr, na.rm = TRUE) * 1.05,
    .groups = "drop"
  )

fig2k <- ggplot(box_df, aes(x = group, y = expr, fill = group)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  scale_fill_manual(values = c(Tumor = "#E64B35", Normal = "#4DBBD5")) +
  facet_wrap(~gene, scales = "free_y", ncol = 4) +
  # 添加P值标注
  geom_text(data = pval_df, aes(x = 1.5, y = y_pos, label = p_label),
            inherit.aes = FALSE, size = 3, vjust = 0) +
  labs(title = "Expression of Prognostic UbRGs (Tumor vs Normal)",
       x = "", y = "Normalized Expression") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

pdf("results/figures/PM12_Boxplot_Genes.pdf",
    width = 4 * min(4, length(final_genes_in_expr)),
    height = 3.5 * ceiling(length(final_genes_in_expr) / 4))
print(fig2k)
dev.off()
cat("  -> PM12 完成\n")

# ================================================================
# Part 11: Hub基因功能验证补充
# ================================================================
cat("\n=== Part11: Hub基因功能验证补充 ===\n")

## --- Fig2L: Hub基因单因素Cox森林图 ---
cat("绘制 PM13: Hub基因单因素Cox森林图\n")

# [Cox结果提取说明]
# summary(coxph)$coefficients[, "exp(coef)"] → HR (已是HR尺度)
# summary(coxph)$conf.int[, "lower .95"]     → CI下限 (已是HR尺度, 勿再exp())
# summary(coxph)$conf.int[, "upper .95"]     → CI上限 (已是HR尺度, 勿再exp())

# 计算8个基因的单因素Cox
unicox_8genes <- lapply(final_genes_in_expr, function(g) {
  expr_g <- as.numeric(expr39582_norm[g, valid_os])
  fit <- coxph(Surv(os_delay[valid_os], os_event[valid_os]) ~ expr_g)
  s <- summary(fit)
  data.frame(
    gene = g,
    HR = s$coefficients["expr_g", "exp(coef)"],
    CI_low = s$conf.int["expr_g", "lower .95"],
    CI_high = s$conf.int["expr_g", "upper .95"],
    p_value = s$coefficients["expr_g", "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
})
unicox_8df <- do.call(rbind, unicox_8genes) %>%
  arrange(HR) %>%
  mutate(
    gene = factor(gene, levels = gene),
    risk = ifelse(HR > 1, "High risk", "Low risk"),
    hr_label = sprintf("%.2f", HR),
    ci_label = sprintf("(%.2f-%.2f)", CI_low, CI_high),
    p_label = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  )

# 森林图
# log尺度森林图
log_lo_8 <- floor(log10(min(unicox_8df$CI_low, na.rm = TRUE) * 0.8))
log_hi_8 <- ceiling(log10(max(unicox_8df$CI_high, na.rm = TRUE) * 1.3))
log_breaks_8 <- 10^seq(log_lo_8, log_hi_8, by = 1)
log_breaks_8 <- log_breaks_8[log_breaks_8 >= min(unicox_8df$CI_low) * 0.8 &
                              log_breaks_8 <= max(unicox_8df$CI_high) * 1.4]
if (length(log_breaks_8) < 3) log_breaks_8 <- scales::log_breaks()(range(c(unicox_8df$CI_low, unicox_8df$CI_high)))

fig2l_forest <- ggplot(unicox_8df, aes(x = HR, y = gene)) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high, color = risk),
                width = 0.3, linewidth = 1, orientation = "y") +
  geom_point(aes(color = risk), size = 4) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  scale_color_manual(values = c("High risk" = "#E64B35", "Low risk" = "#4DBBD5")) +
  scale_x_log10(breaks = log_breaks_8, labels = sprintf("%.2g", log_breaks_8)) +
  labs(x = "Hazard Ratio (95% CI)", y = "") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = c(0.98, 0.5),
    legend.justification = c("right", "center"),
    legend.title = element_blank(),
    axis.text.y = element_text(face = "italic", size = 11),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dotted")
  )

# 表格
table_8 <- unicox_8df %>%
  select(gene, hr_label, ci_label, p_label) %>%
  pivot_longer(cols = c(hr_label, ci_label, p_label), names_to = "col", values_to = "val") %>%
  mutate(x_pos = case_when(col == "hr_label" ~ 1, col == "ci_label" ~ 2, col == "p_label" ~ 3))

fig2l_table <- ggplot(table_8, aes(x = x_pos, y = gene)) +
  geom_text(aes(label = val), size = 3.5, vjust = 0.5) +
  scale_x_continuous(limits = c(0.5, 3.5), breaks = 1:3,
                     labels = c("HR", "95% CI", "P"), expand = c(0, 0)) +
  theme_void(base_size = 11) +
  theme(axis.text.x = element_text(face = "bold", size = 10))

fig2l_combined <- fig2l_forest + fig2l_table + plot_layout(widths = c(2, 1)) +
  plot_annotation(
    title = paste0("Univariate Cox Analysis of ", length(final_genes_in_expr), " Hub Genes (OS)"),
    subtitle = "Individual prognostic value for overall survival",
    theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
  )

pdf("results/figures/PM13_Forest_HubGenes_UniCox.pdf", width = 9, height = 5)
print(fig2l_combined)
dev.off()
cat("  -> PM13 完成\n")

## --- Fig2M: Hub基因表达相关性热图 ---
cat("绘制 PM14: Hub基因表达相关性热图\n")

expr_8genes <- t(expr39582_norm[final_genes_in_expr, valid_os])
colnames(expr_8genes) <- final_genes_in_expr
cor_matrix <- cor(expr_8genes, method = "spearman")

# 计算p值
cor_p <- matrix(NA, ncol(cor_matrix), ncol(cor_matrix))
for(i in 1:ncol(cor_matrix)) {
  for(j in 1:ncol(cor_matrix)) {
    cor_p[i,j] <- cor.test(expr_8genes[,i], expr_8genes[,j], method = "spearman")$p.value
  }
}

# 标记显著性
stars <- ifelse(cor_p < 0.001, "***",
                ifelse(cor_p < 0.01, "**",
                       ifelse(cor_p < 0.05, "*", "")))

pdf("results/figures/PM14_Correlation_HubGenes.pdf", width = 8, height = 7)
pheatmap(cor_matrix,
         display_numbers = stars,
         number_color = "black",
         fontsize_number = 12,
         color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
         breaks = seq(-1, 1, length.out = 101),
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize = 10,
         main = paste0("Spearman Correlation of ", length(final_genes_in_expr), " Hub Genes\n(* p<0.05, ** p<0.01, *** p<0.001)"))
dev.off()
cat("  -> PM14 完成\n")

## --- Fig2N: Hub基因风险评分贡献度 ---
cat("绘制 PM15: Hub基因风险评分贡献度\n")

# 使用多因素Cox系数(与风险评分公式一致)
# 风险评分 = sum(coef * expr),贡献度应该基于同一套系数
hub_coef_for_contrib <- multicox_coef_df$coef[match(final_genes_in_expr, multicox_coef_df$gene)]
contribution_df <- data.frame(
  gene = final_genes_in_expr,
  coef = hub_coef_for_contrib,
  mean_expr = rowMeans(expr39582_norm[final_genes_in_expr, valid_os]),
  stringsAsFactors = FALSE
) %>%
  mutate(
    contribution = coef * mean_expr,
    direction = ifelse(contribution > 0, "Risk-increasing", "Risk-decreasing"),
    gene = factor(gene, levels = gene[order(abs(contribution))])
  )

# 计算贡献度范围,用于调整x轴范围
max_contrib <- max(abs(contribution_df$contribution))
x_limit <- max_contrib * 1.4  # 留出足够空间显示标签

fig2n <- ggplot(contribution_df, aes(x = gene, y = contribution, fill = direction)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", contribution)),
            hjust = ifelse(contribution_df$contribution > 0, -0.1, 1.1),
            size = 3.5) +
  scale_fill_manual(values = c("Risk-increasing" = "#E64B35", "Risk-decreasing" = "#4DBBD5")) +
  coord_flip() +
  scale_y_continuous(limits = c(-x_limit, x_limit), expand = c(0, 0)) +
  labs(title = paste0("Risk Score Contribution of ", length(final_genes_in_expr), " Hub Genes"),
       subtitle = "Contribution = MultiCox Coefficient × Mean Expression",
       x = "", y = "Contribution to Risk Score") +
  theme_classic(base_size = 11) +
  theme(legend.position = c(0.98, 0.02),
        legend.justification = c("right", "bottom"),
        legend.title = element_blank())

pdf("results/figures/PM15_Contribution_HubGenes.pdf", width = 9, height = 6)
print(fig2n)
dev.off()
cat("  -> PM15 完成\n")

## --- Fig2O: Hub基因RFS单因素Cox森林图 ---
cat("绘制 PM16: Hub基因RFS单因素Cox森林图\n")

# 修复: valid_rfs必须同时满足Tumor条件(原代码混入17个Normal样本)
valid_rfs_tumor <- group_39582 == "Tumor" & !is.na(rfs_event) & !is.na(rfs_delay) & rfs_delay > 0
cat("  RFS有效Tumor样本:", sum(valid_rfs_tumor), "\n")

unicox_8genes_rfs <- lapply(final_genes_in_expr, function(g) {
  expr_g <- as.numeric(expr39582_norm[g, valid_rfs_tumor])
  fit <- coxph(Surv(rfs_delay[valid_rfs_tumor], rfs_event[valid_rfs_tumor]) ~ expr_g)
  s <- summary(fit)
  data.frame(
    gene = g,
    HR = s$coefficients["expr_g", "exp(coef)"],
    CI_low = s$conf.int["expr_g", "lower .95"],
    CI_high = s$conf.int["expr_g", "upper .95"],
    p_value = s$coefficients["expr_g", "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
})
unicox_8df_rfs <- do.call(rbind, unicox_8genes_rfs) %>%
  arrange(HR) %>%
  mutate(
    gene = factor(gene, levels = gene),
    risk = ifelse(HR > 1, "High risk", "Low risk"),
    hr_label = sprintf("%.2f", HR),
    ci_label = sprintf("(%.2f-%.2f)", CI_low, CI_high),
    p_label = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  )

# log尺度森林图
log_lo_rfs <- floor(log10(min(unicox_8df_rfs$CI_low, na.rm = TRUE) * 0.8))
log_hi_rfs <- ceiling(log10(max(unicox_8df_rfs$CI_high, na.rm = TRUE) * 1.3))
log_breaks_rfs <- 10^seq(log_lo_rfs, log_hi_rfs, by = 1)
log_breaks_rfs <- log_breaks_rfs[log_breaks_rfs >= min(unicox_8df_rfs$CI_low) * 0.8 &
                                  log_breaks_rfs <= max(unicox_8df_rfs$CI_high) * 1.4]
if (length(log_breaks_rfs) < 3) log_breaks_rfs <- scales::log_breaks()(range(c(unicox_8df_rfs$CI_low, unicox_8df_rfs$CI_high)))

fig2o_forest <- ggplot(unicox_8df_rfs, aes(x = HR, y = gene)) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high, color = risk),
                width = 0.3, linewidth = 1, orientation = "y") +
  geom_point(aes(color = risk), size = 4) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  scale_color_manual(values = c("High risk" = "#E64B35", "Low risk" = "#4DBBD5")) +
  scale_x_log10(breaks = log_breaks_rfs, labels = sprintf("%.2g", log_breaks_rfs)) +
  labs(x = "Hazard Ratio (95% CI)", y = "") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = c(0.98, 0.5),
    legend.justification = c("right", "center"),
    legend.title = element_blank(),
    axis.text.y = element_text(face = "italic", size = 11),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dotted")
  )

table_8_rfs <- unicox_8df_rfs %>%
  select(gene, hr_label, ci_label, p_label) %>%
  pivot_longer(cols = c(hr_label, ci_label, p_label), names_to = "col", values_to = "val") %>%
  mutate(x_pos = case_when(col == "hr_label" ~ 1, col == "ci_label" ~ 2, col == "p_label" ~ 3))

fig2o_table <- ggplot(table_8_rfs, aes(x = x_pos, y = gene)) +
  geom_text(aes(label = val), size = 3.5, vjust = 0.5) +
  scale_x_continuous(limits = c(0.5, 3.5), breaks = 1:3,
                     labels = c("HR", "95% CI", "P"), expand = c(0, 0)) +
  theme_void(base_size = 11) +
  theme(axis.text.x = element_text(face = "bold", size = 10))

fig2o_combined <- fig2o_forest + fig2o_table + plot_layout(widths = c(2, 1)) +
  plot_annotation(
    title = paste0("Univariate Cox Analysis of ", length(final_genes_in_expr), " Hub Genes (RFS)"),
    subtitle = "Individual prognostic value for recurrence-free survival",
    theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
  )

pdf("results/figures/PM16_Forest_HubGenes_RFS.pdf", width = 9, height = 5)
print(fig2o_combined)
dev.off()
cat("  -> PM16 完成\n")

# 保存单基因验证结果
write.csv(unicox_8df, "results/tables/Table13_HubGenes_UniCox_OS.csv", row.names = FALSE)
write.csv(unicox_8df_rfs, "results/tables/Table14_HubGenes_UniCox_RFS.csv", row.names = FALSE)
write.csv(contribution_df, "results/tables/Table15_HubGenes_Contribution.csv", row.names = FALSE)
write.csv(cor_matrix, "results/tables/Table16_HubGenes_Correlation.csv")

# ================================================================
# Part 12: RFS 预后建模系列 (PM17-PM20)
# ================================================================
cat("\n=== Part12: RFS 预后建模系列 ===\n")

# RFS 有效样本
valid_rfs_tumor <- group_39582 == "Tumor" & !is.na(rfs_event) & !is.na(rfs_delay) & rfs_delay > 0
cat("RFS 有效 Tumor 样本:", sum(valid_rfs_tumor), "\n")

if (sum(valid_rfs_tumor) > 0) {

  ## --- PM17: RFS KM 训练集 ---
  cat("绘制 PM17: RFS KM 训练集\n")

  surv_rfs <- Surv(rfs_delay[valid_rfs_tumor], rfs_event[valid_rfs_tumor])
  risk_rfs <- risk_score[valid_rfs_tumor]
  median_rfs_risk <- median(risk_rfs, na.rm = TRUE)
  group_rfs <- ifelse(risk_rfs > median_rfs_risk, "High", "Low")

  fit_rfs_km <- survfit(surv_rfs ~ group_rfs)

  pm17 <- ggsurvplot(
    fit_rfs_km,
    data = data.frame(rfs_delay = rfs_delay[valid_rfs_tumor],
                      rfs_event = rfs_event[valid_rfs_tumor],
                      group_rfs = group_rfs),
    pval = TRUE, pval.size = 4,
    risk.table = TRUE, risk.table.height = 0.22, risk.table.col = "strata",
    palette = c("#E64B35", "#4DBBD5"),
    title = "Kaplan-Meier: Recurrence-Free Survival",
    subtitle = paste0("GSE39582 Training Set (n=", sum(valid_rfs_tumor), ")"),
    xlab = "Time (months)",
    ylab = "Recurrence-Free Survival Probability",
    legend.labs = c("High Risk", "Low Risk"),
    conf.int = FALSE,
    surv.median.line = "hv",
    ggtheme = theme_classic(base_size = 12)
  )

  pdf("results/figures/PM17_KM_RFS_TrainingSet.pdf", width = 7, height = 7)
  print(pm17, newpage = FALSE)
  dev.off()
  cat("  -> PM17 完成\n")

  ## --- PM18: RFS 单基因 KM ---
  cat("绘制 PM18: RFS 单基因 KM\n")

  km_rfs_plot_list <- list()

  for (g in genes_show) {
    expr_g <- as.numeric(expr39582_norm[g, valid_rfs_tumor])
    grp_g <- ifelse(expr_g > median(expr_g, na.rm = TRUE), "High", "Low")

    df_g <- data.frame(
      RFS_time = rfs_delay[valid_rfs_tumor],
      RFS_event = rfs_event[valid_rfs_tumor],
      group = grp_g,
      gene = g
    )

    fit_g <- survfit(Surv(RFS_time, RFS_event) ~ group, data = df_g)
    pv_g <- survdiff(Surv(RFS_time, RFS_event) ~ group, data = df_g)
    pv_n <- 1 - pchisq(pv_g$chisq, df = 1)

    fig_g <- ggsurvplot(
      fit_g, data = df_g,
      pval = TRUE, pval.size = 4,
      risk.table = TRUE, risk.table.height = 0.22, risk.table.col = "strata",
      palette = c("#E64B35", "#4DBBD5"),
      title = paste0("Kaplan-Meier: ", g),
      subtitle = paste0("GSE39582 RFS (n=", nrow(df_g), ", p=", signif(pv_n, 3), ")"),
      xlab = "Time (months)",
      ylab = "Recurrence-Free Survival Probability",
      legend.labs = c("High Risk", "Low Risk"),
      conf.int = FALSE,
      surv.median.line = "hv",
      ggtheme = theme_classic(base_size = 12)
    )

    km_rfs_plot_list[[g]] <- fig_g$plot
  }

  if (length(km_rfs_plot_list) > 0) {
    combined_rfs_km <- wrap_plots(km_rfs_plot_list, ncol = 3)
    pdf("results/figures/PM18_SingleGene_KM_RFS.pdf", width = 14, height = 5 * ceiling(length(km_rfs_plot_list)/3))
    print(combined_rfs_km)
    dev.off()
    cat("  -> PM18 完成\n")
  }

  ## --- PM19: RFS 独立预后森林图 ---
  cat("绘制 PM19: RFS 独立预后森林图\n")

  # 准备临床数据
  age_rfs <- clin_39582$age[valid_rfs_tumor]
  # 【修复P0】clin_39582无gender字段，改用sex
  sex_rfs <- clin_39582$sex[valid_rfs_tumor]
  tnm_rfs <- clin_39582$tnm_stage[valid_rfs_tumor]

  # TNM 分期数值化 (tnm_stage 已是数值: 1-4)
  tnm_numeric_rfs <- tnm_rfs

  # 单因素 Cox
  unicox_vars_rfs <- list(
    Risk_Score = risk_rfs,
    TNM_Stage = tnm_numeric_rfs,
    Age = age_rfs
  )

  # 【P2/P3】重构为 data.frame + formula 方式，避免 var_val[valid_idx] 命名问题
  rfs_delay_vec  <- as.numeric(clin_39582$rfs_delay)
  rfs_event_vec  <- as.numeric(clin_39582$rfs_event)

  # 构建干净的 data.frame（只用 valid_rfs_tumor 过滤，避免 double subsetting）
  rfs_df_for_cox <- data.frame(
      x        = c(risk_rfs, tnm_numeric_rfs, age_rfs),
      time     = c(rfs_delay_vec[valid_rfs_tumor], rfs_delay_vec[valid_rfs_tumor], rfs_delay_vec[valid_rfs_tumor]),
      event    = c(rfs_event_vec[valid_rfs_tumor], rfs_event_vec[valid_rfs_tumor], rfs_event_vec[valid_rfs_tumor]),
      variable = rep(c("Risk_Score", "TNM_Stage", "Age"), each = sum(valid_rfs_tumor)),
      stringsAsFactors = FALSE
  ) %>% filter(!is.na(x), time > 0)

  unicox_results_rfs <- lapply(c("Risk_Score", "TNM_Stage", "Age"), function(var_name) {
      df_var <- filter(rfs_df_for_cox, variable == var_name)
      if (nrow(df_var) < 10) return(NULL)
      fit <- coxph(Surv(time, event) ~ x, data = df_var)
      s   <- summary(fit)
      data.frame(
          Variable = var_name,
          HR       = s$coefficients["x", "exp(coef)"],
          CI_low   = s$conf.int["x", "lower .95"],
          CI_high  = s$conf.int["x", "upper .95"],
          p_value  = s$coefficients["x", "Pr(>|z|)"],
          n        = nrow(df_var),
          stringsAsFactors = FALSE
      )
  })
  unicox_df_rfs <- do.call(rbind, unicox_results_rfs[!sapply(unicox_results_rfs, is.null)])

  if (nrow(unicox_df_rfs) > 0) {
    unicox_df_rfs <- unicox_df_rfs %>%
      mutate(
        Variable = factor(Variable, levels = Variable),
        risk = ifelse(HR > 1, "High risk", "Low risk"),
        hr_label = sprintf("%.2f", HR),
        ci_label = sprintf("(%.2f-%.2f)", CI_low, CI_high),
        p_label = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
      )

    # 森林图 - 与PM13完全一致的风格
    log_lo_19 <- floor(log10(min(unicox_df_rfs$CI_low, na.rm = TRUE) * 0.8))
    log_hi_19 <- ceiling(log10(max(unicox_df_rfs$CI_high, na.rm = TRUE) * 1.3))
    log_breaks_19 <- 10^seq(log_lo_19, log_hi_19, by = 1)
    log_breaks_19 <- log_breaks_19[log_breaks_19 >= min(unicox_df_rfs$CI_low) * 0.8 &
                                     log_breaks_19 <= max(unicox_df_rfs$CI_high) * 1.4]
    if (length(log_breaks_19) < 3) log_breaks_19 <- scales::log_breaks()(range(c(unicox_df_rfs$CI_low, unicox_df_rfs$CI_high)))
    # 确保 Age 的 CI_low 也在 x 轴范围内
    x_lo_val <- min(unicox_df_rfs$CI_low)
    if (x_lo_val < 1) {
      x_lo_val <- max(0.5, x_lo_val * 0.95)
      x_lo_10 <- floor(log10(x_lo_val))
      extra_breaks <- 10^seq(x_lo_10, log10(1), length.out = 3)
      extra_breaks <- extra_breaks[extra_breaks < 1 & extra_breaks >= x_lo_val]
      log_breaks_19 <- sort(unique(c(log_breaks_19, extra_breaks)))
    }

    # 森林图 - 与PM13完全一致
    pm19_forest <- ggplot(unicox_df_rfs, aes(x = HR, y = Variable)) +
      geom_errorbar(aes(xmin = CI_low, xmax = CI_high, color = risk),
                    width = 0.3, linewidth = 1, orientation = "y") +
      geom_point(aes(color = risk), size = 4) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.8) +
      scale_color_manual(values = c("High risk" = "#E64B35", "Low risk" = "#4DBBD5")) +
      scale_x_log10(breaks = log_breaks_19, labels = sprintf("%.2g", log_breaks_19)) +
      labs(x = "Hazard Ratio (95% CI)", y = "") +
      theme_classic(base_size = 11) +
      theme(
        legend.position = "none",
        axis.text.y = element_text(size = 11),
        panel.grid.major.x = element_line(color = "grey90", linetype = "dotted")
      )

    # 表格 - 与PM13完全一致
    table_rfs <- unicox_df_rfs %>%
      select(Variable, hr_label, ci_label, p_label) %>%
      pivot_longer(cols = c(hr_label, ci_label, p_label), names_to = "col", values_to = "val") %>%
      mutate(x_pos = case_when(col == "hr_label" ~ 1, col == "ci_label" ~ 2, col == "p_label" ~ 3))

    table_ann_rfs <- ggplot(table_rfs, aes(x = x_pos, y = Variable, label = val)) +
      geom_text(size = 3.5, vjust = 0.5) +
      scale_x_continuous(limits = c(0.5, 3.5), breaks = 1:3,
                         labels = c("HR", "95% CI", "P"), expand = c(0, 0)) +
      theme_void(base_size = 11) +
      theme(axis.text.x = element_text(face = "bold", size = 10))

    pm19_combined <- pm19_forest + table_ann_rfs + plot_layout(widths = c(2, 1)) +
      plot_annotation(
        title = "Independent Prognostic Factors (RFS)",
        subtitle = "Univariate Cox regression for recurrence-free survival",
        theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
      )

    pdf("results/figures/PM19_Independent_RFS.pdf", width = 9, height = 5)
    print(pm19_combined)
    dev.off()
    cat("  -> PM19 完成\n")
  }

  ## --- PM20: RFS 分层 KM ---
  cat("绘制 PM20: RFS 分层 KM\n")

  # 构建分层数据框(与PM8一致的方式)
  tnm_rfs_for_pm20 <- clin_39582$tnm_stage[valid_rfs_tumor]
  strat_rfs_df <- data.frame(
    risk_score = risk_rfs,
    rfs_delay = rfs_delay[valid_rfs_tumor],
    rfs_event = rfs_event[valid_rfs_tumor],
    tnm_raw = tnm_rfs_for_pm20,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      TNM_num = suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(tnm_raw)))),
      TNM_grp = case_when(
        TNM_num %in% 1:2 ~ "Stage I/II",
        TNM_num %in% 3:4 ~ "Stage III/IV",
        TRUE ~ NA_character_
      ),
      risk_group = ifelse(risk_score > median(risk_score, na.rm = TRUE), "High", "Low")
    )

  groups_rfs_strat <- unique(strat_rfs_df$TNM_grp[!is.na(strat_rfs_df$TNM_grp)])

  pm20_plots <- list()

  for (g in groups_rfs_strat) {
    sub_df <- strat_rfs_df %>% filter(TNM_grp == g, !is.na(TNM_num))
    if (sum(sub_df$risk_group == "High") < 5 || sum(sub_df$risk_group == "Low") < 5) next

    fit_s <- survfit(Surv(rfs_delay, rfs_event) ~ risk_group, data = sub_df)
    pv <- survdiff(Surv(rfs_delay, rfs_event) ~ risk_group, data = sub_df)
    pv_num <- 1 - pchisq(pv$chisq, df = 1)

    fig_strat <- ggsurvplot(
      fit_s, data = sub_df,
      pval = TRUE, pval.size = 4,
      risk.table = TRUE, risk.table.height = 0.22, risk.table.col = "strata",
      palette = c("#E64B35", "#4DBBD5"),
      title = "Kaplan-Meier: Risk Score Groups",
      subtitle = paste0(g, " (n=", nrow(sub_df), ", p=", signif(pv_num, 3), ")"),
      xlab = "Time (months)",
      ylab = "Recurrence-Free Survival Probability",
      legend.labs = c("High Risk", "Low Risk"),
      conf.int = FALSE,
      surv.median.line = "hv",
      ggtheme = theme_classic(base_size = 11)
    )

    pm20_plots[[g]] <- fig_strat$plot
  }

  if (length(pm20_plots) > 0) {
    combined_strat_rfs <- wrap_plots(pm20_plots, ncol = 2)
    pdf("results/figures/PM20_Stratified_KM_RFS.pdf", width = 14, height = 7)
    print(combined_strat_rfs)
    dev.off()
    cat("  -> PM20 完成\n")
  } else {
    cat("  -> PM20 跳过: 分层样本不足\n")
  }

  cat("\n========================================\n")
  cat("RFS 系列完成!PM17-PM20 (4张)\n")
  cat("========================================\n")

} else {
  cat("警告: RFS 有效样本不足,跳过 RFS 系列图表\n")
}

# ================================================================
# Part 13: 保存中间数据供Step4/5使用
# ================================================================
cat("\n=== 保存中间数据 ===\n")
saveRDS(surv_df,            "data/processed/surv_df_train.rds")
saveRDS(final_genes_in_expr,"data/processed/final_genes_in_expr.rds")
saveRDS(gene_coef_used,     "data/processed/gene_coef_used.rds")
saveRDS(median_rs,          "data/processed/median_rs.rds")


cat("Step3 完成!\n")
cat("最终预后基因:", paste(final_genes, collapse = ", "), "\n")
cat("图表: PM1-PM16 (OS系列 16张) + PM17-PM20 (RFS系列 4张)\n")
cat("表格: Table8-16\n")

