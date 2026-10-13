##############################################################
##  CRLM 泛素相关基因研究
##  Step 2: 基因筛选与预后建模
##
##  分析流程：
##    1. 四象限分类（Q2 U Q3 = 转移相关）
##    2. GSE72718 验证 + Tier 分层
##    3. GSE39582: 单因素Cox → LASSO → 多因素Cox → RiskScore
##
##  输出文件：
##    图表：S2A-S2D, PM1-PM2, Supp_S2A-Supp_S2O
##    表格：Table_S2_01~09, Table_S2_S1~S5
##############################################################

setwd("X:/CRLM")

# ============================================================
#  库加载
# ============================================================
library(limma)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(survival)
library(survminer)
library(tibble)
library(tidyr)
library(scales)
library(patchwork)
library(pROC)
library(glmnet)
library(timeROC)

# ============================================================
#  常量定义
# ============================================================
COL <- list(
  UP   = "#E64B35",
  DOWN = "#4DBBD5",
  NS   = "#ADB6B6",
  LM   = "#CC79A7",
  PT   = "#D55E00",
  Q1   = "#E41A1C",
  Q2   = "#377EB8",
  Q3   = "#4DAF4A",
  Q4   = "#984EA3"
)

# ============================================================
#  输出目录
# ============================================================
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ============================================================
#  PART 0: 数据加载
# ============================================================
cat("\n========== PART 0: 数据加载 ==========\n")

# GSE39582（训练集）
expr_39582  <- readRDS("data/processed/expr39582_norm.rds")
pheno_39582 <- readRDS("data/processed/pheno39582.rds")
clin_39582  <- readRDS("data/processed/clin_39582.rds")

# GSE49355（转移发现集）
expr_49355  <- readRDS("data/processed/expr49355_norm.rds")
group_49355 <- readRDS("data/processed/group49355.rds")

# GSE72718（验证集）
expr_72718  <- readRDS("data/processed/expr72718_norm.rds")
pheno_72718 <- readRDS("data/processed/pheno72718.rds")

# 泛素基因库
ub_genes <- readRDS("data/processed/ubiquitin_genes.rds")

cat(sprintf("GSE39582: %d 基因 x %d 样本\n", nrow(expr_39582), ncol(expr_39582)))
cat(sprintf("GSE49355: LM=%d, PT=%d, Normal=%d\n",
    sum(group_49355=="LM"), sum(group_49355=="PT"), sum(group_49355=="Normal")))
cat(sprintf("泛素基因: %d\n", length(ub_genes)))

# 加载 group_39582（Tumor/Normal）用于 PM2 森林图
group_39582  <- readRDS("data/processed/group_39582.rds")
os_event_tmp <- as.numeric(clin_39582$os_event)
os_time_tmp  <- as.numeric(clin_39582$os_delay)
valid_os_tmp <- !is.na(os_event_tmp) & !is.na(os_time_tmp) & os_time_tmp > 0
valid_os_tumor <- valid_os_tmp & (group_39582 == "Tumor")
cat(sprintf("PM2_Forest: OS有效Tumor = %d / 总OS = %d\n", sum(valid_os_tumor), sum(valid_os_tmp)))

# ============================================================
#  PART 1: GSE49355 四象限分类
# ============================================================
cat("\n========== PART 1: GSE49355 四象限分类 ==========\n")

# 筛选到 UbRG
expr_ub <- expr_49355[rownames(expr_49355) %in% ub_genes, ]
cat(sprintf("GSE49355 UbRG: %d 基因 x %d 样本\n", nrow(expr_ub), ncol(expr_ub)))

# Limma 设计矩阵
group_levels <- unique(as.character(group_49355))
group_factor <- factor(group_49355, levels = group_levels)
design <- model.matrix(~ 0 + group_factor)
colnames(design) <- gsub("group_factor", "", colnames(design))

# 重命名 Normal → NC（用于 AUC 计算）
if ("Normal" %in% colnames(design)) {
  colnames(design)[colnames(design) == "Normal"] <- "NC"
}
group_49355_nc <- group_49355
group_49355_nc[group_49355_nc == "Normal"] <- "NC"

cat(sprintf("Design 列: %s\n", paste(colnames(design), collapse = ", ")))
cat(sprintf("分组计数: LM=%d, PT=%d, NC=%d\n",
    sum(group_49355_nc=="LM"), sum(group_49355_nc=="PT"), sum(group_49355_nc=="NC")))

# 构建对比
contrast_list <- list()
if (all(c("LM", "NC") %in% colnames(design))) contrast_list[["LM_vs_NC"]] <- "LM - NC"
if (all(c("PT", "NC") %in% colnames(design))) contrast_list[["PT_vs_NC"]] <- "PT - NC"
if (all(c("LM", "PT") %in% colnames(design))) contrast_list[["LM_vs_PT"]] <- "LM - PT"

contrast <- do.call(makeContrasts, c(contrast_list, list(levels = design)))

# Limma 分析
fit <- lmFit(expr_ub, design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)

# 提取各对比结果
res_LM_NC <- if ("LM_vs_NC" %in% colnames(contrast)) {
  topTable(fit, coef = "LM_vs_NC", number = Inf, adjust.method = "BH")
} else data.frame()

res_PT_NC <- if ("PT_vs_NC" %in% colnames(contrast)) {
  topTable(fit, coef = "PT_vs_NC", number = Inf, adjust.method = "BH")
} else data.frame()

res_LM_PT <- if ("LM_vs_PT" %in% colnames(contrast)) {
  topTable(fit, coef = "LM_vs_PT", number = Inf, adjust.method = "BH")
} else data.frame()

# 构建四象限分类数据框
quadrant <- data.frame(gene = rownames(expr_ub))

if (nrow(res_LM_NC) > 0) {
  res_LM_NC$gene <- rownames(res_LM_NC)
  quadrant <- quadrant %>% left_join(select(res_LM_NC, gene, logFC) %>% rename(logFC_LM_NC = logFC), by = "gene")
} else quadrant$logFC_LM_NC <- NA

if (nrow(res_PT_NC) > 0) {
  res_PT_NC$gene <- rownames(res_PT_NC)
  quadrant <- quadrant %>% left_join(select(res_PT_NC, gene, logFC) %>% rename(logFC_PT_NC = logFC), by = "gene")
} else quadrant$logFC_PT_NC <- NA

if (nrow(res_LM_PT) > 0) {
  res_LM_PT$gene <- rownames(res_LM_PT)
  quadrant <- quadrant %>% left_join(select(res_LM_PT, gene, logFC) %>% rename(logFC_LM_PT = logFC), by = "gene")
} else quadrant$logFC_LM_PT <- NA

# 四象限分类：
# Q1: LM+ PT+（两者在癌症中均上调）
# Q2: LM+ PT-（LM上调，PT下调 - 转移富集）
# Q3: LM- PT-（两者均下调）
# Q4: LM- PT+（PT上调）
quadrant <- quadrant %>% mutate(
  quadrant = case_when(
    logFC_LM_NC > 0 & logFC_PT_NC > 0 ~ "Q1",
    logFC_LM_NC > 0 & logFC_PT_NC < 0 ~ "Q2",
    logFC_LM_NC < 0 & logFC_PT_NC < 0 ~ "Q3",
    logFC_LM_NC < 0 & logFC_PT_NC > 0 ~ "Q4",
    TRUE ~ "NS"
  ),
  is_Q2 = quadrant == "Q2",
  is_Q3 = quadrant == "Q3"
)

print(table(quadrant$quadrant))

focal_ubRG <- quadrant %>% filter(quadrant %in% c("Q2", "Q3")) %>% pull(gene)
cat(sprintf("\nUbRG (Q2+Q3): %d\n", length(focal_ubRG)))

# 保存 DEG 结果
saveRDS(res_LM_NC, "data/processed/DEG_GSE49355_LM_vs_NC.rds")
saveRDS(res_PT_NC, "data/processed/DEG_GSE49355_PT_vs_NC.rds")
saveRDS(res_LM_PT, "data/processed/DEG_GSE49355_LM_vs_PT.rds")

# ============================================================
#  PART 2: GSE72718 验证（LM vs PT）
# ============================================================
cat("\n========== PART 2: GSE72718 验证 ==========\n")

expr_72718_ub <- expr_72718[rownames(expr_72718) %in% ub_genes, ]
cat(sprintf("GSE72718 UbRG: %d 基因 x %d 样本\n", nrow(expr_72718_ub), ncol(expr_72718_ub)))

# 解析组织类型
tissue_col <- 'tissue type:ch1'
pheno_72718$group <- ifelse(
  grepl('colorectal cancer liver metastasis', pheno_72718[[tissue_col]]), 'LM', 'PT'
)
pheno_72718$group <- factor(pheno_72718$group, levels = c('PT', 'LM'))
design_72718 <- model.matrix(~ 0 + pheno_72718$group)
colnames(design_72718) <- levels(pheno_72718$group)

contrast_72718 <- makeContrasts(LM_vs_PT = LM - PT, levels = design_72718)

fit_72718 <- lmFit(expr_72718_ub, design_72718)
fit_72718 <- contrasts.fit(fit_72718, contrast_72718)
fit_72718 <- eBayes(fit_72718)

res_72718 <- topTable(fit_72718, coef = "LM_vs_PT", number = Inf)
res_72718$gene <- rownames(res_72718)

saveRDS(res_72718, "data/processed/DEG_GSE72718_LM_PT.rds")
write.csv(res_72718, "results/tables/Table_S2_01_DEG_GSE72718_LM_vs_PT.csv")

cat(sprintf("GSE72718 LM vs PT DEGs: %d\n", sum(res_72718$adj.P.Val < 0.05 & abs(res_72718$logFC) > 0.5)))

# ============================================================
#  PART 3: Tier 分类 + AUC 筛选
# ============================================================
cat("\n========== PART 3: Tier 分类 + AUC ==========\n")

# 基于 quadrant 和 logFC 大小的 Tier 分层
quadrant <- quadrant %>% mutate(
  Tier = case_when(
    quadrant == "Q2" & abs(logFC_LM_NC) > 1  ~ "Tier1",   # 强 LM 上调
    quadrant == "Q2" & abs(logFC_LM_NC) <= 1 ~ "Tier2",   # 弱 LM 上调
    quadrant == "Q1"                         ~ "Tier3",   # 两者均上调
    quadrant == "Q3"                         ~ "Tier4",   # 两者均下调
    TRUE                                     ~ "Tier5"    # 其他
  )
)

cat("\nTier 分布:\n")
print(table(quadrant$Tier))

# ROC 分析的 AUC 计算
cat("\n计算 AUC...\n")
auc_df <- data.frame(gene = quadrant$gene, AUC_LM_PT = NA, AUC_LM_NC = NA)

for (i in 1:nrow(auc_df)) {
  g <- auc_df$gene[i]
  if (g %in% rownames(expr_ub)) {
    expr_g <- as.numeric(expr_ub[g, ])
    
    # LM vs PT
    idx <- group_49355 %in% c("LM", "PT")
    if (sum(idx) > 10) {
      roc_obj <- tryCatch(roc(group_49355[idx], expr_g[idx], quiet = TRUE), error = function(e) NULL)
      if (!is.null(roc_obj)) auc_df$AUC_LM_PT[i] <- auc(roc_obj)
    }
    
    # LM vs NC
    idx <- group_49355_nc %in% c("LM", "NC")
    if (sum(idx) > 10) {
      roc_obj <- tryCatch(roc(group_49355_nc[idx], expr_g[idx], quiet = TRUE), error = function(e) NULL)
      if (!is.null(roc_obj)) auc_df$AUC_LM_NC[i] <- auc(roc_obj)
    }
  }
}

quadrant <- quadrant %>% left_join(auc_df, by = "gene") %>% mutate(pass_AUC = AUC_LM_PT > 0.7 | AUC_LM_NC > 0.7)

cat(sprintf("Q2Q3 计数: %d\n", sum(quadrant$quadrant %in% c("Q2", "Q3"))))
cat(sprintf("pass_AUC 计数: %d\n", sum(quadrant$pass_AUC, na.rm=TRUE)))

# [决策] 不使用 AUC 过滤 - Q2+Q3 全部进入 Cox 分析（对齐主流 SCI 框架）
final_candidates <- quadrant %>% filter(quadrant %in% c("Q2", "Q3")) %>% pull(gene)
cat(sprintf("\nQ2+Q3（无AUC过滤）: %d\n", length(final_candidates)))

write.csv(auc_df, "results/tables/Table_S2_S2_AUC_screening.csv", row.names = FALSE)

# 四象限散点图
p_quad <- ggplot(quadrant, aes(x = logFC_PT_NC, y = logFC_LM_NC, color = quadrant)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c(Q1 = COL$Q1, Q2 = COL$Q2, Q3 = COL$Q3, Q4 = "#FF7F00", NS = "gray70")) +
  theme_minimal() +
  labs(title = "Quadrant Classification", x = "logFC (PT vs NC)", y = "logFC (LM vs NC)")

ggsave("results/figures/Supp_S2A_quadrant_scatter.pdf", p_quad, width = 10, height = 10)

# ============================================================
#  PART 4: 单因素 Cox 分析
# ============================================================
cat("\n========== PART 4: 单因素 Cox 分析 ==========\n")

# OS 和 RFS 终点
os_event <- as.numeric(clin_39582$os_event)
os_time  <- as.numeric(clin_39582$os_delay)
rfs_event <- as.numeric(clin_39582$rfs_event)
rfs_time  <- as.numeric(clin_39582$rfs_delay)

valid_os  <- !is.na(os_event) & !is.na(os_time)  & os_time > 0
valid_rfs <- !is.na(rfs_event) & !is.na(rfs_time) & rfs_time > 0

cat(sprintf("OS 有效: %d\n", sum(valid_os)))
cat(sprintf("RFS 有效: %d\n", sum(valid_rfs)))

candidates <- intersect(final_candidates, rownames(expr_39582))
cat(sprintf("候选基因在 expr_39582: %d\n", length(candidates)))

# 辅助函数：单因素 Cox
run_unicox <- function(genes, expr, event, time, endpoint) {
  valid <- !is.na(event) & !is.na(time) & time > 0
  event <- event[valid]
  time  <- time[valid]
  expr  <- expr[, valid, drop = FALSE]
  
  results <- list()
  n_success <- 0
  
  for (i in seq_along(genes)) {
    g <- genes[i]
    x <- as.numeric(expr[g, ])
    res <- tryCatch({
      fit <- coxph(Surv(time, event) ~ x)
      s <- summary(fit)
      data.frame(
        gene = g,
        HR = s$coefficients["x", "exp(coef)"],
        CI_low = s$conf.int["x", "lower .95"],
        CI_high = s$conf.int["x", "upper .95"],
        p_value = s$coefficients["x", "Pr(>|z|)"],
        stringsAsFactors = FALSE
      )
    }, error = function(e) NULL)
    
    if (!is.null(res)) { results[[i]] <- res; n_success <- n_success + 1 }
  }
  
  result_df <- do.call(rbind, results)
  cat(sprintf("  run_unicox [%s]: %d/%d 基因成功\n", endpoint, n_success, length(genes)))
  return(result_df)
}

# OS Cox
cat("\n--- OS 单因素 Cox ---\n")
unicox_os <- run_unicox(candidates, expr_39582, os_event, os_time, "OS")
if (is.null(unicox_os) || nrow(unicox_os) == 0) stop("unicox_os 为空")
unicox_os_sig <- unicox_os %>% filter(p_value < 0.05)
cat(sprintf("OS 显著: %d\n", nrow(unicox_os_sig)))

# RFS Cox
cat("\n--- RFS 单因素 Cox ---\n")
unicox_rfs <- run_unicox(candidates, expr_39582, rfs_event, rfs_time, "RFS")
if (is.null(unicox_rfs) || nrow(unicox_rfs) == 0) stop("unicox_rfs 为空")
unicox_rfs_sig <- unicox_rfs %>% filter(p_value < 0.05)
cat(sprintf("RFS 显著: %d\n", nrow(unicox_rfs_sig)))

# 合并用于 LASSO
prognosis_genes <- union(unicox_os_sig$gene, unicox_rfs_sig$gene)
cat(sprintf("\nLASSO 候选: %d\n", length(prognosis_genes)))

write.csv(unicox_os, "results/tables/Table_S2_02_UniCox_OS.csv", row.names = FALSE)
write.csv(unicox_rfs, "results/tables/Table_S2_03_UniCox_RFS.csv", row.names = FALSE)

# ============================================================
#  PART 5: LASSO + Bootstrap 稳定性
# ============================================================
cat("\n========== PART 5: LASSO + Bootstrap ==========\n")

# LASSO 输入
expr_lasso <- t(expr_39582[prognosis_genes, valid_os])
surv_lasso <- Surv(os_time[valid_os], os_event[valid_os])

set.seed(123)
cv_fit <- cv.glmnet(expr_lasso, surv_lasso, family = "cox", alpha = 1, nfolds = 10)

lambda_min <- cv_fit$lambda.min
lambda_1se <- cv_fit$lambda.1se
cat(sprintf("lambda.min: %.6f\n", lambda_min))
cat(sprintf("lambda.1se: %.6f\n", lambda_1se))

lasso_coef <- coef(cv_fit, s = "lambda.min")
lasso_genes <- rownames(lasso_coef)[lasso_coef[, 1] != 0]
cat(sprintf("LASSO 选中 (lambda.min): %d\n", length(lasso_genes)))

# Bootstrap 稳定性分析
cat("\n--- Bootstrap 稳定性 (100 次迭代) ---\n")
n_bootstrap <- 100
bootstrap_selected <- vector("list", n_bootstrap)

for (i in seq_len(n_bootstrap)) {
  set.seed(123 + i)
  boot_idx <- sample(nrow(expr_lasso), replace = TRUE)
  boot_x <- expr_lasso[boot_idx, , drop = FALSE]
  boot_y <- surv_lasso[boot_idx]
  
  tryCatch({
    boot_fit <- glmnet(boot_x, boot_y, family = "cox", lambda = lambda_min, alpha = 1)
    boot_coef <- coef(boot_fit, s = "lambda.min")
    boot_selected <- rownames(boot_coef)[which(boot_coef != 0)]
    bootstrap_selected[[i]] <- boot_selected
  }, error = function(e) {
    bootstrap_selected[[i]] <<- character(0)
  })
  
  if (i %% 20 == 0) cat(sprintf("  Bootstrap: %d/%d\n", i, n_bootstrap))
}

# 频率统计
all_genes <- unique(unlist(bootstrap_selected))
freq_vec <- sapply(all_genes, function(g) sum(sapply(bootstrap_selected, function(x) g %in% x)))
freq_df <- data.frame(
  Gene = names(freq_vec),
  Freq = as.numeric(freq_vec),
  Rate = as.numeric(freq_vec) / n_bootstrap
) %>% arrange(desc(Rate))

cutoff <- 0.80
stable_genes <- freq_df$Gene[freq_df$Rate >= cutoff]
cat(sprintf("稳定基因 (>=%d%% 在 %d 次 bootstrap): %d\n", cutoff*100, n_bootstrap, length(stable_genes)))
cat("前10:\n")
print(head(freq_df, 10))

# Bootstrap 柱状图
if (nrow(freq_df) > 0) {
  p_bootstrap <- ggplot(freq_df[1:min(30, nrow(freq_df)), ], aes(x = reorder(Gene, Rate), y = Rate * 100)) +
    geom_bar(stat = "identity", fill = "#4DBBD5") +
    geom_hline(yintercept = cutoff * 100, color = "#E64B35", linetype = "dashed") +
    labs(x = "基因", y = "选择率 (%)", title = "Bootstrap LASSO 稳定性",
         subtitle = sprintf("%d/%d 基因 >= %.0f%% 在 %d 次 bootstrap", length(stable_genes), length(all_genes), cutoff*100, n_bootstrap)) +
    coord_flip() + theme_bw(base_size = 11)
  ggsave("results/figures/S2D_Bootstrap_Stability.pdf", p_bootstrap, width = 10, height = max(5, nrow(freq_df)*0.3 + 2))
}

saveRDS(freq_df, "data/processed/bootstrap_freq.rds")
saveRDS(stable_genes, "data/processed/stable_genes.rds")

# 使用稳定基因进行多因素 Cox（若 < 3 则回退到 lasso_genes）
multi_cand_genes <- if (length(stable_genes) >= 3) stable_genes else lasso_genes
cat(sprintf("\n多因素 Cox 候选: %d\n", length(multi_cand_genes)))

saveRDS(cv_fit, "data/processed/cv_lasso.rds")
write.csv(data.frame(gene = lasso_genes), "results/tables/Table_S2_04_LASSO_genes.csv", row.names = FALSE)

# CV 图
pdf("results/figures/S2C_LASSO_CV.pdf", width = 8, height = 6)
plot(cv_fit)
abline(v = log(lambda_min), col = "red", lty = 2)
abline(v = log(lambda_1se), col = "blue", lty = 3)
dev.off()

# ============================================================
#  PART 6: 多因素 Cox（Hub 基因选择）
# ============================================================
cat("\n========== PART 6: 多因素 Cox ==========\n")

expr_multi <- t(expr_39582[multi_cand_genes, valid_os])
multi_data <- data.frame(time = os_time[valid_os], event = os_event[valid_os], expr_multi)

formula_str <- paste0("Surv(time, event) ~ ", paste(multi_cand_genes, collapse = " + "))
multi_fit <- coxph(as.formula(formula_str), data = multi_data)
s_multi <- summary(multi_fit)

multi_res <- data.frame(
  gene    = names(coef(multi_fit)),
  HR      = exp(coef(multi_fit)),
  CI_low  = exp(confint(multi_fit))[, 1],
  CI_high = exp(confint(multi_fit))[, 2],
  p_value = s_multi$coefficients[, "Pr(>|z|)"]
) %>% mutate(sig = p_value < 0.05)

hub_genes <- multi_res %>% filter(sig) %>% pull(gene)
cat(sprintf("Hub 基因 (p < 0.05): %d\n", length(hub_genes)))
cat(paste(hub_genes, collapse = ", "), "\n")

write.csv(multi_res, "results/tables/Table_S2_05_MultiCox_results.csv", row.names = FALSE)
saveRDS(multi_fit, "data/processed/multi_cox_fit.rds")
saveRDS(hub_genes, "data/processed/hub_genes.rds")

# ============================================================
#  PART 7: Risk Score 计算
# ============================================================
cat("\n========== PART 7: Risk Score ==========\n")

coef_vec <- as.numeric(coef(multi_fit)[hub_genes])
names(coef_vec) <- hub_genes

if (length(hub_genes) == 1) {
  risk_score <- as.numeric(expr_39582[hub_genes, valid_os]) * coef_vec
} else {
  expr_hub <- t(expr_39582[hub_genes, valid_os])
  risk_score <- as.numeric(expr_hub %*% coef_vec)
}

median_score <- median(risk_score)
risk_group <- ifelse(risk_score >= median_score, "High", "Low")

cat(sprintf("中位 cut-off: %.4f\n", median_score))
cat(sprintf("High: %d, Low: %d\n", sum(risk_group == "High"), sum(risk_group == "Low")))

# KM 生存图
surv_df <- data.frame(time = os_time[valid_os], event = os_event[valid_os], group = risk_group)
km_fit <- survfit(Surv(time, event) ~ group, data = surv_df)

p_km <- ggsurvplot(
  km_fit, data = surv_df,
  pval = TRUE, pval.size = 4,
  risk.table = TRUE, risk.table.height = 0.22, risk.table.col = "strata",
  palette = c("#E64B35", "#4DBBD5"),
  title = "Kaplan-Meier: Risk Score Groups",
  subtitle = paste0("GSE39582 Full Cohort (n=", nrow(surv_df), ")"),
  xlab = "Time (months)",
  ylab = "Overall Survival Probability",
  legend.labs = c("High Risk", "Low Risk"),
  conf.int = FALSE,
  surv.median.line = "hv",
  ggtheme = theme_classic(base_size = 12)
)

pdf("results/figures/S2A_KM_RiskScore.pdf", width = 7, height = 7)
print(p_km, newpage = FALSE)
dev.off()

# 时间依赖 ROC
os_time_days <- os_time[valid_os] * (365.25 / 12)
roc_res <- timeROC(T = os_time_days, delta = os_event[valid_os], marker = risk_score,
                   cause = 1, weighting = "marginal", times = c(365, 3*365, 5*365))

roc_plot <- rbind(
  data.frame(FPR=roc_res$FP[,1], TPR=roc_res$TP[,1], Time=sprintf("1y AUC=%.3f", roc_res$AUC[1])),
  data.frame(FPR=roc_res$FP[,2], TPR=roc_res$TP[,2], Time=sprintf("3y AUC=%.3f", roc_res$AUC[2])),
  data.frame(FPR=roc_res$FP[,3], TPR=roc_res$TP[,3], Time=sprintf("5y AUC=%.3f", roc_res$AUC[3]))
)

p_roc <- ggplot(roc_plot, aes(x = FPR, y = TPR, color = Time)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c("#E64B35", "#0072B2", "#009E73")) +
  labs(title = "Time-dependent ROC (Training Set)", x = "1-Specificity", y = "Sensitivity", color = "") +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.65, 0.25))

pdf("results/figures/S2B_ROC_TimeDependent.pdf", width = 6, height = 6)
print(p_roc)
dev.off()

# 保存 risk score
risk_df <- data.frame(sample = colnames(expr_39582)[valid_os], score = risk_score, group = risk_group)
write.csv(risk_df, "results/tables/Table_S2_06_RiskScore.csv", row.names = FALSE)
saveRDS(risk_score, "data/processed/risk_score.rds")

# ============================================================
#  PART 8: 保存输出
# ============================================================
cat("\n========== PART 8: 保存输出 ==========\n")

saveRDS(final_candidates, "data/processed/focal_ubRG.rds")
saveRDS(prognosis_genes, "data/processed/prognosis_genes.rds")
saveRDS(lasso_genes, "data/processed/lasso_genes.rds")
saveRDS(quadrant, "data/processed/quadrant_classification.rds")
saveRDS(coef_vec, "data/processed/hub_coef.rds")

# group_39582 供 Step3 使用
group_for_step3 <- ifelse(pheno_39582$dataset == "Non Tumoral", "Normal", "Tumor")
names(group_for_step3) <- colnames(expr_39582)
saveRDS(group_for_step3, "data/processed/group_39582.rds")
saveRDS(hub_genes, "data/processed/final_genes.rds")
saveRDS(data.frame(gene = names(coef_vec), coef = coef_vec), "data/processed/final_coef_df.rds")

multi_res_export <- data.frame(gene = multi_res$gene, HR = multi_res$HR, CI_low = multi_res$CI_low, CI_high = multi_res$CI_high, p_value = multi_res$p_value)
saveRDS(multi_res_export, "data/processed/multi_cox_res.rds")

cat(sprintf("  [Step3] %d Hub 基因 -> final_genes.rds\n", length(hub_genes)))
cat(sprintf("  [Step3] %d 系数 -> multicox_coef.rds\n", length(coef_vec)))

# ============================================================
#  PART 9: 单基因 KM 分析（完整队列）
# ============================================================
cat("\n========== PART 9: 单基因 KM（完整队列） ==========\n")

km_full_results <- data.frame(gene = character(), cohort = character(), n = integer(), HR = numeric(), CI_low = numeric(), CI_high = numeric(), p_logrank = numeric(), stringsAsFactors = FALSE)

for (g in hub_genes) {
  expr_g <- as.numeric(expr_39582[g, valid_os])
  median_g <- median(expr_g)
  group_g <- ifelse(expr_g >= median_g, "High", "Low")
  
  surv_g <- data.frame(time = os_time[valid_os], event = os_event[valid_os], group = factor(group_g, levels = c("Low", "High")))
  
  if (nrow(surv_g) == 0 || sum(!is.na(surv_g$time) & !is.na(surv_g$event)) < 10) next
  
  fit_g <- survfit(Surv(time, event) ~ group, data = surv_g)
  p_val <- tryCatch({ sd_obj <- survdiff(Surv(time, event) ~ group, data = surv_g); 1 - pchisq(sd_obj$chisq, df = 1) }, error = function(e) NA)
  
  cox_result <- tryCatch({
    cox_g <- coxph(Surv(time, event) ~ group, data = surv_g)
    list(hr = as.numeric(exp(coef(cox_g)))[1], ci_low = as.numeric(exp(confint(cox_g)))[1], ci_high = as.numeric(exp(confint(cox_g)))[2])
  }, error = function(e) list(hr = NA, ci_low = NA, ci_high = NA))
  
  km_full_results <- rbind(km_full_results, data.frame(gene = g, cohort = "Full", n = sum(valid_os), HR = cox_result$hr, CI_low = cox_result$ci_low, CI_high = cox_result$ci_high, p_logrank = p_val))
  
  # KM 图（字母代码 B-G 用于完整队列）
  gene_letters <- c("B", "C", "D", "E", "F", "G")
  letter_code <- gene_letters[which(hub_genes == g)]
  
  p_km_g <- ggsurvplot(
    fit_g, data = surv_g,
    pval = TRUE, pval.size = 4,
    risk.table = TRUE, risk.table.height = 0.22, risk.table.col = "strata",
    palette = c("#E64B35", "#4DBBD5"),
    title = paste0("Kaplan-Meier: ", g),
    subtitle = paste0("Full Cohort (n=", sum(valid_os), ")"),
    xlab = "Time (months)",
    ylab = "Overall Survival Probability",
    legend.labs = c("Low Risk", "High Risk"),
    conf.int = FALSE,
    surv.median.line = "hv",
    ggtheme = theme_classic(base_size = 12)
  )
  
  pdf(sprintf("results/figures/Supp_S2%s_KM_full_%s.pdf", letter_code, g), width = 7, height = 7)
  print(p_km_g, newpage = FALSE)
  dev.off()
  
  cat(sprintf("  %s: p = %.4f %s\n", g, p_val, ifelse(p_val < 0.05, "*", "")))
}

# ============================================================
#  PART 10: RFS KM 分析（完整队列）
# ============================================================
cat("\n--- RFS KM（完整队列） ---\n")

km_full_rfs_results <- data.frame()

for (g in hub_genes) {
  expr_g_rfs <- as.numeric(expr_39582[g, valid_rfs])
  median_g_rfs <- median(expr_g_rfs)
  group_g_rfs <- ifelse(expr_g_rfs >= median_g_rfs, "High", "Low")
  
  surv_g_rfs <- data.frame(time = rfs_time[valid_rfs], event = rfs_event[valid_rfs], group = factor(group_g_rfs, levels = c("Low", "High")))
  fit_g_rfs <- survfit(Surv(time, event) ~ group, data = surv_g_rfs)
  p_val_rfs <- tryCatch({ sd_obj_rfs <- survdiff(Surv(time, event) ~ group, data = surv_g_rfs); 1 - pchisq(sd_obj_rfs$chisq, df = 1) }, error = function(e) NA)
  
  km_full_rfs_results <- rbind(km_full_rfs_results, data.frame(gene = g, cohort = "Full_RFS", n = sum(valid_rfs), p_logrank = p_val_rfs))
  
  gene_letters_rfs <- c("B", "C", "D", "E", "F", "G")
  letter_code_rfs <- gene_letters_rfs[which(hub_genes == g)]
  
  p_km_g_rfs <- ggsurvplot(
    fit_g_rfs, data = surv_g_rfs,
    pval = TRUE, pval.size = 4,
    risk.table = TRUE, risk.table.height = 0.22, risk.table.col = "strata",
    palette = c("#E64B35", "#4DBBD5"),
    title = paste0("Kaplan-Meier: ", g),
    subtitle = paste0("Full Cohort RFS (n=", sum(valid_rfs), ")"),
    xlab = "Time (months)",
    ylab = "Recurrence-Free Survival Probability",
    legend.labs = c("Low Risk", "High Risk"),
    conf.int = FALSE,
    surv.median.line = "hv",
    ggtheme = theme_classic(base_size = 12)
  )
  
  pdf(sprintf("results/figures/Supp_S2%s_RFS_KM_full_%s.pdf", letter_code_rfs, g), width = 7, height = 7)
  print(p_km_g_rfs, newpage = FALSE)
  dev.off()
  
  cat(sprintf("  %s RFS: p = %.4f %s\n", g, p_val_rfs, ifelse(p_val_rfs < 0.05, "*", "")))
}

# ============================================================
#  PART 11: Stage IV 亚组分析
# ============================================================
cat("\n========== PART 11: Stage IV 亚组分析 ==========\n")

stage_col <- clin_39582$tnm_stage
is_stage4 <- (stage_col == "4") & valid_os
cat(sprintf("Stage IV 样本: %d\n", sum(is_stage4)))

if (sum(is_stage4) >= 30) {
  
  # --- Stage IV OS KM ---
  cat("\n--- OS KM (Stage IV) ---\n")
  km_stage4_results <- data.frame()
  
  for (g in hub_genes) {
    expr_g <- as.numeric(expr_39582[g, is_stage4])
    median_g <- median(expr_g)
    group_g <- ifelse(expr_g >= median_g, "High", "Low")
    
    surv_g <- data.frame(time = os_time[is_stage4], event = os_event[is_stage4], group = factor(group_g, levels = c("Low", "High")))
    fit_g <- survfit(Surv(time, event) ~ group, data = surv_g)
    p_val <- tryCatch({ sd_obj <- survdiff(Surv(time, event) ~ group, data = surv_g); 1 - pchisq(sd_obj$chisq, df = 1) }, error = function(e) NA)
    
    cox_result <- tryCatch({
      cox_g <- coxph(Surv(time, event) ~ group, data = surv_g)
      list(hr = as.numeric(exp(coef(cox_g)))[1], ci_low = as.numeric(exp(confint(cox_g)))[1], ci_high = as.numeric(exp(confint(cox_g)))[2])
    }, error = function(e) list(hr = NA, ci_low = NA, ci_high = NA))
    
    km_stage4_results <- rbind(km_stage4_results, data.frame(gene = g, cohort = "StageIV", n = sum(is_stage4), HR = cox_result$hr, CI_low = cox_result$ci_low, CI_high = cox_result$ci_high, p_logrank = p_val))
    
    # KM 图（字母代码 H-M 用于 Stage IV）
    gene_letters <- c("H", "I", "J", "K", "L", "M")
    letter_code <- gene_letters[which(hub_genes == g)]
    
    p_km_g <- ggsurvplot(
      fit_g, data = surv_g,
      pval = TRUE, pval.size = 4,
      risk.table = TRUE, risk.table.height = 0.22, risk.table.col = "strata",
      palette = c("#E64B35", "#4DBBD5"),
      title = paste0("Kaplan-Meier: ", g),
      subtitle = paste0("Stage IV Subgroup (n=", sum(is_stage4), ")"),
      xlab = "Time (months)",
      ylab = "Overall Survival Probability",
      legend.labs = c("Low Risk", "High Risk"),
      conf.int = FALSE,
      surv.median.line = "hv",
      ggtheme = theme_classic(base_size = 12)
    )
    
    pdf(sprintf("results/figures/Supp_S2%s_KM_stage4_%s.pdf", letter_code, g), width = 7, height = 7)
    print(p_km_g, newpage = FALSE)
    dev.off()
    
    cat(sprintf("  %s: p = %.4f %s\n", g, p_val, ifelse(p_val < 0.05, "*", "")))
  }
  
  # --- RiskScore Stage IV ---
  cat("\n--- RiskScore Stage IV ---\n")
  expr_stage4 <- t(expr_39582[hub_genes, is_stage4])
  risk_stage4 <- as.numeric(expr_stage4 %*% coef_vec)
  median_risk_s4 <- median(risk_stage4)
  group_stage4 <- ifelse(risk_stage4 >= median_risk_s4, "High", "Low")
  
  surv_stage4 <- data.frame(time = os_time[is_stage4], event = os_event[is_stage4], group = factor(group_stage4, levels = c("Low", "High")))
  fit_rs_s4 <- survfit(Surv(time, event) ~ group, data = surv_stage4)
  sd_rs_s4 <- survdiff(Surv(time, event) ~ group, data = surv_stage4)
  p_rs_s4 <- 1 - pchisq(sd_rs_s4$chisq, df = 1)
  cat(sprintf("RiskScore Stage IV p = %.4f %s\n", p_rs_s4, ifelse(p_rs_s4 < 0.05, "*", "")))
  
  p_km_stage4 <- ggsurvplot(
    fit_rs_s4, data = surv_stage4,
    pval = TRUE, pval.size = 4,
    risk.table = TRUE, risk.table.height = 0.22, risk.table.col = "strata",
    palette = c("#E64B35", "#4DBBD5"),
    title = "Kaplan-Meier: Risk Score Groups",
    subtitle = paste0("Stage IV Subgroup (n=", nrow(surv_stage4), ")"),
    xlab = "Time (months)",
    ylab = "Overall Survival Probability",
    legend.labs = c("High Risk", "Low Risk"),
    conf.int = FALSE,
    surv.median.line = "hv",
    ggtheme = theme_classic(base_size = 12)
  )
  
  pdf("results/figures/Supp_S2N_KM_RiskScore_StageIV.pdf", width = 7, height = 7)
  print(p_km_stage4, newpage = FALSE)
  dev.off()
  
  # --- Stage IV RFS KM ---
  is_stage4_rfs <- (stage_col == "4") & valid_rfs
  cat(sprintf("Stage IV RFS 样本: %d\n", sum(is_stage4_rfs)))
  
  if (sum(is_stage4_rfs) >= 15) {
    cat("\n--- RFS KM (Stage IV) ---\n")
    km_s4_rfs_results <- data.frame()
    
    for (g in hub_genes) {
      expr_g_s4_rfs <- as.numeric(expr_39582[g, is_stage4_rfs])
      median_g_s4_rfs <- median(expr_g_s4_rfs)
      group_g_s4_rfs <- ifelse(expr_g_s4_rfs >= median_g_s4_rfs, "High", "Low")
      
      surv_g_s4_rfs <- data.frame(time = rfs_time[is_stage4_rfs], event = rfs_event[is_stage4_rfs], group = factor(group_g_s4_rfs, levels = c("Low", "High")))
      fit_g_s4_rfs <- survfit(Surv(time, event) ~ group, data = surv_g_s4_rfs)
      p_val_s4_rfs <- tryCatch({ sd_obj_s4_rfs <- survdiff(Surv(time, event) ~ group, data = surv_g_s4_rfs); 1 - pchisq(sd_obj_s4_rfs$chisq, df = 1) }, error = function(e) NA)
      
      km_s4_rfs_results <- rbind(km_s4_rfs_results, data.frame(gene = g, cohort = "StageIV_RFS", n = sum(is_stage4_rfs), p_logrank = p_val_s4_rfs))
      
      gene_letters_s4_rfs <- c("H", "I", "J", "K", "L", "M")
      letter_code_s4_rfs <- gene_letters_s4_rfs[which(hub_genes == g)]
      
      p_km_g_s4_rfs <- ggsurvplot(
        fit_g_s4_rfs, data = surv_g_s4_rfs,
        pval = TRUE, pval.size = 4,
        risk.table = TRUE, risk.table.height = 0.22, risk.table.col = "strata",
        palette = c("#E64B35", "#4DBBD5"),
        title = paste0("Kaplan-Meier: ", g),
        subtitle = paste0("Stage IV Subgroup RFS (n=", sum(is_stage4_rfs), ")"),
        xlab = "Time (months)",
        ylab = "Recurrence-Free Survival Probability",
        legend.labs = c("Low Risk", "High Risk"),
        conf.int = FALSE,
        surv.median.line = "hv",
        ggtheme = theme_classic(base_size = 12)
      )
      
      pdf(sprintf("results/figures/Supp_S2%s_RFS_KM_stage4_%s.pdf", letter_code_s4_rfs, g), width = 7, height = 7)
      print(p_km_g_s4_rfs, newpage = FALSE)
      dev.off()
      
      cat(sprintf("  %s Stage IV RFS: p = %.4f %s\n", g, p_val_s4_rfs, ifelse(p_val_s4_rfs < 0.05, "*", "")))
    }
    
    # RiskScore Stage IV RFS
    expr_stage4_rfs <- t(expr_39582[hub_genes, is_stage4_rfs])
    risk_stage4_rfs <- as.numeric(expr_stage4_rfs %*% coef_vec)
    median_risk_s4_rfs <- median(risk_stage4_rfs)
    group_stage4_rfs <- ifelse(risk_stage4_rfs >= median_risk_s4_rfs, "High", "Low")
    
    surv_s4_rfs <- data.frame(time = rfs_time[is_stage4_rfs], event = rfs_event[is_stage4_rfs], group = factor(group_stage4_rfs, levels = c("Low", "High")))
    fit_rs_s4_rfs <- survfit(Surv(time, event) ~ group, data = surv_s4_rfs)
    sd_rs_s4_rfs <- survdiff(Surv(time, event) ~ group, data = surv_s4_rfs)
    p_rs_s4_rfs <- 1 - pchisq(sd_rs_s4_rfs$chisq, df = 1)
    cat(sprintf("RiskScore Stage IV RFS p = %.4f %s\n", p_rs_s4_rfs, ifelse(p_rs_s4_rfs < 0.05, "*", "")))
    
    p_km_s4_rfs <- ggsurvplot(
      fit_rs_s4_rfs, data = surv_s4_rfs,
      pval = TRUE, pval.size = 4,
      risk.table = TRUE, risk.table.height = 0.22, risk.table.col = "strata",
      palette = c("#E64B35", "#4DBBD5"),
      title = "Kaplan-Meier: Risk Score Groups",
      subtitle = paste0("Stage IV Subgroup RFS (n=", nrow(surv_s4_rfs), ")"),
      xlab = "Time (months)",
      ylab = "Recurrence-Free Survival Probability",
      legend.labs = c("High Risk", "Low Risk"),
      conf.int = FALSE,
      surv.median.line = "hv",
      ggtheme = theme_classic(base_size = 12)
    )
    
    pdf("results/figures/Supp_S2N_RFS_KM_RiskScore_StageIV.pdf", width = 7, height = 7)
    print(p_km_s4_rfs, newpage = FALSE)
    dev.off()
  } else {
    cat("  Stage IV RFS 样本 < 15，跳过\n")
  }
  
  # --- 森林图：完整队列 vs Stage IV ---
  unicox_full <- run_unicox(hub_genes, expr_39582[, valid_os_tumor], os_event[valid_os_tumor], os_time[valid_os_tumor], "OS_Full")
  unicox_stage4 <- run_unicox(hub_genes, expr_39582[, is_stage4], os_event[is_stage4], os_time[is_stage4], "OS_Stage4")
  
  n_full <- sum(valid_os)
  n_s4 <- sum(is_stage4)
  
  forest_data <- unicox_full %>%
    mutate(cohort = paste0("Full cohort (n=", n_full, ")")) %>%
    bind_rows(unicox_stage4 %>% mutate(cohort = paste0("Stage IV (n=", n_s4, ")"))) %>%
    mutate(
      cohort = factor(cohort, levels = c(paste0("Full cohort (n=", n_full, ")"), paste0("Stage IV (n=", n_s4, ")"))),
      gene = factor(gene, levels = hub_genes)
    )
  
  lbl_full <- paste0("Full cohort (n=", n_full, ")")
  lbl_s4 <- paste0("Stage IV (n=", n_s4, ")")
  cohort_colors <- setNames(c("#E64B35", "#4DBBD5"), c(lbl_full, lbl_s4))
  x_max_data <- max(forest_data$CI_high, na.rm = TRUE) * 1.5
  
  p_forest <- ggplot(forest_data, aes(x = HR, y = gene, color = cohort)) +
    geom_point(size = 4, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0.3, linewidth = 1, orientation = "y", position = position_dodge(width = 0.5)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.8) +
    scale_x_log10(limits = c(NA, x_max_data)) +
    scale_color_manual(values = cohort_colors, name = NULL) +
    guides(color = guide_legend(byrow = TRUE)) +
    theme_classic(base_size = 11) +
    labs(x = "Hazard Ratio (95% CI)", y = "", title = "Full Cohort vs Stage IV Comparison") +
    theme(
      legend.position = c(0.90, 0.5),
      legend.justification = c(0, 0.5),
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_text(size = 10, color = "#333333"),
      legend.margin = margin(l = 4, r = 4, t = 4, b = 4),
      legend.background = element_rect(fill = "white", color = NA),
      axis.text.y = element_text(face = "italic", size = 11),
      panel.grid.major.x = element_line(color = "grey90", linetype = "dotted"),
      plot.margin = margin(l = 5, r = 80, t = 10, b = 5)
    )
  
  ggsave("results/figures/PM2_Forest_Comparison.pdf", p_forest, width = 10, height = 6)
  
  # --- 保存 Stage IV 结果 ---
  stage4_results <- data.frame(
    gene = hub_genes,
    HR_full = unicox_full$HR,
    p_full = unicox_full$p_value,
    HR_stage4 = unicox_stage4$HR,
    p_stage4 = unicox_stage4$p_value
  )
  write.csv(stage4_results, "results/tables/Table_S2_07_StageIV_validation.csv", row.names = FALSE)
  
  # --- 汇总 KM 结果 ---
  km_all_results <- rbind(km_full_results, km_stage4_results)
  write.csv(km_all_results, "results/tables/Table_S2_S4_single_gene_KM.csv", row.names = FALSE)
  
  cat("  -> PM2, Supp_S2B-Supp_S2O 已生成\n")
  
} else {
  cat("  Stage IV 样本 < 30，跳过\n")
}

# ============================================================
#  最终总结
# ============================================================
cat("#  Step 2 完成！\n")
cat(sprintf("#  Q2+Q3 UbRG: %d\n", length(focal_ubRG)))
cat(sprintf("#  单因素 Cox 显著: %d\n", length(prognosis_genes)))
cat(sprintf("#  LASSO 选中: %d\n", length(lasso_genes)))
cat(sprintf("#  Bootstrap 稳定: %d\n", length(stable_genes)))
cat(sprintf("#  Hub 基因: %d (%s)\n", length(hub_genes), paste(hub_genes, collapse = ", ")))
cat("#  输出：S2A-S2D, PM2, Supp_S2A-Supp_S2O, Table_S2_01~07, Table_S2_S1~S4\n")
