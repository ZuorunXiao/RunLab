##############################################################
##  CRLM 泛素相关基因研究 - Step 5: 治疗反应预测
##
##  【分析逻辑】
##   基于风险评分对治疗反应进行预测，聚焦免疫治疗和靶向治疗相关：
##   1. 免疫治疗响应预测 (IPS - Immunophenoscore)
##   2. TIDE算法代理（T细胞功能失调评分 & TMB代理）
##   3. 化疗/靶向药物敏感性预测（GDSC数据库代理基因集）
##   4. 风险分组治疗决策可视化（四象限图）
##############################################################

setwd("X:/CRLM")
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(ggpubr)
library(scales)
library(GSVA)  # for ssGSEA in TIDE calculation

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)

# ================================================================
# Part 0: 加载数据
# ================================================================
expr39582_norm  <- readRDS("data/processed/expr39582_norm.rds")
surv_df_train   <- readRDS("data/processed/surv_df_train.rds")
ssgsea_df       <- readRDS("data/processed/ssgsea_df.rds")
estimate_df     <- readRDS("data/processed/estimate_df.rds")

tumor_samps <- surv_df_train$sample_id
expr_tumor  <- expr39582_norm[, tumor_samps, drop = FALSE]
risk_group  <- setNames(surv_df_train$risk_group, surv_df_train$sample_id)
risk_score  <- setNames(surv_df_train$risk_score,  surv_df_train$sample_id)

cat("Tumor样本数:", ncol(expr_tumor), "\n")

# ================================================================
# Part 1: IPS (Immunophenoscore) 免疫表型评分
# ================================================================

# IPS 四大模块基因集（来自Charoentong 2017）
ips_genes <- list(
  MHC = c("HLA-A","HLA-B","HLA-C","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1",
           "HLA-DRA","HLA-DRB1","HLA-E","HLA-F","TAPBP","B2M","CALR"),
  CP  = c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT",
           "CD244","CD160","CD276","VTCN1","SIGLEC15"),
  EC  = c("CD8A","CD8B","GZMB","PRF1","EOMES","TBX21"),
  SC  = c("FOXP3","IL2RA","ENTPD1","TNFRSF9")  # CTLA4 only in CP module (not duplicated)
)

# 计算每个模块的z-score均值
calc_ips_module <- function(genes, expr_mat) {
  genes_use <- intersect(genes, rownames(expr_mat))
  if (length(genes_use) < 2) return(rep(0, ncol(expr_mat)))
  module_expr <- expr_mat[genes_use, , drop = FALSE]
  # z-score per gene
  module_z <- t(scale(t(module_expr)))
  colMeans(module_z, na.rm = TRUE)
}

ips_df <- data.frame(
  sample_id = colnames(expr_tumor),
  MHC       = calc_ips_module(ips_genes$MHC, expr_tumor),
  CP        = calc_ips_module(ips_genes$CP,  expr_tumor),
  EC        = calc_ips_module(ips_genes$EC,  expr_tumor),
  SC        = calc_ips_module(ips_genes$SC,  expr_tumor)
) %>%
  mutate(
    # IPS = z_MHC + z_EC + z_SC - z_CP  (Charoentong et al., Cell Reports 2017)
    IPS = MHC + EC + SC - CP,
    risk_group = risk_group[sample_id],
    risk_score = as.numeric(risk_score[sample_id])
  )

# 高低风险组IPS差异
ips_p <- wilcox.test(ips_df$IPS[ips_df$risk_group == "High"],
                      ips_df$IPS[ips_df$risk_group == "Low"], exact = FALSE)$p.value
ips_rho <- cor.test(ips_df$risk_score, ips_df$IPS, method = "spearman")
cat("IPS: High vs Low p=", signif(ips_p, 3), "\n")
cat("IPS-RiskScore Spearman rho=", round(ips_rho$estimate, 3),
    "p=", signif(ips_rho$p.value, 3), "\n")

write.csv(ips_df, "results/tables/Table21_IPS_Scores.csv", row.names = FALSE)

# ================================================================
# Part 2: TIDE-like 功能失调评分
# ================================================================

# 效应T细胞功能失调标志（来自TIDE文献）
dysfunction_genes <- c("PDCD1","HAVCR2","LAG3","TIGIT","CD244","CTLA4","EOMES",
                        "TOX","TOX2","ENTPD1","NR4A1","BATF","PRDM1")
exclusion_genes   <- c("TGFB1","TGFB2","WNT5A","VEGFA","CXCL5","CSF1","MMP9",
                        "CAV1","FAP","PDPN","ACTA2","S100A4","PECAM1")

dys_genes_use <- intersect(dysfunction_genes, rownames(expr_tumor))
exc_genes_use <- intersect(exclusion_genes, rownames(expr_tumor))

# 用ssGSEA计算enrichment score（而非简单colMeans）
tide_gene_sets <- list(Dysfunction = dys_genes_use, Exclusion = exc_genes_use)
tide_ssgsea <- tryCatch({
  param <- ssgseaParam(as.matrix(expr_tumor), tide_gene_sets, normalize = TRUE)
  gsva(param, verbose = FALSE)
}, error = function(e) {
  # fallback to old API
  gsva(as.matrix(expr_tumor), tide_gene_sets, method = "ssgsea",
       ssgsea.norm = TRUE, verbose = FALSE)
})

# tide_ssgsea: 行=基因集, 列=样本
tide_df <- data.frame(
  sample_id      = colnames(expr_tumor),
  Dysfunction    = as.numeric(tide_ssgsea["Dysfunction", ]),
  Exclusion      = as.numeric(tide_ssgsea["Exclusion", ])
) %>%
  mutate(
    # TIDE评分 = 功能障碍 + 免疫排斥（正数=免疫逃逸倾向高，免疫治疗效果差）
    # 参考: Jiang et al., Nature Cancer 2020 - TIDE框架
    TIDE_score     = Dysfunction + Exclusion,
    risk_group     = risk_group[sample_id],
    risk_score     = as.numeric(risk_score[sample_id])
  )

tide_p <- wilcox.test(tide_df$TIDE_score[tide_df$risk_group == "High"],
                       tide_df$TIDE_score[tide_df$risk_group == "Low"], exact = FALSE)$p.value
cat("TIDE: High vs Low p=", signif(tide_p, 3), "\n")

# ================================================================
# Part 3: 药物敏感性预测 (基于已发表的CRLM相关药物基因标志)
# ================================================================

# 基于基因表达代理GDSC化疗敏感性
# 使用文献报道的结直肠癌相关药物作用靶点基因（公开数据）
drug_gene_map <- list(
  Oxaliplatin  = c("ERCC1","ERCC2","XPC","BRCA1"),
  Irinotecan   = c("TOP1","UGT1A1","CES1","CES2"),
  Bevacizumab  = c("VEGFA","VEGFB","KDR","FLT1"),
  Cetuximab    = c("EGFR","KRAS","NRAS","BRAF","PIK3CA"),
  FOLFOX       = c("ERCC1","TYMS","DPYD","UMPS"),
  Fluorouracil = c("TYMS","DPYD","OPRT","UMPS"),
  Regorafenib  = c("KIT","PDGFRA","PDGFRB","KDR","RAF1"),
  Trifluridine = c("TK1","TYMS","DPYD"),
  mTOR_inh     = c("MTOR","RPTOR","RICTOR","AKT1","PIK3CA"),
  Proteasome   = c("PSMA1","PSMA2","PSMA7","PSMB5","PSMB7","UBE2C")
)

# 计算每个药物的预测敏感性评分（基因表达均值）
drug_scores <- lapply(names(drug_gene_map), function(drug) {
  genes_use <- intersect(drug_gene_map[[drug]], rownames(expr_tumor))
  if (length(genes_use) < 2) return(NULL)
  score <- colMeans(expr_tumor[genes_use, , drop = FALSE], na.rm = TRUE)
  data.frame(sample_id = names(score), drug = drug, score = as.numeric(score))
}) %>% bind_rows()

drug_df <- drug_scores %>%
  left_join(data.frame(sample_id = names(risk_group), risk_group = risk_group),
            by = "sample_id")

# 高低风险组药物敏感性差异
drug_diff <- drug_df %>%
  group_by(drug) %>%
  summarise(
    mean_high = mean(score[risk_group == "High"], na.rm = TRUE),
    mean_low  = mean(score[risk_group == "Low"],  na.rm = TRUE),
    p_value   = wilcox.test(score[risk_group == "High"],
                             score[risk_group == "Low"], exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(adj_p = p.adjust(p_value, method = "BH"),
         log2FC = mean_high - mean_low,
         direction = ifelse(log2FC > 0, "Higher in High", "Higher in Low"))

cat("药物敏感性差异结果:\n")
print(drug_diff %>% select(drug, mean_high, mean_low, p_value, direction))
write.csv(drug_diff, "results/tables/Table22_Drug_Sensitivity.csv", row.names = FALSE)

# ================================================================
# Part 4: 治疗建议总结
# ================================================================

treat_df <- ips_df %>%
  select(sample_id, IPS, risk_group, risk_score) %>%
  left_join(tide_df %>% select(sample_id, TIDE_score), by = "sample_id") %>%
  left_join(estimate_df %>% select(sample_id, ImmuneScore, StromalScore), by = "sample_id") %>%
  mutate(
    # ICI_benefit: 高IPS(免疫原性强) + 低TIDE(免疫逃逸弱) = 免疫治疗可能获益
    # 参考: Charoentong et al. Cell 2017 (IPS), Jiang et al. Nature Cancer 2020 (TIDE)
    ICI_benefit = ifelse(IPS > median(IPS, na.rm=TRUE) & TIDE_score < median(TIDE_score, na.rm=TRUE),
                         "ICI_Benefit", "ICI_No_Benefit"),
    Recommend = case_when(
      ICI_benefit == "ICI_Benefit" & risk_group == "High" ~ "ICI + Chemo",
      ICI_benefit == "ICI_Benefit" & risk_group == "Low"  ~ "ICI Monotherapy",
      ICI_benefit != "ICI_Benefit" & risk_group == "High" ~ "Targeted + Chemo",
      TRUE ~ "Standard Chemo"
    )
  )

treatment_summary <- treat_df %>%
  count(risk_group, Recommend) %>%
  group_by(risk_group) %>%
  mutate(pct = n / sum(n) * 100)
cat("治疗建议分布:\n")
print(treatment_summary)

write.csv(treat_df, "results/tables/Table23_Treatment_Decisions.csv", row.names = FALSE)

# ================================================================
# Part 5: 可视化
# ================================================================

## --- Fig4A: IPS箱线图 ---
cat("绘制 TP1: IPS箱线图\n")
fig4a <- ggplot(ips_df, aes(x = risk_group, y = IPS, fill = risk_group)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = c(High = "#E64B35", Low = "#4DBBD5")) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     comparisons = list(c("High", "Low")),
                     hide.ns = TRUE) +
  labs(title = "Immunophenoscore (IPS)",
       subtitle = sprintf("rho = %.3f", ips_rho$estimate),
       caption  = sprintf("p = %.3f (Wilcoxon)", ips_p),
       x = "Risk Group", y = "IPS Score") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

pdf("results/figures/TP1_IPS_Boxplot.pdf", width = 5, height = 5)
print(fig4a)
dev.off()
cat("  -> TP1 完成\n")

## --- Fig4B: TIDE箱线图 ---
cat("绘制 TP2: TIDE箱线图\n")
fig4b <- ggplot(tide_df, aes(x = risk_group, y = TIDE_score, fill = risk_group)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = c(High = "#E64B35", Low = "#4DBBD5")) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     comparisons = list(c("High", "Low")),
                     hide.ns = TRUE) +
  labs(title = "Immune Evasion Score (TIDE-proxy)",
       subtitle = "Wilcoxon test",
       caption  = sprintf("p = %.3f", tide_p),
       x = "Risk Group", y = "TIDE Score") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

pdf("results/figures/TP2_TIDE_Boxplot.pdf", width = 5, height = 5)
print(fig4b)
dev.off()
cat("  -> TP2 完成\n")

## --- Fig4C: 药物敏感性条形图 ---

# 格式化p值避免0.000，并计算标签位置
format_p_drug <- function(p) {
  if (p < 0.001) return("p<0.001")
  else return(sprintf("p=%.3f", p))
}

drug_diff$p_label <- sapply(drug_diff$p_value, format_p_drug)

# 计算合适的x轴范围给标签留空间
x_range <- max(abs(drug_diff$log2FC)) * 1.3

fig4c <- ggplot(drug_diff, aes(x = reorder(drug, log2FC), y = log2FC,
                                fill = direction)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_text(aes(label = p_label,
                hjust = ifelse(log2FC >= 0, -0.1, 1.1)),
            size = 3, nudge_x = 0) +
  coord_flip() +
  scale_fill_manual(values = c("Higher in High" = "#E64B35", "Higher in Low" = "#4DBBD5")) +
  scale_y_continuous(limits = c(-x_range, x_range), expand = c(0.05, 0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(title = "Drug Sensitivity Prediction (High vs Low Risk)",
       subtitle = "Based on drug target gene expression",
       x = "Drug / Treatment", y = "Difference (High - Low)", fill = "") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top")

pdf("results/figures/TP3_Drug_Sensitivity.pdf", width = 9, height = max(5, nrow(drug_diff) * 0.4))
print(fig4c)
dev.off()
cat("  -> TP3 完成\n")

## --- TP4: 四象限治疗决策图 ---
quad_summary <- treat_df %>%
  count(Recommend) %>%
  mutate(pct = n / sum(n) * 100)

fig4d <- ggplot(quad_summary, aes(x = reorder(Recommend, pct), y = pct, fill = Recommend)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct, n)), vjust = 0.5, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c(
    "ICI + Chemo"       = "#E64B35",
    "ICI Monotherapy"   = "#F79646",
    "Targeted + Chemo"  = "#4DBBD5",
    "Standard Chemo"    = "#70AD47"
  )) +
  labs(title = "Treatment Recommendation Distribution",
       subtitle = "Based on IPS, TIDE, and Risk Group",
       x = "", y = "Percentage (%)", fill = "Treatment") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

pdf("results/figures/TP4_Treatment_Decision.pdf", width = 8, height = 5)
print(fig4d)
dev.off()
cat("  -> TP4 完成\n")

cat("Step5 完成！\n")
cat("图表: TP1-TP4 (4张)\n")
cat("表格: Table21-23\n")

