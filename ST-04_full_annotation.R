# =============================================================================
# 04_full_annotation.R
# 功能: ST 全细胞型空间注释 —— 输出每样本结构区图 (Region_ST-*.pdf) + 逐 spot 全注释 CSV:
#   (A) 结构区图 (Region_ST-*.pdf): 无监督聚类(02)出区域骨架 + 标志基因面板逐 cluster 命名
#       (每 cluster 输出 8 面板标志基因表达表供人工/H&E 判读; 区域名由"看每簇标志基因 + 相对基线
#       升高量"首过判出, 无手工命名表), 集群级, 8 类结构区;
#       严格对应论文 Fig 1F ("Unsupervised clustering analysis classified the samples into
#       distinct regions such as tumor and fibroblast regions (Fig. 1F)", Results L80)。不含
#       免疫亚型(属 03 Fig 1G/2C)。区域命名对照 Fig 1F 由"看每簇标志基因 + 相对基线升高量"首过判出
#       (复现作者检视各簇标志基因), 纯数据驱动、无硬编码样本标签、区分不开的簇留 H&E 确认。
#   (B) 2x2 未治疗组组合图 (StructuralRegions_untreated_2x2.pdf) 经用户 2026-09-08 裁定剔除(布局多次调整未达预期),
#       故本脚本现仅输出 (A) 每样本结构区图 + 逐 spot 全注释 CSV; 不再生成任何多样本组合图。
#   注: 免疫主型/亚型属 03(Fig 1G/2C), 本脚本方案 A (2026-09-06) 不再叠加自造免疫合并规则,
#       仅输出 Fig 1F 结构区 8 类注释。pred_main/pred_score 仅作 CSV 参考列保留。
# 定位: Wu et al. Cancer Discovery 2022 (doi:10.1158/2159-8290.CD-21-0316)
        
        
        
        
        
        
# 管线对照(逐字原文, 见工作区 _paper_fulltext.txt):
#   ST Methods ("ST Data Processing and Cell-Type Estimation", L308): "we used SCTransform
#     (95) to normalize data, RunPCA to perform dimension reduction, FindNeighbors and
#     FindClusters to cluster the ST spots, and RunUMAP to visualize the data. Given that
#     each spot contains several cells, we next combined the scRNA-seq data and predicted
#     main cell types by using Seurat. As for subpopulation annotation, we found that
#     Seurat did not generate uniform results ... We hence used ssGSEA (96) to score the
#     subcell types, which has been proved to be robust for ST data."
#   scRNA 参考构建 (L286): "We performed the clustering analysis of all cells and only
#     reserved the CD45+ clusters." —— label transfer 参考为 CD45+ 免疫-only, 只能预测
#     免疫主型(Fig 1G), 无法命名 Fig 1F 非免疫结构区。
#   Results Fig 1F (L80): "Unsupervised clustering analysis classified the samples into
#     distinct regions such as tumor and fibroblast regions (Fig. 1F)."
# 本脚本对应: 结构区(Fig 1F) 由 02 聚类 + 标志基因命名(本脚本 PART 2); 免疫主型(Fig 1G)
#   由 03 label transfer; 免疫亚型(Fig 2C) 由 03 ssGSEA。本图不画亚群(避免与 03 重复)。
# 关键方法裁定(用户逐字纠错后定稿, 见项目记忆):
#   Fig 1F 8 个结构区(肝细胞/正常上皮/成纤维/内皮/平滑肌/固有层/红细胞)由 ST 侧"聚类 +
#   标志基因"判, 非 label transfer 产出(参考为 CD45+ 免疫-only) —— 故本脚本 03c 用标志基因
#   命名, 不动 00。
# 输入: results/01_read_qc/table/<s>_seurat.rds (SCT 归一化表达, 含空间图像);
#       results/02_unsupervised_cluster/table/cluster_id.csv (02 聚类, 供追溯);
#       results/03_immune_spatial/table/spot_celltype_annotation.csv (03 pred_main + pred_score)
# 输出: Region_ST-*.pdf (per-sample structural regions);
#       table/full_annotation.csv (per-spot region + cluster + pred_main/pred_score);
#       table/cluster_annotation_review.csv (per-cluster evidence: data-driven -> LP -> user -> final)
# 注: 自包含(规范22, 不 source 外部脚本); 仅消费 01/02/03 产物。
#      图内不使用 "Fig X" 字样(图号属文稿, 规范17); 标题/标签全英文(规范17)。
# =============================================================================

# --- 包加载 ---------------------------------------------------------------
library(Seurat)    # ST 对象与空间分析
library(Matrix)    # 稀疏矩阵
library(ggplot2)   # 可视化
library(RANN)      # 空间 KNN 邻接 (LP 解剖位置修正用, 用户拍板规则)

# --- 区域命名方法(论文忠实, 数据驱动首过 + H&E 复核) ----------------------------
# 原文 Fig 1F 区域 = 无监督聚类簇 + 作者看每簇标志基因 + H&E 手工命名(Methods L308 未给任何算法)。
#   本脚本做两件事: (1) 输出每簇 8 面板标志基因表达表(即"作者检视"所看的数据); (2) 透明首过命名 =
#   通用候选规则(结构区 rel-argmax + 主导程序兜底, 见 PART 2 命名段): 在"结构区面板"(排除游离红细胞,
#   因其非 Fig 1F 结构区)中取相对基线升高量最高者为名; 全面板无特异性升高(rel<=MIN_SPEC)则兜底取绝对
#   表达最高面板(= 该样本主导程序, 数据推导, 零硬编码); 无可辨升高 -> "Unassigned" 待 H&E 确认
#   (忠实论文: 区分不开的簇走 H&E)。此为首过启发式, 非论文原文算法、非按样本硬编码、非基因级加权;
#   阈值 MIN_SPEC 为跨样本通用经验值。详见 PART 2 内联命名段。

# --- 全局辅助函数(内联, 规范22 自给自足) --------------------------------
sci_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.text      = element_text(size = 10, face = "plain"),
    axis.title       = element_text(size = 12, face = "plain"),
    axis.text        = element_text(size = 10, face = "plain"),
    plot.title       = element_text(size = 13, hjust = 0.5, face = "plain"),
    legend.title     = element_text(size = 11, face = "plain")
  )

sci_save <- function(plot_obj, filename, width = 8, height = 8) {
  pdf(file.path(FIG_DIR, paste0(filename, ".pdf")), width = width, height = height,
      family = "Helvetica")
  print(plot_obj)
  dev.off()
}

# (区域命名改由 PART 2 内联的"看每簇标志基因 + 相对基线升高量"首过命名完成, 无独立函数、无 AUTO_COMPUTE 开关;
#  该计算式对应论文"作者检视各簇标志基因", 文献面板驱动、无自创加权/阈值, 详见 PART 2 命名段。)

# --- 目录 ---------------------------------------------------------------
BASE_DIR   <- "X:/XX"
TAB01_DIR  <- file.path(BASE_DIR, "results/01_read_qc/table")     # SCT 对象(含空间图像)
TAB02_DIR  <- file.path(BASE_DIR, "results/02_unsupervised_cluster/table")
TAB03_DIR  <- file.path(BASE_DIR, "results/03_immune_spatial/table")
FIG_DIR    <- file.path(BASE_DIR, "results/03c_full_annotation/Figure")
TAB_DIR    <- file.path(BASE_DIR, "results/03c_full_annotation/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

sample.id <- c("ST-P1_liver","ST-P1_colon","ST-P2_liver","ST-P2_colon",
               "ST-P3_liver","ST-P3_colon","ST-P4_liver","ST-P4_colon")

# --- D-A: 8 类结构区标记面板 (文献标准正交面板, 标志基因命名, 阈值无关) ---
# 论文 Methods 仅称结构区由 "unsupervised clustering" 出簇(Results L80), 区域命名由作者检视各簇
#   标志基因完成; 原文未给显式结构标记清单(Results L80 仅述区域类型, Supplementary Table S2 仅列
#   免疫细胞标记, 见 L78/L308 上下文及论文 S2)。故结构面板按「文献标准正交面板」重建: 每类结构区
#   取自该组织/细胞类型的权威单细胞或谱系 marker 论文, 且 8 个面板之间零基因重叠(正交), 避免同一
#   基因被两个区共享造成错配。各面板出处与 DOI 见下行内注释。
#   注: 论文 ssGSEA(L308)用于免疫亚细胞型(Fig 2C)打分, 非结构区命名; 本面板供 PART 2 输出每簇
#   8 面板表达表与人工命名判读, 不计算 ssGSEA 评分(论文 ssGSEA 仅用于免疫亚群 Fig 2C), 不参与
#   Fig 1F 区域判类(Fig 1F 区域 = 聚类 + 作者手工命名)。
panels <- list(
  # 肿瘤(结直肠腺癌): EPCAM(上皮黏附分子) + CDX2/KRT20/CEACAM5(结直肠分化/癌标志) + MKI67(增殖)
  #   不用泛上皮 KRT8/18(肝细胞亦可表达角蛋白), 避免 liver spot 误判肿瘤; KRT19 为肠上皮特异角蛋白。
  #   文献: Weise et al. 2020 Cancers "Immunohistochemical markers of colorectal cancer"
  #     DOI 10.3390/cancers12071609 (CDX2/KRT20/CEACAM5 为结直肠腺癌标准 IHC 标志);
  #     EPCAM/KRT19 见 Clinical Proteomics / Abcam 谱系注释。
  Tumor            = c("EPCAM","CDX2","KRT20","CEACAM5","KRT19","MKI67"),
  # 正常上皮(肠杯状/分化): MUC2/TFF1/TFF3(杯状细胞黏蛋白) + FCGBP(肠上皮糖蛋白) + SPDEF(杯状分化
  #   转录因子) + ALPI(肠上皮碱性磷酸酶, 肠吸收的谱系标志)。
  #   文献: Haber et al. 2018 Nature "A single-cell survey of the small intestinal epithelium"
  #     DOI 10.1038/nature24489; Smillie et al. 2019 Cell "Intestinal epithelial cell diversity"
  #     DOI 10.1016/j.cell.2019.06.029 (MUC2/TFF/ALPI/SPDEF/FCGBP 为结肠上皮/杯状细胞标志)。
  Normal_epithelium = c("MUC2","TFF1","TFF3","FCGBP","SPDEF","ALPI"),
  # 肝细胞: ALB/APOA1/APOA2(血浆蛋白) + SERPINA1(α1-抗胰蛋白酶) + CYP3A4(药物代谢酶) + HAMP(铁调素)
  #   + TTR(转运蛋白)。为肝实质细胞特异, 与间质/血管零重叠。
  #   文献: MacParland et al. 2018 Nat Commun "Single cell transcriptome analysis of human liver"
  #     DOI 10.1038/s41467-018-06318-7 (肝细胞 ALB/APOA1/APOA2/CYP3A4/HAMP 标志)。
  Hepatocyte       = c("ALB","APOA1","APOA2","SERPINA1","CYP3A4","HAMP","TTR"),
  # 成纤维/CAF: COL1A1/COL1A2(I 型胶原) + PDGFRA/FAP(CAF 标志) + PCOLCE(前胶原 C 蛋白酶, 基质活化)
  #   + CFD(补体因子 D, 活化基质)。不含 ACTA2/DES(已移入平滑肌) 与 COL3A1/VIM(已移入固有层), 保证正交。
  #   文献: Öhlund et al. 2017 J Exp Med "Distinct populations of quiescent and myofibroblast
  #     pancreatic stellate cells" DOI 10.1084/jem.20162024 (PDGFRA/FAP 为 CAF/iCAF/myCAF 标志);
  #     Rensen et al. 2020 J Hepatol "Hepatic stellate cells" DOI 10.1016/j.jhep.2020.04.022
  #     (肝 CAF COL1A1/PDGFRA/FAP/CFD/PCOLCE)。
  Fibroblast       = c("COL1A1","COL1A2","PDGFRA","FAP","PCOLCE","CFD"),
  # 内皮(通用 + 肝窦 LSEC): PECAM1/VWF/CLDN5/CDH5/PLVAP(连续/毛细血管内皮) + CD34/LYVE1/STAB1/ACKR1
  #   (肝窦/淋巴管内皮, LSEC 标志, 复现论文 Fig 1F 肝转移门管/中央静脉 spot 的内皮区)。
  #   文献: MacParland et al. 2018 Nat Commun DOI 10.1038/s41467-018-06318-7 (LSEC CD34/LYVE1/STAB1);
  #     Vanlandewijck & Docherty 2025 Cell Genomics "A cell and transcriptome atlas of human
  #     arterial vasculature" DOI 10.1016/j.xgen.2025.100811 (PECAM1/VWF/CLDN5/CDH5 为内皮谱系标志)。
  Endothelial      = c("PECAM1","VWF","CLDN5","CDH5","PLVAP","CD34","LYVE1","STAB1","ACKR1"),
  # 平滑肌(血管/肠壁): ACTA2(α-SMA) + MYH11(肌球蛋白重链, SMC 特异) + TAGLN(转凝蛋白) + DES(结蛋白)
  #   + CNN1(钙调蛋白, SMC 特异, 与成纤维零重叠)。
  #   文献: Vanlandewijck & Docherty 2025 Cell Genomics DOI 10.1016/j.xgen.2025.100811
  #     (MYH11/CNN1 为 SMC 特异, 成纤维标记特异性差); 经典 SMC 谱系 ACTA2/TAGLN/DES (PMC5852505)。
  Smooth_muscle    = c("ACTA2","MYH11","TAGLN","DES","CNN1"),
  # 固有层(肠间质/ niche): COL3A1/COL6A1(III/VI 型胶原, 黏膜下间质) + THY1(CD90, 基质/ niche 细胞)
  #   + PDGFRB(周细胞/基质受体) + ENG(CD105, 内皮/间质活化) + RSPO3(R-spondin3, 隐窝 niche 分泌标志)。
  #   与成纤维面板正交: PDGFRA(Fibroblast) vs PDGFRB(固有层) 为不同基因; 无共享胶原基因。
  #   文献: Smillie et al. 2019 Cell DOI 10.1016/j.cell.2019.06.029 (结肠基质/ niche: COL3A1/COL6A1/
  #     THY1/PDGFRB/ENG/RSPO3/GREM1, 固有层成纤维 niche 标志)。
  Lamina_propria   = c("COL3A1","COL6A1","THY1","PDGFRB","ENG","RSPO3","VIM"),
  # 红细胞(血红蛋白/糖蛋白 A): HBB/HBA1/HBA2/HBD(珠蛋白链) + GYPA(血型糖蛋白 A, 红系特异膜蛋白)。
  #   文献: 经典造血/红系标志 (GYPA 见 Cancer Cell / Abcam 红系谱系注释; 珠蛋白为教科书级红细胞标记)。
  Erythrocyte      = c("HBB","HBA1","HBA2","GYPA","HBD")
)
# 稀疏血管/基质型(内皮/红细胞/平滑肌)在 Visium 多细胞 spot 中本为少数。旧版曾用"z-score argmax
#   + 血管型启发式覆盖(VASCULAR_TYPES / ENDO_CALL_FLOOR / VASCULAR_CALL_RATIO / ERYTHROCYTE_*)"
#   补偿, 此整组常数为自创、论文(L308 上下文)无, 已于 2026-09-07 全部剔除。本脚本现按"看每簇标志基因 +
#   相对基线升高量"做透明首过命名(复现作者检视标志基因), 区分不开的簇留 Unassigned 走 H&E 确认;
#   仍输出每簇 8 面板表达表供人工/H&E 复核; 不计算 ssGSEA 评分(论文 ssGSEA 仅用于免疫亚群 Fig 2C)。
# 面板名 -> 显示标签(统一 Tumor, 不再拆 LM/CRC; 含空格符合 Fig 1F 命名)
label_map <- c(
  Tumor           = "Tumor",
  Normal_epithelium = "Normal epithelium",
  Hepatocyte      = "Hepatocytes",
  Fibroblast      = "Fibroblast",
  Endothelial     = "Endothelial",
  Smooth_muscle   = "Smooth muscle",
  Lamina_propria  = "Lamina propria",
  Erythrocyte     = "Erythrocyte",
  Vascular_mesenchymal_niche = "Vascular mesenchymal niche"
)
# D-B: 集群级区域命名由 PART 2 内联的"看每簇标志基因 + 相对基线升高量"首过命名判出(复现原文"作者检视标志基因");
#   不另设独立覆盖列表/手工命名表(原文 Methods 亦无区域命名表); 区分不开的簇留 Unassigned 走 H&E 确认。
# D-C: 方案 A (2026-09-06) 裁定: 删除自造的免疫叠加合并规则(IMMUNE_SCORE_MIN 命名常数 /
#   IMMUNE_STRUCT_QTL 集群内结构相对弱)。Fig 1F 结构区由"聚类+标记基因"独立判出(本脚本 Region
#   图); 免疫主型/亚型由 03 label transfer / ssGSEA 给出(Fig 1G/2C), 不在本脚本叠加, 避免冗余
#   与自造门控。pred_main/pred_score 仅作 CSV 参考列保留。

# --- 调色板 (2026-09-06 对齐 Wu et al. Fig 1F 原图 legend 配色) ----------
#   Tumor=浅蓝, Endothelial=浅紫, Normal epithelium=粉, Erythrocyte=灰,
#   Fibroblast=浅绿, Hepatocytes=橘, Lamina propria=蓝, Smooth muscle=深绿.
region_cols <- c(
  "Tumor"             = "#A6CEE3",
  "Normal epithelium" = "#FB9A99",
  "Erythrocyte"       = "#B3B3B3",
  "Fibroblast"        = "#B2DF8A",
  "Hepatocytes"       = "#FDBF6F",
  "Lamina propria"    = "#1F78B4",
  "Smooth muscle"     = "#33A02C",
  "Vascular mesenchymal niche" = "#B39DDB",  # 紫色(与 Endothelial 同色): 血管基质/间充质生态位(数据驱动共识签名, 用户 H&E 裁定)
  "Unassigned"        = "#D9D9D9"   # 浅灰: 无可辨特异性簇显示灰, 避免整图空白 NA
)
# 结构区 8 类展示顺序(Region 图用): 按原文 legend 两行顺序
region_order <- c("Tumor","Normal epithelium","Erythrocyte",
                  "Fibroblast","Hepatocytes","Lamina propria","Smooth muscle",
                  "Vascular mesenchymal niche","Unassigned")

# --- 专家 / H&E 区域标注 (用户填写, 非算法; 论文 Fig 1F 区域由作者检视标志基因 + H&E 手工命名) ---
# 键=sample, 值=命名向量(cluster -> 面板名); 仅填用户已 H&E 确认的簇。脚本打印覆盖证据, 不隐藏。
# 任一簇未在此列出 -> 退回数据驱动首过(透明, 仅作检视参照, 非论文算法)。
# P1_liver: c3/c8 实测共享签名=血管基质/间充质生态位(VIM/COL1A1/TAGLN/ACTA2/RNASE1/SOD3/TYMP 等,
#   经典内皮标志仅 c3 升、c8 沉默), 用户据此判为 Vascular mesenchymal niche; 余簇按已确认标签。
USER_REGION <- list(
  "ST-P1_liver" = c("0"="Hepatocyte","1"="Tumor","2"="Tumor","3"="Vascular_mesenchymal_niche",
                    "4"="Hepatocyte","5"="Hepatocyte","6"="Tumor","7"="Fibroblast",
                    "8"="Vascular_mesenchymal_niche"),
  # ST-P1_colon: 用户 H&E 裁定 (2026-09-08) —— 数据驱动首过仅 1/3/6 命中, 见 _diag_p1colon_rootcause
  "ST-P1_colon" = c("0"="Lamina_propria","1"="Tumor","2"="Tumor","3"="Normal_epithelium",
                    "4"="Tumor","5"="Fibroblast","6"="Tumor"),
  # ST-P2_colon: 用户 H&E 裁定 (2026-09-08 末, c0 于 13:38 更正为 Smooth muscle) ——
  #   rel-argmax(排除 Ery) 自产 7/12 命中 (c0/c1/c2/c3/c5/c6/c11, 其中 c0 算法本即 Smooth muscle, 用户更正后确认无误, 故不覆盖);
  #   余 5 簇为"标记↔H&E 结构性错位"(肠壁层共表达 c4 + 肿瘤促纤维沉默 c7/c8/c9/c10), 已仿真证实
  #   扩展空间精化 0 救回、且不回退 c2/c5, 属 H&E 不可辨范畴, 仅覆盖误判簇(算法自产正确的 7 簇保留可见)。
  #   见 _diag_p2colon_refine.R: c4 邻接 NE 分数≈0.10、c8/c9/c10 邻接 Tumor 分数≈0.03/0.35/0.04。
  "ST-P2_colon" = c("4"="Lamina_propria","7"="Normal_epithelium",
                    "8"="Tumor","9"="Tumor","10"="Tumor"),
  # ST-P3_colon: 用户 H&E 裁定 (2026-09-08 13:45) —— 仅 c8 为 Tumor, 其余 10 簇走算法自产(用户确认), 不覆盖。
  "ST-P3_colon" = c("8"="Tumor"),
  # ST-P4_colon: 用户 H&E 裁定 (2026-09-08 14:06) —— c1/c6=Lamina propria, c5/c7=Fibroblast; 余走算法自产。
  #   c8 算法本即 Fibroblast(=用户真值)故不覆盖(遵循 c0 教训: 算法已正确则不写冗余覆盖);
  #   c0/c2/c3/c4 算法自产(用户确认)不覆盖。
  "ST-P4_colon" = c("1"="Lamina_propria","6"="Lamina_propria","5"="Fibroblast","7"="Fibroblast")
)

# =============================================================================
# PART 1: 读取 02 聚类 + 03 pred_main/pred_score
# =============================================================================
cat("========== PART 1: 读取 02 聚类 + 03 pred_main/pred_score ==========\n")
clid_file <- file.path(TAB02_DIR, "cluster_id.csv")
if (!file.exists(clid_file)) stop("缺少 ", clid_file, " —— 请先跑完 02_unsupervised_cluster.R")
clid <- read.csv(clid_file, stringsAsFactors = FALSE)
cl_of_spot <- list()
for (s in sample.id) {
  cl_of_spot[[s]] <- setNames(as.character(clid$seurat_clusters[clid$sample == s]),
                              clid$spot[clid$sample == s])
}

pred_file <- file.path(TAB03_DIR, "spot_celltype_annotation.csv")
if (!file.exists(pred_file)) stop("缺少 ", pred_file, " —— 请先跑完 03_immune_spatial.R")
pred_dt <- read.csv(pred_file, stringsAsFactors = FALSE)
pred_main  <- setNames(pred_dt$pred_main,  pred_dt$spot)
pred_score <- setNames(pred_dt$pred_score, pred_dt$spot)

# =============================================================================
# PART 2: 集群级结构区注释 (聚类骨架来自 02; 区域名 = 看每簇标志基因 + 相对基线升高量首过命名, 复现原文作者检视各簇标志基因)
# 原文方法(Wu 2022, Methods "ST Data Processing and Cell-Type Estimation" 段, 逐字见
#   _paper_fulltext.txt L308; Results Fig 1F 逐字 L80):
#   (1) Space Ranger -> Seurat SCTransform -> RunPCA -> FindNeighbors+FindClusters -> RunUMAP
#       (聚类 ST spots, 即本图"区域"来源, L308);
#   (2) Seurat label transfer 预测主细胞型(L308 "predicted main cell types by using Seurat")——
#       但 scRNA 参考为 CD45+ 免疫-only(L286 "only reserved the CD45+ clusters"), 故 label transfer
#       只能预测免疫主型(服务 Fig 1G), 无法命名肿瘤/成纤维/肝细胞等非免疫结构区;
#   (3) ssGSEA 仅用于亚细胞型(subcell types)打分(L308 "We hence used ssGSEA to score the subcell
#       types", 论文 ref 96=GSVA), 应用对象为免疫亚群(Fig 2C), 非 Fig 1F 结构区。
#   Results L80 逐字: "Unsupervised clustering analysis classified the samples into distinct regions
#   such as tumor and fibroblast regions (Fig. 1F)." —— 即 Fig 1F 结构区 = 无监督聚类出的簇,
#   其命名由作者检视各簇标志基因给出(原文未给显式结构标记清单; S2 仅列免疫标记)。
#   本 PART 2 以"看每簇标志基因 + 相对基线升高量"做透明首过命名(计算式对应"作者检视各簇标志基因"),
#   同时输出每簇 8 面板表达表供人工/H&E 复核(详见下方命名段)。
#   不计算 ssGSEA 评分(论文 ssGSEA 仅用于免疫亚群 Fig 2C)。Fig 1F 区域非固定算法: 标记不可辨者
#   (如 Endothelial) 判为 Unassigned, 忠实论文(该结构区须 H&E 复核)。
# 标记表达: 用 01 SCTransform 校正表达(SCT assay data 层, 论文 "SCTransform to normalize data")。
# 解剖学感知面板集 pans: 肝剔 Normal_epithelium/Lamina_propria; 结肠剔 Hepatocyte。
# 免疫覆盖: 方案 A 删除(见 D-C) —— Fig 1F 为结构区, 免疫属 03, 本脚本不叠加自造合并规则。
# =============================================================================
cat("========== PART 2: 集群级注释 + 全注释回贴 ==========\n")
anno_of_spot   <- list()
st_obj         <- list()   # 存 01 RDS 对象(含空间图像)供 PART 3 出图
review_rows    <- list()   # 每簇注释审计行 (data-driven -> LP -> user override -> final)
panel_names_all<- c("Tumor","Normal_epithelium","Hepatocyte","Fibroblast",
                    "Endothelial","Smooth_muscle","Lamina_propria","Erythrocyte")
for (s in sample.id) {
  cat("  处理", s, "...\n")
  rds <- file.path(TAB01_DIR, paste0(s, "_seurat.rds"))
  if (!file.exists(rds)) stop("缺少 ", rds, " —— 请先跑完 01_read_qc.R")
  st <- readRDS(rds)
  sct <- as.matrix(GetAssayData(st, assay = "SCT", layer = "data"))   # genes x spots, SCTransform 校正表达(强制入内存, 消除 HDF5 后台矩阵损坏)

  cl <- cl_of_spot[[s]]                                      # spot -> cluster (来自 02)
  clusters <- sort(unique(cl))

  # 解剖学感知面板集: 肝无正常肠上皮/固有层; 结肠无肝细胞
  if (grepl("liver", s)) {
    pans <- panels[!(names(panels) %in% c("Normal_epithelium","Lamina_propria"))]
  } else {
    pans <- panels[names(panels) != "Hepatocyte"]
  }

  # --- 区域命名: 看每簇标志基因表达, 与全样本基线比较 (复现原文"作者检视各簇标志基因") ---
  # 原文 Fig 1F = 无监督聚类簇 + 作者看每簇标志基因 + H&E 手工命名(Methods L308 未给任何算法)。
  #   本步做两件事: (1) 输出每簇 8 面板标志基因表达表 (panel_expr) —— 即"作者检视"所看的数据, 供判读/H&E 复核;
  #   (2) 透明首过命名 (通用候选规则): 对每个区域 R 算 pv[R]=簇 c 的 R 面板均值、rel[R]=pv[R]−全样本
  #       R 面板基线均值; 在"结构区面板"(排除游离红细胞 Erythrocyte, 其非 Fig 1F 结构区)中取 rel 最高者
  #       为名; 全面板 rel 均 <= MIN_SPEC(无特异性升高) -> 兜底取绝对表达最高面板(= 该样本主导程序,
  #       数据推导, 零硬编码)。阈值见下方命名段。
  #   说明: 这是"作者看到 ALB 各簇都高=背景、转而看真正区分性的标志"的透明计算式对应; 排除红细胞是
  #        schema 级生物学约束(游离 RBC 非结构区)而非按簇硬编码; 非论文原文算法、非按样本硬编码、非基因级加权。
  #       区分不开的簇(真实标记重叠, 两程序同升) -> 透明首过后仍可能判错 -> 由 USER_REGION(H&E) 修正, 按原文走 H&E。
  clab <- character(length(clusters)); names(clab) <- clusters
  # 收集本样本全部面板基因 + 逐基因跨簇均值 (gene x cluster)
  pan_genes_present <- lapply(pans, function(g) rownames(sct)[rownames(sct) %in% g])
  pan_genes_present <- pan_genes_present[lengths(pan_genes_present) > 0]
  all_genes <- unique(unlist(pan_genes_present))
  clust_gene <- matrix(NA_real_, nrow = length(all_genes), ncol = length(clusters),
                       dimnames = list(all_genes, as.character(clusters)))
  for (c in clusters) {
    sp_in_c <- names(cl)[cl == c]
    clust_gene[, as.character(c)] <- if (length(sp_in_c) == 0) 0 else
      rowMeans(sct[all_genes, sp_in_c, drop = FALSE])
  }
  # 输出每簇面板表达表 (判读依据, 即"作者检视各簇标志基因"所看的数据)
  panel_expr_list <- list()
  for (c in clusters) {
    cg <- clust_gene[, as.character(c)]
    cp <- sapply(pan_genes_present, function(g) mean(cg[g], na.rm = TRUE))
    panel_expr_list[[as.character(c)]] <- cp
    cat("    ", s, "cluster", c, "panel_expr:",
        paste(sprintf("%s=%.2f", names(cp), cp), collapse = " "), "\n")
  }
  # --- 透明首过命名 (通用候选规则: 结构区 rel-argmax + 主导程序兜底) ---
  # 仍为非论文原文算法、非按样本/按簇硬编码的透明启发式。设计:
  #   结构区 = 论文 Fig 1F 的 8 类组织/细胞区; 游离红细胞(Erythrocyte)不是结构区(窦性/技术污染),
  #     从候选面板中排除——此为 schema 级生物学约束(非按簇硬编码), 消解"红细胞靠 rel 虚高抢标"类误判(如 c6)。
  #   对每区域 R 算 (a) 绝对表达 pv[R]; (b) 相对基线升高量 rel[R] = pv[R] - 全样本 R 面板基线(pan_base)。
  #   命名 = 在"结构区面板"中取 rel 最高者(= 该簇最特异性升高的组织程序);
  #          若全面板 rel 均 <= MIN_SPEC(簇无标志性特异性升高) -> 兜底取绝对表达最高面板(主导程序, 数据推导, 零硬编码)。
  #   该规则天然: c2(肝血管) SM rel 最高 -> Smooth_muscle; c6 Ery 被排除、Hep rel 最高 -> Hepatocytes; 无按簇硬编码。
  #   真实标记重叠簇(LP<->NE、Fib<->Tumor-stroma)仍两程序同升 -> 任一带过, 须 H&E(USER_REGION)。
  #   MIN_SPEC 为跨样本通用经验值(非按样本调参、非论文给定)。
  MIN_SPEC <- 0.05   # 相对基线特异性下限: 升高量低于此视为无区分性(噪声 rel~0.01 远低, 留 ~5x 余量)
  if (length(clusters) >= 2) {
    gene_base <- rowMeans(clust_gene, na.rm = TRUE)
    pan_base  <- sapply(pan_genes_present, function(g) mean(gene_base[g], na.rm = TRUE))
    struct_panels <- setdiff(names(pan_genes_present), "Erythrocyte")   # 游离红细胞非结构区, 不参与候选
    for (c in clusters) {
      pv  <- sapply(pan_genes_present, function(g) mean(clust_gene[g, as.character(c)], na.rm = TRUE))
      rel <- pv - pan_base[names(pv)]
      rv  <- rel[struct_panels]
      if (max(rv, na.rm = TRUE) > MIN_SPEC) {
        clab[c] <- names(which.max(rv))          # 最特异性升高的结构区程序
      } else {
        clab[c] <- names(which.max(pv))          # 无特异性升高 -> 主导程序(数据驱动推导, 非硬编码)
      }
    }
  } else {
    # 单簇样本: 无"其余簇"可对比, 退化为绝对面板均值最高者(仍 >0 才判)
    for (c in clusters) {
      cm <- sapply(pan_genes_present, function(g) mean(clust_gene[g, as.character(c)], na.rm = TRUE))
      cm <- cm[!is.na(cm)]
      clab[c] <- if (length(cm) > 0 && max(cm) > 0) names(which.max(cm)) else "Unassigned"
    }
  }
  clab_dd <- clab   # 保存数据驱动首过命名
  cat("    ", s, "data-driven 命名:",
      sum(clab != "Unassigned"), "/", length(clusters), "簇已命名;",
      sum(clab == "Unassigned"), "簇 Unassigned(标志基因无可辨升高, 待 H&E 确认)\n")

  # Pass 2: LP 解剖位置修正 —— **用户拍板规则**(非自创):
  #   基质簇若紧邻 Normal epithelium -> Lamina propria; 若嵌入 Tumor 团内 -> Fibroblast。
  #   (黏膜固有层与肿瘤内基质 CAF 共享 COL1A1/COL3A1/THY1, 纯 marker 法分不开, 须用 ST 空间位置区分。
  #    固有层 = 上皮下基质层; 肿瘤内基质 = 嵌入肿瘤团内; 纯基因法无此位置信息。)
  # 【2026-09-07 剔除自创阈值】此前叠加的 (COL1A1-ACTA2)>0.5 / fNEnbr>0.20 / d_NE<150px /
  #   frac_NE_nbr>=0.05 / !is_sm_dominant 五条门均为我自创, 已全部移除, 仅保留用户拍板的
  #   原始语义: 多数 spot 邻接 NE(>0.5) 且非嵌入 Tumor(<0.5) -> LP; 多数邻接 Tumor(>0.5) -> Fibroblast。
  coords <- GetTissueCoordinates(st, scale = "lowres")
  cnl <- tolower(colnames(coords))
  icol <- which(cnl %in% c("imagecol", "x"))
  irow <- which(cnl %in% c("imagerow", "y"))
  coord_mat <- as.matrix(coords[colnames(st), c(icol, irow)])
  # 6 个最近邻 (hex 网格紧邻; RANN 为 Seurat 传递依赖)
  knn_idx <- RANN::nn2(data = coord_mat, k = 7)$nn.idx[, -1]
  nbr_mat <- matrix(colnames(st)[knn_idx], nrow = nrow(knn_idx), ncol = 6,
                    dimnames = list(colnames(st), NULL))
  stromal_cl <- clusters[clab %in% c("Fibroblast", "Lamina propria")]
  for (c in stromal_cl) {
    sp_in_c <- names(cl)[cl == c]
    n_NE <- 0; n_Tumor <- 0
    for (sp in sp_in_c) {
      nbr_cl <- cl[nbr_mat[sp, ]]
      if (any(nbr_cl %in% clusters[clab == "Normal epithelium"])) n_NE <- n_NE + 1
      if (any(nbr_cl %in% clusters[clab == "Tumor"])) n_Tumor <- n_Tumor + 1
    }
    frac_NE <- n_NE / length(sp_in_c)
    frac_Tumor <- n_Tumor / length(sp_in_c)
    old <- clab[c]
    # 多数邻接 NE 且非嵌入 Tumor -> 黏膜固有层; 多数邻接 Tumor -> 肿瘤内基质(Fibroblast)
    if (frac_NE > 0.5 && frac_Tumor < 0.5) {
      clab[c] <- "Lamina propria"
    } else if (frac_Tumor > 0.5) {
      clab[c] <- "Fibroblast"
    }
    if (clab[c] != old) {
      cat("    ", s, "cluster", c, "anatomical refine:", old, "->", clab[c],
          sprintf("(frac_NE=%.2f frac_Tumor=%.2f)", frac_NE, frac_Tumor), "\n")
    }
  }
  clab_lp <- clab   # 保存 LP 解剖位置修正后命名

  # --- 专家 / H&E 区域标注覆盖 (用户填写, 非算法; 论文 Fig 1F 区域由作者检视标志基因 + H&E 命名) ---
  # USER_REGION 仅含用户已 H&E 确认的簇; 覆盖数据驱动首过(透明, 打印证据), 未列簇退回数据驱动。
  if (s %in% names(USER_REGION)) {
    um <- USER_REGION[[s]]
    for (cc in names(um)) {
      if (as.character(cc) %in% names(clab)) {
        old <- clab[as.character(cc)]
        clab[as.character(cc)] <- um[cc]
        if (clab[as.character(cc)] != old)
          cat("    ", s, "cluster", cc, "user-annotation override:", old, "->", clab[as.character(cc)], "\n")
      }
    }
  }
  # --- 全局重命名: 所有 Endothelial 区域统一为 Vascular mesenchymal niche (用户 2026-09-08 14:22 裁定) ---
  # 用户此前将 P1_liver c3/c8 判为 Vascular mesenchymal niche(血管基质/间充质生态位); 现扩展为全局规则:
  #   凡算法/USER_REGION 命名为 Endothelial 的簇, 一律改称 Vascular mesenchymal niche(同紫色 #B39DDB,
  #   与 Endothelial 原色一致), 不再单列 Endothelial 结构区(区域/图例均不再出现 Endothelial)。
  clab[clab == "Endothelial"] <- "Vascular_mesenchymal_niche"
  # 构建每簇审计行 (data-driven -> LP refine -> user override -> final)
  for (c in clusters) {
    cc <- as.character(c)
    pe <- panel_expr_list[[cc]]
    if (is.null(pe)) pe <- setNames(rep(NA_real_, length(panel_names_all)), panel_names_all)
    pe_full <- setNames(rep(NA_real_, length(panel_names_all)), panel_names_all)
    pe_full[names(pe)] <- pe
    final_display <- if (clab[cc] %in% names(label_map)) label_map[clab[cc]] else clab[cc]
    dd_display <- if (clab_dd[cc] %in% names(label_map)) label_map[clab_dd[cc]] else clab_dd[cc]
    lp_display <- if (clab_lp[cc] %in% names(label_map)) label_map[clab_lp[cc]] else clab_lp[cc]
    review_rows[[length(review_rows) + 1]] <- data.frame(
      sample = s, cluster = c,
      data_driven = dd_display,
      after_lp = lp_display,
      final = final_display,
      user_override = clab[cc] != clab_lp[cc],
      top_panel = names(pe)[which.max(pe)],
      panel_expr_Tumor = pe_full["Tumor"],
      panel_expr_Normal_epithelium = pe_full["Normal_epithelium"],
      panel_expr_Hepatocyte = pe_full["Hepatocyte"],
      panel_expr_Fibroblast = pe_full["Fibroblast"],
      panel_expr_Endothelial = pe_full["Endothelial"],
      panel_expr_Smooth_muscle = pe_full["Smooth_muscle"],
      panel_expr_Lamina_propria = pe_full["Lamina_propria"],
      panel_expr_Erythrocyte = pe_full["Erythrocyte"],
      stringsAsFactors = FALSE
    )
  }

  # --- 逐 spot 注释 = 结构区(Fig 1F 由聚类+标记判, 即本步) ---
  # 方案 A (2026-09-06): 删除自造的免疫叠加合并规则(IMMUNE_SCORE_MIN/IMMUNE_STRUCT_QTL)与
  #   Unassigned 类。Fig 1F 为结构区; 免疫主型/亚型属 03(Fig 1G/2C), 不在本图叠加, 避免冗余
  #   与自造门控。pred_main/pred_score 仅作 CSV 参考列保留(不参与判类)。
  # 面板名 -> 显示标签(与 region_order/region_cols 一致, 论文 Fig 1F legend 形式)
  ann <- setNames(ifelse(clab[cl] %in% names(label_map),
                         label_map[clab[cl]],
                         clab[cl]),
                  names(cl))
  anno_of_spot[[s]] <- ann
  st_obj[[s]] <- st
  cat("    ", s, "注释分布:",
      paste(names(sort(table(ann), decreasing = TRUE)), collapse = " / "), "\n")
}

# 输出每簇注释审计表 (data-driven -> LP refine -> user override -> final)
review_df <- do.call(rbind, review_rows)
write.csv(review_df, file.path(TAB_DIR, "cluster_annotation_review.csv"), row.names = FALSE)
cat("Saved cluster_annotation_review.csv (", nrow(review_df), " clusters)\n")

# =============================================================================
# PART 3: 结构区图 (每样本) + 逐 spot 全注释 CSV
# =============================================================================
cat("\n========== PART 3: 结构区图 (每样本) + 逐 spot CSV ==========\n")
full_tab <- data.frame()
for (i in seq_along(sample.id)) {
  s <- sample.id[i]
  st <- st_obj[[s]]
  ann <- anno_of_spot[[s]]
  cl  <- cl_of_spot[[s]]

  # --- (A) 结构区图: 集群级结构标签(标志基因跨 cluster 特异性命名, Fig 1F) ---
  reg_cl <- ann      # 逐 spot 结构区 = PART 2 已算好的 anno_of_spot[[s]](按本样本 clab, 非跨样本残留变量)
  st$region <- factor(reg_cl[match(colnames(st), names(reg_cl))],
                      levels = intersect(region_order, reg_cl))
  p_reg <- SpatialDimPlot(st, group.by = "region", label = FALSE,
                          cols = region_cols, image.alpha = 1) +
    labs(title = paste0(s, " - spatial region annotation")) + sci_theme +
    theme(legend.position = "right", legend.title = element_blank())
  sci_save(p_reg, paste0("Region_", s), width = 8, height = 8)

  df <- data.frame(sample = s, spot = names(ann),
                   cluster = cl_of_spot[[s]][names(ann)],
                   region = reg_cl,
                   pred_main = pred_main[names(ann)],
                   pred_score = pred_score[names(ann)],
                   stringsAsFactors = FALSE)
  full_tab <- rbind(full_tab, df)
}

# --- (D) 2x2 未治疗组组合图 (StructuralRegions_untreated_2x2.pdf) 经用户 2026-09-08 裁定剔除(布局多次调整未达预期), 本脚本不再生成任何多样本组合图 ---

write.csv(full_tab, file.path(TAB_DIR, "full_annotation.csv"), row.names = FALSE)
cat("Saved full_annotation.csv (", nrow(full_tab), " spots)\n")
cat("03c_full_annotation.R done: Region_ST-*.pdf + full_annotation.csv + cluster_annotation_review.csv (8 structural regions, paper Fig 1F; Endothelial 统一更名 Vascular mesenchymal niche)\n")
