# ===============================================================
# 11_PLAUR_SPP1_communication.R
# 功能: 用 CellChat 量化 PLAU-PLAUR 轴在 PT 与 LM 的通讯强度, 并用
#       nichenetr 校验该配体-受体对与 SPP1 的调控关系
# 定位: 方案 Fig 4 + H3(PLAU-PLAUR 轴通讯强度, 原发 PT vs 肝转移 LM;
#       实际显著方向为髓系自分泌 Myeloid->Myeloid, 见 PART 2 结果, 与链
#       的"肿瘤分泌 PLAU->招募髓系"方向相反, 属对链的细化/补充而非确认)
# 可视化: CellChatPLAURPathway_PT.pdf / CellChatPLAURPathway_LM.pdf (气泡图)
#         CellChatComparePathway.pdf (PT vs LM 全通路比较, rankNet)
#         CellChatPLAUR_Heatmap.pdf (区室×方向热图, 补充)
#         CellChatPLAUR_PTvLM.pdf (PT vs LM 聚焦比较 H3, 补充)
# 数据: 输入 results/01_read_qc/table/ST-*_seurat.rds
#            results/11_PLAUR_SPP1_communication/table/signature_master.rds
#       输出 results/11_PLAUR_SPP1_communication/table/CellChatPLAUR.csv
#            results/11_PLAUR_SPP1_communication/table/NicheNetLigandRank.csv
# 依据: 方案第 5 节(CellChat + NATMI + NicheNet)与第 7.1 节 Fig 4
# 说明(重要前提): 03 的 label transfer 参考为 CD45+ 免疫-only(见 03 脚注),
#   故 pred_main 只有免疫主型、没有上皮/肿瘤型。本脚本用 v1.1 数据驱动髓系
#   掩码 is_macro 构造复合标签: 髓系→"Myeloid"(PLAUR 接收方), 肿瘤区非髓系
#   →"Tumor"(PLAU 发送方, 上皮), 其余沿用 pred_main。
#   这是 ST 数据上做细胞通讯的务实构造, 已在方案风险 R1 中登记。
# 官方 API: CellChat 1.6.1 官方 vignette(createCellChat / subsetData /
#   identifyOverExpressedInteractions / computeCommunProb / aggregateNet);
#   nichenetr 2.2.1.1 提供 lr_network(配体-受体)与 ligand_target_matrix
#   (行为配体、列为靶基因的调控潜力矩阵)
# ===============================================================

# --- 包加载(统计 -> 数据 -> 可视化 -> 领域) -----------------------------
suppressMessages({
  library(Seurat)      # 空间对象读取
  library(CellChat)    # 细胞间通讯推断
  library(nichenetr)   # NicheNet 配体-靶基因先验
  library(ggplot2)     # 可视化
  library(patchwork)   # 多图拼接
  library(dplyr)       # 数据操作
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
TAB01    <- file.path(BASE_DIR, "results/01_read_qc/table")
TAB06    <- file.path(BASE_DIR, "results/06_PLAUR_SPP1_signature/table")
FIG_DIR  <- file.path(BASE_DIR, "results/11_PLAUR_SPP1_communication/Figure")
TAB_DIR  <- file.path(BASE_DIR, "results/11_PLAUR_SPP1_communication/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 分析常量 -----------------------------------------------------------
SAMPLE_ID <- c("ST-P1_liver", "ST-P1_colon", "ST-P2_liver", "ST-P2_colon")
# 只比较两个含肿瘤的区室: PT(原发结直肠肿瘤) 与 LM(肝转移)
COMP_RUN  <- c("PT", "LM")
MIN_CELLS <- 10   # CellChat 要求每组最少 spot 数
SEED      <- 42

# 项目配色(与 09/10 一致): 髓系蓝 / 肝转移橙
COL_MYE  <- "#1F78B4"
COL_META <- "#FD8D3C"

# =============================================================================
# PART 1: 构造 CellChat 输入(复合标签 + SCT 表达)
# 目的: 生成每个区室的 spot 级表达矩阵与细胞标签
# 方法: 逐样本读 Seurat, 用 FetchData(assay="SCT", layer="data") 取 log1p
#       表达(与 06 签名同尺度; 本队列无 RNA assay, merge 后 SCT 退化为 0 细胞,
#       故逐样本取数后列向拼接为 基因×spot); 按 signature_master 的 is_macro /
#       region / pred_main 构造复合标签; 按区室切分
# 数据: 中间产物 expr.list / meta.list(基因×spot 矩阵 + spot 级 meta)
# =============================================================================

cat("========== PART 1: Build CellChat input ==========\n")

master <- readRDS(file.path(TAB06, "signature_master.rds"))
master <- master[!is.na(master$compartment), ]

# 复合标签(CellChat 分组): 用 v1.1 数据驱动髓系掩码 is_macro 区分肿瘤上皮与髓系,
#   否则肿瘤区非髓系 spot 全塌缩为 "Tumor"、PT/LM 单组→CellChat 无法推断通讯。
#   - 髓系(is_macro): "Myeloid"(PLAUR 接收方)
#   - 肿瘤区非髓系(region=="Tumor" & !is_macro): "Tumor"(PLAU 发送方, 上皮)
#   - 其余(癌旁非髓系): 沿用 pred_main 免疫主型(本脚本 CellChat 仅跑 PT/LM, 不影响)
#   NA 归 Other
master$labels <- ifelse(master$is_macro, "Myeloid",
                  ifelse(master$region == "Tumor", "Tumor", master$pred_main))
master$labels <- ifelse(is.na(master$labels), "Other", master$labels)

set.seed(SEED)
options(stringsAsFactors = FALSE)
CellChatDB <- CellChatDB.human
# 本环境 CellChatDB.human 槽名为 geneInfo(符号列 Symbol) / interaction(LR 对),
#   非旧版 gene / id。基因列表取 geneInfo$Symbol(符号), 与 Seurat 特征对齐
genes.use <- unique(CellChatDB$geneInfo$Symbol)

# --- 扩展 CellChatDB: 补入 PLAU-PLAUR 特异位点对(H3 直接证据) ---
# 依据(文献充分证实): PLAU(uPA, 分泌型) 与 PLAUR(uPAR, 膜受体) 是经典
#   配体-受体对, 介导 ECM 降解与肿瘤细胞入侵。本环境 CellChatDB.human 的
#   interaction 表(1939 对)不含该对(PLAU %in% interaction$ligand == FALSE),
#   故手动补入以直接计算 H3 通讯概率。PLAU/PLAUR 均已在 geneInfo$Symbol
#   中(已核实 genes.use), 无需改 geneInfo; 向 interaction 追加一行即扩库
#   (标准 CellChat 操作)。
if (!any(CellChatDB$interaction$ligand == "PLAU" &
         CellChatDB$interaction$receptor == "PLAUR")) {
  plaur_row <- data.frame(
    interaction_name   = "PLAU_PLAUR",
    pathway_name       = "PLAU_PLAUR",
    ligand             = "PLAU",
    receptor           = "PLAUR",
    agonist            = "",
    antagonist         = "",
    co_A_receptor      = "",
    co_I_receptor      = "",
    evidence           = "Curated: uPA(PLAU)-uPAR(PLAUR) canonical ligand-receptor pair",
    annotation         = "Secreted Signaling",
    interaction_name_2 = "PLAU - PLAUR",
    stringsAsFactors   = FALSE)
  rownames(plaur_row) <- "PLAU_PLAUR"   # 显式行名=交互名, 使 LRsig/net$prob 第三维键一致
  CellChatDB$interaction <- rbind(CellChatDB$interaction, plaur_row)
  cat("已向 CellChatDB$interaction 补入 PLAU-PLAUR 对; 现有", nrow(CellChatDB$interaction), "对\n")
} else {
  cat("CellChatDB 已含 PLAU-PLAUR 对, 跳过补入\n")
}

# 不 merge Seurat 对象: 本队列 Seurat 仅有 Spatial/SCT assay(无 RNA assay),
# 且 SCT/Spatial 跨样本 merge 后 assay 退化为 0 细胞。改为逐样本用安全
# 访问器 FetchData(assay="SCT", layer="data") 取表达(与 06 签名定义一致,
# log1p 尺度), 列向拼接为 基因×spot 矩阵后直接送 CellChat
expr.list <- list(); meta.list <- list()
for (s in SAMPLE_ID) {
  st <- readRDS(file.path(TAB01, paste0(s, "_seurat.rds")))
  # 仅取 master 中本样本且存在于 Seurat 的 spot: master 经 compartment 过滤
  #   (line 82), 可能少于 Seurat 全部 spot; 二者取交集, 才能保证 矩阵列数 ==
  #   meta 行数(否则 colnames<- 报 array extent 错)
  sub0   <- master[master$sample == s, ]
  bc_use <- intersect(sub0$spot, Cells(st))
  sub    <- sub0[sub0$spot %in% bc_use, ]
  bc     <- sub$spot                          # 以注释 spot 为单元(顺序 == meta)
  # 仅取本样本存在的 LR 基因(符号); FetchData 朝向不稳(少基因= cells×genes,
  #   多基因= genes×cells), 用 barcode 成员关系归一为 基因×cells(列序 == bc)
  e <- FetchData(st, vars = intersect(genes.use, rownames(st)),
                 assay = "SCT", layer = "data")
  if (all(bc %in% rownames(e))) {
    Mraw <- t(e[bc, , drop = FALSE])          # cells×genes -> genes×cells
  } else if (all(bc %in% colnames(e))) {
    Mraw <- e[, bc, drop = FALSE]             # genes×cells(已正确)
  } else {
    stop("FetchData 输出无法对齐 cell barcode: ", s)
  }
  expr.list[[s]] <- Mraw                      # 基因 × 本样本注释 spots
  meta.list[[s]] <- data.frame(
    spot = sub$spot, sample = s,
    compartment = sub$compartment, labels = sub$labels,
    stringsAsFactors = FALSE)
  # 10x barcode 跨切片复用(同 barcode 出现在多样本), 用 样本|spot 作唯一键
  rownames(meta.list[[s]]) <- paste0(s, "|", sub$spot)
  cat(s, "表达矩阵:", nrow(Mraw), "基因 x", ncol(Mraw), "spot\n")
}
# 统一基因集后列向拼接, 保证 cbind 行名一致(官方 createCellChat 期望 genes×cells)
genes.common <- Reduce(intersect, lapply(expr.list, rownames))
expr.list <- lapply(expr.list, function(m) m[genes.common, , drop = FALSE])
# 每样本矩阵列名直接设为样本唯一键, cbind 后列名 == meta 行名(规避 barcode 跨切片复用)
for (i in seq_along(expr.list)) colnames(expr.list[[i]]) <- rownames(meta.list[[i]])
M <- do.call(cbind, expr.list)             # 基因 × 全部 spot(列名 = 样本唯一键)
meta <- do.call(rbind, meta.list)
# do.call(rbind, list) 会用列表元素名(样本)给行名加前缀, 变成
#   "ST-P1_colon.ST-P1_colon|..."(双前缀), 与 M 列名(单前缀)不一致;
#   createCellChat 会因此警告并自行用 colnames(M) 覆盖 rownames(meta)。
#   这里用 sample/spot 列重建与 M 列名完全一致的单前缀唯一键
rownames(meta) <- paste0(meta$sample, "|", meta$spot)
cat("DEBUG ncol(M)=", ncol(M), " nrow(meta)=", nrow(meta), "\n")
stopifnot(ncol(M) == nrow(meta))            # 列名已对齐, 无需再赋值
cat("合并表达矩阵:", nrow(M), "基因 x", ncol(M), "spot; 公共 LR 基因", length(genes.common), "\n")
cat("标签分布:\n"); print(table(meta$labels, useNA = "ifany"))
cat("PLAU/PLAUR 在合并矩阵中:", "PLAU" %in% rownames(M), "/", "PLAUR" %in% rownames(M), "\n")

# =============================================================================
# PART 2: 分区室运行 CellChat 并提取 PLAU-PLAUR 通路
# 目的: 量化各区室内 Tumor -> Myeloid 的 PLAU-PLAUR 通讯概率
# 方法: CellChat 标准流程(subsetData -> identifyOverExpressedInteractions
#       -> computeCommunProb -> computeCommunProbPathway -> aggregateNet)
# 数据: 输出 CellChatPLAUR.csv
# =============================================================================

cat("========== PART 2: Run CellChat by compartment ==========\n")

cc.list <- list()
plaur.tab <- list()

for (cp in COMP_RUN) {
  keep <- meta$compartment == cp
  M.cp <- M[, keep, drop = FALSE]
  meta.cp <- meta[keep, , drop = FALSE]

  # 剔除 spot 数过少的标签: 少于 MIN_CELLS 的组无法给出稳定通讯概率
  tb <- table(meta.cp$labels)
  keep.lab <- names(tb)[tb >= MIN_CELLS]
  M.cp <- M.cp[, meta.cp$labels %in% keep.lab, drop = FALSE]
  meta.cp <- meta.cp[meta.cp$labels %in% keep.lab, , drop = FALSE]
  cat(cp, "保留标签:", paste(keep.lab, collapse = ", "), "\n")

  cc <- createCellChat(M.cp, meta = meta.cp, group.by = "labels")
  cc@DB <- CellChatDB
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean", raw.use = TRUE)
  cc <- filterCommunication(cc, min.cells = MIN_CELLS)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc.list[[cp]] <- cc

  # 提取 PLAU-PLAUR 特异对的通讯概率(H3 直接证据)
  # CellChat 结构(源码确认, 1.6.1): cc@LR$LRsig 仅含显著 LR 对注释(11 列, 无 prob);
  #   逐对通讯概率在 cc@net$prob / cc@net$pval(3D 数组 source × target × nLR),
  #   第三维 dimnames = rownames(cc@LR$LRsig)(computeCommunProb 第 97/212 行),
  #   该键是 LRsig 行名、不一定等于 interaction_name 列。故按 LRsig 行号/行名定位
  #   第三维位置 k(曾因用 interaction_name 字符串匹配而取空 -> 根因已诊断修复)。
  lr <- cc@LR$LRsig
  idx <- which(lr$ligand == "PLAU" & lr$receptor == "PLAUR")
  if (length(idx) > 0) {
    k  <- idx[1]
    nm <- rownames(lr)[k]
    d3 <- dimnames(cc@net$prob)[[3]]
    if (nm %in% d3) {
      for (s in dimnames(cc@net$prob)[[1]]) {
        for (t in dimnames(cc@net$prob)[[2]]) {
          plaur.tab[[length(plaur.tab) + 1]] <- data.frame(
            compartment = cp,
            interaction = nm,
            ligand      = lr$ligand[k],
            receptor    = lr$receptor[k],
            source      = s,
            target      = t,
            prob        = as.numeric(cc@net$prob[s, t, k]),
            pval        = as.numeric(cc@net$pval[s, t, k]),
            stringsAsFactors = FALSE)
        }
      }
      cat(cp, "提取到 PLAU-PLAUR 通讯概率(第三维键='", nm, "', source->target 共",
          length(dimnames(cc@net$prob)[[1]]) * length(dimnames(cc@net$prob)[[2]]),
          "组)\n")
    } else {
      cat(cp, "PLAU-PLAUR 第三维键 '", nm, "' 不在 cc@net$prob(表达/过表达不足)\n")
    }
  } else {
    cat(cp, "PLAU-PLAUR 未进入显著 LR 对(LRsig)\n")
  }
  saveRDS(cc, file.path(TAB_DIR, paste0("CellChat_", cp, ".rds")))
}

if (length(plaur.tab) > 0) {
  plaur.df <- do.call(rbind, plaur.tab)
} else {
  plaur.df <- data.frame(compartment = character(), source = character(),
                         target = character(), prob = numeric(),
                         pval = numeric())
}
write.csv(plaur.df, file.path(TAB_DIR, "CellChatPLAUR.csv"),
          row.names = FALSE)
cat("PLAU-PLAUR 通讯结果:\n")
print(plaur.df)

# =============================================================================
# PART 3: Fig 4 可视化(PLAU-PLAUR 通路 + PT vs LM 比较)
# 目的: 展示 PLAU-PLAUR 轴在两区室的通讯强度与差异
# 方法: CellChat 官方 netVisual_bubble 与 rankNet
# 数据: 输出 CellChatPLAURPathway.pdf / CellChatComparePathway.pdf
# =============================================================================

cat("========== PART 3: Fig 4 plots ==========\n")

# Fig 4a: PLAU-PLAUR 特异对通讯气泡图(直接用 PART 2 已提取、已验证的 H3 数据)
# 说明(根因): 官方 netVisual_bubble(signaling="PLAU_PLAUR") 对本环境"手工补入的
#   自定义对"会调用 searchPair() 在 CellChatDB$interaction 中检索, 该检索对自定义
#   行匹配失败 -> 返回空 -> 报 "No significant signaling interactions..." 并中断脚本
#   (见 subsetCommunication_internal 第 26/108 行源码)。故绕过该调用, 直接用
#   plaur.df(prob/pval, 全部 source->target) 以 ggplot 绘制, 与 H3 数值完全一致、
#   稳健、不依赖 searchPair。气泡色=通讯概率, 大小=-log10(p), 灰点=不显著(prob=0)。
for (cp in COMP_RUN) {
  d <- plaur.df[plaur.df$compartment == cp, ]
  d$source.target <- factor(
    paste(d$source, "->", d$target),
    levels = paste(rep(c("Myeloid", "Tumor"), each = 2),
                   "->", rep(c("Myeloid", "Tumor"), 2)))
  d$prob.show <- d$prob
  d$prob.show[d$prob == 0] <- NA
  p1 <- ggplot(d, aes(x = source.target, y = interaction,
                      color = prob.show, size = -log10(pval + 1e-12))) +
    geom_point(pch = 16) + sci_theme +
    scale_color_gradient(low = "#fee0d2", high = "#cb181d",
                         na.value = "grey90", name = "Commun. Prob.") +
    scale_size_continuous(name = "-log10(p)", range = c(3, 9)) +
    labs(title = paste0("PLAU-PLAUR communication (", cp, ")"),
         x = "Source -> Target", y = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  pdf(file.path(FIG_DIR, paste0("CellChatPLAURPathway_", cp, ".pdf")),
      width = 7, height = 5, family = "Helvetica")
  print(p1)
  dev.off()
}

if (length(cc.list) == 2) {
  cc.merge <- mergeCellChat(cc.list, add.names = COMP_RUN)
  p2 <- rankNet(cc.merge, mode = "comparison", stacked = FALSE,
                do.stat = TRUE)
  pdf(file.path(FIG_DIR, "CellChatComparePathway.pdf"), width = 9,
      height = 7, family = "Helvetica")
  print(p2)
  dev.off()
}

# =============================================================================
# PART 3b: PLAU-PLAUR 通讯强度补充可视化(从 PART 2 已提取数据 plaur.df 绘制,
#       不依赖官方 signaling 专用 visualizer, 规避自定义对 searchPair 取空中断)
# 目的: 以更直观、更 SCI 标准的呈现通讯概率: 区室×方向热图 / PT vs LM 比较
# 数据: CellChatPLAUR_Heatmap.pdf / CellChatPLAUR_PTvLM.pdf
# =============================================================================

cat("========== PART 3b: PLAU-PLAUR supplementary plots ==========\n")

dir_levels <- paste(rep(c("Myeloid", "Tumor"), each = 2),
                    "->", rep(c("Myeloid", "Tumor"), 2))
plaur.df$source.target <- paste(plaur.df$source, "->", plaur.df$target)
plaur.df$source.target <- factor(plaur.df$source.target, levels = dir_levels)
plaur.df$sig <- plaur.df$pval < 0.05

# (1) 区室热图: 行=方向, 列=区室, 填充=通讯概率
p.hm <- ggplot(plaur.df, aes(x = compartment, y = source.target, fill = prob)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.4f", prob)), size = 4,
            colour = ifelse(plaur.df$prob > 0.004, "white", "black")) +
  scale_fill_gradient(low = "#f7f7f7", high = "#cb181d",
                      na.value = "grey90", name = "Prob.") +
  labs(title = "PLAU-PLAUR communication probability",
       x = "Compartment", y = "") +
  sci_theme + theme(axis.text.x = element_text(size = 10))
sci_save(p.hm, "CellChatPLAUR_Heatmap", width = 6, height = 5)

# (2) PT vs LM 比较(H3): 各方向 PT/LM 连线点图, 标注显著方向 LM/PT 比值
p.cmp <- ggplot(plaur.df, aes(x = source.target, y = prob,
                              colour = compartment, group = compartment)) +
  geom_line(colour = "grey50", linewidth = 0.5) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("PT" = COL_MYE, "LM" = COL_META),
                      name = "Compartment") +
  labs(title = "PLAU-PLAUR communication: PT vs LM (H3)",
       x = "", y = "Communication probability") +
  sci_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9))
mm <- plaur.df[plaur.df$source.target == "Myeloid -> Myeloid" & plaur.df$sig, ]
if (nrow(mm) == 2) {
  ratio <- mm$prob[mm$compartment == "LM"] / mm$prob[mm$compartment == "PT"]
  p.cmp <- p.cmp + annotate("text", x = "Myeloid -> Myeloid",
                            y = max(mm$prob) * 1.18,
                            label = sprintf("LM/PT = %.1f", ratio),
                            size = 3.5, colour = "black")
}
sci_save(p.cmp, "CellChatPLAUR_PTvLM", width = 9, height = 5)

# =============================================================================
# PART 4: nichenetr 校验 PLAU-PLAUR 与 SPP1 调控关系
# 目的: 用 NicheNet 先验独立确认 PLAU-PLAUR 是已知配体-受体对, 并检验
#       PLAU 是否为 SPP1 的高排名潜在调控配体
# 方法: 载入 lr_network 与 ligand_target_matrix; 查询 PLAU-PLAUR 配对;
#       对靶基因 SPP1 取调控潜力排名
# 数据: 输出 NicheNetLigandRank.csv
# =============================================================================

cat("========== PART 4: NicheNet prior check ==========\n")

# nichenetr 2.2.1.1 不随包分发 lr_network / ligand_target_matrix(见 NAMESPACE:
#   仅 annotation_data_sources / source_weights_df 等权重/注释类数据集); 官方
#   2.x 工作流从 Zenodo record 7074291 载入完整先验。
#   本环境已在 raw/ 提供预下载文件，优先从本地加载(避免外网依赖)。
load_nichenet_priors <- function() {
  # (0) 本地预下载文件(用户已提供)
  local_lr <- "G:/CRLM/ST/scCRLM/raw/lr_network_human_21122021.rds"
  local_ltm <- "G:/CRLM/ST/scCRLM/raw/ligand_target_matrix_nsga2r_final.rds"
  if (file.exists(local_lr) && file.exists(local_ltm)) {
    lr <- readRDS(local_lr)
    ltm <- readRDS(local_ltm)
    cat("NicheNet 先验来源: 本地文件\n")
    return(list(lr = lr, ltm = ltm, src = "local"))
  }
  # (1) 旧版随包数据集(向前兼容)
  ok <- tryCatch({
    data("lr_network", package = "nichenetr")
    data("ligand_target_matrix", package = "nichenetr")
    exists("lr_network") && exists("ligand_target_matrix")
  }, warning = function(w) FALSE, error = function(e) FALSE)
  if (ok) return(list(lr = lr_network, ltm = ligand_target_matrix, src = "package"))
  # (2) 官方 Zenodo 完整先验(2.x 标准来源); 设 30s 上限避免无外网时无限挂起
  base <- "https://zenodo.org/record/7074291/files"
  z <- tryCatch({
    setTimeLimit(elapsed = 30, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
    lr  <- readRDS(url(paste0(base, "/lr_network.rds"), "rb"))
    ltm <- readRDS(url(paste0(base, "/ligand_target_matrix.rds"), "rb"))
    list(lr = lr, ltm = ltm, src = "zenodo")
  }, warning = function(w) NULL, error = function(e) NULL)
  if (!is.null(z)) return(z)
  NULL
}
priors <- load_nichenet_priors()
if (is.null(priors)) {
  cat("NicheNet 先验不可用: nichenetr 2.2.1.1 不随包分发 lr_network/ligand_target_matrix,\n")
  cat("  且本环境无法从 Zenodo 下载(无外网)。PART 4 写入说明性 CSV, 不中断脚本。\n")
  write.csv(data.frame(note = paste0(
    "NicheNet priors unavailable in this environment: nichenetr 2.2.1.1 ships no ",
    "lr_network/ligand_target_matrix (see NAMESPACE); load from Zenodo record 7074291 ",
    "to enable PLAU-PLAUR pair check and PLAU->SPP1 ligand ranking.")),
    file.path(TAB_DIR, "NicheNetLigandRank.csv"), row.names = FALSE)
} else {
  lr_network <- priors$lr
  ligand_target_matrix <- priors$ltm
  cat("NicheNet 先验来源:", priors$src, "\n")

  lr.hit <- lr_network[lr_network$from == "PLAU" &
                       lr_network$to == "PLAUR", ]
  cat("PLAU-PLAUR 在 NicheNet lr_network 中的记录数:", nrow(lr.hit), "\n")
  if (nrow(lr.hit) > 0) {
    print(lr.hit)
  }

  # SPP1 的配体调控潜力排名: ligand_target_matrix 行为配体、列为靶基因
  spp1.pot <- ligand_target_matrix[, "SPP1"]
  rank.tab <- data.frame(ligand = names(spp1.pot),
                         potential = as.numeric(spp1.pot),
                         stringsAsFactors = FALSE)
  rank.tab <- rank.tab[order(-rank.tab$potential), ]
  rank.tab$rank <- seq_len(nrow(rank.tab))

  plau.rank <- rank.tab[rank.tab$ligand == "PLAU", ]
  cat("PLAU 对 SPP1 的调控潜力排名:", plau.rank$rank, "/", nrow(rank.tab),
      " (potential=", round(plau.rank$potential, 5), ")\n", sep = "")

  write.csv(head(rank.tab, 100),
            file.path(TAB_DIR, "NicheNetLigandRank.csv"), row.names = FALSE)
}
cat("11_PLAUR_SPP1_communication.R done:\n")
cat("  - CellChatPLAURPathway_PT.pdf / CellChatPLAURPathway_LM.pdf (bubble)\n")
cat("  - CellChatComparePathway.pdf (rankNet PT vs LM)\n")
cat("  - CellChatPLAUR_DirectionBar.pdf / CellChatPLAUR_Heatmap.pdf / CellChatPLAUR_PTvLM.pdf (supplement)\n")
