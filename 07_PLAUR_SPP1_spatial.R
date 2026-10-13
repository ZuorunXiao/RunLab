# ===============================================================
# 07_PLAUR_SPP1_spatial.R  
# 功能: 以空间热力图(离散 spot 瓦片)展示 4 个轴基因(SPP1/PLAUR/PLAU/FBXL12)的
#       **全切片 SCT 归一化表达**(每个 spot 都着色, 不限于肿瘤/巨噬子集),
#       颜色 = Min-Max 归一化 (0-1).
# 定位:回答 H2 空间邻近 / H3 轴强度 / H4 共定位
#       用 sct_<gene> 全切片表达(所有 spot 有值), 而非 子集限制签名
#       (SPP1_mac/PLAUR_mac/PLAU_tumor/FBXL12_tumor, 子集外 NA). 子集限制签名
#       仍作为 08/09/10/11/12/13 的统计单位; 本图仅作全切片表达可视化.
#       样本 = 本研究自有分析: 含 colon 原发病灶(P1/P2 colon)与肝转移(P1/P2 liver)共 4 样本; 仅参照原文可视化风格.
# 可视化: 每 (基因 x 样本) = 一张独立 PDF. 无 H&E 背景; 离散瓦片
#        (geom_tile, 按组织坐标铺排, 互不外推); 虚线凹包(concaveman)勾勒肿瘤区
#        (is_tumor=TRUE, 多灶按连通块分别描边); 加 "Tumor" 标签. 一图一pdf, 不合并.
# 数据: signature_master.rds (x,y,is_tumor, sct_* 基因表达列). 无需外部合并.
# 与 Wu et al. 原文 Fig 2C 区别 + 后续说明: 原文画 3 免疫 ssGSEA 亚群(复现锚点),
#   本图按轴画 4 基因全切片表达(回答 H2/H3/H4 空间背景); 后续为**本研究自有分析**,
#   仅参照原文可视化风格(黄-橙-红连续色阶、离散瓦片、虚线肿瘤轮廓), 不受原文"仅 LM"
#   限制 —— 故样本含 colon 原发病灶(P1/P2 colon)与肝转移(P1/P2 liver), 共 4 样本.
# ===============================================================

# --- 包加载 -------------------------------------------------------------
suppressMessages(library(ggplot2))   # 仅依赖 ggplot2
suppressMessages(library(concaveman))  # 凹包肿瘤边界(替代 chull 凸包, 更贴参考图)

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

sci_save <- function(plot_obj, filename, width = 7, height = 7) {
  pdf(file.path(FIG_DIR, paste0(filename, ".pdf")), width = width,
      height = height, family = "Helvetica")
  print(plot_obj)
  dev.off()
}

# 肿瘤 spot 凹轮廓(concaveman, 比 chull 凸包更贴合参考图肿瘤边界; 多灶时按连通块分别描边)
# spacing: 该样本 spot 最近邻间距(px), 用于连通块阈值(2.5*spacing 内视为同一灶)
# min_size: 连通块最小 spot 数, 低于此视为噪声斑点, 不描边(抑制 2-3 spot 假岛, 贴近参考图干净边界)
# 返回 list(polys, members, cluster_ids):
#   polys       : list of data.frame(x,y) 多边形顶点
#   members     : integer 向量(长度=肿瘤 spot 数), 每 spot 的连通块 id (噪声=0)
#   cluster_ids : 与 polys 平行, 每个多边形对应的连通块 id (用于取"最大灶"内真实 spot 质心)
# 无肿瘤 spot 返回 list(polys=list(), members=integer(0), cluster_ids=integer(0))
tumor_concave <- function(df, spacing, concavity = 2, min_size = 10) {
  tum <- df[df$is_tumor, c("x", "y")]
  if (nrow(tum) < 3)
    return(list(polys = list(), members = integer(0), cluster_ids = integer(0)))
  pts <- as.matrix(tum)
  # 连通块: 基于坐标距离 single-linkage 聚类(阈值 2.5*spacing), 避免把不相邻肿瘤灶强行连成一块
  d  <- as.dist(as.matrix(dist(pts)))
  cl <- cutree(hclust(d, method = "single"), h = 2.5 * spacing)
  polys <- list()
  cluster_ids <- integer(0)
  members <- cl                       # 每肿瘤 spot 的连通块 id, 噪声块置 0
  for (k in sort(unique(cl))) {
    sub_pts <- pts[cl == k, , drop = FALSE]
    if (nrow(sub_pts) < min_size) { members[cl == k] <- 0; next }
    poly <- tryCatch(
      concaveman::concaveman(sub_pts, concavity = concavity),
      error = function(e) sub_pts[chull(sub_pts[, "x"], sub_pts[, "y"]), ]
    )
    poly <- as.data.frame(poly)
    if (!all(c("x", "y") %in% colnames(poly))) colnames(poly) <- c("x", "y")
    polys[[length(polys) + 1]] <- poly
    cluster_ids <- c(cluster_ids, k)
  }
  list(polys = polys, members = members, cluster_ids = cluster_ids)
}

# --- 目录结构 -----------------------------------------------------------
BASE_DIR <- "X:/XX"
TAB06    <- file.path(BASE_DIR, "results/06_PLAUR_SPP1_signature/table")
TAB03    <- file.path(BASE_DIR, "results/03_immune_spatial/table")
FIG_DIR  <- file.path(BASE_DIR, "results/07_PLAUR_SPP1_spatial/Figure")
TAB_DIR  <- file.path(BASE_DIR, "results/07_PLAUR_SPP1_spatial/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 分析常量 -----------------------------------------------------------
# Fig 2C 四列 = 2 名未治疗患者 x (肝转移 LM, 结肠原发 CRC)
SAMPLE_ID <- c("ST-P1_colon", "ST-P1_liver", "ST-P2_colon", "ST-P2_liver")   # 本研究自有分析: colon 原发 + liver 转移 共 4 样本
# 亚群: 列名(rds) -> 显示标签 (英文, 规范17)
# 4 个轴基因的**全切片** SCT 表达(sct_* 列, 所有 spot 有值, 无子集限制) ->
#   每个 spot 都着色, 显示全部位表达量(用户指令: 不只肿瘤区域).
SUBSETS <- c(
  "sct_SPP1"    = "SPP1",
  "sct_PLAUR"   = "PLAUR",
  "sct_PLAU"    = "PLAU",
  "sct_FBXL12"  = "FBXL12"
)
SCALED_PREFIX <- "sc_"   # Min-Max 列前缀
# 黄-橙-红 连续色阶 (ColorBrewer YlOrRd 风格): 低=浅黄, 中=橙, 高=深红; 无蓝
HEAT_COLORS <- colorRampPalette(c("#FFFFCC", "#FED976", "#FEB24C",
                                  "#FD8D3C", "#FC4E2A", "#E31A1C",
                                  "#BD0026", "#800026"))(200)

# --- 读取 (4 轴基因 sct_* 表达列已在 signature_master.rds, 无需外部合并) -
cat("========== 07: Fig 2C 4-axis-gene full-section expression (plan §7.1) ==========\n")
master <- readRDS(file.path(TAB06, "signature_master.rds"))
cat("master rows:", nrow(master),
    " samples:", paste(unique(master$sample), collapse = ","), "\n")
cat("axis gene sct_* columns present:",
    all(c("sct_SPP1","sct_PLAUR","sct_PLAU","sct_FBXL12") %in% colnames(master)),
    "\n")

# --- Min-Max 归一化 (每基因全切片表达 0-1) ----------------------------
for (col in names(SUBSETS)) {
  v   <- master[[col]]
  rng <- range(v, na.rm = TRUE)
  master[[paste0(SCALED_PREFIX, col)]] <- (v - rng[1]) / (rng[2] - rng[1])
  cat(sprintf("  %-26s min=%.3f max=%.3f -> scaled 0-1\n", col, rng[1], rng[2]))
}

# 瓦片半宽: 以真实最近邻 spot 间距估计 (中位数 Euclidean 距离).
# 注意: 不能用 diff(unique(x))——ST 对齐后坐标经旋转/取整, 会出现 1 px 假象,
# 导致瓦片只有 1 px 宽、在图上不可见. 实际 Visium spot 间距 ≈ 200 px.
set.seed(20240909)
spot_spacing <- function(df) {
  n <- min(1000, nrow(df))
  idx <- sample.int(nrow(df), n)
  dm <- as.matrix(dist(df[idx, c("x", "y"), drop = FALSE]))
  diag(dm) <- Inf
  median(apply(dm, 1, min))
}
hw_map <- sapply(SAMPLE_ID, function(s) {
  df <- master[master$sample == s & !is.na(master$x), ]
  if (nrow(df) == 0) return(NA)
  spot_spacing(df)
})
cat("median nearest-neighbor spot spacing per sample (px):\n")
print(round(hw_map, 1))
HW_FACTOR <- 0.55   # 瓦片半宽 = 0.55 × 间距 → 瓦片宽度 ≈ 1.1 × 间距, 轻微重叠、填满点阵

# --- 清理旧版 PDF (Round H 授权彻底删除原可视化) ------------------------
old <- list.files(FIG_DIR, pattern = "\\.pdf$", full = TRUE)
if (length(old)) { file.remove(old); cat("removed", length(old), "old PDF(s)\n") }

# --- 主循环: 每 (亚群 x 样本) = 一张独立 PDF ---------------------------
n <- 0
for (s in SAMPLE_ID) {
  dfs <- master[master$sample == s & !is.na(master$x), ]
  if (nrow(dfs) == 0) { cat("WARN: no coords for", s, "\n"); next }
  sp <- if (is.na(hw_map[s])) median(hw_map, na.rm = TRUE) else hw_map[s]
  h <- tumor_concave(dfs, spacing = sp)
  s_lab <- gsub("_", " ", s)   # ST-P1 liver / ST-P1 colon / ...
  HW  <- HW_FACTOR * hw_map[s]

  for (col in names(SUBSETS)) {
    scol  <- paste0(SCALED_PREFIX, col)
    label <- SUBSETS[col]
    sub   <- dfs[!is.na(dfs[[scol]]), ]
    if (nrow(sub) == 0) { cat("WARN: no", col, "in", s, "\n"); next }

    p <- ggplot() +
      geom_tile(data = sub, aes(x = x, y = y, fill = .data[[scol]]),
                width = 2 * HW, height = 2 * HW, color = NA) +
      scale_fill_gradientn(colors = HEAT_COLORS, name = "Scaled proportion",
                           limits = c(0, 1), na.value = "grey90") +
      (if (length(h$polys) > 0) {
         hp <- do.call(rbind, lapply(seq_along(h$polys), function(i)
           cbind(h$polys[[i]], group = i)))
         big_i <- which.max(sapply(h$polys, nrow))   # 最大肿瘤灶
         big_cl <- h$cluster_ids[big_i]
         # 标签置于"最大灶内真实肿瘤 spot 质心"(而非凹多边形顶点均值,
         # 顶点均值对凹形会掉到区域外); 占位点恒在 spot 凸包内, 紧凑灶即在区域内
         cent <- dfs[dfs$is_tumor & h$members == big_cl, ]
         list(
           geom_polygon(data = hp, aes(x = x, y = y, group = group),
                        fill = NA, color = "#111111", linetype = "dashed",
                        linewidth = 1.1, inherit.aes = FALSE),
           annotate("text", x = mean(cent$x), y = mean(cent$y), label = "Tumor",
                    size = 4.5, fontface = "bold", color = "#111111",
                    hjust = 0.5, vjust = 0.5)
         )
       } else geom_blank()) +
      coord_fixed() + scale_y_reverse() + sci_theme +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 13, face = "plain")) +
      labs(title = paste0(label, "  |  ", s_lab))
    sci_save(p, paste0("HeatmapSubset_", col, "_", s))
    n <- n + 1
    cat("saved HeatmapSubset_", col, "_", s, ".pdf\n", sep = "")
  }
}

# --- 审计表: 每 spot 的 scaled proportion (可复核) ----------------------
out_cols <- c("sample", "spot", "x", "y", "is_tumor",
              paste0(SCALED_PREFIX, names(SUBSETS)))
write.csv(master[out_cols], file.path(TAB_DIR, "SubsetScaledProportion.csv"),
          row.names = FALSE)

cat("07 done:", n, "PDFs (4 genes x 4 samples: colon+liver, P1/P2), one plot per PDF.\n")
cat("palette: yellow-orange-red scaled-proportion; tiles discrete (HW=", round(HW, 2),
    "); dashed concaveman tumor outline; no H&E.\n", sep = "")
