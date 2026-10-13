# ===============================================================
# 10_PLAUR_SPP1_colocalization.R
# 功能: 检验 FBXL12-high 肿瘤 spot 与 PLAUR-high 巨噬 spot 的空间共定位,
#       并桥接 PLAUR+SPP1+ 共定位巨噬(效应端)与 FBXL12-high 肿瘤(分泌端)
# 定位: 方案 Fig 3 + H2(招募的空间相关)与 H4(PLAUR/SPP1 共定位)
# 可视化:
#   ColocalizationKcross.pdf         - H2/H4 跨类 L-cross 全局包络曲线
#   ColocalizationNearestNeighbour.pdf - H2/H4 平均最近邻跨类距离置换检验
#   ColocalizationSpatialMap.pdf     - A/B/C 高表达点集在组织坐标的空间叠加
#   MyeloidCoexpression.pdf          - H4a: 髓系内 PLAUR_mac~SPP1_mac Spearman rho
#   DoubleHighBoundary.pdf           - H4b: 双高髓系 spot 到 FBXL12-high 肿瘤距离
# 数据: 输入 results/06_PLAUR_SPP1_signature/table/signature_master.rds
#       输出 table/ColocalizationTests.csv / MyeloidCoexpression.csv
#            table/DoubleHighBoundary.csv
# 依据: 方案第 5 节(spatstat K-cross / L-cross + 全局 envelope + 置换)
# 方法:
#   (1) 距离标定: GetTissueCoordinates() 返回图像像素坐标。Visium 阵列为
#       六边形网格, 相邻 spot 中心距标称 100 um。以切片内最近邻距离的中位
#       数作为阵列 pitch 的像素长度并令其等于 100 um, 得到 um/px(data-
#       driven 标定, 不依赖 scale.factors 的具体定义, 已在 06 实测:
#       P1_liver 中位最近邻 = 200 px, 故 um/px = 0.5)
#   (2) 主检验: spatstat Kcross 的全局 envelope(nsim=999, 随机标签零假设
#       rlabel), 以 L-cross = sqrt(K/pi) - r 呈现; 观测曲线高于包络上界
#       表示两类 spot 在该尺度上显著共定位(聚集), 低于下界表示分离
#   (3) 补充检验: 平均最近邻跨类距离的置换检验(透明、不依赖 spatstat 细节)
# ===============================================================

# --- 包加载(统计 -> 数据 -> 可视化 -> 领域) -----------------------------
suppressMessages({
  library(spatstat.geom)    # ppp 点模式对象与 rlabel
  library(spatstat.explore) # Kcross 与 envelope 全局包络
  library(ggplot2)          # 可视化
  library(reshape2)         # 长宽格式转换
  library(matrixStats)      # rowMins 计算最近邻距离(显式依赖, 规范22)
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

# 显著性星号: p 值 -> */**/***/ns
SigStar <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

# --- 目录结构 -----------------------------------------------------------
BASE_DIR <- "X:/XX"
TAB06    <- file.path(BASE_DIR, "results/06_PLAUR_SPP1_signature/table")
FIG_DIR  <- file.path(BASE_DIR, "results/10_PLAUR_SPP1_colocalization/Figure")
TAB_DIR  <- file.path(BASE_DIR, "results/10_PLAUR_SPP1_colocalization/table")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 分析常量 -----------------------------------------------------------
SAMPLE_ID  <- c("ST-P1_liver", "ST-P1_colon", "ST-P2_liver", "ST-P2_colon")
HIGH_Q     <- 0.75   # 阳性阈值: 子集内签名表达的上四分位(方案三档之一)
PITCH_UM   <- 100    # Visium 相邻 spot 中心距标称值(um)
N_SIM      <- 999    # 全局 envelope 模拟次数(方案第5节: 全局 envelope >=1000, 取 999)
N_PERM     <- 999    # 最近邻距离置换检验次数
SEED       <- 42

# 项目配色(与 09 一致): A=FBXL12-high 肿瘤(分泌端, 红);
# B=PLAUR-high 巨噬(蓝); C=SPP1-high 巨噬(橙)
COL_A <- "#E31A1C"
COL_B <- "#1F78B4"
COL_C <- "#FD8D3C"
TEST_COL <- c("H2_FBXL12_PLAUR" = COL_B, "H4_PLAUR_SPP1" = COL_C)

# =============================================================================
# PART 1: 定义各切片的高表达点集并做距离标定
# 目的: 得到每切片的 A(FBXL12-high 肿瘤)与 B(PLAUR-high 巨噬)点坐标
# 方法: 在子集(肿瘤/巨噬)内按 HIGH_Q 取上四分位; 最近邻中位数标定 um/px
# 数据: 中间产物, 供 PART 2-3 / PART 4b 使用
# =============================================================================

cat("========== PART 1: Define high point sets ==========\n")

set.seed(SEED)
master <- readRDS(file.path(TAB06, "signature_master.rds"))

# UmPerPixel: 以一个坐标数据框的最近邻距离中位数标定物理尺度
UmPerPixel <- function(xy) {
  d <- as.matrix(dist(xy))
  diag(d) <- Inf
  PITCH_UM / median(apply(d, 1, min))
}

pat.list <- list()
for (s in SAMPLE_ID) {
  d <- master[master$sample == s, ]
  d <- d[!is.na(d$x) & !is.na(d$y), ]

  # A: 肿瘤 spot 中 FBXL12 表达上四分位
  dt <- d[d$is_tumor & !is.na(d$sct_FBXL12), ]
  # B: 巨噬 spot 中 PLAUR 表达上四分位
  dm <- d[d$is_macro & !is.na(d$sct_PLAUR), ]
  # C: 巨噬 spot 中 SPP1 表达上四分位(用于 H4 的 PLAUR/SPP1 共定位)
  ds <- d[d$is_macro & !is.na(d$sct_SPP1), ]

  qa <- quantile(dt$sct_FBXL12, HIGH_Q, names = FALSE)
  qb <- quantile(dm$sct_PLAUR, HIGH_Q, names = FALSE)
  qc <- quantile(ds$sct_SPP1, HIGH_Q, names = FALSE)

  pa <- dt[dt$sct_FBXL12 >= qa, c("x", "y")]
  pb <- dm[dm$sct_PLAUR >= qb, c("x", "y")]
  pc <- ds[ds$sct_SPP1 >= qc, c("x", "y")]

  um.px <- UmPerPixel(d[, c("x", "y")])
  pat.list[[s]] <- list(all = d[, c("x", "y")], A = pa, B = pb, C = pc,
                        um_px = um.px, nA = nrow(pa), nB = nrow(pb),
                        nC = nrow(pc))
  cat(sprintf("%s: nA=%d nB=%d nC=%d um/px=%.4f\n", s, nrow(pa),
              nrow(pb), nrow(pc), um.px))
}

# =============================================================================
# PART 2: Kcross 全局 envelope 检验(主检验)
# 目的: 检验 A 类与 B 类点是否存在显著空间聚集(共定位)
# 方法: 合并 A/B 为带标记 ppp, 以 rlabel(随机标签)为零假设做 envelope,
#       统计量取 L-cross = sqrt(Kcross/pi) - r
# 数据: 输出 ColocalizationTests.csv 的 kcross 部分
# =============================================================================

cat("========== PART 2: Kcross global envelope ==========\n")

# MakePpp: 由两类点坐标与总观测窗构造多类型点模式
MakePpp <- function(all.xy, p1, p2, name1, name2) {
  win <- owin(xrange = range(all.xy$x), yrange = range(all.xy$y))
  xy <- rbind(p1, p2)
  mk <- factor(c(rep(name1, nrow(p1)), rep(name2, nrow(p2))),
               levels = c(name1, name2))
  ppp(x = xy$x, y = xy$y, window = win, marks = mk)
}

env.list <- list()
test.rows <- list()

for (s in SAMPLE_ID) {
  pt <- pat.list[[s]]
  # 检验 1: H2 - FBXL12-high 肿瘤 vs PLAUR-high 巨噬
  # 检验 2: H4 - PLAUR-high 巨噬 vs SPP1-high 巨噬
  spec <- list(
    list(id = "H2_FBXL12_PLAUR", p1 = pt$A, p2 = pt$B,
         n1 = "A", n2 = "B"),
    list(id = "H4_PLAUR_SPP1", p1 = pt$B, p2 = pt$C,
         n1 = "B", n2 = "C")
  )

  for (sp in spec) {
    X <- MakePpp(pt$all, sp$p1, sp$p2, sp$n1, sp$n2)
    # 理论 K = pi r^2; 取 L = sqrt(K/pi) - r 使零假设下期望为 0
    env <- envelope(X, Kcross, nsim = N_SIM, i = sp$n1, j = sp$n2,
                    simulate = expression(rlabel(X)),
                    correction = "Ripley",
                    transform = expression(sqrt(./pi) - r),
                    savefuns = FALSE, verbose = FALSE)
    df <- as.data.frame(env)
    df$sample <- s
    df$test <- sp$id
    df$r_um <- df$r * pt$um_px
    # L 值同样按 um/px 缩放, 使 y 轴与 x 轴(r_um)单位一致(均为 um)
    df$obs <- df$obs * pt$um_px
    df$lo  <- df$lo  * pt$um_px
    df$hi  <- df$hi  * pt$um_px
    env.list[[paste(s, sp$id)]] <- df

    # 汇总: 观测 L-cross 在最大尺度处的取值与其是否越出包络
    k <- nrow(df)
    obs <- df$obs[k]
    lo <- df$lo[k]
    hi <- df$hi[k]
    test.rows[[length(test.rows) + 1]] <- data.frame(
      sample = s, test = sp$id, method = "Kcross envelope",
      r_um_max = round(df$r_um[k], 1),
      L_obs = round(obs, 4), L_lo = round(lo, 4), L_hi = round(hi, 4),
      above_upper = obs > hi, below_lower = obs < lo,
      stringsAsFactors = FALSE
    )
    cat(sprintf("%s %s: L_obs=%.3f [%.3f, %.3f] above=%s\n", s, sp$id,
                obs, lo, hi, obs > hi))
  }
}

env.tab <- do.call(rbind, env.list)

# =============================================================================
# PART 3: 平均最近邻跨类距离的置换检验(补充检验)
# 目的: 用一个不依赖 spatstat 细节的透明指标复核共定位
# 方法: 对 A 类每点求到 B 类的最近距离取均值; 随机置换 A/B 标签 999 次
#       得到零分布, 计算双侧 p 值
# 数据: 输出 ColocalizationTests.csv 的 nn 部分 + nn.tab(供 PART 3 图)
# =============================================================================

cat("========== PART 3: Nearest neighbour permutation test ==========\n")

nn.rows <- list()
for (s in SAMPLE_ID) {
  pt <- pat.list[[s]]
  spec <- list(
    list(id = "H2_FBXL12_PLAUR", p1 = pt$A, p2 = pt$B),
    list(id = "H4_PLAUR_SPP1", p1 = pt$B, p2 = pt$C)
  )
  for (sp in spec) {
    pool <- rbind(sp$p1, sp$p2)
    n1 <- nrow(sp$p1); n2 <- nrow(pool) - n1
    if (n1 < 2 || n2 < 2) {
      obs <- NA_real_; null <- NA_real_; p.lo <- NA_real_; p.hi <- NA_real_
    } else {
      # 预先计算完整成对距离矩阵一次, 置换仅做索引(避免每次重算 dist 致超时)
      D <- as.matrix(dist(pool))
      obs <- mean(rowMins(D[1:n1, (n1 + 1):nrow(D), drop = FALSE]))
      null <- replicate(N_PERM, {
        idx <- sample(nrow(pool))
        a <- idx[1:n1]; b <- idx[(n1 + 1):nrow(pool)]
        mean(rowMins(D[a, b, drop = FALSE]))
      })
      p.lo <- (sum(null <= obs, na.rm = TRUE) + 1) / (N_PERM + 1)
      p.hi <- (sum(null >= obs, na.rm = TRUE) + 1) / (N_PERM + 1)
    }
    test.rows[[length(test.rows) + 1]] <- data.frame(
      sample = s, test = sp$id, method = "NN permutation",
      r_um_max = NA, L_obs = round(obs * pt$um_px, 2),
      L_lo = round(mean(null, na.rm = TRUE) * pt$um_px, 2), L_hi = NA,
      above_upper = NA, below_lower = NA,
      stringsAsFactors = FALSE
    )
    nn.rows[[length(nn.rows) + 1]] <- data.frame(
      sample = s, test = sp$id,
      obs_um = if (is.na(obs)) NA_real_ else round(obs * pt$um_px, 2),
      null_mean_um = if (is.na(obs)) NA_real_
        else round(mean(null, na.rm = TRUE) * pt$um_px, 2),
      null_sd_um = if (is.na(obs)) NA_real_
        else round(sd(null) * pt$um_px, 2),
      p_lo = p.lo, p_hi = p.hi,
      stringsAsFactors = FALSE
    )
    cat(sprintf("%s %s: NN obs=%.1f um, null mean=%.1f um, p(closer)=%.4f\n",
                s, sp$id, obs * pt$um_px, mean(null, na.rm = TRUE) * pt$um_px,
                p.lo))
  }
}

nn.tab <- do.call(rbind, nn.rows)
nn.tab$p_two <- ifelse(is.na(nn.tab$p_lo), NA_real_,
                       pmin(1, 2 * pmin(nn.tab$p_lo, nn.tab$p_hi)))
nn.tab$diff_um <- nn.tab$obs_um - nn.tab$null_mean_um
nn.tab$z_score <- nn.tab$diff_um / nn.tab$null_sd_um

test.tab <- do.call(rbind, test.rows)
write.csv(test.tab, file.path(TAB_DIR, "ColocalizationTests.csv"),
          row.names = FALSE)
write.csv(nn.tab, file.path(TAB_DIR, "NearestNeighbourPermutation.csv"),
          row.names = FALSE)

# =============================================================================
# PART 3b: 最近邻置换检验结果可视化
# 目的: 直观呈现观测平均跨类最近邻距离偏离随机标签零分布的程度与方向
#       负值 = 观测距离短于随机期望 = 空间共定位; 正值 = 空间分离
# 方法: 每检验分面; 柱/点 = observed - null mean distance (um); 0 线为随机期望;
#       显著性星号标注双侧 p
# 数据: ColocalizationNearestNeighbour.pdf
# =============================================================================

cat("========== PART 3b: NN permutation plot ==========\n")

lab.map <- c("H2_FBXL12_PLAUR" = "FBXL12 tumor vs PLAUR macrophage",
             "H4_PLAUR_SPP1" = "PLAUR macrophage vs SPP1 macrophage")

nn.tab$sample <- factor(nn.tab$sample, levels = SAMPLE_ID)
nn.tab$test <- factor(nn.tab$test,
                      levels = c("H2_FBXL12_PLAUR", "H4_PLAUR_SPP1"))
nn.tab$test_label <- lab.map[as.character(nn.tab$test)]

p.nn <- ggplot(nn.tab, aes(x = sample, y = diff_um, fill = diff_um < 0)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = SigStar(p_two), y = pmax(diff_um, 0)),
            vjust = -1.0, size = 4, colour = "black") +
  scale_fill_manual(values = c("TRUE" = COL_C, "FALSE" = COL_B),
                    labels = c("closer than random", "farther than random"),
                    name = "Direction") +
  facet_wrap(~ test_label, scales = "free_y") +
  labs(title = "Deviation from random labelling in cross-type nearest-neighbour distance",
       x = "Sample", y = "observed - null mean distance (um)") +
  sci_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
                    strip.text = element_text(size = 9),
                    legend.position = "bottom")
sci_save(p.nn, "ColocalizationNearestNeighbour", width = 9, height = 6)

# =============================================================================
# PART 4: Fig 3 共定位曲线图(L-cross 观测值 + 包络带)
# 目的: 可视化 H2/H4 在各切片的共定位曲线, 观测线按检验上色
# 方法: ggplot2 分面, 阴影为全局包络, 0 线为随机标签零假设
# 数据: ColocalizationKcross.pdf
# =============================================================================

cat("========== PART 4: Fig 3 colocalization plot ==========\n")

env.tab$test <- factor(env.tab$test,
                       levels = c("H2_FBXL12_PLAUR", "H4_PLAUR_SPP1"))
env.tab$test_label <- lab.map[as.character(env.tab$test)]

p.kc <- ggplot(env.tab, aes(x = r_um)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = test), alpha = 0.25,
              colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4,
             colour = "grey40") +
  geom_line(aes(y = obs, colour = test), linewidth = 0.9) +
  scale_colour_manual(values = TEST_COL, name = "Comparison",
                      labels = lab.map) +
  scale_fill_manual(values = TEST_COL, name = "Comparison",
                    labels = lab.map) +
  facet_wrap(~ sample, ncol = 4, scales = "free_y") +
  labs(title = "Cross-type L function with global envelope",
       x = "Distance (um)", y = "L(r) - r (deviation from random labelling)") +
  sci_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 8))
sci_save(p.kc, "ColocalizationKcross", width = 12, height = 4.5)

# =============================================================================
# PART 4b: 空间共定位散点图(A/B/C 高表达点集在组织坐标叠加)
# 目的: 最直观地呈现共定位——高表达点集在组织切片中的空间邻近关系
# 方法: 背景=全部 spot(灰); A=FBXL12-high 肿瘤(红), B=PLAUR-high 巨噬(蓝),
#       C=SPP1-high 巨噬(橙); 每样本分面, 等比坐标
# 数据: ColocalizationSpatialMap.pdf
# =============================================================================

cat("========== PART 4b: Spatial colocalization map ==========\n")

map.rows <- list()
for (s in SAMPLE_ID) {
  pt <- pat.list[[s]]
  bg <- data.frame(x = pt$all$x, y = pt$all$y, cls = "All spots")
  a  <- data.frame(x = pt$A$x, y = pt$A$y, cls = "A FBXL12 tumor")
  b  <- data.frame(x = pt$B$x, y = pt$B$y, cls = "B PLAUR macrophage")
  c  <- data.frame(x = pt$C$x, y = pt$C$y, cls = "C SPP1 macrophage")
  m <- rbind(bg, a, b, c)
  m$sample <- s
  map.rows[[length(map.rows) + 1]] <- m
}
map.df <- do.call(rbind, map.rows)
map.df$cls <- factor(map.df$cls,
  levels = c("All spots", "A FBXL12 tumor", "B PLAUR macrophage",
             "C SPP1 macrophage"))
map.df$sample <- factor(map.df$sample, levels = SAMPLE_ID)
cls_col <- c("All spots" = "#BFBFBF", "A FBXL12 tumor" = COL_A,
             "B PLAUR macrophage" = COL_B, "C SPP1 macrophage" = COL_C)

p.map <- ggplot(map.df, aes(x = x, y = y, colour = cls)) +
  geom_point(size = 0.8, alpha = 0.55) +
  scale_colour_manual(values = cls_col, name = "Spot set") +
  facet_wrap(~ sample, ncol = 2) +
  coord_fixed(ratio = 1) +
  labs(title = "Spatial colocalization of high-expression spot sets",
       x = "Tissue x (pixel)", y = "Tissue y (pixel)", colour = "Spot set") +
  sci_theme + theme(strip.text = element_text(size = 9),
                    legend.text = element_text(size = 9))
sci_save(p.map, "ColocalizationSpatialMap", width = 10, height = 8)

# =============================================================================
# PART 5: H4a 髓系内 PLAUR_mac × SPP1_mac 共表达相关(细胞级共定位证据)
# 目的: 直接检验"PLAUR+SPP1+ 共定位巨噬"——同一髓系 spot 内两签名是否正相关
#       这是对 H4 的核心检验, 区别于 PART2-3 的"两高表达点集空间共现"
# 方法: 每样本 + 合并, 对 is_macro spot 的 PLAUR_mac ~ SPP1_mac 做 Spearman 相关
# 数据: MyeloidCoexpression.csv
# =============================================================================

cat("========== PART 5: H4a myeloid PLAUR~SPP1 co-expression ==========\n")

coex.rows <- list()
for (s in SAMPLE_ID) {
  dm <- master[master$sample == s & master$is_macro, ]
  v1 <- dm$PLAUR_mac; v2 <- dm$SPP1_mac
  keep <- !is.na(v1) & !is.na(v2)
  if (sum(keep) >= 3) {
    r <- cor.test(v1[keep], v2[keep], method = "spearman")
    coex.rows[[length(coex.rows) + 1]] <- data.frame(
      sample = s, n_myeloid = sum(keep),
      rho = round(r$estimate, 4), p = signif(r$p.value, 3),
      stringsAsFactors = FALSE)
    cat(sprintf("%s: rho=%.3f p=%.2e (n=%d)\n", s, r$estimate,
                r$p.value, sum(keep)))
  }
}
dm <- master[master$is_macro, ]
v1 <- dm$PLAUR_mac; v2 <- dm$SPP1_mac; keep <- !is.na(v1) & !is.na(v2)
r <- cor.test(v1[keep], v2[keep], method = "spearman")
coex.rows[[length(coex.rows) + 1]] <- data.frame(
  sample = "ALL", n_myeloid = sum(keep),
  rho = round(r$estimate, 4), p = signif(r$p.value, 3),
  stringsAsFactors = FALSE)
coex.tab <- do.call(rbind, coex.rows)
write.csv(coex.tab, file.path(TAB_DIR, "MyeloidCoexpression.csv"),
          row.names = FALSE)
cat(sprintf("ALL: rho=%.3f p=%.2e (n=%d)\n", r$estimate, r$p.value,
            sum(keep)))

# =============================================================================
# PART 6: H4b 双高(PLAUR+SPP1+)髓系 spot 与 FBXL12-high 肿瘤的空间邻近
# 目的: PLAUR+SPP1+ 共定位巨噬是否比其它髓系细胞更靠近 FBXL12-high 肿瘤(分泌端)
#       此检验将 H2(分泌端 FBXL12-high 肿瘤)与 H4(效应端双高巨噬)直接桥接
# 方法: 每样本 is_macro 内双高(两签名>=样本内 HIGH_Q)分类;
#       计算各髓系 spot 到 PART1 A(FBXL12-high 肿瘤)最近距离(px->um);
#       Mann-Whitney 比较双高 vs 其余髓系 spot 到肿瘤的最近距离
# 数据: DoubleHighBoundary.csv
# =============================================================================

cat("========== PART 6: H4b double-high myeloid proximity to FBXL12 tumor ==========\n")

bh.rows <- list()
for (s in SAMPLE_ID) {
  dm <- master[master$sample == s & master$is_macro, ]
  A  <- pat.list[[s]]$A
  if (nrow(dm) < 5 || nrow(A) < 1) {
    cat(sprintf("%s: skip (n_myeloid=%d nA=%d)\n", s, nrow(dm), nrow(A)))
    next
  }
  q1 <- quantile(dm$PLAUR_mac, HIGH_Q, na.rm = TRUE)
  q2 <- quantile(dm$SPP1_mac, HIGH_Q, na.rm = TRUE)
  dm$dh <- !is.na(dm$PLAUR_mac) & !is.na(dm$SPP1_mac) &
           dm$PLAUR_mac >= q1 & dm$SPP1_mac >= q2
  D <- as.matrix(dist(rbind(dm[, c("x", "y")], A)))
  m <- nrow(dm)
  dA_um <- apply(D[1:m, (m + 1):nrow(D), drop = FALSE], 1, min) *
           pat.list[[s]]$um_px
  # 双高组与其余组均需 >=2 个观测方能做 Mann-Whitney; 否则标记 NA
  # (出现于某签名在髓系内饱和、几乎全部 spot 落入双高阈值的情况)
  if (sum(dm$dh) >= 2 && sum(!dm$dh) >= 2) {
    w <- wilcox.test(dA_um[dm$dh], dA_um[!dm$dh])
    wp <- signif(w$p.value, 3)
  } else {
    wp <- NA_real_
  }
  bh.rows[[length(bh.rows) + 1]] <- data.frame(
    sample = s, n_myeloid = nrow(dm), n_double_high = sum(dm$dh),
    med_dist_double_high_um = if (sum(dm$dh) >= 1) round(median(dA_um[dm$dh]), 1) else NA_real_,
    med_dist_other_um = if (sum(!dm$dh) >= 1) round(median(dA_um[!dm$dh]), 1) else NA_real_,
    W_p = wp, stringsAsFactors = FALSE)
  cat(sprintf("%s: n_dh=%d medDist dh=%.0f vs other=%.0f p=%s\n", s,
              sum(dm$dh),
              if (sum(dm$dh) >= 1) round(median(dA_um[dm$dh]), 0) else NA,
              if (sum(!dm$dh) >= 1) round(median(dA_um[!dm$dh]), 0) else NA,
              if (is.na(wp)) "NA(degenerate)" else sprintf("%.2e", wp)))
}
bh.tab <- do.call(rbind, bh.rows)
write.csv(bh.tab, file.path(TAB_DIR, "DoubleHighBoundary.csv"),
          row.names = FALSE)

# =============================================================================
# PART 7: H4a/H4b 可视化(拆分为独立 PDF, 颜色统一项目配色)
# 目的: H4a 共表达强度(rho)与 H4b 双高巨噬对肿瘤的邻近程度各自成图
# 数据: MyeloidCoexpression.pdf + DoubleHighBoundary.pdf
# =============================================================================

cat("========== PART 7: H4a/H4b plots ==========\n")

coex.tab$star <- SigStar(coex.tab$p)

p_coex <- ggplot(coex.tab, aes(x = sample, y = rho, fill = rho > 0)) +
  geom_col() + geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sprintf("%.2f", rho)), vjust = -0.6, size = 3,
            colour = "black") +
  geom_text(aes(label = star, y = ifelse(rho > 0, rho, 0)),
            vjust = -2.0, size = 4, colour = "black") +
  scale_fill_manual(values = c("TRUE" = COL_A, "FALSE" = COL_B)) +
  labs(title = "PLAUR_mac ~ SPP1_mac Spearman rho in myeloid spots",
       y = "Spearman rho", x = "") +
  sci_theme + theme(legend.position = "none",
                    axis.text.x = element_text(angle = 20, hjust = 1, size = 9))
sci_save(p_coex, "MyeloidCoexpression", width = 8, height = 6)

# 效应量: other myeloid 与 double-high 的中位距离差(正 = double-high 更近)
bh.tab$diff_um <- bh.tab$med_dist_other_um - bh.tab$med_dist_double_high_um
bh.tab$sig <- SigStar(bh.tab$W_p)
bh.tab$sample <- factor(bh.tab$sample, levels = SAMPLE_ID)

p_bh <- ggplot(bh.tab, aes(x = sample, y = diff_um, fill = diff_um > 0)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = sprintf("%.0f", diff_um)), vjust = -0.6, size = 3,
            colour = "black") +
  geom_text(aes(label = sig, y = pmax(diff_um, 0)), vjust = -2.0, size = 4,
            colour = "black") +
  scale_fill_manual(values = c("TRUE" = COL_C, "FALSE" = COL_B)) +
  labs(title = "Proximity advantage of double-high myeloid to FBXL12-high tumor",
       y = "median distance difference (um)", x = "") +
  sci_theme + theme(legend.position = "none",
                    axis.text.x = element_text(angle = 20, hjust = 1, size = 9))
sci_save(p_bh, "DoubleHighBoundary", width = 8, height = 6)

cat("10_PLAUR_SPP1_colocalization.R done:\n")
cat("  - ColocalizationKcross.pdf\n")
cat("  - ColocalizationNearestNeighbour.pdf\n")
cat("  - ColocalizationSpatialMap.pdf\n")
cat("  - MyeloidCoexpression.pdf\n")
cat("  - DoubleHighBoundary.pdf\n")
