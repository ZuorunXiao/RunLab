# =====================================================================
# RARP-21_workflow_figure.R
# 功能：绘制全流程示意图（四个证据层 x 七个判定模块 x 定向综合出口），
#       按期刊版心尺寸输出可直接付印的矢量 PDF。
# 定位：图 1 的唯一图源；图内系列与样本量取自补充表 S1，网络拟合与模块
#       计数取自 S08 的输出表，判定模块的编号与正文 2.1 节的枚举逐条对应。
# 可视化：输出尺寸即最终印刷尺寸（宽 6.0 英寸，与版心同宽），正文不小于
#       6 pt；四色语义（米色=数据集、深蓝=判定模块、浅蓝=支撑分析、绿色=
#       综合出口），判定模块在框顶左角带编号圆标；无网格、无阴影。
# =====================================================================

library(grid)  # 圆角矩形、任意坐标文本与折线箭头；示意图不需要绘图系统

# --- 出口与画布 ---
FIG_DIR <- "X:/RAPA/04_manuscript"  # 图件落盘目录
FIG_FILE <- "RAPA_analysis_flowchart.pdf"  # 图件文件名
PAGE_W_IN <- 6  # 页面宽（英寸），等于版心宽，输出即印刷尺寸
CANVAS_W <- 900  # 画布宽（用户单位）
UNIT_PT <- PAGE_W_IN * 72 / CANVAS_W  # 一用户单位对应的磅值（0.48）

# --- 版式常量 ---
MARGIN_X <- 20  # 左边界
RULE_X2 <- 846  # 层分隔线与内容区右端，留出右侧汇流通道
DAT_X <- 20  # 数据集列的左边界
DAT_W <- 164  # 数据集框宽
DAT_H <- 46  # 数据集框高
DAT_STEP <- 52  # 数据集框行距（框高 + 间隔）
BUS_X <- 194  # 进线汇流竖线，位于数据集列与分析列之间
ANAL_X0 <- 220  # 首个分析框左边界
ANAL_GAP <- 24  # 同行相邻分析框的间隔
BOX_MAX_W <- (RULE_X2 - ANAL_X0 - 2 * ANAL_GAP) / 3  # 单框宽上限，即三框行
BOX_H <- 100  # 分析框高
ROW_GAP <- 40  # 同层两行分析框的间隔，留出编号圆标与出线横线的余量
OUT_X <- 220  # 综合出口框左边界，与分析列左边界对齐
OUT_W <- 626  # 综合出口框宽，右边界与内容区右端对齐
OUT_H <- 108  # 综合出口框高
OUT_GAP <- 70  # 末层底部到综合出口框的距离
LANE_GAP <- 46  # 本层底部到下一层层标题的距离
FIRST_TOP <- 86  # 首层内容顶沿，留出图例行
RULE_UP <- 10  # 层顶沿到层分隔线
HEAD_UP <- 20  # 层顶沿到层标题行中心
BOTTOM_MARGIN <- 26  # 画布底部留白
BOX_R <- 5  # 圆角半径（用户单位）
OUT_R <- 7  # 综合出口框圆角半径
COL_X <- 872  # 右侧汇流竖线的横坐标，距画布右缘 28 单位
BUS_OUT <- 14  # 分析框底到出线横线的垂直距离
ARROW_LEN <- 10  # 箭头头部长度（用户单位）
ARROW_HALF <- 6.8  # 箭头头部半宽，总宽约 13.6 单位
LINE_W <- 3  # 连线粗细（用户单位）
FRAME_W <- 2  # 数据集框描边（用户单位）
ANAL_W <- 2.2  # 分析框描边（用户单位）
OUT_W_LINE <- 2.6  # 综合出口框描边（用户单位）
CORNER_N <- 12  # 每个圆角的插值点数
TEXT_PAD <- 10  # 文字距框左边界

# --- 行内文本位置（自框顶算起，取行中心） ---
TITLE_DY <- 19  # 分析框标题行中心
LINE_DY <- c(42, 62, 82)  # 分析框三行正文的行中心
DAT_TITLE_DY <- 15  # 数据集框标题行中心
DAT_SUB_DY <- 31  # 数据集框副标题行中心
OUT_TITLE_DY <- 23  # 综合出口标题行中心

# --- 编号圆标 ---
BADGE_R <- 8.5  # 圆标半径（用户单位），圆心落在框顶边上
BADGE_DX <- 16  # 圆心距框左边界的水平距离
BADGE_SIZE <- 11  # 圆标数字字号（5.28 pt）

# --- 字号（用户单位；1 单位 = 0.48 pt，故 12.5 单位即 6.0 pt） ---
SIZE_LANE <- 13  # 层标题（6.24 pt）
SIZE_TITLE <- 13.5  # 分析框标题（6.48 pt）
SIZE_BODY <- 12.5  # 框内正文（6.0 pt）
SIZE_OUT_TITLE <- 15  # 综合出口标题（7.2 pt）
SIZE_LEGEND <- 12.5  # 图例（6.0 pt）

# --- 配色 ---
COL_TITLE <- "#1a1a1a"  # 标题与层名
COL_BODY <- "#333333"  # 正文
COL_SUB <- "#4a4a4a"  # 副标题
COL_LINE <- "#2e2e2e"  # 箭头与连线
COL_RULE <- "#cccccc"  # 层分隔线
FILL_DAT <- "#f5e9d3"  # 数据集填充
EDGE_DAT <- "#a8843f"  # 数据集描边
FILL_ANAL <- "#d5e0ea"  # 判定模块填充
EDGE_ANAL <- "#3f6484"  # 判定模块描边
FILL_SUP <- "#eef2f7"  # 支撑分析填充
EDGE_SUP <- "#8fa5b8"  # 支撑分析描边
FILL_OUT <- "#d6efe6"  # 综合出口填充
EDGE_OUT <- "#2e8b74"  # 综合出口描边
COL_WHITE <- "#ffffff"  # 画布底色

# =====================================================================
# PART 1 版面几何
# 目的：由各层数据集框数与各行分析框数递推全部坐标与框宽，不手工累加。
# 方法：层高取数据集列与分析框块的高者，两列各自纵向居中；分析框宽度按
#       该层同行最密的行等分，并设上限；画布高由末层推出。
# 数据：数据集框数 5/4/3/1；分析框行结构 4+4、3、2、1。
# =====================================================================

n.dat <- c(5L, 4L, 3L, 1L)  # 四层各自的数据集框数量
row.n <- list(c(4L, 4L), 3L, 2L, 1L)  # 四层各自逐行的分析框数量
lane.name <- c("Disease bulk",  # 第一层：疾病队列 bulk
               "Rapamycin perturbation bulk",  # 第二层：雷帕霉素扰动 bulk
               "Single-cell",  # 第三层：单细胞
               "Spatial")  # 第四层：空间转录组

dat.h <- DAT_H + (n.dat - 1L) * DAT_STEP  # 各层数据集列总高
box.h <- vapply(row.n, function(r) {
  length(r) * BOX_H + (length(r) - 1L) * ROW_GAP  # 各行分析框加行距
}, 0)
box.w <- vapply(row.n, function(r) {
  n.max <- max(r)  # 本层同行的最多框数
  min((RULE_X2 - ANAL_X0 - (n.max - 1L) * ANAL_GAP) / n.max, BOX_MAX_W)
}, 0)

layer.h <- pmax(dat.h, box.h)  # 层高取两列的高者
layer.top <- FIRST_TOP + c(0, cumsum(layer.h + LANE_GAP)[-length(layer.h)])
dat.top <- layer.top + (layer.h - dat.h) / 2  # 数据集列纵向居中
box.top <- layer.top + (layer.h - box.h) / 2  # 分析框块纵向居中
head.y <- layer.top - HEAD_UP  # 层标题纵坐标
rule.y <- layer.top - RULE_UP  # 层分隔线纵坐标

out.top <- layer.top[length(layer.h)] + layer.h[length(layer.h)] + OUT_GAP
CANVAS_H <- out.top + OUT_H + BOTTOM_MARGIN  # 画布高由版面推出
PAGE_H_IN <- CANVAS_H * UNIT_PT / 72  # 页面高（英寸）

# 各层逐行的框左边界与出线横线纵坐标
box.x <- lapply(seq_along(row.n), function(i) {
  lapply(row.n[[i]], function(n) ANAL_X0 + (seq_len(n) - 1L) * (box.w[i] +
                                                                ANAL_GAP))
})
bus.y <- lapply(seq_along(row.n), function(i) {
  vapply(seq_along(row.n[[i]]), function(r) {
    row.bottom <- box.top[i] + (r - 1L) * (BOX_H + ROW_GAP) + BOX_H
    row.bottom + if (r < length(row.n[[i]])) ROW_GAP / 2 else BUS_OUT
  }, 0)
})

# =====================================================================
# PART 2 图内文字
# 目的：集中登记四层的数据集、分析框与出口文字，便于逐条核对。
# 方法：每个分析框一个列表，含标题、至多三行正文与判定模块编号（支撑分析
#       取 NA）；行宽受框宽约束，故一律用短标签，细节留在图注。
# 数据：样本量与设计取自补充表 S1；网络构建成败取自 S08 拟合结果表；
#       编号与正文 2.1 节的七个判定模块一一对应。
# =====================================================================

dat.lab <- list(
  cbind(c("GSE232911 (n=246)", "GSE57691 (n=68)", "GSE43292 (n=64)",
          "GSE100927 (n=104)", "GSE28829 (n=29)"),
        c("AAA wall, 220 vs 26", "AAA 29 vs 10 control",
          "plaque 32 vs intact", "three arterial beds",
          "early vs advanced")),
  cbind(c("GSE124020 (n=12)", "GSE323983 (n=23)", "GSE6002 (n=4)",
          "GSE65145 (n=6)"),
        c("macrophage anchor", "main perturbation arm", "paired donors",
          "reverse arm (IL-6)")),
  cbind(c("GSE237230 (n=8)", "GSE247238 (n=19)", "GSE212388 (n=8)"),
        c("AAA atlas (10x)", "plaque atlas (10x)", "EC, flow +/- rapa")),
  cbind("GSE277170 (n=56)", "coronary GeoMx")
)

anal.box <- list(
  list(  # 第一层：上行是读数与判定主链，下行是结构性与支撑性分析
    list(
      list(title = "Segment scoring", mod = NA,
           lines = c("five segments", "246 samples", "3 derived indices")),
      list(title = "Panel-level DE", mod = 1L,
           lines = c("limma, per series", "12-gene panel", "BH FDR")),
      list(title = "Stage aggregation", mod = 2L,
           lines = c("per-segment log2FC", "Cohen's d", "30 segment rows")),
      list(title = "Meta-analysis", mod = 3L,
           lines = c("random-effects SMD", "heterogeneity I2",
                     "sign agreement"))
    ),
    list(
      list(title = "Co-expression", mod = NA,
           lines = c("WGCNA (unsigned)", "4 of 6 series",
                     "22/19/12/5 modules")),
      list(title = "Deconvolution", mod = NA,
           lines = c("BayesPrism, EPIC", "quanTIseq", "median consensus")),
      list(title = "Stratification", mod = NA,
           lines = c("region, stage, burden", "18 strata levels",
                     "Wilcoxon, Holm")),
      list(title = "Fraction proxy", mod = 7L,
           lines = c("125 fraction rows", "up/down/indeterm.",
                     "accounting only"))
    )
  ),
  list(
    list(
      list(title = "Panel-level DE", mod = 1L,
           lines = c("limma, per series", "12-gene panel", "BH FDR")),
      list(title = "Direction arbitration", mod = 4L,
           lines = c("six arms, five rapa", "Wilson 95% CI",
                     "up/down/indeterm.")),
      list(title = "Immune enrichment", mod = 6L,
           lines = c("independent arm", "ssGSEA, fgsea, cameraPR",
                     "C7 + GO/KEGG sets"))
    )
  ),
  list(
    list(
      list(title = "Cell annotation", mod = NA,
           lines = c("Seurat v5, marker panel", "12 lineages, 59 genes",
                     "36 clusters")),
      list(title = "Perturbation readout", mod = 5L,
           lines = c("one arm per drug", "sample-level limma", "n=2 vs 2 per arm"))
    )
  ),
  list(
    list(
      list(title = "Spatial segment mapping", mod = NA,
           lines = c("GeoMx DSP, panel QC", "5 segments x ROI",
                     "no drug arm"))
    )
  )
)

out.title <- "Directional synthesis"  # 综合出口标题
out.lines <- c(
  "cross-layer coverage with uncertainty disclosure",
  "five-segment cascade as the common currency",
  "indeterminate reported as indeterminate at every layer"
)

# 版面按行结构排布，文字按层登记；两者不一致即在落图前失败
stopifnot(identical(lapply(anal.box, function(l) {
  vapply(l, length, 0L)
}), lapply(row.n, function(r) as.integer(r))))

# =====================================================================
# PART 3 绘制
# 目的：在统一坐标系（原点左上、y 向下）内绘制四层与综合出口。
# 方法：圆角矩形、文本、线段、圆四类图元；文字非 ASCII 或超框即报错；
#       进线用汇流竖线加逐框短箭头，出线用逐框下垂加逐行横线汇入右竖线。
# 数据：PART 1 的坐标与 PART 2 的文字。
# =====================================================================

Pt <- function(size) size * UNIT_PT  # 用户单位字号换算为磅

AssertAscii <- function(label) {
  # 非 ASCII 字形在 PDF 设备上会被静默替换成字面字节，故设为硬失败
  if (grepl("[^ -~]", label)) {
    stop(sprintf("图内文字含非 ASCII 字符：%s", label))
  }
}

DrawBox <- function(x, y, w, h, fill, edge, lwd, r = BOX_R) {
  # 自建圆角矩形路径：grid.roundrect 在本机 PDF 设备上画不出正确圆角
  # 四角圆心按顺时针排列，每角以 CORNER_N 点逼近四分之一圆弧
  cx <- c(x + r, x + w - r, x + w - r, x + r)  # 左上、右上、右下、左下
  cy <- c(y + r, y + r, y + h - r, y + h - r)
  a0 <- c(pi, 1.5 * pi, 0, 0.5 * pi)  # 各角起始角（y 向下的坐标系）
  px <- numeric(0)  # 路径横坐标累加
  py <- numeric(0)  # 路径纵坐标累加
  for (k in 1:4) {
    ang <- a0[k] + seq(0, pi / 2, length.out = CORNER_N)
    px <- c(px, cx[k] + r * cos(ang))
    py <- c(py, cy[k] + r * sin(ang))
  }
  grid.polygon(px, py, default.units = "native",
               gp = gpar(fill = fill, col = edge, lwd = Pt(lwd)))
}

DrawText <- function(x, y, label, size, col = COL_BODY, bold = FALSE,
                     anchor = "left") {
  # 文本纵向按行中心定位，横向默认左对齐
  AssertAscii(label)
  grid.text(label, x = unit(x, "native"), y = unit(y, "native"),
            just = c(anchor, "centre"),
            gp = gpar(fontsize = Pt(size), col = col,
                      fontface = if (bold) "bold" else "plain"))
}

CheckFit <- function(label, size, avail) {
  # 文本实宽越过框内宽即报错；渲染溢出不会自行暴露，故设为硬失败
  probe <- textGrob(label, gp = gpar(fontsize = Pt(size)))
  w.in <- convertWidth(grobWidth(probe), "inches", valueOnly = TRUE)
  if (w.in > avail * UNIT_PT / 72) {
    stop(sprintf("文字超出框宽：%s（%.1f > %.1f pt）", label, w.in * 72,
                 avail * UNIT_PT))
  }
}

DrawArrow <- function(x0, y0, x1, y1) {
  # 箭头 = 直线段 + 自绘实心三角头；不使用 arrow()，其头部几何在翻转
  # 纵坐标的视口内不可控（头部会溢出到目标框内）
  dx <- x1 - x0  # 方向分量
  dy <- y1 - y0
  len <- sqrt(dx^2 + dy^2)
  ux <- dx / len  # 单位方向
  uy <- dy / len
  bx <- x1 - ARROW_LEN * ux  # 头部底边中点
  by <- y1 - ARROW_LEN * uy
  grid.segments(x0 = unit(x0, "native"), y0 = unit(y0, "native"),
                x1 = unit(bx, "native"), y1 = unit(by, "native"),
                gp = gpar(col = COL_LINE, lwd = Pt(LINE_W)))
  grid.polygon(x = unit(c(bx - ARROW_HALF * (-uy), x1,
                          bx - ARROW_HALF * uy), "native"),
               y = unit(c(by - ARROW_HALF * ux, y1,
                          by + ARROW_HALF * ux), "native"),
               gp = gpar(fill = COL_LINE, col = COL_LINE,
                         lwd = Pt(LINE_W)))
}

DrawLine <- function(x0, y0, x1, y1) {
  # 无箭头连线
  grid.segments(x0 = unit(x0, "native"), y0 = unit(y0, "native"),
                x1 = unit(x1, "native"), y1 = unit(y1, "native"),
                gp = gpar(col = COL_LINE, lwd = Pt(LINE_W)))
}

DrawBadge <- function(x, y, n) {
  # 判定模块编号：白底描边圆压在框顶边上，数字居中
  grid.circle(x = unit(x, "native"), y = unit(y, "native"),
              r = unit(BADGE_R * UNIT_PT / 72, "inches"),
              gp = gpar(fill = COL_WHITE, col = EDGE_ANAL,
                        lwd = Pt(ANAL_W)))
  grid.text(n, x = unit(x, "native"), y = unit(y, "native"),
            gp = gpar(fontsize = Pt(BADGE_SIZE), col = COL_TITLE,
                      fontface = "bold"))
}

cairo_pdf(filename = file.path(FIG_DIR, FIG_FILE),
          width = PAGE_W_IN, height = PAGE_H_IN, family = "serif",
          bg = COL_WHITE)
grid.newpage()
pushViewport(viewport(xscale = c(0, CANVAS_W), yscale = c(CANVAS_H, 0)))

grid.rect(gp = gpar(fill = COL_WHITE, col = NA))  # 画布底色

# --- 图例：四类图元，第二项带编号圆标 ---
legend.x <- c(20, 132, 330, 520)  # 四个色块的左边界
legend.fill <- c(FILL_DAT, FILL_ANAL, FILL_SUP, FILL_OUT)  # 与四类框对应
legend.edge <- c(EDGE_DAT, EDGE_ANAL, EDGE_SUP, EDGE_OUT)  # 图例描边
legend.lab <- c("dataset (GEO)", "verdict module (1-7)",
                "supporting analysis", "integration output")  # 图例文字
for (i in seq_along(legend.x)) {
  DrawBox(legend.x[i], 14, 14, 12, legend.fill[i], legend.edge[i], 1.5, r = 2)
  DrawText(legend.x[i] + 20, 20, legend.lab[i], SIZE_LEGEND)
  CheckFit(legend.lab[i], SIZE_LEGEND, 150)
}
DrawBadge(legend.x[2] + 7, 14, "1")  # 图例示例编号

# --- 四层：标题、分隔线、数据集框、分析框、进线与出线 ---
for (i in seq_along(n.dat)) {
  DrawText(MARGIN_X, head.y[i], lane.name[i], SIZE_LANE, COL_TITLE, bold = TRUE)
  DrawLine(MARGIN_X, rule.y[i], RULE_X2, rule.y[i])  # 层分隔线

  # 数据集框：标题加粗、副标题常规
  for (j in seq_len(n.dat[i])) {
    box.y <- dat.top[i] + DAT_STEP * (j - 1L)  # 本行框顶
    DrawBox(DAT_X, box.y, DAT_W, DAT_H, FILL_DAT, EDGE_DAT, FRAME_W)
    DrawText(DAT_X + TEXT_PAD, box.y + DAT_TITLE_DY, dat.lab[[i]][j, 1],
             SIZE_BODY, COL_TITLE, bold = TRUE)
    DrawText(DAT_X + TEXT_PAD, box.y + DAT_SUB_DY, dat.lab[[i]][j, 2],
             SIZE_BODY, COL_SUB)
    CheckFit(dat.lab[[i]][j, 1], SIZE_BODY, DAT_W - 2 * TEXT_PAD)
    CheckFit(dat.lab[[i]][j, 2], SIZE_BODY, DAT_W - 2 * TEXT_PAD)
  }

  # 进线汇流竖线：纵跨本层首末框心，由数据集列水平引入
  lane.cy <- c(box.top[i] + BOX_H / 2,
               box.top[i] + (length(row.n[[i]]) - 1L) * (BOX_H + ROW_GAP) +
                 BOX_H / 2)
  DrawLine(DAT_X + DAT_W, layer.top[i] + layer.h[i] / 2, BUS_X,
           layer.top[i] + layer.h[i] / 2)
  DrawLine(BUS_X, lane.cy[1], BUS_X, lane.cy[2])

  # 分析框：标题 + 正文 + 编号圆标；逐框由汇流竖线引入
  for (r in seq_along(row.n[[i]])) {
    row.top <- box.top[i] + (r - 1L) * (BOX_H + ROW_GAP)
    for (k in seq_len(row.n[[i]][r])) {
      b <- anal.box[[i]][[r]][[k]]
      bx <- box.x[[i]][[r]][k]  # 本框左边界
      fill <- if (is.na(b$mod)) FILL_SUP else FILL_ANAL  # 支撑分析与模块
      edge <- if (is.na(b$mod)) EDGE_SUP else EDGE_ANAL
      DrawBox(bx, row.top, box.w[i], BOX_H, fill, edge, ANAL_W)
      # 进线只接本行首框：逐框各画一条会横穿同行的前序框
      if (k == 1L) {
        DrawArrow(BUS_X, row.top + BOX_H / 2, bx - 4, row.top + BOX_H / 2)
      }
      DrawText(bx + TEXT_PAD, row.top + TITLE_DY, b$title, SIZE_TITLE,
               COL_TITLE, bold = TRUE)
      CheckFit(b$title, SIZE_TITLE, box.w[i] - 2 * TEXT_PAD)
      for (m in seq_along(b$lines)) {
        DrawText(bx + TEXT_PAD, row.top + LINE_DY[m], b$lines[m],
                 SIZE_BODY, COL_BODY)
        CheckFit(b$lines[m], SIZE_BODY, box.w[i] - 2 * TEXT_PAD)
      }
      if (!is.na(b$mod)) {
        DrawBadge(bx + BADGE_DX, row.top, as.character(b$mod))
      }
      # 出线：本框底部下垂至该行横线
      DrawLine(bx + box.w[i] / 2, row.top + BOX_H, bx + box.w[i] / 2,
               bus.y[[i]][r])
    }
    # 该行横线自首框心右行至右侧汇流竖线
    DrawLine(box.x[[i]][[r]][1] + box.w[i] / 2, bus.y[[i]][r], COL_X,
             bus.y[[i]][r])
  }
}

# --- 右侧汇流竖线与指向综合出口的箭头 ---
DrawLine(COL_X, bus.y[[1]][1], COL_X, out.top + OUT_H / 2)
DrawArrow(COL_X, out.top + OUT_H / 2, OUT_X + OUT_W + 5, out.top + OUT_H / 2)

# --- 综合出口 ---
DrawBox(OUT_X, out.top, OUT_W, OUT_H, FILL_OUT, EDGE_OUT, OUT_W_LINE,
        r = OUT_R)
DrawText(OUT_X + OUT_W / 2, out.top + OUT_TITLE_DY, out.title, SIZE_OUT_TITLE,
         COL_TITLE, bold = TRUE, anchor = "centre")
CheckFit(out.title, SIZE_OUT_TITLE, OUT_W - 2 * TEXT_PAD)
for (m in seq_along(out.lines)) {
  DrawText(OUT_X + OUT_W / 2, out.top + OUT_TITLE_DY + 24 + (m - 1) * 20,
           out.lines[m], SIZE_BODY, COL_BODY, anchor = "centre")
  CheckFit(out.lines[m], SIZE_BODY, OUT_W - 2 * TEXT_PAD)
}

popViewport()
invisible(dev.off())
