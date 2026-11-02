# ===============================================================
# RARP-06_layer1.R
# 功能: 层 1：12 基因面板单基因方向（步骤规格 §8：S05）
# 定位: 层1 · 12 基因面板单基因方向（S05）
# 可视化: Layer1 Panel FoldChange By Contrast；Layer1 EffectSize And Significance Tracks
# ===============================================================
PVAL_CUTOFF <- 0.05  # BH 校正后阈值
LOGFC_MIN <- 0.585  # |log2FC| >= 0.585（= 1.5 倍表达量）
MIN_ARM_N <- 10L  # §9 Q12：n < 10 的臂不进主聚合

vec.panel12 <- c("IL6", "IL1B", "TNF", "CCL2", "CXCL8", "MMP9",
  "SPP1", "NLRP3", "IL10", "TGFB1", "MERTK", "GAS6")

# 每系列报告的对比（**显式声明**，取代旧写法 grep(paste0("^", term), colnames)[1L]）。
# ① 位置取列 = 因子水平**字母序**的副产物，与方案意图无关。GSE57691 的 `~ 0 + group`
# 下第一列是 `groupAOD`，其系数≈AOD 组均值而非任何组间对比；GSE323983 下第一列是
# `IR vs IR+Rapa`，同样不是「扰动 vs 对照」。
# ②对 GSE57691 明文要求「三对比显式声明」，单系数实现无法表达。
# 口径：分子 = 疾病 / 扰动侧，分母 = 对照 / 参考侧；`primary` 标记进 S13 方向投票的那一个。
# 各级别字符串逐字取自实测（S01 的 sample_axis.tsv 与 design_keys.tsv），
# 证据 01_data/verify/check_s05_design_feasibility.R。
dat.contrasts <- rbind(
  data.frame(acc = "GSE232911", term = "tunica",
    contrast = "Media_vs_Adventitia",
    num = "Media", den = "Adventitia", primary = TRUE,
    note = "组织层间对比（该系列无对照侧）", stringsAsFactors = FALSE),
  data.frame(acc = "GSE57691", term = "group",
    contrast = c("largeAAA_vs_control", "smallAAA_vs_control",
      "AOD_vs_control"),
    num = c("large AAA", "small AAA", "AOD"), den = "control",
    primary = c(TRUE, FALSE, FALSE),
    note = "方案 §5 要求三对比显式声明", stringsAsFactors = FALSE),
  data.frame(acc = "GSE43292", term = "tissue", contrast = "Intact_vs_Plaque",
    num = "Macroscopically intact tissue", den = "Atheroma plaque",
    primary = TRUE, note = "", stringsAsFactors = FALSE),
# 即该轴的实质读数；primary 取颈动脉（病变 29 / 对照 12，样本量最大）。
  data.frame(acc = "GSE100927", term = "tissue",
    contrast = c("carotid_lesion_vs_control", "femoral_lesion_vs_control",
      "popliteal_lesion_vs_control"),
    num = c("Atherosclerotic carotid artery",
      "Atherosclerotic femoral artery",
      "Atherosclerotic infra-popliteal artery"),
    den = c("Control carotid artery", "Control femoral artery",
      "Control infra-popliteal artery"),
    primary = c(TRUE, FALSE, FALSE),
    note = "血管床轴：三床各自病变 vs 对照", stringsAsFactors = FALSE),
  data.frame(acc = "GSE28829", term = "phenotype",
    contrast = "advanced_vs_early_plaque",
    num = "advanced atherosclerotic plaque",
    den = "early atherosclerotic plaque", primary = TRUE,
    note = "进展轴：晚期 vs 早期斑块", stringsAsFactors = FALSE),
  data.frame(acc = "GSE323983", term = "arm",
    contrast = c("IR_vs_CTRL", "Rapa_vs_CTRL", "RapaIR_vs_CTRL",
      "RapaIR_vs_IR"),
    num = c("IR", "Rapa", "IR with Rapa", "IR with Rapa"),
    den = c("Ctrl", "Ctrl", "Ctrl", "IR"),
    primary = c(FALSE, FALSE, FALSE, TRUE),
    note = "主读出 = RapaIR_vs_IR（Rapa 在 IR 背景下的额外效应）",
    stringsAsFactors = FALSE),
# vs Cy3=LPS（SOFT 官方字段）—— 每样本**本身已是**一个「Rapa 效应」配对读值，
# 对照已编码在 M 值内，故无独立对照臂。对比 = 4 donor 效应均值 vs 0，
# den 记 "0" 表单侧（对比构造见 PART 2 的 den=="0" 分支）。
# 设计取 `effect` **单水平**（4 样本同为 "RapaEffect"）：`~ 0 + effect` 得
# 1 列设计 → 残差 df = 4 - 1 = 3 可估。若把 4 个 GSM 当 4 个水平则 4 列 = 4
# 样本 → df = 0，eBayes 无自由度（实测 01_data/verify/segfault_bisect.R：
# 全规模 × 全 NA 比例均报 "No residual degrees of freedom"）。
# n=4 功效下限 → primary=FALSE（不进 S13 方向投票、不进 4/7 达标正典）。
  data.frame(acc = "GSE6002", term = "effect", contrast = "RapaEffect_vs_0",
    num = "RapaEffect", den = "0", primary = FALSE,
    note = paste0("4 donor 双色 M 值均值 vs 0（对照已编码在 M 内）；n=4",
      "功效下限，补充臂"),
    stringsAsFactors = FALSE)
)


TAB_DIR <- "03_results/S05/table"
FIG_DIR <- "03_results/S05/Figure"
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

library(limma)  # 层 1 单基因 limma 拟合（eBayes + topTable）
library(ggplot2)  # 层 1 方向图；主题与出口用 §sci_theme / §sci_save

# 主题： §sci_theme 官方定义（theme_bw + 关网格），全脚本统一口径
sci_theme <- ggplot2::theme_bw(base_size = 12) +
ggplot2::theme(
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank(),
  strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
  legend.position  = "right",
  legend.text      = ggplot2::element_text(size = 10, face = "plain"),
  axis.title       = ggplot2::element_text(size = 12, face = "plain"),
  axis.text        = ggplot2::element_text(size = 10, face = "plain"),
  plot.title       = ggplot2::element_text(size = 14, hjust = 0.5,
    face = "plain"),
  legend.title     = ggplot2::element_text(size = 11, face = "plain")
)

# 出口： §sci_save 官方定义，仅 PDF，落本步 FIG_DIR。
# 可读性核验钩子（ · 交付前须核验可读性）：仅当 RAPA_QC_PNG=1 时额外导 PNG
# 到 FIG_DIR/_qc/ 供核验；默认不产出，交付物仍为纯 PDF（符合 §sci_save）。
sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  pdf(pdf_file, width = width, height = height)
  print(plot_obj)
  dev.off()
  if (Sys.getenv("RAPA_QC_PNG") == "1") {
    dir.create(file.path(FIG_DIR, "_qc"), showWarnings = FALSE)
    png(file.path(FIG_DIR, "_qc", paste0(filename, ".png")),
      width = width * 150, height = height * 150, res = 150)
    print(plot_obj)
    dev.off()
  }
}

# PART 1 每系列的设计映射（键名必须与 S01 design_keys.tsv 中的实测键名对应） ---------------

# 注意：GSE232911 分组表里的主体列名为 `subject`（S01 按 title 解析），不是 `patient`。
dat.models <- rbind(
  data.frame(acc = "GSE232911", formula = "~ subject + tunica + thrombus",
    term = "tunica", key.pattern = "^tunica$", factor.vars = "subject",
    stringsAsFactors = FALSE),
  data.frame(acc = "GSE57691", formula = "~ 0 + group", term = "group",
    key.pattern = "^group$|^disease state$", factor.vars = "",
    stringsAsFactors = FALSE),
  data.frame(acc = "GSE43292", formula = "~ patient + tissue", term = "tissue",
    key.pattern = "^tissue$|^disease state$", factor.vars = "patient",
    stringsAsFactors = FALSE),
# 六水平 = 3 血管床 x {病变, 对照}。无主体标识（S04 实测 Q5/Q6 均 N/A），
# 故设计不含配对协变量。
  data.frame(acc = "GSE100927", formula = "~ tissue", term = "tissue",
    key.pattern = "^tissue$", factor.vars = "", stringsAsFactors = FALSE),
  data.frame(acc = "GSE28829", formula = "~ phenotype", term = "phenotype",
    key.pattern = "^phenotype$", factor.vars = "", stringsAsFactors = FALSE),
  data.frame(acc = "GSE323983", formula = "~ donor + arm", term = "arm",
    key.pattern = "^arm$|^treatment$|^condition$", factor.vars = "donor",
    stringsAsFactors = FALSE),
# GSE6002（4-donor 双色）：样本即效应（M=log2(Cy5/Cy3)，Cy5=LPS+Rapa vs Cy3=LPS）。
# `effect` 取**单水平**（4 样本同为 "RapaEffect"）→ `~ 0 + effect` 得 1 列设计，
# 残差 df = 4 - 1 = 3 可估（= 4 个配对值的 one-sample 检验）。
# 对比侧（den = "0"）在 PART 2 单侧分支只给 num 计权，即「均值是否 ≠ 0」。
# 方案 §487 定位「补充臂（非主聚合）」→ primary = FALSE，不进 S13 方向投票。
  data.frame(acc = "GSE6002", formula = "~ 0 + effect", term = "effect",
    key.pattern = "^effect$", factor.vars = "",
    stringsAsFactors = FALSE)
)

dat.contrasts$formula <- dat.models$formula[match(dat.contrasts$acc,
    dat.models$acc)]
if (anyNA(dat.contrasts$formula)) {
  stop(paste0("S05: 对比注册表中有系列未在 dat.models 登记公式：",
      paste(unique(dat.contrasts$acc[is.na(dat.contrasts$formula)]),
        collapse = ", ")))
}

dat.contrasts$factor.vars <- dat.models$factor.vars[match(dat.contrasts$acc,
    dat.models$acc)]
if (anyNA(dat.contrasts$factor.vars)) {
  stop(paste0("S05: 对比注册表中有系列未在 dat.models 登记 factor.vars：",
      paste(unique(dat.contrasts$acc[is.na(dat.contrasts$factor.vars)]),
        collapse = ", ")))
}

dat.keys <- read.delim("03_results/S01/table/design_keys.tsv",
  stringsAsFactors = FALSE)

# 按系列把 characteristics 展开成样本级设计表
BuildDesign <- function(acc) {
  if (acc == "GSE232911") {
    return(read.delim("03_results/S01/table/gse232911_grouping.tsv",
        stringsAsFactors = FALSE))
  }
# 无 S01 设计键（4-donor 配对未登记），设计表由表达矩阵列名直接构造。
# `effect` 取**单水平常量** "RapaEffect"（4 样本同一水平）：`~ 0 + effect`
# 给 1 列设计 → 残差 df = 3 可估。若取样本 ID 作 4 个水平则 df = 0 不可解
# （实测 01_data/verify/segfault_bisect.R）。
# n=4 功效下限如实披露：S11 侧单样本 t + Holm 实测
# （01_data/verify/probe_gse6002_semantics.R：CCL2 p.Holm=0.064 / TNF 0.100 /
# IL6 0.306 / IL10 0.343 / IL1B 0.560，无一 Holm 显著）。
  if (acc == "GSE6002") {
    obj <- readRDS(file.path("03_results", "S03", "table",
      paste0(acc, "_gene.rds")))
    return(data.frame(sample_id = colnames(obj$expr),
      effect = rep("RapaEffect", ncol(obj$expr)),
      stringsAsFactors = FALSE))
  }
  if (acc == "GSE323983") {
    fp <- "03_results/S01/table/sample_axis.tsv"
    if (!file.exists(fp)) stop(paste0("S05: 需要 ", fp, "（请先跑 S01）。"))
    dat.ax <- read.delim(fp, stringsAsFactors = FALSE)
    dat.ax <- dat.ax[dat.ax$acc == acc, c("sample_id", "donor", "arm"),
      drop = FALSE]
    if (nrow(dat.ax) == 0L) {
      stop(paste0("S05: sample_axis.tsv 中缺 ", acc, " 的行。"))
    }
    return(dat.ax)
  }
  m <- dat.models[dat.models$acc == acc, ]
  dat.k <- dat.keys[dat.keys$acc == acc, ]
  vec.key <- unique(dat.k$key)
  hit.name <- grep(paste0("^", m$term, "$"), vec.key, ignore.case = TRUE)
  if (length(hit.name) > 0L) {
    key.use <- vec.key[hit.name[1L]]
  } else {
    hit <- grep(m$key.pattern, vec.key, ignore.case = TRUE)
    if (length(hit) != 1L) {
      stop(paste0("S05: ", acc, " 的 term '", m$term,
          "' 无法唯一定位设计键；key.pattern '",
          m$key.pattern, "' 命中 ", length(hit), " 个键（须恰为 1 个）：",
          paste(vec.key[hit], collapse = " | "),
          " ; 该系列全部键: ", paste(vec.key, collapse = " | ")))
    }
    key.use <- vec.key[hit[1L]]
  }
  dat.long <- dat.k[dat.k$key == key.use, c("sample_id", "value")]
  names(dat.long) <- c("sample_id", m$term)
  vec.extra <- setdiff(all.vars(as.formula(m$formula)), m$term)
  for (v in vec.extra) {
    hit.v <- grep(paste0("^", v, "$"), vec.key, ignore.case = TRUE)
    if (length(hit.v) == 0L) {
      stop(paste0("S05: ", acc, " 的模型变量 '", v,
          "' 在 design_keys.tsv 中没有对应键；实测键: ",
          paste(vec.key, collapse = " | ")))
    }
    dat.pv <- dat.k[dat.k$key == vec.key[hit.v[1L]], c("sample_id", "value")]
    names(dat.pv) <- c("sample_id", v)
    dat.long <- merge(dat.long, dat.pv, by = "sample_id", all = FALSE)
  }
# 合并后行数必须等于该系列的样本数；少行说明有样本缺某个协变量键 —— 那属于
# 静默丢样本，必须报错而不是让 model.frame 悄悄 na.omit。
  n.expect <- length(unique(dat.k$sample_id))
  if (nrow(dat.long) != n.expect) {
    stop(paste0("S05: ", acc, " 的设计表合并后 ", nrow(dat.long),
        " 行 != 样本数 ",
        n.expect, "；有样本缺协变量键，不得静默丢样本。"))
  }
  dat.long$sample_id <- as.character(dat.long$sample_id)
  dat.long
}

# PART 2 逐个系列跑面板 -----------------------------------------------------------------------

dat.layer1 <- data.frame(acc = character(0), term = character(0),
  contrast = character(0),
  gene = character(0), log2FC = numeric(0), adj.P.Val = numeric(0),
  direction = character(0), significant = character(0),
  n_num = integer(0), n_den = integer(0),
  primary = logical(0), exploratory = logical(0),
  stringsAsFactors = FALSE)
dat.panel.miss <- data.frame(acc = character(0), gene = character(0),
  stringsAsFactors = FALSE)
vec.deferred <- character(0)
# df=0 的 4-donor M 值不是「缺数据」而是「统计上不可解」，
# 主读数由 S11 单样本 t + Holm 消费（凭据 probe_gse6002_semantics.R）。
vec.df0 <- character(0)

for (acc in unique(dat.models$acc)) {
# 探针级 _expr.rds 的双轴身份均被替换（p1..pN / s1..sN），面板按名匹配与设计表
# 按 sample_id 连接都会落空（实测见 01_data/verify/check_s03_rownames_contract.R）。
  file.expr <- file.path("03_results", "S03", "table", paste0(acc, "_gene.rds"))
  if (!file.exists(file.expr)) {
    vec.deferred <- c(vec.deferred, acc)
    next
  }
  expr <- readRDS(file.expr)$expr
  dat.design <- BuildDesign(acc)
  dat.design <- dat.design[dat.design$sample_id %in% colnames(expr), ,
    drop = FALSE]
  rownames(dat.design) <- dat.design$sample_id
  m <- dat.models[dat.models$acc == acc, ]
  vec.fv <- strsplit(m$factor.vars, ";", fixed = TRUE)[[1L]]
  vec.fv <- vec.fv[nzchar(vec.fv)]
  for (v in vec.fv) {
    if (!v %in% colnames(dat.design)) {
      stop(paste0("S05: ", acc, " 声明的分类协变量 '", v,
          "' 不在设计表内（实测列: ",
          paste(colnames(dat.design), collapse = ", "), "）。"))
    }
    dat.design[[v]] <- factor(dat.design[[v]])
  }
  dat.cmp <- dat.contrasts[dat.contrasts$acc == acc, , drop = FALSE]
  if (nrow(dat.cmp) == 0L) {
    stop(paste0("S05: ", acc,
        " 在 dat.contrasts 里没有登记任何对比；禁止静默跳过。"))
  }

# 面板基因按行名（基因符号）精确匹配。§9 Q4 的本意是「缺哪个报哪个、禁止同义替换」，
# 不整体剔除该系列（GSE323983 缺 TNF/MMP9/NLRP3/IL10，实测见
# 01_data/verify/check_panel_presence.R）。
  idx.gene <- match(vec.panel12, rownames(expr))
  if (anyNA(idx.gene)) {
    vec.miss <- vec.panel12[is.na(idx.gene)]
    cat(sprintf(
        "[S05] %-11s panel missing %d/%d: %s -- recorded, series continues\n",
        acc, length(vec.miss), length(vec.panel12), paste(vec.miss,
          collapse = ", ")))
    dat.panel.miss <- rbind(dat.panel.miss, data.frame(acc = acc,
        gene = vec.miss,
        stringsAsFactors = FALSE))
  }
  vec.gene.use <- vec.panel12[!is.na(idx.gene)]
  mat.expr <- expr[idx.gene[!is.na(idx.gene)], , drop = FALSE]
  rownames(mat.expr) <- vec.gene.use

  dat.f <- model.frame(as.formula(m$formula), data = dat.design)
  tt.design <- attr(dat.f, "terms")
  dat.f[[m$term]] <- droplevels(as.factor(dat.f[[m$term]]))
  tab.arm <- table(dat.f[[m$term]])
# 单水平说明设计键解析错了（例如键 pattern 命中到 all-constant 的 characteristics），
# 与其让 model.matrix 抛 "contrasts can be applied only to factors with 2 or more
# levels"，不如在此给出可诊断的报错。
# 单侧对比系列（GSE6002：den = "0"，与 0 比）**本就只有 1 个水平**，不在此拦。
  lgl.onesided <- all(dat.cmp$den == "0")
  if (nlevels(dat.f[[m$term]]) < 2L && !lgl.onesided) {
    stop(paste0("S05: ", acc, " 的处理项 '", m$term, "' 只有 ",
        nlevels(dat.f[[m$term]]), " 个水平（取值: ",
        paste(unique(dat.f[[m$term]]), collapse = ", "),
        "），无法做组间对比；检查 PART 1 的键解析。"))
  }
  dat.f <- dat.f[, vapply(dat.f, function(x) !is.numeric(x) || length(unique(
          x)) > 1L, TRUE),
    drop = FALSE]
# GSE6002：单水平因子在 model.matrix 下无法设 contrasts（R 硬限制，实测报
# "contrasts can be applied only to factors with 2 or more levels"），
# 故显式构造 1 列设计矩阵（4 donor 共同均值，等价 one-sample t），
# 列名沿用 `paste0(term, level)` = "effectRapaEffect"，与下方对比构造无缝衔接。
  if (acc == "GSE6002") {
    mat.design <- matrix(1, nrow = nrow(dat.f), ncol = 1L,
      dimnames = list(rownames(dat.f), paste0(m$term, "RapaEffect")))
  } else {
    mat.design <- stats::model.matrix(tt.design, data = dat.f)
  }

  mat.expr <- mat.expr[, rownames(dat.f), drop = FALSE]

# GSE6002（4 donor × `~ 0 + effect`）设计矩阵 4 列 = 4 样本 → df = 0，
# eBayes 无自由度可估（实测 01_data/verify/segfault_bisect.R 全规模/全 NA 比例
# 均报 "No residual degrees of freedom"，属数学不可解而非段错误）。
# 命中即 defer（记 vec.deferred，不静默崩、不阻断其他系列）——
# 把「4 donor 单对比功效下限」变成可追溯的显式降级。
  if (ncol(mat.design) >= ncol(mat.expr)) {
    vec.df0 <- c(vec.df0, acc)
    cat(sprintf("[S05] %-11s df=0 (n_col_design = %d >= n_col_expr = %d) -> deferred; ",
        acc, ncol(mat.design), ncol(mat.expr)))
    cat(paste0("residual df = 0, eBayes cannot estimate; S11",
      "single-sample t + Holm consumes directly.\n"))
    next
  }

# 对比矩阵：按 dat.contrasts 逐条显式声明构造，不再按字母序取第一列。
# 被吸收为基准的水平（有截距的模型里）在设计矩阵中没有对应列，此时该侧不加权。
  vec.lv <- levels(dat.f[[m$term]])
  vec.col <- paste0(m$term, vec.lv)
  mat.con <- matrix(0, nrow = ncol(mat.design), ncol = nrow(dat.cmp),
    dimnames = list(colnames(mat.design), dat.cmp$contrast))
  for (j in seq_len(nrow(dat.cmp))) {
# 单侧对比（GSE6002，den = "0"）：与 0 比，只在 num 侧计权
# —— 语义 = 「4 donor M 值均值是否显著 ≠ 0」（对照已编码在 M 值内）。
    if (identical(dat.cmp$den[j], "0")) {
      vec.col.j <- paste0(m$term, dat.cmp$num[j])
      if (!vec.col.j %in% colnames(mat.design)) {
        stop(paste0("S05: ", acc, " / ", dat.cmp$contrast[j],
            " 的 num 列 '", vec.col.j, "' 不在设计矩阵内（实测列: ",
            paste(colnames(mat.design), collapse = ", "), "）。"))
      }
      mat.con[vec.col.j, j] <- 1
      next
    }
    vec.lv.j <- c(dat.cmp$num[j], dat.cmp$den[j])
    i.j <- match(vec.lv.j, vec.lv)
    if (anyNA(i.j)) {
      stop(paste0("S05: ", acc, " / ", dat.cmp$contrast[j],
          " 的对比水平不存在: ",
          paste(vec.lv.j, collapse = " vs "), " ; 实测水平: ",
          paste(vec.lv, collapse = " | ")))
    }
    lgl.has <- vec.col[i.j] %in% colnames(mat.design)
    if (!any(lgl.has)) {
      stop(paste0("S05: ", acc, " / ", dat.cmp$contrast[j],
          " 的两侧都被吸收为基准，无法估计。"))
    }
    mat.con[vec.col[i.j][lgl.has], j] <- c(1, -1)[lgl.has]
  }

  fit <- limma::lmFit(mat.expr, design = mat.design)
  fit <- limma::contrasts.fit(fit, mat.con)
  fit <- limma::eBayes(fit)

  for (j in seq_len(nrow(dat.cmp))) {
    dat.tt <- limma::topTable(fit, coef = dat.cmp$contrast[j], number = Inf,
      adjust = "BH")
    dat.tt$gene <- rownames(dat.tt)
# 方向分类：不显著 = indeterminate（§4.2，禁止回退为 0）
    dat.tt$direction <- ifelse(dat.tt$adj.P.Val >= PVAL_CUTOFF, "indeterminate",
      ifelse(dat.tt$logFC >= LOGFC_MIN, "up",
        ifelse(dat.tt$logFC <= -LOGFC_MIN, "down",
          "indeterminate")))
    dat.tt$significant <- ifelse(dat.tt$adj.P.Val < PVAL_CUTOFF &
      abs(dat.tt$logFC) >= LOGFC_MIN, "YES", "NO")
    i.num <- match(dat.cmp$num[j], vec.lv)
    n.num <- as.integer(tab.arm[i.num])
# 单侧对比（den = "0"）无对照臂：n_den = 0 → 必然判 exploratory，
# 正合「n=4 功效下限不进主聚合」的语义（若走 match 得 NA，min() 变 NA
# 会让后续 if(NA) 报错，故显式给 0）。
    n.den <- if (identical(dat.cmp$den[j], "0")) 0L
      else as.integer(tab.arm[match(dat.cmp$den[j], vec.lv)])
    lgl.small <- min(n.num, n.den) < MIN_ARM_N
    if (lgl.small) {
      cat(sprintf("[S05] %-11s %-16s arms n = %d vs %d < %d -> exploratory\n",
          acc, dat.cmp$contrast[j], n.num, n.den, MIN_ARM_N))
    }
    dat.layer1 <- rbind(dat.layer1, data.frame(
        acc = acc, term = m$term, contrast = dat.cmp$contrast[j],
        gene = dat.tt$gene,
        log2FC = dat.tt$logFC, adj.P.Val = dat.tt$adj.P.Val,
        direction = dat.tt$direction, significant = dat.tt$significant,
        n_num = n.num, n_den = n.den,
        primary = dat.cmp$primary[j], exploratory = lgl.small,
        stringsAsFactors = FALSE))
  }
  cat(sprintf("[S05] %-11s term=%-7s contrasts=%d genes=%d sig=%d\n", acc,
      m$term,
      nrow(dat.cmp), sum(dat.layer1$acc == acc),
      sum(dat.layer1$significant[dat.layer1$acc == acc] == "YES")))

  dat.design.out <- as.data.frame(dat.f)
  dat.design.out$sample_id <- rownames(dat.f)
  write.csv(dat.design.out, file.path(TAB_DIR, paste0("design_used_", acc,
        ".csv")),
    row.names = FALSE)
}

write.csv(dat.layer1, file.path(TAB_DIR, "layer1_gene_dir.csv"),
  row.names = FALSE)
write.csv(dat.panel.miss, file.path(TAB_DIR, "panel_missing.csv"),
  row.names = FALSE)
# 对比注册表落盘，作为层 2（S06）/层 3（S07）选对比的**单一事实源**：
# 它们原先各自用 grep("^term")[1L] 取设计矩阵第一列，与层 1 的 P0-3 是同一类位置取列缺陷。
write.csv(dat.contrasts, file.path(TAB_DIR, "contrast_spec.csv"),
  row.names = FALSE)

# PART 2.5 图件 -------------------------------------------------------------------------------
# 命名遵守 16（描述性英文，禁 Fig 前缀）；图内文字全英文（ 17）；
# 一律单面板，不引入任何拼图包。
# 原分面柱状图会撑到 9 个面板、且无处安放「该平台未注释该基因」这一事实，肉眼不可判读。
# 改为 图 1 单面板热图（含灰格表示未注释）、图 2 单面板柱（按是否探索性着色）。

dat.l1.plot <- dat.layer1
# 对比标签单行化（ 17 禁多行标签），用 "系列: 对比" 形式
dat.l1.plot$cmp <- paste0(dat.l1.plot$acc, ": ", dat.l1.plot$contrast)
vec.cmp <- unique(dat.l1.plot$cmp)
# 完整网格：缺失的面板基因也要占格（灰底），否则「没画出来」与「测了不显著」无法区分
dat.grid <- expand.grid(cmp = vec.cmp, gene = vec.panel12,
  stringsAsFactors = FALSE)
i.grid <- match(paste(dat.grid$cmp, dat.grid$gene),
  paste(dat.l1.plot$cmp, dat.l1.plot$gene))
dat.grid$log2FC <- dat.l1.plot$log2FC[i.grid]
dat.grid$significant <- dat.l1.plot$significant[i.grid]
dat.grid$exploratory <- dat.l1.plot$exploratory[i.grid]
dat.grid$cmp <- factor(dat.grid$cmp, levels = vec.cmp)
# 基因轴按预注册面板顺序（vec.panel12），非字母序，保证读图顺序与方案一致。
dat.grid$gene <- factor(dat.grid$gene, levels = vec.panel12)

p.fc <- ggplot(dat.grid, aes(x = gene, y = cmp, fill = log2FC)) +
geom_tile(color = "white", linewidth = 0.4) +
# 既不遮挡 log2FC 的双向色，也能一眼看出显著格；判据阈值不变。
geom_text(data = dat.grid[!is.na(dat.grid$significant) &
    dat.grid$significant == "YES", , drop = FALSE],
  label = "*", size = 5.5, vjust = 0.75, colour = "grey10") +
# 缺测格的灰必须与「近零 log2FC 的白」拉开距离：grey85 与中点白几乎同色，
# 读者无法区分「该平台未注释」与「测了但无差异」（ 可判读性）。
scale_fill_gradient2(low = "#3B6BA5", mid = "white", high = "#C44E52",
  midpoint = 0, na.value = "grey60",
  name = expression(log[2]~FC)) +
# y 轴离散扩展：9 个对比标签上下各留一点余量，避免首尾贴边。
scale_y_discrete(limits = rev(vec.cmp), expand = expansion(add = 0.6)) +
# 副标题显式换行：单行 89 字符在 9 in 画布上右端溢出被裁（与 caption 同类根因）。
# 画布加宽到 10.5 in 同时解决 x 轴末三基因标签（TGFB1/MERTK/GAS6）视觉粘连。
labs(title = "Panel-gene fold change by declared contrast",
  subtitle = paste0(
    "Grey = symbol not annotated on that platform (not measured)\n",
    "Asterisk = BH < 0.05 and |log2FC| >= 0.585"),
  x = "Panel gene", y = "Series: contrast") +
sci_theme
sci_save(p.fc, "Layer1_Panel_FoldChange_By_Contrast", width = 10.5, height = 7)

dat.l1.plot$pass.fc <- !is.na(dat.l1.plot$log2FC) & abs(dat.l1.plot$log2FC) >=
LOGFC_MIN
dat.l1.count <- aggregate(list(n_sig = dat.l1.plot$significant == "YES",
    n_fc = dat.l1.plot$pass.fc),
  by = list(cmp = dat.l1.plot$cmp), FUN = sum)
dat.l1.count$n_gene <- as.integer(table(factor(dat.l1.plot$cmp,
      levels = vec.cmp))[
    match(dat.l1.count$cmp, vec.cmp)])
dat.l1.count$n_only <- dat.l1.count$n_fc - dat.l1.count$n_sig
dat.l1.count$exploratory <- dat.l1.plot$exploratory[match(dat.l1.count$cmp,
    dat.l1.plot$cmp)]
dat.l1.count$cmp <- factor(dat.l1.count$cmp, levels = vec.cmp)
# y 轴标签：探索性对比加 † 后缀（ 17 禁多行标签，故用符号 + caption 释义）
vec.lab.y <- as.character(dat.l1.count$cmp)
vec.lab.y[dat.l1.count$exploratory] <- paste0(vec.lab.y[
    dat.l1.count$exploratory],
  " \u2020")
names(vec.lab.y) <- as.character(dat.l1.count$cmp)

# 堆叠长表：两段 = 双达标 / 仅效应量达标
dat.l1.stack <- rbind(
  data.frame(cmp = dat.l1.count$cmp,
    strata = "BH < 0.05 and |log2FC| >= 0.585",
    n = dat.l1.count$n_sig, stringsAsFactors = FALSE),
  data.frame(cmp = dat.l1.count$cmp,
    strata = "|log2FC| >= 0.585 only",
    n = dat.l1.count$n_only, stringsAsFactors = FALSE))
dat.l1.stack$strata <- factor(dat.l1.stack$strata,
  levels = c("BH < 0.05 and |log2FC| >= 0.585",
    "|log2FC| >= 0.585 only"))

sci_save(ggplot(dat.l1.stack, aes(y = cmp, x = n, fill = strata)) +
  geom_col(width = 0.7) +
# 柱右端标「效应量达标数 / 该平台已测数」：分母差异（GSE323983 为 7 而非 12）
# 必须在读数处可见，不与柱长（绝对数）混淆。
  geom_text(data = dat.l1.count,
    aes(y = cmp, x = n_fc, label = sprintf("%d/%d", n_fc,
        n_gene)),
    hjust = -0.22, size = 3.2, colour = "grey20",
    inherit.aes = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  scale_y_discrete(limits = rev(vec.cmp), labels = vec.lab.y,
    expand = expansion(add = 0.6)) +
  scale_fill_manual(values = c("BH < 0.05 and |log2FC| >= 0.585" = "#4C72B0",
      "|log2FC| >= 0.585 only" = "#A8C4E0"),
    name = "Panel genes") +
# 单行时向左溢出画布左缘被裁。故显式换行 + hjust = 0 左对齐。
  labs(title = "Effect-size and significance tracks per declared contrast",
    subtitle = paste0("Bar length = genes with |log2FC| >= 0.585; ",
      "filled segment = those also passing BH < 0.05\n",
      "Labels = effect-size passers / panel genes measured"),
    caption = paste0(
      "Denominator = 12 panel genes except GSE323983, whose platform ",
      "does not annotate 5 of the 12\n",
      "(IL1B, TNF, MMP9, NLRP3, IL10); ",
      "\u2020 = exploratory (arms n < 10)"),
    x = "Panel genes", y = "Series: contrast") +
  sci_theme +
  theme(plot.caption = element_text(hjust = 0)),
  "Layer1_EffectSize_And_Significance_Tracks", width = 10, height = 6.5)

# PART 3 收口 -------------------------------------------------------------------------------

n.indet <- sum(dat.layer1$direction == "indeterminate")
dat.sig.cmp <- aggregate(list(n_sig = dat.layer1$significant == "YES"),
  by = list(acc = dat.layer1$acc,
    contrast = dat.layer1$contrast),
  FUN = sum)
dat.nosig <- dat.sig.cmp[dat.sig.cmp$n_sig == 0L, c("acc", "contrast")]
message(paste0("\n===== S05 layer1 summary =====\nseries modelled = ",
    length(unique(dat.layer1$acc)),
    "  contrasts = ", length(unique(dat.layer1$contrast)),
    "  deferred (no expr from S03) = ",
    if (length(vec.deferred) == 0L) "none" else paste(vec.deferred,
      collapse = ", "),
    "  df=0 (4-donor M-value, S11 consumes) = ",
    if (length(vec.df0) == 0L) "none" else paste(vec.df0,
      collapse = ", "),
    "\nrows = ", nrow(dat.layer1),
    "  indeterminate = ", n.indet,
    "  (" , round(100 * n.indet / max(nrow(dat.layer1), 1L)), "% )",
    paste0("\npanel genes not measured on some platform (recorded,",
      "not substituted) = "),
    nrow(dat.panel.miss),
    "\nacc x contrast with zero significant panel gene = ", nrow(dat.nosig)))
if (nrow(dat.nosig) > 0L) {
  print(dat.nosig)
}
if (length(vec.deferred) > 0L) {
  stop(paste0("S05 STOP: series without an expression matrix from S03: ",
      paste(vec.deferred, collapse = ", "),
      " . Run S03 for these series first; the script will not skip them."))
}
if (length(vec.df0) > 0L) {
  message(paste0("[S05] df=0 series deferred to S11 (4-donor M-value): ",
      paste(vec.df0, collapse = ", "),
      " . Readout: single-sample t-test + Holm, n=4 power floor.",
      " Evidence file: 01_data/verify/probe_gse6002_semantics.R."))
}
message("S05 OK. tables -> 03_results/S05/table/")
