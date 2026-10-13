#===============================================================
## S04_confound.R —— 混杂与配对诊断（步骤规格 §8：S04）
## 目的：对 S03 产出的表达矩阵做「批次 / 块 / 配对结构」诊断：
##       top 可变基因 PCA，按设计因子着色，估计 surrogate variables，
##       并对 §9 Q5（配对因子完整）、Q6（块因子不主导 PCA）给出实测结论。
##
## 定位：分析脚本。本步只估计与诊断，**不改表达值**；显式建模留在 S05 / S06 的模型内，
##       避免二次校正叠加（规则 25：不做无效守卫与重复处理）。
##
## 产出：03_results/S04/table/sv_object.rds    各系列 SV 对象
##       03_results/S04/table/confound_diag.csv 诊断结论（Q5 / Q6）
##       03_results/S04/Figure/PCA_Coloured_By_Design_Factor_<acc>.pdf 每系列单面板 PCA
##       03_results/S04/Figure/Surrogate_Variables_Boxplot_<acc>.pdf    每系列单面板 SV 箱线
##       03_results/S04/Figure/Confound_Diagnosis_By_Series.pdf         单面板热图：系列 x Q5/Q6
##
## 图件清单（规则 35 可视化强制 · SCI 发表导向；动笔前已调研文献确定风格/形式/数量）：
##   可视化: S04 共 3 类（PCA / SV 每系列各一张 + 诊断图一张）
##     1. PCA_Coloured_By_Design_Factor_<acc>.pdf
##        —— 形式：单面板散点（PC1 vs PC2，按设计因子着色，轴标注方差占比）
##        —— 科学问题：批次 / 块 / 配对结构是否主导样本分布（Q6 的前置观察）。
##        —— 参照文献：PCA 按设计因子着色为批次诊断标准首图（Leek 2010, Nat Rev Genet,
##           doi:10.1038/nrg2825, PMID 20084099）。
##     2. Surrogate_Variables_Boxplot_<acc>.pdf
##        —— 形式：单面板箱线（SV1..SVk 的分布）
##        —— 科学问题：替代变量在样本间是否有足够散布（是否值得作协变量）。
##        —— 参照文献：sva 替代变量诊断（Leek & Storey 2007, PLoS Genet,
##           doi:10.1371/journal.pgen.0030161, PMID 18081495）。
##     3. Confound_Diagnosis_By_Series.pdf
##        —— 形式：单面板热图（x = 系列，y = Q5/Q6 断言，tile 着色 PASS/FAIL，格内写判定值）
##        —— 科学问题：各系列配对完整性（Q5）与块因子（Q6）逐项通过情况汇总。
##        —— 参照文献：断言状态矩阵热图为 QC 汇总标准画法。
##
## 图件口径：sci_theme / sci_save 逐字取自规范 §sci_theme / §sci_save 官方定义；
##       不引入任何拼图包，单面板 ggplot2 出图；可读性核验见 sci_save 之 RAPA_QC_PNG 钩子。
##
## 前置：S03 通过。所需包 sva / variancePartition / RUVSeq（未装则显式 stop）。
#===============================================================

## PART 0 常量与主题 ------------------------------------------------------------------

SEED <- 20260925
set.seed(SEED)
FIG_DIR <- "03_results/S04/Figure"
TAB_DIR <- "03_results/S04/table"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

N_PCA_GENES <- 1000L          ## PCA 用 top 可变基因，避免全矩阵协方差爆炸
## Q6 判据：块因子在单个主成分上的最大 R2 >= 0.5 记为「主导」。
Q6_R2_CUT <- 0.5
##   该阈值是本步的**操作定义**，非文献取值 —— 语义来自
##   方案 §9 Q6「块因子不主导 PCA」，「主导」需一个可判定边界，
##   取「解释一半以上方差」为界。未取到支持该具体数值的文献，
##   如实标注为约定值，不得在正文声称有文献出处。

library(sva)  ## 批次与主成分协变量估计（full model 口径，供 Q5/Q6 判定）
## 批次效应方差拆分（与 sva 并列，判是否需要 SV 校正）
library(variancePartition)
library(RUVSeq)  ## RLE 归一化旁证（若可用则作为第二意见证）
library(ggplot2)  ## 诊断图；主题与出口用规范 §sci_theme / §sci_save

## 主题：规范 §sci_theme 官方定义（theme_bw + 关网格），全脚本统一口径
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

## 出口：规范 §sci_save 官方定义，仅 PDF，落本步 FIG_DIR。
## 可读性核验钩子（规则 22 · 交付前须核验可读性）：仅当 RAPA_QC_PNG=1 时额外导 PNG
## 到 FIG_DIR/_qc/ 供核验；默认不产出，交付物仍为纯 PDF（符合 §sci_save）。
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

## PART 1 逐个已归一化系列做诊断 ---------------------------------------------------------
## 待诊断系列从 S03 的落盘目录列举：本步自己的 table 目录里没有表达矩阵。
## 2026-09-29 起改读**基因级**矩阵 `<acc>_gene.rds`（行 = 基因符号、列 = 真实样本 ID）。
## 理由：探针级 `<acc>_expr.rds` 的行名是 p1..pN、列名是 s1..sN，双轴身份均被替换，
## 下游按基因名/样本 ID 的连接全部落空（实测见 01_data/verify/check_s03_rownames_contract.R）。
vec.acc <- gsub("_gene\\.rds$", "",
  list.files(file.path("03_results", "S03", "table"),
    pattern = "_gene\\.rds$", full.names = FALSE))
if (length(vec.acc) == 0L) {
  stop(paste0("S04 FAILED: no gene-level expression matrix found under",
    "03_results/S03/table/ . "),
    "Run S03 first; S04 will not invent a series.")
}
vec.acc <- vec.acc[order(vec.acc)]

SUBJECT_KEY <- c(              ## 系列 -> 主体标识键（逐字取自 S01 实测键名）
  ## gse232911_grouping.tsv 的 subject 列（89 主体 / 246 样本）
  GSE232911 = "subject",
  GSE21545  = "individual",    ## design_keys.tsv
  GSE43292  = "patient",       ## design_keys.tsv
  GSE47472  = "subjects",      ## design_keys.tsv
  GSE57691  = "subjects",      ## design_keys.tsv
  GSE323983 = "donor"          ## sample_axis.tsv 的 donor 列
)
FILE.GROUPING <- "03_results/S01/table/gse232911_grouping.tsv"
FILE.AXIS <- "03_results/S01/table/sample_axis.tsv"
FILE.KEYS <- "03_results/S01/table/design_keys.tsv"
for (fp in c(FILE.GROUPING, FILE.AXIS, FILE.KEYS)) {
  if (!file.exists(fp)) stop(paste0("S04: 需要 ", fp,
      " 才能判定 Q5/Q6（请先跑 S01）。"))
}
## 循环内复用，避免重复 IO
dat.grouping <- read.delim(FILE.GROUPING, stringsAsFactors = FALSE)
dat.axis <- read.delim(FILE.AXIS, stringsAsFactors = FALSE)
dat.keys <- read.delim(FILE.KEYS, stringsAsFactors = FALSE)

## 取某系列的：主体标识（配对）、主体内因子（= 区组，供 Q6）、PCA 着色因子。
##   故无主体登记的系列仍用「第一个可用的设计键」着色（如 GSE100927 的 tissue、
##   GSE21545 的 sample type），而 **Q6 只认主体内因子**，两个语义不混。
BuildPairInfo <- function(acc, vec.sample) {
  out <- data.frame(sample_id = vec.sample, subject = NA_character_,
    block = NA_character_, color = NA_character_,
    stringsAsFactors = FALSE)
  setret <- function(blk, blk.n, col, col.n) {
    out$block <- blk
    out$color <- col
    list(dat = out, block.name = blk.n, color.name = col.n)
  }
  if (acc == "GSE232911") {
    i <- match(out$sample_id, dat.grouping$sample_id)
    out$subject <- as.character(dat.grouping$subject[i])
    v <- as.character(dat.grouping$tunica[i])
    return(setret(v, "tunica", v, "tunica"))
  }
  if (acc == "GSE323983") {
    ax <- dat.axis[dat.axis$acc == acc, , drop = FALSE]
    i <- match(out$sample_id, ax$sample_id)
    out$subject <- as.character(ax$donor[i])
    v <- as.character(ax$arm[i])
    return(setret(v, "arm", v, "arm"))
  }
  sub <- dat.keys[dat.keys$acc == acc, c("sample_id", "key", "value"),
    drop = FALSE]
  if (nrow(sub) == 0L) {
    return(setret(NA_character_, "", NA_character_, ""))
  }
  wide <- reshape(sub, idvar = "sample_id", timevar = "key", direction = "wide")
  names(wide) <- sub("^value\\.", "", names(wide))
  i <- match(out$sample_id, wide$sample_id)
  vec.cand <- setdiff(names(wide), "sample_id")
  if (!acc %in% names(SUBJECT_KEY)) {
    ## 无主体登记：只用第一个非空设计键给 PCA 着色，block 保持 NA（Q6 判 N/A）
    for (k in vec.cand) {
      v <- as.character(wide[[k]][i])
      if (!all(is.na(v))) return(setret(NA_character_, "", v, k))
    }
    return(setret(NA_character_, "", NA_character_, ""))
  }
  key.subj <- SUBJECT_KEY[[acc]]
  out$subject <- as.character(wide[[key.subj]][i])
  v.fb <- NA_character_
  k.fb <- ""
  for (k in setdiff(vec.cand, key.subj)) {
    v <- as.character(wide[[k]][i])
    if (all(is.na(v))) next
    if (!nzchar(k.fb)) {                      ## 第一个可用键作 PCA 着色回退
      v.fb <- v
      k.fb <- k
    }
    tab <- tapply(v, out$subject, function(x) length(unique(x[!is.na(x)])))
    if (any(tab > 1L, na.rm = TRUE)) {
      return(setret(v, k, v, k))  ## 主体内因子：既是区组也是着色因子
    }
  }
  setret(NA_character_, "", v.fb, k.fb)
}

dat.diag <- data.frame(acc = character(0), n_sv = integer(0),
  pca_var_explained = character(0), q5_paired = character(0),
  q6_block = character(0), note = character(0),
  stringsAsFactors = FALSE)
list.sv <- list()

for (acc in vec.acc) {
  obj <- readRDS(file.path(TAB_DIR, "..", "..", "S03", "table",
      paste0(acc, "_gene.rds")))
  expr <- obj$expr
  cat(sprintf("\n[S04] %s   expr: %d symbols x %d samples\n", acc, nrow(expr),
      ncol(expr)))

  info.pair <- BuildPairInfo(acc, colnames(expr))
  dat.design <- data.frame(sample = colnames(expr), stringsAsFactors = FALSE)
  if (acc == "GSE232911") {
    ## SV 模型（下方 ~ patient + tunica + thrombus）需要 thrombus 列，故此处仍 merge 分组表
    dat.design <- merge(dat.design, dat.grouping, by.x = "sample",
      by.y = "sample_id", sort = FALSE)
    stopifnot(all(dat.design$sample %in% colnames(expr)))
    dat.design$patient <- dat.design$subject
  } else {
    dat.design <- cbind(dat.design,
      info.pair$dat[, c("subject", "block"), drop = FALSE])
  }
  ## 防丢样本：设计表必须覆盖表达矩阵的全部样本（否则 PCA 着色出现 NA 而不报错）
  stopifnot(all(colnames(expr) %in% dat.design$sample))
  dat.design$pca.color <- ifelse(is.na(info.pair$dat$color), "not registered",
    info.pair$dat$color)
  dat.design$pca.color.name <- if (nzchar(
      info.pair$color.name)) info.pair$color.name
  else "(none registered in S01)"

  ## top 可变基因 -> PCA
  ## NA sd 必须先剔除：双色系列（GSE6002）的基因级矩阵含 NA（limma 官方语义：
  ## MA.RG 源码 R/G <= 0 置 NA），avereps 归并后部分基因行非 NA 值 < 2 -> sd=NA。
  ## 旧写法 min(N_PCA_GENES, length(vec.sd)) 用**含 NA 的全长度**，
  ## sort 默认剔除 NA 后 names 不足 -> seq_len 越界 -> names 给出 NA -> 表达矩阵
  ## 下标越界（2026-09-29 实测 GSE6002：1042 行中 139 行 sd=NA，903 < 1000）。
  vec.sd <- apply(expr, 1L, stats::sd)
  vec.sd <- vec.sd[!is.na(vec.sd)]
  vec.gene <- names(sort(vec.sd, decreasing = TRUE))[seq_len(min(N_PCA_GENES,
        length(vec.sd)))]
  mat.pca <- prcomp(t(expr[vec.gene, , drop = FALSE]), scale. = FALSE)
  vec.var <- round((mat.pca$sdev^2 / sum(mat.pca$sdev^2))[1:3], 3)
  dat.pca <- as.data.frame(mat.pca$x[, 1:2])
  colnames(dat.pca) <- c("PC1", "PC2")
  dat.pca$pca.color <- as.character(dat.design$pca.color)

  dat.pca$acc <- acc
  ## [user 2026-09-29] 标题与图例长文本按词换行：GSE6002 的设计键值
  ##   "CD14+ isolated monocytes from buffy coat."（42 字符）令标题单行 57 字符
  ##   超出 6.5 in 画布宽，实测标题左端 "GSE6" 被裁、图例条右溢——是裁切不是
  ##   数据错误。strwrap 按词折行且对相同输入确定性输出，离散映射分组不变。
  txt.title <- paste(strwrap(sprintf("PCA of %s coloured by %s", acc,
        dat.design$pca.color.name[1L]), width = 48), collapse = "\n")
  txt.legend <- paste(strwrap(dat.design$pca.color.name[1L],
      width = 26), collapse = "\n")
  dat.pca$pca.color <- vapply(dat.pca$pca.color, function(s) paste(
      strwrap(s, width = 26), collapse = "\n"), character(1L))
  names(dat.pca$pca.color) <- NULL
  sci_save(ggplot(dat.pca, aes(x = PC1, y = PC2, colour = pca.color)) +
    geom_point(size = 1.2, alpha = 0.8) +
    sci_theme +
    labs(title = txt.title,
      x = sprintf("PC1 (%.1f%%)", 100 * vec.var[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * vec.var[2]),
      colour = txt.legend),
    sprintf("PCA_Coloured_By_Design_Factor_%s", acc),
    width = 6.5, height = 5.5)

  ## Q5：配对因子完整性 —— 每个主体是否都有 >= 2 个样本（配对设计的可估计前提）
  vec.subj <- info.pair$dat$subject
  lgl.subj <- !is.na(vec.subj) & nzchar(vec.subj)
  if (sum(lgl.subj) == 0L) {
    q5 <- "NOT_APPLICABLE"
    note5 <- "no subject identifier registered in S01"
  } else {
    tab.subj <- table(vec.subj[lgl.subj])
    n.ok5 <- sum(tab.subj >= 2L)
    q5 <- if (n.ok5 == length(tab.subj)) "PASS" else "FAIL"
    note5 <- sprintf(
      "%d subjects (%d/%d samples carry an id), %d with >=2 samples",
      length(tab.subj), sum(lgl.subj), length(vec.subj), n.ok5)
  }

  ## Q6：块因子是否主导 PCA（方案 §9 原文语义）。
  ##   >= Q6_R2_CUT 视为「主导」→ FAIL。
  if (!nzchar(info.pair$block.name) || all(is.na(info.pair$dat$block))) {
    q6 <- "NOT_APPLICABLE"
    note6 <- "no within-subject block factor found in S01"
  } else {
    vec.blk <- factor(info.pair$dat$block)
    r2.max <- max(vapply(seq_len(min(3L, ncol(mat.pca$x))), function(j) {
          summary(stats::lm(mat.pca$x[, j] ~ vec.blk))$r.squared
        }, numeric(1L)))
    q6 <- if (r2.max >= Q6_R2_CUT) "FAIL" else "PASS"
    note6 <- sprintf("%s: max R2 on PC1-PC3 = %.2f", info.pair$block.name,
      r2.max)
  }

  ## SV 估计（仅对已有设计因子的系列做；其余无法估）
  n.sv <- NA_integer_
  note.sv <- ""
  if (acc == "GSE232911" && length(unique(dat.design$patient)) > 1L) {
    fml.sv <- ~ patient + tunica + thrombus
    mat.model <- stats::model.matrix(fml.sv, data = dat.design)
    ## mod0 显式给出：sva:::irwsva.build 源码
    ##   `if (is.null(mod0)) { mod0 <- mod[, 1] }` —— 官方零模型就是「仅截距」。
    ##   显式传参把该语义固定下来，避免日后换基/换列序时 mod[, 1] 悄悄变义。
    mat.mod0 <- stats::model.matrix(~ 1, data = dat.design)
    obj.sv <- sva::sva(expr[vec.gene, , drop = FALSE], mod = mat.model,
      mod0 = mat.mod0, method = "two-step")
    note.sv <- paste0("SV design = ", paste(deparse(fml.sv), collapse = ""),
      " ; method = two-step (irw numerically singular at ncol(mod)=",
      ncol(mat.model), ") ; n_sv=", ncol(obj.sv$sv),
      " ; aligned to plan 5")
    n.sv <- ncol(obj.sv$sv)
    list.sv[[acc]] <- obj.sv
    dat.sv <- as.data.frame(obj.sv$sv)
    names(dat.sv) <- paste0("SV", seq_len(ncol(obj.sv$sv)))
    dat.sv.long <- data.frame(
      sv = rep(names(dat.sv), each = nrow(dat.sv)),
      value = as.numeric(as.matrix(dat.sv)))
    dat.sv.long$acc <- acc
    sci_save(ggplot(dat.sv.long, aes(x = sv, y = value, fill = sv)) +
      geom_boxplot(width = 0.6, outlier.size = 1.4,
        colour = "grey25") +
      geom_hline(yintercept = 0, linetype = "dashed",
        colour = "grey40") +
      ## 着色按 SV 区分（user 2026-09-29：原图无颜色不可判读）。
      ## 图例省略：fill 与 x 轴一一对应，重复图例只会占版面。
      scale_fill_manual(values = c(SV1 = "#4C72B0", SV2 = "#DD8452",
          SV3 = "#55A868", SV4 = "#C44E52",
          SV5 = "#8172B3", SV6 = "#937860"),
        guide = "none") +
      labs(title = sprintf("Surrogate variables estimated for %s", acc),
        subtitle = "Surrogate variables from sva (two-step method)",
        x = "Surrogate variable", y = "SV value") +
      sci_theme,
      sprintf("Surrogate_Variables_Boxplot_%s", acc),
      width = 6.5, height = 5)
  }
  cat(sprintf("  PCA var(1-3) = %s | SV = %s | Q5=%s Q6=%s\n",
      paste(vec.var, collapse = ", "), n.sv,
      format(q5), format(q6)))

  dat.diag <- rbind(dat.diag, data.frame(
      acc = acc, n_sv = if (is.na(n.sv)) -1L else as.integer(n.sv),
      pca_var_explained = paste(vec.var, collapse = ";"),
      q5_paired = q5, q6_block = q6,
      note = paste0(note5, " | ", note6, " | ",
        note.sv), stringsAsFactors = FALSE))
}

if (length(list.sv) > 0L) {
  saveRDS(list.sv, file.path(TAB_DIR, "sv_object.rds"))
}
write.csv(dat.diag, file.path(TAB_DIR, "confound_diag.csv"), row.names = FALSE)

## PART 1.5 跨系列诊断图 -----------------------------------------------------------------------
dat.diag.plot <- rbind(
  data.frame(acc = dat.diag$acc, check = "Q5 paired",
    result = dat.diag$q5_paired,
    stringsAsFactors = FALSE),
  data.frame(acc = dat.diag$acc, check = "Q6 block", result = dat.diag$q6_block,
    stringsAsFactors = FALSE))
dat.diag.plot$acc <- factor(dat.diag.plot$acc)
## 判定值的取值域显式断言：`factor(levels=)` 会把未登记值静默变成 NA（把「通过」画成空白格），
## 故先断言取值域再映射。当前取值域 = PASS / FAIL / NOT_APPLICABLE。
stopifnot(all(dat.diag.plot$result %in% c("PASS", "FAIL", "NOT_APPLICABLE")))
dat.diag.plot$result3 <- factor(dat.diag.plot$result,
  levels = c("PASS", "FAIL", "NOT_APPLICABLE"),
  labels = c("PASS", "FAIL", "N/A"))
sci_save(ggplot(dat.diag.plot, aes(x = acc, y = check, fill = result3)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c(PASS = "#4C72B0", FAIL = "#C44E52",
      "N/A" = "grey70"),
    name = "Diagnostic result") +
  labs(title = "Pairing and block-factor diagnosis per series",
    subtitle = paste0("Q5 = every subject has >= 2 samples; Q6 = block factor ",
      "explains < ", Q6_R2_CUT, " of PC1-PC3 variance\n",
      "N/A = no subject or block factor registered in S01 ",
      "(not a failure)"),
    x = "Series", y = "Diagnostic") +
  sci_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  "Confound_Diagnosis_By_Series", width = 10, height = 3.5)

## PART 2 收口 -----------------------------------------------------------------------------
vec.fail <- dat.diag$acc[dat.diag$q5_paired == "FAIL" | dat.diag$q6_block ==
  "FAIL"]
message(paste0("\n===== S04 confound summary =====\nseries diagnosed = ", nrow(
      dat.diag),
    "  series with >=1 FAIL = ", length(vec.fail)))
if (length(vec.fail) > 0L) {
  message(paste0("[FAIL] ", paste(vec.fail, collapse = ", "),
      paste0("\n  该判定为诊断结论，见 confound_diag.csv 的 q5_paired",
        "/ q6_block / note 列。")))
}
message("S04 OK. tables + figures -> 03_results/S04/")
