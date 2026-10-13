#===============================================================
## S09_deconv.R —— 免疫微环境去卷积与 M1/M2 极化（步骤规格 §8：S09）
##
## 目的：对组织 bulk 用多个去卷积工具估计免疫细胞比例，取多工具共识，
##       再与层 2 段分做**偏相关**（控制其余细胞类型），
##       检验「段活性是否主要由免疫细胞组成驱动」。
##
## 定位：分析脚本。三工具 = quantiseqr + EPIC + BayesPrism（2026-09-28
##       用户裁定 A/A/A/A：删 timerc；BayesPrism 走三段式；参考用
##       动脉粥样硬化斑块 atlas 聚合 GEP；共识取五大 lineage）。
## 产出：03_results/S09/table/deconv_plan.csv
##       03_results/S09/table/tool_availability.csv
##       03_results/S09/table/immune_consensus_<acc>.csv
##       03_results/S09/table/tool_celltype_fractions_<acc>.csv
##       03_results/S09/table/segment_fraction_corr.csv
##       03_results/S09/Figure/Deconvolution_Composition_StackedBar.pdf
##       03_results/S09/Figure/Deconvolution_Fraction_ByLineage.pdf
##       03_results/S09/Figure/Deconvolution_Fraction_vs_SegmentScore.pdf
##       03_results/S09/Figure/Deconvolution_PartialCorrelation_Heatmap.pdf
##
## 图件清单（规则 35 可视化强制 · SCI 发表导向；动笔前已调研文献确定风格/形式/数量）：
##   可视化: S09 共 4 类
##     1. Deconvolution_Composition_StackedBar.pdf
##        —— 形式：每样本堆叠柱（x = sample，y = 共识比例，fill = lineage，按系列分面）
##        —— 科学问题：各样本免疫细胞组成（去卷积主图）。
##        —— 参照文献：CIBERSORT (Newman et al., Nat Methods 2015) Fig 2A 堆叠柱；
##           Nat Sci Rep 2026 Fig 4A；quanTIseq 官方 vignette 明示 stacked bar。
##     2. Deconvolution_Fraction_ByLineage.pdf
##        —— 形式：分组箱线（x = lineage，y = 共识比例，fill = 系列）
##        —— 科学问题：五大 lineage 比例在三动脉系列间的分布比较。
##        —— 参照文献：Nat Sci Rep 2026 Fig 4B 分组箱线惯例。
##     3. Deconvolution_Fraction_vs_SegmentScore.pdf
##        —— 形式：分面散点网格（行 = segment，列 = 系列，color = lineage）
##        —— 科学问题：段分与免疫细胞比例在样本级的共变形态（覆盖全部系列）。
##        —— 参照文献：Avila Cobos 2020, Nat Commun, PMID 33077803。
##     4. Deconvolution_PartialCorrelation_Heatmap.pdf
##        —— 形式：热图（行 = segment，列 = lineage，fill = 偏相关 r，p<0.05 标 "*"）
##        —— 科学问题：控制其余 lineage 后，各段分与该 lineage 比例的偏相关方向量级。
##        —— 参照文献：CIBERSORT Fig 2B 热图惯例；TIMER2.0, PMID 32442275。
##
## 图件口径：sci_theme / sci_save 逐字取自规范 §sci_theme / §sci_save 官方定义；
##       分面一律 facet_wrap / facet_grid，不引入任何拼图包；可读性核验见
##       sci_save 之 RAPA_QC_PNG 钩子。
##
## 前置：S03 + S06；S09a_build_reference.R 产出的 BayesPrism 参考文件。
##       去卷积工具无中央基线表可查，一律在**调用点硬断言**形参与槽位。
#===============================================================

## PART 0 常量与工具注册表 --------------------------------------------------------------

TAB_DIR <- "03_results/S09/table"
FIG_DIR <- "03_results/S09/Figure"
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

library(ggplot2)  ## 去卷积图；主题与出口用规范 §sci_theme / §sci_save

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

## 出口：规范 §sci_save 官方定义，仅 PDF，落本步 FIG_DIR
sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  pdf(pdf_file, width = width, height = height)
  print(plot_obj)
  dev.off()
  ## 可读性核验钩子（规则 35）：默认关闭，交付物仍为纯 PDF；
  ## 置 RAPA_QC_PNG=1 时额外导 _qc/ PNG 供逐张目视核验，不改变 PDF 产出。
  if (Sys.getenv("RAPA_QC_PNG") == "1") {
    dir.create(file.path(FIG_DIR, "_qc"), showWarnings = FALSE)
    png(file.path(FIG_DIR, "_qc", paste0(filename, ".png")),
      width = width * 150, height = height * 150, res = 150)
    print(plot_obj)
    dev.off()
  }
}

MIN_SAMPLES <- 10L            ## 少于此样本数不做去卷积
MIN_FRACTION <- 0.001         ## 低于该值记 NA（避免把噪声当 0）
CONSENSUS_FN <- "median"      ## 多工具共识：跨工具按 lineage 取中位数
LOG2_MAX <- 50  ## 线性标度判据：表达矩阵最大值低于它视为 log2 尺度
## （log2 表达通常 <= 25，线性 TPM/强度 > 100）

vec.segment <- c("S1_initiation", "S2_amplification", "S3_effect",
  "S4_resolution", "S5_chronic")
vec.acc <- unique(c("GSE232911", "GSE57691", "GSE43292"))

## lineage 共识映射（2026-09-28 裁决 D6-A）：
##   依据 XDec-SM（PMC7049190）区室级 common denominator 与 TIMER2.0
##   （PMID 32442275）「全工具共有类型」策略：共识只取三工具都可估的
##   五大 lineage；亚型（M1/M2、Tregs）按区室口径加总；工具特异类型
##   （Monocytes/Neutrophils/DC、CAFs/Endothelial/SMC/Other）不进共识、
##   单独落盘报告。关键词在「小写 + [._-] 转空格」归一后的名称上做包含
##   匹配，对 BayesPrism 参考的任意命名（如 "B cells"/"CD4 T cells"）
##   同样成立，不需要预知 atlas 的字面类型名。
dat.lineage.kw <- data.frame(
  lineage = c("B.cells", "B.cells", "T.cells.CD4", "T.cells.CD4",
    "T.cells.CD8", "T.cells.CD4", "NK.cells", "NK.cells", "Macrophages"),
  keyword = c("b cell", "bcells", "cd4", "treg", "cd8", "t cell",
    "nk", "natural killer", "macrophage"),
  stringsAsFactors = FALSE)
list.tool <- list(
  list(tool = "quantiseqr", pkg = "quantiseqr",
    ref_file = NA_character_, orient = "CxS",
    steps = list(
      list(fn = "run_quantiseq", matrix_arg = "expression_data",
        chain_arg = NA_character_,
        ## 默认值逐项实测于 quantiseqr 1.20.0 formals（2026-09-28）：
        ## signature_matrix 默认即 "TIL10"；is_arraydata=FALSE 因输入
        ## 是 S03 归一化矩阵而非 quanTIseq 自家 RNA-seq 管道产物；
        ## is_tumordata=FALSE 因动脉组织非肿瘤。
        fixed_args = list(signature_matrix = "TIL10",
          is_arraydata = FALSE,
          is_tumordata = FALSE,
          scale_mRNA = TRUE, method = "lsei",
          column = "gene_symbol", return_se = FALSE),
        ref_args = NULL, out_slot = NA_character_,
        transform = "df_to_cxs"))),
  list(tool = "EPIC", pkg = "EPIC",
    ref_file = NA_character_, orient = "SxC",
    steps = list(
      ## EPIC 1.1.7 实测（2026-09-27）：形参是 bulk 不是 expr；
      ## 返回 list(mRNAProportions, cellFractions, fit.gof)，比例矩阵是
      ## cellFractions。reference 必须写 NULL 而不是 NA：EPIC 内部先判
      ## is.null(reference) 再判 is.character(reference)，传 NA 会被
      ## is.character(NA) 判真，走进"预置参考名比对"分支后 stop。
      list(fn = "EPIC", matrix_arg = "bulk",
        chain_arg = NA_character_,
        fixed_args = list(reference = NULL, scaleExprs = TRUE,
          withOtherCells = TRUE,
          constrainedSum = TRUE,
          rangeBasedOptim = FALSE),
        ref_args = NULL, out_slot = "cellFractions",
        transform = NA_character_))),
  list(tool = "BayesPrism", pkg = "BayesPrism",
    ## 参考由 S09a_build_reference.R 产出：斑块 atlas（Traeuble 2025,
    ## PMID 40931012，同疾病先例且其自身即用 BayesPrism 做 bulk 去卷积）
    ## 按细胞类型聚合的 type-by-gene GEP（summed counts）。
    ref_file = "01_data/references/BayesPrism_reference_atherosclerosis.rds",
    orient = "SxC",
    steps = list(
      ## 步骤 1/3：构造 prism 对象。实测 new.prism 形参（2026-09-27）：
      ## reference, input.type, cell.type.labels, cell.state.labels,
      ## key, mixture, outlier.cut, outlier.fraction, pseudo.min。
      ## mixture 官方要求 sample-by-gene（与 limma 习惯相反），故
      ## transform="t"；无恶性细胞 -> key=NULL。
      list(fn = "new.prism", matrix_arg = "mixture",
        chain_arg = NA_character_,
        fixed_args = list(input.type = "GEP", key = NULL,
          outlier.cut = 0.01, outlier.fraction = 0.1,
          pseudo.min = 1e-8),
        ref_args = c(reference = "gep",
          cell.type.labels = "cell.type.labels",
          cell.state.labels = "cell.state.labels"),
        out_slot = NA_character_, transform = "t"),
      ## 步骤 2/3：Gibbs 采样。n.cores=4 实测加速 2.99x 且不挂（探针
      ## probe_s09_parallel.R：4000 基因 12 样本，1核 306s -> 4核 102s；
      ## snowfall 干净 init/stop，无 PSOCK 死挂）。注：项目笔记的"并行慢 23-38x"
      ## 指 BiocParallel 的 SnowParam（另一后端），与 BayesPrism 的 snowfall 不同，
      ## 此处不适用。seed=123L 显式设定（项目规则，防 MCMC 不可复现；同 n.cores 下
      ## 可复现）。gibbs.control 轻量设置 chain.length=500, burn.in=100（thinning=2），
      ## 三系列 BayesPrism 实测约 0.7h（4核），总 ~1h。
      list(fn = "run.prism", matrix_arg = NA_character_,
        chain_arg = "prism",
        fixed_args = list(n.cores = 4L,
          gibbs.control = list(chain.length = 500L, burn.in = 100L,
            thinning = 2L)),
        ref_args = NULL, out_slot = NA_character_,
        transform = NA_character_, seed = 123L),
      ## 步骤 3/3：取细胞类型比例。实测 get.fraction(bp, which.theta=
      ## "final", state.or.type="type") 返回 样本 x 细胞类型 矩阵，
      ## 且与 bp@posterior.theta_f@theta 全等（2026-09-27 all.equal=TRUE）。
      list(fn = "get.fraction", matrix_arg = NA_character_,
        chain_arg = "bp",
        fixed_args = list(which.theta = "final",
          state.or.type = "type"),
        ref_args = NULL, out_slot = NA_character_,
        transform = NA_character_))))

## 登记表落盘为长表：每工具每步一行，含全部实参（deparse 成字符串），
## 便于复核与复现——登记表与真实调用之间不允许存在看不到的字段
dat.plan <- do.call(rbind, lapply(list.tool, function(t) {
      do.call(rbind, lapply(seq_along(t$steps), function(i) {
            s <- t$steps[[i]]
            data.frame(tool = t$tool, step = i, fn = s$fn,
              matrix_arg = as.character(s$matrix_arg),
              chain_arg = as.character(s$chain_arg),
              out_slot = as.character(s$out_slot),
              transform = as.character(s$transform),
              fixed_args = if (length(s$fixed_args) == 0L) "" else
              paste(vapply(names(s$fixed_args), function(nm)
                  paste0(nm, "=", paste(deparse(s$fixed_args[[nm]]),
                      collapse = " ")), ""),
                collapse = "; "),
              ref_args = if (is.null(s$ref_args)) "" else
              paste(names(s$ref_args), "->", s$ref_args, collapse = "; "),
              ref_file = as.character(t$ref_file),
              orient = t$orient, stringsAsFactors = FALSE)
          }))
    }))
## 必须排在 dat.plan 构造之后（旧写法在构造前写盘，报 object not found）
write.csv(dat.plan, file.path(TAB_DIR, "deconv_plan.csv"), row.names = FALSE)

## PART 1 工具可用性登记 + 外部参考存在性核验 ------------------------------------------
## 签名正确性不在本段判定，而是在每个工具的调用点做形参/槽位硬断言（PART 2）。
## 本段只回答两个问题：包装了没 / 参考文件在不在。

dat.avail <- data.frame(tool = character(0), pkg = character(0),
  status = character(0),
  note = character(0), stringsAsFactors = FALSE)
for (t in list.tool) {
  has.pkg <- requireNamespace(t$pkg, quietly = TRUE)
  ref.status <- if (is.na(t$ref_file)) "NA" else
  if (file.exists(t$ref_file)) "PRESENT" else "REFERENCE_MISSING"
  status <- if (!has.pkg) "PKG_NOT_INSTALLED" else
  if (ref.status == "REFERENCE_MISSING") "REFERENCE_MISSING" else "PRESENT"
  dat.avail <- rbind(dat.avail, data.frame(
      tool = t$tool, pkg = t$pkg, status = status,
      note = paste0("steps=", length(t$steps), "; ref=", ref.status),
      stringsAsFactors = FALSE))
}
write.csv(dat.avail, file.path(TAB_DIR, "tool_availability.csv"),
  row.names = FALSE)

vec.notready <- dat.avail$tool[dat.avail$status != "PRESENT"]
if (length(vec.notready) > 0L) {
  stop(paste0("S09 STOP: deconvolution tool(s) not usable: ",
      paste(vec.notready, collapse = ", "),
      " (status in 03_results/S09/table/tool_availability.csv). ",
      "Either install the package or supply its reference file."))
}

## 参考文件核验：参考缺失时不得开跑（去卷积结果完全由参考决定）
vec.ref.missing <- list.tool[vapply(list.tool, function(t)
    !is.na(t$ref_file) && !file.exists(t$ref_file), TRUE)]
if (length(vec.ref.missing) > 0L) {
  stop("S09 STOP: external reference(s) missing: ",
    paste(vapply(vec.ref.missing, function(t) t$ref_file, ""),
      collapse = "; "),
    ". Deconvolution results are entirely reference-determined; ",
    "run 02_scripts/S09a_build_reference.R first.")
}

## PART 2 共识比例与段分偏相关 -----------------------------------------------------------

## 名称归一 + lineage 关键词匹配；未命中返回 NA（单独报告，不进共识）
MapLineage <- function(vec.ct) {
  vec.norm <- tolower(gsub("[._-]", " ", vec.ct))
  vec.out <- rep(NA_character_, length(vec.ct))
  for (i in seq_len(nrow(dat.lineage.kw))) {
    hit <- is.na(vec.out) &
    grepl(dat.lineage.kw$keyword[i], vec.norm, fixed = TRUE)
    vec.out[hit] <- dat.lineage.kw$lineage[i]
  }
  vec.out
}

## quantiseqr 的返回变换：data.frame（Sample + 10 类 + Other）->
## CxS 矩阵（去 Sample 与 Other 两列后转置；Other 为 quanTIseq 未分类桶，
## 属工具特异类型，走单独报告通道）
FunctionDfToCxS <- function(df.f) {
  vec.drop <- c("Sample", "Other")
  vec.hit <- intersect(vec.drop, colnames(df.f))
  if (!"Sample" %in% vec.hit) {
    stop("S09 STOP: quantiseqr returned a data.frame without 'Sample' ",
      "column; actual colnames: ",
      paste(colnames(df.f), collapse = ", "), ".")
  }
  t(as.matrix(df.f[, setdiff(colnames(df.f), vec.drop), drop = FALSE]))
}

## BayesPrism 基因集过滤（2026-09-30，方案 B）：全 16,035 共享基因在 GSE232911
## 246 样本上 Gibbs 采样会编译级崩溃（无 R 级报错，sample 1 后即 EXIT=1）；
## 4000 基因基准实测可正常跑完且快（6 样本 122s）。故对 BayesPrism 仅保留 bulk
## 中高变基因子集（标准去卷积做法，低变基因主要贡献噪声），参考 GEP 同步按同一
## 基因集裁剪以保证 gene 名对齐。quantiseqr/EPIC 不受影响（自带签名/参考）。
GENE_FILTER_N <- 4000L

## 逐步执行一个工具的全部 steps，返回最终比例矩阵（原始朝向）
FunctionRunTool <- function(t, mat.expr) {
  obj.prev <- NULL        ## 上一步产物
  ref.obj <- NULL         ## 参考对象（按需载入一次）
  gn.keep <- NULL         ## BayesPrism 专用：本系列高变基因集
  for (s in t$steps) {
    fn.f <- get(s$fn, envir = asNamespace(t$pkg))
    ## 硬断言①：matrix_arg / chain_arg 必须真在该函数签名里
    vec.formal <- names(formals(fn.f))
    vec.need <- c(s$matrix_arg, s$chain_arg)
    vec.need <- vec.need[!is.na(vec.need)]
    vec.arg.miss <- setdiff(vec.need, vec.formal)
    if (length(vec.arg.miss) > 0L) {
      stop(paste0("S09 STOP: ", t$tool, " step ", s$fn,
          " lacks required arg(s): ",
          paste(vec.arg.miss, collapse = ", "),
          "; actual formals: ", paste(vec.formal, collapse = ", ")))
    }
    args.call <- s$fixed_args
    if (!is.null(s$ref_args)) {
      if (is.null(ref.obj)) {
        ref.obj <- readRDS(t$ref_file)
        if (t$tool == "BayesPrism") {
          vec.var <- apply(mat.expr, 1, var)
          ## 取 top-N 高变基因的**名称**（非位置）：order() 返回整数位置，
          ## 须用 rownames(mat.expr)[位置] 转回基因名，否则与参考 gep 的字符
          ## 列名求交为空 -> 空矩阵 -> validate.input 报错（与基准探针同陷阱）。
          gn.names <- rownames(mat.expr)[head(order(-vec.var), GENE_FILTER_N)]
          gn.keep <- intersect(gn.names, colnames(ref.obj$gep))
          if (length(gn.keep) == 0L) {
            stop("S09 STOP: BayesPrism gene filter yielded 0 genes; ",
              "check mat.expr rownames vs reference colnames alignment.")
          }
          ref.obj$gep <- ref.obj$gep[, gn.keep, drop = FALSE]
        }
      }
      vec.mem.miss <- setdiff(s$ref_args, names(ref.obj))
      if (length(vec.mem.miss) > 0L) {
        stop(paste0("S09 STOP: reference for ", t$tool, " lacks member(s): ",
            paste(vec.mem.miss, collapse = ", "),
            "; actual members: ",
            paste(names(ref.obj), collapse = ", ")))
      }
      for (nm in names(s$ref_args)) {
        args.call[[nm]] <- ref.obj[[s$ref_args[[nm]]]]
      }
    }
    ## 输入级转置（transform = "t"）：部分工具要求 mixture 为 sample x gene
    ##（如 BayesPrism new.prism，与 limma 习惯的 gene x sample 相反）。
    ## 旧写法把 transpose 误用到函数结果（prism S4 对象）上，导致 mixture
    ## 仍以 gene x sample 传入，colnames 为样本 GSM 号，与参考基因名交集为 0，
    ## 触发 "gene names of reference and mixture do not match"。此处改为转输入。
    mat.use <- mat.expr
    if (t$tool == "BayesPrism" && !is.null(gn.keep)) {
      mat.use <- mat.expr[gn.keep, , drop = FALSE]
    }
    if (!is.na(s$transform) && s$transform == "t") mat.use <- t(mat.use)
    if (!is.na(s$matrix_arg)) args.call[[s$matrix_arg]] <- mat.use
    if (!is.na(s$chain_arg)) args.call[[s$chain_arg]] <- obj.prev
    ## 可复现：若本步声明了 seed，调用前显式设定（MCMC/Gibbs 随机性）
    if (!is.null(s$seed) && !is.na(s$seed)) set.seed(s$seed)
    res <- do.call(fn.f, args.call)
    obj.prev <- if (!is.na(s$out_slot)) {
      ## 硬断言②：结果槽位必须真的存在（取值前验，不猜槽位）
      if (!is.list(res) || !s$out_slot %in% names(res)) {
        stop(paste0("S09 STOP: ", t$tool, " step ", s$fn,
          " return has no slot '", s$out_slot,
          "'; actual top-level names: ",
          if (is.list(res)) paste(names(res), collapse = ", ")
          else "<not a list>",
          ". Fix out_slot in list.tool before rerunning."))
      }
      res[[s$out_slot]]
    } else {
      res
    }
    ## 输出级变换：仅 df_to_cxs（quantiseqr data.frame -> CxS）适用；
    ## "t" 已在输入级处理，此处不再对结果转置，避免对 S4/非矩阵对象误操作。
    if (!is.na(s$transform) && s$transform == "df_to_cxs") {
      obj.prev <- FunctionDfToCxS(obj.prev)
    }
  }
  obj.prev
}

f.consensus <- switch(CONSENSUS_FN,
  median = stats::median,
  mean = mean,
  stop("S09: CONSENSUS_FN must be 'median' or 'mean'"))

dat.corr <- data.frame(acc = character(0), segment = character(0),
  fraction = character(0), n_celltype = integer(0),
  partial_r = numeric(0), t_stat = numeric(0),
  p = numeric(0), stringsAsFactors = FALSE)
dat.consensus.all <- list()

for (acc in vec.acc) {
  ## 2026-09-29 起改读**基因级**矩阵：解卷积工具（EPIC/quantiseq/… 与其签名矩阵）都以
  ## 基因符号为键，探针级 p1..pN 行名对不上任何签名；列名更须与 S06 的段分矩阵同名可比。
  file.expr <- file.path("03_results", "S03", "table", paste0(acc, "_gene.rds"))
  if (!file.exists(file.expr)) {
    stop(paste0("S09 STOP: gene-level expression matrix for ", acc,
        " missing (S03)."))
  }
  obj.expr <- readRDS(file.expr)
  mat.expr <- obj.expr$expr
  if (ncol(mat.expr) < MIN_SAMPLES) {
    stop(paste0("S09 STOP: ", acc, " has ", ncol(mat.expr),
        " samples (< ", MIN_SAMPLES, "); not eligible for deconvolution."))
  }
  ## 线性标度（2026-09-28 裁决 D2-A 附属）：Avila Cobos PMID 33077803
  ## 线性标度最优；BayesPrism 官方文档 log 不可用。三工具统一喂同一
  ## 线性矩阵，保证共识口径一致。
  if (max(mat.expr) < LOG2_MAX) {
    mat.expr <- 2^mat.expr
  }
  ## 段分（S06 产物：段 x 样本）；列须与表达矩阵列同名才可比
  file.score <- file.path("03_results", "S06", "table",
    "segment_score_matrix.rds")
  if (!file.exists(file.score)) {
    stop(paste0("S09 STOP: ", file.score, " missing (run S06 first)."))
  }
  ## S06 产物为 per-acc 列表，须按名取本系列段分矩阵（与 S08 list.score[[acc]] 同口径）
  mat.score <- readRDS(file.score)[[acc]]
  vec.col.miss <- setdiff(colnames(mat.expr), colnames(mat.score))
  if (length(vec.col.miss) > 0L) {
    stop(paste0("S09 STOP: segment score matrix lacks ", length(vec.col.miss),
        " sample column(s) present in ", acc, " expression matrix ",
        "(e.g. ", vec.col.miss[1L], "). Column alignment is mandatory."))
  }

  ## 逐个工具跑（steps 多段式），收集「细胞类型 x 样本」比例矩阵
  ## —— PART 级 checkpoint（规则 33）：本系列产物在盘则复用，避免重复昂贵 BayesPrism。
  file.frac <- file.path(TAB_DIR, paste0("tool_celltype_fractions_", acc,
    ".csv"))
  file.cons <- file.path(TAB_DIR, paste0("immune_consensus_", acc, ".csv"))
  reuse.mode <- if (file.exists(file.cons)) "consensus"
    else if (file.exists(file.frac)) "fractions" else "none"
  list.frac <- list()
  dat.ct.long <- data.frame(tool = character(0), tool_ct = character(0),
    lineage = character(0), sample = character(0),
    fraction = numeric(0), stringsAsFactors = FALSE)
  if (reuse.mode == "none") {
    for (t in list.tool) {
      mat.raw <- FunctionRunTool(t, mat.expr)
      ## 朝向归一化：本脚本下游一律按「细胞类型 x 样本」处理
      if (is.na(t$orient)) {
        stop(paste0("S09 STOP: ", t$tool, " has orient = NA in list.tool; ",
            "measure the returned matrix orientation and register it. ",
            "Guessing here would produce a silently wrong consensus."))
      }
      if (t$orient == "SxC") mat.f <- t(mat.raw)
      else if (t$orient == "CxS") mat.f <- mat.raw
      else {
        stop("S09 STOP: ", t$tool, " orient must be 'CxS' or 'SxC', got '",
          t$orient, "'.")
      }
      mat.f <- as.matrix(mat.f)
      rownames(mat.f) <- gsub("^ *| *$", "", rownames(mat.f))
      ## 朝向校验后才能改写列名：朝向一旦搞反，这步会把样本名写到
      ## 「细胞类型」行上且不报错。
      vec.col.miss <- setdiff(colnames(mat.expr), colnames(mat.f))
      if (length(vec.col.miss) > 0L) {
        stop(paste0("S09 STOP: ", t$tool, " fraction matrix lacks ",
            length(vec.col.miss), " sample column(s) of ", acc,
            " (e.g. ", vec.col.miss[1L], "); dim = ",
            paste(dim(mat.f), collapse = "x"),
            " after orient='", t$orient, "'. Register the real orientation."))
      }
      colnames(mat.f) <- colnames(mat.expr)
      ## lineage 映射：命中的按区室口径加总；未命中的全量进单独报告表
      vec.lin <- MapLineage(rownames(mat.f))
      dat.ct.long <- rbind(dat.ct.long, data.frame(
          tool = t$tool,
          tool_ct = rep(rownames(mat.f), times = ncol(mat.f)),
          lineage = rep(vec.lin, times = ncol(mat.f)),
          sample = rep(colnames(mat.f), each = nrow(mat.f)),
          fraction = as.numeric(mat.f),
          stringsAsFactors = FALSE))
      vec.lin.ok <- unique(vec.lin[!is.na(vec.lin)])
      if (length(vec.lin.ok) == 0L) {
        stop(paste0("S09 STOP: ", t$tool, " mapped zero cell type(s) to the ",
            "five shared lineages; rownames were: ",
            paste(rownames(mat.f), collapse = ", "), "."))
      }
      list.frac[[t$tool]] <- rowsum(mat.f[!is.na(vec.lin), , drop = FALSE],
        group = vec.lin[!is.na(vec.lin)])
    }
    ## 工具特异类型全量落盘（未映射 + 已映射的原始类型都在表里，可追溯）
    write.csv(dat.ct.long,
      file.path(TAB_DIR, paste0("tool_celltype_fractions_", acc, ".csv")),
      row.names = FALSE)
  } else if (reuse.mode == "fractions") {
    ## 仅工具级产物在（偏相关/共识前崩溃）：从 fractions 重建 list.frac，跳过 BayesPrism。
    ## 重建 = 同 lineage 的 cell type 比例按样本求和（与 rowsum 等价），再按样本名对齐。
    dat.ct.long <- read.csv(file.frac)
    for (tl in unique(dat.ct.long$tool)) {
      sub <- dat.ct.long[dat.ct.long$tool == tl, ]
      m <- as.matrix(xtabs(fraction ~ lineage + sample, data = sub))
      m <- m[, colnames(mat.expr),
        drop = FALSE]   # 对齐样本顺序（xtabs 默认字母序）
      list.frac[[tl]] <- m
    }
    cat(sprintf(paste0("\n[S09] %-11s fractions checkpoint reused (BayesPrism",
      "skipped)\n"), acc))
  } else {
    ## reuse.mode == "consensus"：共识已落盘，下方直接载入，跳过工具循环与共识计算
    cat(sprintf("\n[S09] %-11s consensus checkpoint reused (full skip)\n", acc))
  }

  if (reuse.mode == "consensus") {
    mat.cons <- as.matrix(read.csv(file.cons, row.names = 1))
    stopifnot(all(colnames(mat.cons) == colnames(mat.expr)))
    vec.lin.shared <- rownames(mat.cons)
  } else {
    ## 跨工具共识：仅保留**全部工具都报出**的 lineage（2026-09-28 裁决 D6-A）
    vec.lin.shared <- Reduce(intersect, lapply(list.frac, rownames))
    if (length(vec.lin.shared) == 0L) {
      stop(paste0("S09 STOP: no lineage reported by all tools for ", acc,
          "; shared lineages are required before taking a consensus."))
    }
    ## 堆叠成 (lineage, 样本, 工具) 三维数组后用 base 取共识，不引入额外依赖
    mat.stack <- array(
      unlist(lapply(list.frac, function(m) m[vec.lin.shared, , drop = FALSE])),
      dim = c(length(vec.lin.shared), ncol(mat.expr), length(list.frac)),
      dimnames = list(vec.lin.shared, colnames(mat.expr), names(list.frac)))
    mat.cons <- apply(mat.stack, c(1L, 2L), f.consensus)
  }
  mat.cons[mat.cons < MIN_FRACTION] <- NA_real_
  stopifnot(all(rownames(mat.cons) == vec.lin.shared))
  dat.consensus.all[[acc]] <- mat.cons
  ## 共识落盘前置到偏相关之前：即便偏相关崩溃也不丢共识（规则 33 续传）
  write.csv(cbind(lineage = rownames(mat.cons), as.data.frame(mat.cons)),
    file.path(TAB_DIR, paste0("immune_consensus_", acc, ".csv")),
    row.names = FALSE)

  ## 偏相关：控制其余 lineage 后，比例与段分的相关
  vec.seg.used <- intersect(vec.segment, rownames(mat.score))
  for (seg in vec.seg.used) {
    v.seg <- as.numeric(mat.score[seg, colnames(mat.expr)])
    for (lin in vec.lin.shared) {
      ## 统一完整病例掩码：段分、本 lineage 共识、其余 lineage 协变量任一为 NA
      ##（含 MIN_FRACTION 阈值化产生的 NA）即剔除，保证 lm 残差长度一致，
      ## 否则 cor.test 报 'x' and 'y' must have the same length。
      v.frac <- as.numeric(mat.cons[lin, colnames(mat.expr)])
      mat.other <- mat.cons[setdiff(vec.lin.shared, lin), colnames(mat.expr),
        drop = FALSE]
      keep <- !is.na(v.seg) & !is.na(v.frac) &
        !apply(is.na(mat.other), 2L, any)
      if (sum(keep) < 3L) {
        stop("S09 STOP: ", acc, " / ", seg, " / ", lin,
          " usable samples < 3 after NA removal (MIN_FRACTION threshold); ",
          "partial correlation infeasible.")
      }
      y <- v.seg[keep]
      v.frac.k <- v.frac[keep]
      mat.cov <- cbind(1, t(as.matrix(mat.other[, colnames(mat.expr)[keep],
              drop = FALSE])))
      r1 <- residuals(lm(y ~ mat.cov))
      r2 <- residuals(lm(v.frac.k ~ mat.cov))
      r <- cor.test(r1, r2)
      if (is.na(r$p.value)) {
        stop("S09 STOP: ", acc, " / ", seg,
          " 的偏相关不可估（残差零方差），按规则 26 即时报错而非抑制。")
      }
      dat.corr <- rbind(dat.corr, data.frame(
          acc = acc, segment = seg, fraction = lin,
          n_celltype = length(vec.lin.shared),
          partial_r = as.numeric(r$estimate), t_stat = as.numeric(r$statistic),
          p = as.numeric(r$p.value), stringsAsFactors = FALSE))
    }
  }
  cat(sprintf("\n[S09] %-11s lineages = %d  segments = %d\n", acc,
      length(vec.lin.shared), length(vec.seg.used)))
}

write.csv(dat.corr, file.path(TAB_DIR, "segment_fraction_corr.csv"),
  row.names = FALSE)

## PART 3 图件 ---------------------------------------------------------------------------
## 命名遵守规范 16（描述性英文，禁 Fig 前缀）；图内文字全英文（规范 17）；
## 图内一律用半角连字符与冒号，不出破折号；分面一律 facet_wrap / facet_grid。
dat.frac.long <- do.call(rbind, lapply(names(dat.consensus.all), function(a) {
      m <- as.matrix(dat.consensus.all[[a]])
      data.frame(acc = a,
        lineage = rep(rownames(m), times = ncol(m)),
        sample = rep(colnames(m), each = nrow(m)),
        fraction = as.numeric(m), stringsAsFactors = FALSE)
    }))
## 系列间 lineage 轴一致（共识只取五共有 lineage），free_x 属冗余分面自由度，去掉。
## 谱系配色：4 共有 lineage 固定离散调色板（色盲友好），全图一致（规则 28）。
LIN_COLOR <- c("B.cells" = "#4477AA", "Macrophages" = "#EE6677",
  "NK.cells" = "#228833", "T.cells.CD4" = "#CCBB44")
## 系列配色：随 consensus.all 实际名称动态分配，无硬编码（规则 28）
ACC_COLOR <- setNames(c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#882255"),
  names(dat.consensus.all))

## 图 1：每样本免疫细胞组成堆叠柱（去卷积主图；CIBERSORT Fig 2A / Nat Sci Rep
## 2026 Fig 4A / quanTIseq vignette 明示 stacked bar）。每系列一面板，每样本一根堆叠柱。
sci_save(ggplot(dat.frac.long, aes(x = sample, y = fraction, fill = lineage)) +
  geom_col(position = "stack", width = 1) +
  facet_wrap(~ acc, scales = "free_x") +
  scale_fill_manual(values = LIN_COLOR) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Immune cell composition per sample",
    x = "Sample", y = "Fraction", fill = "Lineage") +
  sci_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
  "Deconvolution_Composition_StackedBar", width = 12, height = 5)

## 图 2：各谱系比例按系列分组比较（箱线；Nat Sci Rep 2026 Fig 4B 标准配套）
sci_save(ggplot(dat.frac.long, aes(x = lineage, y = fraction, fill = acc)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 1) +
  scale_fill_manual(values = ACC_COLOR) +
  labs(title = "Consensus lineage fractions by series",
    x = "Lineage", y = "Fraction", fill = "Series") +
  sci_theme,
  "Deconvolution_Fraction_ByLineage", width = 9, height = 5)

## 图 3：谱系比例 vs 段分（散点），遍历全部系列。
## 三元组 (lineage, segment, sample) 显式对齐。
dat.pt <- do.call(rbind, lapply(names(dat.consensus.all), function(a) {
      mat.c <- as.matrix(dat.consensus.all[[a]])
      dat.s <- readRDS(file.path("03_results", "S06", "table",
          "segment_score_matrix.rds"))[[a]]
      do.call(rbind, lapply(rownames(mat.c), function(lin) {
            do.call(rbind, lapply(rownames(dat.s), function(seg) {
                  data.frame(acc = a, lineage = lin, segment = seg,
                    fraction = as.numeric(mat.c[lin, ]),
                    score = as.numeric(dat.s[seg, colnames(mat.c)]),
                    stringsAsFactors = FALSE)
                }))
          }))
    }))
dat.pt$segment <- factor(dat.pt$segment, levels = vec.segment)
sci_save(ggplot(dat.pt, aes(x = fraction, y = score, color = lineage)) +
  geom_point(size = 1.0, alpha = 0.6) +
  facet_grid(segment ~ acc) +
  scale_color_manual(values = LIN_COLOR) +
  labs(title = "Lineage fraction versus segment score",
    x = "Consensus fraction", y = "Segment score", color = "Lineage") +
  sci_theme,
  "Deconvolution_Fraction_vs_SegmentScore", width = 14, height = 8)

## 图 4：段分-谱系偏相关热图（替代原裸点图；CIBERSORT Fig 2B 热图惯例）。
## 双色渐变以 0 为中心，p < 0.05 标 "*"。
dat.hm <- dat.corr[!is.na(dat.corr$partial_r), ]
dat.hm$lineage_name <- dat.hm$fraction
dat.hm$segment <- factor(dat.hm$segment, levels = vec.segment)
sci_save(ggplot(dat.hm, aes(x = lineage_name, y = segment, fill = partial_r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_wrap(~ acc) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, limits = c(-1, 1)) +
  geom_text(aes(label = ifelse(p < 0.05, "*", "")), size = 6, color = "black") +
  labs(title = "Partial correlation of lineage fraction with segment score",
    x = "Lineage", y = "Segment", fill = "Partial r") +
  sci_theme +
  ## 其 axis.text；hjust=1 顶端对齐刻度、vjust=0.5 水平居中，
  ## 长名（Macrophages / T.cells.CD4）垂直展开不截断。
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)),
  "Deconvolution_PartialCorrelation_Heatmap", width = 11, height = 5.5)

## PART 4 收口 ----------------------------------------------------------------------------

n.null <- sum(is.na(dat.corr$p))
message(paste0("\n===== S09 deconvolution summary =====\ntool registry rows = ",
    nrow(dat.avail), "  ready = ",
    sum(dat.avail$status == "PRESENT"),
    "\nseries processed = ", length(dat.consensus.all),
    "\ncorrelation rows = ", nrow(dat.corr),
    "  NA p = ", n.null))
if (n.null > 0L) {
  message("[NOTE] ", n.null, " partial correlation(s) yielded NA p-value ",
    "(e.g. zero variance after residualisation); they are reported as NA, ",
    "NOT imputed as non-significant.")
}
message("S09 OK. tables + figure -> 03_results/S09/")
