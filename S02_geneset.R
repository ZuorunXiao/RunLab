#===============================================================
## S02_geneset.R —— 五段基因集锁定（步骤规格 §8：S02）
## 目的：把已锁定的炎症基因集（MSigDB 派生的精确 gs_name 清单）落到「段 / 亚段 /
##       来源 / 基因集名」四列结构，用 msigdbr 做**全量精确成员检验**，
##       并产出段间重叠矩阵供 S13 披露重叠基因。
## 定位：分析脚本。只读基因集表，不碰表达矩阵。
## 产出：03_results/S02/table/inflammation_geneset_final.csv  段/亚段/来源/gs_name/成员
# （Q3 处置：仅落「富集可用」集 member_count >= MIN_GS_SIZE，下游 S06/S07/S12 直接消费）
# 03_results/S02/table/segment_member_count.csv         全量审计：每段/每集成员数与 pass_min_sz 判定
# 03_results/S02/table/subthreshold_genesets_excluded.csv  Q3 显式登记被排除的 <MIN_GS_SIZE
# 集（留痕）
##       03_results/S02/table/segment_overlap.csv              段间 Jaccard
##       03_results/S02/table/segment_overlap_genes.csv        跨段共享基因
# 03_results/S02/Figure/Geneset_MemberCount_Distribution.pdf
# 五段成员数分布（小提琴+箱线+散点，min.sz=20 参考）
##       03_results/S02/Figure/Geneset_Composition_By_Segment_Source.pdf  单面板堆叠柱：段 x 来源
# 03_results/S02/Figure/Segment_GeneOverlap_UpSet.pdf         段间基因共享 UpSet 图（Conway 2017）
## 图件清单（规则 35 可视化强制 · SCI 发表导向；动笔前已调研已发表文献确定风格/形式/数量）：
##   可视化: S02 共 3 幅
##     1. Geneset_MemberCount_Distribution.pdf
##        —— 形式：小提琴图 + 箱线图 + 抖动散点（纵轴 log10，min.sz=20 红色虚线参考，亚阈值集标灰）
##        —— 回答的科学问题：五段锁定基因集的成员数分布是否足以支撑下游富集（≥20），
##           并直观暴露 Q3 显式排除的亚阈值集（灰色点）。
##        —— 参照文献：基因集规模分布是 MSigDB 富集前标准 QC（Subramanian 2005
##           PMID 16199517；Liberzon 2011 PMID 21773919 之 C2/C5 集合规模惯例）；
##           散点叠加箱线为 RNA-seq 基因集 QC 通用画法。
##     2. Geneset_Composition_By_Segment_Source.pdf
##        —— 形式：单面板堆叠柱状图（x = 段，fill = 来源）
##        —— 回答的科学问题：每个炎症段的锁定基因集分别来自哪些 MSigDB 子库来源，揭示段内来源构成。
##        —— 参照文献：来源构成堆叠柱为基因集注释溯源标准画法（Liberzon 2015
##           PMID 25971742 之 MSigDB v5 集合来源拆解）。
##     3. Segment_GeneOverlap_UpSet.pdf
##        —— 形式：UpSet 图（集合交集矩阵，五段各作为一个集合）
##        —— 回答的科学问题：五段在基因层面是独立信号还是高度重叠（§4.2 重叠折扣的依据），
##           并量化跨段共享基因规模。
##        —— 参照文献：Conway JR et al. Bioinformatics 33(18):2938-2940, 2017
# (doi:10.1093/bioinformatics/btx364, PMID 28645171) —— UpSetR，集合交集的可发表标准画法；
##           基础方法 Lex A et al. IEEE TVCG 20(12):1983-1992, 2014
##           (doi:10.1109/TVCG.2014.2346248)。
##           （运行环境若未装 UpSetR，自动降级为 5×5 段间 Jaccard 热图，同样可发表。）
##
## 图件口径：sci_theme / sci_save 逐字取自规范 §sci_theme / §sci_save 官方定义；
##       不引入任何拼图包，单面板 ggplot2 出图；可读性核验见 sci_save 之 RAPA_QC_PNG 钩子。
## 前置：S01 通过；msigdbr 已安装（未装则本脚本直接停，不降级、不臆造成员）。
## 口径说明（填补方案 §8 S02 的一处留白）：
##   「subsegment」在四列结构中没有现成来源，本脚本按下列确定性规则生成并登记：
##   subsegment = 该 gs_name 所属规律标签 pattern_label（来自
##   inflammation_geneset_lock_v2.csv / _supplement.csv 的 pattern_label 列）；
##   若某 gs_name 未被规律标签覆盖，则取 gs_collection。该规则只做标签化，
##   不增删任何基因集，故不影响 §10 预注册的成员集合。
#===============================================================

## PART 0 常量与目录 ---------------------------------------------------------------

## 表格写出：文件名一律 .csv，内容用逗号分隔（名实相符；下游 read.csv 直接消费）。
## 禁 write.delim（base/utils 无此函数，勿凭记忆调用）。
WriteCsv <- function(dat, path) {
  write.csv(dat, path, quote = FALSE, row.names = FALSE, na = "NA")
}

dir.tab <- "03_results/S02/table"
FIG_DIR <- "03_results/S02/Figure"
dir.create(dir.tab, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

library(ggplot2)  ## 基因集构成图；主题与出口用规范 §sci_theme / §sci_save

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
## 可读性核验钩子（规则 22 三级审核 · 交付前须核验可读性）：仅当显式设置环境变量
## RAPA_QC_PNG=1 时，额外导出一份 PNG 预览到 FIG_DIR/_qc/ 供人工/AI 核验；
## 该 PNG 不进入交付物、默认不产出，保持交付物为纯 PDF（符合规范 §sci_save）。
sci_save <- function(plot_obj, filename, width = 8, height = 6) {
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  pdf(pdf_file, width = width, height = height)
  print(plot_obj)
  dev.off()
  if (Sys.getenv("RAPA_QC_PNG") == "1") {
    dir.create(file.path(FIG_DIR, "_qc"), showWarnings = FALSE)
    png_file <- file.path(FIG_DIR, "_qc", paste0(filename, ".png"))
    png(png_file, width = width * 150, height = height * 150, res = 150)
    print(plot_obj)
    dev.off()
  }
}

MIN_GS_SIZE <- 20L  ## Q3：每段可用基因 >= 20（与 §8 S02 / S06 min.sz 一致）
## §4.2：两段 Jaccard > 0.5 视为同一信号源，其「联合贡献」x 0.5 折扣。
## 折扣的施加位置是 S13（在段 x 数据集 x 三层方向 的汇总处），本脚本只负责
## 在 segment_overlap.csv 里给出 jaccard 与 overlap_discounted 判定，
## 保证「判定证据在测量处、折扣施加在汇总处」两处可分别复核。
JACCARD_MAX <- 0.5
vec.acc.gset <- c("01_data/verify/inflammation_geneset_lock_v2.csv",
  "01_data/verify/inflammation_geneset_supplement.csv")
file.cand <- "01_data/verify/inflammation_geneset_lock_candidates.csv"
file.check <- "01_data/verify/inflammation_geneset_membership_check.csv"

## 五段（§1 / §4.2 命名；与 lock_candidates$segment 的取值一致）
vec.segment <- c("S1_initiation", "S2_amplification", "S3_effect",
  "S4_resolution", "S5_chronic")

## PART 1 读入并合并锁定表 -----------------------------------------------------------

dat.gset <- do.call(rbind, lapply(vec.acc.gset, function(p) {
  transform(read.csv(p, stringsAsFactors = FALSE), source_file = basename(p))
    }))

## 段归属：只认 lock_candidates 中显式标定的 segment；**禁止推断**
dat.cand <- read.csv(file.cand, stringsAsFactors = FALSE)
dat.seg <- unique(dat.cand[, c("gs_name", "segment")])
dat.check <- read.csv(file.check, stringsAsFactors = FALSE)

dat.final <- dat.gset[, c("pattern_label", "collection_tag", "gs_name")]
names(dat.final) <- c("subsegment", "source", "gs_name")
dat.final$subsegment[is.na(
    dat.final$subsegment) | dat.final$subsegment == ""] <-
dat.final$source[is.na(dat.final$subsegment) | dat.final$subsegment == ""]
dat.final <- merge(dat.final, dat.seg, by = "gs_name", all.x = TRUE)
dat.final$segment <- ifelse(is.na(dat.final$segment), "UNMAPPED",
  as.character(dat.final$segment))
dat.final <- unique(dat.final[, c("segment", "subsegment", "source",
      "gs_name")])
dat.final <- dat.final[order(dat.final$segment, dat.final$subsegment,
    dat.final$gs_name),
]

## Q3 / §8：全部 gs_name 必须有显式段归属，不得出现 UNMAPPED
dat.unmapped <- dat.final[dat.final$segment == "UNMAPPED", ]
if (nrow(dat.unmapped) > 0L) {
  stop(paste0("S02 FAILED: ", nrow(dat.unmapped),
      " geneset(s) have no explicit segment ",
      "assignment, e.g. ", paste(head(dat.unmapped$gs_name, 3L),
        collapse = ", "),
      " ; §10 预注册禁止在结果产生后再指派段。"))
}
vec.bad.seg <- setdiff(unique(dat.final$segment), vec.segment)
if (length(vec.bad.seg) > 0L) {
  stop(paste0("S02 FAILED: segment value(s) outside the five-segment scheme: ",
      paste(vec.bad.seg, collapse = ", ")))
}

## PART 2 全量精确成员检验（msigdbr，规则 23：不得凭记忆写基因集名） --------------------

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  stop(paste0("S02 FAILED: package 'msigdbr' not installed ; cannot",
    "fetch exact gene-set members. "),
    "Install it first (Bioconductor) — the script will not guess membership.")
}
dat.msig <- msigdbr::msigdbr(species = "Homo sapiens")
vec.msig.on <- unique(dat.msig$gs_name)

dat.final$in_msigdb <- dat.final$gs_name %in% vec.msig.on

## 交叉核对：与 01_data/verify 中 S01 阶段已完成的 MSigDB 成员核验表比对。
## 若那里记为 FOUND 而本次 msigdbr 匹配不到，说明 MSigDB 版本漂移或基因集名改写；
## 二者混用会让「成员清单」不可复现，故直接 stop，不做静默降级。
vec.check.hit <- intersect(dat.check$candidate[dat.check$verdict == "FOUND"],
  dat.final$gs_name)
dat.false.neg <- setdiff(vec.check.hit, dat.final$gs_name[dat.final$in_msigdb])
if (length(dat.false.neg) > 0L) {
  stop(paste0("S02 FAILED: ", length(dat.false.neg),
      " gene-set(s) recorded FOUND in ",
      basename(file.check), " but not matched in the current MSigDB release: ",
      paste(head(dat.false.neg, 5L), collapse = ", "),
      " . MSigDB version drift — §10 pre-registration requires an exact, ",
      "reproducible member list; resolve before proceeding."))
}
message(paste0("[S02] cross-check vs ", basename(file.check), ": ", length(
      vec.check.hit),
    " FOUND gene-set(s), mismatches = ", length(dat.false.neg)))

dat.final$members <- NA_character_
dat.final$member_count <- NA_integer_

vec.all.member <- list()
for (gs in unique(dat.final$gs_name)) {
  idx <- which(dat.msig$gs_name == gs)
  vec.member <- unique(as.character(dat.msig$gene_symbol[idx]))
  dat.final$members[dat.final$gs_name == gs] <- if (length(idx) > 0L) paste(
    vec.member,
    collapse = ";") else ""
  dat.final$member_count[dat.final$gs_name == gs] <- length(vec.member)
  vec.all.member[[gs]] <- vec.member
}

dat.size <- unique(dat.final[, c("segment", "gs_name", "member_count")])
dat.size$pass_min_sz <- dat.size$member_count >= MIN_GS_SIZE
n.fail.size <- sum(!dat.size$pass_min_sz)

## PART 3 段间重叠 ---------------------------------------------------------------------

## 每段基因并集：以 dat.final 的段归属为准（PART 1 已确保无 UNMAPPED）。
vec.seg.gene <- lapply(vec.segment, function(s)
  sort(unique(unlist(vec.all.member[dat.final$gs_name[
            dat.final$segment == s]]))))

## 五段 x 五段全部有序对均落行（含 n_shared=0 的空对），保证段覆盖可核验；
## overlap_discounted 供 S13 施加 §4.2 的 x 0.5 折扣。
dat.overlap <- data.frame(row.segment = character(0), col.segment = character(
    0),
  n_genes_a = integer(0), n_genes_b = integer(0),
  jaccard = numeric(0), n_shared = integer(0),
  overlap_discounted = logical(0),
  shared_genes = character(0), stringsAsFactors = FALSE)
for (i in seq_along(vec.segment)) {
  for (j in seq_along(vec.segment)) {
    vec.a <- vec.seg.gene[[i]]
    vec.b <- vec.seg.gene[[j]]
    vec.shared <- intersect(vec.a, vec.b)
    num <- length(vec.shared)
    den <- length(union(vec.a, vec.b))
    dat.overlap <- rbind(dat.overlap, data.frame(
        row.segment = vec.segment[i], col.segment = vec.segment[j],
        n_genes_a = length(vec.a), n_genes_b = length(vec.b),
        jaccard = if (den > 0L) num / den else NA_real_,
        n_shared = num,
        overlap_discounted = if (den > 0L) (num / den) > JACCARD_MAX else FALSE,
        shared_genes = paste(vec.shared,
          collapse = ";"), stringsAsFactors = FALSE))
  }
}
dat.ov.hit <- dat.overlap[dat.overlap$n_shared > 0L, ]
dat.overlap.gene <- unique(paste(dat.ov.hit$shared_genes,
    dat.ov.hit$row.segment,
    dat.ov.hit$col.segment, sep = " || "))

## PART 4 写出 ---------------------------------------------------------------------------
## Q3 处置（用户 2026-09-29 决策：显式登记并排除，下游仅用 >= MIN_GS_SIZE 的集）：
##   - inflammation_geneset_final.csv 只落「富集可用」集（member_count >= MIN_GS_SIZE）。
##     下游 S06/S07/S12 直接读此表，无需各自再过滤，从根上避免亚阈值集污染 cameraPR/roast。
##   - segment_member_count.csv 保留全量审计（含 pass_min_sz 判定列）。
##   - subthreshold_genesets_excluded.csv 显式登记被排除的 <MIN_GS_SIZE 集（gs_name/segment/
##     member_count/reason），留痕以满足协议「explicitly recorded」本意，下游不消费。

dat.final$pass_min_sz <- dat.final$member_count >= MIN_GS_SIZE
dat.final.pass <- dat.final[dat.final$pass_min_sz,
  c("segment", "subsegment", "source", "gs_name",
    "in_msigdb", "members", "member_count", "pass_min_sz")]
WriteCsv(dat.final.pass, file.path(dir.tab, "inflammation_geneset_final.csv"))
WriteCsv(dat.size, file.path(dir.tab, "segment_member_count.csv"))

## 亚阈值集显式登记（Q3 排除留痕；不进下游富集）
dat.excl <- dat.size[!dat.size$pass_min_sz, c("gs_name", "segment",
    "member_count")]
dat.excl$reason <- paste0("below Q3 min.sz=", MIN_GS_SIZE,
  " (", n.fail.size, " sets excluded in total)")
WriteCsv(dat.excl, file.path(dir.tab, "subthreshold_genesets_excluded.csv"))

WriteCsv(dat.overlap, file.path(dir.tab, "segment_overlap.csv"))
WriteCsv(data.frame(shared_gene_string = dat.overlap.gene,
    stringsAsFactors = FALSE),
  file.path(dir.tab, "segment_overlap_genes.csv"))

## PART 4.5 图件 ---------------------------------------------------------------------------
## 命名遵守规范（描述性英文，禁 Fig 前缀；图内文字全英文）；不引入拼图包，单面板出图。
## 三幅图分别对应本步三个科学问题（见文件头图件清单 / 规则 35）。

dat.size$segment <- factor(dat.size$segment, levels = vec.segment)
dat.final$segment <- factor(dat.final$segment, levels = vec.segment)

## 图 1：五段成员数分布（小提琴 + 箱线 + 抖动散点，log10 纵轴，min.sz=20 参考线）。
## 直接替代原「1522 个基因集名作 x 轴 + free_x」的不可读柱图；
## 灰点 = Q3 排除的亚阈值集（<20），蓝点 = 富集可用集，红色虚线标注阈值。
p.mem <- ggplot(dat.size, aes(x = segment, y = member_count)) +
geom_violin(fill = "#4C72B0", alpha = 0.22, scale = "width", color = NA) +
geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA,
  color = "black") +
geom_jitter(aes(color = pass_min_sz), width = 0.18, size = 0.7, alpha = 0.5) +
scale_y_log10() +
geom_hline(yintercept = MIN_GS_SIZE, linetype = "dashed", color = "#C44E52") +
scale_color_manual(values = c("FALSE" = "#999999", "TRUE" = "#4C72B0"),
  name = "pass min.sz") +
labs(title = "Gene-set member-count distribution across the five segments",
  subtitle = paste0("Red dashed line: min.sz = ", MIN_GS_SIZE,
    " (Q3 threshold); grey points = excluded sub-threshold sets"),
  x = "Inflammation segment", y = "Members matched in MSigDB (log10 scale)") +
sci_theme
sci_save(p.mem, "Geneset_MemberCount_Distribution", width = 9, height = 5.5)

## 图 2：段内来源构成（单面板堆叠柱，x = 段，fill = 来源）。
## 替代原 facet_grid(segment ~ source) 产生 25 个多为空面板的不可读分面。
dat.comp <- aggregate(gs_name ~ segment + source, dat.final,
  function(x) length(unique(x)))
names(dat.comp)[3L] <- "n_genesets"
## [user 2026-09-29] C8/H 是 msigdbr collection 代号（取值域实测
## msigdbr_collections(): C1-C9 + H），单独出现不可读，图例改全名；
## 其余来源代号（C2:CP:* / C5:GO:BP / C7:IMMUNESIGDB）本身自明，不动。
## 只改 dat.comp（图 2 专用聚合表），dat.final 的 source 原值保留。
vec.lab.src <- c("H" = "H: Hallmark", "C8" = "C8: Cell type signatures")
dat.comp$source <- ifelse(is.na(vec.lab.src[dat.comp$source]),
  dat.comp$source, vec.lab.src[dat.comp$source])
p.src <- ggplot(dat.comp, aes(x = segment, y = n_genesets, fill = source)) +
geom_col() +
labs(title = "Locked gene sets per segment by source",
  x = "Inflammation segment", y = "Number of gene sets", fill = "Source") +
sci_theme
sci_save(p.src, "Geneset_Composition_By_Segment_Source", width = 10,
  height = 5.5)

## 图 3：段间基因共享 UpSet 图（5 段各作为集合；Conway 2017, PMID 28645171）。
## 替代原 5x5 稀疏 Jaccard 热图：UpSet 直接量化「每段独有基因」与「跨段共享基因」规模，
## 是集合交集的可发表标准画法。运行环境未装 UpSetR 时降级为 5x5 Jaccard 热图。
if (requireNamespace("UpSetR", quietly = TRUE)) {
  vec.upset.list <- setNames(lapply(vec.seg.gene, function(g) unique(
        as.character(g))),
    vec.segment)
  p.ov <- UpSetR::upset(
    UpSetR::fromList(vec.upset.list),
    order.by = "freq", nsets = length(vec.segment),
    number.angles = 30, point.size = 3, line.size = 1.2,
    main.bar.color = "#4C72B0", sets.bar.color = "#55A868",
    text.scale = c(1.3, 1.1, 1.1, 1.0, 1.2, 1.0))
  sci_save(p.ov, "Segment_GeneOverlap_UpSet", width = 9, height = 6)
} else {
  message("[S02] UpSetR 未安装，段间重叠降级为 5x5 Jaccard 热图。")
  dat.ov.plot <- dat.overlap
  dat.ov.plot$row.segment <- factor(dat.ov.plot$row.segment,
    levels = vec.segment)
  dat.ov.plot$col.segment <- factor(dat.ov.plot$col.segment,
    levels = vec.segment)
  dat.ov.plot$jaccard[is.na(dat.ov.plot$jaccard)] <- 0
  p.ov <- ggplot(dat.ov.plot, aes(x = col.segment, y = row.segment,
      fill = jaccard)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f\n(n=%d)", jaccard, n_shared)),
    size = 3) +
  scale_fill_distiller(palette = "Blues", direction = 1, name = "Jaccard") +
  labs(title = "Pairwise Jaccard overlap between the five segments",
    x = "Segment", y = "Segment") +
  sci_theme
  sci_save(p.ov, "Segment_GeneOverlap_UpSet", width = 7, height = 6.5)
}

## PART 5 收口 ---------------------------------------------------------------------------

vec.seg.empty <- vec.segment[lengths(vec.seg.gene) == 0L]
if (length(vec.seg.empty) > 0L) {
  message(paste0("[NOTE] segment(s) with 0 genes after MSigDB matching: ",
      paste(vec.seg.empty, collapse = ", "),
      " -> they contribute no geneset to S05/S06/S07 and must be disclosed."))
}
message(paste0("\n===== S02 geneset summary =====\ngenesets locked = ", nrow(
      dat.final),
    "  in MSigDB = ", sum(dat.final$in_msigdb), "/", nrow(dat.final),
    "\nsegments = ", length(vec.segment),
    "  genesets below min.sz=", MIN_GS_SIZE, " = ", n.fail.size,
    "  -> enrichment-eligible (written) = ", nrow(dat.final.pass),
    "\npairs flagged overlap_discounted (Jaccard > ", JACCARD_MAX, ") = ",
    sum(dat.overlap$overlap_discounted)))
dat.size.by.seg <- aggregate(member_count ~ segment, dat.size, function(x)
  c(n = length(x), min = min(x), median = median(x)))
print(dat.size.by.seg)
## Q3 闸门：原 stop 改为「显式登记后继续」（用户 2026-09-29 决策）。
## 326 个 <MIN_GS_SIZE 的集已写入 subthreshold_genesets_excluded.csv，
## 仅作审计留痕，不参与下游富集；下游读 inflammation_geneset_final.csv（已过滤）。
if (n.fail.size > 0L) {
  message(paste0("[Q3] ", n.fail.size, " gene-set(s) with < ", MIN_GS_SIZE,
      " genes explicitly recorded in subthreshold_genesets_excluded.csv",
      " and excluded from downstream enrichment (decision 2026-09-29: ",
      "explicitly recorded, not re-quantified). Enrichment consumes the ",
      nrow(dat.final.pass), " remaining sets (min.sz=", MIN_GS_SIZE, ")."))
}
if (sum(!dat.final$in_msigdb) > 0L) {
  message(paste0("[NOTE] ", sum(!dat.final$in_msigdb),
      " gene-set(s) not found in current MSigDB release; their member list is ",
      "empty and they will be excluded downstream by min.sz=", MIN_GS_SIZE,
      "."))
}
message("S02 OK. tables -> 03_results/S02/table/")
