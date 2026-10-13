# S08b_module_function.R
# 功能：对 S08 最强模块 ME2（528 基因，GSE232911）做 GO BP / Hallmark 功能注释
# 定位：为手稿 "r=0.915 resolution stage reflects biology" 提供可追溯的生物学事实
#       产物：03_results/S08/table/me2_go_enrich.csv / me2_hallmark_enrich.csv
# 依赖：clusterProfiler, org.Hs.eg.db, msigdbr（本机已装，本地注释无需联网）
# 数据：03_results/S08/table/wgcna_modules_GSE232911.csv
suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(dplyr)
})
set.seed(20261009)

MOD_CSV  <- "X:/XX/03_results/S08/table/wgcna_modules_GSE232911.csv"
GO_OUT   <- "X:/XX/03_results/S08/table/me2_go_enrich.csv"
HK_OUT   <- "X:/XX/03_results/S08/table/me2_hallmark_enrich.csv"

mod <- read.csv(MOD_CSV, stringsAsFactors = FALSE)
me2      <- mod$gene[mod$module == "ME2"]
universe <- unique(mod$gene)                       # 进入 WGCNA 的全部基因 = 背景

id <- bitr(me2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
univ_id <- bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

# GO Biological Process
go <- enrichGO(gene = id$ENTREZID, universe = univ_id, OrgDb = org.Hs.eg.db,
               ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
go_df <- as.data.frame(go) %>% arrange(p.adjust) %>% head(15)
write.csv(go_df, GO_OUT, row.names = FALSE)

# Hallmark（生物学状态概览）
h <- msigdbr(species = "Homo sapiens", collection = "H")
h_df <- h[, c("gs_name", "gene_symbol")]
colnames(h_df) <- c("term", "gene")
hk <- enricher(gene = me2, universe = universe, TERM2GENE = h_df,
               pAdjustMethod = "BH", pvalueCutoff = 0.05)
hk_df <- if (!is.null(hk)) as.data.frame(hk) %>% arrange(p.adjust) %>% head(10) else data.frame()

cat("ME2 genes mapped:", nrow(id), "/", length(me2), "\n")
cat("=== GO BP top 10 ===\n")
print(go_df %>% select(Description, Count, p.adjust) %>% head(10))
if (nrow(hk_df) > 0) {
  write.csv(hk_df, HK_OUT, row.names = FALSE)
  cat("=== Hallmark top 8 ===\n")
  print(hk_df %>% select(Description, Count, p.adjust) %>% head(8))
}
