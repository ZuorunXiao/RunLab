
gene_order_path <- file.path(
  "X:/CRLM/scRNA/XX/01_QC/table",
  "gene_order.txt")
if (file.exists(gene_order_path)) file.remove(gene_order_path)

# 重新进入 PART 7 gene_order 逻辑
# 直接执行 rtracklayer 解析（验证本地 GTF 是否正常）
gtf_file <- "X:/CRLM/scRNA/XX/gencode.v44.primary_assembly.annotation.gtf.gz"
cat("解析 GTF...\n")
gtf <- rtracklayer::import(gtf_file, format = "gtf")
cat("GTF 加载成功，gene 数:", sum(gtf$type == "gene"), "\n")