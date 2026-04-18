library(limma)

# ----------------------------
# LOAD DATA
# ----------------------------
setwd("C:/analysis/Transcriptome layer")

expr <- read.table(
  "results/transcriptome/S2/expression_normalized.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

meta <- read.table(
  "results/transcriptome/S2/metadata.tsv",
  header = TRUE,
  sep = "\t"
)

# ----------------------------
# DESIGN MATRIX
# ----------------------------
group <- factor(meta$status)
design <- model.matrix(~ group)

# ----------------------------
# FIT MODEL
# ----------------------------
fit <- lmFit(expr, design)
fit <- eBayes(fit)

# ----------------------------
# RESULTS
# ----------------------------
res <- topTable(fit, coef=2, number=Inf)

res$gene <- rownames(res)

res <- res[, c("gene","logFC","P.Value","adj.P.Val")]

write.table(
  res,
  "results/transcriptome/S2/de_limma.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)