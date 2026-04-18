# ==============================
# WGCNA FULL PIPELINE (ALL GENES)
# ==============================

library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ------------------------------
# STEP 1 — Load data
# ------------------------------
expr <- read.table(
  "results/transcriptome/S2/expression_for_wgcna.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

# transpose → WGCNA needs samples as rows
datExpr <- t(expr)

# ------------------------------
# STEP 2 — Check good genes/samples
# ------------------------------
gsg <- goodSamplesGenes(datExpr, verbose = 3)

if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# ------------------------------
# STEP 3 — Sample clustering (QC)
# ------------------------------
sampleTree <- hclust(dist(datExpr), method = "average")

pdf("results/networks/sample_clustering.pdf")
plot(sampleTree, main = "Sample clustering")
dev.off()

# ------------------------------
# STEP 4 — Soft-threshold selection
# ------------------------------
powers <- c(1:20)

sft <- pickSoftThreshold(datExpr,
                         powerVector = powers,
                         verbose = 5)

pdf("results/networks/soft_threshold.pdf")

par(mfrow = c(1,2))

# scale-free topology
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit",
     type = "n")

text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers,
     col = "red")

abline(h = 0.9, col = "blue")

# mean connectivity
plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n")

text(sft$fitIndices[,1],
     sft$fitIndices[,5],
     labels = powers,
     col = "red")

dev.off()

# ------------------------------
# STEP 5 — Choose power
# ------------------------------
softPower <- 6   # adjust based on plot (usually 5–8)

# ------------------------------
# STEP 6 — Adjacency matrix
# ------------------------------
adjacency <- adjacency(datExpr,
                       power = softPower,
                       type = "signed")

save(adjacency,
     file = "results/networks/adjacency_full.RData")

# ------------------------------
# STEP 7 — TOM matrix
# ------------------------------
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# ------------------------------
# STEP 8 — Gene clustering
# ------------------------------
geneTree <- hclust(as.dist(dissTOM), method = "average")

pdf("results/networks/gene_clustering.pdf")
plot(geneTree, main = "Gene clustering (TOM)")
dev.off()

# ------------------------------
# STEP 9 — Module detection
# ------------------------------
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = 30
)

moduleColors <- labels2colors(dynamicMods)

pdf("results/networks/modules_colors.pdf")
plotDendroAndColors(geneTree,
                    moduleColors,
                    "Modules",
                    dendroLabels = FALSE)
dev.off()

# ------------------------------
# STEP 10 — Module eigengenes
# ------------------------------
MEs <- moduleEigengenes(datExpr,
                        colors = moduleColors)$eigengenes

# ------------------------------
# STEP 11 — Save module assignments
# ------------------------------
modules <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors
)

write.table(modules,
            "results/modules/modules_wgcna_full.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

# ------------------------------
# STEP 12 — Export coexpression edges
# ------------------------------
# threshold for edge export
threshold <- 0.2

edges <- which(adjacency > threshold, arr.ind = TRUE)

edge_list <- data.frame(
  Gene1 = colnames(datExpr)[edges[,1]],
  Gene2 = colnames(datExpr)[edges[,2]],
  Weight = adjacency[edges]
)

# remove self loops
edge_list <- edge_list[edge_list$Gene1 != edge_list$Gene2, ]

write.table(edge_list,
            "results/networks/coexp_wgcna_full.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)