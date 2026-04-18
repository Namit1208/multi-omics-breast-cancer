# Load adjacency matrix
load("results/networks/adjacency.RData")

print("Adjacency loaded")

genes <- colnames(adj)

edges <- data.frame()
count <- 0

print("Extracting edges...")

for (i in 1:(ncol(adj)-1)) {
  for (j in (i+1):ncol(adj)) {
    
    val <- adj[i, j]
    
    # ✅ NA-safe + LOWER threshold (important fix)
    if (!is.na(val) && val > 0.2) {
      edges <- rbind(edges,
                     data.frame(
                       gene1 = genes[i],
                       gene2 = genes[j],
                       weight = val
                     ))
      count <- count + 1
    }
  }
}

print(paste("Total edges:", count))

# Save edge list
write.table(edges,
            file = "results/networks/coexp_wgcna_AD.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

print("Network saved as coexp_wgcna_AD.tsv ✅")