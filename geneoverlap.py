import pandas as pd

tx = pd.read_csv("results/networks/coexp_wgcna_AD.tsv", sep="\t")
ppi = pd.read_csv("results/networks/disease_ppi.tsv", sep="\t")

tx_genes = set(tx['gene1']).union(set(tx['gene2']))
ppi_genes = set(ppi['gene1']).union(set(ppi['gene2']))

overlap = tx_genes.intersection(ppi_genes)

print("TX genes:", len(tx_genes))
print("PPI genes:", len(ppi_genes))
print("Overlap:", len(overlap))