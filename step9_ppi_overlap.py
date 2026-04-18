import pandas as pd

# Load
tx = pd.read_csv("results/networks/coexp_wgcna_AD.tsv", sep="\t")
ppi = pd.read_csv("results/networks/string_ppi.tsv", sep="\t")

# Normalize gene names
for df in [tx, ppi]:
    df['gene1'] = df['gene1'].str.upper().str.strip()
    df['gene2'] = df['gene2'].str.upper().str.strip()

#  KEY STEP — make edges undirected
tx['edge'] = tx.apply(lambda x: tuple(sorted([x['gene1'], x['gene2']])), axis=1)
ppi['edge'] = ppi.apply(lambda x: tuple(sorted([x['gene1'], x['gene2']])), axis=1)

# Keep only edge + weight
tx = tx[['edge', 'weight']]
ppi = ppi[['edge']]

# Merge (intersection)
merged = pd.merge(tx, ppi, on='edge')

# Recover gene1, gene2
merged[['gene1', 'gene2']] = pd.DataFrame(merged['edge'].tolist(), index=merged.index)

# Save
merged[['gene1', 'gene2', 'weight']].to_csv(
    "results/networks/intersection_tx_ppi_AD.tsv",
    sep="\t",
    index=False
)

print("Intersection edges:", len(merged))