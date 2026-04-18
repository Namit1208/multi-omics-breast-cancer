import pandas as pd

ppi = pd.read_csv("data/ppi/string_ppi.tsv", sep="\t")

# Keep strong interactions only
ppi = ppi[ppi['combined_score'] > 700]

# Rename columns
ppi = ppi.rename(columns={
    'protein1': 'gene1',
    'protein2': 'gene2'
})

# Keep only required columns
ppi = ppi[['gene1', 'gene2']]

# Normalize gene names
ppi['gene1'] = ppi['gene1'].str.upper().str.strip()
ppi['gene2'] = ppi['gene2'].str.upper().str.strip()

# Save
ppi.to_csv("results/networks/string_ppi.tsv", sep="\t", index=False)

print("STRING PPI ready ✅")