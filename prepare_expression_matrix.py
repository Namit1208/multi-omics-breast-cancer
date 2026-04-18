import pandas as pd
import os

os.makedirs("results/transcriptome/S1", exist_ok=True)

# load GEO expression matrix
df = pd.read_csv(
    "data/GSE81538_gene_expression_405_transformed.csv"
)

# first column is gene symbol
df = df.rename(columns={df.columns[0]: "gene"})

# save as pipeline format
df.to_csv(
    "results/transcriptome/S1/expression_raw.tsv",
    sep="\t",
    index=False
)

print("Expression matrix prepared.")