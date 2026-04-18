import pandas as pd
import numpy as np
import os

os.makedirs("results/networks", exist_ok=True)

# load normalized expression
expr = pd.read_csv(
    "results/transcriptome/S1/expression_normalized.tsv",
    sep="\t"
)

# gene names
genes = expr["gene"]

# expression values
expr_values = expr.drop(columns=["gene"])

# compute correlation matrix
corr = expr_values.T.corr(method="pearson")

edges = []

threshold = 0.7

for i in range(len(genes)):
    for j in range(i+1, len(genes)):

        weight = corr.iloc[i, j]

        if abs(weight) >= threshold:

            edges.append({
                "nodeA": genes.iloc[i],
                "nodeB": genes.iloc[j],
                "weight": abs(weight),
                "source": "WGCNA"
            })

edges_df = pd.DataFrame(edges)

edges_df.to_csv(
    "results/networks/coexp_S1.tsv",
    sep="\t",
    index=False
)

print("Coexpression network created.")
print("Edges:", len(edges_df))