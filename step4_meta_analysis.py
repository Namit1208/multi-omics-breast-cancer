import pandas as pd
import numpy as np
import os
from scipy.stats import combine_pvalues

os.makedirs("results/transcriptome", exist_ok=True)

# Load DE results
s1 = pd.read_csv("results/transcriptome/S1/de.tsv", sep="\t")
s2 = pd.read_csv("results/transcriptome/S2/de.tsv", sep="\t")

# Rename columns
s1 = s1.rename(columns={
    "logFC": "logFC_S1",
    "pvalue": "pvalue_S1"
})

s2 = s2.rename(columns={
    "logFC": "logFC_S2",
    "pvalue": "pvalue_S2"
})

# Merge by gene
merged = pd.merge(s1, s2, on="gene")

results = []

for _, row in merged.iterrows():

    logfc_values = [row["logFC_S1"], row["logFC_S2"]]
    pvals = [row["pvalue_S1"], row["pvalue_S2"]]

    # meta logFC = mean
    meta_logfc = np.mean(logfc_values)

    # combine p-values (Fisher method)
    meta_pval = combine_pvalues(pvals)[1]

    direction = "up" if meta_logfc > 0 else "down"

    results.append({
        "gene": row["gene"],
        "meta_logFC": meta_logfc,
        "meta_pvalue": meta_pval,
        "direction": direction
    })

meta_df = pd.DataFrame(results)

meta_df.to_csv(
    "results/transcriptome/meta_DEG_BC.tsv",
    sep="\t",
    index=False
)

print("Meta-analysis completed.")
print("Output: results/transcriptome/meta_DEG_BC.tsv")