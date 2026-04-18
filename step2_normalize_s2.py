import pandas as pd
import numpy as np
import os

# create output folder if needed
os.makedirs("results/transcriptome/S2", exist_ok=True)

# load raw expression matrix
expr = pd.read_csv(
    "data/S2/expression_raw.csv",
    low_memory=False
)

# rename first column to gene
expr = expr.rename(columns={expr.columns[0]: "gene"})

# set gene column as index
expr = expr.set_index("gene")

# log transform (log2(x+1))
expr = np.log2(expr.clip(lower=0) + 1)

# z-score normalize each gene across samples
expr_norm = expr.sub(expr.mean(axis=1), axis=0)
expr_norm = expr_norm.div(expr.std(axis=1), axis=0)

# reset index
expr_norm = expr_norm.copy()
expr_norm = expr_norm.reset_index()
expr_norm = expr_norm.reset_index()

# save normalized matrix
expr_norm.to_csv(
    "results/transcriptome/S2/expression_normalized.tsv",
    sep="\t",
    index=False
)

print("S2 normalization complete.")