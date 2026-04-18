import pandas as pd
import numpy as np

# load expression matrix
df = pd.read_csv(
    "results/transcriptome/S1/expression_raw.tsv",
    sep="\t"
)

# first column should be gene names
df = df.set_index(df.columns[0])

# convert everything to numeric
df = df.apply(pd.to_numeric, errors="coerce")

# log transform
df = np.log2(df + 1)

# z-score normalization per gene
df = df.sub(df.mean(axis=1), axis=0)
df = df.div(df.std(axis=1), axis=0)

# remove genes with NaN
df = df.dropna()

# save normalized matrix
df.to_csv(
    "results/transcriptome/S1/expression_normalized.tsv",
    sep="\t"
)

print("Normalization complete")
print("Shape:", df.shape)