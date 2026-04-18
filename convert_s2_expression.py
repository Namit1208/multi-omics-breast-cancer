import pandas as pd

expr = pd.read_csv(
    "data/S2/expression_raw.csv"
)

expr.rename(columns={expr.columns[0]: "gene"}, inplace=True)

expr.to_csv(
    "results/transcriptome/S2/expression_raw.tsv",
    sep="\t",
    index=False
)

print("S2 expression file fixed.")
print("Shape:", expr.shape)