import pandas as pd

expr = pd.read_csv(
    "results/transcriptome/S2/expression_raw.tsv",
    sep="\t"
)

print(expr.shape)
print(expr.columns[:20])