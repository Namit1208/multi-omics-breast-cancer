import pandas as pd
from scipy.stats import ttest_ind

expr = pd.read_csv(
    "results/transcriptome/S1/expression_normalized.tsv",
    sep="\t",
    index_col=0
)

meta = pd.read_csv(
    "results/transcriptome/S1/metadata.tsv",
    sep="\t"
)

cases = meta[meta.status=="case"]["sample"]
controls = meta[meta.status=="control"]["sample"]

results = []

for gene in expr.index:

    case_vals = expr.loc[gene, cases]
    ctrl_vals = expr.loc[gene, controls]

    stat,p = ttest_ind(case_vals, ctrl_vals)

    logfc = case_vals.mean() - ctrl_vals.mean()

    results.append([gene,logfc,p])

df = pd.DataFrame(results, columns=["gene","logFC","pvalue"])

df.to_csv(
    "results/transcriptome/S1/de.tsv",
    sep="\t",
    index=False
)