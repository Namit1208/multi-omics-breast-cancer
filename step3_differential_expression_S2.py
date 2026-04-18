import pandas as pd
import scipy.stats as stats
import os

os.makedirs("results/transcriptome/S2", exist_ok=True)

# load normalized expression
expr = pd.read_csv(
    "results/transcriptome/S2/expression_normalized.tsv",
    sep="\t"
)

# load metadata
meta = pd.read_csv(
    "results/transcriptome/S2/metadata.tsv",
    sep="\t"
)

# identify case/control samples
case_samples = meta[meta["status"] == "case"]["sample"].tolist()
control_samples = meta[meta["status"] == "control"]["sample"].tolist()

results = []

for _, row in expr.iterrows():

    gene = row["gene"]

    case_vals = pd.to_numeric(row[case_samples], errors="coerce").values
    control_vals = pd.to_numeric(row[control_samples], errors="coerce").values

    if len(case_vals) == 0 or len(control_vals) == 0:
        continue

    logfc = case_vals.mean() - control_vals.mean()

    tstat, pval = stats.ttest_ind(case_vals, control_vals)

    # log fold change
    logfc = case_vals.mean() - control_vals.mean()

    # t-test
    tstat, pval = stats.ttest_ind(case_vals, control_vals)

    results.append({
        "gene": gene,
        "logFC": logfc,
        "pvalue": pval
    })

de_df = pd.DataFrame(results)

de_df.to_csv(
    "results/transcriptome/S2/de.tsv",
    sep="\t",
    index=False
)

print("S2 differential expression completed.")