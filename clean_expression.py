import pandas as pd

df = pd.read_csv("results/transcriptome/S2/expression_normalized.tsv", sep="\t")

df = df.drop(columns=["index"])
df = df.set_index("gene")

df.to_csv("results/transcriptome/S2/expression_for_wgcna.tsv", sep="\t")

print("File cleaned ✅")