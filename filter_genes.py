import pandas as pd

df = pd.read_csv("results/transcriptome/S2/de.tsv", sep="\t")

# Keep only valid gene-like names (uppercase letters/numbers)
df = df[df['gene'].str.match(r'^[A-Z0-9]+$')]

# Remove antisense / weird ones again just in case
df = df[~df['gene'].str.contains(r'-AS|^MT-|^5S_|^7SK|^7SL')]

# Drop NA
df = df[df['gene'].notna()]

df.to_csv("results/transcriptome/S2/de_filtered.tsv", sep="\t", index=False)

print("Strict filtering done ✅")