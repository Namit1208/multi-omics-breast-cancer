import pandas as pd
import requests
import os

os.makedirs("results/networks", exist_ok=True)

# Load DEG genes
deg = pd.read_csv("results/transcriptome/meta_DEG_BC.tsv", sep="\t")
genes = deg['gene'].dropna().unique().tolist()

print("Total genes:", len(genes))

# STRING API
string_api_url = "https://string-db.org/api"
output_format = "tsv"
method = "network"

params = {
    "identifiers": "%0d".join(genes[:2000]),  # limit (STRING constraint)
    "species": 9606,
    "required_score": 700,
    "caller_identity": "chatgpt_pipeline"
}

response = requests.post(
    f"{string_api_url}/{output_format}/{method}",
    data=params
)

# Save raw
with open("results/networks/disease_ppi_raw.tsv", "w") as f:
    f.write(response.text)

# Load + clean
ppi = pd.read_csv("results/networks/disease_ppi_raw.tsv", sep="\t")

ppi = ppi[['preferredName_A', 'preferredName_B', 'score']]
ppi.columns = ['gene1', 'gene2', 'weight']

ppi.to_csv("results/networks/disease_ppi.tsv", sep="\t", index=False)

print("Saved disease_ppi.tsv ✅")