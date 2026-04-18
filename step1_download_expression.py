import GEOparse
import pandas as pd
import os

registry = pd.read_csv("manifests/expression_studies.tsv", sep="\t")

for _, row in registry.iterrows():
    
    study = row["study_id"]
    accession = row["accession"]

    print("Downloading", accession)

    gse = GEOparse.get_GEO(accession, destdir="data")

    os.makedirs(f"results/transcriptome/{study}", exist_ok=True)

    for gsm in gse.gsms:
        data = gse.gsms[gsm].table
        data.to_csv(f"results/transcriptome/{study}/expression_raw.tsv", sep="\t")
        break