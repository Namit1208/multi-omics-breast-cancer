import pandas as pd
import gseapy as gp
import os

os.makedirs("results/enrichment", exist_ok=True)

modules = pd.read_csv("results/modules/modules_txppi_AD.tsv", sep="\t")

for m in modules['module'].unique():
    genes = modules[modules['module'] == m]['gene'].tolist()
    
    if len(genes) < 5:
        continue  # skip tiny modules (important)
    
    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2021'],
        organism='human'
    )
    
    enr.results.to_csv(f"results/enrichment/module_{m}.tsv", sep="\t", index=False)