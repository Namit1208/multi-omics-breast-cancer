import pandas as pd
import gseapy as gp
import os

os.makedirs("results/enrichment", exist_ok=True)

# Load modules
modules = pd.read_csv("results/modules/final_tx_ppi_modules.tsv", sep="\t")

# Load disease PPI → background genes
ppi = pd.read_csv("results/networks/disease_ppi.tsv", sep="\t")

background_genes = set(ppi['gene1']).union(set(ppi['gene2']))
background_genes = list(background_genes)

print("Background genes:", len(background_genes))

all_results = []

for m in modules['module'].unique():
    genes = modules[modules['module'] == m]['gene'].tolist()

    # Skip tiny modules
    if len(genes) < 5:
        continue

    try:
        enr = gp.enrichr(
            gene_list=genes,
            gene_sets=['GO_Biological_Process_2021'],
            organism='human',
            background=background_genes   # 🔥 KEY LINE
        )

        res = enr.results

        res['module_id'] = f"TX{m}"
        res = res[['module_id', 'Term', 'P-value']]

        all_results.append(res)

    except Exception as e:
        print(f"Error in module {m}: {e}")

# Combine all modules
if all_results:
    final = pd.concat(all_results, ignore_index=True)
    final.columns = ['module_id', 'term', 'pvalue']

    final.to_csv("results/enrichment/enrich_txppi_AD.tsv", sep="\t", index=False)
    print("Saved enrich_txppi_AD.tsv ✅")
else:
    print("No enrichment results generated ❌")