import pandas as pd
tx = pd.read_csv("results/networks/consensus_tx.tsv", sep="\t")
ppi = pd.read_csv("results/networks/ppi_tx_overlap.tsv", sep="\t")

combined = pd.concat([tx, ppi], ignore_index=True)

combined.to_csv("results/networks/tx_ppi_combined.tsv", sep="\t", index=False)