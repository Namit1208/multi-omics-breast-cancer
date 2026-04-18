import os
import pandas as pd
import networkx as nx
import community as community_louvain

# create folder if missing
os.makedirs("results/modules", exist_ok=True)

df = pd.read_csv("results/networks/consensus_tx.tsv", sep="\t")

G = nx.Graph()

for _, row in df.iterrows():
    G.add_edge(row['gene1'], row['gene2'], weight=row['weight'])

partition = community_louvain.best_partition(G)

modules = pd.DataFrame(list(partition.items()), columns=["gene", "module"])

modules.to_csv("results/modules/tx_modules.tsv", sep="\t", index=False)
