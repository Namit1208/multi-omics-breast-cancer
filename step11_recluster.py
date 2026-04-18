import os
import pandas as pd
import igraph as ig
import leidenalg

os.makedirs("results/modules", exist_ok=True)

df = pd.read_csv("results/networks/intersection_tx_ppi_AD.tsv", sep="\t")

edges = list(zip(df['gene1'], df['gene2']))
g = ig.Graph.TupleList(edges, directed=False)

partition = leidenalg.find_partition(g, leidenalg.ModularityVertexPartition)

modules = []
for i, cluster in enumerate(partition):
    for node in cluster:
        modules.append([g.vs[node]['name'], i])

modules = pd.DataFrame(modules, columns=["gene", "module"])

modules.to_csv("results/modules/modules_txppi_AD.tsv", sep="\t", index=False)

print("Leiden clustering done ✅")