import pandas as pd

file_path = "data/S2/expression_raw.csv"

# read only header
df = pd.read_csv(file_path, nrows=0)

print(df.columns.tolist())