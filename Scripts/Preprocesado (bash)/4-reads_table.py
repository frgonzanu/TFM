import pandas as pd
import glob
from functools import reduce
import os

# buscar en subcarpetas
files = glob.glob("**/*ReadsPerGene.out.tab", recursive=True)

dfs = []

for f in files:
    sample = os.path.basename(f).replace("_ReadsPerGene.out.tab", "")
    
    df = pd.read_csv(f, sep="\t", header=None)
    
    df.columns = ["gene", "unstranded", "strand_plus", "strand_minus"]
    
    # usar solo unstranded
    df = df[["gene", "unstranded"]]
    df = df.rename(columns={"unstranded": sample})
    
    dfs.append(df)

# merge
merged = reduce(lambda left, right: pd.merge(left, right, on="gene"), dfs)

merged.to_csv("counts_matrix.tsv", sep="\t", index=False)
