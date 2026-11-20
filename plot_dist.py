import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

df = pd.read_csv("combinatorial.tsv", sep="\t")

df = df[df["combinatorial"] == True]
df = df[df["Mutagen"].isin(["EMS", "DMS", "MMS"])]
df = df[df["dist"] < 200]
print (df[df["dist"] > 100])
g = sns.FacetGrid(data=df, col="Mutagen", sharex=True)
g.map(sns.histplot, "dist", bins=np.arange(1, 200, 2))
g.set(xscale="log")
g.savefig("genome.dist.png")
