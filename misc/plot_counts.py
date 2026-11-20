import pandas as pd

df = pd.read_csv("csv/mutations/combined.tsv", sep="\t")
df["mutagen"] = df["sample"].apply(lambda s: s.split("_")[0])
print (df.groupby("mutagen").size())