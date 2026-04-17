import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans"]

metadata = pd.read_excel(
    "data/41467_2020_15912_MOESM4_ESM.xlsx", sheet_name="Sample_description", skiprows=1
).fillna({"Normal": "CD0001b", "Mutagen": "none", "Drug.concentration": -1})

metadata = metadata[metadata["Drug.concentration"] != 0]
metadata["Genotoxin"] = metadata.apply(lambda row: row["Mutagen"] if row["Type"] != "mut.acc." else "None (M.A.)", axis=1)

raw_counts = pd.read_excel("data/41467_2020_15912_MOESM4_ESM.xlsx", 
                         sheet_name="Mutation counts for all samples", 
                         skiprows=1).rename(columns={"Unnamed: 0": "Sample"})
raw_counts["count"] = np.sum(raw_counts.values[:, 1:97], axis=1)

raw_counts = raw_counts.merge(metadata)
mut_counts = raw_counts.groupby("Genotoxin").agg(total=("count", "sum")).reset_index()
mut_counts["total"] /= 10**4

f, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 6), sharey=True)

counts = metadata.groupby("Genotoxin").size().reset_index().rename(columns={0: "count"})
sns.barplot(data=counts, y="Genotoxin", x="count", ax=ax1, orient="h", hue="Genotoxin", ec="k", lw=1, palette="colorblind")
sns.barplot(data=mut_counts, y="Genotoxin", x="total", ax=ax2, orient="h", hue="Genotoxin", ec="k", lw=1, palette="colorblind")
# ax2.set_xscale("log")
ax1.set_xlabel("Number of strains")
ax2.set_xlabel("Number of SNV mutations " + r"($\times 10^4$)")
ax1.set_ylabel("Genotoxin treatment")
for ax in (ax1, ax2):
    sns.despine(ax=ax)
f.tight_layout()
f.savefig("meta.png", dpi=200)