import pysam
from collections import Counter, defaultdict, namedtuple
import numpy as np
import pandas as pd
from bx.intervals.intersection import Interval, IntervalTree
import csv
import tqdm
import os
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
import scipy.stats as ss

plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

from snakemake.script import snakemake

# catalog raw counts of mutations
raw_counts = pd.read_excel(
    "data/41467_2020_15912_MOESM4_ESM.xlsx",
    sheet_name="Mutation counts for all samples",
    skiprows=1,
).rename(columns={"Unnamed: 0": "Sample"})
raw_counts["count"] = np.sum(raw_counts.values[:, 1:99], axis=1)


metadata = pd.read_excel(
    snakemake.input.metadata, sheet_name="Sample_description", skiprows=1
).fillna({"Normal": "CD0001b", "Mutagen": "none", "Drug.concentration": -1})

metadata["mutagen"] = metadata.apply(
    lambda row: row["Mutagen"] if row["Type"] != "mut.acc." else "None (M.A.)", axis=1
)

# merge raw counts with metadata
raw_counts = raw_counts.merge(metadata)
raw_counts = raw_counts[
    (raw_counts["mutagen"] == "None (M.A.)")
    | (
        (raw_counts["mutagen"] != "None (M.A.)")
        & (raw_counts["Drug.concentration"] != 0)
    )
]
print(raw_counts[raw_counts["mutagen"] == "Mitomycin"])
raw_counts = raw_counts.groupby("mutagen").agg(raw_count=("count", "sum")).reset_index()

mutations = pd.read_csv(snakemake.input.mutations, sep="\t").fillna(
    {"mutagen": "None (M.A.)"}
)

# remove mutations in mutagen treated strains if dose is 0
mutations = mutations[
    (mutations["mutagen"] == "None (M.A.)")
    | ((mutations["mutagen"] != "None (M.A.)") & (mutations["dose"] != 0))
]

# count number of strains treated with each mutagen
strain_counts = (
    metadata.groupby("mutagen").size().reset_index().rename(columns={0: "strain_count"})
)

# count up mutations in every strain
counts = (
    mutations.groupby(["mutagen"]).size().reset_index().rename(columns={0: "mut_count"})
)
merged = strain_counts.merge(counts)

merged = merged.merge(raw_counts)

merged["mut_count"] /= 10**4
merged["raw_count"] /= 10**4


merged = merged[merged["mutagen"] != "MMS/EMS/DMS"]

f, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 6), sharey=True)

order = merged["mutagen"].unique()
for gi, g in enumerate(order):
    raw = merged[merged["mutagen"] == g]["raw_count"].values[0]
    flt = merged[merged["mutagen"] == g]["mut_count"].values[0]
    if g not in ("EMS", "MMS", "DMS"):
        continue
    ax2.text(
        raw + 0.15,
        gi + 0.15,
        f"{round(flt * 10**4):,} ({round(raw * 10**4):,})",
    )

sns.barplot(
    data=merged,
    order=order,
    y="mutagen",
    x="strain_count",
    ax=ax1,
    orient="h",
    hue="mutagen",
    ec="k",
    lw=1,
    palette="colorblind",
)
sns.barplot(
    data=merged,
    order=order,
    y="mutagen",
    x="raw_count",
    ax=ax2,
    orient="h",
    hue="mutagen",
    alpha=0.75,
    ec="w",
    lw=0,
    palette="colorblind",
)
sns.barplot(
    data=merged,
    order=order,
    y="mutagen",
    x="mut_count",
    ax=ax2,
    orient="h",
    hue="mutagen",
    ec="k",
    lw=1,
    palette="colorblind",
)
ax1.set_xlabel("Number of strains")
ax2.set_xlabel("Number of SNVs " + r"($\times 10^4$)")
ax1.set_ylabel("Genotoxin treatment")
for ax in (ax1, ax2):
    sns.despine(ax=ax)
f.tight_layout()
f.savefig(snakemake.output.fig, dpi=200)
