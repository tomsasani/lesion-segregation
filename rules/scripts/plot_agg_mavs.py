import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from snakemake.script import snakemake
import numpy as np

plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans"]


mav_counts = pd.read_csv(snakemake.input.counts, sep="\t")

metadata = pd.read_excel(
    snakemake.input.metadata,
    sheet_name="Sample_description",
    skiprows=1,
)

mav_counts = mav_counts.merge(metadata, left_on="sample", right_on="Sample")

# filter to the multiallelics
mav_counts = mav_counts[mav_counts["is_multiallelic"]]
# get counts of MAVs in each permutation trial
perm_ma_counts = mav_counts.groupby("permuted").agg(n_mavs=("count", "sum")).reset_index()

f, ax = plt.subplots()

k = "count"
perm_ma_counts["is_empirical"] = perm_ma_counts["permuted"] == 0
_min, _max = perm_ma_counts["n_mavs"].min(), perm_ma_counts["n_mavs"].max()
bins = np.linspace(_min, _max, 100)

ax.hist(
    perm_ma_counts[~perm_ma_counts["is_empirical"]]["n_mavs"].values,
    ec="k",
    lw=1,
    bins=bins,
    label="random strain",
    color="gainsboro",
)
ax.axvline(
    perm_ma_counts[perm_ma_counts["is_empirical"]]["n_mavs"].values[0],
    ls=":",
    c="k",
    label="focal strain",
)
ax.set_xlabel("Number of SNVs with multi-allelic support")
ax.set_ylabel("Count")
ax.legend(title="Source of read evidence", shadow=True, loc="center")
sns.despine(ax=ax)
f.tight_layout()
f.savefig(snakemake.output.agg_mavs, dpi=200)
