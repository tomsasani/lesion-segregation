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

mav_counts = mav_counts[mav_counts["is_multiallelic"]]
_ma_counts = mav_counts.groupby("permuted").agg(ma=("count", "sum")).reset_index()

f, ax = plt.subplots()

k = "count"
_ma_counts["is_empirical"] = _ma_counts["permuted"] == 0
_min, _max = _ma_counts["ma"].min(), _ma_counts["ma"].max()
bins = np.linspace(_min, _ma_counts[_ma_counts["is_empirical"]]["ma"].values[0], 100)
bins = np.linspace(_min, _max, 100)
ax.hist(
    _ma_counts[~_ma_counts["is_empirical"]]["ma"].values,
    ec="k",
    lw=1,
    bins=bins,
    label="random strain",
    color="gainsboro",
)
ax.axvline(
    _ma_counts[_ma_counts["is_empirical"]]["ma"].values[0],
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
