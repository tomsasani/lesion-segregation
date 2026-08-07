
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from snakemake.script import snakemake
import numpy as np

plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans"]

res = pd.read_csv(snakemake.input.mavs, sep="\t")


f, ax = plt.subplots(figsize=(4, 6))
res["empirical"] = res["permuted"].apply(lambda p: "focal strain" if p == 0 else "random strain")
sns.stripplot(
    data=res[res["is_multiallelic"]],
    x="empirical",
    y="ma_reads",
    hue="empirical",
    palette="colorblind",
    ax=ax,
    alpha=0.5,
)
sns.boxplot(
    data=res[res["is_multiallelic"]],
    x="empirical",
    y="ma_reads",
    color="white",
    ax=ax,
    width=0.4,
    fliersize=0,
)

ax.set_xlabel("Source of read evidence")
ax.set_ylabel("# of reads supporting 'extra' allele")
sns.despine(ax=ax)
f.tight_layout()
f.savefig(snakemake.output.n_reads, dpi=200)

f, ax = plt.subplots(figsize=(5, 6))
for cat, cat_df in res[res["is_multiallelic"]].groupby("empirical"):
    if cat == "focal strain":
        ax.axvline(cat_df["ma_reads"].mean(), label="focal strain", ls=":", c="k")
    else:
        vals = cat_df.groupby("permuted")["ma_reads"].mean()
        ax.hist(vals, bins=np.arange(3, 8), ec="w", lw=1, label="random strain", color=sns.color_palette("colorblind", 1)[0])
ax.set_xlabel("Average # of reads supporting\n'extra' allele")
ax.set_ylabel("# of trials")
ax.legend(title="Source of read evidence", shadow=True)
sns.despine(ax=ax)
f.tight_layout()
f.savefig(snakemake.output.n_reads_hist, dpi=200)