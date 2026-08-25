import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from snakemake.script import snakemake
import numpy as np

plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans"]

N_PERM = int(snakemake.wildcards.PERM)

res = pd.read_csv(snakemake.input.mavs, sep="\t")

empirical = res[res["permuted"] == 0]
permuted = res[res["permuted"] != 0]

# plot # of multiallelics remaining after applying various thresholds
f, ax = plt.subplots()
ax2 = ax.twinx()
thresholds = range(1, 6)
o = []

for t in thresholds:
    mas = empirical[empirical["ma_reads"] >= t].shape[0]
    perm_mas = (
        permuted[permuted["ma_reads"] >= t]
        .groupby("permuted")
        .size()
        .reset_index()
        .rename(columns={0: "count"})
    )
    # fill in potential 0 counts
    mux = pd.Index(list(range(1, N_PERM + 1)))
    perm_mas = perm_mas.set_index("permuted").reindex(mux, fill_value=0).reset_index()
    perm_mas = perm_mas["count"].values
    odds_ratio = mas / np.mean(perm_mas)
    o.append(
        {
            "t": t,
            "or": odds_ratio,
            "n": mas,
            "perm_mas": np.mean(perm_mas),
        }
    )
o = pd.DataFrame(o)
sns.lineplot(data=o, x="t", y="n", ax=ax, c="cornflowerblue")
sns.lineplot(data=o, x="t", y="or", ax=ax2, c="firebrick")

# at 3 read threshold
thresh = o[o["t"] == 3]
ax.axvline(x=3, ls=":", c="gainsboro", zorder=-1)
ax.scatter(3, thresh["n"].values[0], c="cornflowerblue", s=100, ec="w", lw=1)
ax2.scatter(3, thresh["or"].values[0], c="firebrick", s=100, ec="w", lw=1)
ax.axhline(y=thresh["n"].values[0], xmin=0, xmax=0.5, c="cornflowerblue", ls=":")
ax2.axhline(y=thresh["or"].values[0], xmin=0.5, xmax=1, c="firebrick", ls=":")
ax.set_ylabel(
    "Number of multi-allelic variants (MAVs)", fontdict={"color": "cornflowerblue"}
)
ax2.set_ylabel(
    "Enrichment of empirical vs. spurious MAVs",
    fontdict={"color": "firebrick"},
)
ax.set_xlabel("Number of reads supporting 'extra' allele")

ax.tick_params(axis="y", labelcolor="cornflowerblue")
ax2.tick_params(axis="y", labelcolor="firebrick")

sns.despine(ax=ax)
sns.despine(ax=ax2, top=True, right=False)
f.tight_layout()
f.savefig(snakemake.output.elbow, dpi=200)