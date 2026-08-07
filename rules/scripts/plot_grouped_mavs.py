import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from snakemake.script import snakemake
import numpy as np

plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans"]


N_PERM = int(snakemake.wildcards.PERM)

mav_counts = pd.read_csv(snakemake.input.counts, sep="\t")


metadata = pd.read_excel(
    snakemake.input.metadata,
    sheet_name="Sample_description",
    skiprows=1,
)

mav_counts = mav_counts.merge(metadata, left_on="sample", right_on="Sample")

group = snakemake.params.group


print (mav_counts)
f, ax = plt.subplots(figsize=(mav_counts[group].nunique() * 2, 5))

palette = sns.color_palette("colorblind", mav_counts[group].nunique())

enrichment = []

rng = np.random.default_rng(42)

order = []
for gi, (g, group_df) in enumerate(mav_counts.groupby(group)):

    tot = mav_counts[(mav_counts["permuted"] == 0) & (mav_counts[group] == g)][
        "total"
    ].sum()

    multiallelics = group_df[group_df["is_multiallelic"]]

    emp = multiallelics[multiallelics["permuted"] == 0].shape[0]
    prm = (
        multiallelics[multiallelics["permuted"] > 0]
        .groupby("permuted")
        .size()
        .reset_index()
        .rename(columns={0: "count"})
    )

    idx = pd.Index(list(range(N_PERM)))
    prm = prm.set_index("permuted").reindex(idx, fill_value=0).reset_index()
    prm_mean = np.mean(prm["count"].values) 
    enrich = round((emp + 1) / (np.mean(prm["count"].values) + 1), 1)

    ax.scatter(
        rng.normal(0, 0.125, size=prm.shape[0]) + gi,
        (prm["count"].values / tot) * 100,
        ec="w",
        label="random strain" if gi == 0 else None,
        color=palette[gi],
    )
    ax.scatter(
        gi,
        (emp / tot) * 100,
        color="k",
        marker="+",
        s=100,
        label="focal strain" if gi == 0 else None,
    )

    order.append(f"{g}\n({tot:,} SNVs)")

ax.set_xticks(range(len(order)))
ax.set_ylabel("% of SNVs with multi-allelic read support")
ax.set_xlabel(group)
ax.set_xticklabels(order)
ax.legend(shadow=True, title="Source of read evidence", loc="center")
sns.despine(ax=ax)

f.tight_layout()
f.savefig(snakemake.output.grouped_mavs, dpi=300)
