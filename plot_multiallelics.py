import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import numpy as np
import tqdm
import scipy.stats as ss
import os
from collections import Counter


def get_ab(summary):
    base_counts = [int(b.split(":")[1]) for b in summary.split("|")]
    min_ab = base_counts[-1] / sum(base_counts)
    return min_ab


plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans"]

THRESHOLD = 3

metadata = pd.read_excel(
    "data/41467_2020_15912_MOESM4_ESM.xlsx",
    sheet_name="Sample_description",
    skiprows=1,
)

N_PERM = 10

MUTAGENS = ["EMS", "DMS", "MMS"]

res = []
samples = []
ma_counts = []
duplicates = []
total_mutations = 0
for i, fh in tqdm.tqdm(enumerate(glob.glob(f"csv/multiallelics/*.{N_PERM}.by_bam.tsv"))):
    sample = fh.split("/")[-1].split(".")[0]
    samples.append(sample)
    try:
        df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    except (pd.errors.EmptyDataError, FileNotFoundError):
        continue

    if df["mutagen"].to_list().pop() not in MUTAGENS:
        continue

    df["is_multiallelic"] = df["ma_reads"] >= THRESHOLD
    total_mutations += df[df["permuted"] == 0].shape[0]
    # store counts
    counts = (
        df.groupby(["sample", "permuted", "is_multiallelic"])
        .size()
        .reset_index()
        .rename(columns={0: "count"})
    )
    mux = pd.MultiIndex.from_product(
        [[sample], list(range(N_PERM + 1)), [True, False]],
        names=["sample", "permuted", "is_multiallelic"],
    )
    counts = (
        counts.set_index(["sample", "permuted", "is_multiallelic"])
        .reindex(mux, fill_value=0)
        .reset_index()
    )
    ma_counts.append(counts)

    df = df[df["ma_reads"] > 0]
    res.append(df)


res = pd.concat(res).fillna("N:0")
res = res.merge(metadata, left_on="sample", right_on="Sample")


ma_counts = pd.concat(ma_counts)
ma_counts = ma_counts.merge(metadata, left_on="sample", right_on="Sample")
ma_totals = (
    ma_counts.groupby(["sample", "permuted"]).agg(total=("count", "sum")).reset_index()
)
ma_counts = ma_counts.merge(ma_totals)

ma_counts[(ma_counts["permuted"] == 0) & (ma_counts["is_multiallelic"])].to_csv("multiallelic.tsv", sep="\t", index=False)
empirical = res[res["permuted"] == 0]
# empirical = empirical.merge(ma_counts[["sample", "total"]])

# empirical[empirical["is_multiallelic"]].to_csv("multiallelic.tsv", sep="\t", index=False)



sb = (
    empirical.groupby(["is_multiallelic", "is_strand_biased"])
    .size()
    .reset_index()
    .rename(columns={0: "count"})
)

res["depth"] = res["summary"].apply(
    lambda s: sum([int(_s.split(":")[1]) for _s in s.split("|")])
)

# plot # of multiallelics remaining after applying various thresholds
# look for elbow!
f, ax = plt.subplots()
ax2 = ax.twinx()
thresholds = range(1, 6)
o = []

permuted = res[res["permuted"] != 0]
for t in thresholds:
    mas = empirical[empirical["ma_reads"] >= t].shape[0]
    perm_mas = (
        permuted[permuted["ma_reads"] >= t]
        .groupby("permuted")
        .size()
        .reset_index()
        .rename(columns={0: "count"})
    )

    mux = pd.Index(list(range(N_PERM + 1)))
    perm_mas = perm_mas.set_index("permuted").reindex(mux, fill_value=0).reset_index()
    perm_mas = perm_mas["count"].values
    odds_ratio = mas / np.mean(perm_mas)
    o.append(
        {
            "t": t,
            "or": odds_ratio,
            "n": mas,
        }
    )
o = pd.DataFrame(o)
sns.lineplot(data=o, x="t", y="n", ax=ax, c="cornflowerblue")
sns.lineplot(data=o, x="t", y="or", ax=ax2, c="firebrick")

# at 3 read threshold
thresh = o[o["t"] == 3]
print(thresh)
ax.axvline(x=3, ls=":", c="gainsboro", zorder=-1)
ax.scatter(3, thresh["n"].values[0], c="cornflowerblue", s=100, ec="w", lw=1)
ax2.scatter(3, thresh["or"].values[0], c="firebrick", s=100, ec="w", lw=1)
ax.axhline(y=thresh["n"].values[0], xmin=0, xmax=0.5, c="cornflowerblue", ls=":")
ax2.axhline(y=thresh["or"].values[0], xmin=0.5, xmax=1, c="firebrick", ls=":")
ax.set_ylabel(
    "Number of multi-allelic variants (MAVs)", fontdict={"color": "cornflowerblue"}
)
ax2.set_ylabel(
    "Enrichment " + r"$\frac{\# MAV}{\mathbb{E}(\# MAV)}$",
    fontdict={"color": "firebrick"},
)
ax.set_xlabel("Number of reads supporting 'extra' allele")

ax.tick_params(axis="y", labelcolor="cornflowerblue")
ax2.tick_params(axis="y", labelcolor="firebrick")

sns.despine(ax=ax)
sns.despine(ax=ax2, top=True, right=False)
f.tight_layout()
f.savefig("elbow.png", dpi=200)

# just ask whether total count of multiallelics (across samples) is more than expected

ma_counts = ma_counts[ma_counts["is_multiallelic"]]
_ma_counts = ma_counts.groupby("permuted").agg(ma=("count", "sum")).reset_index()

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
    label="random",
    color="gainsboro",
)
ax.axvline(
    _ma_counts[_ma_counts["is_empirical"]]["ma"].values[0],
    ls=":",
    c="k",
    label="empirical",
)
ax.set_xlabel("Number of SNVs with multi-allelic support")
ax.set_ylabel("Count")
ax.legend(title="SNVs", shadow=True)
sns.despine(ax=ax)
f.tight_layout()
f.savefig("multi_perm_overall.png", dpi=200)


# support = []
# for p, pdf in res.groupby("permuted"):
#     pdf = pdf[pdf["is_multiallelic"]]
#     if pdf.shape[0] == 0:
#         avg_support = 0
#         continue

#     else:
#         pdf["frac_support"] = pdf.apply(lambda row: int(row["summary"].split("|")[-1].split(":")[1]) / row["depth"], axis=1)

#         avg_support = pdf["ma_reads"].mean()
#         # avg_support = pdf["frac_support"].mean()

#     support.append(
#         {
#             "permuted": p,
#             "status": p == 0,
#             "support": avg_support,
#         }
#     )

# support = pd.DataFrame(support)
# f, ax = plt.subplots()
# ax.hist(support[~support["status"]]["support"], lw=1, ec="w", label="Random")
# ax.axvline(support[support["status"]]["support"].values[0], ls=":", c="k", label="Empirical")
# ax.set_ylabel("Count")
# ax.legend(shadow=True)
# ax.set_xlabel("Average # of reads supporting 'extra' allele")
# sns.despine(ax=ax)
# f.tight_layout()
# f.savefig("read_evidence.png", dpi=200)


###
# for a given "group" (based on genotype or mutagen type), figure out how many
# MAVs were detected using the "true" and randomly sampled SNVs.
###


group = "Mutagen"

prm = res[(res["permuted"] > 0) & (res["is_multiallelic"])]
prm_counts = prm.groupby([group, "permuted"]).size().reset_index().rename(columns={0: "count"})
emp = res[(res["permuted"] == 0) & (res["is_multiallelic"])]
emp_counts = emp.groupby(group).size().reset_index().rename(columns={0: "count"})


f, ax = plt.subplots()
sns.stripplot(data=prm_counts, x=group, y="count", hue=group, ax=ax, alpha=0.25)
sns.pointplot(data=prm_counts, x=group, y="count", hue=group, ax=ax, color="gainsboro")
sns.stripplot(data=emp_counts, x=group, y="count", ax=ax)
f.savefig("mabym.png", dpi=200)


f, ax = plt.subplots(figsize=(6, 5))

palette = sns.color_palette("colorblind", 5)

enrichment = []

rng = np.random.default_rng(42)

order = []
for gi, (g, group_df) in enumerate(res.groupby(group)):

    tot = ma_counts[(ma_counts["permuted"] == 0) & (ma_counts[group] == g)][
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
        label="random" if gi == 0 else None,
        color=palette[gi],
    )
    ax.scatter(
        gi,
        (emp / tot) * 100,
        color="k",
        marker="+",
        s=100,
        label="empirical" if gi == 0 else None,
    )
    # ax.scatter(gi, (prm_mean / tot) * 100, color="gainsboro", marker="_", s=100)

    order.append(f"{g}\n({tot:,} SNVs)")
# ax.set_yscale("log")
ax.set_xticks(range(len(order)))
ax.set_ylabel("% of SNVs with multi-allelic support")
ax.set_xlabel("Genotoxin")
ax.set_xticklabels(order)
ax.legend(shadow=True, title="MAV type", loc="center")
sns.despine(ax=ax)

f.tight_layout()
f.savefig("multiallelic_by_mutagen.png", dpi=200)
