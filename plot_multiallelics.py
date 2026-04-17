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

res = []
for i, fh in tqdm.tqdm(enumerate(glob.glob("csv/multiallelics/*.20.tsv"))):
    try:
        df = pd.read_csv(fh, sep="\t")
    except (pd.errors.EmptyDataError, FileNotFoundError):
        continue
    df["is_multiallelic"] = df["ma_reads"] >= THRESHOLD
    res.append(df)
res = pd.concat(res).fillna("N:0")

res["depth"] = res["summary"].apply(lambda s: sum([int(_s.split(":")[1]) for _s in s.split("|")]))
res = res[res["depth"].between(15, 150)]

###
# are multiallelics more common in certain tech?
# NOTE: this is tough, because later-sequenced samples on Xten might be a certain mutagen/genotype
# that's more likely to have MAVs
###
empirical = res[res["permuted"] == 0]
manifest = pd.read_csv("data/combined_file_manifest.tsv", sep="\t")
smp2tech = dict(zip(manifest["Sample / DONOR ID"], manifest["experiment_title"]))
empirical["tech"] = empirical["sample"].apply(lambda s: smp2tech[s])


empirical["frac"] = empirical.apply(lambda row: int(row["summary"].split("|")[-1].split(":")[1]) / row["depth"], axis=1)

f, ax = plt.subplots()
sns.stripplot(data=empirical[empirical["is_multiallelic"]], x="tech", y="frac", ax=ax, alpha=0.5)
sns.boxplot(data=empirical[empirical["is_multiallelic"]], x="tech", y="frac", ax=ax)
ax.set_ylabel("Fraction of reads supporting third allele")
ax.set_xlabel("Sequencing technology")
f.savefig("ev.png")

# tech = []
# for smp, smp_df in empirical.groupby("sample"):
#     tot = smp_df.shape[0]
#     mav = smp_df[smp_df["is_multiallelic"]].shape[0]
#     tech.append(
#         {
#             "sample": smp,
#             "tech": "HiSeq 2000" if "2000" in smp2tech[smp] else "HiSeq X10",
#             "frac": mav / tot,
#             "tot": tot,
#         }
#     )
# tech = pd.DataFrame(tech)
# value = "ev"
# tech = tech.merge(metadata, left_on="sample", right_on="Sample")
# print (tech)

# replicates = (
#     tech.groupby(["Mutagen", "Genotype", "Drug.concentration"])
#     .agg(dupped=("tech", lambda s: len(set(s)) == 2))
#     .reset_index()
# )
# replicates = replicates[replicates["dupped"]]

# replicate_samples = tech.merge(replicates)

# # figure out the minimum number of samples per replicate in either the hiseq or xten dataset
# min_samples = (
#     replicate_samples.groupby(["Mutagen", "Genotype", "Drug.concentration"])
#     .agg(
#         minim=("tech", lambda t: Counter(t).most_common()[-1][-1]),
#     )
#     .reset_index()
# )
# replicate_samples = replicate_samples.merge(min_samples)
# # make sure sample sizes are equivalent between techs

# downsampled = []
# for (m, g, d), mgdf in replicate_samples.groupby(["Mutagen", "Genotype", "Drug.concentration"]):
#     min_smp = mgdf["minim"].unique()
#     assert len(min_smp) == 1
#     min_smp = min_smp[0]
#     for _, tech_df in mgdf.groupby("tech"):
#         if tech_df.shape[0] > min_smp:
#             tech_df = tech_df.sample(n=min_smp, replace=False)
#             downsampled.append(tech_df)
#         else:
#             downsampled.append(tech_df)
# downsampled = pd.concat(downsampled)

# replicate_samples = downsampled["sample"].unique()
# tech = tech[tech["sample"].isin(replicate_samples)]
# bins = np.linspace(tech[value].min(), tech[value].max(), 100)
# f, ax = plt.subplots()
# trials = 100
# for ti, (t, tdf) in enumerate(tech.groupby("tech")):

#     cumsum = np.zeros((trials, bins.shape[0] - 1))

#     for bs in range(trials):
#         _tdf = tdf.sample(frac=1, replace=True)
#         hist, _ = np.histogram(_tdf[value].values, bins=bins)
#         hist_fracs = hist / np.sum(hist)
#         cumsum[bs, :] = np.cumsum(hist_fracs)
#     ax.plot(bins[:-1], np.mean(cumsum, axis=0), label=t)
#     ax.fill_between(
#         bins[:-1],
#         np.percentile(cumsum, 2.5, axis=0),
#         np.percentile(cumsum, 97.5, axis=0),
#         alpha=0.25,
#     )
# ax.legend()
# # ax.set_yscale("symlog")
# f.savefig("tech.png", dpi=200)


# filter to high-ish depth

res = res.merge(metadata, left_on="sample", right_on="Sample")

# plot # of multiallelics remaining after applying various thresholds
# look for elbow!
ma = res[res["permuted"] == 0]
nma = res[res["permuted"] != 0]
f, ax = plt.subplots()
ax2 = ax.twinx()
thresholds = range(1, 9)
o = []
for t in thresholds:
    mas = ma[ma["ma_reads"] >= t].shape[0]
    perm_mas = nma[nma["ma_reads"] >= t].groupby("permuted").size().reset_index().rename(columns={0: "count"})["count"].values
    odds_ratio = mas / np.mean(perm_mas)
    o.append({"t": t, "or": odds_ratio, "p": np.sum(perm_mas >= mas), "n": mas})
o = pd.DataFrame(o)
sns.lineplot(data=o, x="t", y="n", ax=ax, c="dodgerblue")
sns.lineplot(data=o, x="t", y="or", ax=ax2, c="firebrick")
ax.set_ylabel("Number of multi-allelic mutations", fontdict={"color": "dodgerblue"})
ax2.set_ylabel("Enrichment " + r"$\frac{MA_{emp}}{\mathbb{E}[MA_{rand}]}$", fontdict={"color": "firebrick"})
ax.set_xlabel("Number of reads supporting 'third' allele")
sns.despine(ax=ax)
sns.despine(ax=ax2, top=True, right=False)
f.tight_layout()
f.savefig("elbow.png", dpi=200)

# just ask whether total count of multiallelics (across samples) is more than expected
counts = (
    res.groupby(["permuted", "is_multiallelic"])
    .size()
    .reset_index()
    .rename(columns={0: "count"})
)
totals = (
    counts.groupby(["permuted"]).agg(total=("count", "sum")).reset_index()
)
counts = counts.merge(totals)


counts = counts[~counts["is_multiallelic"]]
counts["true_count"] = counts["total"] - counts["count"]
f, ax = plt.subplots()
k = "true_count"
counts["is_empirical"] = counts["permuted"] == 0
_min, _max = counts[k].min(), counts[k].max()
bins = np.linspace(int(_min * 0.9), int(_max * 1.1), 50)
ax.hist(counts[~counts["is_empirical"]][k].values, ec="w", lw=1, bins=bins, label="random")
ax.axvline(counts[counts["is_empirical"]][k].values[0], ls=":", c="k", label="empirical")
ax.set_xlabel(f"Number of SNVs with multi-allelic support")
ax.set_ylabel("Count")
ax.legend(title="SNVs", shadow=True)
sns.despine(ax=ax)
f.tight_layout()
f.savefig("multi_perm_overall.png", dpi=200)


support = []
for p, pdf in res.groupby("permuted"):
    pdf = pdf[pdf["is_multiallelic"]]
    if pdf.shape[0] == 0:
        avg_support = 0
        continue

    else:
        pdf["frac_support"] = pdf.apply(lambda row: int(row["summary"].split("|")[-1].split(":")[1]) / row["depth"], axis=1)

        avg_support = pdf["ma_reads"].mean()
        avg_support = pdf["frac_support"].mean()

    support.append(
        {
            "permuted": p,
            "status": p == 0,
            "support": avg_support,
        }
    )

support = pd.DataFrame(support)
f, ax = plt.subplots()
ax.hist(support[~support["status"]]["support"], lw=1, ec="w", label="Random")
ax.axvline(support[support["status"]]["support"].values[0], ls=":", c="k", label="Empirical")
ax.set_ylabel("Count")
ax.legend(shadow=True)
ax.set_xlabel("Average # of reads supporting 'extra' allele")
sns.despine(ax=ax)
f.tight_layout()
f.savefig("read_evidence.png", dpi=200)


###
# for a given "group" (based on genotype or mutagen type), figure out how many
# MAVs were detected using the "true" and randomly samples SNVs.
###


# res = res[res["mutagen"] == "MMS"]

group = "mutagen"

n = max([2, res[group].nunique()])

bins = np.linspace(0, 20, 50)

enrichment = []
f, axarr = plt.subplots(n, figsize=(7, 3 * n), sharex=True, sharey=True)
for gi, (group, group_df) in enumerate(res.groupby(group)):

    tot = group_df[group_df["permuted"] == 0].shape[0]
    

    multiallelics = group_df[group_df["is_multiallelic"]]

    emp = multiallelics[multiallelics["permuted"] == 0].shape[0]
    prm = multiallelics[multiallelics["permuted"] > 0].groupby("permuted").size().reset_index().rename(columns={0: "count"})

    idx = pd.Index(list(range(100)))
    prm = prm.set_index("permuted").reindex(idx, fill_value=0).reset_index()

    enrich = emp / np.mean(prm["count"].values)
    enrichment.append({"group": group, "enrichment": enrich, "p": np.sum(prm["count"].values >= emp) / 100})

    ax = axarr[gi]
    
    ax.hist(prm["count"].values, bins=bins, lw=1, ec="w", label="permuted")
    ax.axvline(emp, ls=":", c="k", label=f"empirical")
    sns.despine(ax=ax)
    ax.set_title(f"{group} (n = {tot:,} SNVs)\n({enrich})")
    ax.set_xlabel("Number of SNVs with multi-allelic support")
    ax.set_ylabel("Count")
    if gi == 0:
        ax.legend(shadow=True)
f.tight_layout()
f.savefig("multiallelic_by_mutagen.png", dpi=200)

enrichment = pd.DataFrame(enrichment)
enrichment["sig"] = enrichment["p"] < 0.05
f, ax = plt.subplots(figsize=(4, 12))
sns.barplot(data=enrichment, y="group", x="enrichment", orient="h", hue="sig", ax=ax)
sns.despine(ax=ax)
f.tight_layout()
f.savefig("enrichment.png", dpi=200)
