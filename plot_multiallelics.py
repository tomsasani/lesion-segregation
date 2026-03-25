import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import numpy as np
import tqdm
import scipy.stats as ss
import os


def get_ab(summary):
    base_counts = [int(b.split(":")[1]) for b in summary.split("|")]
    min_ab = base_counts[-1] / sum(base_counts)
    return min_ab


plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]


metadata = pd.read_excel(
    "data/41467_2020_15912_MOESM4_ESM.xlsx",
    sheet_name="Sample_description",
    skiprows=1,
)

res = []
for fh in glob.glob("csv/multiallelics/*.tsv"):
    if not os.path.exists(fh): continue
    if os.path.getsize(fh) == 1: continue
    df = pd.read_csv(fh, sep="\t")
    res.append(df)
res = pd.concat(res).fillna("N:0")

res = res.merge(metadata, left_on="sample", right_on="Sample")

res = res[res["Mutagen"].isin(["EMS", "DMS", "MMS"])]

empirical = res[res["permuted"] == 0]
print (empirical.groupby(["is_multiallelic", "is_strand_biased"]).size())

empirical[empirical["is_multiallelic"] == True].to_csv("multiallelics.tsv", sep="\t", index=False)

# get samples with actual multiallelics
relevant_samples = empirical[empirical["is_multiallelic"] == True]["sample"].unique()
sample_counts = empirical.groupby("sample").size().reset_index().rename(columns={0: "count"})
relevant_samples = sample_counts.query("count >= 100")["sample"].unique()
# res = res[res["sample"].isin(relevant_samples)]


print (res.groupby("Genotype").size())

# bootstrap resample empirical counts

groups = ["Mutagen", "Genotype", "Drug.concentration"]

bs = []

for g, gdf in empirical.groupby(groups):
    # boostrap resample
    for t in range(100):
        _gdf = gdf.sample(frac=1, replace=True)
        n_multiallelic = _gdf["is_multiallelic"].sum()
        n_total = _gdf.shape[0]
        frac_multiallelic = n_multiallelic / n_total
        out_dict = {
            "t": t,
            "total": n_total,
            "frac": frac_multiallelic,
            "count": n_multiallelic,
            "true_frac": 1 - frac_multiallelic,
        }
        for k,v in zip(groups, g):
            out_dict.update({k: v})
        bs.append(out_dict)


bs = pd.DataFrame(bs)

f, axarr = plt.subplots(1, bs["Mutagen"].nunique(), figsize=(16, 5), sharex=True, sharey=True)
for mi, (mutagen, mutagen_df) in enumerate(bs.groupby("Mutagen")):

    max_doses = mutagen_df.groupby("Genotype").agg(max_dose=("Drug.concentration", "max")).reset_index()
    mutagen_df = mutagen_df.merge(max_doses)
    mutagen_df = mutagen_df[mutagen_df["Drug.concentration"] == mutagen_df["max_dose"]]

    ax = axarr[mi]
    sns.pointplot(
        data=mutagen_df,
        x="Genotype",
        y="frac",
        # hue="sample",
        legend=False,
        dodge=True,
        ax=ax,
        linestyle="none",
        errorbar=("ci", 95)
    )
    sns.despine(ax=ax)
    ax.set_title(f"Strains treated with {mutagen}")
    ax.set_xlabel("Strain genotypes")
axarr[0].set_ylabel("Fraction of SNVs with multi-allelic support\n(mean +/- 95% CI)")
f.suptitle("Multi-allelic mutation support in strains treated with the highest dose of each mutagen")
f.tight_layout()
f.savefig("multiallelic_by_mutagen.png", dpi=200)

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


counts = counts[counts["is_multiallelic"]]
counts["frac"] = 1 - counts["count"] / counts["total"]
f, ax = plt.subplots()
counts["status"] = (counts["permuted"] == 0).astype(int)
sns.stripplot(data=counts, x="status", y="count", ax=ax)
ax.set_ylabel("Number of SNVs with multi-allelic support")
ax.set_xlabel("")
ax.set_xticks([0, 1])
ax.set_xticklabels(["Permuted", "Empirical"])
f.savefig("multi_perm_overall.png", dpi=200)


# count the number of multiallelic variants

# NOTE we can't use the mutagen/drug.concentration/genotype columns
# here, becaues sample labels were permuted across conditions during multiallelic cataloguing
counts = (
    res.groupby(["sample", "permuted", "is_multiallelic"])
    .size()
    .reset_index()
    .rename(columns={0: "count"})
)
totals = (
    counts.groupby(["sample", "permuted"]).agg(total=("count", "sum")).reset_index()
)
counts = counts.merge(totals)


counts = counts[~counts["is_multiallelic"]]
counts["frac"] = 1 - counts["count"] / counts["total"]


rng = np.random.default_rng(42)
res = []
f, ax = plt.subplots()
for s, sdf in counts.groupby("sample"):
    perm = sdf[sdf["permuted"] > 0]["frac"].values
    mean_perm = np.mean(perm) + 1e-5
    emp = sdf[sdf["permuted"] == 0]["frac"].values
    print (sdf[sdf["permuted"] == 0])
    assert emp.shape[0] == 1
    emp = emp[0] + 1e-5
    # if emp == 0 and mean_perm == 0: continue
    x1, x2 = rng.normal(0, scale=0.05), rng.normal(1, scale=0.05)
    ax.plot([x1, x2], [emp, mean_perm], c="gainsboro", zorder=-1)
    ax.scatter([x1, x2], [emp, mean_perm], c="darkgrey", zorder=1)
    for kind, val in zip(("emp", "perm"), (emp, mean_perm)):
        res.append({"sample": s, "status": kind, "val": val})

f.savefig('test.png')
res = pd.DataFrame(res)

f, ax = plt.subplots()
sns.stripplot(data=res, x="status", y="val", hue="status", ax=ax, zorder=0, alpha=0.5)

for g, gdf in res.groupby("status"):
    avg = gdf["val"].mean()
    ax.scatter(0 if g == "emp" else 1, avg, zorder=1, c="k")
ax.set_xticks([0, 1])
ax.set_xlim(-0.5, 1.5)
ax.set_xticklabels(["Empirical", "Permuted"])
ax.set_xlabel("")
ax.set_ylabel("Fraction of sites with multi-allelic variation")
ax.set_yscale("log")
sns.despine(ax=ax)
f.tight_layout()
f.savefig("multi_perm.png", dpi=200)
