import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans"]


N_PERM = 0


res = []
for fh in glob.glob(f"csv/multiallelics/*.{N_PERM}.*.tsv"):
    try:
        df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    except (pd.errors.EmptyDataError, FileNotFoundError):
        continue

    df["strategy"] = fh.split("/")[-1].split(".")[-2]
    df["is_multiallelic"] = df["ma_reads"] >= 3
    res.append(df)
res = pd.concat(res)


grouped = res.groupby(["sample", "strategy"]).agg(count=("is_multiallelic", "sum"), total=("is_multiallelic", lambda m: len(m))).reset_index()
grouped["frac"] = grouped["count"] / grouped["total"]
print (grouped)
print (grouped.groupby("strategy").agg(mav=("count", "sum")))
pivoted = pd.pivot(grouped, index="sample", columns=["strategy"], values=["frac"]).reset_index()
pivoted.columns = ["sample", "downsample", "regular"]

f, ax = plt.subplots(figsize=(6, 5))
ax.axline((0, 0), slope=1, ls=":", c="k", zorder=-1)

sns.scatterplot(data=pivoted, x="regular", y="downsample", ax=ax, alpha=0.5, lw=1, ec="w", zorder=1)
sns.despine(ax=ax)

ax.set_xlabel("MAV fraction (calculated using original CRAM)")
ax.set_ylabel("MAV fraction (calculated using downsampled CRAM)")
f.tight_layout()
f.savefig("ds.png", dpi=200)

manifest = pd.read_csv("data/combined_file_manifest.tsv", sep="\t").rename(columns={"Sample / DONOR ID": "sample"})
res = res.merge(manifest)

res = res[res["strategy"] == "downsample"]

tech_comp = []

for (m, g, d), mgdf in res.groupby(["mutagen", "genotype", "dose"]):
    if mgdf["experiment_title"].nunique() > 1:
        tech_comp.append(mgdf)
tech_comp = pd.concat(tech_comp)

# print (tech_comp.groupby(["mutagen", "genotype", "dose", "sample", "experiment_title"]).size())

# test if seq tech has impact
seq_tech_totals = tech_comp.groupby(["experiment_title"]).agg(t=("is_multiallelic", "sum"), tt=("is_multiallelic", lambda m: len(m))).reset_index()
seq_tech_totals["perm"] = 0
seq_tech_totals["emp"] = 1
seq_tech_totals["frac"] = seq_tech_totals["t"] / seq_tech_totals["tt"]

perms = [seq_tech_totals]

for i in range(1, 101):
    tech_comp["experiment_title"] = np.random.permutation(tech_comp["experiment_title"].values)
    seq_tech_totals = tech_comp.groupby(["experiment_title"]).agg(t=("is_multiallelic", "sum"), tt=("is_multiallelic", lambda m: len(m))).reset_index()
    seq_tech_totals["perm"] = i
    seq_tech_totals["emp"] = 0
    seq_tech_totals["frac"] = seq_tech_totals["t"] / seq_tech_totals["tt"]
    perms.append(seq_tech_totals)

perms = pd.concat(perms)
print (perms)
f, ax = plt.subplots()
sns.stripplot(data=perms, x="emp", y="frac", hue="experiment_title", dodge=True, alpha=0.1, ax=ax)
f.savefig("tech.png")

for kind, kind_df in perms.groupby("experiment_title"):
    emp = kind_df[kind_df["emp"] == 1]["frac"].values[0]
    p = kind_df[kind_df["emp"] == 0]["frac"].values
    print (kind, np.sum(p >= emp), np.sum(p <= emp))