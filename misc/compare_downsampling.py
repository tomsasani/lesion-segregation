import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans"]

dp_orig = []
for fh in glob.glob("csv/mosdepth/*.mosdepth.summary.txt"):
    df = pd.read_csv(fh, sep="\t")
    smp = fh.split("/")[-1].split(".")[-4]
    df["sample"] = smp
    df["source"] = "original"

    df = df[df["chrom"] == "total"]
    dp_orig.append(df)
dp_orig = pd.concat(dp_orig)

dp_dwns = []
for fh in glob.glob("csv/mosdepth/downsampled/*.mosdepth.summary.txt"):
    df = pd.read_csv(fh, sep="\t")
    smp = fh.split("/")[-1].split(".")[-4]
    df["sample"] = smp
    df["source"] = "downsampled"
    df = df[df["chrom"] == "total"]

    dp_dwns.append(df)
dp_dwns = pd.concat(dp_dwns)
depths = pd.concat([dp_dwns, dp_orig])
# depths = dp_orig.merge(dp_dwns, on=["sample", "chrom", "length"])[["sample", "mean_x", "mean_y", "source_x", "source_y"]]

N_PERM = 0

res = []
for fh in glob.glob(f"csv/multiallelics/*.{N_PERM}.*.tsv"):
    try:
        df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    except (pd.errors.EmptyDataError, FileNotFoundError):
        continue

    df["source"] = fh.split("/")[-1].split(".")[-2]
    df["is_multiallelic"] = df["ma_reads"] >= 3
    res.append(df)
res = pd.concat(res)

print(res[(res["sample"] == "CD0495d") & (res["is_multiallelic"])])


grouped = (
    res.groupby(["sample", "source"])
    .agg(
        count=("is_multiallelic", "sum"),
        total=("is_multiallelic", len),
    )
    .reset_index()
)
grouped["frac"] = grouped["count"] / grouped["total"]
# grouped["source"] = grouped["source"].apply(lambda s: "downsampled" if s == "downsample" else "original")
merged = grouped.merge(depths)

wide = pd.pivot(data=merged, index="sample", columns=["source"], values=["count", "mean", "frac"]).reset_index()
wide.columns = ["sample", "dwns_count", "orig_count", "dwns_dp", "orig_dp", "dwns_frac", "orig_frac"]

wide["count_diff"] = wide["orig_count"] - wide["dwns_count"]
wide["frac_diff"] = wide["orig_frac"] - wide["dwns_frac"]

print (wide["orig_count"].sum(), wide["dwns_count"].sum())
print (wide[wide["dwns_count"] > wide["orig_count"]])
f, ax = plt.subplots()
wide = wide[wide["orig_frac"] != wide["dwns_frac"]]
ax.plot([wide["dwns_dp"].values, wide["orig_dp"].values], [wide["dwns_frac"].values, wide["orig_frac"]])
f.savefig("t.png")
print (grouped)
print (grouped.groupby("strategy").agg(mav=("count", "sum")))
pivoted = pd.pivot(grouped, index="sample", columns=["strategy"], values=["count"]).reset_index()
pivoted.columns = ["sample", "downsample", "regular"]
pivoted = pivoted[pivoted["regular"] > 0]

f, ax = plt.subplots()
ax.axline((0, 0), slope=1, ls=":", c="k", zorder=-1)

sns.scatterplot(data=pivoted, x="regular", y="downsample", ax=ax, alpha=0.5, lw=1, ec="w", zorder=1)
sns.despine(ax=ax)

ax.set_xlabel("MAV fraction (calculated using original CRAM)")
ax.set_ylabel("MAV fraction (calculated using downsampled CRAM)")
f.tight_layout()
f.savefig("ds.png", dpi=200)

rng = np.random.default_rng(42)
f, ax = plt.subplots()
pivoted = pivoted[pivoted["regular"] > 0]
print (pivoted)
x1 = rng.normal(0, 0.01, size=pivoted.shape[0])
x2 = rng.normal(1, 0.01, size=pivoted.shape[0])
y1 = pivoted["regular"].values
y2 = pivoted["downsample"].values
ax.scatter(x1, y1)
ax.scatter(x2, y2)
diffs = pivoted["downsample"].values / pivoted["regular"].values
print (diffs)
# diffs[np.isnan(diffs)] = 0
diffs[diffs > 1] = 1
alphas = 1 - diffs
# alphas += 0.25
print (alphas)
for i in range(x1.shape[0]):
    ax.plot([x1[i], x2[i]], [y1[i], y2[i]], alpha=alphas[i], c="gainsboro")
# ax.plot(
#     [rng.normal(0, 0.01, size=pivoted.shape[0]), rng.normal(1, 0.01, size=pivoted.shape[0])],
#     [pivoted["regular"].values, pivoted["downsample"].values],
#     alpha=alphas,
# )
f.savefig("dss.png")
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
