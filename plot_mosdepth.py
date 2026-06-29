import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import glob
import numpy as np
import tqdm

plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans"]


mav_counts = pd.read_csv("multiallelic.tsv", sep="\t")

samples = mav_counts["sample"].unique()

colors = sns.color_palette("colorblind", 3)

res = []
for fh in tqdm.tqdm(glob.glob("data/mosdepth/*.summary.txt")):
    smp = fh.split("/")[-1].split(".")[0]
    df = pd.read_csv(fh, sep="\t")
    df = df[df["chrom"] == "total"]
    df["sample"] = smp
    res.append(df)
res = pd.concat(res)

res = res.merge(mav_counts)
res["frac"] = res["count"] / res["total"]

manifest = pd.read_csv("data/combined_file_manifest.tsv", sep="\t").rename(
    columns={"Sample / DONOR ID": "sample", "experiment_title": "Platform"}
)
manifest["Platform"] = manifest["Platform"].apply(lambda t: " ".join(t.split(" ")[:3]))

res = res.merge(manifest)

f, ax = plt.subplots()
sns.stripplot(
    data=res,
    x="Mutagen",
    y="mean",
    hue="Platform",
    dodge=True,
    palette="colorblind",
    alpha=0.5,
    ax=ax,
)
sns.despine(ax=ax)
ax.set_ylabel("Average genome-wide sequencing depth")
f.tight_layout()
f.savefig("mean.png", dpi=200)
