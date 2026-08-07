import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import tqdm

from snakemake.script import snakemake

plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans"]

mutations = pd.read_csv(snakemake.input.mutations, sep="\t")

samples = mutations["sample"].unique()

colors = sns.color_palette("colorblind", 3)

res = []
for fh in tqdm.tqdm(snakemake.input.summaries):
    smp = fh.split("/")[-1].split(".")[0]
    df = pd.read_csv(fh, sep="\t")
    df = df[df["chrom"] == "total"]
    df["sample"] = smp
    res.append(df)
res = pd.concat(res)

manifest = pd.read_csv(snakemake.input.manifest, sep="\t").rename(
    columns={"Sample / DONOR ID": "sample", "experiment_title": "Platform"}
)
manifest["Platform"] = manifest["Platform"].apply(lambda t: " ".join(t.split(" ")[:3]))
metadata = pd.read_excel(
    snakemake.input.metadata,
    sheet_name="Sample_description",
    skiprows=1,
)

res = res.merge(manifest)
res = res.merge(metadata, left_on="sample", right_on="Sample")

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
f.savefig(snakemake.output.mosdepth, dpi=200)
