import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import numpy as np

metadata = pd.read_excel(
    "data/41467_2020_15912_MOESM4_ESM.xlsx",
    sheet_name="Sample_description",
    skiprows=1,
)

res = []
for fh in glob.glob("csv/multiallelics/*.tsv"):
    df = pd.read_csv(fh, sep="\t")
    res.append(df)
res = pd.concat(res).fillna("N:0")


res = res.merge(metadata, left_on="sample", right_on="Sample")

counts = (
    res.groupby(
        [
            "sample",
            "Mutagen",
            "Drug.concentration",
            "Genotype",
            "permuted",
            "is_multiallelic",
        ]
    )
    .size()
    .reset_index()
    .rename(columns={0: "count"})
)
totals = counts.groupby(["sample", "permuted"]).agg(total=("count", "sum")).reset_index()
counts = counts.merge(totals)
counts["frac"] = counts["count"] / counts["total"]

counts = counts[counts["is_multiallelic"] == False]

counts["true_frac"] = 1 - counts["frac"]

res = []
f, ax = plt.subplots()
for s, sdf in counts.groupby("sample"):
    perm = sdf[sdf["permuted"] > 0]["true_frac"].values
    mean_perm = np.mean(perm)
    emp = sdf[sdf["permuted"] == 0]["true_frac"].values
    assert emp.shape[0] == 1
    emp = emp[0]
    ax.plot([0, 1], [emp, mean_perm])
    for kind, val in zip(("emp", "perm"), (emp, mean_perm)):
        res.append({"sample": s, "status": kind, "val": val})
f.savefig('test.png')
res = pd.DataFrame(res)
f, ax = plt.subplots()
sns.stripplot(data=res, x="status", y="val", hue="status", ax=ax)
ax.set_yscale('log')
ax.set_xticks([0, 1])
ax.set_xticklabels(["Empirical", "Permuted"])
ax.set_xlabel("")
ax.set_ylabel("Fraction of sites with multi-allelic variation")
sns.despine(ax=ax)
f.savefig("multi_perm.png", dpi=200)
