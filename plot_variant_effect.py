from cyvcf2 import VCF
import numpy as np
import pandas as pd
import csv
import tabix
import matplotlib.pyplot as plt
import tqdm
import scipy.stats as ss
import seaborn as sns


tb = tabix.open("data/ce11.phyloP135way.bedg.gz")

# These queries are identical. A query returns an iterator over the results.

# res = []

# with open("csv/mutations/combined.tsv", "rt") as infh:
#     csvf = csv.reader(infh, delimiter="\t")
#     header = None
#     for l in tqdm.tqdm(csvf):
#         if header is None:
#             header = l
#             continue
#         d = dict(zip(header, l))
#         records = tb.query(f"chr{d['chrom']}", int(d["start"]), int(d["end"]))
#         r = [_ for _ in records]
#         assert len(r) <= 1
#         p = np.nan
#         if len(r) == 1:
#             p = float(r[0][-1])
#         d.update({"phylop": p})
#         res.append(d)
# res = pd.DataFrame(res)
# res.to_csv("phylop.tsv", sep="\t", index=False)
res = pd.read_csv("phylop.tsv", sep="\t")
res = res[res["mutagen"].isin(["EMS", "DMS", "MMS"])]
_min, _max = res["phylop"].min(), res["phylop"].max()
f, ax = plt.subplots()
bins = np.arange(_min, _max, 0.5)
for mutagen, mutagen_df in res.dropna().groupby("mutagen"):
    scores = mutagen_df["phylop"]
    if scores.shape[0] < 1_000: continue
    bootstrap = np.empty((100, bins.shape[0] - 1))
    for bs in range(100):
        sampled = np.random.choice(scores, replace=True, size=scores.shape[0])
        hist, edges = np.histogram(sampled, bins=bins)
        hist_frac = hist / np.sum(hist)
        bootstrap[bs] = np.cumsum(hist_frac)
    u = np.mean(bootstrap, axis=0)

    # n50
    hist, edges = np.histogram(scores, bins=np.arange(_min, _max, 0.01))
    hist_frac = hist / np.sum(hist)

    lo = np.percentile(bootstrap, q=2.5, axis=0)
    hi = np.percentile(bootstrap, q=97.5, axis=0)
    ax.plot(bins[:-1], u, label=mutagen)
    ax.fill_between(bins[:-1], lo, hi)
ax.legend(title="Mutagen", shadow=True)
sns.despine(ax=ax)
ax.set_xlabel("phyloP score (low " + r"$\to$" + " high conservation)")
ax.set_ylabel("Cumulative fraction of SNVs")
# ax.set_xscale("symlog")
f.savefig("phylop.png", dpi=200)
