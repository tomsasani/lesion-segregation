from cyvcf2 import VCF
import numpy as np
import pandas as pd
import csv
import tabix
import matplotlib.pyplot as plt
import tqdm
import scipy.stats as ss
import seaborn as sns
import itertools

from snakemake.script import snakemake

tb = tabix.open(snakemake.input.phylop_scores)

res = []

with open(snakemake.input.mutations, "rt") as infh:
    csvf = csv.reader(infh, delimiter="\t")
    header = None
    for l in tqdm.tqdm(csvf):
        if header is None:
            header = l
            continue
        d = dict(zip(header, l))
        records = tb.query(f"chr{d['chrom']}", int(d["start"]), int(d["end"]))
        r = [_ for _ in records]
        assert len(r) <= 1
        p = np.nan
        if len(r) == 1:
            p = float(r[0][-1])
        d.update({"phylop": p})
        res.append(d)

res = pd.DataFrame(res)

res.to_csv(snakemake.output.tsv, sep="\t", index=False)

res = res[res["mutagen"].isin(["EMS", "DMS", "MMS"])]

_min, _max = res["phylop"].min(), res["phylop"].max()
f, ax = plt.subplots()
bins = np.arange(_min, _max, 0.01)

out = []

for mutagen, mutagen_df in res.groupby("mutagen"):
    scores = mutagen_df["phylop"]
    tot = scores.shape[0]

    scores = scores[~np.isnan(scores)]
    hist, edges = np.histogram(scores, bins=bins)
    hist_frac = hist / np.sum(hist)
    cumsum = np.cumsum(hist_frac)
    ii = np.where(cumsum >= 0.9)[0].min()
    print (mutagen, edges[ii])

    out.append((mutagen, scores.to_list()))
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
    ax.plot(bins[:-1], u, label=f"{mutagen} ({tot:,} SNVs)")
    ax.fill_between(bins[:-1], lo, hi)
ax.legend(title="Mutagen", shadow=True)
sns.despine(ax=ax)
ax.set_xlabel("phyloP score (low " + r"$\to$" + " high conservation)")
ax.set_ylabel("Cumulative fraction of SNVs")
f.tight_layout()
f.savefig(snakemake.output.phylop, dpi=200)


for i, j in itertools.combinations(range(3), r=2):
    print (out[i][0], out[j][0], ss.ks_2samp(out[i][1], out[j][1]).pvalue)
