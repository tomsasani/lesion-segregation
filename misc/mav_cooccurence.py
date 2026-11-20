import pandas as pd
import itertools
import numpy as np
import tqdm
import matplotlib.pyplot as plt
import seaborn as sns


plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans"]


mavs = pd.read_csv("https://raw.githubusercontent.com/tomsasani/lesion-segregation/refs/heads/main/csv/multiallelics.filtered.tsv", sep="\t")

N = 10

n_co_emp = 0
n_pairs_emp = 0
for sample, sample_df in mavs.groupby("sample"):
    if sample_df.shape[0] < 2: continue
    for chrom_a, chrom_b in itertools.combinations(sample_df["chrom"].to_list(), 2):
        n_pairs_emp += 1
        n_co_emp += int(chrom_a == chrom_b)

mutations = pd.read_csv("csv/mutations/combined.tsv", sep="\t")
res = []
rng = np.random.default_rng(42)
for t in tqdm.tqdm(range(N)):
    n_co = 0
    n_pairs = 0
    for sample, sample_df in mavs.groupby("sample"):
        n_mavs = sample_df.shape[0]
        if n_mavs < 2: continue
        sample_muts = mutations[mutations["sample"] == sample]
        # randomly sample n "mavs"
        _mavs = sample_muts.sample(n=n_mavs, replace=False, random_state=rng)
        for chrom_a, chrom_b in itertools.combinations(_mavs["chrom"].to_list(), 2):
            n_pairs += 1
            n_co += int(chrom_a == chrom_b)

    res.append({"trial": t, "n_pairs": n_pairs, "n_co": n_co})
res = pd.DataFrame(res)

f, ax = plt.subplots(figsize=(7, 5))
ax.hist(res["n_co"], bins=np.arange(res["n_co"].min(), max([n_co_emp, res["n_co"].max()])), lw=1, ec="k", color="cornflowerblue", label="Randomly-assigned MAVs")
ax.axvline(n_co_emp, ls=":", c="k", label="Empirical MAVs")
ax.set_xlabel("Intra-strain MAV pairs on the same chromosome")
ax.set_ylabel("Count")
ax.legend(shadow=True)
sns.despine(ax=ax)
f.savefig("cooc.png", dpi=200)
print (np.sum(res["n_co"].values >= n_co_emp) / N)
print (res["n_co"].mean())


# another approach. one MAV-bearing chromosomes, AFs of nearby variants
mavs = pd.read_csv("csv/multiallelics.100.original.filtered.tsv", sep="\t")

res = []

for (sample, chrom), sdf in mavs.groupby(["sample", "chrom"]):
    rel_muts = mutations[(mutations["sample"] == sample) & ~(mutations["chrom"] == chrom)]
    rel_muts = rel_muts[["sample", "chrom", "af"]]
    rel_muts["is_mav_bearing"] = True
    res.append(rel_muts)

res = pd.concat(res)

f, ax = plt.subplots()
sns.histplot(data=res, x="af", ax=ax)
f.savefig("t.png")