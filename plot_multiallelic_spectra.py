import pandas as pd
import seaborn as sns
from global_vars import REVCOMP
import os
import glob
from mutyper import Ancestor
import tqdm
import scipy.stats as ss
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import logomaker
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import GridSearchCV

def calc_entropy(r):
    counts = Counter(r.upper())
    probabilities = np.array(list(counts.values())) / len(r)
    shannon_entropy = ss.entropy(probabilities, base=2)
    return shannon_entropy


K = 11
ancestor = Ancestor(
    "data/ref/ce10.renamed.fa",
    sequence_always_upper=True,
    k=K,
)

res = []
for fh in glob.glob("csv/multiallelics/*.0.tsv"):
    if os.path.getsize(fh) == 1: continue
    df = pd.read_csv(fh, sep="\t")
    res.append(df)
res = pd.concat(res).dropna(subset="summary")
# res = res[~res["is_strand_biased"]]
# remove sites with no true mutation in the reads
res["n_bases"] = res["summary"].apply(lambda s: s.count(":"))
res = res.query("n_bases >= 2")


## spectra
muts = []
for i, row in res.iterrows():
    ref = row["ref"]

    total = 0
    non_ref_counts = []
    for b in row["summary"].split("|"):
        alt, count = b.split(":")
        count = int(count)
        if alt == ref: 
            # total += count
            continue
        non_ref_counts.append(count)
        total += count
    non_ref_fracs = [c / total for c in non_ref_counts]
    d = row.to_dict()  
    d.update(
        {
            "a1": non_ref_fracs[0],
            "a2": non_ref_fracs[1] if row["is_multiallelic"] else 0.,
        }
    )
    muts.append(d)

muts = pd.DataFrame(muts)
print (muts[muts["is_multiallelic"]])

vals = np.concatenate((muts[muts["is_multiallelic"]]["a1"].values, muts[muts["is_multiallelic"]]["a2"].values))
vals = muts[muts["is_multiallelic"]]["a1"].values

def gmm_bic_score(estimator, X):
    """Callable to pass to GridSearchCV that will use the BIC score."""
    # Make it negative since GridSearchCV expects a score to maximize
    return -estimator.bic(X)

param_grid = {
    "n_components": range(1, 5),
    "covariance_type": ["spherical", "tied", "diag", "full"],
    "init_params": ["random_from_data", "kmeans"]
}
grid_search = GridSearchCV(
    GaussianMixture(), param_grid=param_grid, scoring=gmm_bic_score
)
grid_search.fit(vals.reshape(-1, 1))
print (vals.shape)
print (grid_search.best_params_)

f, ax = plt.subplots()
ax.hist(vals, bins=np.arange(0.5, 1, 0.025))
# sns.stripplot(data=muts[muts["is_multiallelic"]], y="a1", ax=ax)
f.savefig("ratio.png")



# NOTE: to compare the sequence entropy of multiallelics vs non, we should
# probably match on the representation of different mutagens/genotypes

counts = res.groupby(["is_multiallelic", "genotype", "mutagen"]).size().reset_index().rename(columns={0: "count"})
totals = counts.groupby("is_multiallelic").agg(total=("count", "sum")).reset_index()
counts = counts.merge(totals)


resampled = []

nma = res[~res["is_multiallelic"]]
ma = res[res["is_multiallelic"]]
ma["t"] = 0
N_PERM = 200
for t in range(N_PERM):

    for i, row in counts[counts["is_multiallelic"]].iterrows():
        genotype, mutagen, count = row["genotype"], row["mutagen"], row["count"]

        sample = nma[(nma["genotype"] == genotype) & (nma["mutagen"] == mutagen)].sample(n=count, replace=False)
        sample["t"] = t
        resampled.append(sample)

resampled = pd.concat(resampled)

res = pd.concat([resampled, ma])


muts = []
for i, row in tqdm.tqdm(res.iterrows()):
    ref = row["ref"]

    for b in row["summary"].split("|"):
        alt, count = b.split(":")
        if alt == ref: 
            continue

        mutation = ancestor.mutation_type(
            row["chrom"],
            row["start"],
            ref,
            alt,
        )
        refseq = list(mutation[0])      
        d = row.to_dict()  
        d.update(
            {
                "refseq": "".join(refseq),
            }
        )
        muts.append(d)
        break

muts = pd.DataFrame(muts)

muts["entropy"] = muts["refseq"].apply(lambda r: calc_entropy(r))

f, ax = plt.subplots()
sns.boxplot(
    data=muts[muts["is_multiallelic"]],
    x="mutagen",
    y="entropy",
    hue="genotype",
    dodge=True,
    ax=ax,
    fliersize=0,
)
f.savefig("en.png")

# f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
# seqs = muts[~muts["is_multiallelic"]]["refseq"].to_list()

# matrix = logomaker.transform_matrix(logomaker.alignment_to_matrix(seqs), from_type="counts", to_type="probability")
# logomaker.Logo(matrix, fade_probabilities=True,shade_below=.5,
#                           fade_below=.5,
#                          stack_order='small_on_top', ax=ax1)
# seqs = muts[muts["is_multiallelic"]]["refseq"].to_list()

# matrix = logomaker.transform_matrix(logomaker.alignment_to_matrix(seqs), from_type="counts", to_type="probability")
# logomaker.Logo(matrix, fade_probabilities=True,shade_below=.5,
#                           fade_below=.5,
#                          stack_order='small_on_top', ax=ax2)
# f.savefig("o.png")

nma = muts[~muts["is_multiallelic"]]
ma = muts[muts["is_multiallelic"]]


ma_entropy = ma["entropy"].values

stats = []
for t, tdf in nma.groupby("t"):
    nma_entropy = tdf["entropy"].values

    stats.append(nma_entropy.mean())
f, ax = plt.subplots(figsize=(8, 6))
ax.hist(stats, label="Randomly sampled SNVs")
ax.axvline(ma["entropy"].mean(), label="Multiallelic SNVs", ls=":", c="k")
ax.legend(shadow=True)
ax.set_xlabel("Average sequence entropy")
ax.set_ylabel("Count")
p = np.sum(np.array(stats) <= ma["entropy"].mean()) / N_PERM
ax.set_title(f"Multiallelic SNVs have lower average sequence entropy\nthan randomly sampled SNVs (matched on genotype and mutagen)\n{p}")
sns.despine(ax=ax)
f.savefig('entropy.png', dpi=200)

# counts = muts.groupby(["mutagen", "genotype", "mutation"]).size().reset_index().rename(columns={0: "count"})

# g = sns.FacetGrid(data=counts, col="mutagen", row="genotype", sharex=True)
# g.map(sns.barplot, "mutation", "count")
# g.savefig("spec.png")
