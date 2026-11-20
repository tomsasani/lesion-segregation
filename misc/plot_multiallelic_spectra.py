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
import scipy.stats as ss


plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans"]


THRESHOLD = 3
N_PERM = 100
res = []
for fh in tqdm.tqdm(glob.glob(f"csv/multiallelics/*.{N_PERM}.tsv")):
    try:
        df = pd.read_csv(fh, sep="\t").fillna("N:0")
    except pd.errors.EmptyDataError:
        continue
    df = df.query("permuted == 0")
    df["is_multiallelic"] = df["ma_reads"] >= THRESHOLD

    res.append(df)

res = pd.concat(res).fillna("N:0")

res = res[res["mutagen"].isin(["EMS", "DMS", "MMS"])]

# remove sites with no true mutation in the reads
res["n_bases"] = res["summary"].apply(lambda s: s.count(":"))
res = res.query("n_bases >= 2")

print (res[res["is_multiallelic"]].groupby("genotype").size())

revcomp = {"C": "G", "G": "C", "A": "T", "T": "A"}

muts = []
Ts = 10

def get_mutant_allele(row: pd.Series):
    ref = row["ref"]
    exp_alt = row["alt"]
    oth_alt = None
    if row["is_multiallelic"]:
        # get other bases
        summary = row["summary"].split("|")
        for bc in summary:
            b, c = bc.split(":")
            if b not in (ref, exp_alt):
                oth_alt = b

    if ref not in ("C", "A"):
        ref, exp_alt = revcomp[ref], revcomp[exp_alt]
        if oth_alt is not None:
            oth_alt = revcomp[oth_alt]
    if oth_alt is None:
        return ref + r"$\to$" + exp_alt
    else:
        return ref + r"$\to$" + oth_alt


for t in tqdm.tqdm(range(Ts)):

    for mutagen, mdf in res.sample(frac=1, replace=True).groupby("mutagen"):

        mdf["mut"] = mdf.apply(
            lambda row: get_mutant_allele(row),
            axis=1)

        spectra = (
            mdf.groupby(
                [
                    "is_multiallelic",
                    "mut",
                ]
            )
            .size()
            .reset_index()
            .rename(columns={0: "count"})
        )

        totals = mdf.groupby("is_multiallelic").size().reset_index().rename(columns={0: "total"})
        spectra = spectra.merge(totals)
        spectra["frac"] = spectra["count"] / spectra["total"]
        spectra["mutagen"] = mutagen
        spectra["t"] = t
        muts.append(spectra)

muts = pd.concat(muts)
wide = muts.pivot(index=["mutagen", "is_multiallelic", "t"], columns="mut", values="frac").fillna(0).reset_index()

# make sure we have an entry for every combo
mux = pd.MultiIndex.from_product([res["mutagen"].unique(), [True, False], list(range(Ts))], names=["mutagen", "is_multiallelic", "t"])
wide = wide.set_index(["mutagen", "is_multiallelic", "t"]).reindex(mux, fill_value=0).reset_index()

tidy = wide.melt(id_vars=["mutagen", "is_multiallelic", "t"])

g = sns.FacetGrid(data=tidy, col="mutagen", aspect=1.5, col_wrap=2)
g.map(
    sns.barplot,
    "mut",
    "value",
    "is_multiallelic",
    order=sorted(tidy["mut"].unique(), key=lambda m: m[0] + m[-1]),
    palette="colorblind",
)
g.add_legend()
g.savefig("spec.png", dpi=200)


# vals = np.concatenate((muts[muts["is_multiallelic"]]["a1"].values, muts[muts["is_multiallelic"]]["a2"].values))
# vals = muts[muts["is_multiallelic"]]["a1"].values

# def gmm_bic_score(estimator, X):
#     """Callable to pass to GridSearchCV that will use the BIC score."""
#     # Make it negative since GridSearchCV expects a score to maximize
#     return -estimator.bic(X)

# param_grid = {
#     "n_components": range(1, 5),
#     "covariance_type": ["spherical", "tied", "diag", "full"],
#     "init_params": ["random_from_data", "kmeans"]
# }
# grid_search = GridSearchCV(
#     GaussianMixture(),
#     param_grid=param_grid,
#     scoring=gmm_bic_score,
#     cv=3,
#     refit=True,
# )
# grid_search.fit(vals.reshape(-1, 1))
# print (grid_search.best_params_)
# estimator = grid_search.best_estimator_

# weights = estimator.weights_
# means = estimator.means_.flatten()
# stds = np.sqrt(estimator.covariances_.flatten())
# if stds.shape[0] == 1:
#     stds = np.repeat(stds, 2)
# print (weights.shape, means.shape, stds.shape)
# # 3. Generate x-values and 4. Calculate PDFs
# x = np.linspace(0, 1, 1000)
# y_components = [w * ss.norm.pdf(x, m, s) for w, m, s in zip(weights, means, stds)]
# y_gmm_total = np.sum(y_components, axis=0)
# print (len(y_components))

# bins = np.arange(0, 1, 0.02)
# f, ax = plt.subplots()
# hist, edges = np.histogram(vals, bins=bins)
# ax.bar(edges[:-1], hist, edges[1] - edges[0], ec="w", lw=1)
# for i, comp in enumerate(y_components):
#     ax.plot(x, comp, '--', label=f'Component {i+1}')
# ax.legend()
# # ax.hist(grid_search.predict(vals.reshape(-1,1)), bins=bins, histtype="step")
# # sns.stripplot(data=muts[muts["is_multiallelic"]], y="a1", ax=ax)
# f.savefig("ratio.png", dpi=200)
