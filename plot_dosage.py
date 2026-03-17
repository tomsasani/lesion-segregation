import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from global_vars import generate_collapsed_3mers
import numpy as np
from scipy.stats import chi2_contingency


PLOT_IND = False


res = pd.read_csv("mutations.tsv", sep="\t")

all_kmers = generate_collapsed_3mers()
sorted_kmers = sorted(all_kmers, key=lambda k: k[1] + k[5])

base_nucs = [k[1] + k[5] for k in sorted_kmers]
B2C = {
    "CA": "blue",
    "CG": "black",
    "CT": "red",
    "AC": "grey",
    "AG": "green",
    "AT": "pink",
}
colors = [B2C[k] for k in base_nucs]


def calculate_slope(gdf,  plot: bool = False):

    if plot:
        f2, axarr = plt.subplots(
            gdf["dose"].nunique(),
            2,
            figsize=(12, 12),
            sharex=True,
            sharey="col",
        )

    doses = []
    _a = []
    _b = []

    wide = gdf.pivot(index=["dose", "EMS_sig"], columns="kmer", values="avg").fillna(0)
    kmers = list(wide.columns)
    wide = wide.reset_index()
    new_cols = ["dose", "EMS_sig"] + kmers

    wide.columns = new_cols
    for k in all_kmers:
        if k not in wide.columns:
            wide[k] = 0.

    for di, (dose, df) in enumerate(wide.groupby("dose")):
        try:
            a = df[df["EMS_sig"] == True][sorted_kmers].values[0]
            b = df[df["EMS_sig"] == False][sorted_kmers].values[0]
        except IndexError:
            continue

        doses.append(dose)
        _a.append(np.sum(a))
        _b.append(np.sum(b))

        if plot:
            axarr[di, 0].bar(range(len(sorted_kmers)), a, 1, color=colors)
            axarr[di, 1].bar(range(len(sorted_kmers)), b, 1, color=colors)
            axarr[di, 0].set_title(dose)
            axarr[di, 0].set_title(round(np.sum(a)))
            axarr[di, 1].set_title(round(np.sum(b)))
            axarr[di, 0].set_xlabel(f"Spectrum at {dose}")
            for i in (0, 1):
                sns.despine(ax=axarr[di, i])

    if plot:
        f2.tight_layout()
        f2.savefig(f'fig/{g}.png', dpi=200)
        plt.close()

    # keep track of slopes of these lines
    slope = (_b[-1] - _b[0]) / (doses[-1] - doses[0])
    
    return slope, doses, _a, _b


f1, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=False, figsize=(10, 4))

genotype2diff = {}

# group by genotype, mutagen, and dose, and sum over replicates
grouped = (
        (
            res.groupby(
                ["sample", "genotype", "mutagen", "dose", "kmer"]
            )
            .size()
            .reset_index()
            .rename(columns={0: "count"})
        )
        .groupby(["genotype", "mutagen", "dose", "kmer"])
        .agg(
            total=("count", "sum"),
        )
        .reset_index()
    )

# for each combination of mutagen and genotype, we want to compare
# spectra against N2 with that same mutagen

for (genotype, mutagen), gdf in grouped.groupby(["genotype", "mutagen"]):
    if genotype == "N2": continue

    # get N2 for this mutagen
    n2_with_mutagen = grouped[
        (grouped["genotype"] == "N2") & (grouped["mutagen"] == mutagen)
    ]

    n_doses = 1#gdf["dose"].nunique()
    f, axarr = plt.subplots(2)

    # for di, (dose, dose_df) in enumerate(gdf.groupby("dose")):

        # n2_with_dose = n2_with_mutagen[n2_with_mutagen["dose"] == dose]
        # if n2_with_dose.shape[0] == 0: continue

    combined = pd.concat([n2_with_mutagen, gdf])
    wide = combined.pivot(
        index="genotype",
        columns="kmer",
        values="total",
    ).fillna(0)

    kmers = list(wide.columns)
    wide = wide.reset_index()
    new_cols = ["genotype"] + kmers
    wide.columns = new_cols
    for k in all_kmers:
        if k not in wide.columns:
            wide[k] = 0

    # zero out non-CT in EMS
    # for k in sorted_kmers:
    #     if mutagen == "EMS": 
    #         if k[1] == "C" and k[5] == "T":
    #             wide[k] = 0

    n2_spec = wide[sorted_kmers].values[0, :]
    gt_spec = wide[sorted_kmers].values[1, :]

    # idx_with_more = int(np.sum(gt_spec) < np.sum(n2_spec))

    # fc = np.zeros(n2_spec.shape[0])
    # for ki, k in enumerate(sorted_kmers):
    #     a_fore = wide[k].values[idx_with_more] + 1
    #     a_back = np.sum(wide[sorted_kmers].values[idx_with_more, :]) - a_fore
    #     b_fore = wide[k].values[1 - idx_with_more] + 1
    #     b_back = np.sum(wide[sorted_kmers].values[1 - idx_with_more, :]) - b_fore
    #     chi2, p, dof, expected = chi2_contingency([[a_fore, b_fore], [a_back, b_back]])
    #     fc[ki] = chi2

    # if n_doses == 1:
    axarr[0].bar(range(len(sorted_kmers)), n2_spec, 1, color=colors)
    axarr[1].bar(range(len(sorted_kmers)), gt_spec, 1, color=colors)
    # else:
    #     axarr[di, 0].bar(range(len(sorted_kmers)), n2_spec, 1, color=colors)
    #     axarr[di, 1].bar(range(len(sorted_kmers)), gt_spec, 1, color=colors)
    # if n_doses == 1:
    axarr[0].set_title(f"N2")
    axarr[1].set_title(f"{genotype}")
    # else:
    #     axarr[di, 0].set_title(f"N2 {dose}")
    #     axarr[di, 1].set_title(f"{genotype}")

    f.tight_layout()
    f.savefig(f"fig/{genotype}.{mutagen}.png", dpi=200)
