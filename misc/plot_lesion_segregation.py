import pandas as pd
from cyvcf2 import VCF
import matplotlib.pyplot as plt
import tqdm
import numpy as np
import seaborn as sns
import glob
import ruptures as rpt
import os

MUT2EXP = {
    "DMS": "T",
    "EMS": "G",
    "MMS": "A",
    "UV": "C",
    "ENU": "G",
}

REVCOMP = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
}


SYMMETRY_THRESH = 0.3

MUTAGEN = "EMS"

mutations = pd.read_csv("/scratch/ucgd/lustre-labs/quinlan/u1006375/lesion-segregation/csv/mutations/combined.tsv", sep="\t")

metadata = pd.read_excel(
    "/scratch/ucgd/lustre-labs/quinlan/u1006375/lesion-segregation/data/41467_2020_15912_MOESM4_ESM.xlsx",
    sheet_name="Sample_description",
    skiprows=1,
).fillna(
    {
        "Normal": "CD0001b",
        "Mutagen": "none",
        "Drug.concentration": -1,
    }
)

metadata["Mutagen"] = metadata.apply(
    lambda row: row["Mutagen"] if row["Type"] != "mut.acc." else "None (M.A.)", axis=1
)

mutations = mutations.merge(metadata, left_on="sample", right_on="Sample")

mutations = mutations[mutations["sample"] == "CD0566d"]

# mutations = mutations[(mutations["Mutagen"] == MUTAGEN) & (mutations["Genotype"] == "N2")]

counts = mutations.groupby("sample").size().reset_index().rename(columns={0: "count"})

top_n = counts.sort_values("count", ascending=False).head(50)["sample"].to_list()

mutations = mutations[mutations["sample"].isin(top_n)]

# figure out the expected "watson" and "crick" strands.
exp_watson = MUT2EXP[MUTAGEN]
exp_crick = REVCOMP[exp_watson]

mutations["strand"] = np.select(
    [
        mutations["ref"].eq(exp_watson),
        mutations["ref"].eq(exp_crick),
    ],
    [0, 1],
    default=-1,
).astype(np.int8)

# mutations["Genotype"] = "N2"
# mutations["Mutagen"] = mutations["sample"].apply(lambda s: s.split("_")[0])
# mutations["Drug.concentration"] = 1.

chrom2idx = dict(zip(["I", "II", "III", "IV", "V", "X"], range(6)))

for (sample, genotype, mutagen, dose), sdf in mutations.groupby(
    [
        "sample",
        "Genotype",
        "Mutagen",
        "Drug.concentration",
    ]
):
    if mutagen == "Control": continue
    print (sample, genotype, mutagen, dose)

   
    breakpoints = []

    all_chroms_have_data = sdf.groupby("chrom").size().reset_index().rename(columns={0: "count"})
    all_chroms_have_data["suff"] = all_chroms_have_data["count"] >= 1
    if not all_chroms_have_data["suff"].sum() == len(chrom2idx): continue

    f, axarr = plt.subplots(
        3, len(chrom2idx), sharex=True, sharey="row", figsize=(18, 5)
    )

    for ci, (chrom, chrom_df) in enumerate(sdf.groupby("chrom")):

        chrom_df = chrom_df.sort_values("start")

        pos = chrom_df["start"].values
        # dist = chrom_df["dist"].values
        dist = np.array([0] + list(pos[1:] - pos[:-1]))
        strands = chrom_df["strand"].values

        if strands.shape[0] > 1:

            algo = rpt.Pelt(model="rbf").fit(strands)
            cur = 1
            for bkpt in algo.predict(pen=np.log10(strands.shape[0])):
                # calculate symmetry score in region
                region_arr = strands[cur:bkpt]
                n_w, n_c = np.sum(region_arr == 0), np.sum(region_arr == 1)
                symmetry = (n_w - n_c) / (n_w + n_c)
                color = (
                    "darkgoldenrod"
                    if symmetry > SYMMETRY_THRESH
                    else "royalblue" if symmetry < -SYMMETRY_THRESH else "gainsboro"
                )

                n_inf_mutations = region_arr.shape[0]
                if n_inf_mutations < 20: continue

                breakpoints.append(
                    {
                        "chrom": chrom,
                        "start": pos[cur - 1],
                        "end": pos[bkpt - 1],
                        "symmetry": symmetry,
                        "n_inf_mutations": region_arr.shape[0],
                        "strand": (
                                0
                                if symmetry > SYMMETRY_THRESH
                                else 1 if symmetry < -SYMMETRY_THRESH else -1
                            ),
                    }
                )

                axarr[2, ci].plot(
                    [pos[cur - 1], pos[bkpt - 1]],
                    [symmetry, symmetry],
                    c=color,
                    lw=5,
                )
                cur = bkpt

        if chrom_df.shape[0] == 0: continue

        for yi, c in zip([0, 1], ["darkgoldenrod", "royalblue"]):
            idxs = np.where(strands == yi)
            label = (
                exp_watson + r"$\to$" + "N"
                if yi == 0
                else exp_crick + r"$\to$" + "N"
            )
            axarr[yi, ci].scatter(
                pos[idxs],
                dist[idxs],
                c=c,
                label=label,
                alpha=0.5,
            )
        axarr[0, ci].set_title(chrom)
        axarr[-1, ci].set_xlabel("Position (Mbp)")

    for chrom, ci in chrom2idx.items():
        axarr[1, ci].set_xticks([])
        axarr[1, ci].yaxis.set_inverted(True)
        for yi in [0, 1]:
            axarr[yi, ci].set_yscale("log")
        sns.despine(ax=axarr[2, ci])
        axarr[2, ci].set_ylim(-1.05, 1.05)

    for yi in [0, 1]:
        axarr[yi, 0].set_ylabel("Inter-SNP distance (bp)", size=8)
        axarr[yi, 0].legend()
        for ci in range(len(chrom2idx)):
            sns.despine(ax=axarr[yi, ci], top=True, right=True)
            sns.despine(ax=axarr[0, ci], bottom=True)
    axarr[2, 0].set_ylabel("Mutational assymmetry", size=8)
    f.tight_layout()
    f.savefig(f"fig/lesion_segregation/{mutagen}.{genotype}.{dose}.{sample}.png", dpi=200)
    plt.close()

    breakpoints = pd.DataFrame(breakpoints)
    breakpoints["sample"] = sample
    breakpoints.to_csv(
        f"csv/breakpoints/{mutagen}.{genotype}.{dose}.{sample}.breakpoints.tsv", sep="\t", index=False
    )
