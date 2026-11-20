import pandas as pd
from cyvcf2 import VCF
import seaborn as sns
import tqdm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import ruptures as rpt
import numpy as np
from collections import defaultdict
from scipy.stats import norm, binom

# mutagens = ["MMS", "EMS", "DMS", "UV"]

# genotypes = ["N2", "polk-1", "polh-1"]

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


# mutations["Genotype"] = "N2"
# mutations["Mutagen"] = mutations["sample"].apply(lambda s: s.split("_")[0])
# mutations["Drug.concentration"] = 1.

# mutations = mutations[mutations["Mutagen"].isin(mutagens)]
mutations = mutations[mutations["Genotype"].isin(["rcq-5"])]
# mutations = mutations[mutations["sample"].isin(["CD0623c", "CD0540c"])]
rng = np.random.default_rng(42)

for (genotype, mutagen), sdf in tqdm.tqdm(mutations.groupby(["Genotype", "Mutagen"])):

    sdf = sdf[sdf["Drug.concentration"] != 0]

    # print (genotype, mutagen)
    # print (sdf[["chrom", "start", "end", "af"]])

    n = sdf["sample"].nunique()
    nc = 6

    f, axarr = plt.subplots(n, nc, figsize=(16, 2.5 * n), sharex="col", sharey=True)
    # f2, ax2 = plt.subplots()

    res = []
    for si, ((dose, sample), ddf) in enumerate(sdf.groupby(["Drug.concentration", "sample"])):    

        B2C = {"A": "red", "T": "blue", "C": "green", "G": "gray"}

        block_means = defaultdict(list)
        block_sizes = []

        for ci, (chrom, chrom_df) in enumerate(ddf.groupby("chrom")):

            starts = chrom_df["start"].values
            balances = chrom_df["af"].values

            if chrom_df.shape[0] < 2: continue

            folded_balances = np.max(np.vstack((balances, 1 - balances)), axis=0)
            folded_balances = balances

            refs = chrom_df["ref"].to_list()
            refcols = [B2C[r] for r in refs]

            algo = rpt.Pelt(model="rbf").fit(folded_balances)
            cur = 1
    
            axarr[si, ci].scatter(
                starts,
                folded_balances,
                c=refcols,
                # label=ref,
                zorder=1,
            )

            axarr[si, ci].set_xlabel(f"Position on {chrom}")
            axarr[si, ci].set_ylabel("Variant allele freq.")
            # axarr[si, ci].legend()
            axarr[si, ci].set_title(f"{mutagen} @ {dose}\nin {sample}")
            sns.despine(ax=axarr[si, ci])

        
    f.tight_layout()
    f.savefig(f"fig/ab/{mutagen}.{genotype}.ab.png", dpi=200)
    plt.close()
