import pandas as pd
from cyvcf2 import VCF
import seaborn as sns
import tqdm
import matplotlib.pyplot as plt
from global_vars import REVCOMP

import numpy as np

metadata = pd.read_excel(
    "41467_2020_15912_MOESM4_ESM.xlsx",
    sheet_name="Sample_description",
    skiprows=1,
)

mutagens = ["MMS"]

genotypes = [
    "polk-1"
    
]

multi = pd.read_csv("multi.txt", sep="\t", names=["sample", "chrom", "pos", "ref", "alt", "base_counts", "ab"]).fillna("N:0")
multi["is_multiallelic"] = multi["base_counts"].apply(lambda b: len(b.split("|")) > 2 and min([int(_b.split(":")[1]) for _b in b.split("|")]) >= 2)
multi = multi[multi["is_multiallelic"] == True]

metadata = metadata[metadata["Mutagen"].isin(mutagens)]
metadata = metadata[metadata["Genotype"].isin(genotypes)]

for (genotype, mutagen), sdf in tqdm.tqdm(metadata.groupby(["Genotype", "Mutagen"])):

    sdf = sdf[sdf["Drug.concentration"] != 0]

    n = sdf["Sample"].nunique()
    nc = 6

    # doses = sorted(sdf["Drug.concentration"].unique())
    # dose2idx = dict(zip(doses, range(len(doses))))

    f, axarr = plt.subplots(n, nc, figsize=(16, 2.5 * n), sharex="col", sharey=True)

    res = []
    for si, ((dose, sample), ddf) in enumerate(sdf.groupby(["Drug.concentration", "Sample"])):    

        multis = multi[multi["sample"] == sample]

        vcf = VCF(f"data/Filtered_VCFS_recovered/SNV/{sample}_filtered_substitutions.vcf.gz", gts012=True)
        # try:
        #     vcf = VCF(f"data/vcf/{sample}.filtered.vcf.gz", gts012=True)
        # except: continue

        sample_res = []
        for v in tqdm.tqdm(vcf):
            ab = v.format("PM")[1][0]
            gq = v.gt_quals[0]
            # if gq < 20: continue
            
            ref, alt = v.REF, v.ALT[0]
                
            sample_res.append(
                {
                    "chrom": v.CHROM,
                    "pos": v.POS,
                    "ref": ref,
                    "alt": alt,
                    # "dp": dp,
                    "ab": ab,
                    "dose": dose,
                }
            )

        sample_res = pd.DataFrame(sample_res)
        # sample_res = sample_res[sample_res["chrom"].isin(["V", "X"])]
        B2C = {"A": "red", "T": "blue", "C": "green", "G": "gray"}

        for ci, (chrom, chrom_df) in enumerate(sample_res.groupby("chrom")):
            # bins = np.arange(1, 2.5e7, 100_000)
            # hist, edges = np.histogram(chrom_df["pos"], bins=bins)
            # axarr[0, ci].scatter(chrom_df["pos"], chrom_df["dp"], c="gray")

            # multi_chrom = multis[multis["chrom"] == chrom]

            # axarr[1, ci].scatter(multi_chrom["pos"], [1.5] * multi_chrom.shape[0], c="k")

            for ref, ref_df in chrom_df.groupby("ref"):
                axarr[si, ci].scatter(ref_df["pos"], ref_df["ab"], c="gray")#B2C[ref], label=ref)
            # sns.despine(ax=axarr[1, ci])
            # sns.despine(ax=axarr[0, ci])
            # axarr[1, ci].set_ylim(0, 1.1)
            axarr[si, ci].set_ylabel("Allele balance")
            # axarr[1, ci].set_xlabel("Chromosomal position (bp)")
            # axarr[1, ci].set_ylabel("Allele balance")
            axarr[si, ci].set_title(f"chr{chrom} @ {dose}\n{sample}")
            sns.despine(ax=axarr[si, ci])
       
    f.tight_layout()
    f.savefig(f"fig/ab/{mutagen}.{genotype}.ab.png", dpi=200)
    plt.close()

    # res = pd.DataFrame(res)
    # f, ax = plt.subplots()

    # sns.histplot(data=res, x="ab", hue="sample", ax=ax, cumulative=True, common_norm=False, stat="proportion", element="step", fill=False)
    # f.savefig(f"fig/ab/{genotype}.{mutagen}.ab.png")
    # plt.close()
    # multi = pd.read_csv("multi.txt", sep="\t", names=["chrom", "pos", "ref", "alt"])


    

# f, ax = plt.subplots()
# sns.kdeplot(data=res, x="ab", ax=ax)
# f.savefig("abdist.png")
