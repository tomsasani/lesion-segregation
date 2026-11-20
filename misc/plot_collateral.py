## are signatures of collateral mutagenesis similar by damage type?
## is SBS5 a common signature of collateral mutagenesis

import pandas as pd
from cyvcf2 import VCF
import matplotlib.pyplot as plt
import tqdm
import numpy as np
import seaborn as sns
import ruptures as rpt
import os

plt.rc("font", size=12)


# NOTE! need to use all mutations here and look for collateral near the strand resolved ones
mutations = pd.read_csv("strand_resolved.tsv", sep="\t")

# limit to sites where the strand of every focal mutations matches the inferred segment bias,
# OR the mutation is not on watson or crick (colalteral)
# mutations = mutations[(mutations["strand"] == -1) | (mutations["strand"] == mutations["segment_label"])]

res = []

for (sample, genotype, mutagen, dose), sdf in mutations.groupby(["sample", "genotype", "mutagen", "dose"]):

    for ci, (chrom, chrom_df) in enumerate(sdf.groupby("chrom")):
        chrom_df = chrom_df.sort_values("pos").reset_index()

        # ignore first DNM
        cdnms = np.where(chrom_df["dist"].values < 100)[0][1:]
        for i in cdnms:
            a = chrom_df.iloc[i - 1]
            b = chrom_df.iloc[i]

            strand = a["segment_label"]
            # assert strand == b["segment_label"]

            # if we're on the watson strand, first mutation is focal.
            focal = a.to_dict()
            collateral = b.to_dict()
            if strand == b["segment_label"]:
                if strand == 0:
                    focal = a.to_dict()
                    collateral = b.to_dict()
                else:
                    focal = b.to_dict()
                    collateral = a.to_dict()


            for l, k in zip(("focal", "collateral"), (focal, collateral)):
                c = {
                    "chrom": chrom,
                    "interval": f"{a['pos']}_{b['pos']}",
                    "which": l,
                    "pos": k["pos"],
                    "ref": k["ref"],
                    "alt": k["alt"],
                    "kmer": f"{k['ref']}>{k['alt']}",
                    "strand": strand,
                    "dist": b["pos"] - a["pos"],
                    "sample": sample,
                    "mutagen": mutagen,
                    "genotype": genotype,
                    "dose": dose,
                }
                res.append(c)


cdnms = pd.DataFrame(res)
cdnms.to_csv("cdnms.tsv", sep="\t", index=False)


cdnms = cdnms.groupby(["genotype", "mutagen", "dose", "which", "kmer"]).size().reset_index().rename(columns={0: "count"})

cdnms = cdnms[cdnms["which"] == "collateral"]
g = sns.FacetGrid(data=cdnms, row="genotype", col="dose", aspect=3)
g.map(sns.barplot, "kmer", "count")
g.savefig("o.png")
