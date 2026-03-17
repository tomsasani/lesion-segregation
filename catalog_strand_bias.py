import pandas as pd
from cyvcf2 import VCF
import matplotlib.pyplot as plt
import tqdm
import numpy as np
import seaborn as sns
import glob
import ruptures as rpt
import os
from global_vars import MUT2EXP, REVCOMP

plt.rc("font", size=12)

revcomp = {"A": "T", "T": "A", "C": "G", "G": "C"}

metadata = pd.read_excel(
    "41467_2020_15912_MOESM4_ESM.xlsx", sheet_name="Sample_description", skiprows=1
)

mutagens = ["EMS", "MMS", "DMS"]

genotypes = [
    "N2",
    "polk-1",
]

metadata = metadata[metadata["Mutagen"].isin(mutagens)]
metadata = metadata[metadata["Genotype"].isin(genotypes)]

sample2normal = dict(zip(metadata["Sample"], metadata["Normal"]))

# map combinations of genotypes and mutagens to concentrations and get max per treatment
treatment2conc = (
    metadata.groupby(["Genotype", "Mutagen"])
    .agg(conc_max=("Drug.concentration", "max"))
    .to_dict()["conc_max"]
)


n_samples = metadata["Sample"].nunique()

VCF_DIR = "data/Filtered_VCFS_recovered/SNV/{0}_filtered_substitutions.vcf.gz"
# VCF_DIR = "data/vcf/{0}.filtered.vcf.gz"

res = []
# catalog mutations in each sample according to their strand bias
for (sample, genotype, mutagen, dose), sdf in metadata.groupby(["Sample", "Genotype", "Mutagen", "Drug.concentration"]):

    if not os.path.exists(VCF_DIR.format(sample)): continue

    vcf = VCF(
        VCF_DIR.format(sample),
        gts012=True,
    )

    # figure out the expected "watson" and "crick" strands.
    exp_watson = MUT2EXP[mutagen]
    exp_crick = REVCOMP[exp_watson]

    

    # loop over VCF for this sample
    prev_chrom = None
    prev_pos = None
    for v in tqdm.tqdm(vcf):
        ref, alt = v.REF, v.ALT[0]


        dist = 0
        if prev_chrom is None:
            prev_chrom = v.CHROM
            prev_pos = v.POS
        else:
            if v.CHROM == prev_chrom:
                dist = v.POS - prev_pos
                prev_chrom = v.CHROM
                prev_pos = v.POS
            else:
                dist = 0
                prev_chrom = v.CHROM
                prev_pos = v.POS


        # determine the strand of this mutation.
        # 0 if the REF matches the watson starnd, 1 if it matches crick
        strand = 0 if ref == exp_watson else 1 if ref == exp_crick else -1
        
        res.append(
            {
                "sample": sample,
                "genotype": genotype,
                "mutagen": mutagen,
                "dose": float(dose),
                "chrom": v.CHROM,
                "pos": v.start,
                "ref": ref,
                "alt": alt,
                "strand": strand,
                "dist": dist,
            }
        )

res = pd.DataFrame(res)

res.to_csv("strand_bias.tsv", sep="\t", index=False)
