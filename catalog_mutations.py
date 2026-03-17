import pandas as pd
from cyvcf2 import VCF
import matplotlib.pyplot as plt
import tqdm
import numpy as np
import seaborn as sns
import glob
import ruptures as rpt
import os
from global_vars import MUT2EXP, REVCOMP, generate_collapsed_3mers
from mutyper import Ancestor

plt.rc("font", size=12)


ancestor = Ancestor(
    "data/ref/ce11.fa",
    k=3,
    sequence_always_upper=True,
)


metadata = pd.read_excel(
    "41467_2020_15912_MOESM4_ESM.xlsx",
    sheet_name="Sample_description",
    skiprows=1,
)

mutagens = metadata["Mutagen"].unique()
genotypes = metadata["Genotype"].unique()


metadata = metadata[metadata["Mutagen"].isin(mutagens)]
metadata = metadata[metadata["Genotype"].isin(genotypes)]

sample2normal = dict(zip(metadata["Sample"], metadata["Normal"]))


VCF_DIR = "data/Filtered_VCFS_recovered/SNV/{0}_filtered_substitutions.vcf.gz"
# VCF_DIR = "data/vcf/{0}.filtered.vcf.gz"

res = []
# catalog mutations in each sample according to their strand bias
for (sample, genotype, mutagen, dose), sdf in tqdm.tqdm(metadata.groupby(["Sample", "Genotype", "Mutagen", "Drug.concentration"])):

    # if dose == 0: continue

    if not os.path.exists(VCF_DIR.format(sample)): 
        print (sample, genotype, mutagen, dose)
        continue

    vcf = VCF(
        VCF_DIR.format(sample),
        gts012=True,
    )

    # loop over VCF for this sample
    prev_chrom = None
    prev_pos = None
    for v in vcf:
        ref, alt = v.REF, v.ALT[0]

        kmer = ancestor.mutation_type(f"chr{v.CHROM}", v.start, ref, alt)
        if None in kmer: 
            print ('Error!')
            continue
        kmer = ">".join(kmer)

        res.append(
            {
                "sample": sample,
                "genotype": genotype,
                "mutagen": mutagen,
                "dose": float(dose),
                "chrom": v.CHROM,
                "pos": v.POS,
                "kmer": kmer,
            }
        )

res = pd.DataFrame(res)

res.to_csv("mutations.tsv", sep="\t", index=False)

