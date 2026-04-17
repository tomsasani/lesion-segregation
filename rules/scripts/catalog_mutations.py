import pandas as pd
from cyvcf2 import VCF
import matplotlib.pyplot as plt
import tqdm
import os
from mutyper import Ancestor
from bx.intervals.intersection import Interval, IntervalTree
from collections import defaultdict
import csv
import gzip


from snakemake.script import snakemake

plt.rc("font", size=12)


# ancestor = Ancestor(
#     snakemake.input.ref,
#     k=3,
#     sequence_always_upper=True,
# )


metadata = pd.read_excel(
    snakemake.input.metadata,
    sheet_name="Sample_description",
    skiprows=1,
)
sample_metadata = metadata[metadata["Sample"] == snakemake.wildcards.SAMPLE]
assert sample_metadata.shape[0] == 1

res = []
# catalog mutations in each sample according to their strand bias

vcf = VCF(snakemake.input.vcf, gts012=True)

exclude_tree = defaultdict(IntervalTree)
with gzip.open(snakemake.input.exclude, "rt") as infh:
    csvf = csv.reader(infh, delimiter="\t")
    for i, l in enumerate(csvf):
        if i == 0:
            header = ["chrom", "start", "end", "annot", "a", "strand"]
            
        d = dict(zip(header, l))
        exclude_tree[d["chrom"]].insert_interval(Interval(int(d["start"]), int(d["end"])))


# loop over VCF for this sample

for v in vcf:

    if v.FILTER not in ("PASS", None): continue

    if len(exclude_tree[v.CHROM].find(v.start, v.end)) > 0:
        continue

    ref = v.REF
    # NOTE: by counting each site twice at multi-alelics, am i double counting?
    passes = False
    afs = []
    for alt_i, alt in enumerate(v.ALT):

        if len(alt) > 1: continue
        af = v.format("PM")[1][0]
        afs.append(af)
        # print (af)
        # if af < 0.2: continue

        # dp = v.format("DP")[0][0]
        # if not 15 <= dp <= 150: continue

        # ad = v.format("AD")[0][alt_i + 1]
        # if ad < 5: continue

        passes = True

    if not passes: continue
    # just pick one alt for the purposes of logging
    # kmer = ancestor.mutation_type(v.CHROM, v.start, ref, alt)
    # if None in kmer:
    #     print('Error!')
    #     continue
    # kmer = ">".join(kmer)

    res.append(
        {
            "chrom": v.CHROM,
            "start": v.start,
            "end": v.end,
            "ref": ref,
            "alt": ":".join(v.ALT),
            "af": max(afs),
            "sample": snakemake.wildcards.SAMPLE,
            "genotype": sample_metadata["Genotype"].unique()[0],
            "mutagen": sample_metadata["Mutagen"].unique()[0],
            "dose": float(sample_metadata["Drug.concentration"].unique()[0]),
        }
    )

res = pd.DataFrame(res)

res.to_csv(snakemake.output.mutations, sep="\t", index=False)
