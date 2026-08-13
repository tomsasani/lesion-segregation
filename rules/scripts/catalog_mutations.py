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


ancestor = Ancestor(
    snakemake.input.ref,
    k=3,
    sequence_always_upper=False,
)


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


source = snakemake.params.var_source

# loop over VCF for this sample
for v in vcf:

    if v.var_type != "snp":
        continue

    if v.FILTER not in ("PASS", None): 
        continue

    chrom = v.CHROM
    if source == "caveman":
        chrom = f"chr{v.CHROM}"

    if len(exclude_tree[chrom].find(v.start, v.end)) > 0:
        continue

    ref = v.REF
    # NOTE: by counting each site twice at multi-alelics, am i double counting?
    passes = False
    if source == "deepsomatic":
        af = v.format("VAF")[0][0]
    elif source == "caveman":
        af = v.format("PM")[1][0]
    else:
        raise AttributeError
        
    alt = v.ALT[0][0]

    kmer = ancestor.mutation_type(chrom, v.start, ref, alt)
    if None in kmer:
        continue
    kmer = ">".join(kmer)

    if source == "caveman":
        dp = 0
        for direction in ("F", "R"):
            for nuc in ("A", "T", "C", "G"):
                _dp = v.format(f"{direction}{nuc}Z")[-1]
                dp += _dp
    else:
        dp = v.format("DP")[0]

    if source == "deepsomatic":
        gq = v.gt_quals[0]
        if gq < 20: 
            continue

    res.append(
        {
            "chrom": chrom,
            "start": v.start,
            "end": v.end,
            "ref": ref,
            "alt": ":".join(v.ALT),
            "af": af,
            "sample": snakemake.wildcards.SAMPLE,
            "genotype": sample_metadata["Genotype"].unique()[0],
            "mutagen": sample_metadata["Mutagen"].unique()[0],
            "dose": float(sample_metadata["Drug.concentration"].unique()[0]),
        }
    )

res = pd.DataFrame(res)

res.to_csv(snakemake.output.mutations, sep="\t", index=False)
