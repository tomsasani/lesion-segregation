import pandas as pd
from cyvcf2 import VCF
from bx.intervals.intersection import Interval, IntervalTree
from collections import defaultdict
import csv
import gzip
from mutyper import Ancestor

from snakemake.script import snakemake

# read in VCF with sample-private mutations
vcf = VCF(snakemake.input.vcf, gts012=True)

# read in reference genome
ancestor = Ancestor(
    snakemake.input.ref,
    k=3,
    sequence_always_upper=False,
)

# create a binary search tree from the exclude regions
exclude_tree = defaultdict(IntervalTree)
with gzip.open(snakemake.input.exclude, "rt") as infh:
    csvf = csv.reader(infh, delimiter="\t")
    for i, l in enumerate(csvf):
        if i == 0:
            header = ["chrom", "start", "end", "annot"]
        d = dict(zip(header, l))
        exclude_tree[d["chrom"]].insert_interval(
            Interval(int(d["start"]), int(d["end"]))
        )

res = []

source = snakemake.params.var_source

# loop over VCF for this sample
for v in vcf:

    if v.var_type != "snp":
        continue
    if v.FILTER not in ("PASS", None):
        continue

    # chromosome may need chr prefix
    chrom = v.CHROM
    if source == "caveman":
        chrom = f"chr{v.CHROM}"

    # check for exclude region overlap
    if len(exclude_tree[chrom].find(v.start, v.end)) > 0:
        continue

    ref = v.REF

    # use "first" ALT allele for mutyper
    kmer = ancestor.mutation_type(chrom, v.start, ref, v.ALT[0])
    # one reason for None is a low-complexity sequence
    if None in kmer:
        continue
    kmer = ">".join(kmer)

    # allele frequency field differs based on
    # variant caller output
    if source == "deepsomatic":
        af = v.format("VAF")[0][0]
    elif source == "caveman":
        af = v.format("PM")[1][0]
    else:
        raise AttributeError

    # total depth also varies based on caller
    if source == "caveman":
        dp = 0
        for direction in ("F", "R"):
            for nuc in ("A", "T", "C", "G"):
                _dp = v.format(f"{direction}{nuc}Z")[-1]
                dp += _dp
    else:
        dp = v.format("DP")[0]

    # filter on deepsomatic calls
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
        }
    )

res = pd.DataFrame(res)

res.to_csv(snakemake.output.mutations, sep="\t", index=False)
