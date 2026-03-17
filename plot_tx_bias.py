import pandas as pd
import matplotlib.pyplot as plt
import tqdm
import numpy as np
import seaborn as sns
from bx.intervals.intersection import Interval, IntervalTree 
from collections import defaultdict, Counter
import gzip
import csv
import pyfaidx
import os

ref = pyfaidx.Fasta("data/ref/ce11.fa")

plt.rc("font", size=12)

CHROMS = ["I", "II", "III", "IV", "V", "X"]
CHROMS = [f"chr{c}" for c in CHROMS]

cell_type = "Germline"

EXP_THRESH = 0.0001
MIN_MUTATIONS_IN_BLOCK = 20
SYMMETRY_THRESH = 0.8

MUTAGEN = "EMS"

tx_data = pd.read_csv("data/Mitotic germ cells.csv")

# is_expressed = tx_data[tx_data["scaled_TPM"] > EXP_THRESH]
# not_expressed = tx_data[tx_data["scaled_TPM"] <= EXP_THRESH]
# not_expressed["germline_exp_bin"] = 0

n_bins = 5

exp_bins, exp_bin_labels = pd.qcut(
    tx_data["scaled_TPM"].astype(float),
    n_bins,
    duplicates="drop",
    labels=list(range(n_bins)),
    retbins=True,
)
tx_data["germline_exp_bin"] = exp_bins
# tx_data = pd.concat([is_expressed, not_expressed])

genes_to_query = tx_data["gene_ID"].unique()

gene2exp = dict(zip(tx_data["gene_ID"], tx_data["germline_exp_bin"]))

revcomp = {"A": "T", "T": "A", "C": "G", "G": "C"}
mut2exp = {
    "EMS": "G",
    "MMS": "A",
}

metadata = pd.read_excel(
    "41467_2020_15912_MOESM4_ESM.xlsx",
    sheet_name="Sample_description",
    skiprows=1,
)

mutagens = [MUTAGEN]
genotypes = ["N2"]

metadata = metadata[metadata["Mutagen"].isin(mutagens)]
metadata = metadata[metadata["Genotype"].isin(genotypes)]

mutation_tree = defaultdict(IntervalTree)
segment_tree = defaultdict(IntervalTree)

segments = []

### need to calc rate using number of appropriate mut2exp nuceotides in region

mutations = pd.read_csv("strand_bias.tsv", sep="\t")

for (sample, genotype, mutagen, dose), sdf in mutations.groupby(["sample", "genotype", "mutagen", "dose"]):

    if genotype not in genotypes or mutagen not in mutagens: continue


    # get breakpoints
    bp_fh = f"csv/{mutagen}.{genotype}.{dose}.{sample}.breakpoints.tsv"

    if not os.path.exists(bp_fh): continue
    breakpoints = pd.read_csv(bp_fh, sep="\t")

    if breakpoints.shape[0] == 0: continue

    # filter to confident breakpoint blocks
    breakpoints = breakpoints[~breakpoints["symmetry"].between(-SYMMETRY_THRESH, SYMMETRY_THRESH)]
    breakpoints = breakpoints[breakpoints["n_inf_mutations"] >= MIN_MUTATIONS_IN_BLOCK]

    if breakpoints.shape[0] == 0: continue

    # add breakpoints to a GLOBAL interval tree
    for i,row in breakpoints.iterrows():
        segment_tree[row["chrom"]].insert_interval(Interval(row["start"], row["end"], strand=row["strand"], value=sample))
    
     # annotate each mutation with a status corresponding to whether it's in
    # a confidently phased interval
    sdf["interval_annotation"] = sdf.apply(
        lambda row: ",".join(
            [
                o.strand if o.value == sample else "?"
                for o in segment_tree[row["chrom"]].find(row["pos"] - 1, row["pos"])
            ]
        ),
        axis=1,
    )
    sdf["interval_annotation"] = sdf["interval_annotation"].apply(lambda a: 0 if a == "watson" else 1 if a == "crick" else -1)
    # filter to only include mutations that occur in a strand-resolved interval
    # and occur on that strand
    sdf = sdf[sdf["interval_annotation"] == sdf["strand"]]
    sdf = sdf[sdf["interval_annotation"] != -1]

    # add strand-resolved mutations to tree
    for i, row in sdf.iterrows():
        chrom, start, end = row["chrom"], row["pos"] - 1, row["pos"]
        strand = row["strand"]
        mutation_tree[f"chr{chrom}"].insert_interval(
            Interval(start, end, strand=strand, value=sample),
        )


# store lesion resolved mutations in each transcript
mutations_in_transcripts = []
# store the total length of each transcript
transcript_lengths = []

seen_genes = []
with gzip.open("data/caenorhabditis_elegans.PRJNA13758.WBPS19.canonical_geneset.gtf.gz", "rt") as infh:
    csvf = csv.reader(infh, delimiter="\t")
    header = ["chrom", "origin", "kind", "start", "end", "?", "strand", "??", "misc"]

    for i,l in tqdm.tqdm(enumerate(csvf)):
        if i == 0: continue

        d = dict(zip(header, l))
        if d["kind"] != "transcript": continue

        chrom, start, end = f"chr{d['chrom']}", int(d["start"]), int(d["end"])
        if chrom not in CHROMS: continue

        # make sure this gene is within a complete segment
        overlapping_segments = segment_tree[chrom.lstrip("chr")].find(start, end)
        if len(overlapping_segments) == 0: continue

        # NOTE: need to make sure that sample matches?
        # NOTE: we're checking here that the transcript has to be in *any* hi-quality
        # segment, but currently we're not confirming that the mutations we're adding 

        overlaps = False
        for o in overlapping_segments:
            if o.start <= start - 100_000 and o.end >= end + 100_000: overlaps = True
        if not overlaps: continue


        symbol = d["misc"].split(";")[0].split(" ")[1].strip("\"")
        if symbol in seen_genes: continue
        seen_genes.append(symbol)

        # get nucleotide context of this gene
        nucs = ref[chrom][start:end].seq.upper()
        nuc_count = dict(Counter(nucs).most_common())

        for n in ("A", "C", "T", "G"):
            if n not in nuc_count:
                nuc_count[n] = 0

        # check that this is a gene that's expressed in the cell
        # type of interest
        if symbol not in gene2exp: 
            expression_bin = 0
        else:
            expression_bin = gene2exp[symbol]

        # record length and expression of transcript, regardless of whether
        # it contains any mutations
        transcript_lengths.append(
            {
                "symbol": symbol,
                "chrom": chrom,
                "tx_s": start,
                "tx_e": end,
                "length": end - start,
                "expression_bin": expression_bin,
                "As": nuc_count["A"],
                "Ts": nuc_count["T"],
                "Cs": nuc_count["C"],
                "Gs": nuc_count["G"],
            }
        )

        # get overlapping mutations
        overlapping_mutations = mutation_tree[chrom].find(start, end)

        # if the strand of this transcript is (+), then, it's transcribed 5' to 3'
        # using the (-) strand as the template
        template_strand = "-" if d["strand"] == "+" else "+"

        template_strand_lesions = 0
        non_template_strand_lesions = 0

        # if no mutations overlap, move on
        if len(overlapping_mutations) == 0: 
            mutations_in_transcripts.append(
                {
                    "symbol": symbol,
                    "gene_strand": d["strand"],
                    "template_strand_lesions": template_strand_lesions,
                    "non_template_strand_lesions": non_template_strand_lesions,
                    "expression_bin": expression_bin,
                })
            continue

        # otherwise, record the strand of the mutation and the
        # direction of transcription. if the direction of transcription
        # is +, the - strand is used as a template, and vice versa.

        # NOTE: here, do we need to confirm that the mutation is 
        
        for o in overlapping_mutations:
            # ensure that this mutation is present in an overlapping segment from the same
            # sample
            if o.value not in [p.value for p in overlapping_segments]:
                continue
            # check if the lesion is on the template strand
            lesion_strand = o.strand
            # it's a template strand lesion if the mutation is strand 0 and 
            # the transcript is -.
            if (lesion_strand == 0 and template_strand == "-") or (lesion_strand == 1 and template_strand == "+"):
                template_strand_lesions += 1
            else:
                non_template_strand_lesions += 1

        mutations_in_transcripts.append(
            {
                "symbol": symbol,
                "gene_strand": d["strand"],
                "template_strand_lesions": template_strand_lesions,
                "non_template_strand_lesions": non_template_strand_lesions,
                "expression_bin": expression_bin,
            }
        )

# now, we have an entry for every transcript with its number of template
# strand and non template strand mutations.
# we also have an entry for every transcript with its number of A/T/C/G nucs
mutations = pd.DataFrame(mutations_in_transcripts)
lengths = pd.DataFrame(transcript_lengths)
# merge the two dataframes
merged = mutations.merge(lengths)

# EMS lesions happen on Gs. so, our template strand mutation rate is # of template strand
# lesions divided by # of Gs (in genes where with strand -), # of template strand mutations
# # divided by # of Cs (in genes where strand is +).  non-template strand mutation rate is # of non-template strand
# lesions divided by # of Gs (in genes with strand +) and # of non-template lesions divided by # of Cs (in genes with strand -).

exp_watson = mut2exp[MUTAGEN]
exp_crick = revcomp[mut2exp[MUTAGEN]]

print (exp_watson, exp_crick)
print (mutations.query("template_strand_lesions > 0"))

merged["template_denominator"] = merged.apply(
    lambda row: (
        row[f"{exp_crick}s"]
        if row["gene_strand"] == "+"
        else row[f"{exp_watson}s"] if row["gene_strand"] == "-" else -1
    ),
    axis=1,
)
merged["non_template_denominator"] = merged.apply(
    lambda row: (
        row[f"{exp_watson}s"]
        if row["gene_strand"] == "+"
        else row[f"{exp_crick}s"] if row["gene_strand"] == "-" else -1
    ),
    axis=1,
)

# bootstrap
bs = []
for _ in range(100):
    _bs = merged.sample(frac=1, replace=True)

    _merged = _bs.groupby(["expression_bin"]).agg(
        n=("symbol", lambda s: len(s)),
        template_denom=("template_denominator", "sum"),
        template_num=("template_strand_lesions", "sum"),
        non_template_denom=("non_template_denominator", "sum"),
        non_template_num=("non_template_strand_lesions", "sum"),
    ).reset_index()

    _merged["template_rate"] = _merged["template_num"] / _merged["template_denom"]
    _merged["non_template_rate"] = _merged["non_template_num"] / _merged["non_template_denom"]
    bs.append(_merged)
# print(merged.pivot(index="expression_bin", columns=["template_rate", "non_template_rate"]))
bs = pd.concat(bs)
print (bs)
f, ax = plt.subplots(figsize=(8, 6))
# sns.pointplot(data=bs, x="expression_bin", y="template_rate")
# sns.pointplot(data=bs, x="expression_bin", y="non_template_rate")
ax.plot(bs["expression_bin"], bs["template_rate"], label="template", alpha=0.2)
ax.plot(bs["expression_bin"], bs["non_template_rate"], label="non_template", alpha=0.2)
ax.legend()
sns.despine(ax=ax)
f.tight_layout()
f.savefig("o.png", dpi=200)
