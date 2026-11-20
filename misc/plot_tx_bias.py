import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import ruptures as rpt
import tqdm
from bx.intervals.intersection import Interval, IntervalTree
from collections import defaultdict, OrderedDict, Counter
import csv
import gzip
import mutyper
import glob

ancestor = mutyper.Ancestor(
    "data/ref/ce11.fa.gz",
    k=1,
    sequence_always_upper=True,
)


def convert(row: pd.Series):
    revcomp = {"A": "T", "T": "A", "C": "G", "G": "C"}

    ref, alt = row["ref"], row["alt"]
    if len(ref) > 1 or len(alt) > 1:
        mut = "indel"
    else:
        if ref not in ("C", "A"):
            ref, alt = revcomp[ref], revcomp[alt]
        mut = f"{ref}>{alt}"
    return mut

# define segments!
SYMMETRY_THRESH = 0.7

MUT2EXP = {
    "EMS": "G",
    "UV": "C",
    "ENU": "T",
    "MMS": "A",
}

REVCOMP = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
}


MUTAGEN = "MMS"

chrom2idx = dict(zip(["I", "II", "III", "IV", "V", "X"], range(6)))

mutations = pd.read_csv("/scratch/ucgd/lustre-labs/quinlan/u1006375/lesion-segregation/csv/mutations/combined.tsv", sep="\t")

mutations = mutations[mutations["ref"].str.len() == 1]
mutations = mutations[mutations["alt"].str.len() == 1]


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


# figure out the expected "watson" and "crick" strands.
exp_watson = MUT2EXP[MUTAGEN]
exp_crick = REVCOMP[exp_watson]

mutations = mutations[mutations["Mutagen"] == MUTAGEN]
mutations = mutations[mutations["Genotype"] == "polk-1"]

# if lesion affects fwd nucleotide G, it'll look like a G>N.
# if lesion affects fwd nucleotide C (rev nucleotide G), it'll look like a C>N.
mutations["strand"] = np.select(
    [
        mutations["ref"].eq(exp_watson),
        mutations["ref"].eq(exp_crick),
    ],
    [0, 1],
    default=-1,
).astype(np.int8)

counts = mutations.groupby("sample").size().reset_index().rename(columns={0: "count"})

top_n = counts.sort_values("count", ascending=False).head(50)["sample"].to_list()

mutations = mutations[mutations["sample"].isin(top_n)]

breakpoints = []
for fh in glob.glob("csv/breakpoints/*.tsv"):
    df = pd.read_csv(fh, sep="\t")
    breakpoints.append(df)
breakpoints = pd.concat(breakpoints)

breakpoints = breakpoints[(breakpoints["symmetry"] >= SYMMETRY_THRESH) | (breakpoints["symmetry"] <= -SYMMETRY_THRESH)]

bins = np.arange(-1, 1.05, 0.05)
f, ax = plt.subplots()
sns.histplot(data=breakpoints, x="symmetry", ax=ax, bins=bins)
f.savefig("symm.png")


txs = defaultdict(IntervalTree)
seen = []
total_tx = 0
template_nuc_dict = Counter({"A": 0, "T": 0, "C": 0, "G": 0})
non_template_nuc_dict = Counter({"A": 0, "T": 0, "C": 0, "G": 0})

with gzip.open("data/caenorhabditis_elegans.PRJNA13758.WBPS19.canonical_geneset.gtf.gz", "rt") as infh:
    csvf = csv.reader(infh, delimiter="\t")
    for i, l in tqdm.tqdm(enumerate(csvf)):
        if i == 0: continue
        chrom, start, end = l[0], l[3], l[4]
        if chrom not in chrom2idx:
            continue
        kind = l[2]
        if kind != "transcript":
            continue
        biotype = l[8].split(";")[4].split(" ")[2].strip('"')
        if biotype != "protein_coding":
            continue
        if len(txs[chrom].find(int(start), int(end))) > 0:
            continue
        strand = l[6]
        nucs = ancestor[f"chr{chrom}"][int(start):int(end)]
        nuc_cts = Counter(str(nucs))
        # if tx is on fwd, template is on rev
        if strand == "+":
            template_nuc_cts = {REVCOMP[k]: v for k, v in nuc_cts.items()}
            non_template_nuc_cts = {k: v for k, v in nuc_cts.items()}
        elif strand == "-":
            non_template_nuc_cts = {REVCOMP[k]: v for k, v in nuc_cts.items()}
            template_nuc_cts = {k: v for k, v in nuc_cts.items()}
        template_nuc_dict.update(template_nuc_cts)
        non_template_nuc_dict.update(non_template_nuc_cts)
        tx_id = l[8].split(";")[2].split(" ")[2].strip('"')
        # if the transcript is (-) stranded, the forward strand is used
        # as the template.
        txs[chrom].insert_interval(
            Interval(
                int(start),
                int(end),
                strand=0 if strand == "-" else 1 if strand == "+" else -1,
                value={"tx": tx_id}
            )
        )
        # seen.append(f"{chrom}:{start}-{end}")
        total_tx += (int(end) - int(start))

print (template_nuc_dict)
print (non_template_nuc_dict)
def annotate(row, tx_tree, breakpoints):
    strain = row["sample"]
    chrom, start = row["chrom"], row["start"]

    bps = breakpoints[(breakpoints["sample"] == strain) & (breakpoints["chrom"] == chrom)]

    # filter to overlapping bps
    bps = bps[(bps["start"] <= row["start"]) & (bps["end"] >= row["start"])]
    if bps.shape[0] == 0:
        return "UNK"
    else:
        # SNV defines the breakpoint
        if bps.shape[0] > 1:
            return "UNK"
        
        lesion_strand = bps["strand"].values[0]
        if lesion_strand == -1:
            return "UNK"

        tx_overlaps = tx_tree[row["chrom"].strip("chr")].find(row["start"], row["start"] + 1)
        tx_strands = [o.strand for o in tx_overlaps]

        if len(set(tx_strands)) != 1:
            # two genes in opposing directions or no genes at all
            return "UNK"
        else:
            tx_strand = tx_strands[0]
            # if lesion on fwd strand, snv_strand == 0.
            # if fwd strand is the template, tx_strand == "-".
            lesion_on_tx_strand = tx_strand == lesion_strand
            return str(int(lesion_on_tx_strand))

# annotate each mutation with its transcriptional strand

mutations["tx_strand"] = mutations.apply(lambda row: annotate(row, txs, breakpoints), axis=1)

mutations = mutations[mutations["tx_strand"] != "UNK"]
# mutations = mutations.query('chrom == "I"')

cts = mutations.groupby(["ref", "alt"]).size().reset_index().rename(columns={0: "count"})
cts["frac"] = cts["count"] / mutations.shape[0]
cts["k"] = cts.apply(lambda row: row["ref"] + row["alt"], axis=1)
spec = mutations.groupby(["ref", "alt", "tx_strand"]).size().reset_index().rename(columns={0: "count"})
tot = spec.groupby("tx_strand").agg(tot=("count", "sum")).reset_index()
spec = spec.merge(tot)
print (spec)
spec["frac"] = spec["count"] / spec["tot"]
spec["k"] = spec.apply(lambda row: row["ref"] + row["alt"], axis=1)
f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(8, 4))
sns.barplot(data=spec, x="k", y="frac", hue="tx_strand", ax=ax1)
sns.barplot(data=cts, x="k", y="frac", ax=ax2)
f.savefig('o.png')
