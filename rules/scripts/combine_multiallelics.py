import pandas as pd
import tqdm

from snakemake.script import snakemake

# read in metadata
metadata = pd.read_excel(
    snakemake.input.metadata,
    sheet_name="Sample_description",
    skiprows=1,
)

res = []
mav_counts = []
sb = []

N_PERM = int(snakemake.wildcards.PERM)

for i, fh in tqdm.tqdm(enumerate(snakemake.input.mavs)):
    sample = fh.split("/")[-1].split(".")[0]

    # some files may be empty if there weren't any SNVs to begin with
    try:
        df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    except (pd.errors.EmptyDataError, FileNotFoundError):
        continue


    df["is_multiallelic"] = df["ma_reads"] >= snakemake.params.threshold

    strand_bias = (
        df[df["permuted"] == 0]
        .groupby(["is_multiallelic", "is_strand_biased"])
        .size()
        .reset_index()
        .rename(columns={0: "count"})
    )
    sb.append(strand_bias)

    # store counts
    counts = (
        df.groupby(["sample", "permuted", "is_multiallelic"])
        .size()
        .reset_index()
        .rename(columns={0: "count"})
    )
    mux = pd.MultiIndex.from_product(
        [[sample], list(range(N_PERM + 1)), [True, False]],
        names=["sample", "permuted", "is_multiallelic"],
    )
    counts = (
        counts.set_index(["sample", "permuted", "is_multiallelic"])
        .reindex(mux, fill_value=0)
        .reset_index()
    )
    mav_counts.append(counts)
    
    # append MAVs with at least 1 read of evidence 
    df = df[df["ma_reads"] > 0]
    res.append(df)

# get strand bias
sb = pd.concat(sb)

chi2_vals = (
    sb.groupby(["is_multiallelic", "is_strand_biased"])
    .agg(total=("count", "sum"))
    .reset_index()["total"]
    .values.reshape((2, 2))
)
chi2 = ss.chi2_contingency(chi2_vals)
print (chi2_vals, chi2.pvalue)

res = pd.concat(res).fillna("N:0")
res = res.merge(metadata, left_on="sample", right_on="Sample")

res.to_csv(snakemake.output.mavs_raw, sep="\t", index=False)

OUT_COLS = [
    "chrom",
    "start",
    "end",
    "ref",
    "af",
    "sample",
    "ident",
    "bam_smp",
    "is_strand_biased",
    "summary",
    "ma_reads",
    "Genotype",
    "Mutagen",
    "Drug.concentration",
]


res[(res["ma_reads"] >= snakemake.params.threshold) & (res["permuted"] == 0)][
    OUT_COLS
].sort_values(["chrom", "start", "end"]).to_csv(
    snakemake.output.mavs_filt, sep="\t", index=False
)

mav_counts = pd.concat(mav_counts)
mav_counts.to_csv(snakemake.output.counts, sep="\t", index=False)
