import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import numpy as np
import tqdm
import scipy.stats as ss
import os
from collections import Counter

from snakemake.script import snakemake

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

# combine all MAVs with >0 evidence

res = pd.concat(res).fillna("N:0")
res = res.merge(metadata, left_on="sample", right_on="Sample")

res.to_csv(snakemake.output.mavs, sep="\t", index=False)

mav_counts = pd.concat(mav_counts)
mav_totals = (
    mav_counts.groupby(["sample", "permuted"]).agg(total=("count", "sum")).reset_index()
)
mav_counts = mav_counts.merge(mav_totals)

mav_counts.to_csv(snakemake.output.counts, sep="\t", index=False)