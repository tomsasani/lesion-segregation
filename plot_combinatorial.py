import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans"]

df = pd.read_csv("combinatorial.tsv", sep="\t")

res = []

unlinked = 0

f, ax = plt.subplots()

for (sample, chrom, pos), sample_df in df.groupby(["sample", "chrom", "pos_1"]):
    # convert to binary summaries
    sample_df["binary"] = sample_df.apply(lambda row: "-".join([str(1 - int(row["nuc_1"] == row["ref_1"])), str(1 - int(row["nuc_2"] == row["ref_2"]))], ), axis=1)
    if all([s in sample_df["binary"].to_list() for s in ("0-0", "1-0", "0-1")]):
        unlinked += 1
        continue
    else:
        print (sample_df["binary"])
    total = sample_df["count"].sum()
    ref_1 = sample_df[sample_df["nuc_1"] == sample_df["ref_1"]]["count"].sum()
    ref_2 = sample_df[sample_df["nuc_2"] == sample_df["ref_2"]]["count"].sum()

    alt_1 = 1 - (ref_1 / total)
    alt_2 = 1 - (ref_2 / total)
    if alt_1 < alt_2:
        alt_1, alt_2 = alt_2, alt_1

    assert alt_1 >= alt_2


    ax.scatter([0, 1], [alt_1, alt_2])
    ax.plot([0, 1],[alt_1, alt_2] )
    # res.append({"chrom": chrom, "pos": pos, "sample": sample, "alt_1": alt_1,})
    # break

f.savefig("paired.png")


# spectrum?
res = []
for (sample, chrom, pos), sample_df in df.groupby(["sample", "chrom", "pos_1"]):
    # convert to binary summaries
    sample_df["binary"] = sample_df.apply(lambda row: "-".join([str(1 - int(row["nuc_1"] == row["ref_1"])), str(1 - int(row["nuc_2"] == row["ref_2"]))], ), axis=1)
    if all([s in sample_df["binary"].to_list() for s in ("0-0", "1-0", "0-1")]):
        unlinked += 1
        continue
    else:
        print (sample_df["binary"])

    total = sample_df["count"].sum()
    ref_1 = sample_df[sample_df["nuc_1"] == sample_df["ref_1"]]["count"].sum()
    ref_2 = sample_df[sample_df["nuc_2"] == sample_df["ref_2"]]["count"].sum()

    alt_1 = 1 - (ref_1 / total)
    alt_2 = 1 - (ref_2 / total)

    sample_df["mut_1"] = sample_df.apply(lambda row: row["ref_1"] + ">" + row["alt_1"], axis=1)
    sample_df["mut_2"] = sample_df.apply(lambda row: row["ref_2"] + ">" + row["alt_2"], axis=1)

    # treat more frequent mutation as the "first"
    if alt_1 >= alt_2:
        a_mut = sample_df["mut_1"].to_list()[0]
        b_mut = sample_df["mut_2"].to_list()[0]
    else:
        a_mut = sample_df["mut_2"].to_list()[0]
        b_mut = sample_df["mut_1"].to_list()[0]
    print (sample_df)
    for mut, l in zip((a_mut, b_mut), ("a", "b")):
        res.append({"mut": mut, "which": l, "mutagen": sample_df["mutagen"].to_list()[0]})

res = pd.DataFrame(res)
print (res)

d = {"G>A": "C>T", "G>C": "C>G", "G>T": "C>A", "T>A": "A>T", "T>C": "A>G", "T>G": "A>C",}

res["mut"] = res["mut"].apply(lambda m: d[m] if m in d else m)

cts = res.groupby(["mut", "which", "mutagen"]).size().reset_index().rename(columns={0: "count"})

g = sns.FacetGrid(data=cts, row="mutagen")
g.map(sns.barplot, "mut", "count", "which")
g.add_legend()
g.savefig("spec.png")