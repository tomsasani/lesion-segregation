import pandas as pd
import tqdm
import glob
from mutyper.ancestor import Ancestor
from collections import Counter
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
import seaborn as sns
import logomaker
import itertools

from snakemake.script import snakemake

plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans"]


K = 21


def calc_entropy(r):
    # remove central nucleotide
    center = K // 2
    r = r[:center] + r[center + 1 :]
    counts = Counter(r.upper())
    probabilities = np.array(list(counts.values())) / len(r)
    shannon_entropy = ss.entropy(probabilities, base=2)
    return shannon_entropy


ancestor = Ancestor(
    snakemake.input.ref,
    sequence_always_upper=True,
    k=K,
)


multiallelics = []
regular = []

for fh in tqdm.tqdm(snakemake.input.mavs):
    try:
        df = pd.read_csv(fh, sep="\t")
    except pd.errors.EmptyDataError:
        continue
    df["is_multiallelic"] = df["ma_reads"] >= snakemake.params.threshold
    multiallelics.append(df[df["is_multiallelic"]])
    regular.append(df[(~df["is_multiallelic"]) & (df["permuted"] == 0)])

multiallelics = pd.concat(multiallelics).fillna("N:0")
regular = pd.concat(regular).fillna("N:0")

res = pd.concat([multiallelics, regular])

muts = []
for i, row in tqdm.tqdm(res.iterrows()):
    ref = row["ref"]

    for b in row["summary"].split("|"):
        alt, count = b.split(":")
        if alt == ref:
            continue

        mutation = ancestor.mutation_type(
            row["chrom"],
            row["start"],
            ref,
            alt,
        )
        if None in mutation:
            continue
        refseq = list(mutation[0])
        d = row.to_dict()
        d.update(
            {
                "refseq": "".join(refseq),
            }
        )
        muts.append(d)
        break

muts = pd.DataFrame(muts)


muts["entropy"] = muts["refseq"].apply(lambda r: calc_entropy(r))

muts["category"] = muts.apply(lambda row: "MAV (spurious)" if (row["is_multiallelic"] and row["permuted"] > 0) else "MAV (empirical)" if (row["is_multiallelic"] and row["permuted"] == 0) else "SNV (biallelic)", axis=1)

bins = np.arange(muts["entropy"].min(), 2, 0.05)
entropy = np.zeros((3, 100, bins.shape[0] - 1))

categories = ["MAV (empirical)", "SNV (biallelic)", "MAV (spurious)"]
bs = []
for t in range(100):
    _muts = muts.sample(frac=1, replace=True)
    for cat_i, cat in enumerate((categories)):
        vals = _muts[_muts["category"] == cat]["entropy"].values
        hist, edges = np.histogram(vals, bins=bins)
        hist_fracs = hist / np.sum(hist)
        entropy[cat_i, t, :] = np.cumsum(hist_fracs)

cat2ct = muts.groupby("category").size().to_dict()
palette = sns.color_palette("colorblind", len(cat2ct))
f, ax = plt.subplots()
for cat_i, cat in enumerate(categories):
    v = entropy[cat_i]
    lo, hi = np.percentile(v, 2.5, axis=0), np.percentile(v, 97.5, axis=0)
    mean = np.mean(v, axis=0)
    ax.plot(bins[:-1], mean, label=f"{cat} (n = {cat2ct[cat]:,})", color=palette[cat_i])
    ax.fill_between(bins[:-1], lo, hi, alpha=0.1, color=palette[cat_i])

for cat1, cat2 in itertools.combinations(categories, r=2):
    vals1 = muts[muts["category"] == cat1]["entropy"].values
    vals2 = muts[muts["category"] == cat2]["entropy"].values

    print (cat1, cat2, ss.ks_2samp(vals1, vals2).pvalue)

ax.set_xlabel(
    "\nShannon entropy of nucleotide sequence +/- 10bp from MAV\n"
    + "lower "
    + r"$\leftarrow$"
    + " complexity "
    + r"$\rightarrow$"
    + " higher"
)
ax.set_ylabel("Cumulative fraction of MAVs (+/- 95% CI)")
ax.legend(title="Variant type", frameon=True, shadow=True, fancybox=True)
sns.despine(ax=ax)
f.tight_layout()
f.savefig(snakemake.output.entropy, dpi=200)

# f, (ax1, ax2) = plt.subplots(2, figsize=(6, 8), sharex=True, sharey=True)
# seqs = muts[muts["category"] == "MAV (empirical)"]["refseq"].to_list()

# matrix = logomaker.transform_matrix(
#     logomaker.alignment_to_matrix(seqs), from_type="counts", to_type="probability"
# )
# logomaker.Logo(
#     matrix,
#     fade_probabilities=True,
#     # shade_below=0.5,
#     fade_below=0.25,
#     stack_order="small_on_top",
#     ax=ax1,
# )

# seqs = muts[muts["category"] == "MAV (spurious)"]["refseq"].to_list()

# matrix = logomaker.transform_matrix(
#     logomaker.alignment_to_matrix(seqs), from_type="counts", to_type="probability"
# )
# logomaker.Logo(
#     matrix,
#     fade_probabilities=True,
#     # shade_below=0.5,
#     fade_below=0.25,
#     stack_order="small_on_top",
#     ax=ax2,
# )
# ax1.set_title("Sequence context surrounding MAVs\n(at true sites)")
# ax2.set_title("Sequence context surrounding MAVs\n(at random sites)")
# # ax1.axhline(0.64, c="k", ls=":", label="Average GC content genome-wide")
# # ax2.axhline(0.64, c="k", ls=":")
# for ax in (ax1, ax2):
#     sns.despine(ax=ax)
#     ax.set_xticks(np.arange(0, 21, 5))
#     ax.set_xticklabels(["-10", "-5", "0", "+5", "+10"])
#     ax.set_xlabel("Position relative to reference nucleotide")
#     ax.set_ylim(0, 1.1)
# # ax1.legend()
# f.tight_layout()
# f.savefig("fig/mavs.logo.png", dpi=200)

