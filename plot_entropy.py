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

plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans"]


def calc_entropy(r):
    counts = Counter(r.upper())
    probabilities = np.array(list(counts.values())) / len(r)
    shannon_entropy = ss.entropy(probabilities, base=2)
    return shannon_entropy


K = 21
ancestor = Ancestor(
    "data/ref/ce10.renamed.fa",
    sequence_always_upper=True,
    k=K,
)

THRESHOLD = 3
res = []
for fh in tqdm.tqdm(glob.glob("csv/multiallelics/*.100.tsv")):
    try:
        df = pd.read_csv(fh, sep="\t")
    except pd.errors.EmptyDataError:
        continue
    df["is_multiallelic"] = df["ma_reads"] >= THRESHOLD
    df = df[(df["mutagen"] == "EMS") & (df["permuted"] == 0)]
    res.append(df)
res = pd.concat(res).fillna("N:0")


res["depth"] = res["summary"].apply(lambda s: sum([int(_s.split(":")[1]) for _s in s.split("|")]))

res = res[res["depth"].between(15, 150)]


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
        if None in mutation: continue
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

muts["empirical"] = muts["permuted"] == 0


bins = np.arange(0.5, 2, 0.05)
entropy = np.zeros((2, 100, bins.shape[0] - 1))

# bs = []
# for t in range(100):
#     _muts = muts.sample(frac=1, replace=True)
#     for i in range(2):
#         if i == 0:
#             hist, edges = np.histogram(_muts[_muts["empirical"]]["entropy"].values, bins=bins)
#         else:
#             hist, edges = np.histogram(_muts[~_muts["empirical"]]["entropy"].values, bins=bins)

#         hist_fracs = hist / np.sum(hist)
#         entropy[i, t, :] = np.cumsum(hist_fracs)

# f, ax = plt.subplots()
# for i in range(2):
#     v = entropy[i]
#     lo, hi = np.percentile(v, 2.5, axis=0), np.percentile(v, 97.5, axis=0)
#     mean = np.mean(v, axis=0)
#     ax.plot(bins[:-1], mean, label="empirical" if i == 0 else "random")
#     ax.fill_between(bins[:-1], lo, hi, alpha=0.1)

# ax.set_xlabel("\nShannon entropy of nucleotide sequence +/- 10bp from MAV\n" + "lower " + r"$\leftarrow$" + " complexity " + r"$\rightarrow$" + " higher")
# ax.set_ylabel("Cumulative fraction of MAVs (+/- 95% CI)")
# ax.legend(title="MAV source", frameon=True, shadow=True, fancybox=True)
# sns.despine(ax=ax)
# f.tight_layout()
# f.savefig('entropy.png', dpi=200)

f, (ax1, ax2) = plt.subplots(2, figsize=(6, 8), sharex=True, sharey=True)
seqs = muts[(muts["is_multiallelic"]) & (muts["empirical"])]["refseq"].to_list()

matrix = logomaker.transform_matrix(logomaker.alignment_to_matrix(seqs), from_type="counts", to_type="probability")
logomaker.Logo(
    matrix,
    fade_probabilities=True,
    # shade_below=0.5,
    fade_below=0.25,
    stack_order="small_on_top",
    ax=ax1,
)

seqs = muts[(~muts["is_multiallelic"]) & (muts["empirical"])]["refseq"].to_list()

matrix = logomaker.transform_matrix(logomaker.alignment_to_matrix(seqs), from_type="counts", to_type="probability")
logomaker.Logo(
    matrix,
    fade_probabilities=True,
    # shade_below=0.5,
    fade_below=0.25,
    stack_order="small_on_top",
    ax=ax2,
)
ax1.set_title("Sequence context surrounding MAVs\n(at true sites)")
ax2.set_title("Sequence context surrounding MAVs\n(at random sites)")
# ax1.axhline(0.64, c="k", ls=":", label="Average GC content genome-wide")
# ax2.axhline(0.64, c="k", ls=":")
for ax in (ax1, ax2):
    sns.despine(ax=ax)
    ax.set_xticks(np.arange(0, 21, 5))
    ax.set_xticklabels(["-10", "-5", "0", "+5", "+10"])
    ax.set_xlabel("Position relative to reference nucleotide")
    ax.set_ylim(0, 1.1)
# ax1.legend()
f.tight_layout()
f.savefig("logo.png", dpi=200)

# f, ax = plt.subplots(figsize=(7, 6))
# ax.hist(stats, label="Randomly sampled SNVs", ec="w", lw=1)
# ax.axvline(ma["entropy"].mean(), label="Multiallelic SNVs", ls=":", c="k")
# ax.legend(shadow=True)
# ax.set_xlabel("Average sequence entropy\n" + "lower " + r"$\leftarrow$" + " complexity " + r"$\rightarrow$" + " higher")
# ax.set_ylabel("Count")
# p = np.sum(np.array(stats) <= ma["entropy"].mean()) / N_PERM
# ax.set_title(f"Multiallelic SNVs have similar average sequence entropy\ncompared to randomly sampled SNVs (p = {p})")
# sns.despine(ax=ax)
# f.tight_layout()
# f.savefig('entropy.png', dpi=200)
