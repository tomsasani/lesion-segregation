import pysam
from collections import Counter, defaultdict, namedtuple
import numpy as np
import pandas as pd
from bx.intervals.intersection import Interval, IntervalTree
import csv
import tqdm
import os
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
import scipy.stats as ss

plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

from snakemake.script import snakemake


def plot(
    focal_read_info,
    ctrl_bam_fh: str = "data/cram/CD0803c.cram",
    fig_fh: str = "o.png",
):

    f, (ax1, ax2) = plt.subplots(
        2,
        figsize=(8, 10),
        sharex=True,
        sharey=True,
    )

    cmap = dict(
        zip(
            ["A", "C", "G", "T"],
            ["green", "blue", "orange", "red"],
        )
    )

    ctrl_samfile = pysam.AlignmentFile(ctrl_bam_fh, "rc")
    # ctrl_read_info = get_read_info(ctrl_samfile, chrom, start, end)

    for ri, ax in zip((focal_read_info, focal_read_info), (ax1, ax2)):
        base2ct = dict(zip(["A", "C", "G", "T"], [0] * 4))

        n_reads = len(ri)
        if n_reads == 0:
            continue
        thresh = 50 / n_reads

        rng = np.random.default_rng(42)

        read_count = 0
        for read in ri:

            s, e = read.start, read.end
            orient = read.orientation
            base = read.base

            base2ct[base] += 1

            if rng.uniform() > thresh: continue

            ax.plot([s, e], [read_count, read_count],
                    color="gainsboro" if orient else "darkgrey",
                    lw=4,
                    zorder=-1,
                    solid_capstyle='round')
            ax.scatter(start, read_count, color=cmap[base], zorder=1)

            read_count += 1

        ax.set_yticks([])
        sns.despine(ax=ax, left=True)
        ax.set_xlabel(f"Position on {chrom}")
        if ax == ax1:
            ax.set_title(
                f"Read support in control sample\n(n = {n_reads})")
        else:
            ax.set_title(
                f"Read support in {snakemake.wildcards.SAMPLE}\n(n = {n_reads})"
            )

        legend = []
        for allele, color in cmap.items():

            point = Line2D(
                [0],
                [0],
                label=f"{allele} ({base2ct[allele]})",
                marker="o",
                markersize=10,
                markerfacecolor=color,
                markeredgecolor="w",
                linestyle="none",
            )
            legend.append(point)
        ax.legend(handles=legend,
                    frameon=False,
                    fontsize=14,
                    title="Allele support")

    f.tight_layout()
    f.savefig(fig_fh, dpi=200)
    plt.close()


def get_read_info(bam, chrom, start, end):

    res = []

    ReadInfo = namedtuple(
        "read",
        [
            "start",
            "end",
            "orientation",
            "base",
            "qual",
        ],
    )

    ct = 0
    for ri, read in enumerate(bam.fetch(chrom, start, end)):
        # if we've counted over 150, move on
        if ri > 150: break
        if read.mapping_quality < 60:
            continue
        if not read.is_paired:
            continue
        if read.is_duplicate:
            continue
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mate_is_unmapped:
            continue
        if read.reference_id != read.next_reference_id:
            continue
        if read.reference_start <= start < read.reference_end:
            dist_from_start = start - read.reference_start
            dist_from_end = read.reference_end - end  # reference_end is exclusive
            is_near_boundary = (dist_from_start < 5) or (dist_from_end <= 5)
            if is_near_boundary: continue
        # get aligned pairs
        pairs = read.get_aligned_pairs(matches_only=True)
        for query_pos, ref_pos in pairs:
            if ref_pos == start:
                # query_pos is the index in the read's sequence
                base = read.query_sequence[query_pos]
                qual = read.query_qualities[query_pos]
                orient = read.is_forward
                if qual < 20: continue

                read_info = ReadInfo(
                    start=read.reference_start,
                    end=read.reference_end,
                    orientation=orient,
                    base=base,
                    qual=qual,
                )
                res.append(read_info)

    return res


n_perm = int(snakemake.wildcards.PERM)

mutations = pd.read_csv(snakemake.input.mutations, sep="\t")


res = []

# read in the input bam
samfile = pysam.AlignmentFile(snakemake.input.bam, "rb" if snakemake.input.bam.endswith("bam") else "rc")
# NOTE: should we be sampling null mutations
# only from the same mutagen/genotype combinations?
for _perm in tqdm.tqdm(range(n_perm + 1)):
    rng = np.random.default_rng(_perm)

    if _perm > 0:

        # first, figure out how many mutations this sample has
        n = mutations[mutations["sample"] ==
                      snakemake.wildcards.SAMPLE].shape[0]
        # permute sample labels instead
        # mutations["sample"] = rng.permutation(mutations["sample"].values)
        # then, remove this sample from the overall dataframe
        sans = mutations[mutations["sample"] != snakemake.wildcards.SAMPLE]
        # ensure we're looking at samples from same treatment regime
        # sans = sans[sans["sample"].isin(snakemake.params.shared_treatment)]
        # randomly sample mutations from the other strains
        sample_mutations = sans.sample(n=n, replace=False)#, random_state=rng)
        sample_mutations["sample"] = snakemake.wildcards.SAMPLE
        # sample_mutations = mutations[mutations["sample"] == snakemake.wildcards.SAMPLE]
    else:
        sample_mutations = mutations[mutations["sample"] ==
                                     snakemake.wildcards.SAMPLE]

    # loop over mutations

    for it in sample_mutations.itertuples(index=False):

        chrom, start, ref = it.true_chrom, it.start, it.ref

        focal_read_info = get_read_info(samfile, chrom, start, start + 1)

        # total number of reads cannot be outside thresholds
        n_reads = len(focal_read_info)

        # figure out the strand bias of each base's reads
        # bias = []
        # is_biased = False
        # for base in ("A", "T", "C", "G"):
        #     # only measure strand bias on ALT bases
        #     if base == ref: continue
        #     base_read_strand = [
        #         r.orientation for r in focal_read_info if r.base == base
        #     ]
        #     if len(base_read_strand) == 0:
        #         bias.append(f"{base}:0:0")
        #     else:
        #         fwd = sum(base_read_strand)
        #         rev = len(base_read_strand) - fwd

        #         binom = ss.binomtest(fwd, fwd + rev, p=0.5)
        #         if binom.pvalue < 0.05:
        #             is_biased = True
        #         bias.append(f"{base}:{fwd}:{rev}")
        # bias = "|".join(bias)

        base_counts = Counter([r.base for r in focal_read_info]).most_common()

        if (_perm == 0 and snakemake.params.plot is True
                and len(base_counts) > 2 and base_counts[-1][-1] > 1):

            plot(
                focal_read_info,
                fig_fh=
                f"fig/multiallelic/{snakemake.wildcards.SAMPLE}.{chrom}:{start}.png"
            )

        base_count_summary = []
        for b, c in base_counts:
            base_count_summary.append(":".join([b, str(c)]))
        base_count_summary = "|".join(base_count_summary)

        is_multiallelic = False
        reads_supporting_ma = 0
        if _perm == 0:
            if len(base_counts) > 2:
                reads_supporting_ma = min([c for b, c in base_counts])
        else:
            if len(base_counts) > 1:
                reads_supporting_ma = min([c for b, c in base_counts])

        row_dict = it._asdict()
        row_dict.update({
            "summary": base_count_summary,
            "permuted": _perm,
            # "is_strand_biased": is_biased,
            "ma_reads": reads_supporting_ma,
        })
        res.append(row_dict)


res = pd.DataFrame(res)
res.to_csv(snakemake.output.mutations, sep="\t", index=False)
