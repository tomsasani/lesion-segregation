import hashlib
from collections import Counter
import sys

import numpy as np
import pandas as pd
import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
import scipy.stats as ss
from pysam import AlignmentFile

from utils import get_read_info

from snakemake.script import snakemake


plt.rc("font", size=14)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]


def plot(
    focal_read_info,
    ctrl_read_info,
    v_start,
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
            ["A", "C", "G", "T", "N"],
            sns.color_palette("colorblind", 5),
        )
    )

    for ri, ax in zip((ctrl_read_info, focal_read_info), (ax1, ax2)):
        fwd_base2ct = {base: 0 for base in "ACGTN"}
        rev_base2ct = {base: 0 for base in "ACGTN"}

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
            if orient:
                fwd_base2ct[base] += 1
            else:
                rev_base2ct[base] += 1

            if rng.uniform() > thresh:
                continue

            adj_s, adj_e = s - v_start, e - v_start
            ax.plot(
                [adj_s, adj_e],
                [read_count, read_count],
                color="gainsboro" if orient else "dimgrey",
                lw=4,
                zorder=-1,
                solid_capstyle="round",
            )
            ax.scatter(0, read_count, color=cmap[base], zorder=1)

            read_count += 1

        ax.set_yticks([])
        sns.despine(ax=ax, left=True)
        ax.set_xlabel("Position relative to mutated base (bp)")
        if ax == ax1:
            ax.set_title(f"Read support in control sample")
        else:
            ax.set_title(
                f"Read support in {snakemake.wildcards.SAMPLE}"
            )

        legend = []
        for allele, color in cmap.items():
            for di, d in enumerate((fwd_base2ct, rev_base2ct)):
                if d[allele] == 0:
                    continue
                point = Line2D(
                    [0],
                    [0],
                    label=f"{allele} ({d[allele]})",
                    marker="o",
                    markersize=10,
                    markerfacecolor=color,
                    markeredgecolor="w",
                    linewidth=3,
                    color="gainsboro" if di == 0 else "dimgrey",
                )
                legend.append(point)
        ax.legend(
            handles=legend,
            frameon=False,
            fontsize=14,
            title="Allele support",
            loc="upper left",
        )

    f.tight_layout()
    f.savefig(fig_fh, dpi=200)
    plt.close()


n_perm = int(snakemake.wildcards.PERM)

mutations = pd.read_csv(snakemake.input.mutations, sep="\t")

# first, figure out how many mutations this sample has
focal = mutations[mutations["sample"] == snakemake.wildcards.SAMPLE]

n = focal.shape[0]
if n == 0:
    pd.DataFrame().to_csv(snakemake.output.mutations, sep="\t", index=False)
    sys.exit()

genotype = focal["genotype"].iloc[0]
mutagen = focal["mutagen"].iloc[0]
conc = focal["dose"].iloc[0]

res = []

for _perm in tqdm.tqdm(range(n_perm + 1)):

    # get a random number generator that's unique to this
    # sample and perm combination
    hex_hash = hashlib.sha256(snakemake.wildcards.SAMPLE.encode('utf-8')).hexdigest()
    hash_integer = int(hex_hash, 16)
    short_hash = hash_integer % 10**8
    rng = np.random.default_rng(short_hash + _perm)

    # read in control samfile
    ctrl_samfile = AlignmentFile(snakemake.input.ctrl_bam, "rc")

    # NOTE: alternative approach to random sampling. in each permutation,
    # use a different strain's BAM file.
    bam_smp = snakemake.wildcards.SAMPLE
    samfile = AlignmentFile(snakemake.input.bam, "rc")
    if _perm > 0:# and snakemake.wildcards.STRATEGY == "by_bam":
        idx = rng.integers(0, len(snakemake.params.samples))
        bam_smp = snakemake.params.samples[idx]
        bam_fh = f"data/cram/{bam_smp}.cram"
        samfile = AlignmentFile(bam_fh, "rc")
       
    sample_mutations = mutations[mutations["sample"] == snakemake.wildcards.SAMPLE].copy()
    sample_mutations.loc[:, "bam_smp"] = bam_smp

    # loop over mutations
    for it in sample_mutations.itertuples(index=False):

        chrom, start, ref, alt = it.chrom, it.start, it.ref, it.alt

        focal_read_info = get_read_info(samfile, chrom, start, start + 1)

        # total number of reads cannot be outside thresholds
        n_reads = len(focal_read_info)

        # figure out the strand bias of each base's reads
        strand_tallies = {b: {'fwd': 0, 'rev': 0} for b in ("A", "C", "G", "T")}
        base_counts_dict = Counter()

        # Just ONE pass over the reads
        for r in focal_read_info:
            b = r.base
            base_counts_dict[b] += 1
            if b in strand_tallies:
                if r.orientation:
                    strand_tallies[b]['fwd'] += 1
                else:
                    strand_tallies[b]['rev'] += 1

        base_counts = base_counts_dict.most_common()

        # Then calculate bias without looping over the reads again
        bias = []
        is_biased = False
        for base in ("A", "T", "C", "G"):
            if base == ref:
                continue
            fwd = strand_tallies[base]['fwd']
            rev = strand_tallies[base]['rev']
            
            if fwd + rev == 0:
                bias.append(f"{base}:0:0")
            else:
                binom = ss.binomtest(fwd, fwd + rev, p=0.5)
                if binom.pvalue < 0.05:
                    is_biased = True
                bias.append(f"{base}:{fwd}:{rev}")

        base_counts = Counter([r.base for r in focal_read_info]).most_common()

        if (
            _perm == 0
            and snakemake.params.plot is True
            and len(base_counts) > 2
            and base_counts[-1][-1] >= 10
        ):

            ctrl_read_info = get_read_info(ctrl_samfile, chrom, start, start + 1)

            plot(
                focal_read_info,
                ctrl_read_info,
                start,
                fig_fh=f"fig/multiallelic/{snakemake.wildcards.SAMPLE}.{chrom}:{start}.png",
            )

        base_count_summary = []
        for b, c in base_counts:
            base_count_summary.append(":".join([b, str(c)]))
        base_count_summary = "|".join(base_count_summary)

        reads_supporting_ma = 0
        if _perm == 0:
            if len(base_counts) > 2:
                reads_supporting_ma = sum(c for b, c in base_counts if b not in (ref, alt))
        else:
            if len(base_counts) > 1:
                reads_supporting_ma = sum(c for b, c in base_counts if b != ref)

        row_dict = it._asdict()
        row_dict.update(
            {
                "summary": base_count_summary,
                "permuted": _perm,
                "is_strand_biased": is_biased,
                "bias": "|".join(bias),
                "ma_reads": reads_supporting_ma,
            }
        )
        res.append(row_dict)

    samfile.close()
res = pd.DataFrame(res)
res.to_csv(snakemake.output.mutations, sep="\t", index=False)
