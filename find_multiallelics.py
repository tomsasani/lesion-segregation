import pysam
from collections import Counter, defaultdict
import numpy as np
import pandas as pd
from bx.intervals.intersection import Interval, IntervalTree
import csv
import tqdm
import os
from cyvcf2 import VCF

def get_read_info(bam, chrom, start, end):

    read_info = []
    
    ct = 0
    for read in bam.fetch(chrom, start, end):
        if read.mapping_quality < 60:
            continue
        if not read.is_paired:
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
                read_info.append((read.reference_start, read.reference_end, ct, orient, base))
                ct += 1
    return read_info


def main(args):

    exclude_tree = defaultdict(IntervalTree)
    with open(args.exclude, "rt") as infh:
        csvf = csv.reader(infh, delimiter="\t")
        for i, l in enumerate(csvf):
            if i == 0:
                header = l
                continue
            d = dict(zip(header, l))
            exclude_tree[d["genoName"]].insert_interval(Interval(int(d["genoStart"]), int(d["genoEnd"])))

    rng = np.random.default_rng(args.permutation)

    # IGV colors
    # cmap = dict(zip(["A", "C", "G", "T"], ["green", "blue", "orange", "red"]))

    # ctrl_samfile = pysam.AlignmentFile(args.ctrl_bam, "rb")

    mutations = pd.read_csv(args.mutations, sep="\t")[["sample", "chrom", "pos"]]
    if args.permutation > 0:
        mutations["sample"] = rng.permutation(mutations["sample"].values)

    # keep track of how many muts are HOM_ALT vs HET
    res = []

    for sample, sample_df in tqdm.tqdm(mutations.groupby("sample")):

        for i, row in sample_df.iterrows():
            chrom = f"chr{row['chrom']}"
            end = row["pos"]
            start = end - 1

            if len(exclude_tree[chrom].find(start, end)) > 0:
                continue

            if os.path.exists(f"data/bam/{sample}.bam"):
                samfile = pysam.AlignmentFile(f"data/bam/{sample}.bam", "rb")
            else:
                continue

            focal_read_info = get_read_info(samfile, chrom, start, end)

            base_counts = Counter([e[-1] for e in focal_read_info]).most_common()

            multiallelic = False
            if len(base_counts) > 2 and base_counts[-1][-1] > 1:
                multiallelic = True

                # f, (ax1, ax2) = plt.subplots(2, figsize=(8, 9), sharex=True, sharey=False)

                # for ri, ax in zip((ctrl_read_info, focal_read_info), (ax1, ax2)):

                #     n_reads = len(ri)
                #     if n_reads == 0:
                #         print (v.CHROM, v.start, v.end)
                #         continue
                #     thresh = 100 / n_reads

                #     rng = np.random.default_rng(42)

                #     for s, e, ct, orient, base in ri:
                #         if rng.uniform() > thresh: continue
                #         ax.plot([s, e], [ct, ct],
                #                 color="gainsboro" if orient else "darkgrey",
                #                 lw=4,
                #                 zorder=-1,
                #                 solid_capstyle='round')
                #         ax.scatter(v.start,
                #                 ct,
                #                 color=cmap[base],
                #                 zorder=0 if base == v.REF else 1)

                #     ax.set_yticks([])
                #     sns.despine(ax=ax, left=True)
                #     ax.set_xlabel("Position in ce11 (Mbp)")
                #     if ax == ax1:
                #         ax.set_title(f"Read support in control sample\n(n = {n_reads})")
                #     else:
                #         ax.set_title(f"Read support in {sample}\n(n = {n_reads})")
                #     legend = []
                #     for allele, color in cmap.items():

                #         point = Line2D(
                #             [0],
                #             [0],
                #             label=f"{allele} (REF)" if allele == v.REF else allele,
                #             marker="o",
                #             markersize=10,
                #             markerfacecolor=color,
                #             markeredgecolor="w",
                #             linestyle="none",
                #         )
                #         legend.append(point)
                #     ax.legend(handles=legend, frameon=False, fontsize=14)

                # f.tight_layout()
                # f.savefig(f"fig/multiallelic/{sample}.{v.CHROM}:{v.start}.png",
                #         dpi=200)
                # plt.close()

            base_count_summary = []
            for b, c in base_counts:
                base_count_summary.append(":".join([b, str(c)]))
            base_count_summary = "|".join(base_count_summary)

            res.append(
                {
                    "sample": sample,
                    "chrom": chrom,
                    "pos": end,
                    "is_multiallelic": multiallelic,
                    "summary": base_count_summary,
                    "permuted": args.permutation,
                }
            )

    res = pd.DataFrame(res)
    res.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":

    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--mutations")
    p.add_argument("--exclude")
    p.add_argument("--out")
    p.add_argument("-permutation", required=False, default=0, type=int)
    p.add_argument("-plot", action="store_true")
    args = p.parse_args()
    main(args)
