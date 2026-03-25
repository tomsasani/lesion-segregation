import pysam
from collections import Counter, defaultdict
import numpy as np
import pandas as pd
from bx.intervals.intersection import Interval, IntervalTree
import csv
import tqdm
import os


def get_read_info(bam, chrom, start, end):

    read_info = defaultdict(list)

    for ri, read in enumerate(bam.fetch(chrom, start, end)):
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
        # get aligned pairs
        pairs = read.get_aligned_pairs(matches_only=True)
        for query_pos, ref_pos in pairs:
            if ref_pos in (start, end):
                # query_pos is the index in the read's sequence
                base = read.query_sequence[query_pos]
                qual = read.query_qualities[query_pos]
                orient = read.is_forward
                if qual < 20: continue
                read_info[read.query_name].append(
                    (
                        read.reference_start,
                        read.reference_end,
                        orient,
                        base,
                    )
                )

    return read_info



mutations = pd.read_csv("csv/mutations/combined.ce10.tsv", sep="\t")


res = []
for sample, sample_df in tqdm.tqdm(mutations.groupby("sample")):
    try:
        samfile = pysam.AlignmentFile(f"data/cram/{sample}.cram", "rc")
    except:
        continue
    # except FileNotFoundError:
    #     samfile = pysam.AlignmentFile(f"data/cram/{sample}.bam", "rb")

    for chrom, chrom_df in sample_df.groupby("chrom"):
        pos = chrom_df["start"].values
        # filter to low dist
        dist = [0] + list(pos[1:] - pos[:-1])
        chrom_df["dist"] = dist
        chrom_df = chrom_df.reset_index()
        for i, row in chrom_df.iterrows():
            if row["dist"] < 1000 and row["dist"] != 0:
                prev, cur = chrom_df.iloc[i - 1]["start"], row["start"]
                info = get_read_info(samfile, chrom, prev, cur)
                combos = []
                for r, v in info.items():
                    if len(v) == 1: continue
                    combos.append((v[0][-1], v[1][-1]))
                uniq = Counter(combos).most_common()
                d1 = row.to_dict()
                d2 = chrom_df.iloc[i - 1].to_dict()
                if len(uniq) < 3:
                    # d1.update({"combinatorial": False})
                    # d2.update({"combinatorial": False})
                    continue
                else:
                    if uniq[-1][-1] < 2: continue
                    # if sample == "CD0842d":
                    #     print(chrom, prev, cur, uniq)
                    #     print ([k for k,v in info.items() if len(v) > 1])
                    d1.update({"combinatorial": True})
                    d2.update({"combinatorial": True})
                res.append(d1)
                res.append(d2)

res = pd.DataFrame(res).drop(columns=["index"])
res.to_csv("combinatorial.tsv", sep="\t", index=False)
