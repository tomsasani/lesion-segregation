import pandas as pd
from cyvcf2 import VCF
import numpy as np
import ruptures as rpt
from global_vars import MUT2EXP, REVCOMP

SYMMETRY_THRESH = 0.3
N_MUTATIONS = 20

mutations = pd.read_csv("strand_bias.tsv", sep="\t")

chrom2idx = dict(zip(["I", "II", "III", "IV", "V", "X"], range(6)))

res = []
for (sample, genotype, mutagen, dose), sdf in mutations.groupby(["sample", "genotype", "mutagen", "dose"]):

    for ci, (chrom, chrom_df) in enumerate(sdf.groupby("chrom")):

        chrom_df = chrom_df.sort_values("pos")

        pos = chrom_df["pos"].values
        dist = chrom_df["dist"].values
        strands = chrom_df["strand"].values

        if strands.shape[0] > 1:
           
            algo = rpt.Pelt(model="rbf").fit(strands)
            cur = 1
            for bkpt in algo.predict(pen=np.log(strands.shape[0])):

                # calculate symmetry score in region
                # NOTE: here we only use the mutations known to be on
                # the watson (0) or crick (1) strand
                region_arr = strands[cur:bkpt]
                n_w, n_c = np.sum(region_arr == 0), np.sum(region_arr == 1)
                symmetry = (n_w - n_c) / (n_w + n_c)

                # ensure region passes symmetry threshold
                expected_strand = None
                if symmetry > SYMMETRY_THRESH:
                    expected_strand = 0
                elif symmetry < -SYMMETRY_THRESH:
                    expected_strand = 1
                
                if expected_strand is None: continue
                
                # NOTE: num inforamtive muations is combo of nw and nc
                n_inf_mutations = n_w + n_c
                if n_inf_mutations < N_MUTATIONS: continue
               
                seg_start = pos[cur - 1]
                seg_end = pos[bkpt - 1]

                # gather relevant mutations
                mutations_in_seg = chrom_df[chrom_df["pos"].between(seg_start, seg_end, inclusive="both")]
                mutations_in_seg["segment_label"] = expected_strand
                # filter to mutations that match symmetry assignment
                # NOTE: not currently doing this so we can use for collateral mutagenesis
                # mutations_in_seg = mutations_in_seg[mutations_in_seg["strand"] == expected_strand]
                # adjust the distance measure
                dists = [0]
                dists.extend(list(mutations_in_seg["pos"].values[1:] - mutations_in_seg["pos"].values[:-1]))
                mutations_in_seg["dist"] = dists
                res.append(mutations_in_seg)

res = pd.concat(res)
res.to_csv("strand_resolved.tsv", sep="\t", index=False)

                


                