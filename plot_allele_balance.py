import pandas as pd
from cyvcf2 import VCF
import seaborn as sns
import tqdm
import matplotlib.pyplot as plt
from global_vars import REVCOMP
import matplotlib.patches as mpatches
import ruptures as rpt
import numpy as np
from collections import defaultdict
from scipy.stats import norm, binom

mutagens = ["MMS", "EMS", "DMS"]#, "Aflatoxin-B1", "AristolochicAcid"]

genotypes = [
   "N2", "polk-1", "xpa-1", "xpc-1", "xpa-1", "agt-1"
    
]

mutations = pd.read_csv("csv/mutations/combined.tsv", sep="\t")

mutations = mutations[mutations["mutagen"].isin(mutagens)]
mutations = mutations[mutations["genotype"].isin(genotypes)]

RECOMB_BOUNDS = {
    "chrI": ((527, 3858), (11040, 14875)),
    "chrII": ((306, 4879), (12020, 14609)),
    "chrIII": ((494, 3722), (10340, 13217)),
    "chrIV": ((720, 3896), (12970, 16712)),
    "chrV": ((643, 5897), (16550, 20337)),
    "chrX": ((572, 6137), (12480, 16417)),
}

rng = np.random.default_rng(42)

for (genotype, mutagen), sdf in tqdm.tqdm(mutations.groupby(["genotype", "mutagen"])):

    sdf = sdf[sdf["dose"] != 0]

    n = sdf["sample"].nunique()
    nc = 6

    f, axarr = plt.subplots(n, nc, figsize=(16, 2.5 * n), sharex="col", sharey=True)
    f2, ax2 = plt.subplots()

    res = []
    for si, ((dose, sample), ddf) in enumerate(sdf.groupby(["dose", "sample"])):    

        B2C = {"A": "red", "T": "blue", "C": "green", "G": "gray"}

        block_means = defaultdict(list)
        block_sizes = []

        for ci, (chrom, chrom_df) in enumerate(ddf.groupby("chrom")):

            starts = chrom_df["start"].values
            balances = chrom_df["af"].values

            if chrom_df.shape[0] < 10: continue

            folded_balances = np.max(np.vstack((balances, 1 - balances)), axis=0)

            algo = rpt.Pelt(model="rbf").fit(folded_balances)
            cur = 1
            for bkpt in algo.predict(pen=np.log(folded_balances.shape[0])):
                folded_balances_block = folded_balances[cur:bkpt]
                folded_balances_block_het = folded_balances_block[folded_balances_block < 0.9]

                for t in range(100):
                    folded_balances_block_het = rng.choice(
                        folded_balances_block_het,
                        size=folded_balances_block_het.shape[0],
                        replace=True,
                    )

                    block_means[t].append(np.mean(folded_balances_block_het))
                block_sizes.append(folded_balances_block_het.shape[0])
                axarr[si, ci].plot(starts[cur:bkpt], np.ones(bkpt - cur), lw=2)
                # axarr[si, ci].hist(balances[cur:bkpt], bins=np.arange(0, 1, 0.02))
                cur = bkpt
            axarr[si, ci].scatter(
                starts,
                folded_balances,
                c="darkgrey",
                # label=ref,
                zorder=1,
            )

            # bounds = RECOMB_BOUNDS[f"{chrom}"]
            # for coords in bounds:
            #     x1, x2 = coords
            #     x1 *= 1_000
            #     x2 *= 1_000
            #     rect = mpatches.Rectangle(
            #         (x1, 0),
            #         x2 - x1,
            #         1,
            #         color="firebrick",
            #         alpha=0.1,
            #         zorder=-1,
            #         ec=None,
            #     )

            #     axarr[si, ci].add_patch(rect)

            axarr[si, ci].set_xlabel(f"Position on {chrom}")
            axarr[si, ci].set_ylabel("Variant allele freq.")
            # axarr[si, ci].legend()
            axarr[si, ci].set_title(f"{mutagen} @ {dose}\nin {sample}")
            sns.despine(ax=axarr[si, ci])

        # 2. Get Binomial weights for each k
    
    
    # 3. For each block, calculate the probability density across all possible k peaks
    
        
        block_sizes = np.array(block_sizes)
        for t in range(100):

            trial_block_means = np.array(block_means[t])
            if trial_block_means.shape[0] == 0: continue

            num_alleles = 20
            ks = np.arange(1, num_alleles) # k from 1 to 2N-1
            theoretical_vafs = ks / num_alleles
            weights = binom.pmf(ks, num_alleles, 0.5)
            weights /= weights.sum()  # Re-normalize because we excluded k=0 and k=2N

            log_likelihood = 0
            for i, v_obs in enumerate(trial_block_means):
                # Probability of this VAF given the mixture of possible k-states
                p_v = np.sum(weights * norm.pdf(v_obs, loc=theoretical_vafs, scale=0.03))
                
                # Weight the log-likelihood by the block size (SNV count or length)
                log_likelihood += block_sizes[i] * np.log(max(p_v, 1e-12))
                



            # theoretical_vafs = np.array([k / num_alleles for k in range(1, num_alleles)])
            # # For each observed block mean, find the distance to the nearest theoretical peak
            # # This uses broadcasting to find absolute differences
            # distances = np.abs(trial_block_means[:, np.newaxis] - theoretical_vafs)
            # min_sq_distances = np.min(distances, axis=1) ** 2
            # # Calculate a total error score (Mean Squared Error) for this N
            # wmse = np.average(min_sq_distances, weights=block_sizes)
            
            res.append({"trial": t, "survivors": n, "error": log_likelihood, "sample": sample, "dose": float(dose)})

    res = pd.DataFrame(res)
    sns.lineplot(data=res, x = "survivors", y= "error", hue="dose", ax=ax2)
    ax2.legend()
    f2.savefig(f"fig/ab/{mutagen}.{genotype}.survivors.png", dpi=200)
    f.tight_layout()
    f.savefig(f"fig/ab/{mutagen}.{genotype}.ab.png", dpi=200)
    plt.close()
