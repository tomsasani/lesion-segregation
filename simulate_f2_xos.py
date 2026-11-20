import numpy as np
import random
import tqdm
from collections import Counter
import itertools
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def crossover(rng, chrom1, chrom2):
    """
    crossover two chroms
    """

    # get chrom1 mutation array
    a1 = chrom1
    a2 = chrom2

    chrom_len = a1.shape[0]

    # xo position
    xo_probs = np.zeros(chrom_len)
    i1, i2 = int(chrom_len * 0.1), int(chrom_len * 0.3)
    i3, i4 = int(chrom_len * 0.7), int(chrom_len * 0.9)
    xo_probs[i1:i2] = 0.5
    xo_probs[i3:i4] = 0.5
    xo_probs /= np.sum(xo_probs)

    n_xos = rng.choice([1, 2])
    n_xos = 1
    if n_xos == 1:
        i1 = rng.choice(chrom_len, p=xo_probs)

        # generate new reciprocoal chroms
        _a1 = np.concat((a1[:i1], a2[i1:]))
        _a2 = np.concat((a2[:i1], a1[i1:]))

    elif n_xos == 2:
        i1, i2 = sorted(rng.choice(chrom_len, p=xo_probs, size=2, replace=False))

        # generate new reciprocoal chroms
        _a1 = np.concat((a1[:i1], a2[i1:i2], a1[i2:]))
        _a2 = np.concat((a2[:i1], a1[i1:i2], a2[i2:]))

    return _a1, _a2


def make_chromosome(
    lesion_strand,
    rng,
    n_muts=20,
    length=1_000,
):

    # generate double-stranded "mutations"
    m = np.zeros(int(length), dtype=np.int8)
    mut_pos = rng.choice(int(length), size=n_muts, replace=False)
    m[mut_pos] = lesion_strand
    return m


def generate_offspring_chroms(rng, c1, c2, n_offspring: int = 10):

    # store diploid chromosomes of offspring
    chroms = []
    for o in range(n_offspring):
        # in each offspring, we assume that a crossover occurred between
        # the parental chromosomes with prob p
        # NOTE: this is not guaranteed to be the case.
        xo = rng.uniform() <= 0.5
        if xo:
            # make recombinant chromosomes
            _c1, _c2 = crossover(rng, c1, c2)
        else:
            _c1, _c2 = c1.copy(), c2.copy()
        # get a list of all four possible recombinants this parent could donate
        recombinants = [c1, c2, _c1, _c2]

        # randomly select two recombinants (gametes)
        idxs = rng.choice(4, replace=False, size=2)
        for i in idxs:
            chroms.append(recombinants[i])

    chroms = np.vstack(chroms)
   
    return chroms

# generate original chromosomes

rng = np.random.default_rng(42)

CHROMS = 4

f, axarr = plt.subplots(CHROMS, CHROMS, figsize=(3 * CHROMS, 3 * CHROMS), sharex=True, sharey=True)

n_offspring = 3

COL = {0: "grey", 1: "blue", 2: "red"}
vafs = []
for (chrom_i, chrom_j) in tqdm.tqdm(itertools.product(range(CHROMS), range(CHROMS))):
    
    # make two parental chromosomes, which we assume are identical in every
    # gamete progenitor (cells derived from Z2/Z2)

    c1_mut = 1 if rng.uniform() <= 0.5 else 2
    c2_mut = 1 if rng.uniform() <= 0.5 else 2

    c1 = make_chromosome(c1_mut, rng, n_muts=rng.choice([25, 50, 100, 200]), length=1e6)
    c2 = make_chromosome(c2_mut, rng, n_muts=rng.choice([25, 50, 100, 200]), length=1e6)

    offspring_chroms = generate_offspring_chroms(rng, c1, c2, n_offspring=n_offspring)
    print (offspring_chroms)
    axarr[chrom_i, chrom_j].set_title(f"{chrom_i},{chrom_j}")

    if chrom_i == 2 and chrom_j == 0:
        f2, ax = plt.subplots()
        # out = np.hstack((c1, c2, offspring_chroms))
        # print (out.shape)
        # print (np.where(out > 0))
        ii = np.where(np.sum(offspring_chroms, axis=0) > 0)[0]

        for _x in ii:
            for _y in range(n_offspring * 2):
                gt = offspring_chroms[_y, _x]
                if gt == 0: continue
                ax.scatter(_x, _y, c=COL[gt], alpha=0.5)
            
        
        f2.savefig("hm.png")

        f3, (ax1, ax2) = plt.subplots(2)
        dists = [0] + list(ii[1:] - ii[:-1])
        for _i, _x in enumerate(ii):
            gts = offspring_chroms[:, _x]
            gts = gts[gts != 0]
            assert np.unique(gts).shape[0] == 1
            gt = gts[0]
            if gt == 1:
                ax = ax1
            else:
                ax = ax2
            ax.scatter(_x, dists[_i], c=COL[gt])
        f3.savefig("dist.png")

    # subset to mutated sites
    rel_sites = np.where(np.sum(offspring_chroms, axis=0) != 0)[0]
    n_chroms = offspring_chroms.shape[0]

    for i in rel_sites:
        nucs = offspring_chroms[:, i]
        uniq_nucs = np.unique(nucs[nucs != 0])
        assert uniq_nucs.shape[0] <= 1
        vaf = np.sum(offspring_chroms[:, i] != 0) / n_chroms
        succ = rng.binomial(n=100, p=vaf)
        if vaf > 0:
            # add binomial sampling noise assuming 100 reads
            succ = rng.binomial(n=150, p=vaf)
            vafs.append(succ / 150)
            axarr[chrom_i, chrom_j].scatter(i, succ / 150, c=COL[uniq_nucs[0]])
    sns.despine(ax=axarr[chrom_i, chrom_j])
f.tight_layout()
f.savefig("vaf.png")

f, ax = plt.subplots()
ax.hist(vafs)
f.savefig("vafs.png")

