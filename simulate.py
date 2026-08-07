import numpy as np
import random 
import tqdm
from collections import Counter
import itertools
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


class Strand(object):
    def __init__(self, has_lesion: bool = True, base: int = 0,):
        self.has_lesion = has_lesion
        self.base = base

    def set_new_lesion_status(self, has_lesion):
        self.has_lesion = has_lesion

class Chromosome(object):

    def __init__(
        self,
        watson_strand,
        crick_strand,
        p_error: float = 0.01,
        p_repair: float = 0.01,
        
    ):
        self.watson_strand = watson_strand
        self.crick_strand = crick_strand
        self.p_error = p_error
        self.p_repair = p_repair

    def has_lesion(self):
        return self.watson_strand.has_lesion or self.crick_strand.has_lesion

    def synthesize(self):
        """
        create two new chromosomes from this single chromosome.
        """
        new_chroms = []
        for i, parent_strand in enumerate((self.watson_strand, self.crick_strand)):
            # assert not (parent_strand.has_mismatch and parent_strand.has_lesion)

            # if this strand has a lesion, decide whether it
            # persists into the daughter cell.
            if parent_strand.has_lesion:
                retains_lesion = random.random() > self.p_repair
                # if it retains its lesion, synthesize the opposite strand
                # with the possibliity of error-prone polymerase
                if retains_lesion:
                    if random.random() < self.p_error:
                        # assume that polymerase randomly selects a base to incorporate
                        child_strand = Strand(has_lesion=False, base=random.randrange(4))
                    else:
                        child_strand = Strand(has_lesion=False, base=parent_strand.base)
                # if the lesion is repaired, update the parent's lesion status
                else:
                    parent_strand.set_new_lesion_status(False)
                    child_strand = Strand(has_lesion=False, base = parent_strand.base)
            else:
                # if the parent strand doesn't have a lesion,
                child_strand = Strand(has_lesion=False, base=parent_strand.base)

            new_chrom = Chromosome(
                parent_strand,
                child_strand,
                p_error=self.p_error,
                p_repair=self.p_repair,
            )
            new_chroms.append(new_chrom)
        return new_chroms


# define the original chromosome. one strand has a lesion,
# they other doesn't.

N_TRIALS = 1_000
N_LESIONS = 100

r_rates = np.logspace(-2, 0, num=10)
e_rates = np.logspace(-2, 0, num=10)
print (r_rates)

arr = np.zeros((r_rates.shape[0], e_rates.shape[0]))

res = []
for idx, (repair_rate, error_rate) in tqdm.tqdm(
    enumerate(itertools.product(r_rates, e_rates))
):
    i, j = divmod(idx, len(r_rates))
    mav = 0
    for t in range(N_TRIALS):

        chrom = Chromosome(
            Strand(has_lesion=False, base=0),
            Strand(has_lesion=True, base=0),
            p_repair=repair_rate,
            p_error=error_rate,
        )

        # divide four times to reach P4
        for div in range(4):
            new_chroms = chrom.synthesize()
            # pick a random chrom to pass down "germline lineage"
            chrom = random.choice(new_chroms)

        # make it to P4?
        has_lesion_p4 = chrom.has_lesion()
        # save time by skipping if we don't have the lesion at P4
        if not (has_lesion_p4 or chrom.watson_strand.base > 0 or chrom.crick_strand.base > 0):
            # res.append(
            #             {
            #                 "trial": t,
            #                 "repair_rate": repair_rate,
            #                 "error_rate": error_rate,
            #                 "uniq_bases": 0,
            #                 "has_lesion_p4": has_lesion_p4,
            #             }
            #         )
            continue
        # # then, divide to get Z2/Z3 and future germ cell population
        pool = chrom.synthesize()
        for div in range(5):
            new_pool = []
            for c in pool:
                daughters = c.synthesize()
                new_pool.extend(daughters)
            pool = new_pool

        bases = []
        for chrom in pool:
            bases.extend([chrom.watson_strand.base, chrom.crick_strand.base])

        # randomly sample 50X coverage, assuming other 50X comes from mom
        random.shuffle(bases)
        bases = bases[:50]

        n_uniq = len(set([b for b in bases if b != 0]))

        if n_uniq > 2:
            mav += 1

    prob = mav / N_TRIALS
    true_prob = 1 - ((1 - prob) ** N_LESIONS)
    print (prob, true_prob)
    arr[i, j] = true_prob


f, ax = plt.subplots()
sns.heatmap(arr, ax=ax)
ax.set_xlabel("TLS error rate")
ax.set_ylabel("Error-free repair rate")
# ax.set_xticklabels(r_rates, rotation=45)
# ax.set_yticklabels(e_rates)
f.savefig("o.png", dpi=200)

# g = sns.jointplot(data=res, x="repair_rate", y="error_rate", )
# print (res.shape)
# print (res.groupby(["repair_rate", "error_rate", "uniq_bases"]).size())
# print (res.groupby(["repair_rate", "error_rate"]).agg(has_lesion=("has_lesion_p4", "sum")))
# g = sns.FacetGrid(data=res, row="repair_rate", col="error_rate")
# g.map(sns.histplot, "uniq_bases", bins=np.arange(5))
# g.set(yscale="log")
# g.savefig("o.png", dpi=200)
