import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import math


# n trials
T = 1

mutations = pd.read_csv("csv/mutations/combined.tsv", sep="\t")
counts = mutations.groupby("sample").size().reset_index().rename(columns={0: "count"})

mavs = pd.read_csv("csv/multiallelics.filtered.tsv", sep="\t")
print (mavs.groupby(["Genotype"]).size().sort_values())

S = counts["sample"].to_list()
M = counts["count"].to_list()
N = [300]
u = 1e-8

res = []
exp = 0
for s, m in zip(S, M):
    n_mav = mavs[mavs["sample"] == s].shape[0]
    for n in N:
        # actual prob
        p = 1 - ((1 - u) ** (m * n))

        # exp successes
        
        exp_succ = (m * n * u)
        exp += exp_succ
        

        res.append({"sample": s, "p": p, "m": m, "n_mav": n_mav})
res = pd.DataFrame(res)
print (exp)
f, axarr = plt.subplots(3, figsize=(6, 12))
sns.histplot(data=res, x="p", ax=axarr[0])
sns.histplot(data=res, x="m",  ax=axarr[1])
sns.regplot(
    data=res, x="m", y="n_mav", ax=axarr[2], scatter_kws={"alpha": 0.2}, lowess=True
)
axarr[0].set_yscale("log")
axarr[1].set_yscale("log")
f.tight_layout()
f.savefig("p.png")


M = [100, 200, 500, 1_000]
N = [100, 200, 300]
L = [100_000]
u = 1e-6

rng = np.random.default_rng(42)

res = []

for m, l, n in itertools.product(M, L, N):
    # actual prob
    p = 1 - ((1 - u) ** (m * n))

    for t in range(T):
        # keep track of overlaps across offspring
        n_overlaps = 0
        for i in range(n):
            # previously mutated sites
            prev = rng.choice(l, size=m, replace=False)
            # newly mutated sites
            new = rng.choice(l, size=int(l*u), replace=False)
            overlap = np.intersect1d(prev, new)
            n_overlaps += overlap.shape[0]
        res.append(
            {
                "t": t,
                "M": m,
                "L": l,
                "N": n,
                "exp": p,
                "emp": n_overlaps > 0,
            }
        )

res = pd.DataFrame(res)

# actual prob
# print (res)


f, ax = plt.subplots()
# sns.pointplot(data=res, x="M", hue="L", y="emp", ax=ax)
sns.pointplot(data=res, x="M", y="exp", hue="N", palette="colorblind", ax=ax)
ax.set_ylabel("Probability of at least one F2 with a 'second hit'")
ax.set_xlabel("Number of mutations in F1")
sns.despine(ax=ax)
f.savefig('sim.png', dpi=200)
