import pandas as pd
from cyvcf2 import VCF
import seaborn as sns
import tqdm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import ruptures as rpt
import numpy as np
from collections import defaultdict
from scipy.stats import norm, binom

mutations = pd.read_csv("/scratch/ucgd/lustre-labs/quinlan/u1006375/lesion-segregation/csv/mutations/combined.tsv", sep="\t")

metadata = pd.read_excel(
    "/scratch/ucgd/lustre-labs/quinlan/u1006375/lesion-segregation/data/41467_2020_15912_MOESM4_ESM.xlsx",
    sheet_name="Sample_description",
    skiprows=1,
).fillna(
    {
        "Normal": "CD0001b",
        "Mutagen": "none",
        "Drug.concentration": -1,
    }
)

metadata["Mutagen"] = metadata.apply(
    lambda row: row["Mutagen"] if row["Type"] != "mut.acc." else "None (M.A.)", axis=1
)

mutations = mutations.merge(metadata, left_on="sample", right_on="Sample")

mutations = mutations[mutations["Mutagen"].isin(["DMS", "MMS", "EMS"])]

mavs = pd.read_csv("csv/multiallelics.10.original.filtered.tsv", sep="\t")

# mutations = mutations[mutations["sample"].isin(mavs["sample"].to_list())]

# mutations = mutations[mutations["sample"] == "CD0496d"]
f, ax = plt.subplots()
ax.hist(mutations[mutations["Mutagen"] == "EMS"]["af"], bins=np.arange(0, 1, 0.02))
f.savefig("aaf.png")
g = sns.FacetGrid(data=mutations, col="Mutagen", sharey=False)
g.map(sns.histplot, "af", bins=np.arange(0, 1.02, 0.02))
g.savefig("af.png")