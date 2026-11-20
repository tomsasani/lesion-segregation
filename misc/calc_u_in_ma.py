import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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

dnms = pd.read_excel(
    "/scratch/ucgd/lustre-labs/quinlan/u1006375/lesion-segregation/data/41467_2020_15912_MOESM4_ESM.xlsx",
    sheet_name="Mutation counts for all samples",
    skiprows=1,
).rename(columns={"Unnamed: 0": "Sample"})

print (dnms)

dnms["total"] = np.sum(dnms.values[:, 1:97], axis=1)
print (dnms)
dnms = dnms[["Sample", "total"]]
relevant_samples = metadata[(metadata["Type"] == "mut.acc.") & (metadata["Generation"] != 0)]


relevant_samples = relevant_samples.merge(dnms)
print (relevant_samples)
relevant_samples["rate"] = relevant_samples["total"] / relevant_samples["Generation"]

n2 = relevant_samples[relevant_samples["Genotype"] == "N2"]["rate"].mean()

relevant_samples["true_rate"] = relevant_samples["rate"] / n2

f, ax = plt.subplots(figsize=(6, 20))
sns.stripplot(data=relevant_samples, y="Genotype", x="true_rate", hue="Genotype", ax=ax)
ax.set_xlabel("Mutation rate excess above wild-type")
sns.despine(ax=ax)
f.tight_layout()
f.savefig('rates.png', dpi=200)