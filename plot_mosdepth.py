import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import glob
import numpy as np


f, ax = plt.subplots()

for fh in glob.glob("mosdepth_out/*.global.dist.txt"):
    smp = fh.split("/")[-1].split(".")[0]
    df = pd.read_csv(fh, sep="\t", names=["chrom", "depth", "frac"])
    df = df[df["chrom"] == "total"]
    df = df[df["depth"] <= 200]
    ax.plot(df["depth"], df["frac"], label=smp)
ax.legend()
f.savefig('mosdepth.png')