import numpy as np
import pandas as pd
import csv
import tabix
import tqdm

from snakemake.script import snakemake

tb = tabix.open(snakemake.input.phylop_scores)

res = []

with open(snakemake.input.mutations, "rt") as infh:
    csvf = csv.reader(infh, delimiter="\t")
    header = None
    for l in tqdm.tqdm(csvf):
        if header is None:
            header = l
            continue
        d = dict(zip(header, l))
        records = tb.query(d['chrom'], int(d["start"]), int(d["end"]))
        r = [_ for _ in records]
        assert len(r) <= 1
        p = np.nan
        if len(r) == 1:
            p = float(r[0][-1])
        d.update({"phylop": p})
        res.append(d)

res = pd.DataFrame(res)

res.to_csv(snakemake.output.tsv, sep="\t", index=False)
