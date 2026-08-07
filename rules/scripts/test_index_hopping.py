from snakemake.script import snakemake
import pandas as pd
from pysam import AlignmentFile

from utils import get_read_info
from collections import Counter


mutations = pd.read_csv(snakemake.input.mutations, sep="\t")
# get multiallelics
ma = mutations[mutations["is_multiallelic"]]

cram = AlignmentFile(snakemake.input.cram, "rc")

res = []
for i, row in ma.iterrows():
    d = row.to_dict()
    ref = d["ref"]

    read_info = get_read_info(
        cram,
        row["true_chrom"],
        row["start"],
        row["end"],
    )
    bases = Counter(r.base for r in read_info).most_common()
    alt_ct = 0
    if len(bases) > 1:
        alt_ct = sum(b[1] for b in bases if b[0] != ref)
    d = row.to_dict()
    d["alt_sample"] = snakemake.wildcards.SAMPLE
    alt_ev = "|".join([":".join(list(map(str, b))) for b in bases])
    d["alt_evidence"] = alt_ev
    d["alt_ct"] = alt_ct
    res.append(d)

res = pd.DataFrame(res)
res.to_csv(snakemake.output.mutations, sep="\t", index=False)
