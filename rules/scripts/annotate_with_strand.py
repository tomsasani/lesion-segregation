import pandas as pd
from cyvcf2 import VCF

from global_vars import MUT2EXP, REVCOMP

from snakemake.script import snakemake


mutations = pd.read_csv(snakemake.input.mutations, sep="\t")

mutagen = mutations["mutagen"].unique()
assert len(mutagen) == 1
mutagen = mutagen[0]

 # figure out the expected "watson" and "crick" strands.
exp_watson = MUT2EXP[mutagen]
exp_crick = REVCOMP[exp_watson]

res = []
for chrom, chrom_df in mutations.groupby("chrom"):

    pos = chrom_df["pos"].values
    dist = [0] + list(pos[1:] - pos[:-1])

    chrom_df_with_strand = []

    for i, row in chrom_df.iterrows():
        d = row.to_dict()
        # determine the strand of this mutation.
        # 0 if the REF matches the watson starnd, 1 if it matches crick
        strand = 0 if row["ref"] == exp_watson else 1 if row["ref"] == exp_crick else -1
        d.update({"strand": strand})
        chrom_df_with_strand.append(d)
    
    chrom_df_with_strand = pd.DataFrame(chrom_df_with_strand)
    chrom_df_with_strand["dist"] = dist
    res.append(chrom_df_with_strand)

res = pd.concat(res)

res.to_csv(snakemake.output.mutations, sep="\t", index=False)
