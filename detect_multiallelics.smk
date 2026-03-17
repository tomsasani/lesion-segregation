import pandas as pd
import glob

existing_bams = [fh for fh in glob.glob("data/bam/*.bam")]
existing_samples = [fh.split("/")[-1].split(".")[0] for fh in existing_bams]

# filter to relevant samples
metadata = pd.read_excel("data/41467_2020_15912_MOESM4_ESM.xlsx", 
                         sheet_name="Sample_description", 
                         skiprows=1).fillna({"Normal": "CD0001b", "Drug.concentration": 0})

mutagens = ["MMS", "EMS", "DMS"]

genotypes = [
    "N2",
    "xpa-1",
    "xpc-1",
    "xpf-1",
    "polk-1",
    "polh-1",
    "rev-3",
    "polq-1",
]

metadata = metadata[metadata["Mutagen"].isin(mutagens)]
metadata = metadata[metadata["Genotype"].isin(genotypes)]

SAMPLE2NORMAL = {}
for k,v in zip(metadata["Sample"], metadata["Normal"]):
    if "," in v:
        SAMPLE2NORMAL[k] = v.split(",")[0]
    elif " " in v:
        SAMPLE2NORMAL[k] = v.split(" ")[0]
    else:
        SAMPLE2NORMAL[k] = v


metadata = metadata[metadata["Drug.concentration"] != 0]

SAMPLES = metadata["Sample"].to_list()

rule all:
    input:
        expand("csv/multiallelics/evidence.{PERMUTATION}.tsv", PERMUTATION=range(100))

rule find_multiallelics:
    input:
        mutations = "mutations.tsv",
        exclude = "data/repeats.tsv",
    output:
        "csv/multiallelics/evidence.{PERMUTATION}.tsv"
    shell:
        """
        python find_multiallelics.py --mutations {input.mutations} \
                                     --exclude {input.exclude} \
                                     --out {output} \
                                     -permutation {wildcards.PERMUTATION} \
        """