import pandas as pd
import glob
import numpy as np

# get ENA IDs and FTP links
manifest = pd.read_csv("data/combined_file_manifest.tsv", sep="\t")

SAMPLE2ACC = dict(zip(manifest["Sample / DONOR ID"], manifest["run_accession"]))

# filter to relevant samples
metadata = pd.read_excel("data/41467_2020_15912_MOESM4_ESM.xlsx", 
                         sheet_name="Sample_description", 
                         skiprows=1).fillna({"Normal": "CD0001b", "Mutagen": "none", "Drug.concentration": 0})

metadata = metadata[~metadata["Sample"].str.contains("d2")]

ALL_SAMPLES = []
for fh in glob.glob("data/Filtered_VCFS_recovered/SNV/*.vcf.gz"):
    s = fh.split("/")[-1].split("_")[0]
    if "d2" in s: continue
    ALL_SAMPLES.append(s)

# map samples to their corresponding normal sample
SAMPLE2NORMAL = {}
for k,v in zip(metadata["Sample"], metadata["Normal"]):
    if "," in v:
        SAMPLE2NORMAL[k] = v.split(",")[0]
    elif " " in v:
        SAMPLE2NORMAL[k] = v.split(" ")[0]
    else:
        SAMPLE2NORMAL[k] = v

SAMPLE2GT = dict(zip(metadata["Sample"], metadata["Genotype"]))

configfile: 'config.yaml'
mutagens = config["mutagens"]

SAMPLES_MA = metadata[(metadata["Type"] == "mut.acc.") & (metadata["Generation"] > 0)]["Sample"].unique()

# filter to samples with drug concentration > 0
metadata = metadata[metadata["Drug.concentration"] != 0]
metadata = metadata[metadata["Mutagen"].isin(mutagens)]
metadata = metadata[metadata["Generation"] == 1]


CHROMS = ["I", "II", "III", "IV", "V", "X"]
CHROMS = [f"chr{c}" for c in CHROMS]

wildcard_constraints:
    SAMPLE = r"CD\d{4,5}[a-s]",

SAMPLES = metadata["Sample"].to_list()

include: "rules/fastq2bam.smk"
include: "rules/variant_calling.smk"
include: "rules/downsampling.smk"
include: "rules/generate_figures.smk"

rule all:
    input:
        "fig/mutation_counts.png",
        "fig/mosdepth.png",
        "fig/phylop.png",
        expand("fig/{plot}.100.original.png", plot = ["n_reads", "n_reads_hist", "elbow", "mav_agg", "mav_grouped", "entropy"])
        

rule catalog_mutations:
    input:
        vcf = "data/Filtered_VCFS_recovered/SNV/{SAMPLE}_filtered_substitutions.vcf.gz",
        metadata = "data/41467_2020_15912_MOESM4_ESM.xlsx",
        exclude = "data/ce11.sdust.bed.gz",
        ref = "data/ref/ce11.fa.gz"
    output:
        mutations = "csv/per-sample/{SAMPLE}.mutations.tsv"
    params:
        var_source = "caveman"
    resources:
        runtime = 5
    script:
        "rules/scripts/catalog_mutations.py"


rule combine_mutations:
    input:
        mutations = expand("csv/per-sample/{SAMPLE}.mutations.tsv", SAMPLE=ALL_SAMPLES)
    output:
        bed = "csv/mutations/combined.bed",
        mutations = "csv/mutations/combined.tsv"
    run:
        import pandas as pd
        import os

        res = []
        for fh in input.mutations:
            if os.path.getsize(fh) == 1: continue
            df = pd.read_csv(fh, sep="\t").fillna({"dose": 0})
            res.append(df)
        res = pd.concat(res)
        res["ident"] = res.apply(lambda row: f"{row['chrom']}:{row['end']}-{row['end']}:{row['sample']}", axis=1)
        res[["chrom", "start", "end", "sample"]].to_csv(output.bed, sep="\t", index=False, header=False)
        res.to_csv(output.mutations, sep="\t", index=False, header=True)


rule find_multiallelics:
    input:
        mutations = "csv/mutations/combined.tsv",
        bam = lambda wildcards: f"data/cram/downsampled/{wildcards.SAMPLE}.cram" if wildcards.STRATEGY == "downsampled" else f"data/cram/{wildcards.SAMPLE}.cram",
        idx = lambda wildcards: f"data/cram/downsampled/{wildcards.SAMPLE}.cram.crai" if wildcards.STRATEGY == "downsampled" else f"data/cram/{wildcards.SAMPLE}.cram.crai",
        ctrl_bam = lambda wildcards: f"data/cram/{SAMPLE2NORMAL[wildcards.SAMPLE]}.cram",
        ctrl_idx = lambda wildcards: f"data/cram/{SAMPLE2NORMAL[wildcards.SAMPLE]}.cram.crai",
        util = "rules/scripts/utils.py"
    output:
        mutations = "csv/multiallelics/{SAMPLE}.{PERM}.{STRATEGY}.tsv"
    params:
        plot = True,
        samples = lambda wildcards: [s for s in ALL_SAMPLES if SAMPLE2GT[s] != SAMPLE2GT[wildcards.SAMPLE] and s != wildcards.SAMPLE],
    resources:
        runtime = 60
    script:
        "rules/scripts/find_multiallelics.py"


rule test_index_hopping:
    input:
        mutations = "csv/multiallelics.{PERM}.{STRATEGY}.tsv",
        cram = "data/cram/{SAMPLE}.cram",
        idx = "data/cram/{SAMPLE}.cram.crai",
    output:
        mutations = "csv/index_hopping/{SAMPLE}.{PERM}.{STRATEGY}.tsv"
    params:
        threshold = 3
    script:
        "rules/scripts/test_index_hopping.py"
