import pandas as pd
import glob
import numpy as np

# get ENA IDs and FTP links
manifest = pd.read_csv("data/combined_file_manifest.tsv", sep="\t")
manifest["filetype"] = manifest["is_bam"].apply(lambda b: "bam" if b else "cram")
manifest["fq_ftp"] = manifest["fastq_ftp"].apply(lambda f: f.split(";")[0].split("_")[0])

SAMPLE2FQPREF = dict(zip(manifest["Sample / DONOR ID"], manifest["fq_ftp"]))
SAMPLE2ACC = dict(zip(manifest["Sample / DONOR ID"], manifest["run_accession"]))
SAMPLE2CRAM = dict(zip(manifest["Sample / DONOR ID"], manifest["aligned_ftp"]))


SAMPLE2FILETYPE = dict(zip(manifest["Sample / DONOR ID"], manifest["filetype"]))

# filter to relevant samples
metadata = pd.read_excel("data/41467_2020_15912_MOESM4_ESM.xlsx", 
                         sheet_name="Sample_description", 
                         skiprows=1).fillna({"Normal": "CD0001b", "Mutagen": "none", "Drug.concentration": 0})

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

raw_counts = pd.read_excel("data/41467_2020_15912_MOESM4_ESM.xlsx", 
                         sheet_name="Mutation counts for all samples", 
                         skiprows=1).rename(columns={"Unnamed: 0": "Sample"})
raw_counts["total"] = np.sum(raw_counts.values[:, 1:97], axis=1)

mutacc = metadata[metadata["Type"].str.contains("acc")]["Sample"].unique()
# get mut. acc samples with relatively high counts
SAMPLES_MA = raw_counts[raw_counts["Sample"].isin(mutacc)].query("total >= 25")["Sample"].to_list()

SAMPLE2GT = dict(zip(metadata["Sample"], metadata["Genotype"]))

configfile: 'config.yaml'
mutagens = config["mutagens"]

# filter to samples with drug concentration > 0
metadata = metadata[metadata["Drug.concentration"] != 0]
metadata = metadata[metadata["Mutagen"].isin(mutagens)]
metadata = metadata[metadata["Generation"] == 1]

include: "rules/fastq2bam.smk"
include: "rules/variant_calling.smk"
include: "rules/test_index_hopping.smk"

CHROMS = ["I", "II", "III", "IV", "V", "X"]
CHROMS = [f"chr{c}" for c in CHROMS]

wildcard_constraints:
    SAMPLE = r"CD\d{4,5}[a-s]",

SAMPLES = metadata["Sample"].to_list()


rule all:
    input:
       
        expand("csv/multiallelics/{SAMPLE}.{PERM}.{STRATEGY}.tsv", SAMPLE=SAMPLES, PERM=[100], STRATEGY=["regular"]),
        

rule index_cram:
    input:
        cram = "data/cram/{SAMPLE}.cram"
    output:
        idx = "data/cram/{SAMPLE}.cram.crai"
    threads: 8
    resources:
        mem_mb = 8_000,
        runtime = 10,
    shell:
        """
        module load samtools/1.16
        
        samtools index -@ {threads} {input.cram}
        """

rule run_mosdepth:
    input:
        cram = "data/cram/{SAMPLE}.cram",
        idx = "data/cram/{SAMPLE}.cram.crai",
        ref = "data/ref/ce11.fa.gz"
    threads: 8
    output: "csv/mosdepth/{SAMPLE}.mosdepth.summary.txt"
    params:
        prefix = lambda wildcards: f"csv/mosdepth/{wildcards.SAMPLE}"
    shell:
        """
        mosdepth -t {threads} \
                 -nx \
                 -f {input.ref} \
                 {params.prefix} \
                 {input.cram}
        """

rule catalog_mutations:
    input:
        vcf = "data/Filtered_VCFS_recovered/SNV/{SAMPLE}_filtered_substitutions.vcf.gz",
        metadata = "data/41467_2020_15912_MOESM4_ESM.xlsx",
        exclude = "data/ce11.sdust.bed.gz",
        ref = "data/ref/ce11.fa.gz"
    output:
        mutations = "csv/per-sample/{SAMPLE}.mutations.tsv"
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
        res["ident"] = res.apply(lambda row: f"chr{row['chrom']}:{row['end']}-{row['end']}:{row['sample']}", axis=1)
        res["true_chrom"] = res["chrom"].apply(lambda c: f"chr{c}")
        res[["true_chrom", "start", "end", "sample"]].to_csv(output.bed, sep="\t", index=False, header=False)
        res.to_csv(output.mutations, sep="\t", index=False, header=True)


rule plot_mutations:
    input:
        mutations = "csv/mutations/combined.tsv",
        metadata = "data/41467_2020_15912_MOESM4_ESM.xlsx"
    output:
        fig = "fig/mutation_counts.png"
    script:
        "rules/scripts/plot_mutation_counts.py"


def get_downsampling_frac(wildcards):
    df = pd.read_csv(f"data/mosdepth/{wildcards.SAMPLE}.mosdepth.summary.txt", sep="\t")
    df = df[df["chrom"] == "total"]
    dp = df["mean"].values[0]

    adj = 30 / dp
    if adj > 1:
        adj = dp
    return adj


rule downsample:
    input:
        cram = "data/cram/{SAMPLE}.cram",
        idx = "data/cram/{SAMPLE}.cram.crai",
        ref = "data/ref/ce11.fa.gz",
        depth = "data/mosdepth/{SAMPLE}.mosdepth.summary.txt"
    output:
        cram = "data/cram/downsampled/{SAMPLE}.cram",
        idx = "data/cram/downsampled/{SAMPLE}.cram.crai"
    params:
        downsampling_frac = lambda wildcards: get_downsampling_frac(wildcards)
    threads: 8
    shell:
        """
        module load samtools/1.16
        
        samtools view -O CRAM \
                      --write-index \
                      -C \
                      -s {params.downsampling_frac} \
                      -T {input.ref} \
                      -@ {threads} \
                      -o {output.cram} \
                      {input.cram}
        """


rule find_multiallelics:
    input:
        mutations = "csv/mutations/combined.tsv",
        bam = lambda wildcards: f"data/cram/downsampled/{wildcards.SAMPLE}.cram" if wildcards.STRATEGY == "downsample" else f"data/cram/{wildcards.SAMPLE}.cram",
        idx = lambda wildcards: f"data/cram/downsampled/{wildcards.SAMPLE}.cram.crai" if wildcards.STRATEGY == "downsample" else f"data/cram/{wildcards.SAMPLE}.cram.crai",
        ctrl_bam = lambda wildcards: f"data/cram/{SAMPLE2NORMAL[wildcards.SAMPLE]}.cram",
        ctrl_idx = lambda wildcards: f"data/cram/{SAMPLE2NORMAL[wildcards.SAMPLE]}.cram.crai"
    output:
        mutations = "csv/multiallelics/{SAMPLE}.{PERM}.{STRATEGY}.tsv"
    params:
        plot = False,
        samples = lambda wildcards: [s for s in ALL_SAMPLES if SAMPLE2GT[s] != SAMPLE2GT[wildcards.SAMPLE]],
    resources:
        runtime = 20
    script:
        "rules/scripts/find_multiallelics.py"


rule combine_multiallelics:
    input:
        mutations = expand("csv/multiallelics/{SAMPLE}.{{PERM}}.{{STRATEGY}}.tsv", SAMPLE=SAMPLES)
    output: 
        mutations = "csv/multiallelics.{PERM}.{STRATEGY}.tsv"
    params:
        threshold = 3
    run:
        import pandas as pd

        res = []
        for fh in input.mutations:
            try:
                df = pd.read_csv(fh, sep="\t")
            except (pd.errors.EmptyDataError, FileNotFoundError):
                continue
            df = df[df["ma_reads"] > 1]
            df = df[df["permuted"] == 0]
            df["is_multiallelic"] = df["ma_reads"] >= params.threshold
            res.append(df)

        res = pd.concat(res)
        res.to_csv(output.mutations, sep="\t", index=False)


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
   