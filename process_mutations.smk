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
# SAMPLE2FILETYPE["CD0048a"] = "bam"
# SAMPLE2FILETYPE["CD0001b"] = "bam"

# filter to relevant samples
metadata = pd.read_excel("data/41467_2020_15912_MOESM4_ESM.xlsx", 
                         sheet_name="Sample_description", 
                         skiprows=1).fillna({"Normal": "CD0001b", "Mutagen": "none", "Drug.concentration": 0})

# only consider samples with at least 2 replicates per treatment
# metadata["treatment"] = metadata.apply(lambda row: f'{row["Genotype"]}_{row["Mutagen"]}', axis=1)
# sample_per_treatment = metadata.groupby(["treatment"]).size().reset_index().rename(columns={0: "count"})
# multi_treatments = sample_per_treatment.query("count > 1")["treatment"].to_list()

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


configfile: 'config.yaml'
mutagens = config["mutagens"]
genotypes = metadata["Genotype"].unique()

# filter to samples with drug concentration > 0
metadata = metadata[metadata["Drug.concentration"] != 0]
metadata = metadata[metadata["Mutagen"].isin(mutagens)]
metadata = metadata[metadata["Genotype"].isin(genotypes)]

# get genotoxin samples
SAMPLES_GT = metadata["Sample"].to_list()

include: "rules/fastq2bam.smk"
include: "rules/variant_calling.smk"

CHROMS = ["I", "II", "III", "IV", "V", "X"]
CHROMS = [f"chr{c}" for c in CHROMS]

wildcard_constraints:
    SAMPLE = r"CD\d{4,5}[a-s]",
    # FILETYPE = f"bam|cram"

SAMPLES = SAMPLES_GT # + SAMPLES_MA
print (len(SAMPLES))
# SAMPLES = [s for s,t in sample2treat.items() if t in multi_treatments]
# manifest = pd.read_csv("data/combined_file_manifest.tsv", sep="\t")
# smp2tech = dict(zip(manifest["Sample / DONOR ID"], manifest["experiment_title"]))

# FILETYPES = [SAMPLE2FILETYPE[s] for s in SAMPLES]

# existing = [fh.split("/")[-1].split(".")[0] for fh in glob.glob("data/cram/*.txt")]
# SAMPLES = [s for s in SAMPLES if s in existing]

rule all:
    input:
        # "csv/mutations/combined.ce10.tsv"
        # expand("csv/multiallelics/{SAMPLE}.{PERM}.tsv", SAMPLE=SAMPLES_GT, PERM=[11]),
        expand("data/bam/{SAMPLE}.bam.bai", SAMPLE=SAMPLES[:500])
        

rule download_cram:
    input:
    output: "data/cram/{SAMPLE}.{FILETYPE}"
    params:
        ftp = lambda wildcards: SAMPLE2CRAM[wildcards.SAMPLE],
    resources:
        runtime = 240,
        slurm_extra = "--exclude=rw160,rw159,rw157,rw158"
    shell:
        """
        wget -O {output} http://{params.ftp}
        sleep 20
        """


rule index:
    input:
        bam = "data/cram/{SAMPLE}.{FILETYPE}"
    output:
        meta = "data/cram/{SAMPLE}.{FILETYPE}.txt"
    threads: 8
    resources:
        mem_mb = 8_000,
        runtime = 20,
        slurm_extra = "--exclude=rw160,rw159"

    shell:
        """
        module load samtools
        
        samtools index -@ {threads} {input.bam}
        touch {output}
        """


rule get_doubletons:
    input:
        vcfs = expand("data/vcf/filtered/{SAMPLE}.vcf.gz", SAMPLE=SAMPLES),
        idxs = expand("data/vcf/filtered/{SAMPLE}.vcf.gz.tbi", SAMPLE=SAMPLES)
    output:
        "data/vcf/merged/to_exclude.vcf.gz"
    threads: 8
    params: threshold = len(SAMPLES) - 2
    shell:
        """
        module load bcftools
        
        bcftools merge --threads {threads} {input.vcfs} | bcftools view -i 'N_MISSING <= {params.threshold}' -Oz -o {output}
        """


rule remove_doubletons:
    input:
        exclude = "data/vcf/merged/to_exclude.vcf.gz",
        vcf = "data/vcf/filtered/{SAMPLE}.vcf.gz"
    output: "data/vcf/filtered/{SAMPLE}.filtered.vcf.gz"
    shell:
        """
        module load bedtools

        bedtools intersect -v -header -a {input.vcf} -b {input.exclude} | bgzip > {output}
        """


rule catalog_mutations:
    input:
        vcf = "data/Filtered_VCFS_recovered/SNV/{SAMPLE}_filtered_substitutions.vcf.gz",
        metadata = "data/41467_2020_15912_MOESM4_ESM.xlsx",
        exclude = "data/ce11.rmsk.bed.gz"
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
            df = pd.read_csv(fh, sep="\t")
            res.append(df)
        res = pd.concat(res)
        res["ident"] = res.apply(lambda row: f"chr{row['chrom']}:{row['end']}-{row['end']}:{row['sample']}", axis=1)
        res["true_chrom"] = res["chrom"].apply(lambda c: f"chr{c}")
        res[["true_chrom", "start", "end", "sample"]].to_csv(output.bed, sep="\t", index=False, header=False)
        res.to_csv(output.mutations, sep="\t", index=False, header=True)


rule liftover_mutations:
    input:
        mutations = "csv/mutations/combined.bed",
        chain = "data/ce11ToCe10.over.chain.gz"
    output:
        mutations = "csv/mutations/liftover.ce10.bed",
        unmapped = "csv/mutations/unmapped.ce10.bed"
    shell:
        """
        liftOver -preserveInput {input.mutations} {input.chain} {output.mutations} {output.unmapped}
        """


rule merge_mutations:
    input:
        mutations = "csv/mutations/liftover.ce10.bed",
        old_mutations = "csv/mutations/combined.tsv"
    output:
        mutations = "csv/mutations/combined.ce10.tsv"
    run:
        import pandas as pd
        old = pd.read_csv(input.old_mutations, sep="\t")
        # old["ident"] = old.apply(lambda row: f"chr{row['chrom']}:{row['end']}-{row['end']}:{row['sample']}", axis=1)
        old = old.drop(columns=["chrom", "start", "end"])
        new = pd.read_csv(input.mutations, sep="\t", names=["chrom", "start", "end", "ident"])

        print (new)
        print (old)
        new = new.merge(old, on="ident")
        print (new)
        new["chrom"] = new["chrom"].apply(lambda c: c.replace("chr", "CHROMOSOME_"))
        new.to_csv(output.mutations, sep="\t", index=False)



rule find_multiallelics:
    input:
        mutations = "csv/mutations/combined.ce10.tsv",
        # mutations = "csv/mutations/combined.tsv",
        # bam = "data/bam/{SAMPLE}.bam",
        # idx = "data/bam/{SAMPLE}.bam.bai"
        bam = lambda wildcards: f"data/cram/{wildcards.SAMPLE}.{SAMPLE2FILETYPE[wildcards.SAMPLE]}",
        meta = lambda wildcards: f"data/cram/{wildcards.SAMPLE}.{SAMPLE2FILETYPE[wildcards.SAMPLE]}.txt"
    output:
        mutations = "csv/multiallelics/{SAMPLE}.{PERM}.tsv"
    params:
        plot = False,
        shared_treatment = lambda wildcards: [k for k,v in sample2treat.items() if v == sample2treat[wildcards.SAMPLE]]
    resources:
        runtime = 20
    script:
        "rules/scripts/find_multiallelics.py"


rule run_cnvkit:
    input:
        bam = "data/bam/{SAMPLE}.bam",
        ctrl_bam = lambda wildcards: "data/bam/" + SAMPLE2NORMAL[wildcards.SAMPLE] + ".bam",
        # reference = "cnvs/CD0001b.targetcoverage.cnn"
    output: "cnvs/{SAMPLE}/{SAMPLE}.cnr"
    threads: 4
    resources:
        runtime = 60
    shell:
        """
        cnvkit.py batch {input.bam} \
                        -m wgs \
                        -f data/ref/ce11.fa \
                        --segment-method none \
                        -d cnvs/{wildcards.SAMPLE}/ \
                        -n {input.ctrl_bam} \
                        -p {threads}
        """

rule plot_cnvkit:
    input: "cnvs/{SAMPLE}/{SAMPLE}.cnr"
    output: "fig/cnvkit/{SAMPLE}.{CHROM}.pdf"
    shell:
        """
        cnvkit.py scatter -c {wildcards.CHROM} \
                          -o {output} \
                          {input}
        """