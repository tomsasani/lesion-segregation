import pandas as pd
import glob

include: "rules/fastq2bam.smk"
include: "rules/variant_calling.smk"


def check_max_conc(row: pd.Series, treatment2conc):
    drug_conc = row["Drug.concentration"]
    mutagen = row["Mutagen"]
    genotype = row["Genotype"]
    treatment = (genotype, mutagen)
    if treatment in treatment2conc:
        if drug_conc == treatment2conc[treatment]:
            return 1
        else:
            return 0
    else:
        return 0

# get ENA IDs and FTP links
manifest = pd.read_csv("data/combined_file_manifest.tsv", sep="\t")

SAMPLE2CRAM = dict(zip(manifest["Sample / DONOR ID"], manifest["aligned_ftp"]))
SAMPLE2FILETYPE = dict(zip(manifest["Sample / DONOR ID"], manifest["is_bam"]))


# filter to relevant samples
metadata = pd.read_excel("data/41467_2020_15912_MOESM4_ESM.xlsx", 
                         sheet_name="Sample_description", 
                         skiprows=1).fillna({"Normal": "CD0001b", "Drug.concentration": 0})

metadata = metadata[metadata["Drug.concentration"] != 0]

# map combinations of genotypes and mutagens to concentrations and get max per treatment
treatment2conc = metadata.groupby(["Genotype", "Mutagen"]).agg(conc_max=("Drug.concentration", "max")).to_dict()["conc_max"]
metadata["is_max_conc"] = metadata.apply(lambda row: check_max_conc(row, treatment2conc), axis=1)

configfile: 'config.yaml'
mutagens: list[str] = config["mutagens"]
genotypes: list[str] = config["genotypes"]
genotypes = metadata["Genotype"].unique()

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

CHROMS = ["I", "II", "III", "IV", "V", "X"]
CHROMS = [f"chr{c}" for c in CHROMS]

wildcard_constraints:
    SAMPLE = r"CD\d{4,5}[a-j]",
    FILETYPE = f"bam|cram"

SAMPLES = metadata["Sample"].to_list()

SAMPLES = [s for s in SAMPLES if s not in ("CD0013a", "CD0117a", "CD0014d", "CD0806b", "CD0110d")]

FILETYPES = ["bam" if SAMPLE2FILETYPE[s] else "cram" for s in SAMPLES]

print (metadata[metadata["Sample"].isin(SAMPLES)].groupby("Genotype").size())
rule all:
    input:
        # expand("csv/multiallelics/{SAMPLE}.{PERMUTATION}.tsv", SAMPLE=SAMPLES, PERMUTATION=range(10)),
        # expand("csv/{SAMPLE}.mutations.tsv", SAMPLE=SAMPLES),
        expand("data/cram/{SAMPLE}.{FILETYPE}", zip, SAMPLE=SAMPLES, FILETYPE=FILETYPES),
        # "csv/mutations/combined.ce10.tsv"
        # expand("fig/cnvkit/{SAMPLE}.{CHROM}.pdf", SAMPLE=SAMPLES, CHROM=CHROMS)


rule download_cram:
    input:
    output: "data/cram/{SAMPLE}.{FILETYPE}"
    params:
        ftp = lambda wildcards: SAMPLE2CRAM[wildcards.SAMPLE],
    resources:
        runtime = 240
    shell:
        """
        curl -o {output} {params.ftp}
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
    shell:
        """
        module load samtools
        
        samtools index -@ {threads} {input.bam}
        touch {output}
        """



rule catalog_mutations:
    input:
        vcf = "data/Filtered_VCFS_recovered/SNV/{SAMPLE}_filtered_substitutions.vcf.gz",
        # vcf = "data/vcf/{SAMPLE}.vcf.gz",
        metadata = "data/41467_2020_15912_MOESM4_ESM.xlsx",
        ref = "data/ref/ce11.fa",
        exclude = "data/repeats.tsv",
    output:
        mutations = "csv/per-sample/{SAMPLE}.mutations.tsv"
    resources:
        runtime = 5
    script:
        "rules/scripts/catalog_mutations.py"


rule combine_mutations:
    input:
        mutations = expand("csv/per-sample/{SAMPLE}.mutations.tsv", SAMPLE=SAMPLES)
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
        res["chrom"] = res["chrom"].apply(lambda c: f"chr{c}")
        res[["chrom", "start", "end", "sample"]].to_csv(output.bed, sep="\t", index=False, header=False)
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
        old["ident"] = old.apply(lambda row: f"{row['chrom']}:{row['end']}-{row['end']}:{row['sample']}", axis=1)
        old = old.drop(columns=["chrom", "start", "end"])
        new = pd.read_csv(input.mutations, sep="\t", names=["chrom", "start", "end", "ident"])
        new = new.merge(old, on="ident")

        new["chrom"] = new["chrom"].apply(lambda c: c.replace("chr", "CHROMOSOME_"))
        new.to_csv(output.mutations, sep="\t", index=False)



rule find_multiallelics:
    input:
        mutations = "csv/mutations/combined.ce10.tsv",
        bam = lambda wildcards: f"data/cram/{wildcards.SAMPLE}.{'bam' if SAMPLE2FILETYPE[wildcards.SAMPLE] else 'cram'}",
        meta = lambda wildcards: f"data/cram/{wildcards.SAMPLE}.{'bam.txt' if SAMPLE2FILETYPE[wildcards.SAMPLE] else 'cram.txt'}"
    output:
        mutations = "csv/multiallelics/{SAMPLE}.{PERMUTATION}.tsv"
    params:
        thresh = 2,
        plot = True
    resources:
        runtime = 5
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