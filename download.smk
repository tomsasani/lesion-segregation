import pandas as pd

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
mapping = pd.read_csv("data/CDnr_to_ENA_ms_Volkova_2020.csv")
file_report_p1 = pd.read_csv("data/filereport_read_run_PRJEB2730_tsv.txt", sep="\t")
file_report_p2 = pd.read_csv("data/filereport_read_run_PRJEB4769_tsv.txt", sep="\t")
file_report = pd.concat([file_report_p1, file_report_p2])

merged = mapping.merge(file_report, left_on="Accession number", right_on="secondary_sample_accession")

merged["fq1"] = merged["fastq_ftp"].apply(lambda f: f.split(";")[0])

# map sample names to FQ names
SAMPLE_NAMES = merged["Sample / DONOR ID"].values
SAMPLE_ACCS = merged["run_accession"].values
FQ_PREF = [fh.split("_")[0] for fh in merged["fq1"].values]

sample2fq = dict(zip(SAMPLE_NAMES, FQ_PREF))
sample2acc = dict(zip(SAMPLE_NAMES, SAMPLE_ACCS))


# filter to relevant samples
metadata = pd.read_excel("data/41467_2020_15912_MOESM4_ESM.xlsx", sheet_name="Sample_description", skiprows=1).fillna({"Normal": "CD0001b", "Drug.concentration": 0})

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

# map combinations of genotypes and mutagens to concentrations and get max per treatment
treatment2conc = metadata.groupby(["Genotype", "Mutagen"]).agg(conc_max=("Drug.concentration", "max")).to_dict()["conc_max"]
metadata["is_max_conc"] = metadata.apply(lambda row: check_max_conc(row, treatment2conc), axis=1)
metadata = metadata[metadata["is_max_conc"] == 1]

metadata = metadata[metadata["Drug.concentration"] != 0]


SAMPLES = metadata["Sample"].to_list()

CHROMS = ["I", "II", "III", "IV", "V", "X"]
CHROMS = [f"chr{c}" for c in CHROMS]

wildcard_constraints:
    SAMPLE = r"[A-Z]{2}[0-9]{4}[a-z]{1}",
    CHROM = r"chr(I|II|III|IV|V|X)"

include: "rules/fastq2bam.smk"
include: "rules/variant_calling.smk"

rule all:
    input:
        vcfs = expand("data/vcf/{SAMPLE}.vcf.gz", SAMPLE=SAMPLES)
        # bams = expand("data/bam/{SAMPLE}.bam", SAMPLE=["CD0496d"])
        # vcfs = expand("data/vcf/per_chrom/{CHROM}.{SAMPLE}.filtered.vcf.gz", CHROM=CHROMS, SAMPLE=["CD0495a"])

rule download_fq:
    input:
    output: temp("data/fq/{SAMPLE}_{FQ_NUM}.fq.gz")
    params:
        accession = lambda wildcards: sample2acc[wildcards.SAMPLE],
        fq_ftp = lambda wildcards: sample2fq[wildcards.SAMPLE] + "_" + wildcards.FQ_NUM + ".fastq.gz"
    resources:
        runtime = 60
    shell:
        """
        wget -O {output} ftp://{params.fq_ftp}
        """

rule combine_bcf:
    input:
        vcf = "data/vcf/{SAMPLE}.vcf.gz",
        chrom_conv = "data/chromosome_conversion.txt"
    output: "data/vcf/{SAMPLE}.filtered.vcf.gz"
    resources:
        runtime = 60
    shell:
        """
        module load bcftools
        
        bcftools annotate --rename-chrs {input.chrom_conv} {input.vcf} | bgzip > {output}
        """

