import pandas as pd

wildcard_constraints:
    SAMPLE = r"CD\d{4,5}[a-s]",

mavs = pd.read_csv("csv/multiallelics.10.original.tsv", sep="\t")

# samples with MAVs
SAMPLES = list(mavs[(mavs["permuted"] == 0) & (mavs["is_multiallelic"])]["sample"].unique())

# filter to relevant samples
metadata = pd.read_excel("data/41467_2020_15912_MOESM4_ESM.xlsx", 
                         sheet_name="Sample_description", 
                         skiprows=1).fillna({"Normal": "CD0001b", "Mutagen": "none", "Drug.concentration": 0})

metadata = metadata[~metadata["Sample"].str.contains("d2")]

# map samples to their corresponding normal sample
SAMPLE2NORMAL = {}
for k,v in zip(metadata["Sample"], metadata["Normal"]):
    if "," in v:
        SAMPLE2NORMAL[k] = v.split(",")[0]
    elif " " in v:
        SAMPLE2NORMAL[k] = v.split(" ")[0]
    else:
        SAMPLE2NORMAL[k] = v

rule all:
    input:
        expand("igv_reports/html/{SAMPLE}.html", SAMPLE=SAMPLES)


rule create_mav_bed:
    input:
        mavs = "csv/multiallelics.10.original.tsv"
    output:
        bed = "igv_reports/bed/{SAMPLE}.multiallelics.bed",
        ann = "igv_reports/bed/{SAMPLE}.status.bed"
    params:
        threshold = 3
    run:
        import pandas as pd
        
        df = pd.read_csv(input.mavs, sep="\t", dtype={"sample_id": str})
        df = df[df["ma_reads"] >= params.threshold]
        df = df[df["sample"] == wildcards.SAMPLE]
        df = df[df["permuted"] == 0]
        df["is_mav"] = True
        df["adj_start"] = df["start"] - 1_000
        df["adj_end"] = df["end"] + 1_000
        df[["chrom", "adj_start", "adj_end"]].to_csv(output.bed, sep="\t", index=False, header=False)
        df[["chrom", "start", "end", "is_mav"]].to_csv(output.ann, sep="\t", index=False, header=False)


rule make_report:
    input:
        bed = "igv_reports/bed/{SAMPLE}.multiallelics.bed",
        ann = "igv_reports/bed/{SAMPLE}.status.bed",
        cram = "data/cram/{SAMPLE}.cram",
        ctrl_cram = lambda wildcards: f"data/cram/{SAMPLE2NORMAL[wildcards.SAMPLE]}.cram",
        ref = "data/ref/ce11.fa.gz"
    output:
        html = "igv_reports/html/{SAMPLE}.html"
    shell:
        """
        create_report {input.bed} \
                      --fasta {input.ref} \
                      --tracks {input.ctrl_cram} {input.cram} {input.ann} \
                      --output {output.html}
        """