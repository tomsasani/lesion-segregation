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

def get_downsampling_frac(wildcards):
    df = pd.read_csv(f"data/mosdepth/{wildcards.SAMPLE}.mosdepth.summary.txt", sep="\t")
    df = df[df["chrom"] == "total"]
    dp = df["mean"].values[0]

    adj = 30 / dp
    if adj > 1:
        adj = 0.999
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


rule run_mosdepth_downsampled:
    input:
        cram = "data/cram/downsampled/{SAMPLE}.cram",
        idx = "data/cram/downsampled/{SAMPLE}.cram.crai",
        ref = "data/ref/ce11.fa.gz"
    threads: 8
    output: "csv/mosdepth/downsampled/{SAMPLE}.mosdepth.summary.txt"
    params:
        prefix = lambda wildcards: f"csv/mosdepth/downsampled/{wildcards.SAMPLE}"
    shell:
        """
        mosdepth -t {threads} \
                 -nx \
                 -f {input.ref} \
                 {params.prefix} \
                 {input.cram}
        """

