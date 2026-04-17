rule index_ref:
    input:
        "data/ref/ce11.fa.gz"
    output:
        "data/ref/ce11.fa.gz.bwt"
    resources:
        runtime = 120
    shell:
        """
        module load bwa
        
        bwa index {input}
        """

rule download_fastq:
    input:
    output: 
        fq1 = temp("data/fastq/{SAMPLE}_1.fastq"),
        fq2 = temp("data/fastq/{SAMPLE}_2.fastq")
    params:
        acc = lambda wildcards: SAMPLE2ACC[wildcards.SAMPLE],
    resources:
        runtime = 240,
    threads: 8
    shell:
        """
        set -euo pipefail
        fasterq-dump -t fasterq_tmp/ -e {threads} -p -O data/fastq/ -o {wildcards.SAMPLE} {params.acc}
        """


rule compress_fastq:
    input: "data/fastq/{SAMPLE}_{FQ}.fastq"
    output: temp("data/fastq/{SAMPLE}_{FQ}.fastq.gz")
    threads: 8
    priority: 100
    resources:
        runtime = 30
    shell:
        """
        bgzip --threads 8 {input}
        """


rule fastq2bam_parabricks:
    input:
        reference = "data/ref/ce11.fa.gz",
        reference_meta = "data/ref/ce11.fa.gz.bwt",
        fq1 = "data/fastq/{SAMPLE}_1.fastq.gz",
        fq2 = "data/fastq/{SAMPLE}_2.fastq.gz",
        sif = "clara-parabricks_4.6.0-1.sif",
    output:
        bam = "data/bam/{SAMPLE}.bam",
        idx = "data/bam/{SAMPLE}.bam.bai"
    priority: 100
    resources:
        mem_mb = 64_000,
        runtime = 30,
        gpu = 1,
        cpus_per_gpu = 64,
        slurm_account = "quinlan-gpu-rw",
        slurm_partition = "quinlan-gpu-rw",
        slurm_extra = "--exclusive"
    script:
        "bash_scripts/fastq2bam.sh"


# rule index_bam:
#     input:
#         bam = "data/bam/{SAMPLE}.bam"
#     output:
#         idx = "data/bam/{SAMPLE}.bam.bai"
#     threads: 8
#     resources:
#         mem_mb = 8_000,
#         runtime = 10,
#     shell:
#         """
#         module load samtools
        
#         samtools index -@ {threads} {input.bam}
#         """
