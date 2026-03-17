rule fastq2bam_parabricks:
    input:
        reference = "data/ref/ce11.fa.gz",
        fq1 = "data/fq/{SAMPLE}_1.fq.gz",
        fq2 = "data/fq/{SAMPLE}_2.fq.gz",
        sif = "clara-parabricks_4.6.0-1.sif",
    output:
        bam = "data/bam/{SAMPLE}.bam"
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


rule index_bam:
    input:
        bam = "data/bam/{SAMPLE}.bam"
    output:
        idx = "data/bam/{SAMPLE}.bam.bai"
    threads: 8
    resources:
        mem_mb = 8_000,
        runtime = 10,
    shell:
        """
        module load samtools
        
        samtools index -@ {threads} {input.bam}
        """
