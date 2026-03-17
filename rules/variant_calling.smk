rule run_deepsomatic:
    input:
        tumor = "data/bam/{SAMPLE}.bam",
        tumor_idx = "data/bam/{SAMPLE}.bam.bai",
        normal = lambda wildcards: f"data/bam/{SAMPLE2NORMAL[wildcards.SAMPLE]}.bam",
        normal_idx = lambda wildcards: f"data/bam/{SAMPLE2NORMAL[wildcards.SAMPLE]}.bam.bai",
        ref = "data/ref/ce11.fa",
        # ref_idx = "data/ref/ce11.fa.gz.fai",
        sif = "clara-parabricks_4.6.0-1.sif"
    output:
        vcf = "data/vcf/{SAMPLE}.vcf.gz",
    resources:
        mem_mb = 64_000,
        runtime = 60,
        gpu = 1,
        cpus_per_gpu = 64,
        slurm_account = "quinlan-gpu-rw",
        slurm_partition = "quinlan-gpu-rw",
        slurm_extra = "--exclusive"
    script:
        "bash_scripts/deepsomatic.sh"