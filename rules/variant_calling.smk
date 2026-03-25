rule run_deepsomatic:
    input:
        tumor = "data/cram/{SAMPLE}.{FILETYPE}",
        tumor_idx = "data/cram/{SAMPLE}.{FILETYPE}.txt",
        normal = "data/cram/CD0001b.cram",
        normal_idx = "data/cram/CD0001b.cram.txt",
        ref = "data/ref/ce10.renamed.fa",
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