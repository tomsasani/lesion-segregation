# rule run_deepsomatic:
#     input:
#         normal = lambda wildcards: f"data/cram/{SAMPLE2NORMAL[wildcards.SAMPLE]}.cram",
#         tumor = "data/cram/{SAMPLE}.cram",
#         ref = "data/ref/ce11.fa.gz",
#         sif = "clara-parabricks_4.6.0-1.sif"
#     output:
#         vcf = "data/vcf/{SAMPLE}.vcf.gz",
#     threads: 24
#     resources:
#         mem_mb = 64_000,
#         runtime = 30,
#         gpu = 1,
#         cpus_per_gpu = 64,
#         slurm_account = "quinlan-gpu-rw",
#         slurm_partition = "quinlan-gpu-rw",
#         slurm_extra = "--exclusive"
#     script:
#         "bash_scripts/deepsomatic.sh"


rule run_deepsomatic_cpu:
    input:
        normal = lambda wildcards: f"data/cram/{SAMPLE2NORMAL[wildcards.SAMPLE]}.cram",
        tumor = "data/cram/{SAMPLE}.cram",
        ref = "data/ref/ce11.fa.gz",
        sif = "deepsomatic_1.10.0.sif"
    output:
        vcf = "data/vcf/{SAMPLE}.vcf.gz",
    threads: 16
    params:
        tmpdir = lambda wildcards: f"deepsomatic_tmp/{wildcards.SAMPLE}"
    resources:
        mem_mb = 64_000,
        runtime = 480,
    script:
        "bash_scripts/deepsomatic_cpu.sh"


rule filter_vcf:
    input: "data/vcf/{SAMPLE}.vcf.gz"
    output: "data/vcf/filtered/{SAMPLE}.vcf.gz"
    shell:
        """
        module load bcftools
        
        bcftools view -f 'PASS' -v snps,mnps -Oz -o {output} {input}
        """


rule index_vcf:
    input: "data/vcf/filtered/{SAMPLE}.vcf.gz"
    output: "data/vcf/filtered/{SAMPLE}.vcf.gz.tbi"
    shell:
        """
        module load bcftools
        
        bcftools index --tbi {input}
        """



