# rule run_deepsomatic:
#     input:
#         tumor = lambda wildcards: f"data/bam/{SAMPLE2NORMAL[wildcards.SAMPLE]}.bam",
#         normal = "data/bam/{SAMPLE}.bam",
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
        tumor = lambda wildcards: f"data/bam/{SAMPLE2NORMAL[wildcards.SAMPLE]}.bam",
        normal = "data/bam/{SAMPLE}.bam",
        ref = "data/ref/ce11.fa.gz",
        sif = "deepsomatic_1.10.0.sif"
    output:
        vcf = "data/vcf/{SAMPLE}.vcf.gz",
    threads: 24
    params:
        tmpdir = "deepsomatic_tmp/"
    resources:
        mem_mb = 64_000,
        runtime = 240,
    script:
        "bash_scripts/deepsomatic_cpu.sh"


rule filter_vcf:
    input:
        "data/vcf/{SAMPLE}.vcf.gz"
    output:
        "data/vcf/prefiltered/{SAMPLE}.vcf.gz"
    shell:
        """
        module load bcftools
        
        bcftools view -f 'PASS' -v snps,mnps -Oz -o {output} {input}
        """

rule remove_exclude:
    input:
        vcf = "data/vcf/prefiltered/{SAMPLE}.vcf.gz",
        exclude = "data/ce10.rmsk.renamed.bed.gz"
    output:
        vcf = "data/vcf/filtered/{SAMPLE}.vcf.gz"
    shell:
        """
        module load bedtools
        
        bedtools intersect -a {input.vcf} -b {input.exclude} -v -header | bgzip > {output}
        """


rule index_vcf:
    input:
        "data/vcf/filtered/{SAMPLE}.vcf.gz"
    output:
        "data/vcf/filtered/{SAMPLE}.vcf.gz.tbi"
    shell:
        """
        module load bcftools
        
        bcftools index --tbi {input}
        """



