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
    priority: 80
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
    priority: 90
    resources:
        runtime = 30,

    shell:
        """
        bgzip --threads 8 {input}
        """


rule clean_with_fastp:
    input:
        fastp_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/fastp",
        fq1 = "data/fastq/{SAMPLE}_1.fastq.gz",
        fq2 = "data/fastq/{SAMPLE}_2.fastq.gz",
    output:
        fq1_clean = temp("data/fastq/{SAMPLE}.1.clean.fastq.gz"),
        fq2_clean = temp("data/fastq/{SAMPLE}.2.clean.fastq.gz"),
        fq1_unpaired = temp("data/fastq/{SAMPLE}.1.unpaired.fastq.gz"),
        fq2_unpaired = temp("data/fastq/{SAMPLE}.2.unpaired.fastq.gz"),
        html_report = "reports/{SAMPLE}-fastp-report.html",
        json_report = "reports/{SAMPLE}-fastp-report.json"
    threads: 4
    resources:
        runtime = 30
    log: "logs/fastp/{SAMPLE}.log"
    shell:
        """
        {input.fastp_binary} --in1 {input.fq1} \
                --in2 {input.fq2} \
                --out1 {output.fq1_clean} \
                --out2 {output.fq2_clean} \
                --unpaired1 {output.fq1_unpaired} \
                --unpaired2 {output.fq2_unpaired} \
                --thread {threads} \
                --html {output.html_report} \
                --json {output.json_report} \
                -V
        """


rule fastq2bam:
    input:
        reference = "data/ref/ce11.fa.gz",
        fq1 = "data/fastq/{SAMPLE}.1.clean.fastq.gz",
        fq2 = "data/fastq/{SAMPLE}.2.clean.fastq.gz",
        fq1_unpaired = "data/fastq/{SAMPLE}.1.unpaired.fastq.gz",
        fq2_unpaired = "data/fastq/{SAMPLE}.2.unpaired.fastq.gz",
        bwa_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/src/bwa-mem2-2.3_x64-linux/bwa-mem2",
    output:
        bam = temp("data/bam/{SAMPLE}.bam"),
    priority: 100
    threads: 16
    resources:
        mem_mb = 32_000
    params:
        K = 96_000_000,
        max_records_in_ram = 5_000_000,
        read_group = lambda wildcards: f"@RG\\tID:{wildcards.SAMPLE}\\tLB:lib1\\tPL:ILLUMINA\\tSM:{wildcards.SAMPLE}\\tPU:{wildcards.SAMPLE}"
    script:
        "bash_scripts/align.sh"


rule mark_dups:
    input:
        bam = "data/bam/{SAMPLE}.bam"
    output:
        bam = temp("data/bam/{SAMPLE}.dupmarked.bam"),
        metrics = "data/dupmetrics/{SAMPLE}.txt"
    resources:
        mem_mb = 32_000
    shell:
        """
        module load gatk/4.6

        gatk MarkDuplicates \
            --java-options -Xmx24g \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics}
        """


rule convert_to_cram:
    input:
       bam = "data/bam/{SAMPLE}.dupmarked.bam",
       reference = "data/ref/ce11.fa.gz",
    output:
        cram = "data/cram/{SAMPLE}.cram"
    threads: 8
    shell:
        """
        module load samtools/1.16

        samtools view -@ {threads} \
                      -C \
                      -o {output.cram} \
                      -T {input.reference} \
                      --output-fmt-option embed_ref=1 \
                      {input.bam}
        """