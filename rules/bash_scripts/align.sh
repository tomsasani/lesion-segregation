#!/bin/bash
set -euo pipefail

module load bwa gatk/4.6 samtools/1.16

export SAMTOOLS_TMPDIR=/scratch/ucgd/lustre-labs/quinlan/u1006375/samtools_tmp

bwa mem \
    -t ${snakemake[threads]} \
    -K ${snakemake_params[K]} \
    -R ${snakemake_params[read_group]} \
    ${snakemake_input[reference]} \
    ${snakemake_input[fq1]} \
    ${snakemake_input[fq2]} \
    | \
  gatk SortSam \
    --java-options -Xmx24g \
    --MAX_RECORDS_IN_RAM ${snakemake_params[max_records_in_ram]} \
    -I /dev/stdin \
    -O ${snakemake_output[bam]} \
    --SORT_ORDER coordinate \
    --TMP_DIR $SAMTOOLS_TMPDIR \
    --CREATE_INDEX true


