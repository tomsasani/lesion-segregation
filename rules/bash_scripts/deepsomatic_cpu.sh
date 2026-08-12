#!/bin/bash
set -e

module load singularity
export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-labs/quinlan/u1006375/deepsomatic_tmp/

singularity exec --cleanenv \
        -H $SINGULARITYENV_TMPDIR \
        -B /usr/lib/locale/:/usr/lib/locale/ \
            ${snakemake_input[sif]} \
            run_deepsomatic \
                --model_type WGS \
                --num_shards ${snakemake[threads]} \
                --reads_normal ${snakemake_input[normal]} \
                --reads_tumor ${snakemake_input[tumor]} \
                --output_vcf ${snakemake_output[vcf]} \
                --sample_name_tumor ${snakemake_wildcards[SAMPLE]} \
                --sample_name_normal 'control' \
                --ref ${snakemake_input[ref]} \
                --intermediate_results_dir ${snakemake_params[tmpdir]}