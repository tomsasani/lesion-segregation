#!/bin/bash
set -e

pwd
nvidia-smi

module load singularity

export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-labs/quinlan/parabricks_tmp/
export SAMTOOLS_TMPDIR=/scratch/ucgd/lustre-labs/quinlan/u1006375/samtools_tmp

singularity exec --cleanenv \
        --nv \
        -H $SINGULARITYENV_TMPDIR \
        --bind $(pwd):/workdir \
        --bind $(pwd):/outputdir \
        --pwd /workdir \
        ${snakemake_input[sif]} \
        pbrun deepsomatic \
            --ref ${snakemake_input[ref]} \
            --in-tumor-bam ${snakemake_input[tumor]} \
            --in-normal-bam ${snakemake_input[normal]} \
            --out-variants ${snakemake_output[vcf]} \
            --tmp-dir $SAMTOOLS_TMPDIR \
            --num-cpu-threads-per-stream 2
