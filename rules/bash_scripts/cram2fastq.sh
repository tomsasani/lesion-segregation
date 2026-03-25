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
        pbrun bam2fq \
            --ref ${snakemake_input[reference]} \
            --in-bam ${snakemake_input[bam]} \
            --out-prefix ${snakemake_params[out_prefix]}