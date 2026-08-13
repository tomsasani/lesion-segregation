## Transgenerational effects of DNA damage in *C. elegans*

authors: Thomas A. Sasani

> This repository contains the necessary code to reproduce figures and analyses in "Transgenerational effects of DNA damage in *C. elegans*," Sasani et al. (2026) *bioRxiv*.

### Prerequisites

#### Clone the repo

```
git clone https://github.com/quinlan-lab/lesion-segregation
```

#### `python` dependencies

Using your favorite `python` environment manager (e.g., `pyenv` or `conda`), install the necessary dependencies in `requirements.txt`.

#### Metadata and variants from Volkova *et al.* (2020)

Our analysis uses data from a [prior study](https://www.nature.com/articles/s41467-020-15912-7) of mutagenesis in *C. elegans*.

1. Download filtered SNV VCF files for each strain in [Supplementary Data 6](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-020-15912-7/MediaObjects/41467_2020_15912_MOESM9_ESM.zip).
Move the folder of variants (`Filtered_VCFS_recovered`) to a top-level directory called `data`.
2. Download strain metadata in [Supplementary Data 1](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-020-15912-7/MediaObjects/41467_2020_15912_MOESM4_ESM.xlsx).
Move the Excel file (`41467_2020_15912_MOESM4_ESM.xlsx`) to a top-level directory called `data`.




### Directory structure

```
run_pipeline.smk # main snakemake pipeline

rules/
    fastq2bam.smk # download FASTQ, align to ce11, convert to CRAM
    variant_calling.smk # use DeepSomatic to call variants (not used for publication analyses)
    generate_figures.smk # use precomputed mutation data to make figures
    downsampling.smk # randomly downsample CRAMs to determine effects of sequencing depth on MAV detection

    scripts/ # python scripts for analysis and figure generation
    bash_scripts/ # bash scripts for alignment/variant calling/etc.

csv/
    per-sample/ # filtered SNVs in each strain
    multiallelics/ # MAVs identified in each strain
```

### Instructions for reproducing figures

We've included "precomputed" files (filtered SNVs in each strain, multi-allelic variants in each strain, etc.) to make figure generation easier.
It's also possible to run the entire pipeline from scratch (downloading FASTQ, aligning to the reference, etc.), though this will require access to a high-performance computing environment.

```
snakemake -j $n_processes \
          -s run.smk
```
