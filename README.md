## Inherited DNA damage in *C. elegans*

authors: Thomas A. Sasani

> This repository contains the necessary code to reproduce figures and analyses in "Inherited DNA damage generates multi-allelic mutations in *C. elegans*," Sasani et al. (2026) *bioRxiv*.

### Prerequisites

#### Clone the repo

```
git clone https://github.com/tomsasani/lesion-segregation
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
    mutations/combined.tsv # filtered SNVs in each strain
    multiallelics.filtered.tsv # MAVs identified in each strain

igv_reports/
    html/ # interactive IGV HTML for every MAV, including reads aligned to the MAV +/- 1kbp in the treated sample and a control

```

### Instructions for reproducing figures

We've included "precomputed" files (filtered SNVs in each strain, multi-allelic variants in each strain, etc.) to make figure generation easier.
It's also possible to run the entire pipeline from scratch (downloading FASTQ, aligning to the reference, etc.), though this will require access to a high-performance computing environment.

#### Mutation file format

Filtered per-strain SNVs are provided in TSV format as follows:

```
chrom	start	end	ref	alt	af	sample
chrIV	13008199	13008200	G	T	0.5	CD0001a
chrV	19417290	19417291	C	T	1.0	CD0001a
```

Per-strain multi-allelic evidence is provided in TSV format as follows:

```
chrom	start	end	ref	af	sample	ident	bam_smp	summary	permuted	is_strand_biased	bias	ma_reads
chrI	14506346	14506347	G	0.53	CD0230d	chrI:14506347-14506347:CD0230d	CD0230d	G:25|A:22	0	False	A:15:7|T:0:0|C:0:0	0
chrIV	5181072	5181073	A	0.51	CD0230d	chrIV:5181073-5181073:CD0230d	CD0230d	T:30|A:26	0	False	T:20:10|C:0:0|G:0:0	0
chrX	515429	515430	T	0.51	CD0230d	chrX:515430-515430:CD0230d	CD0230d	G:31|T:30	0	False	A:0:0|C:0:0|G:20:11	0
```

Most columns are self-explanatory, but descriptions for the others are:

`af`: variant allele frequency of the "called" alternate allele in the VCF from Volkova et al. (2020).

`ident`: unique identifier for the original SNV (`chrom`:`start`:`end`:`sample`)

`bam_smp`: the strain we used to calculate per-base evidence at this SNV. in random trials, this strain will *not* match the strain listed in `ident`.

`summary`: pipe-separated evidence for each nucleotide at this SNV in `bam_smp`'s alignments. 

`permuted`: `0` if we used the focal strain's alignments to count up evidence for each nucleotide. a trial identifier (`1`, `2`, `3`, and so on) if we used a random strain's alignments.

`is_strand_biased`: whether a non-reference allele was strand biased at a binomial p-value threshold of 0.05.

`bias`: for every non-reference allele, the count of reads aligned to the forward or reverse strand.

`ma_reads`: number of reads supporting an 'extra' allele at this SNV in `bam_smp`. if `bam_smp` is the focal sample (i.e., the sample in whom we originally identified the SNV), this is the count of reads supporting a 'third' allele. if `bam_smp` is a random sample, this is the count of reads supporting any non-reference allele.

If you want to generate figures, just run the following:

```
snakemake -j $n_processes \
          -s run.smk
```
