rule plot_mutations:
    input:
        mutations = "csv/mutations/combined.tsv",
        metadata = "data/41467_2020_15912_MOESM4_ESM.xlsx"
    output:
        fig = "fig/mutation_counts.png"
    script:
        "scripts/plot_mutation_counts.py"


rule combine_candidate_multiallelics:
    input:
        mavs = expand("csv/multiallelics/{SAMPLE}.{{PERM}}.{{STRATEGY}}.tsv", SAMPLE=SAMPLES),
        metadata = "data/41467_2020_15912_MOESM4_ESM.xlsx"
    output:
        mavs = "csv/multiallelics.{PERM}.{STRATEGY}.tsv",
        counts = "csv/counts.{PERM}.{STRATEGY}.tsv"
    params:
        threshold = 3
    script:
        "scripts/combine_multiallelics.py"


rule plot_read_support:
    input:
        mavs = "csv/multiallelics.{PERM}.{STRATEGY}.tsv",
        
    output:
        n_reads = "fig/n_reads.{PERM}.{STRATEGY}.png",
        n_reads_hist = "fig/n_reads_hist.{PERM}.{STRATEGY}.png"
    script:
        "scripts/plot_read_support.py"


rule plot_elbow:
    input:
        mavs = "csv/multiallelics.{PERM}.{STRATEGY}.tsv",
    output:
        elbow = "fig/elbow.{PERM}.{STRATEGY}.png",
    script:
        "scripts/plot_elbow.py"


rule plot_mav_agg:
    input:
        # mavs = "csv/multiallelics.{PERM}.{STRATEGY}.tsv",
        counts = "csv/counts.{PERM}.{STRATEGY}.tsv",
        metadata = "data/41467_2020_15912_MOESM4_ESM.xlsx"

    output:
        agg_mavs = "fig/mav_agg.{PERM}.{STRATEGY}.png",
    script:
        "scripts/plot_agg_mavs.py"


rule plot_grouped_mavs:
    input:
        # mavs = "csv/multiallelics.{PERM}.{STRATEGY}.tsv",
        counts = "csv/counts.{PERM}.{STRATEGY}.tsv",
        metadata = "data/41467_2020_15912_MOESM4_ESM.xlsx"
    output:
        grouped_mavs = "fig/mav_grouped.{PERM}.{STRATEGY}.png",
    params:
        group = "Mutagen"
    script:
        "scripts/plot_grouped_mavs.py"


rule plot_entropy:
    input:
        mavs = expand("csv/multiallelics/{SAMPLE}.{{PERM}}.{{STRATEGY}}.tsv", SAMPLE=SAMPLES),
        ref = "data/ref/ce11.fa.gz"
    output:
        entropy = "fig/entropy.{PERM}.{STRATEGY}.png"
    params:
        threshold = 3
    script:
        "scripts/plot_entropy.py"


rule annotate_with_phylop:
    input:
        mutations = "csv/mutations/combined.tsv",
        phylop_scores = "data/ce11.phyloP135way.bedg.gz"
    output:
        phylop = "csv/mutations/combined.phylop.tsv"
    script:
        "scripts/annotate_with_phylop.py"


rule plot_phylop:
    input:
        mutations = "csv/mutations/combined.phylop.tsv"
    output:
        phylop = "fig/phylop.png"
    params:
        mutagens = ["DMS", "EMS", "MMS"]
    script:
        "scripts/plot_phylop.py"


rule plot_mosdepth:
    input:
        summaries = expand("data/mosdepth/{SAMPLE}.mosdepth.summary.txt", SAMPLE=SAMPLES),
        manifest = "data/combined_file_manifest.tsv",
        mutations = "csv/mutations/combined.tsv",
        metadata =  "data/41467_2020_15912_MOESM4_ESM.xlsx"
    output:
        mosdepth = "fig/mosdepth.png"
    script:
        "scripts/plot_mosdepth.py"

