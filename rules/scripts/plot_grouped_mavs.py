import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import argparse


from snakemake_argparse_bridge import snakemake_compatible

@snakemake_compatible(mapping={
    'perm': 'wildcards.sample',
    'counts': 'input.counts',
    'metadata': 'input.metadata',
    'group': 'params.group',
    'out': 'output.grouped_mavs',
})
def main():
    p = argparse.ArgumentParser()
    p.add_argument("-perm", default=100)
    p.add_argument("--counts")
    p.add_argument("--metadata")
    p.add_argument("--out")
    p.add_argument("-group", default="Mutagen")

    args = p.parse_args()

    plt.rc("font", size=13)
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Liberation Sans"]

    N_PERM = int(args.perm)

    mav_counts = pd.read_csv(args.counts, sep="\t")

    metadata = pd.read_excel(
        args.metadata,
        sheet_name="Sample_description",
        skiprows=1,
    )

    mav_counts = mav_counts.merge(metadata, left_on="sample", right_on="Sample")

    print (mav_counts)

    group = args.group

    f, ax = plt.subplots(figsize=(mav_counts[group].nunique() * 2, 5))

    palette = sns.color_palette("colorblind", mav_counts[group].nunique())

    rng = np.random.default_rng(42)

    order = []
    for gi, (g, group_df) in enumerate(mav_counts.groupby(group)):

        tot = mav_counts[(mav_counts["permuted"] == 0) & (mav_counts[group] == g)][
            "total"
        ].sum()
        tot = group_df[(group_df["permuted"] == 0) & (group_df["is_multiallelic"])]["total"].sum()

        multiallelics = group_df[group_df["is_multiallelic"]]

        emp = multiallelics[multiallelics["permuted"] == 0]["count"].sum()
        prm = (
            multiallelics[multiallelics["permuted"] > 0]
            .groupby("permuted")["count"]
            .sum()
        )

        prm_mean = np.mean(prm.values) 
        enrich = round((emp + 1) / (np.mean(prm.values) + 1), 1)

        ax.scatter(
            rng.normal(0, 0.125, size=prm.shape[0]) + gi,
            (prm.values / tot) * 100,
            ec="w",
            label="random strain" if gi == 0 else None,
            color=palette[gi],
        )
        ax.scatter(
            gi,
            (emp / tot) * 100,
            color="k",
            marker="+",
            s=100,
            label="focal strain" if gi == 0 else None,
        )

        order.append(f"{g}\n({tot:,} SNVs)")

    ax.set_xticks(range(len(order)))
    ax.set_ylabel("% of SNVs with multi-allelic read support")
    ax.set_xlabel(group)
    ax.set_xticklabels(order)
    ax.legend(shadow=True, title="Source of read evidence", loc="center")
    sns.despine(ax=ax)

    f.tight_layout()
    f.savefig(args.out, dpi=300)

if __name__ == "__main__":
    main()
