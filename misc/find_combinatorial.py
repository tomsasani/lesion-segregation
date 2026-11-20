from collections import Counter, defaultdict, namedtuple
import numpy as np
import pandas as pd
import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from pysam import AlignmentFile

plt.rc("font", size=12)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

ReadInfo = namedtuple(
        "read",
        [
            "name",
            "start",
            "end",
            "orientation",
            "base",
            "base_pos",
            "qual",
        ],
    )

def get_read_info(bam, chrom, start, end):

    res = []
    seen_read_names = set()
    coords = set()

    for ri, read in enumerate(bam.fetch(chrom, start, end)):
        if read.mapping_quality < 60:
            continue
        if not read.is_paired:
            continue
        if read.is_duplicate:
            continue
        if (
            read.is_unmapped
            or read.is_secondary
            or read.is_supplementary
            or read.mate_is_unmapped
        ):
            continue
        if read.reference_id != read.next_reference_id:
            continue

        # 5. Prevent double-counting if Read 1 and Read 2 overlap at this position
        if read.query_name in seen_read_names:
            continue

        r_start, r_end = read.reference_start, read.reference_end

        # Ensure the read actually overlaps the specific target `start` position
        if not (r_start <= start < r_end):
            continue

        dist_from_start = start - r_start
        dist_from_end = r_end - end  # reference_end is exclusive
        is_near_boundary = (dist_from_start < 5) or (dist_from_end <= 5)
        if is_near_boundary:
            continue
        
        start_end = (r_start, r_end)
        if start_end in coords:
            continue
        coords.add(start_end)

        # get aligned pairs
        pairs = read.get_aligned_pairs(matches_only=True)
        for query_pos, ref_pos in pairs:
            if ref_pos == start:
                # query_pos is the index in the read's sequence
                base = read.query_sequence[query_pos]
                qual = read.query_qualities[query_pos]
                orient = read.is_forward
                if qual < 20:
                    break
                read_info = ReadInfo(
                    name=read.query_name,
                    start=read.reference_start,
                    end=read.reference_end,
                    orientation=orient,
                    base=base,
                    base_pos=ref_pos,
                    qual=qual,
                )
                res.append(read_info)
                seen_read_names.add(read.query_name)
                break

    return res

IGV_PALETTE = [
    "#FCA35C",
    "#A5D2E3",
    "#CF1E1F",
    "#2871AC",
    "#000000",
]

def plot(
    read_info,
    v_start,
    v_end,
    ref_1 = "N",
    ref_2 = "N",
    fig_fh: str = "o.png",
):

    f, ax = plt.subplots(
        # 2,
        # figsize=(8, 10),
        # sharex=True,
        # sharey=True,
    )

    cmap = dict(
        zip(
            ["A", "C", "G", "T", "N"],
            IGV_PALETTE,
        )
    )

    thresh = 100 / len(read_info)

    rng = np.random.default_rng(42)
    read_count = 0

    fwd_base_combos = []
    rev_base_combos = []

    mid_point = (v_end + v_start) // 2
    for read_name in read_info.keys():
        ri = read_info[read_name]

        if rng.uniform() > thresh:
            continue

        if len(ri) != 2:
            continue

        # get start and end of this read
        starts = [read.start for read in ri]
        ends = [read.end for read in ri]
        orients = [read.orientation for read in ri]
        
        plotted = []
        for s, e, o in zip(starts, ends, orients):
            
            if (s, e, o) in plotted: continue

            adj_s, adj_e = s - mid_point, e - mid_point

            ax.plot(
                [adj_s, adj_e],
                [read_count, read_count],
                color="gainsboro" if o else "dimgrey",
                lw=4,
                zorder=-1,
                solid_capstyle="round",
            )
            base_combos = tuple(read.base for read in ri)
            if o:
                fwd_base_combos.append(base_combos)
            else:
                rev_base_combos.append(base_combos)
            
            for read in ri:

                base = read.base
                base_pos = read.base_pos
                base_qual = read.qual

                # calculate relative position
                rel_pos = base_pos - mid_point
                # ask whether this base is REF or ALT
                relevant_ref = ref_1 if rel_pos < 0 else ref_2
                is_ref = base == relevant_ref
                if not is_ref:
                    ax.scatter(base_pos - mid_point, read_count, color="#1d75bc", zorder=1, alpha=base_qual / 41, s=75)

            
            plotted.append((s, e, o))
        read_count += 1
    ax.set_yticks([])
    sns.despine(ax=ax, left=True)
    ax.set_xlabel("Position relative to mutated base (bp)")

    
    # if read_count < 15: return None
    legend = []
    legend_labels = []

    

    for li, l in enumerate((fwd_base_combos, rev_base_combos)):
        set_l = Counter(l).most_common()

        for combo, count in set_l:
            if len(combo) != 2:
                continue

            connection_line = Line2D(
                            [-0.5, 0.5],
                            [0],
                            color=(
                                "gainsboro" if li == 0 else "dimgrey"
                            ),  # Match your read track colors!
                            linewidth=3,
                            linestyle="-",
                            zorder=1,
                        )


            # 2. Create the first dot (positioned on the left)
            dot1 = Line2D(
                [-0.5],
                [0],
                marker="o",
                markersize=10,
                markerfacecolor="blue" if combo[0] == ref_1 else "red",#cmap[combo[0]],
                markeredgecolor="w",
                linestyle="none",
                zorder=2,
            )

            # 3. Create the second dot (positioned on the right)
            dot2 = Line2D(
                [0.5],
                [0],
                marker="o",
                markersize=10,
                markerfacecolor="blue" if combo[1] == ref_2 else "red",#cmap[combo[1]],
                markeredgecolor="w",
                linestyle="none",
                zorder=2,
            )

            # Bundle all three into the tuple
            legend.append((dot1, connection_line, dot2))
            legend_labels.append(f"{'-'.join(combo)} (n={count})")


    ax.legend(
        handles=legend,
        labels=legend_labels,
        handler_map={tuple: HandlerTuple(ndivide=None)},
        frameon=False,
        title="Allele support",
    )

    f.tight_layout()
    f.savefig(fig_fh, dpi=200)
    plt.close()


mutations = pd.read_csv("csv/mutations/combined.tsv", sep="\t")
metadata = pd.read_excel(
    "/scratch/ucgd/lustre-labs/quinlan/u1006375/lesion-segregation/data/41467_2020_15912_MOESM4_ESM.xlsx",
    sheet_name="Sample_description",
    skiprows=1,
).fillna(
    {
        "Normal": "CD0001b",
        "Mutagen": "none",
        "Drug.concentration": -1,
    }
)

metadata["Mutagen"] = metadata.apply(
    lambda row: row["Mutagen"] if row["Type"] != "mut.acc." else "None (M.A.)", axis=1
)

mutations = mutations.merge(metadata, left_on="sample", right_on="Sample")

dists = []
res = []
for (sample, mut), sample_df in tqdm.tqdm(mutations.groupby(["sample", "Mutagen"])):

    samfile = AlignmentFile(f"data/cram/{sample}.cram", "rb")

    for chrom, chrom_df in sample_df.groupby("chrom"):
        if chrom_df.shape[0] < 2: continue
        pos = chrom_df["start"].values
        # filter to low dist
        dist = [0] + list(pos[1:] - pos[:-1])
        chrom_df["dist"] = dist
        chrom_df = chrom_df.reset_index()
        for i, row in chrom_df.iterrows():

            combo = {}

            all_reads = defaultdict(list)

            # ignore probable collaterals
            # if row["dist"] < 20: continue
            if row["dist"] == 0: continue
            if row["dist"] > 750: continue
            dists.append(row["dist"])

            site_1, site_2 = chrom_df.iloc[i - 1].to_dict(), row.to_dict()
            ref_1, ref_2 = site_1["ref"], site_2["ref"]

            for site in (site_1["start"], site_2["start"]):
                info = get_read_info(samfile, chrom, site, site + 1)
                for r in info:
                    all_reads[r.name].append(r)
        
           
            combos = []
            for read in all_reads:
                read_info = all_reads[read]
                # make sure reads span both sites
                if len(read_info) < 2:
                    continue

                # read pairs overlap and sequence the same base twice
                read_info_filt = []
                seen_bases = []
                for r in read_info:
                    if (r.base_pos, r.base) in seen_bases:
                        continue
                    read_info_filt.append(r)
                    seen_bases.append((r.base_pos, r.base))
                # if len(read_info_filt) < 2: 
                #     continue
                combos.append(tuple(_v.base for _v in read_info_filt))
            uniq = Counter(combos).most_common()

            # keep track of all possible combinatorial
            # pairs to calculate fractions
            if len(uniq) < 3:
                continue
            else:
                if uniq[2][-1] < 2:
                    continue
                
                # print (uniq)
                plot(
                        all_reads, site_1["start"], site_2["start"], ref_1=ref_1, ref_2=ref_2, fig_fh=f"fig/combinatorial/{mut}.{sample}.{chrom}_{site_1['start']}.{ref_1}.{ref_2}.png"
                    )

            print (uniq)
            for (n1, n2), ct in uniq:
                res.append({"sample": sample, "chrom": site_1["chrom"], 
                          "pos_1": site_1["start"], "ref_1": site_1["ref"], "alt_1": site_1["alt"],
                          "pos_2": site_2["start"], "ref_2": site_2["ref"], "alt_2": site_2["alt"],
                          "nuc_1": n1, "nuc_2": n2, "count": ct, "mutagen": row["Mutagen"]})
            

res = pd.DataFrame(res)
res.to_csv("combinatorial.tsv", sep="\t", index=False)


