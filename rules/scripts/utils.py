from collections import namedtuple


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
