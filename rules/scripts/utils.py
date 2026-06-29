from collections import namedtuple


def get_read_info(bam, chrom, start, end):

    res = []

    ReadInfo = namedtuple(
        "read",
        [
            "start",
            "end",
            "orientation",
            "base",
            "qual",
        ],
    )

    coords = []
    for ri, read in enumerate(bam.fetch(chrom, start, end)):
        # if we've counted over 150, move on
        # if ri > 150: break
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
        if read.reference_start <= start < read.reference_end:
            dist_from_start = start - read.reference_start
            dist_from_end = read.reference_end - end  # reference_end is exclusive
            is_near_boundary = (dist_from_start < 5) or (dist_from_end <= 5)
            if is_near_boundary:
                continue

        start_end = f"{read.reference_start}_{read.reference_end}"
        if start_end in coords:
            continue

        coords.append(start_end)
        # get aligned pairs
        pairs = read.get_aligned_pairs(matches_only=True)
        for query_pos, ref_pos in pairs:
            if ref_pos == start:
                # query_pos is the index in the read's sequence
                base = read.query_sequence[query_pos]
                qual = read.query_qualities[query_pos]
                orient = read.is_forward
                if qual < 20:
                    continue

                read_info = ReadInfo(
                    start=read.reference_start,
                    end=read.reference_end,
                    orientation=orient,
                    base=base,
                    qual=qual,
                )
                res.append(read_info)

    return res
