#!/usr/bin/env python3
"""
read_and_base_filters.py

Applies the read filters and the base filters (Methods, Variant calling and
filtering) to the strand consensus reads of a SAM stream, by setting base
qualities to 0. Reads SAM on stdin, writes SAM to stdout. variant_calling.py
then uses only the bases with quality >= its -q (40).

Read filters (every base of the read is set to 0):
    read length outside [--min-length, --max-length]
    passes (np tag) below --min-np
    divergence 1 - M / (M + X + I + D) above --threshold, from the CIGAR
    (= and X operations, as written by pbmm2 --preset CCS)

Base filters (the base is set to 0):
    within --padding bases of an insertion or deletion
    within --edge bases of either read end: 0 for the untrimmed reads, 450 for
        the trimmed reads (Supplementary Figs. 3-4)
    sm (subreads matching the base) below --min-sm
    sm / (sm + sx) at or below --min-sm-ratio

Usage:
    samtools view -q 60 -h <sample>.sorted.bam T2T1 \
    | python read_and_base_filters.py -i 3 -d 0.01 -e 450 -np 5 -sm 5 -r 0.8 \
          -lmin <min length> -lmax <max length> \
    | python variant_calling.py -i - ...
"""
import sys, argparse, pysam  # type: ignore


def read_divergence(aln):
    """1 - (= bases) / (= + X + I + D bases) of the alignment, in [0, 1]."""
    stats   = aln.get_cigar_stats()[0]
    n_match = stats[7]   # = exact matches
    n_mis   = stats[8]   # X mismatches
    n_ins   = stats[1]   # I insertions (read bases)
    n_del   = stats[2]   # D deletions  (ref bases)
    aln_len = n_match + n_mis + n_ins + n_del
    if aln_len == 0:
        return 0.0
    return 1.0 - (n_match / aln_len)


def mask_alignment(aln, pad, div_threshold, edge, np_min, sm_min,
                   sm_ratio_min, min_length, max_length):
    if aln.is_unmapped or not aln.cigartuples:
        return aln

    qlen = aln.query_length
    quals = aln.query_qualities
    if quals is None or qlen == 0:
        return aln

    # ------------------------------------------------------------------
    # Read filters
    # ------------------------------------------------------------------
    # Read length
    if min_length > 0 or max_length > 0:
        rlen = aln.infer_read_length() or qlen
        if (min_length > 0 and rlen < min_length) or \
           (max_length > 0 and rlen > max_length):
            aln.query_qualities = [0] * qlen
            return aln

    # Number of passes (np tag)
    if np_min > 0:
        try:
            np_passes = int(aln.get_tag("np"))
        except KeyError:
            np_passes = 0
        if np_passes < np_min:
            aln.query_qualities = [0] * qlen
            return aln

    # Divergence from the reference
    if read_divergence(aln) > div_threshold:
        aln.query_qualities = [0] * qlen
        return aln

    # ------------------------------------------------------------------
    # Base filters
    # ------------------------------------------------------------------
    mask = [False] * qlen
    qpos = 0

    # Distance to an insertion or deletion
    for op, length in aln.cigartuples:
        if op in (0, 7, 8):  # M, =, X
            qpos += length
        elif op == 1:  # I
            start = max(0, qpos - pad)
            end = min(qlen, qpos + length + pad)
            for i in range(start, end):
                mask[i] = True
            qpos += length
        elif op == 2:  # D
            start = max(0, qpos - pad)
            end = min(qlen, qpos + pad)
            for i in range(start, end):
                mask[i] = True
        elif op == 4:  # S (soft clip)
            qpos += length

    # Read ends
    for i in range(min(edge, qlen)):
        mask[i] = True
        mask[-(i + 1)] = True

    # Subreads matching the base: sm and sm / (sm + sx)
    if sm_min > 0 or sm_ratio_min > 0:
        sms = aln.get_tag("sm")
        sxs = aln.get_tag("sx")
        # sm and sx in the order of SEQ for reads aligned to the minus strand
        if aln.is_reverse:
            sms = sms[::-1]
            sxs = sxs[::-1]

        for i in range(qlen):
            sm_val = sms[i]
            sx_val = sxs[i]

            # Minimum sm
            if sm_min > 0 and sm_val < sm_min:
                mask[i] = True
                continue

            # Minimum sm / (sm + sx), strictly above the threshold
            if sm_ratio_min > 0:
                total = sm_val + sx_val
                if total == 0 or (sm_val / total) <= sm_ratio_min:
                    mask[i] = True

    # Set the filtered bases to quality 0
    if any(mask):
        for i, flag in enumerate(mask):
            if flag:
                quals[i] = 0
        aln.query_qualities = quals

    return aln


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Set base qualities to 0 by the read filters (length, passes, "
            "divergence) and the base filters (indels, read ends, sm, sm/(sm+sx))."
        )
    )
    parser.add_argument("-i", "--padding", type=int, default=3,
                        help="bases on each side of an insertion or deletion "
                             "(default: 3)")
    parser.add_argument("-d", "--threshold", type=float, default=0.01,
                        help="maximum divergence of a read (default: 0.01)")
    parser.add_argument("-e", "--edge", type=int, default=0,
                        help="bases at each read end (default: 0)")
    parser.add_argument("-np", "--min-np", dest="np_min", type=int, default=5,
                        help="minimum number of passes (np tag); 0 = off "
                             "(default: 5)")
    parser.add_argument("-sm", "--min-sm", dest="sm_min", type=int, default=5,
                        help="minimum sm of a base; 0 = off (default: 5)")
    parser.add_argument("-r", "--min-sm-ratio", dest="sm_ratio_min",
                        type=float, default=0.8,
                        help="sm / (sm + sx) of a base must be above this; "
                             "0 = off (default: 0.8)")
    parser.add_argument("-lmin", "--min-length", dest="min_length", type=int,
                        default=0,
                        help="minimum read length; 0 = off (default: 0)")
    parser.add_argument("-lmax", "--max-length", dest="max_length", type=int,
                        default=0,
                        help="maximum read length; 0 = off (default: 0)")
    args = parser.parse_args()

    with pysam.AlignmentFile(sys.stdin, "r") as infile, \
         pysam.AlignmentFile(sys.stdout, "w", header=infile.header) as outfile:
        for aln in infile.fetch(until_eof=True):
            aln = mask_alignment(
                aln,
                args.padding, args.threshold, args.edge,
                args.np_min, args.sm_min, args.sm_ratio_min,
                args.min_length, args.max_length,
            )
            outfile.write(aln)


if __name__ == "__main__":
    main()
