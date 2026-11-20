#!/usr/bin/env python3
"""
reference_context_counts.py

Trinucleotide counts and lengths of the T2T autosomes, over the whole
autosomes and outside the masked sites (Supplementary Table 6).

A site is counted when it and both neighbours are A, C, G or T. Its context is
the reference triplet centred on it, read from the + strand, as in the
`triplet` column of the cntx and variant tables. All 64 triplets are reported.

Output
    --counts-out   one row per triplet: triplet, autosomes, outside_mask
    --lengths-out  one row per autosome plus a row `autosomes` with the totals:
                   chrom, length_bp, acgt_bp, non_acgt_bp, mask_bp,
                   outside_mask_bp, outside_mask_acgt_bp

Usage
    python reference_context_counts.py --ref b6_t2t_mm39XYM.fa --chroms T2T1 ... T2T19 \\
        --mask excluded_regions.bed --counts-out t2t_reference_context_counts.tsv \\
        --lengths-out t2t_autosome_lengths.tsv
"""
import argparse
import itertools
import sys

import numpy as np
import pandas as pd
import pysam

BASES = "ACGT"
TRIPLETS = ["".join(t) for t in itertools.product(BASES, repeat=3)]
CODE = np.full(256, 4, dtype=np.uint8)
for i, b in enumerate(BASES):
    CODE[ord(b)] = i
    CODE[ord(b.lower())] = i


def load_mask(path, chroms):
    bed = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"],
                      dtype={"chrom": str})
    return {c: g[["start", "end"]].to_numpy(np.int64) for c, g in bed[bed.chrom.isin(chroms)].groupby("chrom")}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref", required=True)
    p.add_argument("--chroms", nargs="+", required=True)
    p.add_argument("--mask", required=True, help="masked sites BED (0-based half-open)")
    p.add_argument("--counts-out", required=True)
    p.add_argument("--lengths-out", required=True)
    args = p.parse_args()

    mask = load_mask(args.mask, args.chroms)
    whole = np.zeros(64, np.int64)
    outside = np.zeros(64, np.int64)
    rows = []
    with pysam.FastaFile(args.ref) as fa:
        for c in args.chroms:
            code = CODE[np.frombuffer(fa.fetch(c).encode(), dtype=np.uint8)]
            L = len(code)
            masked = np.zeros(L, bool)
            for s, e in mask.get(c, np.empty((0, 2), np.int64)):
                masked[s:e] = True
            acgt = code < 4

            ok = acgt[:-2] & acgt[1:-1] & acgt[2:]                 # centre = position i + 1
            idx = (code[:-2] << 4) | (code[1:-1] << 2) | code[2:]   # uint8; valid where ok
            whole += np.bincount(idx[ok], minlength=64)
            outside += np.bincount(idx[ok & ~masked[1:-1]], minlength=64)

            rows.append(dict(chrom=c, length_bp=L, acgt_bp=int(acgt.sum()), non_acgt_bp=int((~acgt).sum()),
                             mask_bp=int(masked.sum()), outside_mask_bp=int((~masked).sum()),
                             outside_mask_acgt_bp=int((acgt & ~masked).sum())))
            sys.stderr.write(f"{c}: {L:,} bp\n")

    lengths = pd.DataFrame(rows)
    total = lengths.drop(columns="chrom").sum()
    lengths = pd.concat([lengths, pd.DataFrame([dict(chrom="autosomes", **total.to_dict())])], ignore_index=True)
    lengths.to_csv(args.lengths_out, sep="\t", index=False)

    counts = pd.DataFrame(dict(triplet=TRIPLETS, autosomes=whole, outside_mask=outside))
    counts.to_csv(args.counts_out, sep="\t", index=False)
    sys.stderr.write(f"autosomes {int(total.length_bp):,} bp; contexts counted: {int(whole.sum()):,} whole, "
                     f"{int(outside.sum()):,} outside the mask\n")


if __name__ == "__main__":
    main()
