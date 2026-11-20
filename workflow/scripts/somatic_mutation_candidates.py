#!/usr/bin/env python3
"""
somatic_mutation_candidates.py

Somatic mutation candidates and their callable bps, from the trimmed mutation
table of mutation_tables.py (Methods, Somatic and germline mutations).

A candidate is a variant of the trimmed table that
    - passes MULTIALLELIC_SITE, IN_REPEATMASKER_MASK and SHARED_WITH_OTHER_LIBRARIES
    - is a singleton (VARIANT_MOLECULES = 1)
    - in the untrimmed reads: is a singleton (UNTRIMMED_VARIANT_MOLECULES = 1)
      and passes UNTRIMMED_MULTIALLELIC_SITE, UNTRIMMED_SHARED_WITH_OTHER_LIBRARIES
      and UNTRIMMED_ABOVE_LIBRARY_COVERAGE_P99
    - has UNTRIMMED_COVERAGE >= --min-cov
The somatic mutations are the candidates with SINGLETONS_ON_MOLECULE <= 2; the
column is kept and the cutoff is not applied here.

Output
    --candidates-out  one row per candidate, columns in COLUMNS; SBS96_CHANNEL
                      is the 96-type label, e.g. A[C>T]G
    --callable-out    the rows of --denominator (callable_bps.py) at --min-cov

Usage
    python somatic_mutation_candidates.py \\
        --mutation-table mutations_450b_trimmed_reads.tsv \\
        --denominator <results>/callable_bps/trimmed/denominator.tsv --min-cov 10 \\
        --candidates-out somatic_mutation_candidates.tsv \\
        --callable-out callable_bases.tsv
"""
import argparse
import sys

import numpy as np
import pandas as pd

FILTERS = ["MULTIALLELIC_SITE", "IN_REPEATMASKER_MASK", "SHARED_WITH_OTHER_LIBRARIES",
           "UNTRIMMED_MULTIALLELIC_SITE", "UNTRIMMED_SHARED_WITH_OTHER_LIBRARIES",
           "UNTRIMMED_ABOVE_LIBRARY_COVERAGE_P99"]
COLUMNS = ["LIBRARY", "ANIMAL", "CHROM", "POS", "REF", "ALT", "TRINUCLEOTIDE_CONTEXT", "SBS96_CHANNEL",
           "COVERAGE", "REPEATMASKER_CLASS", "ZMW", "SINGLETONS_ON_MOLECULE",
           "UNTRIMMED_SINGLETONS_ON_MOLECULE", "UNTRIMMED_COVERAGE",
           "UNTRIMMED_ABOVE_LIBRARY_COVERAGE_P99",
           "FWD_STRAND_SM", "REV_STRAND_SM", "FWD_STRAND_SX", "REV_STRAND_SX",
           "FWD_STRAND_NP", "REV_STRAND_NP", "FWD_POSITION_IN_READ", "REV_POSITION_IN_READ",
           "FWD_READ_LENGTH", "REV_READ_LENGTH"]
CUTOFF = {"4kb": 2, "10kb": 2}            # used only for the printed check
_COMP = str.maketrans("ACGT", "TGCA")


def sbs96(triplet, ref, alt):
    """96-type label, e.g. A[C>T]G, written from the pyrimidine of the base pair."""
    if ref in "AG":
        triplet, ref, alt = triplet.translate(_COMP)[::-1], ref.translate(_COMP), alt.translate(_COMP)
    return f"{triplet[0]}[{ref}>{alt}]{triplet[2]}"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mutation-table", required=True)
    p.add_argument("--denominator", required=True)
    p.add_argument("--min-cov", type=int, default=10)
    p.add_argument("--candidates-out", required=True)
    p.add_argument("--callable-out", required=True)
    args = p.parse_args()

    t = pd.read_csv(args.mutation_table, sep="\t", dtype=str, keep_default_na=False)
    num = lambda c: t[c].astype(int)
    keep = ((t[FILTERS].astype(int).sum(axis=1) == 0) & (num("VARIANT_MOLECULES") == 1)
            & (num("UNTRIMMED_VARIANT_MOLECULES") == 1) & (num("UNTRIMMED_COVERAGE") >= args.min_cov))
    c = t[keep].copy()
    if (c.TRINUCLEOTIDE_CONTEXT.str[1] != c.REF).any():
        sys.exit("[error] reference base is not the centre of the trinucleotide context")
    c["SBS96_CHANNEL"] = [sbs96(x, r, a) for x, r, a in zip(c.TRINUCLEOTIDE_CONTEXT, c.REF, c.ALT)]
    c[COLUMNS].to_csv(args.candidates_out, sep="\t", index=False)

    d = pd.read_csv(args.denominator, sep="\t")
    d = d[d.floor == args.min_cov]
    if d.empty:
        sys.exit(f"[error] {args.denominator} has no rows at minimum coverage {args.min_cov}")
    d.to_csv(args.callable_out, sep="\t", index=False)

    cut = np.where(c.LIBRARY.str.endswith("10kb"), CUTOFF["10kb"], CUTOFF["4kb"])
    n_cut = int((c.SINGLETONS_ON_MOLECULE.astype(int) <= cut).sum())
    print(f"{len(c):,} candidates at minimum coverage {args.min_cov} -> {args.candidates_out}\n"
          f"{n_cut:,} with SINGLETONS_ON_MOLECULE <= {CUTOFF['4kb']} (4 kb) / {CUTOFF['10kb']} (10 kb)\n"
          f"callable bases at minimum coverage {args.min_cov} -> {args.callable_out}")


if __name__ == "__main__":
    main()
