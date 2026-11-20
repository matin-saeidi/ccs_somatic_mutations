#!/usr/bin/env python3
"""
strand_length_difference_threshold.py

Threshold of the strand length difference filter (Methods, Molecule filters;
Supplementary Table 3): a percentile of |forward read length - reverse read
length| over the molecules of one library. Passed to variant_calling.py as
--max-len-div.

Input:  <sample>.np_len.txt (columns zmw, strand, np, read_length, with a
        header). Molecules with exactly one forward and one reverse strand
        consensus read are used.
Output: line 1 of --out (or stdout): the threshold, in bases.
        --relative: the threshold as a fraction of the mean read length of the
        molecule (for variant_calling.py --max-len-div-frac).
        --stats: the percentile table as comment lines.

Usage:  python strand_length_difference_threshold.py -i <sample>.np_len.txt \
            -o <sample>.strand_length_difference.txt --pct 99 --stats
"""
import argparse
import sys

import numpy as np
import pandas as pd

QS = [50, 75, 90, 95, 97.5, 99, 99.5, 99.9]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--np-len", required=True,
                   help="<sample>.np_len.txt with columns zmw/strand/np/read_length")
    p.add_argument("-o", "--out", default="-", help="output path, or '-' for stdout")
    p.add_argument("--pct", type=float, default=99.0,
                   help="percentile of |flen-rlen| to report (default: 99)")
    p.add_argument("--relative", action="store_true",
                   help="report the threshold as a fraction of the mean read length")
    p.add_argument("--stats", action="store_true",
                   help="also write the percentile table as comment lines")
    args = p.parse_args()

    df = pd.read_csv(args.np_len, sep="\t",
                     usecols=["zmw", "strand", "read_length"],
                     dtype={"zmw": np.int64, "strand": "category",
                            "read_length": np.int32})
    n_rows = len(df)

    # Molecules (ZMWs) with exactly one forward and one reverse read
    cnt = df.groupby(["zmw", "strand"], observed=True).size().unstack(fill_value=0)
    for s in ("fwd", "rev"):
        if s not in cnt:
            sys.exit(f"[error] no '{s}' rows in {args.np_len}")
    good = cnt.index[(cnt["fwd"] == 1) & (cnt["rev"] == 1)]
    df = df[df.zmw.isin(good)]

    f = df[df.strand == "fwd"].set_index("zmw")["read_length"]
    r = df[df.strand == "rev"].set_index("zmw")["read_length"]
    j = pd.concat([f.rename("f"), r.rename("r")], axis=1, join="inner")
    if j.empty:
        sys.exit(f"[error] no paired molecules in {args.np_len}")

    dl = (j.f - j.r).abs().to_numpy(dtype=np.float64)
    ln = ((j.f + j.r) / 2).to_numpy(dtype=np.float64)
    vals = dl / ln if args.relative else dl
    thr = float(np.percentile(vals, args.pct))

    sys.stderr.write(
        f"[{args.np_len}] rows={n_rows:,} paired_molecules={len(j):,} "
        f"median_len={np.median(ln):,.0f}\n"
        f"  p{args.pct:g} of |flen-rlen| = "
        f"{thr:.6f}{'' if args.relative else ' bp'}"
        f"{' (fraction of mean read length)' if args.relative else ''}\n")

    lines = [f"{thr:.6f}" if args.relative else f"{thr:.0f}"]
    if args.stats:
        lines.append("# percentile\tabsolute_bp\trelative_frac")
        for q in QS:
            lines.append(f"# {q}\t{np.percentile(dl, q):.0f}\t"
                         f"{np.percentile(dl / ln, q):.6f}")
        lines.append(f"# paired_molecules\t{len(j)}\t")
        lines.append(f"# median_read_length\t{np.median(ln):.0f}\t")

    text = "\n".join(lines) + "\n"
    if args.out == "-":
        sys.stdout.write(text)
    else:
        with open(args.out, "w") as fh:
            fh.write(text)


if __name__ == "__main__":
    main()
