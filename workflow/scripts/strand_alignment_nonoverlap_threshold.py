#!/usr/bin/env python3
"""
strand_alignment_nonoverlap_threshold.py

Threshold of the strand alignment non-overlap filter (Methods, Molecule
filters; Supplementary Table 3): a percentile, over the molecules of one
library, of the number of sites covered by only one of the two alignments of
its strand consensus reads,

    non-overlap = |union of the two reference spans| - |intersection|

Passed to variant_calling.py as --max-strand-nonoverlap.

Molecules used: the forward and reverse reads of a ZMW, primary alignments with
MAPQ >= --min-mapq on the same contig of --regions, that pass the read filters
(np >= --np-min, length within --bounds, divergence <= --div-max), the strand
length difference filter (--lendiv) and have no supplementary alignment.

Output: line 1 of --out (or stdout): the threshold, in bases, rounded down.
        --stats: the percentile table and molecule counts as comment lines.

Usage:  python strand_alignment_nonoverlap_threshold.py -b <sample>.sorted.bam \
            --bounds <sample>.read_length_bounds.txt \
            --lendiv <sample>.strand_length_difference.txt \
            --regions T2T1,T2T2,...,T2T19 -o <sample>.strand_alignment_nonoverlap.txt \
            --pct 99 --stats
"""
import argparse
import subprocess
import sys

import numpy as np
import pandas as pd

QS = [50, 90, 95, 99, 99.5, 99.9]

# One line per alignment: zmw, strand, contig, start, end (0-based), read
# length, np, CIGAR = bases, CIGAR = + X + I + D, SA tag present
AWK = r'''
BEGIN { OFS = "\t" }
{
    split($1, q, "/")
    np = -1; sa = 0
    for (i = 12; i <= NF; i++) {
        if ($i ~ /^np:i:/) np = substr($i, 6)
        else if ($i ~ /^SA:Z:/) sa = 1
    }
    eq = 0; aln = 0; span = 0
    nn = split($6, n, /[^0-9]+/)
    split($6, op, /[0-9]+/)
    for (k = 1; k < nn; k++) {
        o = op[k + 1]
        if (o == "=") { eq += n[k]; aln += n[k]; span += n[k] }
        else if (o == "X" || o == "D") { aln += n[k]; span += n[k] }
        else if (o == "I") aln += n[k]
        else if (o == "M" || o == "N") span += n[k]
    }
    print q[2], q[4], $3, $4 - 1, $4 - 1 + span, length($10), np, eq, aln, sa
}
'''
COLS = ["zmw", "strand", "contig", "start", "end", "length", "np", "n_eq", "aln_len", "sa"]
DTYPES = {"zmw": np.int64, "strand": "category", "contig": "category", "start": np.int64,
          "end": np.int64, "length": np.int32, "np": np.int16, "n_eq": np.int32,
          "aln_len": np.int32, "sa": np.int8}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-b", "--bam", required=True)
    p.add_argument("--bounds", required=True, help="file with 'lmin lmax' on line 1")
    p.add_argument("--lendiv", required=True, help="file with max |flen-rlen| on line 1")
    p.add_argument("--regions", required=True, help="comma-separated contigs")
    p.add_argument("-o", "--out", default="-", help="output path, or '-' for stdout")
    p.add_argument("--pct", type=float, default=99.0,
                   help="percentile of non-overlap to report (default: 99)")
    p.add_argument("--np-min", type=int, default=5)
    p.add_argument("--div-max", type=float, default=0.01)
    p.add_argument("--min-mapq", type=int, default=60)
    p.add_argument("--threads", type=int, default=3)
    p.add_argument("--stats", action="store_true",
                   help="also write the percentile table and counts as comment lines")
    args = p.parse_args()

    lmin, lmax = map(int, open(args.bounds).readline().split()[:2])
    lendiv = int(open(args.lendiv).readline().split()[0])

    view = ["samtools", "view", "-@", str(args.threads), "-F", "0x904",
            "-q", str(args.min_mapq), "--keep-tag", "np,SA", args.bam] + args.regions.split(",")
    sam = subprocess.Popen(view, stdout=subprocess.PIPE)
    awk = subprocess.Popen(["awk", AWK], stdin=sam.stdout, stdout=subprocess.PIPE)
    sam.stdout.close()
    r = pd.read_csv(awk.stdout, sep="\t", header=None, names=COLS, dtype=DTYPES)
    awk.wait()
    sam.wait()
    if sam.returncode != 0 or awk.returncode != 0:
        sys.exit(f"[error] samtools/awk exited with {sam.returncode}/{awk.returncode}")
    if (r.np < 0).any():
        sys.exit(f"[error] read without an np tag in {args.bam}")

    # Read filters
    r["ok"] = ((r.np >= args.np_min) & (r.length >= lmin) & (r.length <= lmax)
               & (1 - r.n_eq / r.aln_len.clip(lower=1) <= args.div_max))
    keep = ["zmw", "contig", "start", "end", "length", "sa", "ok"]
    m = r.loc[r.strand == "fwd", keep].merge(r.loc[r.strand == "rev", keep],
                                             on=["zmw", "contig"], suffixes=("_f", "_r"))
    n_reads, n_pairs = len(r), len(m)
    del r
    # Molecule filters: strand length difference, chimeric molecules
    lend = (m.length_f - m.length_r).abs() <= lendiv
    whole = m.ok_f & m.ok_r
    sa = (m.sa_f > 0) | (m.sa_r > 0)
    u = m[lend & whole & ~sa]
    if u.empty:
        sys.exit(f"[error] no usable molecules in {args.bam}")
    nonov = ((np.maximum(u.end_f, u.end_r) - np.minimum(u.start_f, u.start_r))
             - np.maximum(0, np.minimum(u.end_f, u.end_r) - np.maximum(u.start_f, u.start_r)))
    nonov = nonov.to_numpy(dtype=np.float64)
    thr = int(np.floor(np.percentile(nonov, args.pct)))

    sys.stderr.write(
        f"[{args.bam}] reads={n_reads:,} pairs={n_pairs:,} usable={len(u):,}\n"
        f"  p{args.pct:g} of strand non-overlap = {thr} bp; "
        f"drops {(nonov > thr).sum():,} of {len(u):,} usable molecules "
        f"({100 * (nonov > thr).mean():.2f}%)\n")

    lines = [f"{thr}"]
    if args.stats:
        lines.append("# percentile\tnonoverlap_bp")
        for q in QS:
            lines.append(f"# {q:g}\t{np.percentile(nonov, q):g}")
        lines.append(f"# reads\t{n_reads}")
        lines.append(f"# pairs_same_contig\t{n_pairs}")
        lines.append(f"# pairs_failing_length_difference\t{(~lend).sum()}")
        lines.append(f"# pairs_failing_whole_read_filters\t{(lend & ~whole).sum()}")
        lines.append(f"# pairs_with_supplementary\t{(lend & whole & sa).sum()}")
        lines.append(f"# usable_molecules\t{len(u)}")
        lines.append(f"# zero_nonoverlap_pct\t{100 * (nonov == 0).mean():.4f}")
        lines.append(f"# dropped_at_threshold\t{(nonov > thr).sum()}")

    text = "\n".join(lines) + "\n"
    if args.out == "-":
        sys.stdout.write(text)
    else:
        with open(args.out, "w") as fh:
            fh.write(text)


if __name__ == "__main__":
    main()
