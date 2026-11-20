#!/usr/bin/env python3
"""
read_table.py

One row per HiFi strand consensus read of a library, from the primary records
of the aligned BAM (unmapped reads included), saved as a compressed .npz. Used
for the sequencing summary and the read and molecule filters of each library
(Supplementary Tables 2 and 3).

Output: .npz with one entry per read
        zmw      int64   ZMW number
        strand   int8    0 = fwd, 1 = rev
        np       int16   number of passes (np tag)
        length   int32   read length (bases in SEQ)
        mapped   bool    has a primary alignment
        mapq     uint8   mapping quality (0 if unmapped)
        contig   int8    0 = unmapped, 1 = T2T autosome (T2T1-T2T19), 2 = other contig
        chrom    int8    autosome number (1-19), 0 if not on a T2T autosome
        n_eq     int32   matched bases in the alignment (CIGAR =)
        aln_len  int32   CIGAR = + X + I + D; divergence = 1 - n_eq / aln_len
        start    int32   alignment start on the reference, 0-based (-1 if unmapped)
        end      int32   alignment end on the reference, 0-based exclusive (-1 if unmapped)
        sa       bool    has a supplementary alignment (SA tag)
        plus: sample, bam, n_reads

With --flagstat, the number of reads is checked against the flagstat "primary"
line. --region restricts the run to one contig (no flagstat check).

Usage:  python read_table.py --bam <sample>.sorted.bam --out <sample>.reads.npz \
            --flagstat <sample>.flagstat.txt --threads 3
"""
import argparse
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

# One line per primary record: zmw, strand, flag, contig, mapq, read length,
# np, CIGAR = bases, CIGAR = + X + I + D, position, reference span, SA tag present
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
    print q[2], q[4], $2, $3, $5, length($10), np, eq, aln, $4, span, sa
}
'''
COLS = ["zmw", "strand", "flag", "rname", "mapq", "length", "np", "n_eq", "aln_len", "pos", "span", "sa"]
CHUNK = 4_000_000


def flagstat_primary(path):
    with open(path) as f:
        m = re.search(r"^(\d+) \+ \d+ primary$", f.read(), re.M)
    if m is None:
        sys.exit(f"[error] no 'primary' line in {path}")
    return int(m.group(1))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--flagstat", default=None)
    p.add_argument("--threads", type=int, default=3)
    p.add_argument("--region", default=None)
    args = p.parse_args()

    sample = os.path.basename(args.bam).split(".")[0]
    view = ["samtools", "view", "-@", str(args.threads), "-F", "0x900",
            "--keep-tag", "np,SA", args.bam]
    if args.region:
        view.append(args.region)

    sam = subprocess.Popen(view, stdout=subprocess.PIPE)
    awk = subprocess.Popen(["awk", AWK], stdin=sam.stdout, stdout=subprocess.PIPE)
    sam.stdout.close()

    autosome_num = {f"T2T{i}": i for i in range(1, 20)}
    parts = {k: [] for k in ("zmw", "strand", "np", "length", "mapped", "mapq", "contig",
                             "chrom", "n_eq", "aln_len", "start", "end", "sa")}
    reader = pd.read_csv(awk.stdout, sep="\t", header=None, names=COLS,
                         dtype={"zmw": np.int64, "strand": str, "flag": np.int32,
                                "rname": str, "mapq": np.int16, "length": np.int32,
                                "np": np.int16, "n_eq": np.int32, "aln_len": np.int32,
                                "pos": np.int64, "span": np.int64, "sa": np.int8},
                         chunksize=CHUNK)
    n = 0
    for ch in reader:
        strand = ch["strand"].to_numpy()
        if not np.isin(strand, ["fwd", "rev"]).all():
            sys.exit(f"[error] unexpected strand field in {args.bam}")
        if (ch["np"].to_numpy() < 0).any():
            sys.exit(f"[error] read without an np tag in {args.bam}")
        mapped = (ch["flag"].to_numpy() & 4) == 0
        chrom  = np.where(mapped, ch["rname"].map(autosome_num).fillna(0).to_numpy(), 0)
        contig = np.where(~mapped, 0, np.where(chrom > 0, 1, 2))
        parts["zmw"].append(ch["zmw"].to_numpy())
        parts["strand"].append((strand == "rev").astype(np.int8))
        parts["np"].append(ch["np"].to_numpy())
        parts["length"].append(ch["length"].to_numpy())
        parts["mapped"].append(mapped)
        parts["mapq"].append(np.where(mapped, ch["mapq"].to_numpy(), 0).astype(np.uint8))
        parts["contig"].append(contig.astype(np.int8))
        parts["chrom"].append(chrom.astype(np.int8))
        parts["n_eq"].append(ch["n_eq"].to_numpy())
        parts["aln_len"].append(ch["aln_len"].to_numpy())
        start = np.where(mapped, ch["pos"].to_numpy() - 1, -1)
        parts["start"].append(start.astype(np.int32))
        parts["end"].append(np.where(mapped, start + ch["span"].to_numpy(), -1).astype(np.int32))
        parts["sa"].append(ch["sa"].to_numpy().astype(bool))
        n += len(ch)
        print(f"  {sample}: {n:,} reads", file=sys.stderr, flush=True)

    awk.wait()
    sam.wait()
    if sam.returncode != 0 or awk.returncode != 0:
        sys.exit(f"[error] samtools/awk exited with {sam.returncode}/{awk.returncode}")

    data = {k: np.concatenate(v) if v else np.array([]) for k, v in parts.items()}

    if args.flagstat and not args.region:
        expected = flagstat_primary(args.flagstat)
        if n != expected:
            sys.exit(f"[error] {sample}: {n:,} reads vs flagstat primary {expected:,}")
        print(f"  {sample}: matches flagstat primary ({expected:,})", file=sys.stderr)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    np.savez_compressed(args.out, sample=sample, bam=args.bam, n_reads=n, **data)
    print(f"wrote {args.out}  ({n:,} reads)", file=sys.stderr)


if __name__ == "__main__":
    main()
