#!/usr/bin/env python
"""
subread_error_sample.py

Random sample of the per-base sm (subreads matching the consensus base) and sx
(subreads mismatching it) of one library's strand consensus reads, for the
error rate per subread (Methods, Estimation of the sequencing error rate):
sum(sx) / sum(sm + sx).

Sampling: a random site of the reference (uniform over contigs >= MIN_CONTIG_LEN,
weighted by length) is looked up in the BAM index and the next READS_PER_LOCUS
records are read, until --n-reads reads are kept. Primary alignments only, each
read at most once, reads of at least 2 * EDGE_TRIM + --k bases; --k positions
per read, drawn without replacement, more than EDGE_TRIM bases from either end.

Output: --out (.npz) with arrays sm, sx, np (one entry per sampled base).

Usage:  python subread_error_sample.py --bam <sample>.sorted.bam \
            --out <sample>.sm_sx.npz --n-reads 100000 --k 100 --seed 0
"""

import argparse
import os
import random

import numpy as np
import pysam

MIN_CONTIG_LEN  = 1_000_000   # contigs sampled: at least this long
READS_PER_LOCUS = 100         # records read per index lookup
EDGE_TRIM       = 450         # bases excluded at each end of the read, as in the trimmed reads


def sample_bam(bam_path, n_reads, k, seed=0,
               reads_per_locus=READS_PER_LOCUS, edge_trim=EDGE_TRIM):
    """Draw `n_reads` primary alignments at random and `k` per-base sm/sx values from each."""
    rng = random.Random(seed)

    bam = pysam.AlignmentFile(bam_path, "rb")
    contigs = [(c, l) for c, l in zip(bam.references, bam.lengths) if l >= MIN_CONTIG_LEN]
    if not contigs:
        raise RuntimeError(f"no contigs >= {MIN_CONTIG_LEN} bp in {bam_path}")
    cum   = np.cumsum([l for _, l in contigs])
    total = int(cum[-1])

    min_len = 2 * edge_trim + k

    seen = set()
    sm_chunks, sx_chunks, np_vals = [], [], []
    n_kept = n_loci = 0

    while n_kept < n_reads:
        x = rng.randrange(total)
        i = int(np.searchsorted(cum, x, side="right"))
        contig = contigs[i][0]
        pos    = x - (int(cum[i - 1]) if i else 0)
        n_loci += 1

        consumed = 0
        for read in bam.fetch(contig, pos):
            if consumed >= reads_per_locus or n_kept >= n_reads:
                break
            consumed += 1

            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            name = read.query_name
            if name in seen:
                continue
            try:
                sm_tag = read.get_tag("sm")
                sx_tag = read.get_tag("sx")
                np_val = read.get_tag("np")
            except KeyError:
                continue

            length = len(sm_tag)
            if length < min_len:
                continue

            seen.add(name)
            idx = np.fromiter(rng.sample(range(edge_trim, length - edge_trim), k),
                              dtype=np.int64, count=k)
            sm_chunks.append(np.frombuffer(sm_tag, dtype=np.uint8)[idx])
            sx_chunks.append(np.frombuffer(sx_tag, dtype=np.uint8)[idx])
            np_vals.append(np_val)
            n_kept += 1

    bam.close()
    return {
        "sm": np.concatenate(sm_chunks),
        "sx": np.concatenate(sx_chunks),
        "np": np.repeat(np.asarray(np_vals, dtype=np.uint16), k),
        "n_reads": n_kept,
        "n_loci":  n_loci,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam",     required=True, help="aligned, indexed BAM of the library")
    ap.add_argument("--out",     required=True, help="output .npz")
    ap.add_argument("--n-reads", type=int, default=100_000)
    ap.add_argument("--k",       type=int, default=100, help="positions sampled per read")
    ap.add_argument("--seed",    type=int, default=0)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    res = sample_bam(args.bam, args.n_reads, args.k, seed=args.seed)
    np.savez_compressed(args.out, sm=res["sm"], sx=res["sx"], np=res["np"])
    eps = res["sx"].sum() / (res["sm"].astype(np.int64).sum() + res["sx"].sum())
    print(f"{args.bam}: {res['n_reads']} reads from {res['n_loci']} loci, "
          f"{res['sm'].size} bases, eps={eps:.6f} -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
