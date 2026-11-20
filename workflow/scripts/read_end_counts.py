#!/usr/bin/env python3
"""
read_end_counts.py

Callable bps and variants by distance from the read end, for the errors at the
edges of consensus reads (Methods, Errors at the edges of consensus reads;
Supplementary Figs. 3-4), in the untrimmed reads.

A callable bp is a site of a molecule where both strand consensus reads pass
the read, molecule and base filters and carry the same base (as in
variant_calling.py, whose get_ssdna_info is used here). Each is counted at its
distance d from an end of a strand consensus read, from that read's own 5' end:
    fwd_5 = fwd qpos          fwd_3 = fwd_len - 1 - fwd qpos
    rev_5 = rev qpos          rev_3 = rev_len - 1 - rev qpos
and by site:
    repeat bit   site in --repeat-bed
    shared bit   site in --shared-bed
    coverage     coverage of the site (cntx informative_cov), 0 .. --cov-max,
                 the last bin holding everything >= --cov-max

Subcommands
    count   one chromosome of one library, from the SAM stream of
            read_and_base_filters.py:
                samtools view -q 60 -h <sample>.sorted.bam T2T1 \\
                | python read_and_base_filters.py -i 3 -d 0.01 -e 0 \\
                      -np 5 -sm 5 -r 0.8 -lmin LMIN -lmax LMAX \\
                | python read_end_counts.py count -i - --chrom T2T1 \\
                      --ref b6_t2t_mm39XYM.fa --cntx T2T1.ext0.qual40.cntx.bed.gz \\
                      --repeat-bed repeat_regions.bed --shared-bed shared_sites.bed \\
                      --max-len-div N --drop-supplementary --max-strand-nonoverlap M \\
                      -o T2T1.npz
    merge   sum the chromosomes of one library:
                python read_end_counts.py merge -o <sample>.npz T2T1.npz T2T2.npz ...

Output arrays (npz)
    ends        ["fwd_5", "fwd_3", "rev_5", "rev_3"]
    win         int64 (4 ends, NCLS, 2 repeat, 2 shared, NCOV, --window)
                d < --window, by read-length class of that read:
                [e0, e1), [e1, e2), ..., [e_last, inf) from --win-len-edges
    win_long    int64 (4 ends, 2 repeat, 2 shared, NCOV, --window-long)
                reads with length >= --min-len-long, d < --window-long
    folded      int64 (2 reads fwd/rev, 2 repeat, 2 shared, NCOV, --fold-max + 1)
                d = min(d5, d3); last bin holds d >= --fold-max
    folded_len  int64 (2 reads, 2 repeat, 2 shared, NCOV, NLEN, NFD)
                by read length in --len-bin bins (last bin >= --len-max) and
                d = min(d5, d3) in --fd-bin bins (last bin >= --fold-max)
    win_len_edges, window, min_len_long, window_long, fold_max, cov_max,
    len_bin, len_max, fd_bin   the counting parameters
    n_informative, n_pairs, n_pairs_len_dropped, n_pairs_sa_dropped,
    n_pairs_nonov_dropped      totals (merge: summed)
    dcall_*     every variant (chrom, pos0, alt, fwd/rev qpos and read length)
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import pysam  # type: ignore

PIPELINE_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PIPELINE_SCRIPTS)
from variant_calling import get_ssdna_info  # noqa: E402

ENDS = ["fwd_5", "fwd_3", "rev_5", "rev_3"]
FLUSH = 5_000_000


def bed_to_bool(path, chrom, length):
    """Boolean array over the chromosome, True inside any BED interval."""
    out = np.zeros(length, dtype=bool)
    if not path:
        return out
    d = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                    names=["chrom", "start", "end"], dtype={"chrom": str})
    d = d[d.chrom == chrom]
    for s, e in zip(d.start.to_numpy(), d.end.to_numpy()):
        out[s:e] = True
    return out


def cntx_to_cov(path, length, cov_max):
    """Coverage per site from a cntx.bed.gz, capped at cov_max."""
    cov = np.zeros(length, dtype=np.uint8)
    for ch in pd.read_csv(path, sep="\t", header=None, usecols=[1, 4],
                          names=["start", "inf"], dtype={"start": np.int64, "inf": np.int64},
                          chunksize=20_000_000):
        cov[ch.start.to_numpy()] = np.minimum(ch.inf.to_numpy(), cov_max)
    return cov


class Hist:
    """Flat histogram filled from buffered index arrays."""

    def __init__(self, shape):
        self.shape = shape
        self.size = int(np.prod(shape))
        self.counts = np.zeros(self.size, dtype=np.int64)
        self.buf, self.n = [], 0

    def add(self, idx):
        if len(idx):
            self.buf.append(idx)
            self.n += len(idx)
            if self.n >= FLUSH:
                self.flush()

    def flush(self):
        if self.buf:
            self.counts += np.bincount(np.concatenate(self.buf), minlength=self.size)
            self.buf, self.n = [], 0

    def array(self):
        self.flush()
        return self.counts.reshape(self.shape)


def count(args):
    with pysam.FastaFile(args.ref) as fa:
        seq = fa.fetch(args.chrom).upper()
    genome = {args.chrom: seq}
    clen = len(seq)

    ncov = args.cov_max + 1
    rep = bed_to_bool(args.repeat_bed, args.chrom, clen)
    sh = bed_to_bool(args.shared_bed, args.chrom, clen)
    cov = cntx_to_cov(args.cntx, clen, args.cov_max)
    strat = ((rep.astype(np.uint16) * 2 + sh.astype(np.uint16)) * ncov
             + cov.astype(np.uint16))
    del rep, sh, cov
    nstrat = 4 * ncov
    sys.stderr.write(f"[edge_effect_counts] {args.chrom}: {clen:,} bp, strata ready\n")

    W, WL, FM = args.window, args.window_long, args.fold_max
    edges = np.asarray(args.win_len_edges, dtype=np.int64)
    ncls = len(edges)
    LB, NLEN = args.len_bin, args.len_max // args.len_bin + 1
    FB, NFD = args.fd_bin, args.fold_max // args.fd_bin + 1
    win = Hist((4, ncls, nstrat, W))
    folded_len = Hist((2, nstrat, NLEN, NFD))
    win_long = Hist((4, nstrat, WL))
    folded = Hist((2, nstrat, FM + 1))

    n_inf = n_pairs = n_len_drop = n_sa_drop = n_nonov_drop = 0
    dcalls = []
    pairs = {}

    bam = pysam.AlignmentFile(sys.stdin if args.input == "-" else args.input,
                              "r" if args.input == "-" else "rb")
    for i, aln in enumerate(bam.fetch(until_eof=True)):
        if i % 100_000 == 0 and i > 0:
            sys.stderr.write(f"{i} reads processed...\n")
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        elems = aln.query_name.split("/")
        key, strand = elems[1], elems[-1]
        if key not in pairs:
            pairs[key] = {strand: aln}
            continue
        pairs[key][strand] = aln
        if len(pairs[key]) != 2:
            continue
        fwd, rev = pairs[key].get("fwd"), pairs[key].get("rev")
        del pairs[key]
        if fwd is None or rev is None:
            continue
        # Molecule filters: strand length difference, chimeric molecules,
        # strand alignment non-overlap
        lf, lr = fwd.query_length, rev.query_length
        if args.max_len_div and abs(lf - lr) > args.max_len_div:
            n_len_drop += 1
            continue
        if args.drop_supplementary and (fwd.has_tag("SA") or rev.has_tag("SA")):
            n_sa_drop += 1
            continue
        if args.max_strand_nonoverlap:
            fs, fe = fwd.reference_start, fwd.reference_end
            rs, re_ = rev.reference_start, rev.reference_end
            nonov = (max(fe, re_) - min(fs, rs)) - max(0, min(fe, re_) - max(fs, rs))
            if nonov > args.max_strand_nonoverlap:
                n_nonov_drop += 1
                continue
        n_pairs += 1

        _, fmm, _ = get_ssdna_info(fwd, genome, args.minqual)
        _, rmm, _ = get_ssdna_info(rev, genome, args.minqual)
        pos, fq, rq = [], [], []
        for p in fmm.keys() & rmm.keys():
            f, r = fmm[p], rmm[p]
            if f[1] == r[1]:
                pos.append(p)
                fq.append(f[5])
                rq.append(r[5])
                if f[1] != f[0]:
                    dcalls.append((p, f[1], f[5], r[5], lf, lr))
        if not pos:
            continue
        n_inf += len(pos)

        P = np.asarray(pos, dtype=np.int64)
        s = strat[P].astype(np.int64)
        for frame, (q, L) in enumerate(((np.asarray(fq, dtype=np.int64), lf),
                                        (np.asarray(rq, dtype=np.int64), lr))):
            d5, d3 = q, L - 1 - q
            fd = np.minimum(np.minimum(d5, d3), FM)
            folded.add((frame * nstrat + s) * (FM + 1) + fd)
            lb = min(L // LB, NLEN - 1)
            folded_len.add(((frame * nstrat + s) * NLEN + lb) * NFD + np.minimum(fd // FB, NFD - 1))
            cls = int(np.searchsorted(edges, L, side="right")) - 1
            for k, d in ((2 * frame, d5), (2 * frame + 1, d3)):
                if cls >= 0:
                    m = d < W
                    win.add(((k * ncls + cls) * nstrat + s[m]) * W + d[m])
                if L >= args.min_len_long:
                    m = d < WL
                    win_long.add((k * nstrat + s[m]) * WL + d[m])
    bam.close()

    sh5 = (4, 2, 2, ncov)
    dc = np.array(dcalls, dtype=object).reshape(-1, 6)
    np.savez_compressed(
        args.output,
        ends=np.array(ENDS),
        win=win.array().reshape((4, ncls, 2, 2, ncov, W)),
        win_long=win_long.array().reshape(sh5 + (WL,)),
        folded=folded.array().reshape((2, 2, 2, ncov, FM + 1)),
        folded_len=folded_len.array().reshape((2, 2, 2, ncov, NLEN, NFD)),
        n_informative=np.int64(n_inf),
        n_pairs=np.int64(n_pairs),
        n_pairs_len_dropped=np.int64(n_len_drop),
        n_pairs_sa_dropped=np.int64(n_sa_drop),
        n_pairs_nonov_dropped=np.int64(n_nonov_drop),
        win_len_edges=edges, window=np.int64(W), min_len_long=np.int64(args.min_len_long),
        window_long=np.int64(WL), fold_max=np.int64(FM), cov_max=np.int64(args.cov_max),
        len_bin=np.int64(LB), len_max=np.int64(args.len_max), fd_bin=np.int64(FB),
        dcall_chrom=np.array([args.chrom] * len(dc)),
        dcall_pos0=dc[:, 0].astype(np.int64),
        dcall_alt=dc[:, 1].astype(str),
        dcall_fwd_qpos=dc[:, 2].astype(np.int64),
        dcall_rev_qpos=dc[:, 3].astype(np.int64),
        dcall_fwd_len=dc[:, 4].astype(np.int64),
        dcall_rev_len=dc[:, 5].astype(np.int64),
    )
    sys.stderr.write(f"Seen {n_inf} informative bases, {n_pairs} pairs, "
                     f"{len(dc)} d-calls; dropped {n_len_drop} pairs by length divergence, "
                     f"{n_sa_drop} by supplementary alignment, {n_nonov_drop} by strand non-overlap\n")


def merge(args):
    sums = ["win", "win_long", "folded", "folded_len",
            "n_informative", "n_pairs", "n_pairs_len_dropped",
            "n_pairs_sa_dropped", "n_pairs_nonov_dropped"]
    cats = ["dcall_chrom", "dcall_pos0", "dcall_alt", "dcall_fwd_qpos",
            "dcall_rev_qpos", "dcall_fwd_len", "dcall_rev_len"]
    out, fixed = {}, None
    for f in args.inputs:
        z = np.load(f)
        if fixed is None:
            fixed = {k: z[k] for k in z.files if k not in sums + cats}
            out = {k: z[k].astype(np.int64) for k in sums}
            out.update({k: [z[k]] for k in cats})
        else:
            for k, v in fixed.items():
                if not np.array_equal(v, z[k]):
                    sys.exit(f"[merge] {k} differs in {f}")
            for k in sums:
                out[k] = out[k] + z[k]
            for k in cats:
                out[k].append(z[k])
    for k in cats:
        out[k] = np.concatenate(out[k])
    np.savez_compressed(args.output, **fixed,
                        chroms=np.array([os.path.basename(f).split(".")[0] for f in args.inputs]),
                        **out)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("count", help="count one chromosome of one library")
    c.add_argument("-i", "--input", required=True, help="SAM/BAM from read_and_base_filters.py, '-' for stdin")
    c.add_argument("--chrom", required=True)
    c.add_argument("--ref", required=True, help="reference FASTA (indexed)")
    c.add_argument("--cntx", required=True, help="cntx.bed.gz of variant_calling.py for this chromosome")
    c.add_argument("--repeat-bed", default=None)
    c.add_argument("--shared-bed", default=None)
    c.add_argument("-q", "--minqual", type=int, default=40)
    c.add_argument("--max-len-div", type=int, default=0,
                   help="strand length difference filter: maximum |fwd_len - rev_len| (0 = off)")
    c.add_argument("--drop-supplementary", action="store_true", default=False,
                   help="chimeric molecule filter: drop molecules with an SA tag on either read")
    c.add_argument("--max-strand-nonoverlap", type=int, default=0,
                   help="strand alignment non-overlap filter: maximum sites covered by "
                        "only one of the two alignments (0 = off)")
    c.add_argument("--win-len-edges", type=int, nargs="+", default=[3000, 4000, 5000],
                   help="lower edges of the read-length classes of `win`")
    c.add_argument("--window", type=int, default=1500)
    c.add_argument("--min-len-long", type=int, default=5000)
    c.add_argument("--window-long", type=int, default=2500)
    c.add_argument("--fold-max", type=int, default=3000)
    c.add_argument("--cov-max", type=int, default=40)
    c.add_argument("--len-bin", type=int, default=250)
    c.add_argument("--len-max", type=int, default=16000)
    c.add_argument("--fd-bin", type=int, default=100)
    c.add_argument("-o", "--output", required=True)

    m = sub.add_parser("merge", help="sum per-chromosome files")
    m.add_argument("-o", "--output", required=True)
    m.add_argument("inputs", nargs="+")

    args = p.parse_args()
    count(args) if args.cmd == "count" else merge(args)


if __name__ == "__main__":
    main()
