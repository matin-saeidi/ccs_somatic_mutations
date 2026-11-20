#!/usr/bin/env python3
"""
strand_mismatch_counts.py

Opportunities and mismatches between the two strand consensus reads of each
molecule (Methods, Estimation of the sequencing error rate; Supplementary
Table 5, Supplementary Figs. 5, 6 and 14), by the number of passes (np) and the
sm and sx of each read.

A compared base is a site of a molecule where both strand consensus reads pass
the read, molecule and base filters (total_cov of the cntx files of
variant_calling.py). A mismatch is a compared base where exactly one of the two
reads carries the reference base (types f and r). The filters, as in
read_and_base_filters.py and variant_calling.py:
    reads       primary, MAPQ >= --mapq; length (with hard clips) in
                [--lmin, --lmax]; np >= --np-min; divergence
                1 - (=) / (= + X + I + D) <= --div-max
    bases       base quality >= --minqual; more than --pad bases from an
                insertion or deletion; more than --ext bases from a read end;
                sm >= --sm-min; sm / (sm + sx) > --sm-ratio; no N in the
                reference triplet
    molecules   both reads present; strand length difference <= --max-len-div;
                no SA tag (--drop-supplementary); strand alignment non-overlap
                <= --max-strand-nonoverlap

Every compared base is counted once for each read, at that read's np, sm and
sx, and by site:
    bit 0   site in --mask-bed (masked sites)
    bit 1   site in --mutation-tables (variants of the same animal)
Every mismatch is counted for the read carrying the non-reference base.

Every compared base is also counted once for each read by that read's np (bins
from --context-np-edges) and by its reference context as that read shows it,
5' to 3' (reads aligned to the minus strand reverse-complemented): the
trinucleotide centred on the base, and the lengths of three runs of identical
bases in the reference, capped at --run-max:
    own      the run containing the base
    before   the run ending just 5' of the base on the read; 0 if the base
             before it is part of its own run
    after    the run starting just 3' of the base on the read; 0 likewise
Mismatch records carry the same context for the read with the mismatch.

Subcommands
    count   one chromosome of one library:
                python strand_mismatch_counts.py count --bam <sample>.sorted.bam \\
                    --chrom T2T1 --ref b6_t2t_mm39XYM.fa --ext 450 \\
                    --lmin LMIN --lmax LMAX --max-len-div N --drop-supplementary \\
                    --max-strand-nonoverlap M --mask-bed excluded_regions.bed \\
                    --mutation-tables <library>/mutations_with_vaf.tsv [...] -o T2T1.npz
            --region START-END (0-based, end exclusive) reads only the
            alignments overlapping that part of the chromosome.
    merge   sum the chromosomes of one library:
                python strand_mismatch_counts.py merge -o <sample>.npz T2T1.npz T2T2.npz ...

Output arrays (npz)
    compared    int64 (2 reads fwd/rev, 4 site classes, NP, SM, SX)  opportunities
    mismatch    int64 (2 reads fwd/rev, 4 site classes, NP, SM, SX)  mismatches,
                by the read with the non-reference base
                np, sm and sx bins are 0 .. --np-max / --sm-max / --sx-max; the
                last bin holds everything at or above the cap
    np_max, sm_max, sx_max      the caps
    context     int64 (NP bins, 4 site classes, 64 trinucleotides, RUN, RUN, RUN)
                opportunities x reads by np bin and context as the read shows
                it: np bin k holds context_np_edges[k] <= np < context_np_edges[k + 1];
                trinucleotide index 16 * b(5') + 4 * b(centre) + b(3') with
                A, C, G, T = 0, 1, 2, 3, then own, before and after run lengths
                0 .. --run-max (last bin: at or above the cap)
    run_max     the run-length cap
    context_np_edges  lower edges of the np bins of `context`
    n_compared, n_pairs, n_pairs_len_dropped, n_pairs_sa_dropped,
    n_pairs_nonov_dropped       totals (merge: summed)
    mm_*        every mismatch: chrom, pos0, ref, triplet, call_type (f/r), err_base,
                err_reverse / other_reverse (read aligned to the reference
                minus strand), err_np/sm/sx, other_np/sm/sx, err_qpos and
                other_qpos (from the read's own 5' end), err_len, stratum, zmw,
                and the context of the read with the mismatch: read_trinuc
                (index as in `context`), own_run, before_run, after_run
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import pysam  # type: ignore

FLUSH = 5_000_000
MM_INT = ["pos0", "err_np", "err_sm", "err_sx", "other_np", "other_sm", "other_sx",
          "err_qpos", "other_qpos", "err_len", "stratum", "zmw",
          "read_trinuc", "own_run", "before_run", "after_run"]
MM_STR = ["chrom", "ref", "triplet", "call_type", "err_base"]
MM_BOOL = ["err_reverse", "other_reverse"]
SUMS = ["compared", "mismatch", "context", "n_compared", "n_pairs", "n_pairs_len_dropped",
        "n_pairs_sa_dropped", "n_pairs_nonov_dropped"]


def bed_to_bool(path, chrom, length):
    """Boolean array over the chromosome, True inside any BED interval."""
    out = np.zeros(length, dtype=bool)
    d = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                    names=["chrom", "start", "end"], dtype={"chrom": str})
    d = d[d.chrom == chrom]
    for s, e in zip(d.start.to_numpy(), d.end.to_numpy()):
        out[s:e] = True
    return out


def positions_to_bool(paths, chrom, length):
    """Boolean array over the chromosome, True at the positions of the variant tables."""
    out = np.zeros(length, dtype=bool)
    for p in paths:
        d = pd.read_csv(p, sep="\t", usecols=["chrom", "pos_1based"], dtype={"chrom": str})
        out[d.pos_1based[d.chrom == chrom].to_numpy() - 1] = True
    return out


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


def strand_view(a, args, n_triplet):
    """Bases of one strand consensus read that pass the read and base filters.

    Returns (rpos, base, sm, sx, qpos from the read's 5' end, np) with rpos
    sorted, or None if the whole read is filtered."""
    seq = a.query_sequence
    qlen = len(seq)
    npass = a.get_tag("np") if a.has_tag("np") else 0
    rlen = a.infer_read_length() or qlen
    if npass < args.np_min or (args.lmin and rlen < args.lmin) or (args.lmax and rlen > args.lmax):
        return None

    mask = np.zeros(qlen, dtype=bool)
    qi, ri = [], []
    qp, rp = 0, a.reference_start
    n_eq = n_x = n_i = n_d = 0
    for op, L in a.cigartuples:
        if op in (0, 7, 8):
            qi.append(np.arange(qp, qp + L)); ri.append(np.arange(rp, rp + L))
            if op == 7:
                n_eq += L
            elif op == 8:
                n_x += L
            qp += L; rp += L
        elif op == 1:
            mask[max(0, qp - args.pad):min(qlen, qp + L + args.pad)] = True
            n_i += L; qp += L
        elif op == 2:
            mask[max(0, qp - args.pad):min(qlen, qp + args.pad)] = True
            n_d += L; rp += L
        elif op == 3:
            rp += L
        elif op == 4:
            qp += L
    tot = n_eq + n_x + n_i + n_d
    if tot and 1.0 - n_eq / tot > args.div_max:
        return None
    if args.ext:
        mask[:args.ext] = True
        mask[max(qlen - args.ext, 0):] = True
    sm = np.asarray(a.get_tag("sm"), dtype=np.int64)
    sx = np.asarray(a.get_tag("sx"), dtype=np.int64)
    if a.is_reverse:
        sm, sx = sm[::-1], sx[::-1]
    smsx = sm + sx
    with np.errstate(divide="ignore", invalid="ignore"):
        mask |= (sm < args.sm_min) | (smsx == 0) | (sm / np.where(smsx == 0, 1, smsx) <= args.sm_ratio)
    q = np.asarray(a.query_qualities, dtype=np.int64)
    qi = np.concatenate(qi) if qi else np.zeros(0, np.int64)
    ri = np.concatenate(ri) if ri else np.zeros(0, np.int64)
    keep = ~mask[qi] & (q[qi] >= args.minqual) & ~n_triplet[ri]
    qi, ri = qi[keep], ri[keep]
    base = np.frombuffer(seq.encode(), dtype=np.uint8)[qi]
    nat = (qlen - 1 - qi) if a.is_reverse else qi
    return ri, base, sm[qi], sx[qi], nat, npass


def run_lengths(code, cap):
    """Length of the run of identical bases containing each position, capped at cap."""
    starts = np.r_[0, np.flatnonzero(code[1:] != code[:-1]) + 1]
    lens = np.diff(np.r_[starts, len(code)])
    return np.repeat(np.minimum(lens, cap).astype(np.uint8), lens)


def read_context(P, rev, code, rl, cap):
    """(trinucleotide index, own, before, after run lengths) of positions P as a
    read aligned to the minus (rev) or plus strand shows them."""
    l, c, r = code[P - 1].astype(np.int64), code[P].astype(np.int64), code[P + 1].astype(np.int64)
    own = rl[P].astype(np.int64)
    left = np.where(l != c, rl[P - 1], 0).astype(np.int64)
    right = np.where(r != c, rl[P + 1], 0).astype(np.int64)
    if rev:
        return 16 * (3 - r) + 4 * (3 - c) + (3 - l), own, right, left
    return 16 * l + 4 * c + r, own, left, right


def count(args):
    with pysam.FastaFile(args.ref) as fa:
        seq = fa.fetch(args.chrom).upper()
    clen = len(seq)
    refb = np.frombuffer(seq.encode(), dtype=np.uint8)
    code = np.full(clen + 1, 4, dtype=np.uint8)            # one N past the end
    for k, b in enumerate(b"ACGT"):
        code[:clen][refb == b] = k
    RUN = args.run_max + 1
    rl = run_lengths(code, args.run_max)
    isn = refb == ord("N")
    n_triplet = isn.copy()                       # N anywhere in the triplet centred on the position
    n_triplet[1:] |= isn[:-1]
    n_triplet[:-1] |= isn[1:]
    stratum = (bed_to_bool(args.mask_bed, args.chrom, clen).astype(np.int64)
               + 2 * positions_to_bool(args.mutation_tables, args.chrom, clen).astype(np.int64))
    sys.stderr.write(f"[strand_mismatch_counts] {args.chrom}: {clen:,} bp, strata ready\n")

    NP, SM, SX = args.np_max + 1, args.sm_max + 1, args.sx_max + 1
    shape = (2, 4, NP, SM, SX)
    compared, mismatch = Hist(shape), Hist(shape)
    np_edges = np.asarray(args.context_np_edges, dtype=np.int64)
    NPB = len(np_edges)
    context = Hist((NPB, 4, 64, RUN, RUN, RUN))

    def cell(strand, s, npv, sm, sx):
        return (((strand * 4 + s) * NP + min(npv, NP - 1)) * SM + np.minimum(sm, SM - 1)) * SX + np.minimum(sx, SX - 1)

    n_compared = n_pairs = n_len = n_sa = n_nonov = 0
    mms = {k: [] for k in MM_INT + MM_STR + MM_BOOL}
    pairs = {}
    bam = pysam.AlignmentFile(args.bam)
    region = dict(contig=args.chrom)
    if args.region:
        s0, e0 = (int(x) for x in args.region.split("-"))
        region.update(start=s0, stop=e0)
    for i, aln in enumerate(bam.fetch(**region)):
        if i % 100_000 == 0 and i > 0:
            sys.stderr.write(f"{i} reads processed...\n")
        if (aln.is_unmapped or aln.is_secondary or aln.is_supplementary
                or aln.mapping_quality < args.mapq):
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
        if args.max_len_div and abs(fwd.query_length - rev.query_length) > args.max_len_div:
            n_len += 1
            continue
        if args.drop_supplementary and (fwd.has_tag("SA") or rev.has_tag("SA")):
            n_sa += 1
            continue
        if args.max_strand_nonoverlap:
            fs, fe, rs, re_ = fwd.reference_start, fwd.reference_end, rev.reference_start, rev.reference_end
            nonov = (max(fe, re_) - min(fs, rs)) - max(0, min(fe, re_) - max(fs, rs))
            if nonov > args.max_strand_nonoverlap:
                n_nonov += 1
                continue
        n_pairs += 1

        F, R = strand_view(fwd, args, n_triplet), strand_view(rev, args, n_triplet)
        if F is None or R is None:
            continue
        P, jf, jr = np.intersect1d(F[0], R[0], assume_unique=True, return_indices=True)
        if not len(P):
            continue
        n_compared += len(P)
        s = stratum[P]
        fb, rb, ref = F[1][jf], R[1][jr], refb[P]
        compared.add(cell(0, s, F[5], F[2][jf], F[3][jf]))
        compared.add(cell(1, s, R[5], R[2][jr], R[3][jr]))
        ctx = {}
        for rd, rnp in ((fwd, F[5]), (rev, R[5])):
            tri, own, bef, aft = read_context(P, rd.is_reverse, code, rl, args.run_max)
            ctx[rd.is_reverse] = (tri, own, bef, aft)
            nb = max(int(np.searchsorted(np_edges, rnp, side="right")) - 1, 0)
            context.add(((((nb * 4 + s) * 64 + tri) * RUN + own) * RUN + bef) * RUN + aft)

        mis = (fb == ref) != (rb == ref)
        if not mis.any():
            continue
        for k in np.nonzero(mis)[0]:
            fe_ = fb[k] != ref[k]
            (E, je, ea), (O, jo, oa) = (((F, jf[k], fwd), (R, jr[k], rev)) if fe_ else
                                        ((R, jr[k], rev), (F, jf[k], fwd)))
            st = 0 if fe_ else 1
            mismatch.add(np.atleast_1d(cell(st, s[k], E[5], E[2][je], E[3][je])))
            p = int(P[k])
            tri, own, bef, aft = (v[k] for v in ctx[ea.is_reverse])
            for name, v in (("read_trinuc", tri), ("own_run", own), ("before_run", bef), ("after_run", aft),
                            ("pos0", p), ("err_np", E[5]), ("err_sm", E[2][je]), ("err_sx", E[3][je]),
                            ("other_np", O[5]), ("other_sm", O[2][jo]), ("other_sx", O[3][jo]),
                            ("err_qpos", E[4][je]), ("other_qpos", O[4][jo]), ("err_len", ea.query_length),
                            ("stratum", s[k]), ("zmw", int(key)),
                            ("chrom", args.chrom), ("ref", seq[p]), ("triplet", seq[p - 1:p + 2]),
                            ("call_type", "f" if fe_ else "r"), ("err_base", chr(E[1][je])),
                            ("err_reverse", ea.is_reverse), ("other_reverse", oa.is_reverse)):
                mms[name].append(v)
    bam.close()

    np.savez_compressed(
        args.output,
        compared=compared.array(), mismatch=mismatch.array(), context=context.array(),
        np_max=np.int64(args.np_max), sm_max=np.int64(args.sm_max), sx_max=np.int64(args.sx_max),
        run_max=np.int64(args.run_max), context_np_edges=np_edges,
        n_compared=np.int64(n_compared), n_pairs=np.int64(n_pairs),
        n_pairs_len_dropped=np.int64(n_len), n_pairs_sa_dropped=np.int64(n_sa),
        n_pairs_nonov_dropped=np.int64(n_nonov),
        **{f"mm_{k}": np.asarray(mms[k], dtype=np.int64) for k in MM_INT},
        **{f"mm_{k}": np.asarray(mms[k], dtype=str) for k in MM_STR},
        **{f"mm_{k}": np.asarray(mms[k], dtype=bool) for k in MM_BOOL},
    )
    sys.stderr.write(f"Seen {n_compared} opportunities, {n_pairs} pairs, {len(mms['pos0'])} mismatches; "
                     f"dropped pairs: {n_len} length divergence, {n_sa} supplementary, "
                     f"{n_nonov} strand non-overlap\n")


def merge(args):
    cats = [f"mm_{k}" for k in MM_INT + MM_STR + MM_BOOL]
    out, fixed = {}, None
    for f in args.inputs:
        z = np.load(f)
        if fixed is None:
            fixed = {k: z[k] for k in z.files if k not in SUMS + cats}
            out = {k: z[k].astype(np.int64) for k in SUMS}
            out.update({k: [z[k]] for k in cats})
        else:
            for k, v in fixed.items():
                if not np.array_equal(v, z[k]):
                    sys.exit(f"[merge] {k} differs in {f}")
            for k in SUMS:
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
    c.add_argument("--bam", required=True, help="aligned, indexed BAM of the library")
    c.add_argument("--chrom", required=True)
    c.add_argument("--region", default=None, help="START-END on --chrom (0-based, end exclusive)")
    c.add_argument("--ref", required=True, help="reference FASTA (indexed)")
    c.add_argument("--mask-bed", required=True, help="masked sites (site class bit 0)")
    c.add_argument("--mutation-tables", nargs="+", required=True,
                   help="mutations_with_vaf.tsv of the animal's libraries (site class bit 1)")
    c.add_argument("--ext", type=int, default=0, help="bases filtered at each read end")
    c.add_argument("--lmin", type=int, default=0, help="minimum read length (0 = off)")
    c.add_argument("--lmax", type=int, default=0, help="maximum read length (0 = off)")
    c.add_argument("--mapq", type=int, default=60)
    c.add_argument("--np-min", type=int, default=5)
    c.add_argument("--sm-min", type=int, default=5)
    c.add_argument("--sm-ratio", type=float, default=0.8, help="sm / (sm + sx) must be above this")
    c.add_argument("--pad", type=int, default=3, help="bases filtered on each side of an indel")
    c.add_argument("--div-max", type=float, default=0.01)
    c.add_argument("-q", "--minqual", type=int, default=40)
    c.add_argument("--max-len-div", type=int, default=0,
                   help="strand length difference filter: maximum |fwd_len - rev_len| (0 = off)")
    c.add_argument("--drop-supplementary", action="store_true",
                   help="chimeric molecule filter: drop molecules with an SA tag on either read")
    c.add_argument("--max-strand-nonoverlap", type=int, default=0,
                   help="strand alignment non-overlap filter: maximum sites covered by "
                        "only one of the two alignments (0 = off)")
    c.add_argument("--np-max", type=int, default=60)
    c.add_argument("--sm-max", type=int, default=80)
    c.add_argument("--sx-max", type=int, default=30)
    c.add_argument("--run-max", type=int, default=8)
    c.add_argument("--context-np-edges", type=int, nargs="+", default=[5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 35],
                   help="lower edges of the np bins of the context counts")
    c.add_argument("-o", "--output", required=True)

    m = sub.add_parser("merge", help="sum per-chromosome files")
    m.add_argument("-o", "--output", required=True)
    m.add_argument("inputs", nargs="+")

    args = p.parse_args()
    count(args) if args.cmd == "count" else merge(args)


if __name__ == "__main__":
    main()
