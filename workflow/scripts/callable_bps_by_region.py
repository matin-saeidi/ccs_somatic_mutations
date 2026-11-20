#!/usr/bin/env python3
"""
callable_bps_by_region.py

Callable bps split by genomic region (accessible and inaccessible chromatin;
Methods, Chromatin accessibility): callable_sites and callable_depth per
(region x trinucleotide x minimum coverage), per library.

A site is callable at minimum coverage F when
    F <= selecting coverage <= its 99th percentile  AND  it is not a masked site
    AND its coverage in --cntx-dir is at least 1
The selecting coverage is the site's coverage in --cntx-dir, or, with
--untrimmed-cntx-dir, its coverage in the untrimmed reads. callable_depth sums
the coverage in --cntx-dir. With --cntx-dir ext450 and --untrimmed-cntx-dir
ext0 these are the callable bps of callable_bps.py; --check-against stops the
run unless the totals over regions and trinucleotides equal that table at
every minimum coverage.

Regions: each --scheme NAME=BED names a BED whose 4th column is a region label
(0-based half-open, intervals non-overlapping). Sites outside every interval
of a scheme get the label "none". Several schemes are cross-tabulated.

Trinucleotide: the cntx triplet (column 6), folded to the 32 pyrimidine-centred
contexts; triplets with a non-ACGT base go to "NNN".

Output: --out, one row per sample x floor x region x triplet, zeros included:
    sample floor p99 <scheme ...> triplet callable_sites callable_depth
    (p99 is named untrimmed_p99 with --untrimmed-cntx-dir)

Usage
    python callable_bps_by_region.py --cntx-dir <results>/variant_calling/ext450 \\
        --untrimmed-cntx-dir <results>/variant_calling/ext0 \\
        --mask excluded_regions.bed --fai b6_t2t_mm39XYM.fa.fai \\
        --samples <library> --chroms T2T1 ... T2T19 \\
        --scheme accessibility=accessibility_regions.bed.gz \\
        --floors 5 6 7 8 9 10 11 12 --ext 450 --qual 40 \\
        --check-against <results>/callable_bps/trimmed/denominator.tsv --out <library>.tsv
"""
import argparse
import itertools
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from callable_bps_untrimmed import (
    MAX_COV_PCT, all_triplets, cntx_chunks, fold_triplet, in_mask, load_mask,
    weighted_percentile,
)

CTX = all_triplets(fold=True) + ["NNN"]     # 32 folded contexts + non-ACGT bin
CTX_IDX = {c: i for i, c in enumerate(CTX)}


def read_scheme(spec):
    """'name=path.bed[.gz]' -> (name, labels, {chrom: DataFrame(start, end, code)})."""
    name, path = spec.split("=", 1)
    d = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2, 3],
                    names=["chrom", "start", "end", "label"],
                    dtype={"chrom": str, "start": np.int64, "end": np.int64, "label": str})
    labels = ["none"] + sorted(d.label.unique())
    code = {l: i for i, l in enumerate(labels)}
    d["code"] = d.label.map(code).astype(np.uint8)
    return name, labels, {c: g[["start", "end", "code"]] for c, g in d.groupby("chrom")}


def label_array(iv, L):
    """Per-base label codes for one contig (0 = none). Overlaps are an error."""
    lab = np.zeros(L, dtype=np.uint8)
    if iv is None:
        return lab
    s, e, v = iv.start.to_numpy(), iv.end.to_numpy(), iv.code.to_numpy()
    o = np.argsort(s)
    s, e, v = s[o], e[o], v[o]
    if (s[1:] < e[:-1]).any():
        sys.exit("[error] a region BED has overlapping intervals")
    for a, b, c in zip(s, e, v):
        lab[a:b] = c
    return lab


def dense_coverage(path, L):
    """informative_cov of one contig's cntx file as a per-base array (0 where absent)."""
    cov = np.zeros(L, dtype=np.int16)
    for ch in cntx_chunks(path, usecols=("start", "informative_cov")):
        cov[ch["start"].to_numpy()] = np.minimum(ch["informative_cov"].to_numpy(), 32767)
    return cov


def process_sample(sample, a, mask, schemes, clen, p99):
    """Callable sites and depth by (region combination, context, selecting coverage 0..p99)."""
    fmin = min(a.floors)
    ncv = int(p99) + 1                                   # selecting-coverage bins 0..p99
    sizes = [len(lbl) for _, lbl, _ in schemes]
    stride = np.cumprod([1] + sizes[::-1])[:-1][::-1]    # row-major combo index
    ncombo = int(np.prod(sizes))
    n_bins = ncombo * len(CTX) * ncv
    sites_h = np.zeros(n_bins, dtype=np.int64)
    depth_h = np.zeros(n_bins, dtype=np.int64)
    n_seen = 0
    for chrom in a.chroms:
        path = os.path.join(a.cntx_dir, sample, f"{chrom}.ext{a.ext}.qual{a.qual}.cntx.bed.gz")
        if not os.path.exists(path):
            sys.exit(f"[{sample}] missing cntx: {path}")
        sel_cov = None
        if a.untrimmed_cntx_dir:
            upath = os.path.join(a.untrimmed_cntx_dir, sample, f"{chrom}.ext0.qual{a.qual}.cntx.bed.gz")
            if not os.path.exists(upath):
                sys.exit(f"[{sample}] missing cntx: {upath}")
            sel_cov = dense_coverage(upath, clen[chrom])
        labs = [label_array(ivs.get(chrom), clen[chrom]) for _, _, ivs in schemes]
        mc = mask.get(chrom)
        for ch in cntx_chunks(path):
            n_seen += len(ch)
            ic = ch["informative_cov"].to_numpy()
            st = ch["start"].to_numpy()
            sel = ic if sel_cov is None else sel_cov[st].astype(np.int32)
            keep = (sel >= fmin) & (sel <= p99) & (ic >= 1)
            if not keep.any():
                continue
            keep &= ~in_mask(mc, st)
            if not keep.any():
                continue
            st, icf, self_ = st[keep], ic[keep].astype(np.int64), sel[keep].astype(np.int64)
            tri = ch["triplet"]
            lut = np.array([CTX_IDX.get(fold_triplet(t) or "NNN", CTX_IDX["NNN"])
                            for t in tri.cat.categories], dtype=np.int64)
            codes = tri.cat.codes.to_numpy()[keep]
            ctx = np.where(codes >= 0, lut[np.maximum(codes, 0)], CTX_IDX["NNN"])
            combo = np.zeros(len(st), dtype=np.int64)
            for lab, sd in zip(labs, stride):
                combo += lab[st].astype(np.int64) * sd
            flat = (combo * len(CTX) + ctx) * ncv + self_
            sites_h += np.bincount(flat, minlength=n_bins)
            depth_h += np.rint(np.bincount(flat, weights=icf, minlength=n_bins)).astype(np.int64)
        sys.stderr.write(f"  [{sample}] {chrom} done ({n_seen:,} positions scanned)\n")
    shape = (ncombo, len(CTX), ncv)
    return sites_h.reshape(shape), depth_h.reshape(shape), sizes


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cntx-dir", required=True, help="cntx files whose informative_cov is counted")
    p.add_argument("--untrimmed-cntx-dir", default=None,
                   help="choose sites by the untrimmed informative_cov (ext0 cntx files) and its p99")
    p.add_argument("--mask", required=True)
    p.add_argument("--fai", required=True, help="reference .fai (contig lengths)")
    p.add_argument("--samples", nargs="+", required=True)
    p.add_argument("--chroms", nargs="+", required=True,
                   help="contigs to read, e.g. T2T1 .. T2T19")
    p.add_argument("--scheme", action="append", required=True, metavar="NAME=BED")
    p.add_argument("--ext", type=int, default=450)
    p.add_argument("--qual", type=int, default=40)
    p.add_argument("--floors", type=int, nargs="+", default=[5, 6, 7, 8, 9, 10, 11, 12])
    p.add_argument("--check-against", default=None, help="denominator.tsv to reproduce")
    p.add_argument("--out", required=True)
    a = p.parse_args()

    fai = pd.read_csv(a.fai, sep="\t", header=None, usecols=[0, 1], names=["n", "L"])
    clen = dict(zip(fai.n, fai.L))
    mask = load_mask(a.mask, a.chroms)
    schemes = [read_scheme(s) for s in a.scheme]
    for name, labels, _ in schemes:
        sys.stderr.write(f"scheme {name}: {labels}\n")
    ref = pd.read_csv(a.check_against, sep="\t") if a.check_against else None
    p99_col = "untrimmed_p99" if a.untrimmed_cntx_dir else "p99"

    out = []
    for s in a.samples:
        t0 = time.time()
        sel_dir = a.untrimmed_cntx_dir or a.cntx_dir
        dist = pd.read_csv(os.path.join(sel_dir, s, "informative_coverage_dist.tsv"), sep="\t")
        p99 = weighted_percentile(dist["coverage"], dist["n_sites"], MAX_COV_PCT)
        sys.stderr.write(f"[{s}] {p99_col}={p99:.0f} scanning ...\n")
        sites_h, depth_h, sizes = process_sample(s, a, mask, schemes, clen, p99)
        combos = list(itertools.product(*[lbl for _, lbl, _ in schemes]))
        for F in a.floors:
            sites = sites_h[:, :, F:].sum(2)
            depth = depth_h[:, :, F:].sum(2)
            if ref is not None:
                r = ref[(ref["sample"] == s) & (ref.floor == F)]
                if len(r) != 1:
                    sys.exit(f"[error] {s} floor {F}: no row in {a.check_against}")
                if (int(sites.sum()) != int(r.callable_sites.iloc[0])
                        or int(depth.sum()) != int(r.callable_depth.iloc[0])):
                    sys.exit(f"[error] {s} floor {F}: totals {int(sites.sum()):,} / "
                             f"{int(depth.sum()):,} != denominator.tsv "
                             f"{int(r.callable_sites.iloc[0]):,} / "
                             f"{int(r.callable_depth.iloc[0]):,}")
            for ci, combo in enumerate(combos):
                for k, ctx in enumerate(CTX):
                    row = {"sample": s, "floor": F, p99_col: p99}
                    row.update({name: lab for (name, _, _), lab in zip(schemes, combo)})
                    row.update(triplet=ctx, callable_sites=int(sites[ci, k]),
                               callable_depth=int(depth[ci, k]))
                    out.append(row)
        msg = " (totals match denominator.tsv at every floor)" if ref is not None else ""
        sys.stderr.write(f"[{s}] done in {time.time()-t0:.0f}s{msg}\n")
    pd.DataFrame(out).to_csv(a.out, sep="\t", index=False)
    sys.stderr.write(f"-> {a.out}\n")


if __name__ == "__main__":
    main()
