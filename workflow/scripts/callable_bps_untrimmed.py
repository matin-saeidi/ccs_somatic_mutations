#!/usr/bin/env python3
"""
callable_bps_untrimmed.py

Callable bps of the untrimmed reads (Methods, Coverage filter; Estimation of
the somatic mutation rate), per library and minimum coverage F. A site is
callable when

    F <= coverage <= the library's 99th percentile of coverage
    AND it is not a masked site

    callable_depth(F) = coverage summed over the callable sites  (callable bps)
    callable_sites(F) = number of callable sites

Coverage is informative_cov of the cntx files of variant_calling.py; the 99th
percentile is that of the library's informative_coverage_dist.tsv. Also holds
the functions that callable_bps.py, callable_bps_by_region.py and
mutation_tables.py import (masked-site lookup, 99th percentile, trinucleotide
folding, cntx reader).

Output: --out, columns sample floor p99 callable_depth callable_sites
        (--fai adds autosome_bp and pct_autosome_callable)

Usage:
    python callable_bps_untrimmed.py --cntx-dir <results>/variant_calling/ext0 \\
        --mask <results>/masked_sites/excluded_regions.bed \\
        --samples <library> --chroms T2T1 ... T2T19 --ext 0 --qual 40 \\
        --floors 5 6 7 8 9 10 --out <library>.tsv
"""
import argparse
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

CNTX_COLS = ["chrom", "start", "end", "total_cov", "informative_cov", "triplet"]
CNTX_DTYPES = {"start": np.int64, "informative_cov": np.int32, "triplet": "category"}
MAX_COV_PCT, CHUNK = 99, 5_000_000
BASES = "ACGT"
_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def fold_triplet(t):
    """Reference triplet -> its pyrimidine-centred form (TGT -> ACA, GGC -> GCC).
    None if the triplet is not pure ACGT."""
    if len(t) != 3 or any(c not in BASES for c in t):
        return None
    return t if t[1] in "CT" else revcomp(t)


def all_triplets(fold=True):
    """The 32 pyrimidine-centred contexts (fold=True) or all 64, sorted."""
    raw = [a + b + c for a in BASES for b in BASES for c in BASES]
    return sorted({fold_triplet(t) for t in raw} if fold else set(raw))


def cntx_chunks(path, usecols=("start", "informative_cov", "triplet"), chunksize=CHUNK):
    """Stream the given columns of a cntx.bed.gz in chunks; the triplet column is
    read as a category. Uses pigz when it is installed, else Python's gzip."""
    kw = dict(sep="\t", header=None, names=CNTX_COLS, usecols=list(usecols),
              dtype={c: CNTX_DTYPES[c] for c in usecols}, chunksize=chunksize)
    pigz = shutil.which("pigz")
    if pigz:
        proc = subprocess.Popen([pigz, "-dc", path], stdout=subprocess.PIPE, bufsize=1 << 24)
        yield from pd.read_csv(proc.stdout, **kw)
        proc.stdout.close()
        if proc.wait() != 0:
            sys.exit(f"[error] pigz failed on {path}")
    else:
        yield from pd.read_csv(path, compression="gzip", **kw)


def weighted_percentile(coverage, site_counts, pct):
    """pct-th percentile of coverage over sites, from a coverage histogram."""
    coverage = np.asarray(coverage, dtype=float)
    site_counts = np.asarray(site_counts, dtype=float)
    order = np.argsort(coverage)
    coverage, site_counts = coverage[order], site_counts[order]
    cumfrac = np.cumsum(site_counts) / site_counts.sum()
    idx = min(np.searchsorted(cumfrac, pct / 100.0), len(coverage) - 1)
    return float(coverage[idx])


def load_mask(path, chroms):
    """BED -> {chrom: (starts, ends)} as sorted, merged numpy arrays."""
    df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                     names=["chrom", "start", "end"], dtype={"chrom": str})
    df = df[df["chrom"].isin(chroms)]
    out = {}
    for c, g in df.groupby("chrom"):
        g = g.sort_values("start")
        s, e = g["start"].to_numpy(np.int64), g["end"].to_numpy(np.int64)
        ms, me = [s[0]], [e[0]]
        for i in range(1, len(s)):
            if s[i] <= me[-1]:
                me[-1] = max(me[-1], e[i])
            else:
                ms.append(s[i]); me.append(e[i])
        out[c] = (np.array(ms, np.int64), np.array(me, np.int64))
    return out


def in_mask(mask_c, starts):
    """True where a 0-based position is inside one of the sorted, non-overlapping
    [s, e) intervals."""
    if mask_c is None:
        return np.zeros(len(starts), dtype=bool)
    ms, me = mask_c
    i = np.searchsorted(ms, starts, side="right") - 1
    ok = i >= 0
    out = np.zeros(len(starts), dtype=bool)
    if ok.any():
        j = i[ok]
        out[ok] = starts[ok] < me[j]
    return out


def add_bincount(hist, values):
    if len(values) == 0:
        return hist
    b = np.bincount(values.astype(np.int64))
    if len(b) > len(hist):
        hist = np.concatenate([hist, np.zeros(len(b) - len(hist), dtype=np.int64)])
    hist[:len(b)] += b
    return hist


def band_sum(hist, lo, hi):
    hi = min(int(hi), len(hist) - 1)
    return 0 if hi < lo else int(hist[lo:hi + 1].sum())


def band_sum_depth(hist, lo, hi):
    hi = min(int(hi), len(hist) - 1)
    if hi < lo:
        return 0
    c = np.arange(lo, hi + 1)
    return int((c * hist[lo:hi + 1]).sum())


def process_sample(sample, cntx_dir, mask, ext, qual, chroms):
    all_h = np.zeros(1, dtype=np.int64)
    msk_h = np.zeros(1, dtype=np.int64)
    for chrom in chroms:
        path = os.path.join(cntx_dir, sample, f"{chrom}.ext{ext}.qual{qual}.cntx.bed.gz")
        if not os.path.exists(path):
            sys.exit(f"[{sample}] missing cntx: {path}")
        mc = mask.get(chrom)
        for chunk in pd.read_csv(path, sep="\t", header=None, names=CNTX_COLS,
                                 usecols=["start", "informative_cov"],
                                 dtype={"start": np.int64, "informative_cov": np.int64},
                                 compression="gzip", chunksize=CHUNK):
            ic = chunk["informative_cov"].to_numpy()
            all_h = add_bincount(all_h, ic)
            m = in_mask(mc, chunk["start"].to_numpy())
            if m.any():
                msk_h = add_bincount(msk_h, ic[m])
    return all_h, msk_h


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cntx-dir", required=True)
    p.add_argument("--mask", required=True, help="masked sites BED (0-based half-open)")
    p.add_argument("--out", required=True)
    p.add_argument("--ext", type=int, default=0)
    p.add_argument("--qual", type=int, default=40)
    p.add_argument("--floors", type=int, nargs="+", default=[5, 6, 7, 8, 9, 10])
    p.add_argument("--samples", nargs="+", required=True)
    p.add_argument("--chroms", nargs="+", required=True,
                   help="contigs to read, e.g. T2T1 .. T2T19")
    p.add_argument("--fai", default=None,
                   help="reference .fai; adds autosome_bp and pct_autosome_callable "
                        "(callable_sites / summed length of --chroms)")
    args = p.parse_args()

    autosome_bp = None
    if args.fai:
        want = set(args.chroms)
        autosome_bp = 0
        with open(args.fai) as fh:
            for line in fh:
                f = line.split("\t")
                if f[0] in want:
                    autosome_bp += int(f[1])
        if autosome_bp == 0:
            sys.exit(f"[error] --fai {args.fai} matched none of --chroms")
        sys.stderr.write(f"callable-footprint denominator: {autosome_bp:,} bp "
                         f"over {len(want)} contigs\n")

    mask = load_mask(args.mask, args.chroms)
    n_iv = sum(len(v[0]) for v in mask.values())
    n_bp = sum(int((v[1] - v[0]).sum()) for v in mask.values())
    sys.stderr.write(f"mask: {n_iv:,} merged intervals, {n_bp:,} bp\n")

    rows = []
    for s in args.samples:
        sys.stderr.write(f"[{s}] scanning cntx ...\n")
        all_h, msk_h = process_sample(s, args.cntx_dir, mask, args.ext,
                                      args.qual, args.chroms)

        dist = pd.read_csv(os.path.join(args.cntx_dir, s,
                                        "informative_coverage_dist.tsv"), sep="\t")
        p99 = weighted_percentile(dist["coverage"], dist["n_sites"], MAX_COV_PCT)

        scan_total, dist_total = int(all_h.sum()), int(dist["n_sites"].sum())
        if scan_total != dist_total:
            sys.stderr.write(f"  [WARN] scan sites {scan_total:,} != dist sites "
                             f"{dist_total:,} -> cntx and "
                             f"informative_coverage_dist.tsv out of sync\n")

        for F in args.floors:
            sites = band_sum(all_h, F, p99) - band_sum(msk_h, F, p99)
            depth = band_sum_depth(all_h, F, p99) - band_sum_depth(msk_h, F, p99)
            row = dict(sample=s, floor=F, p99=p99,
                       callable_depth=depth, callable_sites=sites)
            if autosome_bp:
                row["autosome_bp"] = autosome_bp
                row["pct_autosome_callable"] = round(100.0 * sites / autosome_bp, 4)
            rows.append(row)
        sys.stderr.write(f"  p99={p99:.0f}  covered_positions={scan_total:,}  "
                         f"masked_positions={int(msk_h.sum()):,}\n")
        if autosome_bp:
            for r in rows:
                if r["sample"] == s:
                    sys.stderr.write(f"    floor {r['floor']:>2}: callable_sites="
                                     f"{r['callable_sites']:>13,}  "
                                     f"{r['pct_autosome_callable']:>6.2f}% of autosomes\n")

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    wide = out.pivot(index="sample", columns="floor", values="callable_depth")
    wide.columns = [f"depth_floor{c}" for c in wide.columns]
    print("\ncallable_depth (SUM informative_cov over kept sites):")
    print(pd.concat([out.groupby("sample")["p99"].first().rename("p99"), wide],
                    axis=1).to_string())
    print(f"\nfull table -> {args.out}")


if __name__ == "__main__":
    main()
