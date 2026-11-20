#!/usr/bin/env python3
"""
callable_bps.py

Callable bps of the trimmed reads (Methods, Coverage filter; Estimation of the
somatic mutation rate; Supplementary Tables 4 and 6). At minimum coverage F a
site is callable when

    F <= untrimmed coverage <= the library's 99th percentile of untrimmed coverage
    AND it is not a masked site
    AND it is covered by at least one molecule in the trimmed reads

    callable_depth = trimmed coverage summed over the callable sites (callable bps)
    callable_sites = number of callable sites

Coverage is informative_cov of the cntx files of variant_calling.py: untrimmed
from <cntx-dir>/ext0, trimmed from <cntx-dir>/ext450. The 99th percentile is
that of the untrimmed informative_coverage_dist.tsv.

Modes
    table   one library: every site of its untrimmed cntx files, counted by
            masked or not, trinucleotide, trimmed and untrimmed coverage; stops
            if the counts differ from either informative_coverage_dist.tsv
    derive  the tables of all libraries -> callable bps per library and minimum
            coverage, and callable bps per library and trinucleotide at
            --context-floor; stops if the untrimmed callable sites differ from
            --check-untrimmed (callable_bps_untrimmed.py)

Table <sample>.npz
    hist      int64 (2, 65, CAP+1, CAP+1): hist[m, k, t, u] = sites with mask
              bit m (1 = masked), trinucleotide CONTEXTS[k] (64 = non-ACGT),
              trimmed coverage t and untrimmed coverage u; coverage above CAP
              is counted at CAP
    contexts  the 64 trinucleotides, then "non-ACGT"
    cap       CAP

Output
    --denominator-out  sample floor untrimmed_p99 callable_depth callable_sites
                       autosome_bp pct_autosome_callable
    --context-out      sample triplet callable_bp callable_depth floor untrimmed_p99
                       (64 trinucleotides, + strand of the reference)

Usage
    python callable_bps.py table --sample S --chroms T2T1 ... T2T19 \\
        --cntx-dir <results>/variant_calling --mask excluded_regions.bed \\
        --fai b6_t2t_mm39XYM.fa.fai --out coverage_tables/S.npz
    python callable_bps.py derive --tables coverage_tables/*.npz \\
        --cntx-dir <results>/variant_calling --fai b6_t2t_mm39XYM.fa.fai \\
        --chroms T2T1 ... T2T19 --floors 5 6 7 8 9 10 11 12 --context-floor 10 \\
        --check-untrimmed <results>/callable_bps/untrimmed/denominator.tsv \\
        --denominator-out denominator.tsv --context-out context_opportunity.tsv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from callable_bps_untrimmed import (
    MAX_COV_PCT, all_triplets, cntx_chunks, in_mask, load_mask, weighted_percentile,
)

QUAL = 40
TRIM_EXT = 450
CAP = 64
NC = CAP + 1
CONTEXTS = all_triplets(fold=False) + ["non-ACGT"]
NCTX = len(CONTEXTS)
NON_ACGT = NCTX - 1
CTX_IDX = {c: i for i, c in enumerate(CONTEXTS[:-1])}
HIST_SHAPE = (2, NCTX, NC, NC)


def cntx_path(cntx_dir, ext, sample, chrom):
    path = os.path.join(cntx_dir, f"ext{ext}", sample, f"{chrom}.ext{ext}.qual{QUAL}.cntx.bed.gz")
    if not os.path.exists(path):
        sys.exit(f"[error] missing cntx: {path}")
    return path


def coverage_dist(cntx_dir, ext, sample):
    return pd.read_csv(os.path.join(cntx_dir, f"ext{ext}", sample, "informative_coverage_dist.tsv"), sep="\t")


def p99_of(dist):
    return weighted_percentile(dist["coverage"], dist["n_sites"], MAX_COV_PCT)


def contig_lengths(fai):
    d = pd.read_csv(fai, sep="\t", header=None, usecols=[0, 1], names=["name", "length"], dtype={"name": str})
    return dict(zip(d.name, d.length))


def chrom_hist(sample, chrom, cntx_dir, mc, L):
    """Sites of one contig's untrimmed cntx file by (masked, context, trimmed, untrimmed coverage)."""
    trimmed = np.zeros(L, np.int16)
    for ch in cntx_chunks(cntx_path(cntx_dir, TRIM_EXT, sample, chrom), usecols=("start", "informative_cov")):
        trimmed[ch["start"].to_numpy()] = np.minimum(ch["informative_cov"].to_numpy(), 32767)

    seen = np.zeros(L, bool)
    hist = np.zeros(int(np.prod(HIST_SHAPE)), np.int64)
    n_rows = n_above = 0
    for ch in cntx_chunks(cntx_path(cntx_dir, 0, sample, chrom)):
        st = ch["start"].to_numpy()
        u = ch["informative_cov"].to_numpy().astype(np.int64)
        t = trimmed[st].astype(np.int64)
        seen[st] = True
        n_rows += len(st)
        n_above += int((t > u).sum())
        m = in_mask(mc, st).astype(np.int64)
        tri = ch["triplet"]
        lut = np.array([CTX_IDX.get(c, NON_ACGT) for c in tri.cat.categories] + [NON_ACGT], np.int64)
        codes = tri.cat.codes.to_numpy().astype(np.int64)
        k = lut[np.where(codes >= 0, codes, len(lut) - 1)]
        flat = ((m * NCTX + k) * NC + np.minimum(t, CAP)) * NC + np.minimum(u, CAP)
        hist += np.bincount(flat, minlength=len(hist))
    n_missing = int((trimmed[~seen] > 0).sum())
    print(f"{sample} {chrom}: {n_rows:,} untrimmed cntx rows; trimmed > untrimmed coverage at {n_above:,}; "
          f"trimmed coverage at positions absent from the untrimmed cntx: {n_missing:,}", flush=True)
    if n_above or n_missing:
        sys.exit("[error] the trimmed and untrimmed cntx files do not pair position by position")
    return hist.reshape(HIST_SHAPE)


def run_table(a):
    lengths = contig_lengths(a.fai)
    mask = load_mask(a.mask, a.chroms)
    hist = np.zeros(HIST_SHAPE, np.int64)
    for chrom in a.chroms:
        hist += chrom_hist(a.sample, chrom, a.cntx_dir, mask.get(chrom), lengths[chrom])

    # sites by coverage must equal each read set's informative_coverage_dist.tsv
    for ext, axis_keep, lo in [(0, 3, 0), (TRIM_EXT, 2, 1)]:
        got = hist.sum(axis=tuple(i for i in range(4) if i != axis_keep))
        d = coverage_dist(a.cntx_dir, ext, a.sample)
        exp = np.zeros(NC, np.int64)
        for c, n in zip(d["coverage"].astype(int), d["n_sites"].astype(np.int64)):
            exp[min(c, CAP)] += n
        ok = np.array_equal(got[lo:], exp[lo:])
        print(f"{a.sample} ext{ext}: positions by coverage {'match' if ok else 'DIFFER FROM'} "
              f"informative_coverage_dist.tsv (coverage >= {lo}; {int(got[lo:].sum()):,} positions)")
        if not ok:
            sys.exit(f"[error] coverage table and ext{ext} informative_coverage_dist.tsv differ")

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez(a.out, hist=hist, contexts=np.array(CONTEXTS), cap=CAP, sample=a.sample)
    print(f"{a.sample}: {len(a.chroms)} chromosomes -> {a.out}")


def run_derive(a):
    want = set(a.chroms)
    autosome_bp = sum(L for c, L in contig_lengths(a.fai).items() if c in want)
    cov = np.arange(NC)
    ref = pd.read_csv(a.check_untrimmed, sep="\t").set_index(["sample", "floor"])

    den_rows, ctx_rows = [], []
    for path in a.tables:
        z = np.load(path)
        s = str(z["sample"])
        h = z["hist"][0].astype(np.int64)                     # unmasked sites: (context, t, u)
        p99u = p99_of(coverage_dist(a.cntx_dir, 0, s))
        if p99u >= CAP:
            sys.exit(f"[error] {s}: untrimmed p99 {p99u} is not below CAP {CAP}")
        for F in sorted(set(a.floors) | {a.context_floor}):
            win = (cov >= F) & (cov <= p99u)
            hw = h[:, :, win]                                  # (context, t, u in window)
            depth = int((cov[None, :, None] * hw).sum())
            sites = int(hw[:, 1:, :].sum())
            if (s, F) in ref.index:
                r = ref.loc[(s, F)]
                u_depth = int((cov[win] * h.sum(axis=(0, 1))[win]).sum())
                u_sites = int(h.sum(axis=(0, 1))[win].sum())
                if u_depth != int(r.callable_depth) or u_sites != int(r.callable_sites):
                    sys.exit(f"[error] {s} floor {F}: untrimmed sites {u_sites:,} / depth {u_depth:,} "
                             f"!= {a.check_untrimmed} {int(r.callable_sites):,} / {int(r.callable_depth):,}")
            if F in a.floors:
                den_rows.append(dict(sample=s, floor=F, untrimmed_p99=p99u, callable_depth=depth,
                                     callable_sites=sites, autosome_bp=autosome_bp,
                                     pct_autosome_callable=round(100.0 * sites / autosome_bp, 4)))
            if F == a.context_floor:
                bp = hw[:, 1:, :].sum(axis=(1, 2))
                dp = (cov[None, :, None] * hw).sum(axis=(1, 2))
                if int(bp.sum()) != sites or int(dp.sum()) != depth:
                    sys.exit(f"[error] {s}: context totals differ from the denominator")
                for k, c in enumerate(CONTEXTS[:-1]):
                    ctx_rows.append(dict(sample=s, triplet=c, callable_bp=int(bp[k]), callable_depth=int(dp[k]),
                                         floor=F, untrimmed_p99=p99u))
                print(f"{s}: untrimmed p99 {p99u:.0f}; at minimum coverage {F}: callable sites {sites:,}, "
                      f"callable depth {depth:,} (non-ACGT sites {int(bp[NON_ACGT]):,})")
        print(f"{s}: untrimmed callable sites and depth match {a.check_untrimmed} at every floor it has")

    den = pd.DataFrame(den_rows)
    for out, df in [(a.denominator_out, den), (a.context_out, pd.DataFrame(ctx_rows))]:
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
        df.to_csv(out, sep="\t", index=False)
        print(f"-> {out}")
    wide = den.pivot(index="sample", columns="floor", values="callable_depth")
    wide.columns = [f"depth_floor{c}" for c in wide.columns]
    print("\ncallable_depth (trimmed bases, sites by untrimmed coverage):")
    print(wide.to_string())


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="mode", required=True)
    q = sub.add_parser("table")
    q.add_argument("--sample", required=True)
    q.add_argument("--chroms", nargs="+", required=True)
    q.add_argument("--cntx-dir", required=True, help="parent of ext0/ and ext450/")
    q.add_argument("--mask", required=True, help="masked sites BED (0-based half-open)")
    q.add_argument("--fai", required=True)
    q.add_argument("--out", required=True)
    q = sub.add_parser("derive")
    q.add_argument("--tables", nargs="+", required=True)
    q.add_argument("--cntx-dir", required=True)
    q.add_argument("--fai", required=True)
    q.add_argument("--chroms", nargs="+", required=True)
    q.add_argument("--floors", type=int, nargs="+", required=True)
    q.add_argument("--context-floor", type=int, required=True)
    q.add_argument("--check-untrimmed", required=True, help="denominator.tsv of callable_bps_untrimmed.py")
    q.add_argument("--denominator-out", required=True)
    q.add_argument("--context-out", required=True)
    a = p.parse_args()
    {"table": run_table, "derive": run_derive}[a.mode](a)


if __name__ == "__main__":
    main()
