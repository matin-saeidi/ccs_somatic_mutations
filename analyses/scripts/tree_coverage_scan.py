#!/usr/bin/env python3
"""
tree_coverage_scan.py

Untrimmed coverage of every library on one autosome, for the tree of the
animals (Methods, Relationships among the animals; Supplementary Fig. 15).

For each minimum coverage C, counts the sites where every library has

    C <= untrimmed coverage <= that library's 99th percentile of untrimmed coverage

inside and outside the RepeatMasker repeats. Also records every library's
untrimmed coverage at each site with a variant in any library (untrimmed
mutation table), and the sites passing in each library on its own outside the
masked sites.

Output  <out>.npz
    libraries     the libraries
    p99           99th percentile of untrimmed coverage of each library
    min_covs      the minimum coverages
    L             int64 (n_min_covs, 2): sites passing in every library,
                  column 0 outside the repeats, column 1 inside
    single        int64 (n_min_covs, n_libraries): sites passing in each library,
                  outside the masked sites
    pos           int64 (n_pos,): sites with a variant, 1-based
    cov           int32 (n_pos, n_libraries): untrimmed coverage at those sites,
                  0 where the library has no molecule

Usage
    python tree_coverage_scan.py --chrom T2T19 --libraries <libraries> \\
        --cntx-dir <results>/variant_calling/ext0 --fai b6_t2t_mm39XYM.fa.fai \\
        --repeats <results>/masked_sites/repeat_regions.bed \\
        --mask <results>/masked_sites/excluded_regions.bed \\
        --denominator <results>/callable_bps/trimmed/denominator.tsv \\
        --calls <results>/mutation_tables/untrimmed_reads/mutations_untrimmed_reads.tsv \\
        --out T2T19.npz
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "pipeline", "scripts"))
from callable_bps_untrimmed import cntx_chunks  # noqa: E402


def bed_bool(path, chrom, n):
    """Boolean array over the chromosome, True inside the BED intervals."""
    b = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                    names=["c", "s", "e"], dtype={"c": str})
    b = b[b.c == chrom]
    m = np.zeros(n, bool)
    for s, e in zip(b.s.to_numpy(), b.e.to_numpy()):
        m[s:e] = True
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-covs", type=int, nargs="+", default=[8, 9, 10])
    ap.add_argument("--libraries", nargs="+", required=True)
    ap.add_argument("--cntx-dir", required=True, help="ext0/ of variant_calling.smk")
    ap.add_argument("--fai", required=True)
    ap.add_argument("--repeats", required=True, help="repeat_regions.bed of masked_sites.py")
    ap.add_argument("--mask", required=True, help="excluded_regions.bed of masked_sites.py")
    ap.add_argument("--denominator", required=True, help="denominator.tsv of callable_bps.py (untrimmed_p99)")
    ap.add_argument("--calls", required=True, help="untrimmed mutation table of mutation_tables.py")
    a = ap.parse_args()
    LIBS = list(a.libraries)

    fai = pd.read_csv(a.fai, sep="\t", header=None, usecols=[0, 1], names=["c", "n"])
    n = int(fai.set_index("c").n[a.chrom])

    den = pd.read_csv(a.denominator, sep="\t")
    p99 = den.groupby("sample").untrimmed_p99.agg(["min", "max"])
    if (p99["min"] != p99["max"]).any():
        sys.exit("[error] untrimmed_p99 differs between minimum coverages")
    p99 = p99["min"].astype(int).reindex(LIBS).to_numpy()

    calls = pd.read_csv(a.calls, sep="\t", usecols=["CHROM", "POS"], dtype={"CHROM": str})
    pos = np.unique(calls.loc[calls.CHROM == a.chrom, "POS"].to_numpy(np.int64))
    start = pos - 1
    cov = np.zeros((len(pos), len(LIBS)), np.int32)

    rep = bed_bool(a.repeats, a.chrom, n)
    full = bed_bool(a.mask, a.chrom, n)
    cnt = np.zeros((len(a.min_covs), n), np.uint8)
    single = np.zeros((len(a.min_covs), len(LIBS)), np.int64)

    for li, lib in enumerate(LIBS):
        path = f"{a.cntx_dir}/{lib}/{a.chrom}.ext0.qual40.cntx.bed.gz"
        for ch in cntx_chunks(path, usecols=("start", "informative_cov")):
            s = ch["start"].to_numpy()
            u = ch["informative_cov"].to_numpy()
            ok = u <= p99[li]
            for k, C in enumerate(a.min_covs):
                sel = s[ok & (u >= C)]
                cnt[k, sel] += 1
                single[k, li] += int((~full[sel]).sum())
            j = np.searchsorted(start, s)
            j[j == len(start)] = 0
            hit = start[j] == s if len(start) else np.zeros(len(s), bool)
            cov[j[hit], li] = u[hit]
        print(f"[{a.chrom}] {lib} done", flush=True)

    L = np.zeros((len(a.min_covs), 2), np.int64)
    for k in range(len(a.min_covs)):
        all8 = cnt[k] == len(LIBS)
        L[k, 0] = int((all8 & ~rep).sum())
        L[k, 1] = int((all8 & rep).sum())

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez(a.out, libraries=np.array(LIBS), p99=p99, min_covs=np.array(a.min_covs),
             L=L, single=single, pos=pos, cov=cov)
    print(f"[{a.chrom}] L outside repeats: "
          + "  ".join(f"C{C} {L[k, 0]:,}" for k, C in enumerate(a.min_covs)), flush=True)
    print(f"[write] {a.out}")


if __name__ == "__main__":
    main()
