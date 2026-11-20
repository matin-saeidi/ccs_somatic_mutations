#!/usr/bin/env python3
"""
tree_sites.py

Site table and callable length of the tree of the animals (Methods,
Relationships among the animals; Supplementary Fig. 15, Supplementary
Table 7), from the per-chromosome scans of tree_coverage_scan.py and the
untrimmed mutation table.

Sites: every site with a variant in any library, on T2T1-T2T19, in the
untrimmed reads, except sites with more than one non-reference base within a
library or between libraries.

Output (--out-dir)
    nj_sites.tsv     one row per site
        CHROM, POS, REF, ALT
        TRINUCLEOTIDE_CONTEXT         reference trinucleotide centred on POS
        IN_REPEATMASKER_MASK          1 = in the repeat component of the masked sites
        <library>_COVERAGE            untrimmed coverage; 0 = no molecule
        <library>_VARIANT_MOLECULES   untrimmed molecules with ALT; 0 = no variant
    callable_L.tsv   one row per chromosome and minimum coverage, plus CHROM = ALL
        L_OUTSIDE_REPEATS   sites where every library has
                            minimum coverage <= untrimmed coverage <= its 99th percentile,
                            outside the repeats
        L_INSIDE_REPEATS    the same, inside them

Stops if the COVERAGE of a variant in the table differs from the scanned
coverage, or if the sites passing in each library differ from callable_sites
of --untrimmed-denominator (callable_bps_untrimmed.py).

Usage
    python tree_sites.py --scans <dir>/scan/T2T*.npz \\
        --calls <results>/mutation_tables/untrimmed_reads/mutations_untrimmed_reads.tsv \\
        --untrimmed-denominator <results>/callable_bps/untrimmed/denominator.tsv --out-dir <dir>
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

CHROMS = [f"T2T{i}" for i in range(1, 20)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scans", nargs="+", required=True, help="<chrom>.npz of tree_coverage_scan.py")
    ap.add_argument("--calls", required=True, help="untrimmed mutation table of mutation_tables.py")
    ap.add_argument("--untrimmed-denominator", required=True, help="denominator.tsv of callable_bps_untrimmed.py")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    CALLS, DENOM_UNTRIMMED = a.calls, a.untrimmed_denominator

    scans = {os.path.basename(p).split(".")[0]: np.load(p) for p in a.scans}
    scans = {c: scans[c] for c in CHROMS}
    s0 = scans[CHROMS[0]]
    libs = [str(x) for x in s0["libraries"]]
    min_covs = [int(x) for x in s0["min_covs"]]

    # ---- callable length -----------------------------------------------------
    rows = []
    for c in CHROMS:
        for k, C in enumerate(min_covs):
            rows.append(dict(CHROM=c, MIN_COVERAGE=C,
                             L_OUTSIDE_REPEATS=int(scans[c]["L"][k, 0]),
                             L_INSIDE_REPEATS=int(scans[c]["L"][k, 1])))
    L = pd.DataFrame(rows)
    tot = (L.groupby("MIN_COVERAGE", as_index=False)[["L_OUTSIDE_REPEATS", "L_INSIDE_REPEATS"]]
             .sum().assign(CHROM="ALL"))
    L = pd.concat([L, tot[L.columns]], ignore_index=True)

    # ---- check: sites passing in each library = the untrimmed callable sites ---
    den = pd.read_csv(DENOM_UNTRIMMED, sep="\t").set_index(["sample", "floor"]).callable_sites
    single = sum(scans[c]["single"] for c in CHROMS)
    bad = 0
    for k, C in enumerate(min_covs):
        for li, lib in enumerate(libs):
            if (lib, C) in den.index and int(den[(lib, C)]) != int(single[k, li]):
                print(f"[check] {lib} min cov {C}: scan {single[k, li]:,} "
                      f"vs denominator {int(den[(lib, C)]):,}")
                bad += 1
    if bad:
        sys.exit("[error] scan does not reproduce the untrimmed denominator.tsv")
    print("[check] per-library callable positions = the untrimmed denominator.tsv: OK")

    # ---- sites ---------------------------------------------------------------
    t = pd.read_csv(CALLS, sep="\t", dtype={"CHROM": str},
                    usecols=["CHROM", "POS", "REF", "ALT", "TRINUCLEOTIDE_CONTEXT", "LIBRARY", "VARIANT_MOLECULES",
                             "COVERAGE", "MULTIALLELIC_SITE", "IN_REPEATMASKER_MASK"])
    t = t[t.CHROM.isin(CHROMS)]
    n_alt = t.groupby(["CHROM", "POS"]).ALT.nunique()
    multi = t.groupby(["CHROM", "POS"]).MULTIALLELIC_SITE.max()
    drop = (n_alt > 1) | (multi == 1)
    key = pd.MultiIndex.from_frame(t[["CHROM", "POS"]])
    t = t[~key.isin(drop[drop].index)]
    print(f"positions called in any library: {len(n_alt):,}; dropped with more than "
          f"one alternative allele: {int(drop.sum()):,}")

    site = t.drop_duplicates(["CHROM", "POS"])[["CHROM", "POS", "REF", "ALT", "TRINUCLEOTIDE_CONTEXT",
                                              "IN_REPEATMASKER_MASK"]]
    site = (site.assign(_o=site.CHROM.map({c: i for i, c in enumerate(CHROMS)}))
                .sort_values(["_o", "POS"]).drop(columns="_o").reset_index(drop=True))

    cov = np.zeros((len(site), len(libs)), np.int64)
    for c in CHROMS:
        m = (site.CHROM == c).to_numpy()
        p = scans[c]["pos"]
        j = np.searchsorted(p, site.POS.to_numpy()[m])
        if not (p[j] == site.POS.to_numpy()[m]).all():
            sys.exit(f"[error] {c}: site missing from the scan")
        cov[m] = scans[c]["cov"][j]

    alt = (t.pivot_table(index=["CHROM", "POS"], columns="LIBRARY",
                         values="VARIANT_MOLECULES", aggfunc="first")
             .reindex(columns=libs).reindex(pd.MultiIndex.from_frame(site[["CHROM", "POS"]]))
             .fillna(0).astype(int).to_numpy())

    # ---- check: table COVERAGE against the scan --------------------------------
    tc = (t.pivot_table(index=["CHROM", "POS"], columns="LIBRARY",
                        values="COVERAGE", aggfunc="first")
            .reindex(columns=libs).reindex(pd.MultiIndex.from_frame(site[["CHROM", "POS"]]))
            .to_numpy())
    called = ~np.isnan(tc)
    n_bad = int((tc[called] != cov[called]).sum())
    print(f"[check] COVERAGE of {int(called.sum()):,} called rows = scanned coverage: "
          f"{'OK' if n_bad == 0 else f'{n_bad:,} differ'}")
    if n_bad:
        sys.exit("[error] coverage mismatch")

    out = site.copy()
    for li, lib in enumerate(libs):
        out[f"{lib}_COVERAGE"] = cov[:, li]
    for li, lib in enumerate(libs):
        out[f"{lib}_VARIANT_MOLECULES"] = alt[:, li]
    os.makedirs(a.out_dir, exist_ok=True)
    out.to_csv(f"{a.out_dir}/nj_sites.tsv", sep="\t", index=False)
    L.to_csv(f"{a.out_dir}/callable_L.tsv", sep="\t", index=False)

    print(f"\nsites: {len(out):,}")
    print(L[L.CHROM == "ALL"].to_string(index=False))
    print(f"[write] {a.out_dir}/nj_sites.tsv\n[write] {a.out_dir}/callable_L.tsv")


if __name__ == "__main__":
    main()
