#!/usr/bin/env python3
"""
masked_sites.py

Masked sites (Methods, Masked sites): the union of
    repeats   RepeatMasker classes --classes (Simple_repeat, Low_complexity,
              Satellite*; a class ending in * matches every class with that
              prefix, e.g. Satellite/centr)
    shared    sites with a variant called in >= --min-shared animals, in the
              untrimmed reads (Results, Identifying somatic mutations)

The RepeatMasker .out names contigs by GenBank accession (OZ077242.1 ...);
they are matched to T2T1..T2T19 by sequence length (--gca-fai, --hybrid-fai).
Libraries of the same animal (--animals) count as one animal.

Output (--out-dir; BED, 0-based half-open):
    excluded_regions.bed   the masked sites
    repeat_regions.bed     the repeat component
    shared_sites.bed       the shared component
    mask_summary.tsv       sizes of each component; shared_unmasked is 0

Usage:
    python masked_sites.py --repeatmasker GCA_964188535.1.repeatMasker.out.gz \\
        --gca-fai GCA_964188535.1_C57BL_6J_T2T_v1_genomic.fna.fai \\
        --hybrid-fai b6_t2t_mm39XYM.fa.fai --rates-dir <results>/call_tables/ext0 \\
        --samples <libraries> --animals <library>=<animal> ... --chroms T2T1 ... T2T19 \\
        --classes Simple_repeat Low_complexity 'Satellite*' --min-shared 2 \\
        --out-dir <results>/masked_sites
"""
import argparse
import gzip
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

DEFAULT_CLASSES = ["Simple_repeat", "Low_complexity", "Satellite*"]


def contig_map(gca_fai, hyb_fai, chroms):
    """GCA accession -> hybrid contig name, matched on sequence length."""
    g = pd.read_csv(gca_fai, sep="\t", header=None, usecols=[0, 1], names=["n", "L"])
    h = pd.read_csv(hyb_fai, sep="\t", header=None, usecols=[0, 1], names=["n", "L"])
    h = h[h.n.isin(chroms)]
    m = h.merge(g, on="L", suffixes=("_hyb", "_gca"))
    if len(m) != len(chroms):
        sys.exit(f"[error] matched {len(m)}/{len(chroms)} contigs by length; "
                 "cannot map RepeatMasker coordinates safely")
    dup = m.L.duplicated().any()
    if dup:
        sys.exit("[error] two contigs share a length; length matching is ambiguous")
    return dict(zip(m.n_gca, m.n_hyb))


def read_repeats(rm_path, cmap, classes):
    """RepeatMasker .out -> {contig: [(start0, end)]} in BED half-open coords.
    A class ending in '*' is a prefix match."""
    want = {c for c in classes if not c.endswith("*")}
    prefixes = tuple(c[:-1] for c in classes if c.endswith("*"))
    out = defaultdict(list)
    kept = total = 0
    with gzip.open(rm_path, "rt") as fh:
        for i, line in enumerate(fh):
            if i < 3:
                continue
            f = line.split()
            if len(f) < 11 or not f[0].isdigit():
                continue
            total += 1
            c = cmap.get(f[4])
            cls_ = f[10]
            if c is None or not (cls_ in want or
                                 (prefixes and cls_.startswith(prefixes))):
                continue
            # RepeatMasker begin/end are 1-based inclusive -> BED half-open
            out[c].append((int(f[5]) - 1, int(f[6])))
            kept += 1
    shown = sorted(want) + [p_ + "*" for p_ in prefixes]
    sys.stderr.write(f"  repeat records: {kept:,} kept of {total:,} "
                     f"(classes: {', '.join(shown)})\n")
    return out


def merge(intervals):
    """Sort and merge overlapping/adjacent [start, end) intervals."""
    if not intervals:
        return []
    intervals.sort()
    out = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def shared_positions(rates_dir, samples, chroms, min_shared, animal_of):
    """Positions called in >= min_shared distinct animals, from the variant tables."""
    frames = []
    for s in samples:
        p = os.path.join(rates_dir, s, "mutations_with_vaf.tsv")
        if not os.path.exists(p):
            sys.exit(f"[error] missing call table: {p}")
        d = pd.read_csv(p, sep="\t", dtype={"chrom": str}, low_memory=False,
                        usecols=["chrom", "pos_1based"])
        d = d[d.chrom.isin(chroms)].drop_duplicates()
        d["animal"] = animal_of(s)
        frames.append(d)
    a = pd.concat(frames).drop_duplicates()
    k = a.groupby(["chrom", "pos_1based"])["animal"].nunique()
    return k[k >= min_shared].index.to_frame(index=False), k


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repeatmasker", required=True, help="*.repeatMasker.out.gz")
    ap.add_argument("--gca-fai", required=True)
    ap.add_argument("--hybrid-fai", required=True)
    ap.add_argument("--rates-dir", required=True,
                    help="dir of <sample>/mutations_with_vaf.tsv")
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--animals", nargs="+", default=None, metavar="LIBRARY=ANIMAL",
                    help="animal of each library; default: each library is its own animal")
    ap.add_argument("--chroms", nargs="+", required=True)
    ap.add_argument("--classes", nargs="+", default=DEFAULT_CLASSES)
    ap.add_argument("--min-shared", type=int, default=2,
                    help="a position is masked if called in >= this many animals")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    chroms = list(args.chroms)

    # animal of each library
    animals = dict(a.split("=", 1) for a in args.animals) if args.animals else {}
    animal_of = lambda s: animals.get(s, s)

    cmap = contig_map(args.gca_fai, args.hybrid_fai, chroms)
    rep = read_repeats(args.repeatmasker, cmap, args.classes)
    rep = {c: merge(v) for c, v in rep.items()}
    rep_bp = sum(e - s for v in rep.values() for s, e in v)
    rep_n = sum(len(v) for v in rep.values())
    sys.stderr.write(f"  repeat component: {rep_n:,} merged intervals, {rep_bp:,} bp\n")

    shared, kcnt = shared_positions(args.rates_dir, args.samples, chroms,
                                    args.min_shared, animal_of)
    sys.stderr.write(f"  recurrence component: {len(shared):,} positions "
                     f"(>= {args.min_shared} animals of "
                     f"{len(set(map(animal_of, args.samples)))})\n")

    # union: repeats + 1 bp per shared position
    allv = {c: list(v) for c, v in rep.items()}
    for c, p in zip(shared.chrom, shared.pos_1based):
        allv.setdefault(c, []).append((p - 1, p))       # 1-based P -> [P-1, P)
    allv = {c: merge(v) for c, v in allv.items()}

    def write_bed(path, d):
        with open(path, "w") as fh:
            for c in chroms:
                for s, e in d.get(c, []):
                    fh.write(f"{c}\t{s}\t{e}\n")

    write_bed(os.path.join(args.out_dir, "repeat_regions.bed"), rep)
    write_bed(os.path.join(args.out_dir, "excluded_regions.bed"), allv)
    with open(os.path.join(args.out_dir, "shared_sites.bed"), "w") as fh:
        for c, p in zip(shared.chrom, shared.pos_1based):
            fh.write(f"{c}\t{p-1}\t{p}\n")

    # check: every shared position is inside the masked sites
    unmasked = 0
    for c, grp in shared.groupby("chrom"):
        iv = allv.get(c, [])
        if not iv:
            unmasked += len(grp); continue
        st = np.array([x[0] for x in iv]); en = np.array([x[1] for x in iv])
        pos = grp.pos_1based.to_numpy() - 1
        i = np.searchsorted(st, pos, side="right") - 1
        ok = (i >= 0) & (pos < en[np.clip(i, 0, None)])
        unmasked += int((~ok).sum())

    tot_bp = sum(e - s for v in allv.values() for s, e in v)
    tot_n = sum(len(v) for v in allv.values())
    summ = pd.DataFrame([dict(
        repeat_classes=",".join(args.classes),
        repeat_intervals=rep_n, repeat_bp=rep_bp,
        shared_positions=len(shared), min_shared=args.min_shared,
        shared_already_in_repeats=int(len(shared) - (tot_bp - rep_bp)),
        total_intervals=tot_n, total_bp=tot_bp,
        shared_unmasked=unmasked)])
    summ.to_csv(os.path.join(args.out_dir, "mask_summary.tsv"), sep="\t", index=False)
    sys.stderr.write(f"  MASK: {tot_n:,} intervals, {tot_bp:,} bp\n")
    sys.stderr.write(f"  shared_unmasked = {unmasked} (must be 0)\n")
    if unmasked:
        sys.exit("[error] some shared positions are not inside the mask")


if __name__ == "__main__":
    main()
