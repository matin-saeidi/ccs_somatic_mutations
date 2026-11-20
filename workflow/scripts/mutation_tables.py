#!/usr/bin/env python
# =============================================================================
# mutation_tables.py — every variant of all libraries in one table per read
# set, with the filters as 0/1 columns (1 = fails) instead of applied:
#
#     ext0    ->  <out-dir>/untrimmed_reads/mutations_untrimmed_reads.tsv
#     ext450  ->  <out-dir>/450b_trimmed_reads/mutations_450b_trimmed_reads.tsv
#
# One row per library x site x non-reference base, sorted by chromosome,
# position, library and base. Per-molecule columns hold one value per variant
# molecule, comma-separated, in the order of ZMW.
#
# Filter columns
#   ABOVE_LIBRARY_COVERAGE_P99   coverage above the library's 99th percentile
#   MULTIALLELIC_SITE            another non-reference base at the site in the library
#   IN_REPEATMASKER_MASK         in the repeat component of the masked sites
#   SHARED_WITH_OTHER_LIBRARIES  a variant at the site in another library
# The trimmed table also holds each variant's values in the untrimmed reads
# (UNTRIMMED_* columns). SINGLETONS_ON_MOLECULE is the number of singletons on
# the molecule that pass every filter; the somatic mutations are the
# singletons on molecules with SINGLETONS_ON_MOLECULE <= 2 (Methods, Errors at
# the edges of consensus reads).
#
# Inputs
#   --calls-dir     ext{0,450}/<library>/mutations_with_vaf.tsv (variant_table.py)
#   --calling-dir   ext{0,450}/<library>/informative_coverage_dist.tsv (variant_calling.smk)
#   --repeat-bed    repeat_regions.bed (masked_sites.py)
#   --repeatmasker  RepeatMasker .out of the T2T assembly; REPEATMASKER_CLASS
#                   lists the classes of RM_CLASSES that hold the site
#   --gca-fai, --reference-fai   contig lengths, to match the RepeatMasker
#                   contig names to T2T1..T2T19
#   --animals       animal of each library (ANIMAL column)
#
# Usage
#   python mutation_tables.py --calls-dir <results>/call_tables \
#       --calling-dir <results>/variant_calling \
#       --repeat-bed <results>/masked_sites/repeat_regions.bed \
#       --repeatmasker GCA_964188535.1.repeatMasker.out.gz \
#       --gca-fai GCA_964188535.1_C57BL_6J_T2T_v1_genomic.fna.fai \
#       --reference-fai b6_t2t_mm39XYM.fa.fai \
#       --samples <libraries> --animals <library>=<animal> ... \
#       --chroms T2T1 ... T2T19 --out-dir <results>/mutation_tables
# =============================================================================
import argparse
import gzip
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from callable_bps_untrimmed import weighted_percentile, load_mask, in_mask
from masked_sites import contig_map, merge

# Inputs, set from the command line in main()
CALLS_DIR = CALLING_DIR = RM_OUT = GCA_FAI = REF_FAI = None
CHROMS, SAMPLES = [], []

MAX_COV_PCT = 99
CUTOFF = {"4kb": 2, "10kb": 2}                 # used only for the printed check
RM_CLASSES = ["Simple_repeat", "Low_complexity", "Satellite", "LINE", "SINE", "LTR"]

DATASETS = {
    0:   dict(name="untrimmed_reads",    repeat=None),
    450: dict(name="450b_trimmed_reads", repeat=None),
}

# per-molecule columns taken from the variant table: (variant table column, output column)
MOLECULE_FIELDS = [
    ("fwd_sm",      "FWD_STRAND_SM"),
    ("rev_sm",      "REV_STRAND_SM"),
    ("fwd_sx",      "FWD_STRAND_SX"),
    ("rev_sx",      "REV_STRAND_SX"),
    ("fwd_np",      "FWD_STRAND_NP"),
    ("rev_np",      "REV_STRAND_NP"),
    ("fwd_qpos",    "FWD_POSITION_IN_READ"),
    ("rev_qpos",    "REV_POSITION_IN_READ"),
    ("fwd_readlen", "FWD_READ_LENGTH"),
    ("rev_readlen", "REV_READ_LENGTH"),
]

FILTERS = ["ABOVE_LIBRARY_COVERAGE_P99", "MULTIALLELIC_SITE",
           "IN_REPEATMASKER_MASK", "SHARED_WITH_OTHER_LIBRARIES"]
UNTRIMMED_FILTERS = ["UNTRIMMED_MULTIALLELIC_SITE", "UNTRIMMED_SHARED_WITH_OTHER_LIBRARIES"]
KEY = ["sample", "chrom", "pos_1based", "alt_base"]


ANIMALS = {}                                    # library -> animal, set from --animals


def animal_of(sample):
    """Animal of a library (the library itself when --animals does not list it)."""
    return ANIMALS.get(sample, sample)


def load_repeatmasker():
    """RepeatMasker .out -> {class: {chrom: (starts, ends)}} for RM_CLASSES, merged,
    0-based half-open. The class is the part of class/family before the '/'."""
    cmap = contig_map(GCA_FAI, REF_FAI, CHROMS)
    iv = defaultdict(lambda: defaultdict(list))
    with gzip.open(RM_OUT, "rt") as fh:
        for i, line in enumerate(fh):
            if i < 3:
                continue
            f = line.split()
            if len(f) < 11 or not f[0].isdigit():
                continue
            c, cls = cmap.get(f[4]), f[10].split("/")[0]
            if c is not None and cls in RM_CLASSES:
                iv[cls][c].append((int(f[5]) - 1, int(f[6])))
    out = {}
    for cls, by_chrom in iv.items():
        out[cls] = {}
        for c, v in by_chrom.items():
            m = merge(v)
            out[cls][c] = (np.array([s for s, _ in m], np.int64), np.array([e for _, e in m], np.int64))
    return out


def in_intervals(idx, chrom, pos1):
    """True where the 1-based position is inside an interval of idx (tests pos - 1)."""
    out = np.zeros(len(pos1), bool)
    for c in np.unique(chrom):
        m = chrom == c
        out[m] = in_mask(idx.get(c), pos1[m] - 1)
    return out


def load_calls(ext):
    num = {"pos_1based": int, "n_informative_molecules": int}
    raw = pd.concat([pd.read_csv(f"{CALLS_DIR}/ext{ext}/{s}/mutations_with_vaf.tsv",
                                 sep="\t", dtype=str, keep_default_na=False).astype(num)
                     for s in SAMPLES], ignore_index=True)
    bad = set(raw.chrom) - set(CHROMS)
    if bad:
        sys.exit(f"[error] calls on non-autosomal contigs: {sorted(bad)}")
    return raw


def build_records(ext, rm_classes):
    """One record per library x position x ALT with the site filters of this read set."""
    cfg = DATASETS[ext]
    raw = load_calls(ext)
    agg = dict(ref=("ref_base", "first"), triplet=("triplet", "first"),
               cov=("n_informative_molecules", "first"),
               nvar=("zmw", "size"), zmw=("zmw", list))
    agg.update({c: (c, ",".join) for c, _ in MOLECULE_FIELDS})
    rec = raw.groupby(KEY, sort=False).agg(**agg).reset_index()
    chrom, pos = rec.chrom.to_numpy(), rec.pos_1based.to_numpy()

    # coverage p99 per library
    p99 = {}
    for s in SAMPLES:
        d = pd.read_csv(f"{CALLING_DIR}/ext{ext}/{s}/informative_coverage_dist.tsv", sep="\t")
        p99[s] = weighted_percentile(d["coverage"], d["n_sites"], MAX_COV_PCT)
    rec["p99"] = rec["sample"].map(p99).astype(int)
    rec["above"] = (rec["cov"] > rec["p99"]).astype(int)

    # multiallelic within the library
    rec["multi"] = (rec.groupby(["sample", "chrom", "pos_1based"]).alt_base
                       .transform("nunique") > 1).astype(int)

    # RepeatMasker mask and classes
    rec["repeat"] = in_intervals(load_mask(cfg["repeat"], CHROMS), chrom, pos).astype(int)
    hits = {cls: in_intervals(rm_classes[cls], chrom, pos) for cls in RM_CLASSES}
    rec["rm_class"] = [",".join(cls for cls in hits if hits[cls][i]) for i in range(len(rec))]

    # other libraries with a call at the same position (any ALT)
    libs = rec.groupby(["chrom", "pos_1based"])["sample"].agg(set).rename("libs")
    rec = rec.join(libs, on=["chrom", "pos_1based"])
    rec["others"] = [",".join(l for l in SAMPLES if l in L and l != s)
                     for L, s in zip(rec.libs, rec["sample"])]
    rec["shared"] = (rec.others != "").astype(int)
    return rec, raw


def singletons_per_molecule(rec, passed):
    """{(sample, zmw): calls on the molecule, at sites where it is the only molecule carrying the
    non-reference allele, that pass `passed`}."""
    sing = rec[passed & (rec.nvar == 1)]
    return sing.groupby([sing["sample"], sing.zmw.str[0]]).size().to_dict()


def per_zmw(rec, per_mol):
    return [",".join(str(per_mol.get((s, z), 0)) for z in zs) for s, zs in zip(rec["sample"], rec.zmw)]


def add_untrimmed(rec, urec, u_per_mol):
    """Untrimmed values of each trimmed call. Every trimmed call is also an untrimmed call."""
    u = urec[KEY + ["nvar", "cov", "p99", "above", "multi", "shared", "others"]].rename(
        columns={c: f"u_{c}" for c in ["nvar", "cov", "p99", "above", "multi", "shared", "others"]})
    rec = rec.merge(u, on=KEY, how="left", validate="one_to_one")
    if rec.u_nvar.isna().any():
        sys.exit(f"[error] {int(rec.u_nvar.isna().sum())} trimmed calls have no untrimmed call")
    rec["u_spm"] = per_zmw(rec, u_per_mol)
    return rec


def to_table(rec, untrimmed):
    """Output table of one read set from its records."""
    out = pd.DataFrame({
        "CHROM": rec.chrom, "POS": rec.pos_1based, "REF": rec.ref, "ALT": rec.alt_base,
        "LIBRARY": rec["sample"], "ANIMAL": rec["sample"].map(animal_of),
        "TRINUCLEOTIDE_CONTEXT": rec.triplet,
        "VARIANT_MOLECULES": rec.nvar, "COVERAGE": rec["cov"],
        "VAF": (rec.nvar / rec["cov"]).round(6),
        "LIBRARY_COVERAGE_P99": rec.p99, "ABOVE_LIBRARY_COVERAGE_P99": rec.above,
        "MULTIALLELIC_SITE": rec.multi,
        "IN_REPEATMASKER_MASK": rec.repeat, "REPEATMASKER_CLASS": rec.rm_class,
        "SHARED_WITH_OTHER_LIBRARIES": rec.shared, "SHARED_WITH_LIBRARIES": rec.others,
        **({"UNTRIMMED_VARIANT_MOLECULES": rec.u_nvar.astype(int),
            "UNTRIMMED_COVERAGE": rec.u_cov.astype(int),
            "UNTRIMMED_LIBRARY_COVERAGE_P99": rec.u_p99.astype(int),
            "UNTRIMMED_ABOVE_LIBRARY_COVERAGE_P99": rec.u_above.astype(int),
            "UNTRIMMED_MULTIALLELIC_SITE": rec.u_multi.astype(int),
            "UNTRIMMED_SHARED_WITH_OTHER_LIBRARIES": rec.u_shared.astype(int),
            "UNTRIMMED_SHARED_WITH_LIBRARIES": rec.u_others} if untrimmed else {}),
        "ZMW": rec.zmw.str.join(","), "SINGLETONS_ON_MOLECULE": rec.spm,
        **({"UNTRIMMED_SINGLETONS_ON_MOLECULE": rec.u_spm} if untrimmed else {}),
        **{name: rec[c] for c, name in MOLECULE_FIELDS},
    })
    order = {c: i for i, c in enumerate(CHROMS)}
    lib_order = {s: i for i, s in enumerate(SAMPLES)}
    out = (out.assign(_c=out.CHROM.map(order), _l=out.LIBRARY.map(lib_order))
              .sort_values(["_c", "POS", "_l", "ALT"])
              .drop(columns=["_c", "_l"])
              .reset_index(drop=True))
    return out


def site_passed(rec, untrimmed):
    """Calls that pass every site filter (SINGLETONS_ON_MOLECULE counts these).
    Trimmed table: the coverage filter uses UNTRIMMED_ABOVE_LIBRARY_COVERAGE_P99."""
    cap = "u_above" if untrimmed else "above"
    ok = rec[[cap, "multi", "repeat", "shared"]].sum(axis=1) == 0
    if untrimmed:
        ok &= (rec.u_nvar == 1) & (rec.u_multi == 0) & (rec.u_shared == 0)
    return ok


def write(tab, raw, ext, out_root):
    cfg = DATASETS[ext]
    out_dir = f"{out_root}/{cfg['name']}"
    os.makedirs(out_dir, exist_ok=True)
    path = f"{out_dir}/mutations_{cfg['name']}.tsv"
    tab.to_csv(path, sep="\t", index=False)

    trimmed = "UNTRIMMED_COVERAGE" in tab
    if trimmed:
        filters = ([f for f in FILTERS if f != "ABOVE_LIBRARY_COVERAGE_P99"]
                   + ["UNTRIMMED_ABOVE_LIBRARY_COVERAGE_P99"] + UNTRIMMED_FILTERS)
    else:
        filters = FILTERS
    s = tab[(tab[filters].sum(axis=1) == 0) & (tab.VARIANT_MOLECULES == 1)]
    if trimmed:
        s = s[s.UNTRIMMED_VARIANT_MOLECULES == 1]
    cut = np.where(s.LIBRARY.str.endswith("10kb"), CUTOFF["10kb"], CUTOFF["4kb"])
    cov = s.UNTRIMMED_COVERAGE if trimmed else s.COVERAGE
    final = s[(s.SINGLETONS_ON_MOLECULE.astype(int) <= cut) & (cov >= 10)]
    print(f"ext{ext}: {len(raw):,} molecules -> {len(tab):,} rows, "
          f"final set {len(final):,} -> {path}")


def main():
    global CALLS_DIR, CALLING_DIR, RM_OUT, GCA_FAI, REF_FAI, CHROMS, SAMPLES
    p = argparse.ArgumentParser(description="Write the two mutation tables.")
    p.add_argument("--calls-dir", required=True, help="parent of ext0/ and ext450/ of variant_table.py")
    p.add_argument("--calling-dir", required=True, help="parent of ext0/ and ext450/ of variant_calling.smk")
    p.add_argument("--repeat-bed", required=True, help="repeat_regions.bed of masked_sites.py")
    p.add_argument("--repeatmasker", required=True, help="RepeatMasker .out.gz of the T2T assembly")
    p.add_argument("--gca-fai", required=True, help=".fai of the T2T assembly FASTA")
    p.add_argument("--reference-fai", required=True, help=".fai of the reference FASTA")
    p.add_argument("--samples", nargs="+", required=True, help="libraries, in the order of the rows")
    p.add_argument("--animals", nargs="+", default=None, metavar="LIBRARY=ANIMAL",
                   help="animal of each library (ANIMAL column); default: the library")
    p.add_argument("--chroms", nargs="+", required=True)
    p.add_argument("--out-dir", required=True)
    args = p.parse_args()
    if args.animals:
        ANIMALS.update(a.split("=", 1) for a in args.animals)
    CALLS_DIR, CALLING_DIR = args.calls_dir, args.calling_dir
    RM_OUT, GCA_FAI, REF_FAI = args.repeatmasker, args.gca_fai, args.reference_fai
    CHROMS, SAMPLES = list(args.chroms), list(args.samples)
    for ext in DATASETS:
        DATASETS[ext]["repeat"] = args.repeat_bed

    rm_classes = load_repeatmasker()

    urec, uraw = build_records(0, rm_classes)
    u_per_mol = singletons_per_molecule(urec, site_passed(urec, False))
    urec["spm"] = per_zmw(urec, u_per_mol)
    write(to_table(urec, False), uraw, 0, args.out_dir)

    rec, raw = build_records(450, rm_classes)
    rec = add_untrimmed(rec, urec, u_per_mol)
    rec["spm"] = per_zmw(rec, singletons_per_molecule(rec, site_passed(rec, True)))
    write(to_table(rec, True), raw, 450, args.out_dir)


if __name__ == "__main__":
    main()
