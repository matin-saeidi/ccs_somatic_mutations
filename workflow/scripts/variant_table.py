#!/usr/bin/env python3
"""
variant_table.py

Table of the variants (type d) of one library and read set, one row per
molecule carrying the non-reference allele, with the coverage and VAF of the site: mutations_with_vaf.tsv.

Input:  <mut-dir>/<sample>/<chrom>.ext<ext>.qual<qual>.muts.bed.gz and
        .cntx.bed.gz from variant_calling.py, for every chromosome of --chroms
        cntx columns: chrom  start  end  total_cov  informative_cov  triplet
Output: --vaf-out (default stdout), columns:
        sample chrom pos_1based ref_base alt_base triplet
        fwd_/rev_ base qual sm sx qpos readlen divergence np
        zmw n_variant_molecules n_informative_molecules n_total_molecules vaf
        n_informative_molecules is the coverage of the site; vaf =
        n_variant_molecules / n_informative_molecules.

Usage:
    python variant_table.py --mut-dir <results>/variant_calling/ext0 \\
        --samples 60wk_male_ctx_2_4kb --chroms T2T1 T2T2 ... T2T19 \\
        --ext 0 --qual 40 --vaf-out mutations_with_vaf.tsv
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict


def load_mutations(muts_path):
    """
    Returns:
        d_records : list of dicts, one per d-type line
        d_counts  : dict (chrom, start) -> count of d lines at that position
    """
    d_records = []
    d_counts  = defaultdict(int)

    if not os.path.exists(muts_path):
        return d_records, d_counts

    with gzip.open(muts_path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            if fields[7] != "d":
                continue

            key = (fields[0], int(fields[1]))
            d_counts[key] += 1
            d_records.append({
                "chrom":    fields[0],
                "start":    int(fields[1]),
                "triplet":  fields[3],
                "fwd_info": fields[4],
                "rev_info": fields[5],
                "zmw":      fields[6],
            })

    return d_records, d_counts


def scan_coverage(cntx_path, needed):
    """(total_cov, informative_cov) of the positions in `needed`, from one
    cntx.bed.gz."""
    needed_cov = {}
    if not os.path.exists(cntx_path):
        return needed_cov

    with gzip.open(cntx_path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            key = (fields[0], int(fields[1]))
            if key in needed:
                needed_cov[key] = (int(fields[3]), int(fields[4]))

    return needed_cov


def compute_sample(sample, mut_dir, ext, qual, chroms):
    all_d_mutations = []
    missing         = []

    for chrom in chroms:
        prefix    = f"{chrom}.ext{ext}.qual{qual}"
        muts_path = os.path.join(mut_dir, sample, f"{prefix}.muts.bed.gz")
        cntx_path = os.path.join(mut_dir, sample, f"{prefix}.cntx.bed.gz")

        if not os.path.exists(muts_path) or not os.path.exists(cntx_path):
            missing.append(chrom)
            continue

        d_records, d_counts = load_mutations(muts_path)
        needed              = {(r["chrom"], r["start"]) for r in d_records}
        cov                 = scan_coverage(cntx_path, needed)

        for rec in d_records:
            key            = (rec["chrom"], rec["start"])
            n_total, n_inf = cov.get(key, (0, 0))
            n_var          = d_counts[key]

            rec["n_total_molecules"]       = n_total
            rec["n_informative_molecules"] = n_inf
            rec["n_variant_molecules"]     = n_var
            rec["vaf"] = n_var / n_inf if n_inf > 0 else float("nan")

            all_d_mutations.append(rec)

    return all_d_mutations, missing


def write_vaf_table(fh, records):
    """
    fwd_info / rev_info from variant_calling.py:
        base : qual : sm : sx : qpos : readlen : divergence : np
    """
    header = [
        "sample", "chrom", "pos_1based", "ref_base", "alt_base", "triplet",
        "fwd_base", "fwd_qual", "fwd_sm", "fwd_sx",
        "fwd_qpos", "fwd_readlen", "fwd_divergence", "fwd_np",
        "rev_base", "rev_qual", "rev_sm", "rev_sx",
        "rev_qpos", "rev_readlen", "rev_divergence", "rev_np",
        "zmw",
        "n_variant_molecules",
        "n_informative_molecules",   # coverage of the site (m + d)
        "n_total_molecules",         # cntx total_cov
        "vaf",
    ]
    fh.write("\t".join(header) + "\n")

    for rec in sorted(records,
                      key=lambda r: (r["sample"], r["chrom"], r["start"])):
        fwd      = rec["fwd_info"].split(":")
        rev      = rec["rev_info"].split(":")
        ref_base = rec["triplet"][1]
        alt_base = fwd[0]

        fh.write("\t".join([
            rec["sample"],
            rec["chrom"],
            str(rec["start"] + 1),
            ref_base,
            alt_base,
            rec["triplet"],
            fwd[0], fwd[1], fwd[2], fwd[3], fwd[4], fwd[5], fwd[6], fwd[7],
            rev[0], rev[1], rev[2], rev[3], rev[4], rev[5], rev[6], rev[7],
            rec["zmw"],
            str(rec["n_variant_molecules"]),
            str(rec["n_informative_molecules"]),
            str(rec["n_total_molecules"]),
            f"{rec['vaf']:.6f}",
        ]) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mut-dir", required=True)
    p.add_argument("--samples", nargs="+", required=True)
    p.add_argument("--chroms", nargs="+", required=True,
                   help="contigs to read, e.g. T2T1 .. T2T19")
    p.add_argument("--ext",  type=int, default=0)
    p.add_argument("--qual", type=int, default=40)
    p.add_argument("--vaf-out", default=None,
                   help="output path for mutations_with_vaf.tsv (default: stdout)")
    p.add_argument("--allow-missing-chroms", action="store_true",
                   help="warn instead of stopping when a chromosome's files are missing")
    args = p.parse_args()

    all_vaf_records = []
    any_missing     = False

    for sample in args.samples:
        print(f"Processing {sample}...", file=sys.stderr)

        mut_records, missing = compute_sample(
            sample, args.mut_dir, args.ext, args.qual, args.chroms
        )

        if missing:
            any_missing = True
            print(f"  {'WARNING' if args.allow_missing_chroms else 'ERROR'} "
                  f"[{sample}]: missing inputs for chromosomes: {missing}",
                  file=sys.stderr)

        print(f"  {len(mut_records):,} d-type mutations loaded", file=sys.stderr)

        for rec in mut_records:
            rec["sample"] = sample
        all_vaf_records.extend(mut_records)

        del mut_records

    if any_missing and not args.allow_missing_chroms:
        sys.exit("[error] missing per-chromosome inputs; refusing to write a "
                 "partial table. Pass --allow-missing-chroms to override.")

    if args.vaf_out:
        with open(args.vaf_out, "w") as fh:
            write_vaf_table(fh, all_vaf_records)
        print(f"VAF table written to: {args.vaf_out}", file=sys.stderr)
    else:
        write_vaf_table(sys.stdout, all_vaf_records)


if __name__ == "__main__":
    main()
