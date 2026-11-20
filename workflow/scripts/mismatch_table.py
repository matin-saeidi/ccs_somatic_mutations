#!/usr/bin/env python3
"""
mismatch_table.py

Table of the mismatches between the two strand consensus reads of each
molecule (types f and r of variant_calling.py: only one read carries a
non-reference base; Methods, Estimation of the sequencing error rate), one row
per mismatch, for one library and read set. With --discordant-out, also the
sites where both reads differ from the reference and from each other (type d_c).

Input:  <mut-dir>/<sample>/<chrom>.ext<ext>.qual<qual>.muts.bed.gz from
        variant_calling.py, for every chromosome of --chroms
Output: --out (gzipped if it ends in .gz), columns in HEADER below:
        err_*   the read with the non-reference base (qpos from its own 5' end,
                err_dist_5p / err_dist_3p from its two ends)
        other_* the read with the reference base
        --discordant-out: both reads' fields in full (fwd_*, rev_*)

Usage:
    python mismatch_table.py --mut-dir <results>/variant_calling/ext450 \\
        --samples 60wk_male_ctx_2_4kb --chroms T2T1 ... T2T19 \\
        --ext 450 --qual 40 --out single_strand_calls.tsv.gz \\
        --discordant-out discordant_calls.tsv.gz
"""

import argparse
import gzip
import os
import sys

HEADER = [
    "sample", "chrom", "pos_1based", "ref_base", "triplet",
    "call_type",        # f or r
    "err_strand",       # fwd or rev: the read with the non-reference base
    "err_base", "err_qual", "err_sm", "err_sx",
    "err_qpos",         # position in that read, from its own 5' end
    "err_readlen",
    "err_dist_5p",      # err_qpos
    "err_dist_3p",      # err_readlen - 1 - err_qpos
    "err_divergence", "err_np",
    "other_base", "other_qual", "other_qpos", "other_readlen", "other_np",
    "zmw",
    "other_sm", "other_sx", "other_divergence",
]
STRAND_FIELDS = ["base", "qual", "sm", "sx", "qpos", "readlen", "divergence", "np"]
DISCORDANT_HEADER = (["sample", "chrom", "pos_1based", "ref_base", "triplet"]
                     + [f"fwd_{x}" for x in STRAND_FIELDS]
                     + [f"rev_{x}" for x in STRAND_FIELDS]
                     + ["zmw"])


def iter_calls(muts_path, want):
    """Yield (fields, mtype) for lines whose call type is in `want`."""
    if not os.path.exists(muts_path):
        return
    with gzip.open(muts_path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            if len(f) < 8:
                continue
            if f[7] in want:
                yield f, f[7]


def main():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mut-dir", required=True)
    p.add_argument("--samples", nargs="+", required=True)
    p.add_argument("--chroms", nargs="+", required=True)
    p.add_argument("--ext",  type=int, default=0)
    p.add_argument("--qual", type=int, default=40)
    p.add_argument("--call-types", nargs="+", default=["f", "r"],
                   choices=["f", "r"],
                   help="mismatch types to write (default: f r)")
    p.add_argument("--out", required=True,
                   help="output path; gzipped if it ends in .gz")
    p.add_argument("--discordant-out", default=None,
                   help="also write the d_c calls to this path; gzipped if it ends in .gz")
    p.add_argument("--allow-missing-chroms", action="store_true")
    args = p.parse_args()

    want    = set(args.call_types) | ({"d_c"} if args.discordant_out else set())
    missing = []
    n_out   = 0
    counts  = {t: 0 for t in want}

    opener = lambda path: (gzip.open if path.endswith(".gz") else open)(path, "wt")
    dc_out = opener(args.discordant_out) if args.discordant_out else None
    if dc_out:
        dc_out.write("\t".join(DISCORDANT_HEADER) + "\n")
    with opener(args.out) as out:
        out.write("\t".join(HEADER) + "\n")

        for sample in args.samples:
            sys.stderr.write(f"Processing {sample}...\n")
            for chrom in args.chroms:
                prefix = f"{chrom}.ext{args.ext}.qual{args.qual}"
                path   = os.path.join(args.mut_dir, sample, f"{prefix}.muts.bed.gz")
                if not os.path.exists(path):
                    missing.append(f"{sample}/{chrom}")
                    continue

                for f, mtype in iter_calls(path, want):
                    fwd = f[4].split(":")     # base:qual:sm:sx:qpos:readlen:div:np
                    rev = f[5].split(":")
                    if mtype == "d_c":
                        dc_out.write("\t".join([sample, f[0], str(int(f[1]) + 1), f[3][1], f[3]]
                                                + fwd + rev + [f[6]]) + "\n")
                        counts[mtype] += 1
                        continue
                    # the read whose base differs from the reference
                    if mtype == "f":
                        err, other, strand = fwd, rev, "fwd"
                    else:
                        err, other, strand = rev, fwd, "rev"

                    qpos    = int(err[4])
                    readlen = int(err[5])
                    out.write("\t".join([
                        sample, f[0], str(int(f[1]) + 1),
                        f[3][1], f[3],                    # ref base, triplet
                        mtype, strand,
                        err[0], err[1], err[2], err[3],
                        str(qpos), str(readlen),
                        str(qpos), str(readlen - 1 - qpos),
                        err[6], err[7],
                        other[0], other[1], other[4], other[5], other[7],
                        f[6],
                        other[2], other[3], other[6],
                    ]) + "\n")
                    n_out += 1
                    counts[mtype] += 1

    if dc_out:
        dc_out.close()
    sys.stderr.write(f"  wrote {n_out:,} calls  ({', '.join(f'{t}={c:,}' for t, c in sorted(counts.items()))})"
                     f" -> {args.out}\n")

    if missing:
        sys.stderr.write(f"  {'WARNING' if args.allow_missing_chroms else 'ERROR'}: "
                         f"missing inputs: {missing[:10]}"
                         f"{' ...' if len(missing) > 10 else ''}\n")
        if not args.allow_missing_chroms:
            sys.exit("[error] missing per-chromosome inputs; refusing to write a "
                     "partial table. Pass --allow-missing-chroms to override.")


if __name__ == "__main__":
    main()
