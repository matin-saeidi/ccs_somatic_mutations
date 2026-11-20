#!/usr/bin/env python3
"""
joint_coverage.py

The two libraries of one animal (the short and the long library of the
60-week-old mouse), site by site, outside the masked sites: both libraries'
untrimmed coverage, with the trimmed coverage of each library there (Methods,
Concordance between the short and long libraries of the 60-week-old mouse;
Supplementary Fig. 8).

Modes
    positions  every site with a variant of either library in the two mutation
               tables -> <out-dir>/call_positions.tsv (chrom, pos_1based)
    chrom      one chromosome -> <out-dir>/parts/<chrom>.npz
    merge      all parts -> <out-dir>/joint_untrimmed_coverage.npz, then removes parts/

Output  joint_untrimmed_coverage.npz
    n_pos    int64 (CAP+1, CAP+1): n_pos[a, b] = unmasked sites with untrimmed
             coverage a in the short library and b in the long library (above CAP
             counted at CAP; sites with 0 in both libraries not counted)
    t4_sum   int64 (CAP+1, CAP+1): the short library's trimmed coverage summed over those sites
    t10_sum  int64 (CAP+1, CAP+1): the same for the long library
    lookup_* at each site of call_positions.tsv: chrom, pos, t4, t10 (trimmed
             coverage), u4, u10 (untrimmed coverage), masked
Trimmed coverage at the sites chosen by untrimmed coverage, e.g. 10..P4 in the
short library and 10..P10 in the long one:
    t4_sum[10:P4+1, 10:P10+1].sum(), t10_sum[10:P4+1, 10:P10+1].sum()

Usage
    python joint_coverage.py positions --mutation-tables <untrimmed table> <trimmed table> \\
        --libraries 60wk_male_ctx_2_4kb 60wk_male_ctx_2_10kb --out-dir <dir>
    python joint_coverage.py chrom --chrom T2T1 --positions <dir>/call_positions.tsv \\
        --cntx-dir <results>/variant_calling --mask excluded_regions.bed \\
        --fai b6_t2t_mm39XYM.fa.fai --libraries 60wk_male_ctx_2_4kb 60wk_male_ctx_2_10kb --out-dir <dir>
    python joint_coverage.py merge --libraries 60wk_male_ctx_2_4kb 60wk_male_ctx_2_10kb --out-dir <dir>
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

# Inputs, set from the command line
CNTX = MASK_BED = FAI = None
LIBS = []                         # [short library, long library]
CAP = 64
CHUNK = 10_000_000


def chrom_len(chrom):
    for line in open(FAI):
        f = line.split("\t")
        if f[0] == chrom:
            return int(f[1])
    sys.exit(f"{chrom} not in {FAI}")


def dense_mask(chrom, L):
    d = pd.read_csv(MASK_BED, sep="\t", header=None, usecols=[0, 1, 2],
                    names=["chrom", "start", "end"], dtype={"chrom": str})
    d = d[d.chrom == chrom]
    edge = np.zeros(L + 1, np.int32)
    np.add.at(edge, d.start.to_numpy(np.int64), 1)
    np.add.at(edge, np.minimum(d.end.to_numpy(np.int64), L), -1)
    return np.cumsum(edge[:L]) > 0


def dense_cov(lib, chrom, ext, L):
    path = f"{CNTX}/ext{ext}/{lib}/{chrom}.ext{ext}.qual40.cntx.bed.gz"
    cov = np.zeros(L, np.int16)
    for ch in pd.read_csv(path, sep="\t", header=None, usecols=[1, 4], names=["start", "ic"],
                          dtype={"start": np.int64, "ic": np.int32}, chunksize=CHUNK, compression="gzip"):
        cov[ch.start.to_numpy()] = np.minimum(ch.ic.to_numpy(), 32767)
    print(f"  {lib} ext{ext} {chrom} read", flush=True)
    return cov


def run_positions(a):
    parts = []
    for t in a.mutation_tables:
        d = pd.read_csv(t, sep="\t", usecols=["CHROM", "POS", "LIBRARY"], dtype={"CHROM": str})
        parts.append(d[d.LIBRARY.isin(LIBS)][["CHROM", "POS"]])
    p = pd.concat(parts).drop_duplicates().sort_values(["CHROM", "POS"])
    os.makedirs(a.out_dir, exist_ok=True)
    p.rename(columns={"CHROM": "chrom", "POS": "pos_1based"}).to_csv(
        f"{a.out_dir}/call_positions.tsv", sep="\t", index=False)
    print(f"{len(p):,} call positions -> {a.out_dir}/call_positions.tsv")


def run_chrom(a):
    chrom = a.chrom
    L = chrom_len(chrom)
    masked = dense_mask(chrom, L)
    t4, t10 = dense_cov(LIBS[0], chrom, 450, L), dense_cov(LIBS[1], chrom, 450, L)
    u4, u10 = dense_cov(LIBS[0], chrom, 0, L), dense_cov(LIBS[1], chrom, 0, L)

    sel = ~masked & ((u4 > 0) | (u10 > 0))
    idx = np.minimum(u4[sel], CAP).astype(np.int64) * (CAP + 1) + np.minimum(u10[sel], CAP).astype(np.int64)
    size = (CAP + 1) ** 2
    n_pos = np.bincount(idx, minlength=size).reshape(CAP + 1, CAP + 1)
    t4_sum = np.bincount(idx, weights=t4[sel].astype(np.float64), minlength=size).reshape(CAP + 1, CAP + 1)
    t10_sum = np.bincount(idx, weights=t10[sel].astype(np.float64), minlength=size).reshape(CAP + 1, CAP + 1)
    n_bad = int(((t4 > u4) | (t10 > u10)).sum())

    pos = pd.read_csv(a.positions, sep="\t", dtype={"chrom": str})
    pos = pos[pos.chrom == chrom]
    p0 = pos.pos_1based.to_numpy(np.int64) - 1
    os.makedirs(f"{a.out_dir}/parts", exist_ok=True)
    np.savez(f"{a.out_dir}/parts/{chrom}.npz", n_pos=n_pos, t4_sum=np.rint(t4_sum).astype(np.int64),
             t10_sum=np.rint(t10_sum).astype(np.int64), n_trim_above_untrim=n_bad,
             lookup_chrom=pos.chrom.to_numpy(str), lookup_pos=p0 + 1,
             lookup_t4=t4[p0], lookup_t10=t10[p0], lookup_u4=u4[p0], lookup_u10=u10[p0], lookup_masked=masked[p0])
    print(f"{chrom}: {int(n_pos.sum()):,} unmasked covered positions; trimmed > untrimmed at {n_bad}", flush=True)


def run_merge(a):
    parts = sorted(glob.glob(f"{a.out_dir}/parts/T2T*.npz"))
    if len(parts) != 19:
        sys.exit(f"expected 19 parts, found {len(parts)}")
    z = [np.load(p) for p in parts]
    out = dict(libs=np.array(LIBS), cap=CAP, n_trim_above_untrim=sum(int(x["n_trim_above_untrim"]) for x in z))
    for k in ["n_pos", "t4_sum", "t10_sum"]:
        out[k] = sum(x[k] for x in z)
    for k in ["chrom", "pos", "t4", "t10", "u4", "u10", "masked"]:
        out[f"lookup_{k}"] = np.concatenate([x[f"lookup_{k}"] for x in z])
    np.savez(f"{a.out_dir}/joint_untrimmed_coverage.npz", **out)
    for p in parts:
        os.remove(p)
    os.rmdir(f"{a.out_dir}/parts")
    print(f"merged {len(parts)} parts -> {a.out_dir}/joint_untrimmed_coverage.npz; "
          f"positions with trimmed > untrimmed coverage: {out['n_trim_above_untrim']}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["positions", "chrom", "merge"])
    ap.add_argument("--libraries", nargs=2, required=True, metavar=("SHORT", "LONG"))
    ap.add_argument("--mutation-tables", nargs="+", help="positions: the mutation tables of mutation_tables.py")
    ap.add_argument("--chrom")
    ap.add_argument("--positions", help="chrom: TSV with columns chrom, pos_1based")
    ap.add_argument("--cntx-dir", help="chrom: parent of ext0/ and ext450/ of variant_calling.smk")
    ap.add_argument("--mask", help="chrom: masked sites BED")
    ap.add_argument("--fai", help="chrom: reference .fai")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    LIBS = list(a.libraries)
    CNTX, MASK_BED, FAI = a.cntx_dir, a.mask, a.fai
    {"positions": run_positions, "chrom": run_chrom, "merge": run_merge}[a.mode](a)
