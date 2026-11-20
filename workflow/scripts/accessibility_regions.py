#!/usr/bin/env python3
"""
accessibility_regions.py

Accessible and inaccessible chromatin of cortical neurons on the T2T
autosomes (Methods, Chromatin accessibility), following Gronska-Peski et al.
2026:

    accessible    snATAC-seq peaks of Sziraki et al. 2023 in which >= --min-frac
                  (2%) of the pooled cells have at least one read
    inaccessible  every other site with a counterpart in GRCm39

The pooled cells are the groups --celltypes x --ages of sziraki_pseudobulk.py
(cortical projection neurons 1-3 and interneurons 1-2, wild type); a peak's
fraction is the cells with >= 1 read in it over the cells of those groups.

The peaks are on GRCm39. A T2T site is accessible if its GRCm39 counterpart,
through the liftOver chain from the T2T assembly to GRCm39 (--chain, UCSC
GenArk), lies in an accessible peak. T2T sites with no counterpart in GRCm39
are left out. Chain contig names (GenBank accessions) are matched to
T2T1..T2T19 by sequence length (--fai).

--validate-fasta compares the T2T and GRCm39 sequence of a random sample of
the lifted accessible intervals.

Output (--out-dir)
  accessibility_regions.bed.gz        chrom start end region (0-based, half-open;
                                      region = accessible | inaccessible)
  accessible_peaks_mm39.bed.gz        the accessible peaks, merged, on GRCm39
  accessibility_regions_summary.tsv   peak counts, bp per region, lift yield

Usage
  python accessibility_regions.py --atac-npz sziraki_atac_pseudobulk.npz \\
      --chain GCA_964188535.1ToMm39.over.chain.gz --fai b6_t2t_mm39XYM.fa.fai \\
      --chroms T2T1 ... T2T19 --validate-fasta b6_t2t_mm39XYM.fa <GRCm39 FASTA> \\
      --out-dir <results>/chromatin_accessibility/regions
"""
import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd

DEFAULT_CELLTYPES = ["Cortical projection neurons 1", "Cortical projection neurons 2",
                     "Cortical projection neurons 3", "Interneurons 1", "Interneurons 2"]
DEFAULT_AGES = ["Young", "WT", "Aged"]
REGIONS = ["inaccessible", "accessible"]                  # labels 1, 2


def merge(s, e):
    """Sort and merge [s, e) intervals; returns merged numpy arrays."""
    if len(s) == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    o = np.lexsort((e, s))
    s, e = np.asarray(s, np.int64)[o], np.asarray(e, np.int64)[o]
    ce = np.maximum.accumulate(e)
    brk = np.flatnonzero(s[1:] > ce[:-1]) + 1
    st = np.concatenate([[0], brk])
    en = np.concatenate([brk, [len(s)]]) - 1
    return s[st], ce[en]


def accessible_peaks(npz, celltypes, ages, min_frac):
    z = np.load(npz)
    for want, have, what in ((celltypes, z["group_celltype"], "cell type"),
                             (ages, z["group_age"], "age")):
        missing = sorted(set(want) - set(have))
        if missing:
            sys.exit(f"[error] {what}(s) not in the reference: {missing}")
    sel = np.isin(z["group_celltype"], celltypes) & np.isin(z["group_age"], ages)
    n_cells = int(z["n_cells"][sel].sum())
    frac = z["n_nonzero"][:, sel].sum(1) / n_cells
    pk = pd.Series(z["peak"]).str.split("-", expand=True)
    peaks = pd.DataFrame({"chrom": pk[0], "start": pk[1].astype(np.int64),
                          "end": pk[2].astype(np.int64), "frac": frac})
    groups = [f"{g} (n={n:,})" for g, n in zip(z["group"][sel], z["n_cells"][sel])]
    return peaks, n_cells, groups


def read_chain(path):
    """Yield (tName, tSize, qName, qSize, qStrand, blocks) with blocks an int64
    array of rows (tStart, qStart_on_qStrand, size)."""
    with gzip.open(path, "rt") as fh:
        hdr, blocks = None, []
        for line in fh:
            f = line.split()
            if not f:
                continue
            if f[0] == "chain":
                if hdr:
                    yield (*hdr, np.array(blocks, dtype=np.int64))
                tName, tSize, tStrand, tStart = f[2], int(f[3]), f[4], int(f[5])
                qName, qSize, qStrand, qStart = f[7], int(f[8]), f[9], int(f[10])
                if tStrand != "+":
                    sys.exit("[error] chain with target strand '-' (not a liftOver chain)")
                hdr = (tName, tSize, qName, qSize, qStrand)
                t, q, blocks = tStart, qStart, []
                continue
            size = int(f[0])
            blocks.append((t, q, size))
            if len(f) == 3:
                t += size + int(f[1])
                q += size + int(f[2])
        if hdr:
            yield (*hdr, np.array(blocks, dtype=np.int64))


def lift(chain_path, tmap, acc):
    """Map accessible mm39 intervals onto T2T through the chain's blocks.
    Returns {t2t_chrom: (image_s, image_e, acc_s, acc_e, blocks_bp)} and the
    per-block table used for validation."""
    img = {c: ([], []) for c in tmap.values()}
    hit = {c: ([], []) for c in tmap.values()}
    val = []                                    # (tchrom, t0, t1, qname, q0, q1, strand)
    for tName, tSize, qName, qSize, qStrand, b in read_chain(chain_path):
        tc = tmap.get(tName)
        if tc is None or len(b) == 0:
            continue
        ts, qs, sz = b[:, 0], b[:, 1], b[:, 2]
        img[tc][0].append(ts)
        img[tc][1].append(ts + sz)
        # forward-strand mm39 interval of each block
        if qStrand == "+":
            qa, qb = qs, qs + sz
        else:
            qa, qb = qSize - (qs + sz), qSize - qs
        ps, pe = acc.get(qName, (None, None))
        if ps is None or len(ps) == 0:
            continue
        i0 = np.searchsorted(pe, qa, side="right")      # first peak ending after qa
        i1 = np.searchsorted(ps, qb, side="left")       # peaks starting before qb
        n = np.maximum(i1 - i0, 0)
        if not n.any():
            continue
        bi = np.repeat(np.arange(len(b)), n)
        pi = np.concatenate([np.arange(a, a + k) for a, k in zip(i0[n > 0], n[n > 0])])
        x = np.maximum(qa[bi], ps[pi])
        y = np.minimum(qb[bi], pe[pi])
        if qStrand == "+":
            t0 = ts[bi] + (x - qa[bi])
            t1 = ts[bi] + (y - qa[bi])
        else:
            t0 = ts[bi] + (qb[bi] - y)
            t1 = ts[bi] + (qb[bi] - x)
        hit[tc][0].append(t0)
        hit[tc][1].append(t1)
        val.append(pd.DataFrame({"tchrom": tc, "t0": t0, "t1": t1, "qname": qName,
                                 "q0": x, "q1": y, "strand": qStrand}))
    out = {}
    for c in tmap.values():
        s = np.concatenate(img[c][0]) if img[c][0] else np.zeros(0, np.int64)
        e = np.concatenate(img[c][1]) if img[c][1] else np.zeros(0, np.int64)
        raw_bp = int((e - s).sum())
        ms, me = merge(s, e)
        if int((me - ms).sum()) != raw_bp:
            sys.exit(f"[error] chain blocks overlap on {c}: a T2T base would have "
                     f"two mm39 images. Not a single-coverage liftOver chain.")
        a_s = np.concatenate(hit[c][0]) if hit[c][0] else np.zeros(0, np.int64)
        a_e = np.concatenate(hit[c][1]) if hit[c][1] else np.zeros(0, np.int64)
        out[c] = (ms, me) + merge(a_s, a_e)
    return out, (pd.concat(val, ignore_index=True) if val else pd.DataFrame())


def validate(val, t2t_fa, mm39_fa, n=2000, seed=0):
    import pysam
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    ft, fm = pysam.FastaFile(t2t_fa), pysam.FastaFile(mm39_fa)
    mnames = set(fm.references)
    v = val[(val.t1 - val.t0) >= 50]
    v = v.sample(min(n, len(v)), random_state=seed)
    same = tot = 0
    for r in v.itertuples():
        qn = r.qname if r.qname in mnames else r.qname.replace("chr", "")
        a = ft.fetch(r.tchrom, int(r.t0), int(r.t1)).upper()
        b = fm.fetch(qn, int(r.q0), int(r.q1)).upper()
        if r.strand == "-":
            b = b.translate(comp)[::-1]
        same += sum(x == y for x, y in zip(a, b))
        tot += len(a)
    return len(v), same / max(tot, 1)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--atac-npz", required=True)
    p.add_argument("--chain", required=True, help="GCA_964188535.1ToMm39.over.chain.gz")
    p.add_argument("--fai", required=True, help="reference .fai (T2T1..T2T19)")
    p.add_argument("--chroms", nargs="+", required=True)
    p.add_argument("--celltypes", nargs="+", default=DEFAULT_CELLTYPES)
    p.add_argument("--ages", nargs="+", default=DEFAULT_AGES)
    p.add_argument("--min-frac", type=float, default=0.02)
    p.add_argument("--validate-fasta", nargs=2, metavar=("T2T_FA", "MM39_FA"), default=None)
    p.add_argument("--out-dir", required=True)
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    # chain target names -> T2T names by length
    h = pd.read_csv(a.fai, sep="\t", header=None, usecols=[0, 1], names=["n", "L"])
    h = h[h.n.isin(a.chroms)]
    tsz = {}
    for tName, tSize, *_ in read_chain(a.chain):
        tsz[tName] = tSize
    by_len = {}
    for n, L in tsz.items():
        by_len.setdefault(L, []).append(n)
    tmap = {}
    for n, L in zip(h.n, h.L):
        cands = by_len.get(L, [])
        if len(cands) != 1:
            sys.exit(f"[error] {n} (length {L:,}) matches {len(cands)} chain targets")
        tmap[cands[0]] = n
    clen = dict(zip(h.n, h.L))

    peaks, n_cells, groups = accessible_peaks(a.atac_npz, a.celltypes, a.ages, a.min_frac)
    sys.stderr.write(f"reference: {n_cells:,} cells pooled from {len(groups)} groups\n"
                     + "".join(f"    {x}\n" for x in groups))
    accp = peaks[peaks.frac >= a.min_frac]
    acc = {c: merge(g.start.to_numpy(), g.end.to_numpy()) for c, g in accp.groupby("chrom")}
    acc_mm39_bp = sum(int((e - s).sum()) for s, e in acc.values())
    pd.concat([pd.DataFrame({"chrom": c, "start": s, "end": e}) for c, (s, e) in acc.items()]
              ).sort_values(["chrom", "start"]).to_csv(
        os.path.join(a.out_dir, "accessible_peaks_mm39.bed.gz"), sep="\t", header=False,
        index=False, compression="gzip")
    auto_mm39 = [f"chr{i}" for i in range(1, 20)]
    acc_auto_bp = sum(int((e - s).sum()) for c, (s, e) in acc.items() if c in auto_mm39)
    sys.stderr.write(f"peaks: {len(peaks):,} total, {len(accp):,} with >= "
                     f"{100*a.min_frac:g}% of cells -> {acc_mm39_bp:,} bp merged "
                     f"(autosomes {acc_auto_bp:,})\n")

    lifted, val = lift(a.chain, tmap, acc)
    bed, tot = [], {"accessible": 0, "inaccessible": 0, "no_mm39_image": 0}
    for c in a.chroms:
        L = clen[c]
        ims, ime, acs, ace = lifted[c]
        lab = np.zeros(L, dtype=np.uint8)
        for s, e in zip(ims, ime):
            lab[s:e] = 1
        for s, e in zip(acs, ace):
            lab[s:e] = 2
        cut = np.flatnonzero(np.diff(lab)) + 1
        s = np.concatenate([[0], cut])
        e = np.concatenate([cut, [L]])
        v = lab[s]
        k = v > 0
        bed.append(pd.DataFrame({"chrom": c, "start": s[k], "end": e[k],
                                 "region": np.array(REGIONS)[v[k] - 1]}))
        cnt = np.bincount(lab, minlength=3)
        tot["no_mm39_image"] += int(cnt[0])
        tot["inaccessible"] += int(cnt[1])
        tot["accessible"] += int(cnt[2])
    bed = pd.concat(bed, ignore_index=True)
    bed.to_csv(os.path.join(a.out_dir, "accessibility_regions.bed.gz"), sep="\t",
               header=False, index=False, compression="gzip")

    auto_bp = sum(clen[c] for c in a.chroms)
    rows = [dict(item="peaks_total", value=len(peaks)),
            dict(item=f"peaks_accessible(frac>={a.min_frac:g})", value=len(accp)),
            dict(item="accessible_bp_mm39_all_chroms", value=acc_mm39_bp),
            dict(item="accessible_bp_mm39_autosomes", value=acc_auto_bp)]
    for k, v in tot.items():
        rows.append(dict(item=f"t2t_{k}_bp", value=v, pct_autosomes=round(100 * v / auto_bp, 3)))
    rows.append(dict(item="lift_yield(t2t_acc/mm39_acc_autosomes)",
                     value=round(tot["accessible"] / max(acc_auto_bp, 1), 5)))
    if a.validate_fasta:
        nv, ident = validate(val, *a.validate_fasta)
        rows.append(dict(item=f"validation_identity({nv}_intervals)", value=round(ident, 6)))
        sys.stderr.write(f"validation: {nv} lifted intervals, T2T vs mm39 identity {ident:.4%}\n")
    rows.append(dict(item=f"reference: {n_cells} cells; types={'|'.join(a.celltypes)}; "
                          f"ages={'|'.join(a.ages)}; min_frac={a.min_frac:g}"))
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(a.out_dir, "accessibility_regions_summary.tsv"), sep="\t", index=False)
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
