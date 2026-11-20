#!/usr/bin/env python3
"""
sziraki_pseudobulk.py

Pseudo-bulk of the Sziraki et al. 2023 mouse-brain atlas (GEO GSE212606) per
(cell type x age), for the accessible chromatin of accessibility_regions.py
(Methods, Chromatin accessibility). The pipeline runs it with --only atac.

Each MatrixMarket file is read once, in chunks:
  RNA   GSM6538356_RNA_gene_count.txt.gz    genes x cells, UMI counts
        -> per (gene, group): sum of UMIs, sum of per-cell CP10k
           (UMIs / cell total * 1e4), number of cells with >= 1 UMI
  ATAC  GSM6538357_ATAC_peak_count.txt.gz   peaks x cells, read counts
        -> per (peak, group): number of cells with >= 1 read

A group is (Main_cluster_name, age). Wild-type cells only: RNA `Type` 3mo,
6mo, 21mo; ATAC `Conditions` prefix Young, WT, Aged. The matrix columns are
checked against the *_cell_annotation.csv rows: the column sums equal UMI_count
(RNA) or reads_in_peaks (ATAC) for every cell.

Output (--out)
  sziraki_rna_pseudobulk.npz
      gene_id, gene_type, gene_name          [n_genes]
      group, group_celltype, group_age       [G]
      n_cells, umi_total                     [G]
      umi_sum    int64   [n_genes, G]
      cp10k_sum  float64 [n_genes, G]   mean CP10k of a pool = sum / n_cells
      n_expr     int32   [n_genes, G]
  sziraki_atac_pseudobulk.npz
      peak ("chr-start-end", as published)   [n_peaks]
      group, group_celltype, group_age, n_cells, reads_total
      n_nonzero  int32   [n_peaks, G]   cells with >= 1 read in the peak

Usage:  python sziraki_pseudobulk.py --indir <dir with the GEO files> --out <dir> --only atac
"""
import argparse
import gzip
import os
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd

CHUNK = 20_000_000
RNA_WT_TYPES = ["3mo", "6mo", "21mo"]
ATAC_WT_PREFIX = ["Young", "WT", "Aged"]


def mtx_header(path):
    """Return (n_rows, n_cols, nnz, n_header_lines) of a gzipped MatrixMarket file."""
    n_hdr = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            n_hdr += 1
            if line.startswith("%"):
                continue
            r, c, nnz = (int(x) for x in line.split())
            return r, c, nnz, n_hdr
    sys.exit(f"[error] no dimension line in {path}")


def stream_mtx(path, n_hdr, chunk=CHUNK):
    """Yield (i, j, v) int32 chunks, 1-based i/j as in the file. Uses pigz when
    it is installed, else Python's gzip."""
    pigz = shutil.which("pigz")
    if pigz:
        proc = subprocess.Popen([pigz, "-dc", path], stdout=subprocess.PIPE,
                                bufsize=1 << 24)
        src = proc.stdout
    else:
        proc, src = None, path
    reader = pd.read_csv(src, sep=" ", header=None, names=["i", "j", "v"],
                         skiprows=n_hdr, dtype=np.int32, engine="c",
                         chunksize=chunk,
                         compression=None if proc else "gzip")
    for ch in reader:
        yield ch["i"].to_numpy(), ch["j"].to_numpy(), ch["v"].to_numpy()
    if proc:
        proc.stdout.close()
        if proc.wait() != 0:
            sys.exit(f"[error] pigz failed on {path}")


def group_index(celltype, age):
    """Cells -> group ids. Groups are sorted (celltype, age) pairs; excluded
    cells (age None) get -1."""
    keep = age.notna().to_numpy()
    pairs = sorted(set(zip(celltype[keep], age[keep])))
    gid = {p: k for k, p in enumerate(pairs)}
    idx = np.full(len(celltype), -1, dtype=np.int64)
    idx[keep] = [gid[p] for p in zip(celltype[keep], age[keep])]
    return idx, pairs


def run_rna(indir, out):
    t0 = time.time()
    ann = pd.read_csv(os.path.join(indir, "GSM6538356_RNA_cell_annotation.csv.gz"),
                      usecols=["UMI_count", "Type", "Main_cluster_name"])
    genes = pd.read_csv(os.path.join(indir, "GSM6538356_RNA_gene_annotation.csv.gz"))
    path = os.path.join(indir, "GSM6538356_RNA_gene_count.txt.gz")
    n_g, n_c, nnz, n_hdr = mtx_header(path)
    if n_g != len(genes) or n_c != len(ann):
        sys.exit(f"[error] RNA dims {n_g}x{n_c} vs annotations {len(genes)}x{len(ann)}")

    age = ann.Type.where(ann.Type.isin(RNA_WT_TYPES))
    grp, pairs = group_index(ann.Main_cluster_name.to_numpy(), age)
    G = len(pairs)
    # 1-based column lookup tables
    grp1 = np.concatenate([[-1], grp])
    inv_umi1 = np.concatenate([[0.0], 1e4 / ann.UMI_count.to_numpy(float)])
    sys.stderr.write(f"RNA: {n_g:,} genes x {n_c:,} cells, {nnz:,} nonzeros; "
                     f"{int((grp >= 0).sum()):,} WT cells in {G} groups\n")

    umi = np.zeros(n_g * G, dtype=np.float64)       # exact for integers < 2^53
    cpk = np.zeros(n_g * G, dtype=np.float64)
    nex = np.zeros(n_g * G, dtype=np.float64)
    colsum = np.zeros(n_c + 1, dtype=np.float64)
    seen = 0
    for i, j, v in stream_mtx(path, n_hdr):
        seen += len(v)
        colsum += np.bincount(j, weights=v, minlength=n_c + 1)
        g = grp1[j]
        k = g >= 0
        flat = (i[k].astype(np.int64) - 1) * G + g[k]
        vk = v[k].astype(np.float64)
        umi += np.bincount(flat, weights=vk, minlength=n_g * G)
        cpk += np.bincount(flat, weights=vk * inv_umi1[j[k]], minlength=n_g * G)
        nex += np.bincount(flat, minlength=n_g * G)
        sys.stderr.write(f"\r  RNA {seen:,}/{nnz:,} ({100*seen/nnz:.0f}%) "
                         f"{time.time()-t0:.0f}s")
    sys.stderr.write("\n")
    if seen != nnz:
        sys.exit(f"[error] RNA read {seen:,} entries, header says {nnz:,}")
    bad = np.flatnonzero(colsum[1:] != ann.UMI_count.to_numpy(float))
    if len(bad):
        sys.exit(f"[error] RNA column sums != UMI_count for {len(bad):,} cells "
                 f"(first: col {bad[0]+1}) -> matrix columns are NOT in "
                 f"annotation order")
    sys.stderr.write("  column sums == UMI_count for every cell: column order OK\n")

    n_cells = np.bincount(grp[grp >= 0], minlength=G)
    umi_tot = np.bincount(grp[grp >= 0],
                          weights=ann.UMI_count.to_numpy(float)[grp >= 0], minlength=G)
    np.savez_compressed(
        os.path.join(out, "sziraki_rna_pseudobulk.npz"),
        gene_id=genes.gene_id.to_numpy(str), gene_type=genes.gene_type.to_numpy(str),
        gene_name=genes.gene_short_name.to_numpy(str),
        group=np.array([f"{c}|{a}" for c, a in pairs]),
        group_celltype=np.array([c for c, _ in pairs]),
        group_age=np.array([a for _, a in pairs]),
        n_cells=n_cells.astype(np.int64), umi_total=umi_tot.astype(np.int64),
        umi_sum=umi.reshape(n_g, G).astype(np.int64),
        cp10k_sum=cpk.reshape(n_g, G),
        n_expr=nex.reshape(n_g, G).astype(np.int32))
    sys.stderr.write(f"  wrote sziraki_rna_pseudobulk.npz ({time.time()-t0:.0f}s)\n")


def run_atac(indir, out):
    t0 = time.time()
    ann = pd.read_csv(os.path.join(indir, "GSM6538357_ATAC_cell_annotation.csv.gz"),
                      usecols=["reads_in_peaks", "Main_cluster_name", "Conditions"])
    peaks = pd.read_csv(os.path.join(indir, "GSM6538357_ATAC_peak_annotation.csv.gz"))
    path = os.path.join(indir, "GSM6538357_ATAC_peak_count.txt.gz")
    n_p, n_c, nnz, n_hdr = mtx_header(path)
    if n_p != len(peaks) or n_c != len(ann):
        sys.exit(f"[error] ATAC dims {n_p}x{n_c} vs annotations {len(peaks)}x{len(ann)}")

    prefix = ann.Conditions.str.split("_").str[0]
    age = prefix.where(prefix.isin(ATAC_WT_PREFIX))
    grp, pairs = group_index(ann.Main_cluster_name.to_numpy(), age)
    G = len(pairs)
    grp1 = np.concatenate([[-1], grp])
    sys.stderr.write(f"ATAC: {n_p:,} peaks x {n_c:,} cells, {nnz:,} nonzeros; "
                     f"{int((grp >= 0).sum()):,} WT cells in {G} groups\n")

    nz = np.zeros(n_p * G, dtype=np.float64)
    colsum = np.zeros(n_c + 1, dtype=np.float64)
    seen = 0
    for i, j, v in stream_mtx(path, n_hdr):
        seen += len(v)
        colsum += np.bincount(j, weights=v, minlength=n_c + 1)
        g = grp1[j]
        k = (g >= 0) & (v > 0)
        flat = (i[k].astype(np.int64) - 1) * G + g[k]
        nz += np.bincount(flat, minlength=n_p * G)
        sys.stderr.write(f"\r  ATAC {seen:,}/{nnz:,} ({100*seen/nnz:.0f}%) "
                         f"{time.time()-t0:.0f}s")
    sys.stderr.write("\n")
    if seen != nnz:
        sys.exit(f"[error] ATAC read {seen:,} entries, header says {nnz:,}")
    bad = np.flatnonzero(colsum[1:] != ann.reads_in_peaks.to_numpy(float))
    if len(bad):
        sys.exit(f"[error] ATAC column sums != reads_in_peaks for {len(bad):,} "
                 f"cells (first: col {bad[0]+1}) -> matrix columns are NOT in "
                 f"annotation order")
    sys.stderr.write("  column sums == reads_in_peaks for every cell: column order OK\n")

    n_cells = np.bincount(grp[grp >= 0], minlength=G)
    rd_tot = np.bincount(grp[grp >= 0],
                         weights=ann.reads_in_peaks.to_numpy(float)[grp >= 0],
                         minlength=G)
    np.savez_compressed(
        os.path.join(out, "sziraki_atac_pseudobulk.npz"),
        peak=peaks.peak.to_numpy(str),
        group=np.array([f"{c}|{a}" for c, a in pairs]),
        group_celltype=np.array([c for c, _ in pairs]),
        group_age=np.array([a for _, a in pairs]),
        n_cells=n_cells.astype(np.int64), reads_total=rd_tot.astype(np.int64),
        n_nonzero=nz.reshape(n_p, G).astype(np.int32))
    sys.stderr.write(f"  wrote sziraki_atac_pseudobulk.npz ({time.time()-t0:.0f}s)\n")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--indir", required=True, help="directory holding the GEO files")
    p.add_argument("--out", required=True)
    p.add_argument("--only", choices=["rna", "atac"], default=None)
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.only in (None, "rna"):
        run_rna(a.indir, a.out)
    if a.only in (None, "atac"):
        run_atac(a.indir, a.out)


if __name__ == "__main__":
    main()
