#!/usr/bin/env bash
# =============================================================================
# build_reference.sh — reference genome and its pbmm2 index (Methods, Alignment)
#
# Contigs of the reference:
#   T2T1..T2T19              the 19 autosomes of the C57BL/6J T2T assembly
#                            (GCA_964188535.1), renamed; variants are called here only
#   CAXLPS0100000xx.1 (x218) the unplaced scaffolds of the T2T assembly (decoys)
#   X, Y, MT                 from GRCm39 (decoys)
#
# Usage:  bash build_reference.sh <T2T FASTA> <T2T assembly report> <GRCm39 FASTA> <output FASTA>
#         FASTA files uncompressed; samtools and pbmm2 on PATH.
# Writes  <output FASTA> and its .fai, <output FASTA without .fa>.mmi (pbmm2,
#         preset CCS), and a .fai next to each input FASTA.
# =============================================================================

set -euo pipefail

# ── Paths / params ───────────────────────────────────────────────────────────
T2T_FASTA=$(realpath "$1")
REPORT=$(realpath "$2")
MM39=$(realpath "$3")
HYBRID=$(realpath -m "$4")
OUT_DIR=$(dirname "$HYBRID")
PRESET="CCS"          # as the --preset of pbmm2 align (alignment.smk)

command -v samtools >/dev/null || { echo "ERROR: samtools not found" >&2; exit 1; }
command -v pbmm2    >/dev/null || { echo "ERROR: pbmm2 not found" >&2; exit 1; }
[[ -s "$MM39"      ]] || { echo "ERROR: GRCm39 FASTA not found: $MM39" >&2; exit 1; }
[[ -s "$T2T_FASTA" ]] || { echo "ERROR: T2T FASTA not found: $T2T_FASTA" >&2; exit 1; }

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# temporary folder, removed on exit
WORK=$(mktemp -d "${OUT_DIR}/.build_hybrid.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

# ── 0. Contig names: autosomes T2T1..T2T19; X, Y and MT of GRCm39 ────────────
samtools faidx "$MM39"
if cut -f1 "$MM39.fai" | grep -qx 'chr1'; then PFX="chr"; else PFX=""; fi
PFX="T2T${PFX}"       # -> T2T1..T2T19 (Ensembl GRCm39 names)
mapfile -t DECOY_CONTIGS < <(cut -f1 "$MM39.fai" | grep -Eix "(chr)?(X|Y|M|MT)")
echo ">>> autosome naming: '${PFX}N'   borrowed decoys: ${DECOY_CONTIGS[*]}"
[[ ${#DECOY_CONTIGS[@]} -ge 3 ]] || { echo "ERROR: expected X, Y, M(T) in mm39; found: ${DECOY_CONTIGS[*]}" >&2; exit 1; }

# ── 1. Index the T2T FASTA ───────────────────────────────────────────────────
echo ">>> [1/5] Indexing source T2T FASTA ..."
samtools faidx "$T2T_FASTA"

# ── 2. Autosome names and unplaced scaffolds, from the assembly report ───────
echo ">>> [2/5] Building autosome rename map + contig list from report ..."
if [[ -s "$REPORT" ]]; then
    # autosomes: GenBank accession (column 5) -> T2T<N>
    awk -F'\t' -v pfx="$PFX" '!/^#/ && $2=="assembled-molecule"{
        n=$1; sub(/.*#chr/,"",n); print $5"\t"pfx n }' "$REPORT" \
        | sort -k2 -V > "$WORK/autosome_rename.tsv"
    awk -F'\t' '!/^#/ && $2=="unplaced-scaffold"{print $5}' "$REPORT" > "$WORK/unplaced_accessions.txt"
    # autosome length from the report, checked in step 4
    EXP_AUTO_BP=$(awk -F'\t' '!/^#/ && $2=="assembled-molecule"{s+=$9} END{print s}' "$REPORT")
else
    echo "    (assembly_report not found; falling back to defline parsing)" >&2
    cut -f1 "$T2T_FASTA.fai" > "$WORK/all_names.txt"
    grep '^>' "$T2T_FASTA" \
      | sed -n "s/^>\([^[:space:]]*\).*chromosome:\{0,1\} \([0-9][0-9]*\).*/\1\t${PFX}\2/p" \
      | sort -k2 -V > "$WORK/autosome_rename.tsv"
    grep '^>' "$T2T_FASTA" | grep -v 'chromosome:' | sed 's/^>//; s/ .*//' > "$WORK/unplaced_accessions.txt"
    EXP_AUTO_BP=""
fi

n_auto=$(wc -l < "$WORK/autosome_rename.tsv")
n_unpl=$(wc -l < "$WORK/unplaced_accessions.txt")
echo "    autosomes: $n_auto   unplaced contigs (decoys): $n_unpl"
cat "$WORK/autosome_rename.tsv"
[[ "$n_auto" -eq 19 ]] || { echo "ERROR: expected 19 autosomes, found $n_auto." >&2; exit 1; }
[[ "$n_unpl" -ge 1  ]] || { echo "ERROR: no unplaced contigs found." >&2; exit 1; }

# ── 3. Autosomes, renamed ────────────────────────────────────────────────────
echo ">>> [3/5] Extracting + renaming autosomes ..."
: > "$WORK/t2t_autosomes.fa"
while IFS=$'\t' read -r acc chrom; do
    samtools faidx "$T2T_FASTA" "$acc" | sed "1s/.*/>${chrom}/" >> "$WORK/t2t_autosomes.fa"
done < "$WORK/autosome_rename.tsv"

# ── 3b. Unplaced scaffolds, names unchanged ──────────────────────────────────
echo ">>> [3b/5] Extracting $n_unpl unplaced contigs as decoys ..."
samtools faidx -r "$WORK/unplaced_accessions.txt" "$T2T_FASTA" > "$WORK/t2t_unplaced.fa"

# ── 4. X, Y and MT of GRCm39; concatenate and check ──────────────────────────
echo ">>> [4/5] Borrowing ${DECOY_CONTIGS[*]} from mm39 + concatenating ..."
samtools faidx "$MM39" "${DECOY_CONTIGS[@]}" > "$WORK/mm39_XYM.fa"
cat "$WORK/t2t_autosomes.fa" "$WORK/t2t_unplaced.fa" "$WORK/mm39_XYM.fa" > "$HYBRID.tmp"
mv -f "$HYBRID.tmp" "$HYBRID"
samtools faidx "$HYBRID"

n_total=$(wc -l < "${HYBRID}.fai")
exp=$(( n_auto + n_unpl + ${#DECOY_CONTIGS[@]} ))
echo "    total contigs: $n_total (expect $exp)"
[[ "$n_total" -eq "$exp" ]] || { echo "ERROR: $n_total contigs, expected $exp" >&2; exit 1; }

dups=$(cut -f1 "${HYBRID}.fai" | sort | uniq -d)
[[ -z "$dups" ]] || { echo "ERROR: duplicate contig names: $dups" >&2; exit 1; }

if [[ -n "$EXP_AUTO_BP" ]]; then
    got=$(awk -v pfx="$PFX" '$1 ~ ("^"pfx"[0-9]+$"){s+=$2} END{print s+0}' "${HYBRID}.fai")
    [[ "$got" -eq "$EXP_AUTO_BP" ]] || { echo "ERROR: autosome bp $got != $EXP_AUTO_BP" >&2; exit 1; }
    echo "    autosome bp verified: $got"
fi
echo "    contigs / lengths:"; cut -f1,2 "${HYBRID}.fai"

# ── 5. pbmm2 index ───────────────────────────────────────────────────────────
echo ">>> [5/5] Indexing hybrid with pbmm2 (--preset $PRESET) ..."
pbmm2 index "$HYBRID" "${HYBRID%.fa}.mmi" --preset "$PRESET"

echo "=== done ==="
echo "FASTA: $HYBRID   ($n_total contigs)"
echo "index: ${HYBRID%.fa}.mmi"
