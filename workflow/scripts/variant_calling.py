#!/usr/bin/env python3
"""
variant_calling.py

Variant calling (Methods, Variant calling and filtering): compares the two
strand consensus reads of each molecule with each other and with the
reference, at every site where the bases of both reads pass the base filters
(quality >= -q, set by read_and_base_filters.py). Applies the molecule filters
before comparing.

Input:  a SAM/BAM stream of one chromosome from read_and_base_filters.py (-i -)
Output: stdout       one line per site of a molecule where either read differs
                     from the reference (BED; sort with sort-bed):
                       chrom  start  end  triplet  fwd_info  rev_info  zmw  type
                     fwd_info, rev_info = base:qual:sm:sx:qpos:readlen:divergence:np
                     (qpos from the read's own 5' end)
                     type  d    variant: both reads carry the same non-reference base
                           d_c  both reads differ from the reference and from each other
                           f    mismatch: only the forward read differs
                           r    mismatch: only the reverse read differs
                           m    both reads carry the reference base (with -m only)
        <prefix>.cntx.bed.gz  one line per site covered by at least one molecule:
                       chrom  start  end  total_cov  informative_cov  triplet
                     total_cov        molecules compared at the site
                     informative_cov  molecules whose two reads carry the same base
                                      (the coverage of the paper)
        --qpos-hist  bases by distance from each read end

Molecule filters:
    --max-len-div N           strand length difference above N bases
    --drop-supplementary      chimeric molecule (SA tag on either read)
    --max-strand-nonoverlap N strand alignment non-overlap above N sites

Usage:
    ... | python variant_calling.py -i - -r b6_t2t_mm39XYM.fa -q 40 \
              -p <out prefix> --max-len-div N --drop-supplementary \
              --max-strand-nonoverlap M | sort-bed - | bgzip > <prefix>.muts.bed.gz
"""
import sys, pysam, argparse, os, pickle  # type: ignore
import numpy as np

def safediv(n, d):
    if d > 0:
        return float(n) / float(d)
    else:
        return np.nan

def load_genome(fasta_file_path):
    """Reference sequences as {contig: sequence}, cached in <fasta>.pkl."""
    cache_file = f"{fasta_file_path}.pkl"

    if os.path.exists(cache_file):
        sys.stderr.write(f"Loading cached reference genome from {cache_file}...\n")
        with open(cache_file, "rb") as f:
            try:
                genome = pickle.load(f)
                sys.stderr.write("...genome loading complete (from cache).\n")
                return genome
            except Exception as e:
                sys.stderr.write(f"Warning: Could not load cache file {cache_file}. Re-parsing. Error: {e}\n")

    sys.stderr.write(f"Loading reference genome from {fasta_file_path} into memory...\n")
    genome = {}
    with pysam.FastaFile(fasta_file_path) as ref:
        for chrom in ref.references:
            genome[chrom] = ref.fetch(chrom).upper()
    sys.stderr.write("...genome loading complete (from fasta).\n")

    try:
        sys.stderr.write(f"Saving genome to cache file: {cache_file} ...\n")
        with open(cache_file, "wb") as f:
            pickle.dump(genome, f, protocol=pickle.HIGHEST_PROTOCOL)
        sys.stderr.write("...cache saving complete.\n")
    except Exception as e:
        sys.stderr.write(f"Warning: Could not save cache file. {e}\n")

    return genome

def get_ssdna_info(aln, genome, minqual):
    """Bases of one strand consensus read with quality >= minqual and no N in
    the reference triplet: {ref position: (ref base, base, sm, sx, qual, qpos,
    read length, NM divergence, triplet)}, with the read's mismatch-only
    divergence and its number of passes."""
    info = {}
    chrom = aln.reference_name
    chrom_seq = genome[chrom]
    length = aln.query_length
    aln_length = aln.query_alignment_length
    tot_div = safediv(aln.get_tag("NM"), aln_length)
    tot_div = f"{tot_div:.3f}"
    read, quals = aln.query_sequence, aln.query_qualities
    sms, sxs = aln.get_tag("sm"), aln.get_tag("sx")
    # sm and sx in the order of SEQ for reads aligned to the minus strand
    if aln.is_reverse:
        sms = sms[::-1]
        sxs = sxs[::-1]

    try:
        np_passes = int(aln.get_tag("np"))
    except KeyError:
        np_passes = 0

    n_mis = 0

    for qpos, rpos in aln.get_aligned_pairs(matches_only=True):
        triplet = chrom_seq[rpos - 1 : rpos + 2]
        if "N" in triplet:
            continue
        ref_base = triplet[1]
        base, sm, sx, qual = read[qpos], sms[qpos], sxs[qpos], quals[qpos]
        # Position in the read counted from the read's own 5' end
        native_qpos = (length - 1 - qpos) if aln.is_reverse else qpos
        if quals[qpos] >= minqual:
            info[rpos] = (ref_base, base, sm, sx, qual, native_qpos, length, tot_div, triplet)
        if ref_base != base:
            n_mis += 1

    sbs_div = safediv(n_mis, aln_length)
    sbs_div = f"{sbs_div:.3f}"

    return sbs_div, info, np_passes


def process_pair(key, fwd, rev, genome, args, collect_qpos=False):
    """Compares the two strand consensus reads of one molecule. Returns the
    compared sites, the sites where both reads carry the same base, the output
    lines and the two read lengths; with collect_qpos, also the qpos of each
    informative site on each read."""
    fwd_sbs_div, fwd_mm, fwd_np = get_ssdna_info(fwd, genome, args.minqual)
    rev_sbs_div, rev_mm, rev_np = get_ssdna_info(rev, genome, args.minqual)

    fwd_n_qual = len(fwd_mm)
    rev_n_qual = len(rev_mm)
    n_qual_ratio = safediv(fwd_n_qual, rev_n_qual)
    n_qual_ratio = f"{n_qual_ratio:.2f}"
    chrom = fwd.reference_name

    seq_lines        = []
    common_bases     = fwd_mm.keys() & rev_mm.keys()
    informative_bases = set()   # positions where fbase == rbase (m or d)
    inf_fq, inf_rq   = ([], []) if collect_qpos else (None, None)
    pair_len         = (fwd.query_length, rev.query_length)

    for pos in common_bases:
        fref, fbase, fsm, fsx, fq, fqpos, flen, ftotdiv, ftrip = fwd_mm[pos]
        _,    rbase, rsm, rsx, rq, rqpos, rlen, rtotdiv, rtrip = rev_mm[pos]

        if fbase == rbase:
            informative_bases.add(pos)
            if collect_qpos:
                inf_fq.append(fqpos)
                inf_rq.append(rqpos)

        if fref == fbase == rbase:
            if ftrip != rtrip:
                sys.exit("Something is off with triplets!")
            if not args.mono:
                continue
            mtype = "m"

        elif any(b != fref for b in (fbase, rbase)):
            if fbase != fref and rbase != fref:
                mtype = "d" if rbase == fbase else "d_c"
            elif fbase != fref:
                mtype = "f"
            elif rbase != fref:
                mtype = "r"
            else:
                sys.exit("Unhandled base config")

        seq_lines.append("\t".join([
            chrom,
            str(pos), str(pos+1),
            ftrip,
            f"{fbase}:{fq}:{fsm}:{fsx}:{fqpos}:{flen}:{fwd_sbs_div}:{fwd_np}",
            f"{rbase}:{rq}:{rsm}:{rsx}:{rqpos}:{rlen}:{rev_sbs_div}:{rev_np}",
            f"r{key}",
            mtype
        ]))

    return common_bases, informative_bases, seq_lines, pair_len, inf_fq, inf_rq


def main():
    p = argparse.ArgumentParser(description="Call variants and mismatches between the two strand consensus reads of each molecule")
    p.add_argument("-i", "--bam", required=True, help="Input BAM/SAM file or '-' for stdin")
    p.add_argument("-r", "--ref", required=True, help="Reference FASTA (indexed)")
    p.add_argument("-p", "--prefix", type=str, required=True, help="Prefix of the output files")
    p.add_argument("-q", "--minqual", type=int, default=0, help="Minimum base quality")
    p.add_argument("-b", "--buffer", type=float, help="Progress message every this many reads")
    p.add_argument("-m", "--mono", action="store_true", default=False, help="Also write the sites where both reads carry the reference base")
    p.add_argument("--max-d-per-pair", type=int, default=0,
                    help="Drop molecules with more than this many variants; "
                    "0 = off. (default: 0)")
    # ---- molecule filters ----------------------------------------------------
    p.add_argument("--max-len-div", type=int, default=0,
                    help="Strand length difference filter: drop the molecule if "
                    "|fwd_readlen - rev_readlen| exceeds this many bases "
                    "(strand_length_difference_threshold.py). 0 = off. (default: 0)")
    p.add_argument("--max-len-div-frac", type=float, default=0.0,
                    help="As --max-len-div, as a fraction of the mean read "
                    "length of the molecule. 0 = off. (default: 0)")
    p.add_argument("--drop-supplementary", action="store_true", default=False,
                    help="Chimeric molecule filter: drop the molecule if either "
                    "strand consensus read has a supplementary alignment (SA tag).")
    p.add_argument("--max-strand-nonoverlap", type=int, default=0,
                    help="Strand alignment non-overlap filter: drop the molecule "
                    "if the sites covered by only one of its two alignments exceed "
                    "this number (strand_alignment_nonoverlap_threshold.py). "
                    "0 = off. (default: 0)")
    # ---- bases by distance from the read end --------------------------------
    p.add_argument("--qpos-hist", default=None,
                    help="Write the informative bases by distance from each read "
                    "end (bins of --qpos-bin) to this path.")
    p.add_argument("--qpos-bin", type=int, default=25,
                    help="Bin width in bases for --qpos-hist (default: 25)")
    p.add_argument("--qpos-sample", type=int, default=100,
                    help="Count --qpos-hist from 1 molecule in N; 1 = every "
                    "molecule. (default: 100)")
    args = p.parse_args()

    if args.max_len_div and args.max_len_div_frac:
        sys.exit("[error] use --max-len-div OR --max-len-div-frac, not both")
    if args.qpos_sample < 1:
        sys.exit("[error] --qpos-sample must be >= 1")
    if args.qpos_bin < 1:
        sys.exit("[error] --qpos-bin must be >= 1")

    bam_mode = "r" if args.bam == "-" else "rb"
    bam = pysam.AlignmentFile(sys.stdin if args.bam == "-" else args.bam, bam_mode)
    genome = load_genome(args.ref)
    pairs = {}

    n_pairs = 0                 # molecules that pass the molecule filters
    n_pairs_len_dropped = 0     # dropped by --max-len-div / --max-len-div-frac
    n_pairs_sa_dropped = 0      # dropped by --drop-supplementary
    n_pairs_nonov_dropped = 0   # dropped by --max-strand-nonoverlap
    # informative bases by distance from the 5' and 3' ends of each read, up to 200 kb
    QPOS_NBINS = max(1, 200_000 // max(args.qpos_bin, 1))
    qpos_hist = ({k: np.zeros(QPOS_NBINS, dtype=np.int64)
                  for k in ("fwd_5", "fwd_3", "rev_5", "rev_3")}
                 if args.qpos_hist else None)

    for i, aln in enumerate(bam.fetch(until_eof=True)):

        if i % int(args.buffer) == 0 and i > 0:
            sys.stderr.write(f"{i} reads processed...\n")

        if i == 0:
            chrom = aln.reference_name
            chrom_len = len(genome[chrom])
            coverage         = np.zeros(chrom_len, dtype=np.uint32)  # m + d + d_c + f + r
            informative_cov  = np.zeros(chrom_len, dtype=np.uint32)  # m + d only

        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue

        # Read name <movie>/<zmw>/ccs/<fwd|rev>: the two reads of a molecule share the ZMW
        qname = aln.query_name
        elems = qname.split("/")
        key, strand = elems[1], elems[-1]

        if key not in pairs:
            pairs[key] = {strand: aln}
            continue

        pairs[key][strand] = aln

        # both strand consensus reads of the molecule
        if len(pairs[key]) == 2:
            fwd, rev = pairs[key].get("fwd"), pairs[key].get("rev")
            if fwd is None or rev is None:
                del pairs[key]
                continue

            # ---- molecule filter: strand length difference -----------------
            if args.max_len_div or args.max_len_div_frac:
                flen_q, rlen_q = fwd.query_length, rev.query_length
                dlen = abs(flen_q - rlen_q)
                too_far = (dlen > args.max_len_div) if args.max_len_div else \
                          (dlen > args.max_len_div_frac * (flen_q + rlen_q) / 2.0)
                if too_far:
                    n_pairs_len_dropped += 1
                    del pairs[key]
                    continue

            # ---- molecule filter: chimeric molecules -----------------------
            if args.drop_supplementary and (fwd.has_tag("SA") or rev.has_tag("SA")):
                n_pairs_sa_dropped += 1
                del pairs[key]
                continue

            # ---- molecule filter: strand alignment non-overlap -------------
            if args.max_strand_nonoverlap:
                fs, fe = fwd.reference_start, fwd.reference_end
                rs, re_ = rev.reference_start, rev.reference_end
                nonov = (max(fe, re_) - min(fs, rs)) - max(0, min(fe, re_) - max(fs, rs))
                if nonov > args.max_strand_nonoverlap:
                    n_pairs_nonov_dropped += 1
                    del pairs[key]
                    continue

            n_pairs += 1
            want_qpos = (qpos_hist is not None and n_pairs % args.qpos_sample == 0)
            bases, info_bases, lines, pair_len, inf_fq, inf_rq = process_pair(
                key, fwd, rev, genome, args, collect_qpos=want_qpos)

            # ---- molecules with more than --max-d-per-pair variants --------
            if args.max_d_per_pair > 0:
                n_d = sum(1 for line in lines if line.endswith("\td"))
                if n_d > args.max_d_per_pair:
                    del pairs[key]
                    continue

            # ---- informative bases by distance from the read ends ----------
            if want_qpos and inf_fq:
                flen_q, rlen_q = pair_len
                fq = np.asarray(inf_fq, dtype=np.int64)
                rq = np.asarray(inf_rq, dtype=np.int64)
                for arr, d in ((qpos_hist["fwd_5"], fq),
                               (qpos_hist["fwd_3"], flen_q - 1 - fq),
                               (qpos_hist["rev_5"], rq),
                               (qpos_hist["rev_3"], rlen_q - 1 - rq)):
                    b = np.bincount(np.clip(d, 0, None) // args.qpos_bin,
                                    minlength=len(arr))
                    arr[:len(b)] += b[:len(arr)]

            np.add.at(coverage,        np.fromiter(bases,      dtype=np.int64), 1)
            np.add.at(informative_cov, np.fromiter(info_bases, dtype=np.int64), 1)

            for line in lines:
                sys.stdout.write(line + "\n")
            del pairs[key]

    with pysam.BGZFile(f"{args.prefix}.cntx.bed.gz", "w") as fh:
        for i, b in enumerate(coverage):
            if b > 0 and i > 0:
                trip = genome[chrom][i - 1 : i + 2]
                ib = informative_cov[i]
                # columns: chrom  start  end  total_cov  informative_cov  triplet
                fh.write(f"{chrom}\t{i}\t{i+1}\t{b}\t{ib}\t{trip}\n".encode())

    if qpos_hist is not None:
        # Counts from 1 molecule in --qpos-sample (recorded in the header), not rescaled
        with open(args.qpos_hist, "w") as fh:
            fh.write(f"# qpos_bin\t{args.qpos_bin}\n")
            fh.write(f"# qpos_sample\t{args.qpos_sample}\n")
            fh.write(f"# pairs_used\t{n_pairs // max(args.qpos_sample, 1)}\n")
            fh.write("bin_start\tfwd_5\tfwd_3\trev_5\trev_3\n")
            last = max((int(np.max(np.nonzero(v)[0])) if v.any() else 0)
                       for v in qpos_hist.values())
            for b in range(last + 1):
                fh.write(f"{b * args.qpos_bin}\t"
                         + "\t".join(str(int(qpos_hist[k][b]))
                                     for k in ("fwd_5", "fwd_3", "rev_5", "rev_3"))
                         + "\n")

    sys.stderr.write(
        f"Seen {coverage.sum()} total bases, "
        f"{informative_cov.sum()} informative bases\n"
    )
    if args.max_len_div or args.max_len_div_frac:
        thr = (f"{args.max_len_div} bp" if args.max_len_div
               else f"{args.max_len_div_frac:.4f} x readlen")
        tot = n_pairs + n_pairs_len_dropped + n_pairs_sa_dropped + n_pairs_nonov_dropped
        sys.stderr.write(
            f"Length-divergence filter ({thr}): dropped {n_pairs_len_dropped} "
            f"of {tot} pairs ({safediv(100 * n_pairs_len_dropped, tot):.2f}%)\n")
    if args.drop_supplementary:
        tot = n_pairs + n_pairs_sa_dropped + n_pairs_nonov_dropped
        sys.stderr.write(
            f"Supplementary-alignment filter: dropped {n_pairs_sa_dropped} "
            f"of {tot} pairs ({safediv(100 * n_pairs_sa_dropped, tot):.2f}%)\n")
    if args.max_strand_nonoverlap:
        tot = n_pairs + n_pairs_nonov_dropped
        sys.stderr.write(
            f"Strand non-overlap filter ({args.max_strand_nonoverlap} bp): dropped "
            f"{n_pairs_nonov_dropped} of {tot} pairs "
            f"({safediv(100 * n_pairs_nonov_dropped, tot):.2f}%)\n")
    bam.close()

if __name__ == "__main__":
    main()
