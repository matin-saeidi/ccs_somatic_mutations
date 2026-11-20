# =============================================================================
# sequencing_errors.smk — mismatches between the two strand consensus reads of
# each molecule, and errors of individual subreads [Methods: Estimation of the
# sequencing error rate; Supplementary Table 5, Supplementary Figs. 5, 6 and 14]
#
# Inputs   the aligned BAMs; <results>/filter_thresholds/;
#          <results>/masked_sites/excluded_regions.bed;
#          <results>/call_tables/ext450/<library>/mutations_with_vaf.tsv
# Outputs  <results>/strand_mismatches/
#            counts/ext450/<sample>/<chrom>.npz   opportunities and mismatches, one per
#                                                 chromosome (strand_mismatch_counts.py)
#            ext450/<sample>.counts.npz           summed over the chromosomes of the library
#          <results>/subread_errors/
#            <sample>.sm_sx.npz                   sm and sx of 100 bases from each of 100,000
#                                                 random reads (subread_error_sample.py)
# =============================================================================

MISM_EXT = 450          # read end trim: the trimmed reads

# histogram caps: the last bin holds everything at or above the cap
NP_MAX, SM_MAX, SX_MAX = 60, 80, 30
RUN_MAX = 8                                        # run lengths of the context counts, last bin >= RUN_MAX
CONTEXT_NP_EDGES = "5 6 7 8 9 10 12 15 20 25 35"   # lower edges of the np bins of the context counts

# libraries of the same animal (config.yaml samples: animal)
ANIMAL_LIBS = {s: [t for t in SAMPLES if ANIMAL[t] == ANIMAL[s]] for s in SAMPLES}

localrules: strand_mismatch_counts_merge


rule sequencing_errors:
    input:
        expand(MISM_DIR + f"/ext{MISM_EXT}/{{sample}}.counts.npz", sample=SAMPLES),
        expand(SUBR_DIR + "/{sample}.sm_sx.npz", sample=SAMPLES),


rule strand_mismatch_counts:
    input:
        bam    = BAM_DIR + "/{sample}.sorted.bam",
        bai    = BAM_DIR + "/{sample}.sorted.bam.bai",
        ref    = REF,
        bounds = THR_DIR + "/{sample}.read_length_bounds.txt",
        lendiv = THR_DIR + "/{sample}.strand_length_difference.txt",
        nonov  = THR_DIR + "/{sample}.strand_alignment_nonoverlap.txt",
        mask   = MASK,
        muts   = lambda w: [TABLE_DIR + f"/ext{MISM_EXT}/{s}/mutations_with_vaf.tsv" for s in ANIMAL_LIBS[w.sample]],
    output:
        MISM_DIR + f"/counts/ext{MISM_EXT}/{{sample}}/{{chrom}}.npz",
    params:
        mapq     = config["min_mapping_quality"],
        indel    = config["indel_distance"],
        div_max  = config["max_divergence"],
        np_min   = config["min_passes"],
        sm_min   = config["min_subread_matches"],
        sm_ratio = config["min_subread_match_fraction"],
    resources:
        time          = "'05:59:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "8000M",
    log:
        LOGS + f"/strand_mismatches/ext{MISM_EXT}/{{sample}}/{{chrom}}.log",
    conda:
        ENV
    shell:
        """
        read LMIN LMAX < {input.bounds}
        read MAXLENDIV < {input.lendiv}
        read MAXNONOV < {input.nonov}
        echo "[filters] lmin=$LMIN lmax=$LMAX max_len_div=$MAXLENDIV max_strand_nonoverlap=$MAXNONOV" > {log}

        python {SCRIPTS}/strand_mismatch_counts.py count \
            --bam {input.bam} --chrom {wildcards.chrom} --ref {input.ref} \
            --ext {MISM_EXT} --lmin $LMIN --lmax $LMAX --mapq {params.mapq} \
            --np-min {params.np_min} --sm-min {params.sm_min} --sm-ratio {params.sm_ratio} \
            --pad {params.indel} --div-max {params.div_max} -q {QUAL} \
            --max-len-div $MAXLENDIV --drop-supplementary --max-strand-nonoverlap $MAXNONOV \
            --mask-bed {input.mask} --mutation-tables {input.muts} \
            --np-max {NP_MAX} --sm-max {SM_MAX} --sx-max {SX_MAX} --run-max {RUN_MAX} \
            --context-np-edges {CONTEXT_NP_EDGES} \
            -o {output} 2>> {log}
        """


rule strand_mismatch_counts_merge:
    input:
        expand(MISM_DIR + f"/counts/ext{MISM_EXT}/{{{{sample}}}}/{{chrom}}.npz", chrom=CHROMS),
    output:
        MISM_DIR + f"/ext{MISM_EXT}/{{sample}}.counts.npz",
    conda:
        ENV
    shell:
        "python {SCRIPTS}/strand_mismatch_counts.py merge -o {output} {input}"


# Seed of each library: its position in config.yaml samples (0, 1, ...)
rule subread_error_sample:
    input:
        bam = BAM_DIR + "/{sample}.sorted.bam",
        bai = BAM_DIR + "/{sample}.sorted.bam.bai",
    output:
        SUBR_DIR + "/{sample}.sm_sx.npz",
    params:
        seed = lambda w: SAMPLES.index(w.sample),
    resources:
        time          = "'02:00:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "4000M",
    log:
        LOGS + "/subread_errors/{sample}.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/subread_error_sample.py \
            --bam {input.bam} --out {output} \
            --n-reads 100000 --k 100 --seed {params.seed} > {log} 2>&1
        """
