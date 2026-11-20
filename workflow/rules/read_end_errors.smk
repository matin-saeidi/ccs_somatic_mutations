# =============================================================================
# read_end_errors.smk — callable bps and variants by distance from the read
# end, in the untrimmed reads [Methods: Errors at the edges of consensus reads;
# Supplementary Figs. 3-4]
#
# The reads, molecules and bases of variant_calling.smk (read end trim 0),
# counted by read_end_counts.py.
#
# Inputs   the aligned BAMs; <results>/filter_thresholds/;
#          <results>/variant_calling/ext0/ (cntx files);
#          <results>/masked_sites/repeat_regions.bed and shared_sites.bed
# Outputs  <results>/read_end_errors/
#            counts/<sample>/<chrom>.npz   one per chromosome
#            <sample>.edge_counts.npz      summed over the chromosomes of the library
# =============================================================================

rule read_end_errors:
    input:
        expand(EDGE_DIR + "/{sample}.edge_counts.npz", sample=SAMPLES),


# Counting windows of read_end_counts.py
WIN_LEN_EDGES = "3000 4000 5000"   # read-length classes of `win`: [3000,4000), [4000,5000), [5000,inf)
WINDOW        = 1500   # d < WINDOW from each end
MIN_LEN_LONG  = 5000   # reads counted in `win_long`
WINDOW_LONG   = 2500
FOLD_MAX      = 3000   # distance from the nearer end, min(d5, d3); last bin >= FOLD_MAX
COV_MAX       = 40     # coverage classes 0..COV_MAX, last bin >= COV_MAX
LEN_BIN       = 250    # read-length bins of `folded_len`, last bin >= LEN_MAX
LEN_MAX       = 16000
FD_BIN        = 100    # distance bins of `folded_len`

localrules: read_end_counts_merge


rule read_end_counts:
    input:
        bam    = BAM_DIR + "/{sample}.sorted.bam",
        bai    = BAM_DIR + "/{sample}.sorted.bam.bai",
        ref    = REF,
        cntx   = CALL_DIR + "/ext0/{sample}/{chrom}.ext0.qual" + str(QUAL) + ".cntx.bed.gz",
        bounds = THR_DIR + "/{sample}.read_length_bounds.txt",
        lendiv = THR_DIR + "/{sample}.strand_length_difference.txt",
        nonov  = THR_DIR + "/{sample}.strand_alignment_nonoverlap.txt",
        rep    = MASK_DIR + "/repeat_regions.bed",
        shared = MASK_DIR + "/shared_sites.bed",
    output:
        EDGE_DIR + "/counts/{sample}/{chrom}.npz",
    params:
        mapq      = config["min_mapping_quality"],
        indel     = config["indel_distance"],
        div_max   = config["max_divergence"],
        np_min    = config["min_passes"],
        sm_min    = config["min_subread_matches"],
        sm_ratio  = config["min_subread_match_fraction"],
        drop_supp = "--drop-supplementary" if config["drop_chimeric_molecules"] else "",
    resources:
        time          = "'11:59:00'",
        ntasks        = 1,
        cpus_per_task = 4,
        mem_per_cpu   = "4000M",
    log:
        LOGS + "/read_end_errors/{sample}/{chrom}.log",
    conda:
        ENV
    shell:
        """
        read LMIN LMAX < {input.bounds}
        read MAXLENDIV < {input.lendiv}
        read MAXNONOV < {input.nonov}
        echo "[filters] lmin=$LMIN lmax=$LMAX max_len_div=$MAXLENDIV max_strand_nonoverlap=$MAXNONOV {params.drop_supp}" > {log}

        (samtools view -q {params.mapq} -h {input.bam} {wildcards.chrom} \
        | python {SCRIPTS}/read_and_base_filters.py \
            -i {params.indel} -d {params.div_max} -e 0 \
            -np {params.np_min} -sm {params.sm_min} -r {params.sm_ratio} \
            -lmin $LMIN -lmax $LMAX \
        | python {SCRIPTS}/read_end_counts.py count \
            -i - --chrom {wildcards.chrom} --ref {input.ref} \
            --cntx {input.cntx} --repeat-bed {input.rep} --shared-bed {input.shared} \
            -q {QUAL} --max-len-div $MAXLENDIV \
            {params.drop_supp} --max-strand-nonoverlap $MAXNONOV \
            --win-len-edges {WIN_LEN_EDGES} --window {WINDOW} \
            --min-len-long {MIN_LEN_LONG} --window-long {WINDOW_LONG} \
            --fold-max {FOLD_MAX} --cov-max {COV_MAX} \
            --len-bin {LEN_BIN} --len-max {LEN_MAX} --fd-bin {FD_BIN} \
            -o {output}) 2>> {log}
        """


rule read_end_counts_merge:
    input:
        expand(EDGE_DIR + "/counts/{{sample}}/{chrom}.npz", chrom=CHROMS),
    output:
        EDGE_DIR + "/{sample}.edge_counts.npz",
    conda:
        ENV
    shell:
        "python {SCRIPTS}/read_end_counts.py merge -o {output} {input}"
