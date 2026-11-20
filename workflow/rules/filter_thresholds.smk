# =============================================================================
# filter_thresholds.smk — library thresholds of the read and molecule filters
# [Methods: Read filters, Molecule filters; Supplementary Table 3]
#
# Inputs   <results>/read_tables/<sample>.np_len.txt; the aligned BAMs
# Outputs  <results>/filter_thresholds/
#            <sample>.read_length_bounds.txt           read length filter: "<min>\t<max>",
#                                                      the mean +/- 2 s.d. (Supplementary Fig. 2)
#            <sample>.strand_length_difference.txt     strand length difference filter:
#                                                      the library's 99th percentile, in bases
#            <sample>.strand_alignment_nonoverlap.txt  strand alignment non-overlap filter:
#                                                      the library's 99th percentile, in sites
#          Line 1 holds the threshold; the other lines are comments.
# =============================================================================

rule filter_thresholds:
    input:
        expand(THR_DIR + "/{sample}.read_length_bounds.txt",          sample=SAMPLES),
        expand(THR_DIR + "/{sample}.strand_length_difference.txt",    sample=SAMPLES),
        expand(THR_DIR + "/{sample}.strand_alignment_nonoverlap.txt", sample=SAMPLES),


rule read_length_bounds:
    input:
        READ_DIR + "/{sample}.np_len.txt",
    output:
        THR_DIR + "/{sample}.read_length_bounds.txt",
    resources:
        time          = "'00:30:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "2000M",
    conda:
        ENV
    shell:
        "python {SCRIPTS}/read_length_bounds.py {input} {output}"


rule strand_length_difference_threshold:
    input:
        READ_DIR + "/{sample}.np_len.txt",
    output:
        THR_DIR + "/{sample}.strand_length_difference.txt",
    params:
        pct = config["strand_length_difference_percentile"],
    resources:
        time          = "'00:30:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "8000M",
    log:
        LOGS + "/filter_thresholds/{sample}.strand_length_difference.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/strand_length_difference_threshold.py \
            -i {input} -o {output} --pct {params.pct} --stats 2>{log}
        """


# Molecules that pass the read filters, the strand length difference filter
# and the chimeric molecule filter
rule strand_alignment_nonoverlap_threshold:
    input:
        bam    = BAM_DIR + "/{sample}.sorted.bam",
        bai    = BAM_DIR + "/{sample}.sorted.bam.bai",
        bounds = THR_DIR + "/{sample}.read_length_bounds.txt",
        lendiv = THR_DIR + "/{sample}.strand_length_difference.txt",
    output:
        THR_DIR + "/{sample}.strand_alignment_nonoverlap.txt",
    params:
        pct     = config["strand_alignment_nonoverlap_percentile"],
        regions = ",".join(CHROMS),
        np_min  = config["min_passes"],
        div_max = config["max_divergence"],
        mapq    = config["min_mapping_quality"],
    resources:
        time          = "'02:00:00'",
        ntasks        = 1,
        cpus_per_task = 4,
        mem_per_cpu   = "8000M",
    log:
        LOGS + "/filter_thresholds/{sample}.strand_alignment_nonoverlap.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/strand_alignment_nonoverlap_threshold.py \
            -b {input.bam} --bounds {input.bounds} --lendiv {input.lendiv} \
            --regions {params.regions} --np-min {params.np_min} \
            --div-max {params.div_max} --min-mapq {params.mapq} --pct {params.pct} \
            --threads 3 -o {output} --stats 2>{log}
        """
