# =============================================================================
# mutation_tables.smk — every variant with the filters as columns; somatic
# mutation candidates [Methods: Variant calling and filtering; Somatic and
# germline mutations; Supplementary Table 7]
#
# Inputs   <results>/call_tables/, <results>/variant_calling/ (coverage
#          distributions), <results>/masked_sites/, <results>/callable_bps/
# Outputs  <results>/mutation_tables/
#            untrimmed_reads/mutations_untrimmed_reads.tsv
#            450b_trimmed_reads/mutations_450b_trimmed_reads.tsv
#                one row per library x site x non-reference base; filters as
#                0/1 columns (mutation_tables.py)
#          <results>/somatic_mutation_candidates/
#            somatic_mutation_candidates.tsv   singletons of the trimmed table at
#                                              callable sites; the somatic mutations
#                                              are the rows with SINGLETONS_ON_MOLECULE <= 2
#            callable_bases.tsv                callable bps at config.yaml min_coverage
#            t2t_reference_context_counts.tsv  trinucleotides of the T2T autosomes, all
#                                              and outside the masked sites
#            t2t_autosome_lengths.tsv          length, ACGT and masked bp per autosome
# =============================================================================

rule mutation_tables:
    input:
        MUT_TABLES + "/untrimmed_reads/mutations_untrimmed_reads.tsv",
        MUT_TABLES + "/450b_trimmed_reads/mutations_450b_trimmed_reads.tsv",
        CAND_DIR + "/somatic_mutation_candidates.tsv",
        CAND_DIR + "/callable_bases.tsv",
        CAND_DIR + "/t2t_reference_context_counts.tsv",
        CAND_DIR + "/t2t_autosome_lengths.tsv",


rule build_mutation_tables:
    input:
        rates   = expand(TABLE_DIR + "/ext{ext}/{sample}/mutations_with_vaf.tsv", ext=EXTS, sample=SAMPLES),
        dists   = expand(CALL_DIR + "/ext{ext}/{sample}/informative_coverage_dist.tsv", ext=EXTS, sample=SAMPLES),
        repeats = MASK_DIR + "/repeat_regions.bed",
        rm      = RM_OUT,
        gfai    = GCA_FAI,
        hfai    = REF_FAI,
    output:
        untrimmed = MUT_TABLES + "/untrimmed_reads/mutations_untrimmed_reads.tsv",
        trimmed   = MUT_TABLES + "/450b_trimmed_reads/mutations_450b_trimmed_reads.tsv",
    resources:
        time          = "'01:00:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "16000M",
    log:
        LOGS + "/mutation_tables/mutation_tables.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/mutation_tables.py \
            --calls-dir     {TABLE_DIR} \
            --calling-dir   {CALL_DIR} \
            --repeat-bed    {input.repeats} \
            --repeatmasker  {input.rm} \
            --gca-fai       {input.gfai} \
            --reference-fai {input.hfai} \
            --samples       {SAMPLES} \
            --animals       {ANIMALS_ARG} \
            --chroms        {CHROMS} \
            --out-dir       {MUT_TABLES} \
        > {log} 2>&1
        """


rule somatic_mutation_candidates:
    input:
        table = MUT_TABLES + "/450b_trimmed_reads/mutations_450b_trimmed_reads.tsv",
        denom = CALLABLE + "/trimmed/denominator.tsv",
    output:
        candidates = CAND_DIR + "/somatic_mutation_candidates.tsv",
        callable   = CAND_DIR + "/callable_bases.tsv",
    params:
        min_cov = config["min_coverage"],
    resources:
        time          = "'00:30:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "4000M",
    log:
        LOGS + "/mutation_tables/somatic_mutation_candidates.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/somatic_mutation_candidates.py \
            --mutation-table {input.table} \
            --denominator    {input.denom} \
            --min-cov        {params.min_cov} \
            --candidates-out {output.candidates} \
            --callable-out   {output.callable} \
        > {log} 2>&1
        """


rule reference_context_counts:
    input:
        ref  = REF,
        mask = MASK,
    output:
        counts  = CAND_DIR + "/t2t_reference_context_counts.tsv",
        lengths = CAND_DIR + "/t2t_autosome_lengths.tsv",
    resources:
        time          = "'01:00:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "8000M",
    log:
        LOGS + "/mutation_tables/reference_context_counts.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/reference_context_counts.py \
            --ref {input.ref} --chroms {CHROMS} --mask {input.mask} \
            --counts-out {output.counts} --lengths-out {output.lengths} \
        2> {log}
        """
