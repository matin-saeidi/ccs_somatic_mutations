# =============================================================================
# callable_bps.smk — coverage filter and callable bps [Methods: Coverage filter;
# Estimation of the somatic mutation rate; Supplementary Tables 4 and 6]
#
# Callable sites: unmasked sites with an untrimmed coverage from the minimum
# coverage to the library's 99th percentile of untrimmed coverage. Callable bps:
# coverage summed over the callable sites.
#
# Inputs   <results>/variant_calling/ (cntx files, coverage distributions);
#          <results>/masked_sites/excluded_regions.bed
# Outputs  <results>/callable_bps/
#            untrimmed/denominator.tsv            untrimmed reads, per library x
#                                                 minimum coverage (callable_bps_untrimmed.py)
#            trimmed/coverage_tables/<sample>.npz sites by masked or not, trinucleotide,
#                                                 trimmed and untrimmed coverage (callable_bps.py)
#            trimmed/denominator.tsv              trimmed reads at the sites chosen by
#                                                 untrimmed coverage, per library x minimum
#                                                 coverage: the denominator of the somatic
#                                                 mutation rate
#          <results>/somatic_mutation_candidates/context_opportunity.tsv
#                                                 trimmed callable bps per library x 64
#                                                 trinucleotides at config.yaml min_coverage
# =============================================================================

rule callable_bps:
    input:
        CALLABLE + "/untrimmed/denominator.tsv",
        expand(CALLABLE + "/trimmed/coverage_tables/{sample}.npz", sample=SAMPLES),
        CALLABLE + "/trimmed/denominator.tsv",
        CAND_DIR + "/context_opportunity.tsv",


# ── Untrimmed reads ──────────────────────────────────────────────────────────
rule callable_bps_untrimmed_library:
    input:
        mask = MASK,
        cntx = expand(CALL_DIR + "/ext0/{{sample}}/{chrom}.ext0.qual{qual}.cntx.bed.gz",
                      chrom=CHROMS, qual=QUAL),
        dist = CALL_DIR + "/ext0/{sample}/informative_coverage_dist.tsv",
    output:
        CALLABLE + "/untrimmed/parts/{sample}.tsv",
    params:
        floors = " ".join(map(str, config["min_coverages_untrimmed"])),
    resources:
        time          = "'02:59:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "8000M",
    log:
        LOGS + "/callable_bps/untrimmed_{sample}.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/callable_bps_untrimmed.py \
            --cntx-dir {CALL_DIR}/ext0 \
            --mask     {input.mask} \
            --samples  {wildcards.sample} \
            --chroms   {CHROMS} \
            --ext      0 \
            --qual     {QUAL} \
            --floors   {params.floors} \
            --out      {output} \
        2>{log}
        """


rule callable_bps_untrimmed:
    input:
        expand(CALLABLE + "/untrimmed/parts/{sample}.tsv", sample=SAMPLES),
    output:
        CALLABLE + "/untrimmed/denominator.tsv",
    shell:
        r"""
        {{ head -1 {input[0]}; tail -q -n +2 {input}; }} > {output}
        """


# ── Trimmed reads, sites chosen by untrimmed coverage ────────────────────────
# Stops if the trimmed and untrimmed cntx files do not pair site by site, or
# if the counts differ from either informative_coverage_dist.tsv
rule coverage_table:
    input:
        mask  = MASK,
        cntx0 = expand(CALL_DIR + "/ext0/{{sample}}/{chrom}.ext0.qual{qual}.cntx.bed.gz",
                       chrom=CHROMS, qual=QUAL),
        cntx  = expand(CALL_DIR + "/ext450/{{sample}}/{chrom}.ext450.qual{qual}.cntx.bed.gz",
                       chrom=CHROMS, qual=QUAL),
        dist0 = CALL_DIR + "/ext0/{sample}/informative_coverage_dist.tsv",
        dist  = CALL_DIR + "/ext450/{sample}/informative_coverage_dist.tsv",
    output:
        CALLABLE + "/trimmed/coverage_tables/{sample}.npz",
    resources:
        time          = "'02:59:00'",
        ntasks        = 1,
        cpus_per_task = 2,
        mem_per_cpu   = "4000M",
    log:
        LOGS + "/callable_bps/coverage_table_{sample}.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/callable_bps.py table \
            --sample   {wildcards.sample} \
            --chroms   {CHROMS} \
            --cntx-dir {CALL_DIR} \
            --mask     {input.mask} \
            --fai      {REF_FAI} \
            --out      {output} \
        > {log} 2>&1
        """


# Stops if the untrimmed callable sites differ from untrimmed/denominator.tsv
rule callable_bps_trimmed:
    input:
        tables = expand(CALLABLE + "/trimmed/coverage_tables/{sample}.npz", sample=SAMPLES),
        dists0 = expand(CALL_DIR + "/ext0/{sample}/informative_coverage_dist.tsv", sample=SAMPLES),
        check  = CALLABLE + "/untrimmed/denominator.tsv",
    output:
        denom   = CALLABLE + "/trimmed/denominator.tsv",
        context = CAND_DIR + "/context_opportunity.tsv",
    params:
        floors = " ".join(map(str, config["min_coverages"])),
        floor  = config["min_coverage"],
    resources:
        time          = "'00:30:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "4000M",
    log:
        LOGS + "/callable_bps/trimmed.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/callable_bps.py derive \
            --tables          {input.tables} \
            --cntx-dir        {CALL_DIR} \
            --fai             {REF_FAI} \
            --chroms          {CHROMS} \
            --floors          {params.floors} \
            --context-floor   {params.floor} \
            --check-untrimmed {input.check} \
            --denominator-out {output.denom} \
            --context-out     {output.context} \
        > {log} 2>&1
        """
