# =============================================================================
# call_tables.smk — variant and mismatch tables per library and read set
#
# Inputs   <results>/variant_calling/ext<ext>/<sample>/ (variant_calling.smk)
# Outputs  <results>/call_tables/ext<ext>/<sample>/
#            mutations_with_vaf.tsv       variants (d), one row per molecule carrying the non-reference allele, with
#                                         the coverage and VAF of the site (variant_table.py)
#            single_strand_calls.tsv.gz   mismatches between the two strand consensus
#                                         reads (f, r) (mismatch_table.py)
#            discordant_calls.tsv.gz      sites where the two reads carry different
#                                         non-reference bases (d_c)
# =============================================================================

rule call_tables:
    input:
        expand(TABLE_DIR + "/ext{ext}/{sample}/mutations_with_vaf.tsv",     sample=SAMPLES, ext=EXTS),
        expand(TABLE_DIR + "/ext{ext}/{sample}/single_strand_calls.tsv.gz", sample=SAMPLES, ext=EXTS),
        expand(TABLE_DIR + "/ext{ext}/{sample}/discordant_calls.tsv.gz",    sample=SAMPLES, ext=EXTS),


rule variant_table:
    input:
        muts = expand(CALL_DIR + "/ext{{ext}}/{{sample}}/{chrom}.ext{{ext}}.qual{qual}.muts.bed.gz",
                      chrom=CHROMS, qual=QUAL),
        cntx = expand(CALL_DIR + "/ext{{ext}}/{{sample}}/{chrom}.ext{{ext}}.qual{qual}.cntx.bed.gz",
                      chrom=CHROMS, qual=QUAL),
    output:
        TABLE_DIR + "/ext{ext}/{sample}/mutations_with_vaf.tsv",
    params:
        mut_dir = CALL_DIR + "/ext{ext}",
        chroms  = " ".join(CHROMS),
    resources:
        time          = "'02:59:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "8000M",
    log:
        LOGS + "/call_tables/ext{ext}/{sample}/variant_table.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/variant_table.py \
            --mut-dir {params.mut_dir} \
            --samples {wildcards.sample} \
            --chroms  {params.chroms} \
            --ext     {wildcards.ext} \
            --qual    {QUAL} \
            --vaf-out {output} \
        2>{log}
        """


rule mismatch_table:
    input:
        muts = expand(CALL_DIR + "/ext{{ext}}/{{sample}}/{chrom}.ext{{ext}}.qual{qual}.muts.bed.gz",
                      chrom=CHROMS, qual=QUAL),
    output:
        tsv = TABLE_DIR + "/ext{ext}/{sample}/single_strand_calls.tsv.gz",
        dc  = TABLE_DIR + "/ext{ext}/{sample}/discordant_calls.tsv.gz",
    params:
        mut_dir = CALL_DIR + "/ext{ext}",
        chroms  = " ".join(CHROMS),
    resources:
        time          = "'01:00:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "4000M",
    log:
        LOGS + "/call_tables/ext{ext}/{sample}/mismatch_table.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/mismatch_table.py \
            --mut-dir {params.mut_dir} \
            --samples {wildcards.sample} \
            --chroms  {params.chroms} \
            --ext     {wildcards.ext} \
            --qual    {QUAL} \
            --out     {output.tsv} \
            --discordant-out {output.dc} \
        2>{log}
        """
