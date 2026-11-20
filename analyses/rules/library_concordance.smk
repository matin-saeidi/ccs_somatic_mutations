# =============================================================================
# library_concordance.smk — the short and the long library of the 60-week-old
# mouse (config.yaml concordance_libraries), site by site [Methods: Concordance
# between the short and long libraries of the 60-week-old mouse; Supplementary Fig. 8]
#
# Inputs   <results>/variant_calling/ (cntx files of both read sets);
#          <results>/masked_sites/excluded_regions.bed; <results>/mutation_tables/
# Outputs  <results>/library_concordance/
#            call_positions.tsv             sites with a variant of either library
#            joint_untrimmed_coverage.npz   unmasked sites by the untrimmed coverage of
#                                           both libraries, with their trimmed coverage;
#                                           both coverages at every site of call_positions.tsv
#                                           (joint_coverage.py)
# =============================================================================

CONC_LIBS = " ".join(config["concordance_libraries"])

rule library_concordance:
    input:
        CONC_DIR + "/joint_untrimmed_coverage.npz",


rule concordance_call_positions:
    input:
        UNTRIMMED_TABLE, TRIMMED_TABLE,
    output:
        CONC_DIR + "/call_positions.tsv",
    resources:
        time          = "'00:30:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "8000M",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/joint_coverage.py positions --mutation-tables {input} \
            --libraries {CONC_LIBS} --out-dir {CONC_DIR}
        """


rule joint_coverage_chrom:
    input:
        positions = CONC_DIR + "/call_positions.tsv",
        mask      = MASK_DIR + "/excluded_regions.bed",
        cntx      = lambda w: [CALL_DIR + f"/ext{e}/{lib}/{w.chrom}.ext{e}.qual40.cntx.bed.gz"
                               for e in (0, 450) for lib in config["concordance_libraries"]],
    output:
        temp(CONC_DIR + "/parts/{chrom}.npz"),
    resources:
        time          = "'01:00:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "16000M",
    log:
        LOGS + "/library_concordance/{chrom}.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/joint_coverage.py chrom --chrom {wildcards.chrom} \
            --positions {input.positions} --cntx-dir {CALL_DIR} --mask {input.mask} \
            --fai {REF_FAI} --libraries {CONC_LIBS} --out-dir {CONC_DIR} > {log} 2>&1
        """


rule joint_coverage:
    input:
        expand(CONC_DIR + "/parts/{chrom}.npz", chrom=CHROMS),
    output:
        CONC_DIR + "/joint_untrimmed_coverage.npz",
    resources:
        time          = "'00:30:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "4000M",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/joint_coverage.py merge --libraries {CONC_LIBS} --out-dir {CONC_DIR}
        """
