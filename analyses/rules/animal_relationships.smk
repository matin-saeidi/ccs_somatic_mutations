# =============================================================================
# animal_relationships.smk — sites and callable length of the tree of the
# animals [Methods: Relationships among the animals; Supplementary Fig. 15,
# Supplementary Table 7]
#
# Inputs   <results>/variant_calling/ext0/ (cntx files); <results>/masked_sites/;
#          <results>/callable_bps/; <results>/mutation_tables/untrimmed_reads/
# Outputs  <results>/animal_relationships/
#            scan/<chrom>.npz   untrimmed coverage of every library at the sites with a
#                               variant; sites passing the coverage filter in every
#                               library (tree_coverage_scan.py)
#            nj_sites.tsv       one row per site with a variant: coverage and variant
#                               molecules of every library (tree_sites.py)
#            callable_L.tsv     sites passing the coverage filter in every library, per
#                               chromosome and minimum coverage
# =============================================================================

rule animal_relationships:
    input:
        TREE_DIR + "/nj_sites.tsv",
        TREE_DIR + "/callable_L.tsv",


rule tree_coverage_scan:
    input:
        cntx    = expand(CALL_DIR + "/ext0/{sample}/{{chrom}}.ext0.qual40.cntx.bed.gz", sample=SAMPLES),
        repeats = MASK_DIR + "/repeat_regions.bed",
        mask    = MASK_DIR + "/excluded_regions.bed",
        denom   = CALLABLE + "/trimmed/denominator.tsv",
        calls   = UNTRIMMED_TABLE,
    output:
        TREE_DIR + "/scan/{chrom}.npz",
    resources:
        time          = "'03:00:00'",
        ntasks        = 1,
        cpus_per_task = 2,
        mem_per_cpu   = "2000M",
    log:
        LOGS + "/animal_relationships/scan_{chrom}.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/tree_coverage_scan.py --chrom {wildcards.chrom} --libraries {SAMPLES} \
            --cntx-dir {CALL_DIR}/ext0 --fai {REF_FAI} --repeats {input.repeats} --mask {input.mask} \
            --denominator {input.denom} --calls {input.calls} --out {output} > {log} 2>&1
        """


rule tree_sites:
    input:
        scans = expand(TREE_DIR + "/scan/{chrom}.npz", chrom=CHROMS),
        calls = UNTRIMMED_TABLE,
        denom = CALLABLE + "/untrimmed/denominator.tsv",
    output:
        sites = TREE_DIR + "/nj_sites.tsv",
        L     = TREE_DIR + "/callable_L.tsv",
    resources:
        time          = "'00:30:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "8000M",
    log:
        LOGS + "/animal_relationships/tree_sites.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/tree_sites.py --scans {input.scans} --calls {input.calls} \
            --untrimmed-denominator {input.denom} --out-dir {TREE_DIR} > {log} 2>&1
        """
