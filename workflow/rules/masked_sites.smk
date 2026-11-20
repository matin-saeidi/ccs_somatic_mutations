# =============================================================================
# masked_sites.smk — masked sites [Methods: Masked sites]
#
# RepeatMasker classes config.yaml repeat_classes (Simple_repeat,
# Low_complexity, Satellite*), plus the sites with a variant called in
# config.yaml min_animals_shared (2) or more animals (config.yaml samples:
# animal) in the untrimmed reads.
# Used for both read sets.
#
# Inputs   config.yaml reference: repeatmasker; the reference and T2T .fai;
#          <results>/call_tables/ext0/<sample>/mutations_with_vaf.tsv
# Outputs  <results>/masked_sites/ (BED, 0-based half-open)
#            excluded_regions.bed   the masked sites
#            repeat_regions.bed     the RepeatMasker component
#            shared_sites.bed       the shared component
#            mask_summary.tsv       sizes of each component
# =============================================================================

rule masked_sites:
    input:
        MASK,
        MASK_DIR + "/repeat_regions.bed",
        MASK_DIR + "/shared_sites.bed",
        MASK_DIR + "/mask_summary.tsv",


rule build_mask:
    input:
        rm    = RM_OUT,
        gfai  = GCA_FAI,
        hfai  = REF_FAI,
        rates = expand(TABLE_DIR + "/ext0/{sample}/mutations_with_vaf.tsv", sample=SAMPLES),
    output:
        mask    = MASK,
        repeats = MASK_DIR + "/repeat_regions.bed",
        shared  = MASK_DIR + "/shared_sites.bed",
        summary = MASK_DIR + "/mask_summary.tsv",
    params:
        classes    = " ".join(f"'{c}'" for c in config["repeat_classes"]),
        min_shared = config["min_animals_shared"],
    resources:
        time          = "'02:00:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "16000M",
    log:
        LOGS + "/masked_sites/build_mask.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/masked_sites.py \
            --repeatmasker {input.rm} \
            --gca-fai      {input.gfai} \
            --hybrid-fai   {input.hfai} \
            --rates-dir    {TABLE_DIR}/ext0 \
            --samples      {SAMPLES} \
            --animals      {ANIMALS_ARG} \
            --chroms       {CHROMS} \
            --classes      {params.classes} \
            --min-shared   {params.min_shared} \
            --out-dir      {MASK_DIR} \
        > {log} 2>&1
        """
