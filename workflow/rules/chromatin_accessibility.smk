# =============================================================================
# chromatin_accessibility.smk — accessible chromatin of cortical neurons and
# its callable bps [Methods: Chromatin accessibility; Fig. 3c]
#
# Accessible: snATAC-seq peaks of Sziraki et al. 2023 in which >= 2% of the
# pooled cortical neurons (cortical projection neurons 1-3, interneurons 1-2,
# wild type) have a read, carried from GRCm39 to the T2T autosomes through the
# liftOver chain; inaccessible: every other site with a GRCm39 counterpart.
#
# Inputs   config.yaml sziraki_atac_dir, liftover_chain, reference: grcm39_fasta;
#          <results>/variant_calling/, <results>/masked_sites/, <results>/callable_bps/
# Outputs  <results>/chromatin_accessibility/
#            sziraki/sziraki_atac_pseudobulk.npz        cells with a read, per peak x
#                                                       (cell type x age) (sziraki_pseudobulk.py)
#            regions/accessibility_regions.bed.gz       accessible / inaccessible sites
#            regions/accessible_peaks_mm39.bed.gz       accessible peaks on GRCm39
#            regions/accessibility_regions_summary.tsv  sizes (accessibility_regions.py)
#            callable_bps_by_region_ext450.tsv.gz       trimmed callable bps per library x
#                                                       minimum coverage x region x 32
#                                                       trinucleotides (callable_bps_by_region.py)
# =============================================================================

rule chromatin_accessibility:
    input:
        ACC_DIR + "/regions/accessibility_regions.bed.gz",
        ACC_DIR + "/callable_bps_by_region_ext450.tsv.gz",


rule sziraki_pseudobulk:
    input:
        cells = config["sziraki_atac_dir"] + "/GSM6538357_ATAC_cell_annotation.csv.gz",
        peaks = config["sziraki_atac_dir"] + "/GSM6538357_ATAC_peak_annotation.csv.gz",
        counts = config["sziraki_atac_dir"] + "/GSM6538357_ATAC_peak_count.txt.gz",
    output:
        ACC_DIR + "/sziraki/sziraki_atac_pseudobulk.npz",
    resources:
        time          = "'02:00:00'",
        ntasks        = 1,
        cpus_per_task = 2,
        mem_per_cpu   = "8000M",
    log:
        LOGS + "/chromatin_accessibility/sziraki_pseudobulk.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/sziraki_pseudobulk.py --indir {config[sziraki_atac_dir]} \
            --out {ACC_DIR}/sziraki --only atac > {log} 2>&1
        """


rule accessibility_regions:
    input:
        npz    = ACC_DIR + "/sziraki/sziraki_atac_pseudobulk.npz",
        chain  = config["liftover_chain"],
        ref    = REF,
        fai    = REF_FAI,
        grcm39 = config["reference"]["grcm39_fasta"],
    output:
        bed   = ACC_DIR + "/regions/accessibility_regions.bed.gz",
        peaks = ACC_DIR + "/regions/accessible_peaks_mm39.bed.gz",
        summ  = ACC_DIR + "/regions/accessibility_regions_summary.tsv",
    resources:
        time          = "'01:00:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "16000M",
    log:
        LOGS + "/chromatin_accessibility/accessibility_regions.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/accessibility_regions.py \
            --atac-npz {input.npz} --chain {input.chain} \
            --fai {input.fai} --chroms {CHROMS} \
            --validate-fasta {input.ref} {input.grcm39} \
            --out-dir {ACC_DIR}/regions 2> {log}
        """


# Stops if the totals over regions and trinucleotides differ from
# callable_bps/trimmed/denominator.tsv at any minimum coverage
rule callable_bps_by_region_library:
    input:
        mask    = MASK,
        denom   = CALLABLE + "/trimmed/denominator.tsv",
        regions = ACC_DIR + "/regions/accessibility_regions.bed.gz",
        cntx0   = expand(CALL_DIR + "/ext0/{{sample}}/{chrom}.ext0.qual{qual}.cntx.bed.gz",
                         chrom=CHROMS, qual=QUAL),
        cntx    = expand(CALL_DIR + "/ext450/{{sample}}/{chrom}.ext450.qual{qual}.cntx.bed.gz",
                         chrom=CHROMS, qual=QUAL),
        dist0   = CALL_DIR + "/ext0/{sample}/informative_coverage_dist.tsv",
    output:
        temp(ACC_DIR + "/callable_bps_by_region/{sample}.tsv"),
    params:
        floors = " ".join(map(str, config["min_coverages"])),
    resources:
        time          = "'02:59:00'",
        ntasks        = 1,
        cpus_per_task = 2,
        mem_per_cpu   = "6000M",
    log:
        LOGS + "/chromatin_accessibility/callable_bps_by_region_{sample}.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/callable_bps_by_region.py \
            --cntx-dir           {CALL_DIR}/ext450 \
            --untrimmed-cntx-dir {CALL_DIR}/ext0 \
            --mask               {input.mask} \
            --fai                {REF_FAI} \
            --samples            {wildcards.sample} \
            --chroms             {CHROMS} \
            --scheme             accessibility={input.regions} \
            --floors             {params.floors} \
            --ext 450 --qual {QUAL} \
            --check-against      {input.denom} \
            --out                {output} \
        2> {log}
        """


rule callable_bps_by_region:
    input:
        expand(ACC_DIR + "/callable_bps_by_region/{sample}.tsv", sample=SAMPLES),
    output:
        ACC_DIR + "/callable_bps_by_region_ext450.tsv.gz",
    shell:
        r"""
        {{ head -1 {input[0]}; tail -q -n +2 {input}; }} | gzip > {output}
        """
