# =============================================================================
# alignment.smk — strand consensus reads aligned to the reference [Methods: Alignment]
#
# pbmm2 (minimap2, preset CCS) aligns the forward and reverse strand consensus
# reads of each library independently; samtools indexes the sorted BAM.
#
# Inputs   config.yaml samples: bam (HiFi BAM per library); the reference (reference.smk)
# Outputs  <results>/alignment/
#            <sample>.sorted.bam, .sorted.bam.bai
#            <sample>.flagstat.txt    samtools flagstat
#            <sample>.coverage.txt    samtools coverage
# =============================================================================

rule alignment:
    input:
        expand(BAM_DIR + "/{sample}.sorted.bam",     sample=SAMPLES),
        expand(BAM_DIR + "/{sample}.sorted.bam.bai", sample=SAMPLES),
        expand(BAM_DIR + "/{sample}.flagstat.txt",   sample=SAMPLES),
        expand(BAM_DIR + "/{sample}.coverage.txt",   sample=SAMPLES),


# Wall time requested: 1 day for a short library (6.5-8 h measured), 2 days
# for the long library
ALIGN_TIME = {s: "'2-00:00:00'" if "10kb" in s else "'1-00:00:00'" for s in SAMPLES}


# Memory: 16 x 4 GB requested, ~48 GB peak measured. --sort-memory is per
# sort thread (6 x 6G = 36 GB in total).
rule align:
    input:
        bam = lambda wildcards: config["samples"][wildcards.sample]["bam"],
        ref = REF_MMI,
    output:
        bam = BAM_DIR + "/{sample}.sorted.bam",
        bai = BAM_DIR + "/{sample}.sorted.bam.bai",
    threads: 16
    resources:
        time          = lambda wildcards: ALIGN_TIME[wildcards.sample],
        ntasks        = 1,
        cpus_per_task = 16,
        mem_per_cpu   = "4000M",
    log:
        LOGS + "/alignment/align/{sample}.log",
    conda:
        ENV
    shell:
        """
        export TMPDIR={BAM_DIR}/tmp
        mkdir -p $TMPDIR

        pbmm2 align \
            --preset        CCS \
            --sort \
            --unmapped \
            --num-threads   {threads} \
            --sort-threads  6 \
            --sort-memory   6G \
            --log-level     INFO \
            --sample        {wildcards.sample} \
            {input.ref} \
            {input.bam} \
            {output.bam} \
        2>{log}

        samtools index {output.bam} 2>>{log}
        """


rule flagstat:
    input:
        bam = BAM_DIR + "/{sample}.sorted.bam",
        bai = BAM_DIR + "/{sample}.sorted.bam.bai",
    output:
        BAM_DIR + "/{sample}.flagstat.txt",
    threads: 4
    resources:
        time          = "'03:30:00'",
        ntasks        = 1,
        cpus_per_task = 4,
        mem_per_cpu   = "2000M",
    log:
        LOGS + "/alignment/flagstat/{sample}.log",
    conda:
        ENV
    shell:
        "samtools flagstat -@ {threads} {input.bam} > {output} 2>{log}"


rule coverage:
    input:
        bam = BAM_DIR + "/{sample}.sorted.bam",
        bai = BAM_DIR + "/{sample}.sorted.bam.bai",
    output:
        BAM_DIR + "/{sample}.coverage.txt",
    threads: 1
    resources:
        time          = "'04:00:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "8000M",
    log:
        LOGS + "/alignment/coverage/{sample}.log",
    conda:
        ENV
    shell:
        "samtools coverage {input.bam} > {output} 2>{log}"
