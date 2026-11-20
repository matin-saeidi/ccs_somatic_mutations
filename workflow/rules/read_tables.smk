# =============================================================================
# read_tables.smk — passes and length of every strand consensus read
#
# Inputs   the aligned BAMs (alignment.smk)
# Outputs  <results>/read_tables/
#            <sample>.np_len.txt   one row per alignment record: zmw, strand
#                                  (fwd/rev), np (passes), read_length (qe - qs);
#                                  read by filter_thresholds.smk
#            <sample>.reads.npz    one entry per primary record (read_table.py);
#                                  Supplementary Tables 2-3, Supplementary Figs. 2 and 13
# =============================================================================

rule read_tables:
    input:
        expand(READ_DIR + "/{sample}.np_len.txt", sample=SAMPLES),
        expand(READ_DIR + "/{sample}.reads.npz",  sample=SAMPLES),


# Read names are <movie>/<zmw>/ccs/<strand>; tags start at SAM column 12
rule read_lengths_and_passes:
    input:
        bam = BAM_DIR + "/{sample}.sorted.bam",
        bai = BAM_DIR + "/{sample}.sorted.bam.bai",
    output:
        READ_DIR + "/{sample}.np_len.txt",
    threads: 8
    resources:
        time          = "'03:00:00'",
        ntasks        = 1,
        cpus_per_task = 8,
        mem_per_cpu   = "2000M",
    log:
        LOGS + "/read_tables/{sample}.np_len.log",
    conda:
        ENV
    shell:
        """
        samtools view -@ {threads} {input.bam} \
            | awk 'BEGIN{{OFS="\\t"; print "zmw","strand","np","read_length"}} {{
                split($1, qname, "/")
                zmw = qname[2]
                strand = qname[4]
                np=""; qs=""; qe=""
                for(i=12;i<=NF;i++) {{
                    if($i~/^np:i:/) np=substr($i,6)
                    if($i~/^qs:i:/) qs=substr($i,6)
                    if($i~/^qe:i:/) qe=substr($i,6)
                }}
                if(np!="" && qs!="" && qe!="") print zmw, strand, np, qe-qs
            }}' \
            > {output} 2>{log}
        """


rule read_table:
    input:
        bam      = BAM_DIR + "/{sample}.sorted.bam",
        bai      = BAM_DIR + "/{sample}.sorted.bam.bai",
        flagstat = BAM_DIR + "/{sample}.flagstat.txt",
    output:
        READ_DIR + "/{sample}.reads.npz",
    resources:
        time          = "'02:00:00'",
        ntasks        = 1,
        cpus_per_task = 4,
        mem_per_cpu   = "4000M",
    log:
        LOGS + "/read_tables/{sample}.reads.log",
    conda:
        ENV
    shell:
        """
        python {SCRIPTS}/read_table.py \
            --bam {input.bam} --flagstat {input.flagstat} --threads 3 \
            --out {output} 2>{log}
        """
