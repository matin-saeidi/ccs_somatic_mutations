# =============================================================================
# reference.smk — reference genome [Methods: Alignment]
#
# The C57BL/6J T2T assembly (GCA_964188535.1) with chromosomes X, Y and MT of
# GRCm39 added. Its autosomes are named T2T1..T2T19 and variants are called on
# these only; the unplaced T2T scaffolds and X, Y and MT are decoys.
#
# Inputs   config.yaml reference: t2t_fasta, t2t_assembly_report, grcm39_fasta
# Outputs  config.yaml reference: out_dir/
#            b6_t2t_mm39XYM.fa, .fa.fai   reference FASTA and its index
#            b6_t2t_mm39XYM.mmi           pbmm2 index (preset CCS)
#            b6_t2t_mm39XYM.fa.pkl        reference sequences read by variant_calling.py
#
# Resources in this file and the others are those we requested on our SLURM
# cluster, with the measured run times or memory where noted; adjust them to
# your data and cluster.
# =============================================================================

rule reference:
    input:
        REF, REF_FAI, REF_MMI, REF + ".pkl",


rule build_reference:
    input:
        t2t    = config["reference"]["t2t_fasta"],
        report = config["reference"]["t2t_assembly_report"],
        grcm39 = config["reference"]["grcm39_fasta"],
    output:
        fasta   = REF,
        fai     = REF_FAI,
        mmi     = REF_MMI,
        t2t_fai = GCA_FAI,
    resources:
        time          = "'04:00:00'",
        ntasks        = 1,
        cpus_per_task = 8,
        mem_per_cpu   = "4000M",
    log:
        LOGS + "/reference/build_reference.log",
    conda:
        ENV
    shell:
        """
        bash {SCRIPTS}/build_reference.sh {input.t2t} {input.report} {input.grcm39} {output.fasta} > {log} 2>&1
        """


# Reference sequences as one Python dictionary, written once before the
# variant-calling jobs read it
rule reference_pickle:
    input:
        ref = REF,
    output:
        pkl = REF + ".pkl",
    resources:
        time          = "'02:00:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "32000M",
    log:
        LOGS + "/reference/reference_pickle.log",
    conda:
        ENV
    shell:
        """
        python -c "import sys; sys.path.insert(0, '{SCRIPTS}'); import variant_calling; variant_calling.load_genome('{input.ref}')" 2>{log}
        """
