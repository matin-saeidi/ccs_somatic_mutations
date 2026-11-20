# =============================================================================
# variant_calling.smk — variants, mismatches and coverage, per library x
# autosome x read set [Methods: Variant calling and filtering]
#
# samtools view (mapping quality filter)
#   | read_and_base_filters.py (read filters and base filters)
#   | variant_calling.py (molecule filters, then the comparison of the two
#     strand consensus reads of each molecule with each other and with the reference)
#
# Read sets (ext, config.yaml read_end_trim): ext0 = untrimmed reads,
# ext450 = trimmed reads (450 bases at each read end not used).
#
# Inputs   the aligned BAMs; the reference; <results>/filter_thresholds/
# Outputs  <results>/variant_calling/ext<ext>/<sample>/
#            <chrom>.ext<ext>.qual40.muts.bed.gz (+ .tbi)  variants (d), mismatches (f, r)
#                                                          and discordant sites (d_c)
#            <chrom>.ext<ext>.qual40.cntx.bed.gz           coverage of every covered site:
#                                                          chrom start end total_cov
#                                                          informative_cov triplet
#            <chrom>.ext<ext>.qual40.qpos.tsv              informative bases by distance
#                                                          from each read end (1 molecule in 100)
#            informative_coverage_dist.tsv                 sites by coverage (informative_cov)
#            total_coverage_dist.tsv                       sites by total_cov
# =============================================================================

rule variant_calling:
    input:
        expand(CALL_DIR + "/ext{ext}/{sample}/{chrom}.ext{ext}.qual{qual}.muts.bed.gz",
               sample=SAMPLES, chrom=CHROMS, ext=EXTS, qual=QUAL),
        expand(CALL_DIR + "/ext{ext}/{sample}/{chrom}.ext{ext}.qual{qual}.cntx.bed.gz",
               sample=SAMPLES, chrom=CHROMS, ext=EXTS, qual=QUAL),
        expand(CALL_DIR + "/ext{ext}/{sample}/{chrom}.ext{ext}.qual{qual}.qpos.tsv",
               sample=SAMPLES, chrom=CHROMS, ext=EXTS, qual=QUAL),
        expand(CALL_DIR + "/ext{ext}/{sample}/total_coverage_dist.tsv",       sample=SAMPLES, ext=EXTS),
        expand(CALL_DIR + "/ext{ext}/{sample}/informative_coverage_dist.tsv", sample=SAMPLES, ext=EXTS),


# Memory: ~9 GB peak per job; wall time: 0.5-2 h measured per job
rule call_variants:
    input:
        bam    = BAM_DIR + "/{sample}.sorted.bam",
        bai    = BAM_DIR + "/{sample}.sorted.bam.bai",
        ref    = REF,
        pkl    = REF + ".pkl",
        bounds = THR_DIR + "/{sample}.read_length_bounds.txt",
        lendiv = THR_DIR + "/{sample}.strand_length_difference.txt",
        nonov  = THR_DIR + "/{sample}.strand_alignment_nonoverlap.txt",
    output:
        muts = CALL_DIR + "/ext{ext}/{sample}/{chrom}.ext{ext}.qual{qual}.muts.bed.gz",
        cntx = CALL_DIR + "/ext{ext}/{sample}/{chrom}.ext{ext}.qual{qual}.cntx.bed.gz",
        qpos = CALL_DIR + "/ext{ext}/{sample}/{chrom}.ext{ext}.qual{qual}.qpos.tsv",
    params:
        pref      = CALL_DIR + "/ext{ext}/{sample}/{chrom}.ext{ext}.qual{qual}",
        mapq      = config["min_mapping_quality"],
        indel     = config["indel_distance"],
        div_max   = config["max_divergence"],
        np_min    = config["min_passes"],
        sm_min    = config["min_subread_matches"],
        sm_ratio  = config["min_subread_match_fraction"],
        drop_supp = "--drop-supplementary" if config["drop_chimeric_molecules"] else "",
    resources:
        time          = "'5:59:00'",
        ntasks        = 1,
        cpus_per_task = 4,
        mem_per_cpu   = "8000M",
    log:
        LOGS + "/variant_calling/ext{ext}/{sample}/{chrom}.ext{ext}.qual{qual}.log",
    conda:
        ENV
    shell:
        """
        read LMIN LMAX < {input.bounds}
        read MAXLENDIV < {input.lendiv}
        read MAXNONOV < {input.nonov}
        echo "[filters] read length $LMIN-$LMAX, strand length difference $MAXLENDIV," \
             "strand alignment non-overlap $MAXNONOV, {params.drop_supp}" > {log}

        (samtools view -q {params.mapq} -h {input.bam} {wildcards.chrom} \
        | python {SCRIPTS}/read_and_base_filters.py \
            -i {params.indel} \
            -d {params.div_max} \
            -e {wildcards.ext} \
            -np {params.np_min} \
            -sm {params.sm_min} \
            -r {params.sm_ratio} \
            -lmin $LMIN \
            -lmax $LMAX \
        | python {SCRIPTS}/variant_calling.py \
            -i - \
            -r {input.ref} \
            -b 10000 \
            -q {wildcards.qual} \
            -p {params.pref} \
            --max-d-per-pair 0 \
            --max-len-div $MAXLENDIV \
            {params.drop_supp} \
            --max-strand-nonoverlap $MAXNONOV \
            --qpos-hist {output.qpos} \
            --qpos-bin 25 \
            --qpos-sample 100 \
        | sort-bed --max-mem 16G - \
        | bgzip > {output.muts}) \
        2>>{log}

        tabix -p bed {output.muts} 2>>{log}
        """


# Sites by coverage, per chromosome: <kind>\t<coverage>\t<count>
# (kind = total | informative)
rule coverage_dist_chrom:
    input:
        cntx = CALL_DIR + "/ext{ext}/{sample}/{chrom}.ext{ext}.qual" + str(QUAL) + ".cntx.bed.gz",
    output:
        partial = temp(CALL_DIR + "/ext{ext}/{sample}/cov_partials/{chrom}.cov.partial.tsv"),
    resources:
        time          = "'00:45:00'",
        ntasks        = 1,
        cpus_per_task = 4,
        mem_per_cpu   = "2000M",
    log:
        LOGS + "/variant_calling/ext{ext}/{sample}/coverage_dist.{chrom}.log",
    conda:
        ENV
    shell:
        """
        bgzip -c -d -@ {resources.cpus_per_task} {input.cntx} \
        | awk -F'\\t' '
            NF < 6 {{
                print "ERROR: <6 cols in {input.cntx} at line " NR ": " $0 > "/dev/stderr";
                err = 1; exit 1
            }}
            {{ t[$4]++; i[$5]++ }}
            END {{
                if (err) exit 1;
                for (c in t) print "total\\t"       c "\\t" t[c];
                for (c in i) print "informative\\t" c "\\t" i[c];
            }}' \
        > {output.partial} 2>>{log}
        """


# Sites by coverage, summed over the autosomes: coverage\tn_sites, by coverage
rule coverage_dist_merge:
    input:
        partials = expand(
            CALL_DIR + "/ext{{ext}}/{{sample}}/cov_partials/{chrom}.cov.partial.tsv",
            chrom=CHROMS,
        ),
    output:
        total       = CALL_DIR + "/ext{ext}/{sample}/total_coverage_dist.tsv",
        informative = CALL_DIR + "/ext{ext}/{sample}/informative_coverage_dist.tsv",
    resources:
        time          = "'00:15:00'",
        ntasks        = 1,
        cpus_per_task = 1,
        mem_per_cpu   = "2000M",
    log:
        LOGS + "/variant_calling/ext{ext}/{sample}/coverage_dist.merge.log",
    shell:
        """
        {{
            printf 'coverage\\tn_sites\\n'
            awk -F'\\t' '$1=="total"       {{ a[$2]+=$3 }} END {{ for (c in a) print c "\\t" a[c] }}' {input.partials} \
            | sort -k1,1n
        }} > {output.total} 2>>{log}

        {{
            printf 'coverage\\tn_sites\\n'
            awk -F'\\t' '$1=="informative" {{ a[$2]+=$3 }} END {{ for (c in a) print c "\\t" a[c] }}' {input.partials} \
            | sort -k1,1n
        }} > {output.informative} 2>>{log}
        """
