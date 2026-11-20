# workflow

Snakemake pipeline from the PacBio HiFi reads (CCS by strand) to the somatic mutations,
their callable bps and the error counts. [Snakefile](Snakefile) runs the steps of
[rules/](rules/); each step calls the scripts of [scripts/](scripts/). Inputs and
parameters are in [../config/config.yaml](../config/config.yaml).

## Running

From the repository root:

```
snakemake --profile profiles/slurm -n            # dry run
snakemake --profile profiles/slurm               # every step
snakemake --profile profiles/slurm callable_bps  # one step, with the steps it needs
```

Each step is a rule named after its file (e.g. `variant_calling` in `variant_calling.smk`).

## Steps

In the order they run. Outputs are in `<results_dir>/` (config.yaml); the header of each
`.smk` file lists its outputs and their columns.

| Step | Contents | Outputs | Paper |
|---|---|---|---|
| [reference.smk](rules/reference.smk) | T2T autosomes + decoys, pbmm2 index | `config.yaml reference: out_dir` | Methods, Alignment |
| [alignment.smk](rules/alignment.smk) | pbmm2 alignment of the strand consensus reads | `alignment/` | Methods, Alignment |
| [read_tables.smk](rules/read_tables.smk) | passes and length of every strand consensus read | `read_tables/` | Supp. Tables 2-3; Supp. Figs. 2, 13 |
| [filter_thresholds.smk](rules/filter_thresholds.smk) | library thresholds: read length, strand length difference, strand alignment non-overlap | `filter_thresholds/` | Methods, Read filters, Molecule filters; Supp. Table 3 |
| [variant_calling.smk](rules/variant_calling.smk) | read, molecule and base filters; variants, mismatches and coverage per library x autosome, untrimmed (ext0) and trimmed (ext450) reads | `variant_calling/` | Methods, Variant calling and filtering |
| [call_tables.smk](rules/call_tables.smk) | variant and mismatch tables per library | `call_tables/` | |
| [masked_sites.smk](rules/masked_sites.smk) | RepeatMasker repeats + sites with a variant in two or more animals | `masked_sites/` | Methods, Masked sites |
| [callable_bps.smk](rules/callable_bps.smk) | coverage filter; callable bps per library and minimum coverage | `callable_bps/` | Methods, Coverage filter; Supp. Tables 4, 6 |
| [mutation_tables.smk](rules/mutation_tables.smk) | every variant with the filters as columns; somatic mutation candidates | `mutation_tables/`, `somatic_mutation_candidates/` | Methods, Somatic and germline mutations; Supp. Table 7 |
| [read_end_errors.smk](rules/read_end_errors.smk) | callable bps and variants by distance from the read end | `read_end_errors/` | Methods, Errors at the edges of consensus reads; Supp. Figs. 3-4 |
| [sequencing_errors.smk](rules/sequencing_errors.smk) | mismatches between strand consensus reads; errors per subread | `strand_mismatches/`, `subread_errors/` | Methods, Estimation of the sequencing error rate; Supp. Table 5; Supp. Figs. 5, 6, 14 |
| [chromatin_accessibility.smk](rules/chromatin_accessibility.smk) | accessible chromatin of cortical neurons; its callable bps | `chromatin_accessibility/` | Methods, Chromatin accessibility; Fig. 3c |

## Main outputs

| File | Contents |
|---|---|
| `somatic_mutation_candidates/somatic_mutation_candidates.tsv` | singletons at callable sites; the somatic mutations are the rows with `SINGLETONS_ON_MOLECULE <= 2` |
| `somatic_mutation_candidates/callable_bases.tsv` | callable bps per library at the minimum coverage (10): the denominator of the somatic mutation rate |
| `somatic_mutation_candidates/context_opportunity.tsv` | callable bps per library and trinucleotide |
| `somatic_mutation_candidates/t2t_reference_context_counts.tsv` | trinucleotides of the T2T autosomes |
| `callable_bps/trimmed/denominator.tsv` | callable bps per library at minimum coverages 5-12 |
| `mutation_tables/untrimmed_reads/`, `mutation_tables/450b_trimmed_reads/` | every variant of every library, filters as 0/1 columns |
| `call_tables/ext450/<library>/single_strand_calls.tsv.gz` | mismatches between the two strand consensus reads of each molecule |
| `strand_mismatches/ext450/<library>.counts.npz` | opportunities and mismatches by passes, sm and sx |
| `read_end_errors/<library>.edge_counts.npz` | callable bps and variants by distance from the read end |
| `chromatin_accessibility/callable_bps_by_region_ext450.tsv.gz` | callable bps in accessible and inaccessible chromatin |

## Resources

Each rule sets the time, CPUs and memory we requested on our cluster; adjust them to yours.
The aligned BAMs take ~250 GB per short library and ~400 GB for the long library.
