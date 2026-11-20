# scripts

Scripts of the pipeline steps. Each one prints its usage with `--help` (except
`read_length_bounds.py` and `build_reference.sh`, whose usage is in their header).

| Script | Contents | Step |
|---|---|---|
| [build_reference.sh](build_reference.sh) | reference FASTA (T2T autosomes, unplaced scaffolds, GRCm39 X/Y/MT) and pbmm2 index | reference |
| [read_table.py](read_table.py) | one entry per strand consensus read: passes, length, alignment | read_tables |
| [read_length_bounds.py](read_length_bounds.py) | read length filter: mean +/- 2 s.d. of the library | filter_thresholds |
| [strand_length_difference_threshold.py](strand_length_difference_threshold.py) | strand length difference filter: 99th percentile of the library | filter_thresholds |
| [strand_alignment_nonoverlap_threshold.py](strand_alignment_nonoverlap_threshold.py) | strand alignment non-overlap filter: 99th percentile of the library | filter_thresholds |
| [read_and_base_filters.py](read_and_base_filters.py) | read filters and base filters, as base quality 0 | variant_calling, read_end_errors |
| [variant_calling.py](variant_calling.py) | molecule filters; variants, mismatches and coverage | variant_calling |
| [variant_table.py](variant_table.py) | variants of a library, with coverage and VAF | call_tables |
| [mismatch_table.py](mismatch_table.py) | mismatches between the two strand consensus reads | call_tables |
| [masked_sites.py](masked_sites.py) | masked sites | masked_sites |
| [callable_bps_untrimmed.py](callable_bps_untrimmed.py) | callable bps of the untrimmed reads; functions imported by the other callable_bps scripts | callable_bps |
| [callable_bps.py](callable_bps.py) | callable bps of the trimmed reads, sites chosen by untrimmed coverage | callable_bps |
| [mutation_tables.py](mutation_tables.py) | every variant with the filters as columns | mutation_tables |
| [somatic_mutation_candidates.py](somatic_mutation_candidates.py) | singletons at callable sites | mutation_tables |
| [reference_context_counts.py](reference_context_counts.py) | trinucleotides of the T2T autosomes | mutation_tables |
| [read_end_counts.py](read_end_counts.py) | callable bps and variants by distance from the read end | read_end_errors |
| [strand_mismatch_counts.py](strand_mismatch_counts.py) | opportunities and mismatches by passes, sm and sx | sequencing_errors |
| [subread_error_sample.py](subread_error_sample.py) | random sample of per-base sm and sx | sequencing_errors |
| [sziraki_pseudobulk.py](sziraki_pseudobulk.py) | Sziraki et al. 2023 snATAC-seq per cell type and age | chromatin_accessibility |
| [accessibility_regions.py](accessibility_regions.py) | accessible and inaccessible chromatin on the T2T autosomes | chromatin_accessibility |
| [callable_bps_by_region.py](callable_bps_by_region.py) | callable bps per region | chromatin_accessibility |
