# analyses

Analyses of this study that use the outputs of [../workflow/](../workflow/) and are not part
of it. Run them after the pipeline, from the repository root:

```
snakemake -s analyses/Snakefile --profile profiles/slurm
```

| Step | Contents | Outputs | Paper |
|---|---|---|---|
| [library_concordance.smk](rules/library_concordance.smk) | the short and the long library of the 60-week-old mouse, site by site (`concordance_libraries` in [config.yaml](../config/config.yaml)) | `library_concordance/` | Methods, Concordance between the short and long libraries; Supp. Fig. 8 |
| [animal_relationships.smk](rules/animal_relationships.smk) | coverage of every library at the sites with a variant, for the tree of the animals | `animal_relationships/` | Methods, Relationships among the animals; Supp. Fig. 15; Supp. Table 7 |

Scripts: [joint_coverage.py](scripts/joint_coverage.py), [tree_coverage_scan.py](scripts/tree_coverage_scan.py),
[tree_sites.py](scripts/tree_sites.py) (`--help` for their usage).
