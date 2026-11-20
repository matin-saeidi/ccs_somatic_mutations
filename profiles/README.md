# profiles

[slurm/config.yaml](slurm/config.yaml): Snakemake profile for a SLURM cluster (the one we used).
Set `account` and `partition` to yours, then run from the repository root with
`--profile profiles/slurm`. The resources of each job are set in its rule
([../workflow/rules/](../workflow/rules/)).
