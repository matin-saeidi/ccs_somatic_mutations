# config

| File | Contents | To edit |
|---|---|---|
| [config.yaml](config.yaml) | input files, output folder and the parameters of every filter, named as in the paper's Methods | the paths |
| [../profiles/slurm/config.yaml](../profiles/slurm/config.yaml) | Snakemake profile for a SLURM cluster: submit command and default resources | `account`, `partition` |

The resources of each job (time, CPUs, memory) are set in its rule, in
[../workflow/rules/](../workflow/rules/).

On another scheduler, or on one machine, replace the profile, e.g.
`snakemake --cores 16 --use-conda`.
