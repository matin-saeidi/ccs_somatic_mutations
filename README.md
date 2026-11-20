# ccs_somatic_mutations

Code for *Cortical neurons accumulate point mutations at similar rates in mice and humans*
(Saeidi, Yoo, Stolyarova, Polleux, Przeworski and Spisak). The methods are described in the paper.

## Installation

```
conda env create -f environment.yml
conda activate snakemake
```

## Pipeline

From PacBio HiFi reads (CCS by strand) to somatic point mutations, callable bps and error counts.

1. In [config/config.yaml](config/config.yaml), set the paths of the input files (the file lists
   where to download each public one), the HiFi BAM and the animal of each library, and the output
   folder. The parameters in the file are those of the paper.
2. In [profiles/slurm/config.yaml](profiles/slurm/config.yaml), set your SLURM account and partition.
3. From the repository root:

```
snakemake --profile profiles/slurm -n                         # dry run
snakemake --profile profiles/slurm                            # pipeline
snakemake -s analyses/Snakefile --profile profiles/slurm      # analyses of this study, after the pipeline
```

The steps and their outputs are listed in [workflow/README.md](workflow/README.md) and
[analyses/README.md](analyses/README.md).

## Figures and tables

| Notebook | Makes |
|---|---|
| [notebooks/main_figures.ipynb](notebooks/main_figures.ipynb) | Figs. 2 and 3 |
| [notebooks/supplementary.ipynb](notebooks/supplementary.ipynb) | Supplementary Tables 2–7 and Supplementary Figs. 1–15 |

In the first cell of each notebook, set `RESULTS` to the output folder of the pipeline; each
figure or table can then be made on its own. Figures are saved to [figures/](figures/), which
holds the figures of the paper, and tables to `tables/`. See [notebooks/README.md](notebooks/README.md).

## Contact

Matin Saeidi (ss6917@columbia.edu).
