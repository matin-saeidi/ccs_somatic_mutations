# notebooks

Notebooks that make the figures and tables of the paper from the outputs of
[../workflow/](../workflow/) and [../analyses/](../analyses/).

| Notebook | Makes | Saved to |
|---|---|---|
| [main_figures.ipynb](main_figures.ipynb) | Figs. 2 and 3 | `../figures/` |
| [supplementary.ipynb](supplementary.ipynb) | Supplementary Tables 2–7, Supplementary Figs. 1–15 | `../tables/`, `../figures/` |

[../figures/](../figures/) holds the figures of the paper. Fig. 1 is a schematic and has no code.

## Running

```
conda env create -f environment.yml
conda activate ccs_figures
jupyter lab
```

In the first cell of each notebook, set `RESULTS` to `results_dir` of
[../config/config.yaml](../config/config.yaml) and run it. Then run the cells of any figure or table
on their own; in supplementary.ipynb, run the Setup cells first.

## Data of other studies

[data/](data/) holds the published spectra and signatures that Fig. 3b compares ours with:

| File | Contents | Source |
|---|---|---|
| `gronska_peski_2026_al_cerebellar_spectra.tsv` | 96-type counts of the autosomal point mutations of each of the 12 *ad libitum* mice, cerebellar neurons (`<mouse>__uncorrected`), and the same counts divided by the trinucleotide ratio of that mouse (`<mouse>`) | Grońska-Pęski et al. 2026, Supplementary Tables 5 and 9 |
| `riva_2020_mSBS_mm10.tsv` | the 11 mouse signatures (mSBS), 96 types | Riva et al. 2020, `figure1/mSBSs.rds` of github.com/team113sanger/mouse-mutatation-signatures |

The COSMIC v3.4 SBS signatures (GRCh37) and the trinucleotide counts of the mm10 and
GRCh37 genomes are read from the SigProfilerAssignment and SigProfilerMatrixGenerator
packages (installed by `environment.yml`).
