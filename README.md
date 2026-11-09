# nanoPhos

Phosphoproteomics workflow for nanogram-scale inputs. Analysis code and figure
notebooks for the nanoPhos manuscript (Oliinyk et al.).

## Quick start

```bash
git clone https://github.com/MannLabs/nanoPhos-figures.git
cd nanoPhos-figures
python -m venv .venv && source .venv/bin/activate     # Python 3.11; Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Figure 4i (limma) also needs R (tested with 4.5.2) and the Bioconductor package `limma`
(`BiocManager::install("limma")`), which the code calls through `rpy2`. On Windows, start
Jupyter from cmd or PowerShell, not Git Bash, because rpy2 cannot locate R from Git Bash.

Then set up the data layout (all paths in the notebooks are relative to the repo root):

```
pride_data/     <- the analysis/ and metadata/ folders of the MassIVE deposit, unchanged
                   (pride_data/analysis/<figure2|figure3|figure4|figure5|suppl_figure2>/,
                    pride_data/metadata/fasta/); gitignored. The deposit also holds
                   analysis/figure2/selectivity_salt_optimization.tsv as a separate file
                   (Supplementary Figure 1a); put it into pride_data/analysis/figure2/
data/           <- committed reference data (funscores, mouse_kegg_annotation,
                   *_conditionSetup.tsv) — already in the repo
tables/         <- committed manuscript tables (Supplementary_Table_3, Supplementary Data 1-2, …)
figures/        <- figure outputs, created automatically on run
source_data/    <- per-panel source data (MetaInfo workbook), created automatically on run
```

Open a figure notebook and **Restart → Run All** (launch Jupyter from the repo root
so the relative data paths resolve; the notebooks also self-correct their working
directory if opened from inside `notebooks/`):

```bash
pip install jupyterlab      # or open the notebooks in any IDE with a Jupyter kernel
jupyter lab notebooks/nanoPhos_Figure2.ipynb
```

### Notebooks

All reproduction notebooks live in [`notebooks/`](notebooks/):

- **`nanoPhos_Figure{2–5}.ipynb`** and **`nanoPhos_Suppl_Figure{1–5}.ipynb`** —
  the analysis reported in the manuscript, with strict per-run Class I localization
  filtering (probability ≥ 0.75). These reproduce the main and supplementary figures.

The DIA-NN vs Spectronaut comparison (Supplementary Figure 3i–l) is part of
`nanoPhos_Suppl_Figure3.ipynb`. It reads the DIA-NN reports from
`pride_data/analysis/figure3/diann/`; the deposit README describes the search settings,
and the complete DIA-NN output and predicted library are in
`search/figure3/diann_output_sorted_hela_cells/` of the deposit.

Exploratory and superseded notebooks, including the original-submission (`v00`)
versions, are kept locally under `archive/` and are not version-controlled.

## Repository layout

| Path | Contents |
|---|---|
| `src/` | Production Python modules: `core` (Class I collapse), `analytics_core_V04`, `limma_utils`, `PeptideCollapse_v4`, `alphaPhosHelperFunctions`, `decoupling_preprocess` (Suppl. Fig. 2), `metainfo_export` (source data), `supplementary_data` (Supplementary Data 1–2) |
| `docs/` | `REPO_STRUCTURE.md`, changelogs, manuscript context |
| `data/` | Small reference CSVs cited by individual figures |
| `pride_data/` | Search-engine reports from the data deposit (gitignored, large) |
| `archive/` | Older versions, exploratory notebooks, deprecated work (gitignored) |
| `notebooks/` | Reproduction notebooks: one per main figure (`nanoPhos_Figure{2–5}.ipynb`) and per supplementary figure (`nanoPhos_Suppl_Figure{1–5}.ipynb`) |
| `tables/` | Committed manuscript tables (`Supplementary_Table_3.xlsx`, hydrophobic-peptide supplementary table, `Supplementary_Data_1.xlsx`, `Supplementary_Data_2.xlsx`) |

> **Note:** `src/alphaquant/` and `docs/phosphonetworks/` are large vendored/reference
> trees and are gitignored. AlphaQuant installs from PyPI via `requirements.txt`.

Full description and tracking policy in [`docs/REPO_STRUCTURE.md`](docs/REPO_STRUCTURE.md).

## Site-level collapse pipeline

The phosphosite collapse implementation lives in
[`src/PeptideCollapse_v4.py`](src/PeptideCollapse_v4.py). It converts Spectronaut
precursor-level BGS exports into a sites × samples quantification matrix.

The default applies **strict per-(site, run) localization filtering** at probability
≥ 0.75: for each phosphosite, the intensity in a given run is retained only if the
localization probability for that site in that specific run meets the cutoff. This
produces Class I localized phosphosite counts in the sense of Olsen et al. (2006) and
matches the per-file pre-filtering convention used by MaxQuant evidence tables and
the µPhos workflow.

A legacy `localization_strategy="global_max"` option preserves the dataset-wide-max
filtering used by the original Hogrebe et al. (2018) PeptideCollapse Perseus plugin —
available for direct A/B comparison.

## Reproducing the figures

Each figure notebook is self-contained. `Restart → Run All` regenerates panel-level
outputs given a populated `pride_data/` and the dependencies in `requirements.txt`.
The first cells of each notebook set the working directory to the repo root and
prepend `src/` to `sys.path`, so the bare module imports below them resolve to the
production code under `src/`. The last cell of each notebook writes the plotted data of
every panel to `source_data/MetaInfo_figures_checklist_v05.xlsx` (one sheet per panel).

To run all notebooks without opening them, from the repo root:

```bash
for nb in notebooks/*.ipynb; do python -m nbconvert --to notebook --execute --output-dir executed "$nb"; done
```

Some steps read precomputed results from the deposit by default. A `REGENERATE` /
`RERUN_AQ` flag at the top of the cell recomputes them: AlphaQuant results (Suppl.
Fig. 1o and 5d) and the decoupling summary tables in `data/decoupling_*.csv` (Suppl.
Fig. 2).

Supplementary Data 1 and 2 are rebuilt with `python src/supplementary_data.py`, and
Supplementary Table 2 (SIGNOR EGFR-pathway coverage) with `python src/signor_coverage.py`.

## Citation

Citation will be added on publication.

## Contact

Denys Oliinyk · <oliinyk@biochem.mpg.de> · Max Planck Institute of Biochemistry
