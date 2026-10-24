# Temporal Hotspots of Distributional Change in Stress and Passive mHealth Data

This repository contains R scripts for the illustrations, change-point detection simulations, hotspot simulations, sensitivity analyses, and real-data application reported in the manuscript *Identification of Temporal Hotspots of Distributional Change in Stress and Passive mHealth Data*.

## Files

**All 11 R scripts and this README are kept in one directory.** The `.R` extensions may be hidden by Windows File Explorer.

```text
README.md
Library.R
Hotspot_simulation_shared.R
Illustrations.R
Simulations_univariate.R
Simulations_bivariate.R
Simulations_hotspot.R
Simulations_hotspot_spurious.R
Simulations_hotspot_dependence.R
Simulations_hotspot_likert.R
Applications_local_updated.R
Applications_shift_local.R
```

| File | Purpose |
| --- | --- |
| `Library.R` | Core change-point detection and hotspot estimation functions used by the analyses. |
| `Hotspot_simulation_shared.R` | Shared functions for the three supplementary hotspot simulations; automatically sources `Library.R`. |
| `Illustrations.R` | Generates simulated-series and illustrative algorithm figures. |
| `Simulations_univariate.R` | Univariate change-point detection comparisons, including no-change simulations. |
| `Simulations_bivariate.R` | Bivariate change-point detection comparisons. |
| `Simulations_hotspot.R` | Main hotspot simulation (hit rate and hotspot length). |
| `Simulations_hotspot_spurious.R` | No-change and one-series-change negative-control simulations. |
| `Simulations_hotspot_dependence.R` | Sensitivity analysis under AR(1) temporal dependence. |
| `Simulations_hotspot_likert.R` | Sensitivity analysis of randomized Likert transformations using Jaccard similarity. |
| `Applications_local_updated.R` | Four-participant hotspot analysis, manuscript tables, and participant figures. |
| `Applications_shift_local.R` | Circular-shift negative-control analysis for the real-data application. |

## How to run

Open a terminal (or RStudio) with the working directory set to the folder containing these scripts. The scripts use relative filenames such as `source("Library.R")`; no separate code folders are required.

### 1. Illustrations

```bash
Rscript Illustrations.R
```

This produces `example_samplepath1.pdf`, `example_result1.pdf`, `example_samplepath2.pdf`, and `example_result2.pdf` in the working directory.

### 2. Main simulations

Each simulation is independent and can be run separately; the commands below do **not** have to be executed in a particular order.

```bash
Rscript Simulations_univariate.R
Rscript Simulations_bivariate.R
Rscript Simulations_hotspot.R
```

The main summary outputs are:

- `univariate_null_results.csv` and `univariate_change_results.csv`
- `bivariate_change_results.csv`
- `hotspot_simulation_summary.csv`

The scripts also save `.RData` results, and the hotspot simulation saves scenario-level raw results.

### 3. Supplementary hotspot simulations

These scripts automatically source `Hotspot_simulation_shared.R`, which in turn sources `Library.R`. You do **not** need to execute either helper script separately.

```bash
Rscript Simulations_hotspot_spurious.R
Rscript Simulations_hotspot_dependence.R
Rscript Simulations_hotspot_likert.R
```

The main summary outputs are `spurious_hotspot_summary.csv`, `dependence_hotspot_summary.csv`, and `likert_hotspot_summary.csv`, respectively. Corresponding raw `.RData` files are also generated.

### 4. Real-data application (restricted input required)

The application scripts require an existing file named `hotspot_4subjects_data.RData` in the **same working directory**, containing the `df_list` object with the four participants' time series. **This input dataset and the data-preparation script are not included in the public repository.** Access to participant-level data is subject to the applicable study permissions and privacy requirements.

With the authorized input available locally, run:

```bash
Rscript Applications_local_updated.R
Rscript Applications_shift_local.R
```

Run `Applications_local_updated.R` first. It produces `data_appl.RData`, `data_appl_table.csv`, `data_appl_table.tex`, and four PDFs under the automatically created `figure/` directory (`hotspot_<subject>.pdf`). The figures use the window sizes selected in `main_figure_G` and format four-digit axis labels without thousands separators.

Then `Applications_shift_local.R` uses the saved observed analyses in `data_appl.RData` to compute circular-shift controls and generate `data_appl_shift_summary.csv`, `data_appl_shift_table.tex`, and `data_appl_shift_table_with_coverage.tex`. It also writes a resumable raw CSV (`data_appl_shift_raw.csv`). Circular shifts preserve the values and within-series order of the passive-sensing data while changing their temporal alignment with self-reported stress; the resulting comparisons are descriptive, not clinical validation.

## Computing environment and reproducibility notes

- R and the packages called by the scripts are required. These include `doParallel`, `doRNG`, `foreach`, `ggplot2`, `patchwork`, and the change-point detection packages used by the comparison scripts (including `breakfast`, `changepoint`, `stepR`, `mosum`, `ecp`, `hdbinseg`, and `changepoints`). Additional packages are referenced within the scripts.
- The univariate and bivariate comparison scripts are configured for **500** Monte Carlo replications. The main hotspot and supplementary simulation scripts are configured for **200** replications per setting; the randomized Likert analysis uses **10** transformations per simulated dataset.
- Simulation scripts use `SLURM_CPUS_PER_TASK` when available and otherwise default to one worker. The simulations, particularly confidence-interval calculations, can be computationally intensive.
- The standalone scripts describe the corresponding analyses, but some manuscript results were computed through separate HPC job-array runs and merged outputs. Re-running a standalone script may not reproduce the exact reported Monte Carlo values to the final decimal if random streams or execution settings differ.
- Outputs are written into the working directory unless otherwise noted. Generated results, logs, and restricted participant-level data should be kept out of public commits.
