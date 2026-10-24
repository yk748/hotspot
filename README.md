# hotspot

This repository contains the R code accompanying the manuscript *Identification of Temporal Hotspots of Distributional Change in Stress and Passive mHealth Data* by Younghoon Kim, Sumanta Basu, and Samprit Banerjee. An earlier version of the manuscript is available on [arXiv](https://arxiv.org/abs/2511.05887) under the title *Identification of Emotionally Stressful Periods Through Tracking Changes in Statistical Features of mHealth Data*.

## Overview

The proposed framework identifies **temporal hotspots**: intervals in which statistical features of passive sensing variables and self-reported stress indicators show changes. We extend moving-sum (MOSUM) change-point detection to identify distributional changes within individual series and across paired series, incorporating mean- and variance-related features using Mahalanobis-distance-based statistics.

Hotspots are constructed using two rules: **thresholding**, which combines time points where the relevant statistics exceed calibrated thresholds, and **confidence intervals**, which use bootstrap-based pointwise confidence intervals for estimated change points. These intervals indicate candidate periods of distributional change; they are not independently validated clinical stress events.

The repository provides the research scripts used for the illustrations, numerical experiments, sensitivity analyses, and real-data application. An R package implementing the framework is under development.

## Repository contents

All **11 R files** are stored in the **same directory** as this README. No subfolders are needed for the source scripts. Some scripts create output directories automatically.

```text
hotspot/
├── README.md
├── Library.R
├── Hotspot_simulation_shared.R
├── Illustrations.R
├── Simulations_univariate.R
├── Simulations_bivariate.R
├── Simulations_hotspot.R
├── Simulations_hotspot_spurious.R
├── Simulations_hotspot_dependence.R
├── Simulations_hotspot_likert.R
├── Applications_local_updated.R
└── Applications_shift_local.R
```

| File | Purpose |
| --- | --- |
| `Library.R` | Core change-point detection, refinement, and hotspot estimation functions. |
| `Hotspot_simulation_shared.R` | Shared utilities for supplementary simulations; sources `Library.R`. |
| `Illustrations.R` | Synthetic examples and algorithm illustrations. |
| `Simulations_univariate.R` | Univariate change-point detection comparisons and no-change experiments. |
| `Simulations_bivariate.R` | Bivariate change-point detection comparisons. |
| `Simulations_hotspot.R` | Main hotspot simulation evaluating identification and interval characteristics. |
| `Simulations_hotspot_spurious.R` | Negative-control scenarios with no change or changes in only one series. |
| `Simulations_hotspot_dependence.R` | Sensitivity to AR(1) temporal dependence. |
| `Simulations_hotspot_likert.R` | Stability under randomized transformations of ordinal stress responses. |
| `Applications_local_updated.R` | Observed-data hotspot analysis, tables, and participant figures. |
| `Applications_shift_local.R` | Circular-shift negative-control analysis of the observed data. |

## Running the code

Set your R working directory to the directory containing these files before executing the scripts. The analyses use relative paths (for example, `source("Library.R")`). Scripts for distinct simulation studies may be run independently; they do not need to be run in the order shown below.

### Illustrations

```bash
Rscript Illustrations.R
```

Creates `example_samplepath1.pdf`, `example_result1.pdf`, `example_samplepath2.pdf`, and `example_result2.pdf`.

### Main numerical experiments

```bash
Rscript Simulations_univariate.R
Rscript Simulations_bivariate.R
Rscript Simulations_hotspot.R
```

Main summary outputs:

- `univariate_null_results.csv` and `univariate_change_results.csv`
- `bivariate_change_results.csv`
- `hotspot_simulation_summary.csv`

Raw simulation results are also saved as `.RData` files.

### Supplementary experiments

```bash
Rscript Simulations_hotspot_spurious.R
Rscript Simulations_hotspot_dependence.R
Rscript Simulations_hotspot_likert.R
```

These scripts source `Hotspot_simulation_shared.R`, which in turn sources `Library.R`. Their main summary outputs are `spurious_hotspot_summary.csv`, `dependence_hotspot_summary.csv`, and `likert_hotspot_summary.csv`, respectively.

### Real-data application (restricted input)

These scripts require `hotspot_4subjects_data.RData`, containing a `df_list` object with the four analyzed participants' time series. **Neither this dataset nor the data-preparation code is included in this repository.** Authorized users must place the input data file in the same working directory before running the application.

```bash
Rscript Applications_local_updated.R
Rscript Applications_shift_local.R
```

Run `Applications_local_updated.R` **before** `Applications_shift_local.R`. The first script generates `data_appl.RData`, `data_appl_table.csv`, `data_appl_table.tex`, and participant plots in `figure/hotspot_<subject>.pdf`. The second uses `data_appl.RData` to calculate circular-shift comparisons and generates `data_appl_shift_summary.csv`, `data_appl_shift_table.tex`, and `data_appl_shift_table_with_coverage.tex`, as well as a resumable `data_appl_shift_raw.csv`.

The stress-only and circular-shift analyses are descriptive diagnostics. They do not establish clinical validity or causal coupling between self-reported stress and passive sensing measures.

## Reproducibility notes

- The scripts require R and the R packages imported in the respective files. Major dependencies include `foreach`, `doParallel`, `doRNG`, `ggplot2`, `patchwork`, `breakfast`, `changepoint`, `stepR`, `mosum`, `ecp`, `hdbinseg`, and `changepoints`; inspect the scripts for additional requirements.
- Univariate and bivariate simulation scripts specify **500 Monte Carlo replications**. The main hotspot and supplementary simulation scripts specify **200 replications** per setting. The randomized Likert experiment applies **10 transformations** per simulated dataset.
- Simulation scripts read `SLURM_CPUS_PER_TASK` when available and otherwise default to one worker. Bootstrap-based calculations can be computationally expensive.
- Some reported results were produced using separate HPC job-array runs and merged outputs. Standalone scripts reproduce the analysis design, but their exact final Monte Carlo digits may differ if execution settings or random-number streams differ.
- Generated `.RData`, `.csv`, `.pdf`, `.tex`, logs, and participant-level input files are outputs or restricted data rather than source files for the public repository. Review all files before committing to a public repository.

## Reference

**Accepted journal article**

Kim, Y., Basu, S., & Banerjee, S. (in press). Identification of Temporal Hotspots of Distributional Change in Stress and Passive mHealth Data. *Sensors*.

**Earlier preprint**

Kim, Y., Basu, S., & Banerjee, S. (2025). *Identification of Emotionally Stressful Periods Through Tracking Changes in Statistical Features of mHealth Data*. arXiv. https://arxiv.org/abs/2511.05887
