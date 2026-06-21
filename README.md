# Ped-CMU ExFPS — External Validation

Reproducible analysis code for the manuscript:

> **Predicting Extubation Failure in Pediatric Cardiac Patients: External
> Validation and Practice Adherence to the Ped-CMU ExFPS**

This repository contains the R script that reproduces the external-validation
analyses of the **Paediatric Chiang Mai University Extubation Failure Predictive
Score (Ped-CMU ExFPS)**.

## The score

The Ped-CMU ExFPS is a points-based score for the risk of extubation failure in
paediatric cardiac patients. On the logit scale:

```
score = 10 * history_of_reintubation
      +  4 * pneumonia
      +  1 * acyanosis
      +  6 * cyanosis (SpO2 > 85%)

linear predictor = -3.286 + 0.231 * score
```

A score of **≥ 5** defines the high-risk group.

## What the script does

`Supplementary_analysis_code.R` runs top to bottom and produces:

1. Data import and harmonisation of the development and external-validation
   datasets into a single pooled ("membership") dataset.
2. **Table 1** — baseline characteristics, development vs. validation.
3. **Table 2** — distribution of the score and its linear predictor.
4. **Table 3** — predictor–outcome associations in each dataset.
5. Case-mix dissimilarity (membership model).
6. Discrimination — AUROC, with single-cohort and development-vs-validation ROC
   curves.
7. Calibration — calibration curve, calibration slope and
   calibration-in-the-large (CITL).
8. **Table 4** — diagnostic indices at the pre-defined cut-off (score ≥ 5).
9. Decision curve analysis.
10. Exploratory analysis of clinician adherence to score-based management.

## Repository layout

```
.
├── Supplementary_analysis_code.R   # the analysis (run this)
├── README.md
└── .gitignore                      # blocks data/ and all output from being committed
```

Two folders are expected locally but are **not** included in this repository
(see *Data availability*):

```
data/      # input .xlsx datasets (de-identified)
output/    # tables, figures and exported data written by the script
```

## How to run

1. Install R (tested with **R 4.5.2**) and the required packages:

   ```r
   install.packages(c(
     "readxl", "dplyr", "tidyr", "tibble", "forcats", "stringr",
     "janitor", "haven", "ggplot2", "gtsummary", "gt", "broom",
     "pROC", "epiR", "binom"))
   ```

2. Place the two de-identified datasets in a `data/` folder:

   - `data/external_validation_dataset.xlsx`
   - `data/development_dataset.xlsx`

3. From the repository root, run the script. It reads from `data/` and writes
   tables and figures to `output/`:

   ```r
   source("Supplementary_analysis_code.R")
   ```

   The `data_dir` and `output_dir` paths can be edited at the top of the script.

## Data availability

The datasets contain de-identified patient information and are **not** included
in this repository. They are available from the corresponding author on
reasonable request.

## Citation

If you use this code, please cite the manuscript above. (Full citation to be
added on publication.)
