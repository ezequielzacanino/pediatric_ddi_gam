# GAM Validation

Methodological validation pipeline of the GAM method for detecting DDI signals in pediatric population. 
The method uses a semi-syntethic data simulation approach to inject DDI signals with developmental dynamics.
Interaction metrics in the additive and multiplicative scales are computed from the GAM and compared against their stratified counterparts (IOR and AC)


## Pipeline

From the `gam_validation/` root:

```powershell
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' 10_augmentation.R   # data augmentation and generation of control sets
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' 02_descriptive.R    # descriptive statistics and graphs
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' 20_null.R           # empirical null distribtion generation
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' 30_metrics.R        # performance evaluation
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' 40_graphs.R         # faceted figures and other graphs
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' 50_posthoc_validation.R  # post-hoc validation 
```
Each script does `source("00_functions.R")`, which in turn sources `01_theme.R`.

## Structure

- `00_functions.R`: global configuration, ATC mapping from the vocabulary and all
  pipeline functions (signal injection, GAM, IOR/AC, null, power, metrics).
- `01_theme.R`: ggplot theme and per-method color palette.
- `02_descriptive.R`: descriptives and figures for the original dataset and semi-syntethic augmented set 
- `10_augmentation.R`: builds the semi-synthetic positive/negative set and the sensitivity analysis.
- `20_null.R`: generation of the empirical null distribution by permutation. computes null thresholds for signal detection.
- `30_metrics.R`: computes performance metrics for the general, power-calibrated and intersection set evaluations.
- `40_graphs.R`: faceted metric figures and global ROC curves from the `30_metrics.R` CSVs.
- `50_posthoc_validation.R`: additional tests aimed to evaluate the methodology.

## Inputs 

- `../faers_parsing/data/processed/ade_raw.csv` - provides the curated dataset parsed by `faers_parsing/`.
- `../data/vocabulary/vocabulary_SNOMED_MEDDRA_RxNorm_ATC/CONCEPT.csv` - shared OMOP vocabulary for mapping.

## Outputs

- `results/<suffix>/augmentation_results/` - control triplets fitted by the GAM and stratified methods, co-administration counts by stage, null-pool metadata for permutation, and descriptive summaries.
- `results/<suffix>/descriptive_results/` - descriptive tables and figures for the original dataset and the control sets (produced by `02_descriptive.R`).
- `results/<suffix>/null_distribution_results/` - null distribution and `null_thresholds*.csv`.
- `results/<suffix>/metrics_results/` - performance metrics, power surfaces and figures.
- `results/<suffix>/posthoc_validation_results/` - tables and figures from the post-hoc validation.

The `suffix` encodes the GAM parameterization (default `sics`).
