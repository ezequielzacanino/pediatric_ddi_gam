# Pediatric DDI GAM Workspace

Workspace with four subprojects for building and validating a pediatric
drug-drug-event (DDI) interaction benchmark using FAERS/openFDA data.

## Structure

```
FAERS_Project/
├── data/
│   └── vocabulary/         # Shared OMOP vocabulary (not versioned, ~3 GB)
├── faers_parsing/          # Download, parsing and construction of ade_raw.csv
├── ddi_reference_set/      # Curation of the pediatric reference set (upstream)
├── gam_benchmark/          # GAM validation of the curated set against ade_raw (downstream)
└── gam_validation/         # Experimental GAM validation pipeline
```

The OMOP vocabulary is centralized in `data/vocabulary/` and referenced by
relative path (`../data/vocabulary/...`) from the projects that need it.

> This repository is **code-only**: data, vocabularies, results and credentials
> are not versioned. Each pipeline regenerates its own outputs locally.

---

## faers_parsing/

**Goal:** Download FAERS data via the openFDA API, transform it with the OMOP
vocabulary and produce `ade_raw.csv` as a model-ready table of pediatric reports.

**Pipeline (in order):**

```powershell
python scripts\python\01_download_openfda_drug_event.py   # Download openFDA JSON
python scripts\python\02_build_openfda_er_tables.py        # Build ER tables with OMOP
python scripts\python\03_parse_pediatric_faers.py          # Filter pediatric patients
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' scripts\R\04_build_ade_raw.R  # Assemble ade_raw
```

**External dependencies:**

- `.openFDA.params` — API credentials (not versioned)
- `data/vocabulary/vocabulary_SNOMED_MEDDRA_RxNorm_ATC/` — OMOP vocabulary (not versioned)

**Main output:** `data/processed/ade_raw.csv`

The Python download/parsing scripts are adapted from Nicholas Giangreco,
*openFDA_drug_event_parsing* (v1.0.0), Zenodo, 2021, DOI
[10.5281/zenodo.4464544](https://doi.org/10.5281/zenodo.4464544).

---

## ddi_reference_set/

**Goal:** Build the curated set of pediatric reference triplets
(drug1 x drug2 x MedDRA event) with explicit ATC and MedDRA codes.
Upstream project: produces the curated set consumed by `gam_benchmark/`.

**Pipeline (from the `ddi_reference_set/` root):**

```powershell
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' scripts\R\02_curate_pediatric_ddi_reference_set.R
```

**Key files:**

| File | Role |
|---|---|
| `00_functions.R` | MedDRA mapping helpers and paths to the shared vocabulary |
| `scripts/R/01_curate_...` | Defines triplets, maps ATC and MedDRA PT/HLT/HLGT |

**Expected inputs:**

- `../data/vocabulary/vocabulary_SNOMED_MEDDRA_RxNorm_ATC/CONCEPT.csv`
- `../data/vocabulary/vocabulary_SNOMED_MEDDRA_RxNorm_ATC/CONCEPT_ANCESTOR.csv`

**Main outputs (deliverables to gam_benchmark):**

- `results/curated_pediatric_ddi_reference_set/curated_pediatric_ddi_triplets.csv`
- `results/curated_pediatric_ddi_reference_set/curated_pediatric_ddi_sources.csv`

---

## gam_benchmark/

**Goal:** Validate the curated set against `ade_raw.csv` using GAM models and
stratified estimators (IOR/AC). Downstream project: consumes the curated set
from `ddi_reference_set/`. The set comprises positive and negative controls; the
benchmark estimates specificity/PPV/NPV/AUC from both (a positive-only set leaves
those metrics undefined).

**Pipeline (from the `gam_benchmark/` root):**

```powershell
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' scripts\R\01_validate_pediatric_ddi_benchmark.R
```

**Key files:**

| File | Role |
|---|---|
| `00_functions.R` | GAM/stratified modeling functions, data loading and mapping |
| `scripts/R/01_validate_...` | Maps the curated set, fits models and produces metrics |

**Expected inputs (by relative path to other projects):**

- `data/processed/ade_raw.csv` (generated in `faers_parsing/`)
- `../data/vocabulary/vocabulary_SNOMED_MEDDRA_RxNorm_ATC/CONCEPT.csv` and `CONCEPT_ANCESTOR.csv`
- `../ddi_reference_set/results/curated_pediatric_ddi_reference_set/curated_pediatric_ddi_triplets.csv`
- `results/sics/null_distribution_results/null_thresholds*.csv` (from the `gam_validation/` pipeline)

**Main outputs:**

- `results/sics/benchmark_validation/benchmark_metrics.csv`
- `results/sics/benchmark_validation/benchmark_triplets_ready_for_modeling.csv`

---

## gam_validation/

Methodological validation of the GAM on semi-synthetic data (injected positive
set) vs stratified estimators (IOR/AC). Produces the null distribution and the
thresholds consumed by `gam_benchmark/`.

**Pipeline (from the `gam_validation/` root):**

```powershell
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' 10_augmentation.R
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' 20_null.R
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' 30_metrics.R
```

**Expected inputs:** `../faers_parsing/data/processed/ade_raw.csv` and
`../data/vocabulary/.../CONCEPT.csv` (ATC mapping).

Scripts: `00_functions.R`, `01_theme.R`, `02_descriptive.R`, `10_augmentation.R`, `20_null.R`, `30_metrics.R`, `40_graphs.R`, `50_posthoc_validation.R`.

---

## End-to-end reproduction

The four subprojects form a single pipeline. Run them in this order — each step's
output feeds the next:

1. **`faers_parsing/`** — download and assemble `data/processed/ade_raw.csv`.
2. **`gam_validation/`** — produce the null distribution and detection thresholds from `ade_raw.csv`.
3. **`ddi_reference_set/`** — curate the pediatric reference triplets.
4. **`gam_benchmark/`** — validate the curated set against `ade_raw.csv`, using the thresholds from step 2 and the triplets from step 3.

See each subproject's README for the exact scripts and per-step commands. All paths
are relative to the project root; nothing needs editing if the layout above is preserved.

---

## Software environment

- **R 4.4.2** (ucrt) — packages: `data.table` 1.16.2, `mgcv` 1.9.1, `MASS` 7.3-61, `akima` 0.6-3.6, `ggplot2` 4.0.1, `scales` 1.4.0, `openxlsx` 4.2.8.1, `pacman` 0.5.1. R scripts load packages via `pacman::p_load(...)`, which installs any missing package on first run.
- **Python 3.12.5** — see [`faers_parsing/requirements.txt`](faers_parsing/requirements.txt).

---

## Data availability and provenance

No data are stored in this repository; every input comes from public sources and is
regenerated locally.

- **FAERS (openFDA):** downloaded via the openFDA drug-event API in `faers_parsing/scripts/python/01_download_openfda_drug_event.py`. 
Requires an openFDA API key in `faers_parsing/.openFDA.params` (not versioned).
- **OMOP vocabulary:** OHDSI Athena export placed in `data/vocabulary/` (SNOMED, MedDRA, RxNorm, ATC). 
- **MedDRA:** MedDRA is a licensed terminology and is **not** redistributed here. Reproducing the MedDRA mapping requires a valid MedDRA license and the corresponding OMOP vocabulary export.


---

## Conventions

- Code and comments in English
- **Exception — `ddi_reference_set/`:** its curation documents (`*.md`: README,
  `INCLUSION_CRITERIA.md`, `CURATION_WORKFLOW.md`, `HUMAN_CURATION_GUIDE.md`) and
  the curation skills are written in **Spanish**, as the independent human curator is a Spanish speaker and follows them directly. All
  code (comments, variable and function names) stays in English as everywhere else.
- Each R project (`ddi_reference_set/`, `gam_benchmark/`, `gam_validation/`) has its own `00_functions.R`; `faers_parsing/` is independent
- OMOP vocabulary centralized in `data/vocabulary/`, referenced by `../data/vocabulary/...`
- Data, vocabularies, results and temporaries are **not versioned**
- R version: 4.4.2 at `C:\Program Files\R\R-4.4.2\`

---

## License

Code released under the MIT License — see [`LICENSE`](LICENSE).
