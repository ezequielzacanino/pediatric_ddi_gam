# FAERS Parsing

Downloads and transforms a pediatric curated version of FAERS from openFDA.

## Execution order

```powershell
python scripts\python\01_download_openfda_drug_event.py
python scripts\python\02_build_openfda_er_tables.py
python scripts\python\03_parse_pediatric_faers.py
R      scripts\R\04_build_ade_raw.R
```

The Python notebooks and scripts are adapted from Nicholas Giangreco,
*ngiangre/openFDA_drug_event_parsing: First release* (v1.0.0), Zenodo, 2021,
DOI [10.5281/zenodo.4464544](https://doi.org/10.5281/zenodo.4464544).

## OMOP Vocabulary 

Scripts `02` and `03` require the following files in
`data/vocabulary/vocabulary_SNOMED_MEDDRA_RxNorm_ATC/`:

| File | Used by |
|---|---|
| `CONCEPT.csv` | `02`, `03` |
| `CONCEPT_RELATIONSHIP.csv` | `02`, `03` |
| `CONCEPT_ANCESTOR.csv` | `03` |

## Handoff 

`data/processed/ade_raw.csv` is used by the following sub-projects:

```powershell
# gam_benchmark
# gam_validation 
```

## Notes

- `.openFDA.params` contains credentials.
- Raw data in `data/raw/`.
- Main output in `data/processed/ade_raw.csv`.
