# Dual-curation validation (human vs agent)

Mini-subproject inside `ddi_reference_set/` 
Estimates the agreement between the **agentic** curation workflow and the **independent human** curator on the same
task

This validation is **self-contained**: it samples its own pairs and compares the
two curators' output.

## Workflow

1. A **stratified, enriched** sample of **50 drug-drug pairs** (25 + 25 by
   default) is drawn. All pairs come from the same base frame:
   > Pairs co-administered in pediatric FAERS (co-reported in
   `>=` `min_pair_coreport` distinct pediatric cases)
   
   They differ only in a **hidden** expected label:
   - **positive** pairs are those also listed as positive controls in candidates pairs (taken from **CRESCENDDI** (Kontsioti
     et al., *Sci Data* 2022;9:72, doi:10.1038/s41597-022-01159-y)).
   - **negative** pairs are co-reported but absent from CRESCENDDI (presumed
     co-administered without a documented interaction). 

2. Two **identical blind workbooks** are emitted, one for the human and one for the agent. Each curator independently searches for drug-drug-event triplets
   that meet `../INCLUSION_CRITERIA.md` for the assigned pairs, documenting each
   found triplet, or marking a pair as `no_triplet_found`.
3. A comparison script maps both sets to a common key and reports the agreement **in both directions** plus the union view

The curators stay blind to each other and to the expected label (agent can not access neither CRESCENDDI nor the human's workbook)

## Order of execution

Run from this directory (`dual_curation_validation/`):

```powershell
# 0. Download "Data Record 1 - Positive Controls.xlsx" from the
#    CRESCENDDI repo (https://github.com/elpidakon/CRESCENDDI/tree/main/data_records)
#    into input/. 
# Drug columns DRUG_1_CONCEPT_NAME / DRUG_2_CONCEPT_NAME are the script defaults
#
# 1. Draw the stratified sample (CRESCENDDI positives + presumed negatives) and
# Build the two blind workbooks. Reads ade_raw once to build the co-reporting frame (cached).

# 2. input/dual_curation_human.xlsx to the human reviewer
#    input/dual_curation_agent.xlsx to the agent.

# 3. Compare and compute the agreement metrics.
R scripts\R\02_compare_dual_curation.R
```

## How each curator fills the workbook

The workbook has **one sheet per assigned pair** (`pair_01` … `pair_NN`). Each pair
sheet shows its two drugs as read-only context in rows 1-2 (they are fixed by the pair). riplet table starting at row 3, with the same
vocabulary-backed dropdowns as the parent input template. 
The full step-by-step procedure both curators follow is in
[`CURATION_GUIDE.md`]. 
For each pair:

- **Found a qualifying triplet** -> add one row per triplet below the header:
  complete the event (at its finest MedDRA level) and the documentation columns. An inline citation can go in the `source_*`/`rationale` columns.
- **No qualifying triplet** -> set `no_triplet_found = yes` on the first row.

The sheet name encodes `pair_id`.
The `pairs` sheet is a read-only worklist of the assigned pairs.

## Matching and metrics

Triplets are keyed on **`(pair_id, MedDRA PT)`**: `pair_id` is shared across both workbooks, and the event is rolled up to PT (via the parent
`resolve_meddra_event_levels()`)
Each key is `matched`, `human_only` or `agent_only`.

- `comparison_triplets.csv`: one row per `(pair, PT)` key with `in_human`,
  `in_agent` and `status` (`matched`/`human_only`/`agent_only`)
- `comparison_summary.csv`: triplet counts (`human`, `agent`, `matched`,
  `human_only`, `agent_only`, `union`), `jaccard_triplets = matched / union`, and
  the **pair-level observed agreement** (each curator: did the pair yield >= 1
  triplet?).

The pair-level observed agreement is the primary statistic

The CRESCENDDI `expected_label` still seeds the sampling strata (script 01) but is
**not** a pediatric gold standard. 
CRESCENDDI is built from adult compendia (BNF / Micromedex) and labels drug-drug-*event* triplets, while this exercise labels a
*pair* and lets the curator find the event, so a "no triplet" on a
CRESCENDDI-positive pair is not necessarily a miss.

## Structure

- `scripts/R/01_sample_pairs_and_build_worksheets.R`: builds the pediatric
  co-reporting frame, draws the stratified (CRESCENDDI-positive + presumed-
  negative) sample and writes the two blind workbooks.
- `scripts/R/02_compare_dual_curation.R`: reads both completed workbooks (one
  sheet per pair) and writes the per-triplet comparison and a small summary
  (counts, Jaccard, pair-level observed agreement).
- `HUMAN_CURATION_GUIDE.md`: step-by-step procedure for the human curator, the
  human counterpart of the agent's `curar-control-positivo` skill.
- `input/dual_curation_human.xlsx`, `input/dual_curation_agent.xlsx`: the blind
  worksheets, one sheet per assigned pair (manual entry; not versioned).
- `input/Data Record 1 - Positive Controls.xlsx`: CRESCENDDI positive controls
  used to label the positive stratum (downloaded from the CRESCENDDI repo; not
  versioned; never shown to the curators).
- `results/coreported_pairs_frame.csv`: cached pediatric co-reporting frame
  (auditable).
- `results/sampled_pairs_key.csv`: the sampled pairs + co-report counts + the
  hidden `expected_label` (positive/negative). The matching key for script 02
  (via `pair_id`); the `expected_label` only seeds the sampling strata. Keep away
  from the curators.
- `results/comparison_triplets.csv`: per-`(pair, PT)` detail
  (`matched`/`human_only`/`agent_only`).
- `results/comparison_summary.csv`: triplet counts, `jaccard_triplets` and
  pair-level observed agreement.

## Configuration (top of script 01)

- `sampling_seed` (default 12345): fixed so the draw is reproducible.
- `n_positive` / `n_negative` (default 25 / 25): CRESCENDDI-matched pairs and
  presumed-negative pairs handed to each curator; their sum is the total.
- `min_pair_coreport` (default 20): plausibility floor for a pair to enter the
  base frame (applied to both strata).
- `crescenddi_file`, `crescenddi_sheet`, `crescenddi_drug1_col`,
  `crescenddi_drug2_col`: the CRESCENDDI positive controls file (default
  `input/Data Record 1 - Positive Controls.xlsx`), its sheet, and the two columns
  holding the interacting drugs (default `DRUG_1_CONCEPT_NAME` /
  `DRUG_2_CONCEPT_NAME`; RxNorm Ingredient names, or OMOP concept_ids; both
  resolved). `.xlsx` and `.csv` are both read.

## Notes

- Reuses the parent `../00_functions.R` (picklists, MedDRA resolver, NICHD and
  evidence vocabularies); every data path (vocabulary, `ade_raw`) is set
  explicitly because the parent path constants are relative to `ddi_reference_set/`.
- Requires the `openxlsx` R package (as the parent input-template scripts do).
