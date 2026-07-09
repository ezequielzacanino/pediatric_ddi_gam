################################################################################
# Negative control candidate generation 
# Script 03_generate_negative_candidates
################################################################################
#
# Builds matched 1:1 negative-control *candidates* 
# from the curated positive set and pediatric FAERS co-reporting. 
#
# matched 1:1,
#   - event_swap: keep a curated pair (matched on the drug pair) and attach a different curated event.
#   - drug_swap : keep a curated event (matched on the event) and replace one drug with another universe drug
# Each positive yields one suggested candidate.

source("00_functions.R", local = TRUE)
library(openxlsx)

################################################################################
# Configuration
################################################################################

input_xlsx <- "./input/ddi_reference_input.xlsx"
curated_triplets_file <-
  "./results/curated_pediatric_ddi_reference_set/curated_pediatric_ddi_triplets.csv"
# FAERS case-level table is produced upstream and shared with gam_benchmark
ade_raw_file <- "../gam_benchmark/data/processed/ade_raw.csv"

output_dir <- "./results/negative_control_candidates/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
candidates_file <- file.path(output_dir, "negative_control_candidates.csv")

# Minimum distinct pediatric co-reports of the pair for a candidate
min_pair_coreport <- 1L

# Minimum distinct pediatric co-reports of the full drug-drug-event triplet
min_triplet_coreport <- 1L

if (!file.exists(input_xlsx)) {
  stop(sprintf("%s was not found. Run scripts/R/00_build_input_template.R first.", input_xlsx))
}
if (!file.exists(curated_triplets_file)) {
  stop("Curated triplets not found. Run scripts/R/02_curate_pediatric_ddi_reference_set.R first.")
}
if (!file.exists(ade_raw_file)) {
  stop(sprintf("%s was not found (FAERS case-level table from gam_benchmark).", ade_raw_file))
}

# Stable unordered drug-pair key so {A,B} and {B,A} collapse to one identifier.
pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "_")

################################################################################
# 1. Positive set: workbook cells + ids resolved by script 02
################################################################################

# Dropdown cell strings come from the workbook
# candidates are copy-paste ready; the resolved ATC/MedDRA concept ids come from the curated output.
wb_triplets <- read_workbook_sheet(input_xlsx, "triplets")
if (!"control_type" %in% names(wb_triplets)) wb_triplets[, control_type := "positive"]
wb_triplets[is.na(control_type) | !nzchar(trimws(control_type)), control_type := "positive"]
for (event_col in c("event_llt", "event_pt", "event_hlt", "event_hlgt")) {
  if (!event_col %in% names(wb_triplets)) wb_triplets[, (event_col) := NA_character_]
  wb_triplets[, (event_col) := as.character(get(event_col))]
}
wb_pos <- wb_triplets[trimws(control_type) == "positive", .(
  triplet_id,
  drug1 = trimws(as.character(drug1)),
  drug2 = trimws(as.character(drug2)),
  event_llt, event_pt, event_hlt, event_hlgt
)]

curated <- fread(curated_triplets_file, select = c(
  "triplet_id", "drug1_atc_concept_id", "drug2_atc_concept_id",
  "meddra_concept_id", "meddra_pt"
))

positives <- merge(wb_pos, curated, by = "triplet_id")
positives[, `:=`(
  drug1_id = as.integer(drug1_atc_concept_id),
  drug2_id = as.integer(drug2_atc_concept_id),
  event_id = as.integer(meddra_concept_id)
)]
positives[, pair_id := pair_key(drug1_id, drug2_id)]
positives[, triplet_key := paste(pair_id, event_id, sep = "_")]

positive_pairs <- unique(positives$pair_id)
positive_triplets <- unique(positives$triplet_key)

# Drug and event universes 
# A pair is represented once, keeping the first curated triplet that uses it as its matched anchor.
drug_universe <- unique(rbindlist(list(
  positives[, .(drug = drug1, drug_id = drug1_id)],
  positives[, .(drug = drug2, drug_id = drug2_id)]
)))
event_universe <- unique(
  positives[, .(event_id, meddra_pt, event_llt, event_pt, event_hlt, event_hlgt)]
)
pair_rep <- unique(positives, by = "pair_id")[, .(
  pair_id, drug1, drug2, drug1_id, drug2_id, matched_triplet_id = triplet_id
)]

################################################################################
# 2. Pediatric FAERS evidence 
################################################################################

universe_ids <- unique(drug_universe$drug_id)
event_ids <- unique(event_universe$event_id)
ade <- fread(ade_raw_file,
             select = c("safetyreportid", "atc_concept_id", "meddra_concept_id", "nichd"))
ade <- ade[nichd %in% niveles_nichd & atc_concept_id %in% universe_ids]

# 2a. Pair co-reporting: distinct pediatric reports containing both drugs
report_drugs <- unique(ade[, .(safetyreportid, atc_concept_id)])
pair_co <- merge(report_drugs, report_drugs, by = "safetyreportid", allow.cartesian = TRUE)
pair_co <- pair_co[atc_concept_id.x < atc_concept_id.y]
pair_coreport <- pair_co[, .(pair_coreport = uniqueN(safetyreportid)),
                         by = .(drug_a = atc_concept_id.x, drug_b = atc_concept_id.y)]
pair_coreport[, pair_id := pair_key(drug_a, drug_b)]

# 2b. Single-drug + event support: distinct pediatric reports where one drug alone co-occurs with a candidate event.
drug_event_support <- unique(
  ade[meddra_concept_id %in% event_ids, .(safetyreportid, atc_concept_id, meddra_concept_id)]
)
drug_event_support <- drug_event_support[, .(evt_support = uniqueN(safetyreportid)),
                                         by = .(drug_id = atc_concept_id, event_id = meddra_concept_id)]

# 2c. Triplet co-reporting: distinct pediatric reports containing both drugs AND the candidate event (at MedDRA PT) 
event_reports <- fread(ade_raw_file, select = c("safetyreportid", "meddra_concept_id", "nichd"))
event_reports <- unique(event_reports[nichd %in% niveles_nichd & meddra_concept_id %in% event_ids,
                                      .(safetyreportid, event_id = meddra_concept_id)])
triplet_co <- merge(pair_co, event_reports, by = "safetyreportid", allow.cartesian = TRUE)
triplet_co[, pair_id := pair_key(atc_concept_id.x, atc_concept_id.y)]
triplet_co_counts <- triplet_co[, .(triplet_coreport = uniqueN(safetyreportid)),
                                by = .(pair_id, event_id)]

################################################################################
# 3. Generation of candidates (matched 1:1)
################################################################################

# 3a. event_swap
cross_join <- function(left, right) {
  left <- copy(left)[, join_key := 1L]
  right <- copy(right)[, join_key := 1L]
  merge(left, right, by = "join_key", allow.cartesian = TRUE)[, join_key := NULL]
}
event_swap <- cross_join(pair_rep, event_universe)
event_swap[, `:=`(strategy = "event_swap", known_interacting_pair = TRUE)]
event_swap[, triplet_key := paste(pair_id, event_id, sep = "_")]
event_swap <- event_swap[!(triplet_key %in% positive_triplets)]

# 3b. drug_swap
anchors <- rbindlist(list(
  positives[, .(matched_triplet_id = triplet_id, event_id, meddra_pt,
                event_llt, event_pt, event_hlt, event_hlgt,
                anchor = drug1, anchor_id = drug1_id)],
  positives[, .(matched_triplet_id = triplet_id, event_id, meddra_pt,
                event_llt, event_pt, event_hlt, event_hlgt,
                anchor = drug2, anchor_id = drug2_id)]
))
drug_swap <- cross_join(anchors, drug_universe)
drug_swap <- drug_swap[anchor_id != drug_id]
drug_swap[, `:=`(
  drug1 = anchor, drug1_id = anchor_id,
  drug2 = drug, drug2_id = drug_id,
  strategy = "drug_swap", known_interacting_pair = FALSE
)]
drug_swap[, pair_id := pair_key(drug1_id, drug2_id)]
drug_swap <- drug_swap[!(pair_id %in% positive_pairs)]
drug_swap[, triplet_key := paste(pair_id, event_id, sep = "_")]
drug_swap <- drug_swap[!(triplet_key %in% positive_triplets)]

keep_cols <- c("matched_triplet_id", "strategy", "known_interacting_pair",
               "pair_id", "triplet_key", "drug1", "drug2", "drug1_id", "drug2_id",
               "event_id", "meddra_pt", "event_llt", "event_pt", "event_hlt", "event_hlgt")
candidates <- rbindlist(list(event_swap[, ..keep_cols], drug_swap[, ..keep_cols]))

################################################################################
# 4. Evidence annotation and plausibility filter
################################################################################

candidates[pair_coreport, on = "pair_id", pair_coreport := i.pair_coreport]
candidates[is.na(pair_coreport), pair_coreport := 0L]
candidates[triplet_co_counts, on = .(pair_id, event_id), triplet_coreport := i.triplet_coreport]
candidates[is.na(triplet_coreport), triplet_coreport := 0L]
candidates[drug_event_support, on = .(drug1_id = drug_id, event_id), evt_support_drug1 := i.evt_support]
candidates[drug_event_support, on = .(drug2_id = drug_id, event_id), evt_support_drug2 := i.evt_support]
candidates[is.na(evt_support_drug1), evt_support_drug1 := 0L]
candidates[is.na(evt_support_drug2), evt_support_drug2 := 0L]
candidates[, single_drug_event_max := pmax(evt_support_drug1, evt_support_drug2)]

# Plausibility: keep only pairs with real pediatric co-reporting.
candidates <- candidates[pair_coreport >= min_pair_coreport]
# Detectability: keep only triplets the benchmark can actually see
candidates <- candidates[triplet_coreport >= min_triplet_coreport]

# Collapse duplicate triplets
setorder(candidates, single_drug_event_max, -pair_coreport, triplet_key)
candidates <- unique(candidates, by = "triplet_key")

# Suggested matched 1:1 pick: the top candidate per matched positive.
candidates[, suggested := seq_len(.N) == 1L, by = matched_triplet_id]

################################################################################
# 5. Candidate form 
################################################################################

setorder(candidates, -suggested, matched_triplet_id, single_drug_event_max, -pair_coreport)
candidates[, triplet_id := sprintf("N%03d", seq_len(.N))]

candidate_note <- paste0(
  "Auto-generated negative candidate (", candidates$strategy,
  "); verify documented absence in compendia/labels and check single-drug attribution."
)

# Workbook `triplets` columns in order; curator-filled fields stay blank.
candidate_sheet <- candidates[, .(
  triplet_id,
  control_type = "negative",
  drug1, drug2,
  event_llt, event_pt, event_hlt, event_hlgt,
  interaction_type = "none",
  mechanism = "",
  evidence_type = "",
  evidence_level = "",
  pediatric_population = "",
  age_range = "",
  ontogenic_modulation = "",
  higher_risk_stages = "",
  ontogeny_evidence = "",
  source_title = "",
  source_year = "",
  source_type = "",
  confidence_level = "",
  rationale = "",
  comments = candidate_note,
  # Decision-support columns (not part of the workbook schema).
  match_strategy = strategy,
  matched_positive = matched_triplet_id,
  known_interacting_pair,
  meddra_pt,
  pair_coreport,
  triplet_coreport,
  evt_support_drug1,
  evt_support_drug2,
  single_drug_event_max,
  suggested
)]

fwrite(candidate_sheet, candidates_file)

cat("\nNegative-control candidates:", candidates_file, "\n")
cat("Positives matched:", uniqueN(positives$triplet_id), "\n")
cat(sprintf("Viable candidates (pair_coreport >= %d and triplet_coreport >= %d): %d\n",
            min_pair_coreport, min_triplet_coreport, nrow(candidate_sheet)))
cat("Suggested (matched 1:1):", sum(candidate_sheet$suggested), "\n")
print(candidate_sheet[, .N, by = match_strategy][order(-N)])
