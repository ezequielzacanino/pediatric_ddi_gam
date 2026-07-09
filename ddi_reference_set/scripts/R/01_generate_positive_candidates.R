################################################################################
# Generation of positive triplets candidates (from the CRESCENDDI dataset)
# Script 01_generate_positive_candidates
################################################################################
#
# Builds externally-anchored list of positive-control *candidates*
# The candidate frame is the intersection of:
#   - CRESCENDDI positive controls 
# (Kontsioti et al., Sci Data 2022;9:72, doi:10.1038/s41597-022-01159-y)
#     restricted to the strongest Micromedex evidence tiers
#   - pediatric FAERS co-reporting -> the pair is actually co-administered
#

source("00_functions.R", local = TRUE)
library(openxlsx)

################################################################################
# Configuration
################################################################################

# FAERS case-level table is produced upstream and shared with gam_benchmark
# (not versioned); referenced by relative path like the OMOP vocabulary.
ade_raw_file <- "../gam_benchmark/data/processed/ade_raw.csv"

output_dir <- "./results/positive_control_candidates/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
candidates_file <- file.path(output_dir, "positive_control_candidates.csv")

# Target size of the curated positive set
target_positive_n <- 50L

# Minimum number of NICHD stages with >= 1 pediatric coadmin for a candidate to be flagged eligible
min_stages_covered <- 2L

# Minimum distinct pediatric co-reports of the pair (any stage)
min_pair_coreport <- 1L

# CRESCENDDI positive controls: one row per positive control
# interacting drugs are RxNorm Ingredient *names* in DRUG_1/2_CONCEPT_NAME. 
# if absent, download
# "Data Record 1 - Positive Controls.xlsx" from the CRESCENDDI repo
# (https://github.com/elpidakon/CRESCENDDI/tree/main/data_records). Both .xlsx and
# .csv are read.
crescenddi_file <- "./dual_curation_validation/input/Data Record 1 - Positive Controls.xlsx"
crescenddi_sheet <- 1
crescenddi_drug1_col <- "DRUG_1_CONCEPT_NAME"
crescenddi_drug2_col <- "DRUG_2_CONCEPT_NAME"
crescenddi_evid_col <- "MICROMEDEX_EVID_LEVEL"
# Strongest Micromedex tiers. "Established" alone is the highest-confidence tier
# "Probable" is added so the high-confidence pool is large enough to reach
# target_positive_n. Tighten to c("Established") if the pool is already ample.
crescenddi_evid_keep <- c("Established", "Probable")
# Optional adult interaction event, carried only as a curation hint. Left blank if
# the column is absent (its name varies across Data Record exports).
crescenddi_event_col <- "EVENT_CONCEPT_NAME"

concept_path <- ruta_concept  # shared OMOP CONCEPT.csv (from 00_functions.R)

if (!file.exists(ade_raw_file)) {
  stop(sprintf("%s was not found (FAERS case-level table from gam_benchmark).", ade_raw_file))
}
if (!file.exists(crescenddi_file)) {
  stop(sprintf(paste0(
    "CRESCENDDI positive controls not found at %s.\n",
    "  Download 'Data Record 1 - Positive Controls.xlsx' from the CRESCENDDI repo\n",
    "  (https://github.com/elpidakon/CRESCENDDI/tree/main/data_records), place it at\n",
    "  that path (the same file dual_curation_validation uses), and set\n",
    "  crescenddi_drug1_col / crescenddi_drug2_col / crescenddi_sheet to its headers."),
    crescenddi_file))
}

# Stable unordered key so {A,B} and {B,A} collapse to one identifier.
substance_pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "||")

################################################################################
# 1. CRESCENDDI positive controls
################################################################################

# RxNorm Ingredient name <-> concept_id, and the ATC 5th dropdown names
concept_dt <- fread(
  concept_path, quote = "",
  select = c("concept_id", "concept_name", "vocabulary_id", "concept_class_id", "invalid_reason")
)

rxnorm_ingredients <- concept_dt[
  tolower(trimws(vocabulary_id)) == "rxnorm" &
    trimws(concept_class_id) == "Ingredient",
  .(concept_id = as.character(concept_id), name = tolower(trimws(concept_name)))
]
ingredient_name_by_id <- setNames(rxnorm_ingredients$name, rxnorm_ingredients$concept_id)
valid_ingredient_names <- unique(rxnorm_ingredients$name)

# Resolve a vector of CRESCENDDI drug values (ids or names) to ingredient names.
resolve_ingredient_name <- function(values) {
  values <- trimws(as.character(values))
  resolved <- unname(ingredient_name_by_id[values])     # try as concept_id
  by_name <- is.na(resolved) & tolower(values) %in% valid_ingredient_names
  resolved[by_name] <- tolower(values[by_name])          # fall back to name
  resolved
}

crescenddi <- if (grepl("\\.xlsx?$", crescenddi_file, ignore.case = TRUE)) {
  read_workbook_sheet(crescenddi_file, crescenddi_sheet)
} else {
  fread(crescenddi_file)
}
for (col in c(crescenddi_drug1_col, crescenddi_drug2_col, crescenddi_evid_col)) {
  if (!col %in% names(crescenddi)) {
    stop(sprintf("Column '%s' not in %s. Available columns: %s",
                 col, crescenddi_file, paste(names(crescenddi), collapse = ", ")))
  }
}
has_event_col <- crescenddi_event_col %in% names(crescenddi)
if (!has_event_col) {
  cat(sprintf("Note: event column '%s' not found in CRESCENDDI; the event hint stays blank.\n",
              crescenddi_event_col))
}

crescenddi <- crescenddi[trimws(get(crescenddi_evid_col)) %in% crescenddi_evid_keep]

cres <- data.table(
  sub1 = resolve_ingredient_name(crescenddi[[crescenddi_drug1_col]]),
  sub2 = resolve_ingredient_name(crescenddi[[crescenddi_drug2_col]]),
  evid = trimws(as.character(crescenddi[[crescenddi_evid_col]])),
  event_hint = if (has_event_col) trimws(as.character(crescenddi[[crescenddi_event_col]])) else NA_character_
)
cres <- cres[!is.na(sub1) & !is.na(sub2) & sub1 != sub2]
cres[, pair_key := substance_pair_key(sub1, sub2)]

# One row per CRESCENDDI substance pair: keep the strongest tier seen
# semicolon-joined list of distinct adult events as the curation hint.
evid_rank <- setNames(seq_along(crescenddi_evid_keep), crescenddi_evid_keep)
cres_pair <- cres[, .(
  crescenddi_evid_level = crescenddi_evid_keep[min(evid_rank[evid])],
  crescenddi_event_hint = paste(unique(event_hint[!is.na(event_hint) & nzchar(event_hint)]),
                                collapse = "; ")
), by = pair_key]
known_pair_keys <- cres_pair$pair_key
known_substances <- unique(c(cres$sub1, cres$sub2))

################################################################################
# 2. Map CRESCENDDI substances to ATC 5th 
################################################################################
# Only ATC-5th concepts whose substance appears in a CRESCENDDI positive pair

atc5_map <- unique(concept_dt[
  tolower(trimws(vocabulary_id)) == "atc" &
    toupper(trimws(concept_class_id)) == "ATC 5TH" &
    (is.na(invalid_reason) | trimws(invalid_reason) == ""),
  .(drug_id = as.integer(concept_id), drug_name = concept_name)
], by = "drug_id")
# Substance = ATC concept_name up to the first ";" (route), lowercased.
atc5_map[, substance := tolower(trimws(sub(";.*$", "", drug_name)))]
atc5_target <- atc5_map[substance %in% known_substances]
target_ids <- unique(atc5_target$drug_id)

if (length(target_ids) == 0L) {
  stop("No ATC 5th concept matched a CRESCENDDI substance; check the drug columns.")
}

################################################################################
# 3. Pediatric FAERS co-reporting of the target pairs, per NICHD stage
################################################################################
# Each safetyreportid carries one NICHD stage
# per-stage counts partition the total co-report count

ade <- fread(ade_raw_file, select = c("safetyreportid", "atc_concept_id", "nichd"))
ade <- ade[nichd %in% niveles_nichd & atc_concept_id %in% target_ids]
report_drugs <- unique(ade[, .(safetyreportid, atc_concept_id, nichd)])

pair_co <- merge(report_drugs, report_drugs, by = c("safetyreportid", "nichd"),
                 allow.cartesian = TRUE)
pair_co <- pair_co[atc_concept_id.x < atc_concept_id.y]
stage_co <- pair_co[, .(coadmin = uniqueN(safetyreportid)),
                    by = .(drug1_id = atc_concept_id.x, drug2_id = atc_concept_id.y, nichd)]

# Resolve drug ids to ATC 5th names/substances and keep only CRESCENDDI pairs.
stage_co <- merge(stage_co,
                  atc5_target[, .(drug1_id = drug_id, drug1 = drug_name, substance1 = substance)],
                  by = "drug1_id")
stage_co <- merge(stage_co,
                  atc5_target[, .(drug2_id = drug_id, drug2 = drug_name, substance2 = substance)],
                  by = "drug2_id")
stage_co[, pair_key := substance_pair_key(substance1, substance2)]
stage_co <- stage_co[pair_key %in% known_pair_keys]

if (nrow(stage_co) == 0L) {
  stop("No CRESCENDDI positive pair is co-reported in pediatric FAERS; nothing to curate.")
}

################################################################################
# 4. Per-route candidate pairs
################################################################################
# An ATC 5th name carries the route ("amikacin; parenteral" vs "; oral"); 
# the same CRESCENDDI substance pair can appear under several routes
# Build one row per route combination, then collapse to the most co-reported route per substance
# the curator can switch route from the dossier.

stage_wide <- dcast(stage_co,
                    pair_key + drug1_id + drug2_id + drug1 + drug2 ~ nichd,
                    value.var = "coadmin", fill = 0L)
# Guarantee one column per canonical stage even if a stage never appears.
for (stage in niveles_nichd) {
  if (!stage %in% names(stage_wide)) stage_wide[, (stage) := 0L]
}
stage_cols <- niveles_nichd
stage_wide[, pair_coreport := rowSums(.SD), .SDcols = stage_cols]
stage_wide[, stages_covered := rowSums(.SD > 0L), .SDcols = stage_cols]
stage_wide[, meets_coverage := stages_covered >= min_stages_covered]
stage_wide <- stage_wide[pair_coreport >= min_pair_coreport]

# Collapse to one representative route per substance pair
setorder(stage_wide, pair_key, -meets_coverage, -pair_coreport, -stages_covered)
candidates <- unique(stage_wide, by = "pair_key")

# Attach CRESCENDDI evidence tier and adult event hint.
candidates <- merge(candidates, cres_pair, by = "pair_key", all.x = TRUE)

################################################################################
# 5. Rank so the agent meets the target with the fewest dead ends
################################################################################
# Stage coverage is the binding acceptance constraint

candidates[, evid_rank := match(crescenddi_evid_level, crescenddi_evid_keep)]
candidates[is.na(evid_rank), evid_rank := length(crescenddi_evid_keep) + 1L]
setorder(candidates, -meets_coverage, -stages_covered, -pair_coreport, evid_rank, pair_key)
candidates[, triplet_id := sprintf("P%03d", seq_len(.N))]

candidate_note <- paste0(
  "Auto-generated positive candidate from CRESCENDDI (", candidates$crescenddi_evid_level,
  "); confirm the pediatric event and attribution to the interaction, then assign a stable T-id.")

################################################################################
# 6. Candidates
################################################################################
# Workbook `triplets` columns first 

candidate_sheet <- data.table(
  triplet_id = candidates$triplet_id,
  control_type = "positive",
  drug1 = candidates$drug1, drug2 = candidates$drug2,
  event_llt = "", event_pt = "", event_hlt = "", event_hlgt = "",
  interaction_type = "",
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
  crescenddi_evid_level = candidates$crescenddi_evid_level,
  crescenddi_event_hint = candidates$crescenddi_event_hint,
  pair_coreport = candidates$pair_coreport,
  stages_covered = candidates$stages_covered,
  meets_coverage = candidates$meets_coverage
)
# Per-stage co-administration columns (coadmin_<stage>), in canonical order.
stage_out <- candidates[, ..stage_cols]
setnames(stage_out, paste0("coadmin_", stage_cols))
candidate_sheet <- cbind(candidate_sheet, stage_out)

fwrite(candidate_sheet, candidates_file)

################################################################################
# 7. Report
################################################################################

n_meets <- sum(candidate_sheet$meets_coverage)
cat("\nPositive-control candidates:", candidates_file, "\n")
cat("CRESCENDDI positive pairs (tiers", paste(crescenddi_evid_keep, collapse = "/"), "):",
    length(known_pair_keys), "\n")
cat("Co-reported in pediatric FAERS:", nrow(candidate_sheet), "\n")
cat(sprintf("Meeting coverage (>= %d NICHD stages): %d (target %d)\n",
            min_stages_covered, n_meets, target_positive_n))
if (n_meets < target_positive_n) {
  cat("  WARNING: fewer coverage-meeting candidates than the target. Widen\n",
      "  crescenddi_evid_keep or lower min_stages_covered.\n")
}
print(candidate_sheet[, .N, by = crescenddi_evid_level][order(crescenddi_evid_level)])
