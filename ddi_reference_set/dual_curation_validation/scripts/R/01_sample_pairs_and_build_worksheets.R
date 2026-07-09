################################################################################
# Dual-curation validation - sampling and blind worksheet builder
# Script 01_sample_pairs_and_build_worksheets
################################################################################
#
# Draws a sample of drug-drug pairs and emits two identical blind workbooks 
#
# Generates two identical sheets (one for human curator and other for the agent)
#  
# Script 02 estimates agreement between curators
#
# Owing to the low prevalence of DDI, the sample is not completly random
# Presumably positive: random sample from CRESCENDDI (Kontsioti et al.)
# Contains DDI from adult population
# 
# Presumably negative: random sample from real coadministrations in FAERS
# Probabilistically negative 
# 
# CRESCENDDI drugs are RxNorm Ingredient concepts 
# the workbook pairs are ATC 5th ("substance; route") 
# There is no clean id-level ATC-5th -> RxNorm-Ingredient. the two are matched on the normalized substance

source("../00_functions.R")
library(openxlsx)

################################################################################
# Configuration
################################################################################

# Refuse to clobber worksheets a curator may already have filled in.
overwrite_existing <- TRUE

sampling_seed <- 12345        # fixed so the draw is reproducible and auditable
n_positive <- 50L                 # CRESCENDDI-backed pairs (presumed positive)
n_negative <- 50L                 # co-reported pairs absent from CRESCENDDI and in FAERS

n_pairs <- n_positive + n_negative  # total pairs handed to each curator
min_pair_coreport <- 20L          # plausibility floor: pair co-reported in >= this

voc_dir <- "../../data/vocabulary/vocabulary_SNOMED_MEDDRA_RxNorm_ATC"
concept_path <- file.path(voc_dir, "CONCEPT.csv")
ade_raw_file <- "../../gam_benchmark/data/processed/ade_raw.csv"

input_dir <- "./input"
results_dir <- "./results"
dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# CRESCENDDI positive controls 
# One row per positive control
# drugs are identified by RxNorm Ingredient *name* in DRUG_1_CONCEPT_NAME /DRUG_2_CONCEPT_NAME  
# Downloaded "Data Record 1 - Positive Controls.xlsx"
# from the CRESCENDDI repo and place it at the path below:
# https://github.com/elpidakon/CRESCENDDI/tree/main/data_records

crescenddi_file <- file.path(input_dir, "Data Record 1 - Positive Controls.xlsx")
crescenddi_sheet <- 1                 # sheet holding the positive-control table
crescenddi_drug1_col <- "DRUG_1_CONCEPT_NAME"
crescenddi_drug2_col <- "DRUG_2_CONCEPT_NAME"

# Restrict positives to the strongest CRESCENDDI evidence tier
crescenddi_evid_col <- "MICROMEDEX_EVID_LEVEL"
crescenddi_evid_keep <- "Established"

frame_file <- file.path(results_dir, "coreported_pairs_frame.csv")
sample_key_file <- file.path(results_dir, "sampled_pairs_key.csv")
human_xlsx <- file.path(input_dir, "dual_curation_human.xlsx")
agent_xlsx <- file.path(input_dir, "dual_curation_agent.xlsx")

if ((file.exists(human_xlsx) || file.exists(agent_xlsx)) && !overwrite_existing) {
  stop(sprintf(
    "%s / %s already exist. Set overwrite_existing <- TRUE to rebuild (this discards manual curation).",
    human_xlsx, agent_xlsx
  ))
}
if (!file.exists(ade_raw_file)) {
  stop(sprintf("%s was not found (FAERS case-level table from gam_benchmark).", ade_raw_file))
}
if (!file.exists(crescenddi_file)) {
  stop(sprintf(paste0(
    "CRESCENDDI positive controls not found at %s.\n",
    "  Download 'Data Record 1 - Positive Controls.xlsx' from the CRESCENDDI repo\n",
    "  (https://github.com/elpidakon/CRESCENDDI/tree/main/data_records), place it at\n",
    "  that path, and set crescenddi_drug1_col / crescenddi_drug2_col / crescenddi_sheet\n",
    "  to match its headers."),
    crescenddi_file))
}

################################################################################
# 1. Build the pediatric co-reporting frame: pediatric co-reported ATC-5th pairs
################################################################################
# Self-join over the full pediatric drug universe

ade <- fread(ade_raw_file, select = c("safetyreportid", "atc_concept_id", "nichd"))
ade <- ade[nichd %in% niveles_nichd & !is.na(atc_concept_id)]

# Distinct (report, drug) membership; a self-join counts every co-exposed pair.
report_drugs <- unique(ade[, .(safetyreportid, atc_concept_id)])
pair_co <- merge(report_drugs, report_drugs, by = "safetyreportid", allow.cartesian = TRUE)
pair_co <- pair_co[atc_concept_id.x < atc_concept_id.y]
pair_frame <- pair_co[, .(pair_coreport = uniqueN(safetyreportid)),
                      by = .(drug1_id = atc_concept_id.x, drug2_id = atc_concept_id.y)]
rm(pair_co); gc()

# Keep only plausibly co-administered pairs.
pair_frame <- pair_frame[pair_coreport >= min_pair_coreport]

################################################################################
# 2. Resolve drug ids to the ATC 5th dropdown value ("substance; route")
################################################################################

# Mapping to the same concept_name the workbook dropdowns use keeps the sampled pairs copy-paste compatible with the curation schema

concept_dt <- fread(
  concept_path, quote = "",
  select = c("concept_id", "concept_name", "vocabulary_id", "concept_class_id", "invalid_reason")
)
atc5_map <- unique(concept_dt[
  tolower(trimws(vocabulary_id)) == "atc" &
    toupper(trimws(concept_class_id)) == "ATC 5TH" &
    (is.na(invalid_reason) | trimws(invalid_reason) == ""),
  .(drug_id = as.integer(concept_id), drug_name = concept_name)
], by = "drug_id")
# Substance = ATC concept_name up to the first ";" (route), lowercased.
atc5_map[, substance := tolower(trimws(sub(";.*$", "", drug_name)))]

# Inner joins drop any non-ATC-5th concept id, leaving a clean drug-level frame.
pair_frame <- merge(pair_frame,
                    atc5_map[, .(drug1_id = drug_id, drug1 = drug_name, substance1 = substance)],
                    by = "drug1_id")
pair_frame <- merge(pair_frame,
                    atc5_map[, .(drug2_id = drug_id, drug2 = drug_name, substance2 = substance)],
                    by = "drug2_id")

setorder(pair_frame, -pair_coreport, drug1_id, drug2_id)
fwrite(pair_frame[, .(drug1, drug2, drug1_id, drug2_id, pair_coreport)], frame_file)

################################################################################
# 3. Flag pairs that are CRESCENDDI positive controls
################################################################################
# CRESCENDDI drugs are resolved to RxNorm Ingredient substance names 
# A frame pair is a known positive when its two substances form one of those pairs.

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

# Data Record 1 ships as .xlsx; a .csv export is read too.
crescenddi <- if (grepl("\\.xlsx?$", crescenddi_file, ignore.case = TRUE)) {
  as.data.table(read.xlsx(crescenddi_file, sheet = crescenddi_sheet))
} else {
  fread(crescenddi_file)
}
for (col in c(crescenddi_drug1_col, crescenddi_drug2_col, crescenddi_evid_col)) {
  if (!col %in% names(crescenddi)) {
    stop(sprintf("Column '%s' not in %s. Available columns: %s",
                 col, crescenddi_file, paste(names(crescenddi), collapse = ", ")))
  }
}

# Keep only Micromedex "Established" positive controls.
crescenddi <- crescenddi[trimws(get(crescenddi_evid_col)) == crescenddi_evid_keep]

cres_drug1 <- resolve_ingredient_name(crescenddi[[crescenddi_drug1_col]])
cres_drug2 <- resolve_ingredient_name(crescenddi[[crescenddi_drug2_col]])
keep <- !is.na(cres_drug1) & !is.na(cres_drug2) & cres_drug1 != cres_drug2
known_pair_keys <- unique(paste(pmin(cres_drug1[keep], cres_drug2[keep]),
                                pmax(cres_drug1[keep], cres_drug2[keep]), sep = "||"))

pair_frame[, substance_pair_key :=
             paste(pmin(substance1, substance2), pmax(substance1, substance2), sep = "||")]
pair_frame[, is_known_ddi := substance_pair_key %in% known_pair_keys]

# An ATC-5th name carries the route ("substance; oral" vs "; parenteral")
# Collapses to one row per unordered substance pair (keeping the most co-reported route)
# curator never sees the same pair twice
pair_frame <- pair_frame[!duplicated(substance_pair_key)]

################################################################################
# 4. Draw the stratified sample and blind the order
################################################################################

positive_frame <- pair_frame[is_known_ddi == TRUE]
negative_frame <- pair_frame[is_known_ddi == FALSE]

if (nrow(positive_frame) < n_positive) {
  stop(sprintf(paste0(
    "Only %d co-reported pairs match a CRESCENDDI positive control, fewer than the\n",
    "  requested %d positives. Lower n_positive or min_pair_coreport, or check that the\n",
    "  CRESCENDDI drug columns resolve to ingredient names."),
    nrow(positive_frame), n_positive))
}
if (nrow(negative_frame) < n_negative) {
  stop(sprintf("Only %d presumed-negative pairs available, fewer than the requested %d.",
               nrow(negative_frame), n_negative))
}

set.seed(sampling_seed)
positive_sample <- positive_frame[sample(.N, n_positive)]
negative_sample <- negative_frame[sample(.N, n_negative)]
positive_sample[, expected_label := "positive"]
negative_sample[, expected_label := "negative"]

# Mix and shuffle so pair_id order leaks nothing about the expected label.
sampled <- rbind(positive_sample, negative_sample)
sampled <- sampled[sample(.N)]
sampled[, pair_id := seq_len(.N)]

# Audit key (NOT given to the curators): carries the hidden truth label
# script 02 / downstream analysis can score against it
fwrite(sampled[, .(pair_id, drug1, drug2, drug1_id, drug2_id,
                   pair_coreport, expected_label)],
       sample_key_file)

################################################################################
# 5. Build the two identical blind workbooks
################################################################################

picklists <- build_vocabulary_picklists(concept_path = concept_path)

# NICHD dropdown options mirror the parent input template.
nichd_options <- c(
  niveles_nichd,
  "term_neonatal,infancy",
  "term_neonatal,infancy,toddler",
  "infancy,toddler",
  "early_childhood,middle_childhood",
  "early_adolescence,late_adolescence"
)

ref_lists <- list(
  ref_llt = picklists$llt,
  ref_pt = picklists$pt,
  ref_hlt = picklists$hlt,
  ref_hlgt = picklists$hlgt,
  ref_nichd = nichd_options
)

# each pair gets its own sheet (pair_NN)
# no_triplet_found lets the curator explicitly record a pair that yields no qualifying triplet
# empty result is unambiguous (vs "not done yet")
pair_triplet_cols <- c(
  "no_triplet_found", "triplet_id",
  "event_llt", "event_pt", "event_hlt", "event_hlgt",
  "interaction_type", "mechanism", "evidence_type", "evidence_level",
  "pediatric_population", "age_range",
  "ontogenic_modulation", "higher_risk_stages", "ontogeny_evidence",
  "source_doi", "source_year", "source_type", "confidence_level",
  "rationale", "comments"
)

# The triplet table starts at this sheet row (rows 1-2 hold the drug context)
# script 02 reads each pair sheet with the same startRow.
header_row <- 3L

# Empty (header-only) table; the curator appends one row per triplet found.
empty_table <- do.call(
  data.table,
  setNames(rep(list(character()), length(pair_triplet_cols)), pair_triplet_cols)
)

# Column positions within the per-pair table, kept in sync with pair_triplet_cols.
col_no_triplet <- match("no_triplet_found", pair_triplet_cols)
col_llt <- match("event_llt", pair_triplet_cols)
col_pt <- match("event_pt", pair_triplet_cols)
col_hlt <- match("event_hlt", pair_triplet_cols)
col_hlgt <- match("event_hlgt", pair_triplet_cols)
col_interaction <- match("interaction_type", pair_triplet_cols)
col_evidence_level <- match("evidence_level", pair_triplet_cols)
col_ontogeny <- match("ontogenic_modulation", pair_triplet_cols)
col_higher_risk <- match("higher_risk_stages", pair_triplet_cols)
col_confidence <- match("confidence_level", pair_triplet_cols)

build_curation_workbook <- function(path, role_label) {
  wb <- createWorkbook()
  addWorksheet(wb, "pairs")

  # Hidden reference sheets backing the list validations.
  for (ref_name in names(ref_lists)) {
    addWorksheet(wb, ref_name, visible = FALSE)
    writeData(wb, ref_name, data.frame(term = ref_lists[[ref_name]]))
  }

  # pairs sheet: the read-only worklist of the assigned pairs
  pairs_sheet <- sampled[, .(pair_id, drug1, drug2)]
  writeData(wb, "pairs", pairs_sheet, withFilter = TRUE)
  freezePane(wb, "pairs", firstRow = TRUE)

  validation_rows <- (header_row + 1L):(header_row + 300L)
  range_ref <- function(ref_name) {
    sprintf("'%s'!$A$2:$A$%d", ref_name, length(ref_lists[[ref_name]]) + 1L)
  }

  # One sheet per assigned pair: the curator fills one row per triplet found
  # sets no_triplet_found = yes if the pair yields none.
  for (i in seq_len(nrow(sampled))) {
    sheet <- sprintf("pair_%02d", sampled$pair_id[i])
    addWorksheet(wb, sheet)

    # Rows 1-2: the pair this sheet is for (read-only context; not entered
    # the curator cannot pick a drug other than the pair's).
    writeData(wb, sheet, "drug1:", startCol = 1, startRow = 1, colNames = FALSE)
    writeData(wb, sheet, sampled$drug1[i], startCol = 2, startRow = 1, colNames = FALSE)
    writeData(wb, sheet, "drug2:", startCol = 1, startRow = 2, colNames = FALSE)
    writeData(wb, sheet, sampled$drug2[i], startCol = 2, startRow = 2, colNames = FALSE)

    # Header at header_row, data below it; freeze the context + header rows.
    writeData(wb, sheet, empty_table, startRow = header_row, withFilter = TRUE)
    freezePane(wb, sheet, firstActiveRow = header_row + 1L, firstActiveCol = 1)

    add_list <- function(col, value) {
      dataValidation(wb, sheet, cols = col, rows = validation_rows,
                     type = "list", value = value, allowBlank = TRUE)
    }
    add_list(col_no_triplet, '"yes"')
    add_list(col_llt, range_ref("ref_llt"))
    add_list(col_pt, range_ref("ref_pt"))
    add_list(col_hlt, range_ref("ref_hlt"))
    add_list(col_hlgt, range_ref("ref_hlgt"))
    add_list(col_higher_risk, range_ref("ref_nichd"))
    add_list(col_interaction, '"pharmacokinetic,pharmacodynamic,mixed,unknown,pharmaceutical"')
    add_list(col_evidence_level, sprintf('"%s"', paste(evidence_levels, collapse = ",")))
    add_list(col_ontogeny, '"yes,no,unknown"')
    add_list(col_confidence, '"high,moderate"')
  }

  saveWorkbook(wb, path, overwrite = TRUE)
  cat(sprintf("  %-6s workbook: %s\n", role_label, path))
}

build_curation_workbook(human_xlsx, "human")
build_curation_workbook(agent_xlsx, "agent")

################################################################################
# 6. Report
################################################################################

cat("\nPediatric co-reporting frame (ATC-5th pairs, co-report >=", min_pair_coreport, "):\n")
cat("  pairs in frame:", nrow(pair_frame), "\n")
cat("  CRESCENDDI-matched (positive) pairs available:", nrow(positive_frame), "\n")
cat("  presumed-negative pairs available:", nrow(negative_frame), "\n")
cat("  frame cached:", frame_file, "\n")
cat("Sampled pairs:", n_pairs, sprintf("(%d positive + %d negative, seed %d)\n",
                                       n_positive, n_negative, sampling_seed))
cat("  sample key (hidden truth labels):", sample_key_file, "\n")
cat("\nHand dual_curation_human.xlsx to the human and dual_curation_agent.xlsx to\n")
cat("the agent; each fills one sheet per pair (pair_01..pair_NN) following\n")
cat("HUMAN_CURATION_GUIDE.md, then run script 02 to compare.\n")
