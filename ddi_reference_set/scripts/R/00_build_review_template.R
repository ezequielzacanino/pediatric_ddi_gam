################################################################################
# Builder of the single review workbook for auditing the agent's curation work.
# Script 00_build_review_template
################################################################################
#
# Produces reviews/agent_review.xlsx, 
# Human review workbook of all triplets the agent curated. 

# The conceptual rubric and the literature mapping in REVIEW_GUIDE.md.
#
# Data sheets are seeded from the input workbook 
# Re-run only to (re)generate the empty/seeded workbook
# Set overwrite_existing <- TRUE for a clean rebuild.
#

library(openxlsx)
library(data.table)

# Guard against destroying manual review. Flip to TRUE only for a clean rebuild.
overwrite_existing <- TRUE

input_xlsx  <- "./input/ddi_reference_input.xlsx"
reviews_dir <- "./reviews"
output_xlsx <- file.path(reviews_dir, "agent_review.xlsx")
dir.create(reviews_dir, showWarnings = FALSE, recursive = TRUE)

if (file.exists(output_xlsx) && !overwrite_existing) {
  stop(sprintf(
    "%s already exists. Set overwrite_existing <- TRUE to rebuild.",
    output_xlsx
  ))
}

################################################################################
# 1. Seed the worklist from the triplets registered in the input workbook
################################################################################
# The agent records every curated triplet as a row in the 'triplets' sheet
# For each triplet_id, locate its dossier under agent/workspace/{positivos,negativos} by id prefix
# A missing dossier leaves the path blank, which itself flags the omission for the reviewer.

locate_dossier <- function(triplet_id) {
  hits <- list.files(
    c("agent/workspace/positivos", "agent/workspace/negativos"),
    pattern = sprintf("^%s_.*\\.md$", triplet_id),
    full.names = TRUE
  )
  if (length(hits)) hits[1] else NA_character_
}

triplets <- as.data.table(read.xlsx(input_xlsx, sheet = "triplets"))
triplets <- triplets[!is.na(triplet_id) & nzchar(triplet_id)]

seed <- data.table(
  triplet_id = triplets$triplet_id,
  dossier = vapply(triplets$triplet_id, locate_dossier, character(1))
)
setorder(seed, triplet_id)

################################################################################
# 2. review sheet (one row per triplet)
################################################################################
# Column order is fixed; the dropdown column indices below must stay in sync.

review_cols <- c(
  "triplet_id", "dossier", "reviewer", "review_date", "review_minutes",
  "verdict", "verdict_reason",
  # performance rubric (see REVIEW_GUIDE.md)
  "sources_all_real",        # A: no fabricated citations (yes/no)
  "metadata_correct",        # A: authors/year/journal/title match the real source
  "all_claims_supported",    # B: each attributed claim is actually in the source
  "sources_relevant",        # C: cited sources are pertinent to the pair+event
  "best_evidence_used",      # C: used the strongest evidence available
  "pediatric_valid",         # D: direct pediatric evidence of the pair (<21)
  "mechanism_sound",         # E: declared mechanism correct, coherent with type
  "interaction_attributed",  # E: event attributed to the interaction (not 1 drug)
  "mapping_correct",         # F: ATC drugs and MedDRA event mapped correctly
  "complete",                # G: no key evidence/events omitted; FAERS coverage ok
  "calibrated",              # H: evidence_level / confidence_level appropriate
  # citation-level counts (can also be derived from the sources sheet)
  "n_sources", "n_fabricated", "n_metadata_errors", "n_unsupported", "n_irrelevant",
  "failure_tags", "comments"
)

review <- as.data.table(matrix(NA_character_, nrow = nrow(seed),
                               ncol = length(review_cols),
                               dimnames = list(NULL, review_cols)))
if (nrow(seed)) {
  review[, triplet_id := seed$triplet_id]
  review[, dossier := seed$dossier]
}

################################################################################
# 3. sources sheet (one row per cited source, seeded from the input workbook)
################################################################################
# Pre-fill triplet_id + citation + pmid_or_doi from the sources sheet of the workbook  
# source_n is a per-triplet running index in the workbook's order
# Rows added by hand still get the same dropdowns and a triplet_id picklist

sources_cols <- c(
  "triplet_id", "source_n", "citation", "pmid_or_doi",
  "exists",         # yes/no  : the citation resolves to a real reference
  "metadata_ok",    # yes/partial/no : authors/year/journal/title match
  "supports_claim", # yes/partial/no : the source states what the agent claims
  "relevant",       # yes/no  : pertinent to this pair+event (not filler)
  "pediatric",      # yes/no/na : provides direct pediatric evidence
  "note"
)

input_sources <- as.data.table(read.xlsx(input_xlsx, sheet = "sources"))
input_sources <- input_sources[!is.na(triplet_id) & nzchar(triplet_id)]
# Number the sources within each triplet in their workbook order
# order the sheet by (triplet_id, source_n) so it lines up with the review worklist.
input_sources[, source_n := rowid(triplet_id)]

sources <- data.table(
  triplet_id     = input_sources$triplet_id,
  source_n       = as.character(input_sources$source_n),
  citation       = input_sources$citation,
  pmid_or_doi    = input_sources$PMID_or_DOI,
  exists         = NA_character_,
  metadata_ok    = NA_character_,
  supports_claim = NA_character_,
  relevant       = NA_character_,
  pediatric      = NA_character_,
  note           = NA_character_
)
setcolorder(sources, sources_cols)
setorder(sources, triplet_id, source_n)

################################################################################
# 4. legend sheet
################################################################################

legend <- data.table(
  field = c(
    "(scale) yes/partial/no",
    "(scale) yes/no",
    "verdict",
    "review.* rubric columns",
    "review.n_*",
    "sources sheet",
    "failure_tags",
    "rubric definitions & literature"
  ),
  meaning = c(
    "partial = verifiable incomplete compliance (e.g. 2 of 3 claims supported); justify in comments when not 'yes'",
    "binary check",
    "accepted_as_is / accepted_with_edits / rejected (verdict is independent of performance)",
    "agent-performance rubric per triplet; A=source veracity, B=support/faithfulness, C=relevance/precision, D=pediatric, E=mechanism, F=mapping, G=completeness, H=calibration",
    "citation-level counts; aggregate as sum(n_fabricated)/sum(n_sources) = fabrication rate, etc. (ALCE-style citation precision/recall over the set)",
    "one row per cited source; keyed by triplet_id; lets you compute citation-level rates without trusting the per-triplet counts",
    paste(
      "comma list from: fabricated_source, wrong_metadata, unsupported_claim,",
      "irrelevant_source, weak_evidence_used, adult_only_evidence, wrong_mechanism,",
      "single_drug_attribution, wrong_meddra_mapping, wrong_atc, omitted_evidence,",
      "faers_error, miscalibrated_confidence, misclassified_negative, none"
    ),
    "see REVIEW_GUIDE.md (ALCE citation precision/recall; AIS; Tam 2024 npj Digit Med; BMC Med Inform Decis Mak 2025)"
  )
)

################################################################################
# 5. Assemble the workbook with dropdowns
################################################################################

wb <- createWorkbook()
addWorksheet(wb, "review")
addWorksheet(wb, "sources")
addWorksheet(wb, "legend")

writeData(wb, "review", review, withFilter = TRUE)
writeData(wb, "sources", sources, withFilter = TRUE)
writeData(wb, "legend", legend, withFilter = TRUE)

# Dropdowns applied to a buffer of rows past the seeded data 
# newly added rows keep them
buffer_rows <- 100L
validation_rows <- 2:(max(nrow(review), nrow(sources)) + buffer_rows + 1L)
yn <- '"yes,no"'
ypn <- '"yes,partial,no"'
ynna <- '"yes,no,na"'

add_list <- function(sheet, col, value) {
  dataValidation(wb, sheet, cols = col, rows = validation_rows,
                 type = "list", value = value, allowBlank = TRUE)
}

# review sheet dropdowns (indices match review_cols order)
add_list("review", 6L, '"accepted_as_is,accepted_with_edits,rejected"')  # verdict
add_list("review", 8L, yn)    # sources_all_real
for (col in c(9L, 10L, 11L, 12L, 14L, 15L, 16L, 17L, 18L)) {
  add_list("review", col, ypn)
}
add_list("review", 13L, yn)   # pediatric_valid

# sources sheet dropdowns
# triplet_id picklist points at the review worklist column
# hand-added source rows can only reference a triplet that exists in the review sheet.
if (nrow(review)) {
  dataValidation(wb, "sources", cols = 1L, rows = validation_rows, type = "list",
                 value = sprintf("'review'!$A$2:$A$%d", nrow(review) + 1L),
                 allowBlank = TRUE)
}
add_list("sources", 5L, yn)    # exists
add_list("sources", 6L, ypn)   # metadata_ok
add_list("sources", 7L, ypn)   # supports_claim
add_list("sources", 8L, yn)    # relevant
add_list("sources", 9L, ynna)  # pediatric

freezePane(wb, "review", firstRow = TRUE)
freezePane(wb, "sources", firstRow = TRUE)
setColWidths(wb, "review", cols = 1:length(review_cols), widths = "auto")
setColWidths(wb, "legend", cols = 1:2, widths = c(28, 90))

saveWorkbook(wb, output_xlsx, overwrite = TRUE)

cat("\nReview workbook:", output_xlsx, "\n")
cat("Seeded triplets (worklist rows):", nrow(seed), "\n")
if (nrow(seed)) cat("  ", paste(seed$triplet_id, collapse = ", "), "\n")
missing_dossier <- seed[is.na(dossier), triplet_id]
if (length(missing_dossier)) {
  cat("Triplets without a dossier found under agent/workspace/:",
      paste(missing_dossier, collapse = ", "), "\n")
}
