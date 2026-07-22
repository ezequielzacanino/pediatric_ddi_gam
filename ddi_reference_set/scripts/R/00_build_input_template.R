################################################################################
# Builds the spreadsheet-driven input template for the curated pediatric set.
#
# Script 00_build_input_template
################################################################################
#
# Produces input/ddi_reference_input.xlsx
# Every drug/event/category column has a vocabulary-backed dropdown
# Re-run this script only to (re)generate the empty workbook
# overwrite_existing <- TRUE to rebuild from scratch.

source("00_functions.R", local = TRUE)
library(openxlsx)

# Guard against destroying manual edits. Flip to TRUE only for a clean rebuild.
overwrite_existing <- TRUE

input_dir <- "./input"
output_xlsx <- file.path(input_dir, "ddi_reference_input.xlsx")
dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)

if (file.exists(output_xlsx) && !overwrite_existing) {
  stop(sprintf(
    "%s already exists. Set overwrite_existing <- TRUE to rebuild empty.",
    output_xlsx
  ))
}

################################################################################
# Empty input schema
################################################################################
# Only the header row plus the vocabulary-backed dropdowns

# triplets sheet schema. event terms are split across four optional MedDRA levels
input_triplets <- data.table(
  triplet_id = character(),
  control_type = character(),
  drug1 = character(),
  drug2 = character(),
  event_llt = character(),
  event_pt = character(),
  event_hlt = character(),
  event_hlgt = character(),
  interaction_type = character(),
  mechanism = character(),
  evidence_type = character(),
  evidence_level = character(),
  pediatric_population = character(),
  age_range = character(),
  ontogenic_modulation = character(),
  higher_risk_stages = character(),
  ontogeny_evidence = character(),
  source_title = character(),
  source_year = character(),
  source_type = character(),
  confidence_level = character(),
  rationale = character(),
  comments = character()
)

input_sources <- data.table(
  triplet_id = character(),
  PMID_or_DOI = character(),
  URL = character(),
  citation = character(),
  source_type = character(),
  notes = character()
)

################################################################################
# Build the workbook with vocabulary-backed dropdowns
################################################################################

picklists <- build_vocabulary_picklists()

wb <- createWorkbook()
addWorksheet(wb, "triplets")
addWorksheet(wb, "sources")

# NICHD dropdown options
nichd_options <- c(
  niveles_nichd,
  "term_neonatal,infancy",
  "term_neonatal,infancy,toddler",
  "infancy,toddler",
  "early_childhood,middle_childhood",
  "early_adolescence,late_adolescence"
)

# Reference sheets hold the dropdown lists
# Each list validation points at the matching reference column.
ref_lists <- list(
  ref_atc = picklists$atc,
  ref_llt = picklists$llt,
  ref_pt = picklists$pt,
  ref_hlt = picklists$hlt,
  ref_hlgt = picklists$hlgt,
  ref_nichd = nichd_options
)
for (ref_name in names(ref_lists)) {
  addWorksheet(wb, ref_name, visible = FALSE)
  writeData(wb, ref_name, data.frame(term = ref_lists[[ref_name]]))
}

writeData(wb, "triplets", input_triplets, withFilter = TRUE)
writeData(wb, "sources", input_sources, withFilter = TRUE)

# Apply list validations to a buffer of empty rows
validation_rows <- 2:202L

range_ref <- function(ref_name) {
  sprintf("'%s'!$A$2:$A$%d", ref_name, length(ref_lists[[ref_name]]) + 1L)
}

# triplets sheet column positions (kept in sync with the input_triplets layout)
col_control <- 2L
col_drug1 <- 3L; col_drug2 <- 4L
col_llt <- 5L; col_pt <- 6L; col_hlt <- 7L; col_hlgt <- 8L
col_interaction <- 9L
col_evidence_level <- 12L
col_ontogeny <- 15L; col_higher_risk_stages <- 16L
col_confidence <- 21L

add_list <- function(col, value) {
  dataValidation(wb, "triplets", cols = col, rows = validation_rows,
                 type = "list", value = value, allowBlank = TRUE)
}

add_list(col_drug1, range_ref("ref_atc"))
add_list(col_drug2, range_ref("ref_atc"))
add_list(col_llt, range_ref("ref_llt"))
add_list(col_pt, range_ref("ref_pt"))
add_list(col_hlt, range_ref("ref_hlt"))
add_list(col_hlgt, range_ref("ref_hlgt"))
add_list(col_higher_risk_stages, range_ref("ref_nichd"))
# Small categorical fields use inline lists. 
# interaction_type carries "none" for negative controls
add_list(col_control, '"positive,negative"')
add_list(col_interaction, '"pharmacokinetic,pharmacodynamic,mixed,unknown,pharmaceutical,none"')
add_list(col_evidence_level, sprintf('"%s"', paste(evidence_levels, collapse = ",")))
add_list(col_ontogeny, '"yes,no,unknown"')
add_list(col_confidence, '"high,moderate"')

freezePane(wb, "triplets", firstRow = TRUE)
freezePane(wb, "sources", firstRow = TRUE)

saveWorkbook(wb, output_xlsx, overwrite = TRUE)

cat("\nInput workbook:", output_xlsx, "\n")
cat("Seed triplets:", nrow(input_triplets), "\n")
cat("Dropdown sizes -> ATC:", length(picklists$atc),
    "| LLT:", length(picklists$llt),
    "| PT:", length(picklists$pt),
    "| HLT:", length(picklists$hlt),
    "| HLGT:", length(picklists$hlgt), "\n")
