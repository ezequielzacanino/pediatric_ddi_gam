################################################################################
# Construction script for the pediatric curated reference set
#
# Script curate_pediatric_ddi_reference_set
################################################################################
#
# Reads the curator-maintained workbook (input/ddi_reference_input.xlsx)
# produces the analysis-ready triplets and sources consumed by gam_benchmark.
#
# Consolidation step: runs after every workbook edit, both after positives are
# curated (script 01) and after negatives are curated (script 03).

source("00_functions.R", local = TRUE)
library(openxlsx)

################################################################################
# Configuration
################################################################################

input_xlsx <- "./input/ddi_reference_input.xlsx"

output_dir <- "./results/curated_pediatric_ddi_reference_set/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

triplets_file <- file.path(output_dir, "curated_pediatric_ddi_triplets.csv")
sources_file <- file.path(output_dir, "curated_pediatric_ddi_sources.csv")

if (!file.exists(input_xlsx)) {
  stop(sprintf("%s was not found. Run scripts/R/00_build_input_template.R first.", input_xlsx))
}
if (!file.exists(ruta_concept)) {
  stop("CONCEPT.csv was not found in data/vocabulary/vocabulary_SNOMED_MEDDRA_RxNorm_ATC.")
}

################################################################################
# Read workbook
################################################################################

triplets <- read_workbook_sheet(input_xlsx, "triplets")
sources <- read_workbook_sheet(input_xlsx, "sources")

# Ensure the four optional event columns exist
for (event_col in c("event_llt", "event_pt", "event_hlt", "event_hlgt")) {
  if (!event_col %in% names(triplets)) triplets[, (event_col) := NA_character_]
  triplets[, (event_col) := as.character(get(event_col))]
}
triplets[, `:=`(drug1 = trimws(as.character(drug1)), drug2 = trimws(as.character(drug2)))]

# Control type: every triplet is a positive control or a negative 
if (!"control_type" %in% names(triplets)) triplets[, control_type := NA_character_]
triplets[, control_type := trimws(as.character(control_type))]
triplets[is.na(control_type) | !nzchar(control_type), control_type := "positive"]

# Evidence level: use the curator-set controlled value when present
# otherwise derive the strongest tier across the triplet's own source_type descriptor
if (!"evidence_level" %in% names(triplets)) triplets[, evidence_level := NA_character_]
triplets[, evidence_level := trimws(as.character(evidence_level))]
evidence_inputs <- rbindlist(list(
  triplets[, .(triplet_id, source_type)],
  sources[, .(triplet_id, source_type)]
), use.names = TRUE)
triplets[derive_triplet_evidence_level(evidence_inputs), on = "triplet_id", evidence_level_derived := i.evidence_level]
triplets[is.na(evidence_level) | !nzchar(evidence_level), evidence_level := evidence_level_derived]
triplets[, evidence_level_derived := NULL]


# Ontogeny columns: record whether the interaction risk is age-modulated
for (ontogeny_col in c("ontogenic_modulation", "higher_risk_stages", "ontogeny_evidence")) {
  if (!ontogeny_col %in% names(triplets)) triplets[, (ontogeny_col) := NA_character_]
  triplets[, (ontogeny_col) := as.character(get(ontogeny_col))]
}
# Normalize casing/whitespace, then default anything missing or outside the controlled vocabulary to "unknown".
triplets[, ontogenic_modulation := tolower(trimws(ontogenic_modulation))]
triplets[is.na(ontogenic_modulation) | !ontogenic_modulation %in% c("yes", "no"),
         ontogenic_modulation := "unknown"]

# higher_risk_stages is a comma-separated subset of niveles_nichd

# Split the comma-separated list: normalize non-breaking spaces, accepts comma or semicolon as separators. 
# then trim/lowercase each token and drop empties.
split_stage_tokens <- function(x) {
  x <- fifelse(is.na(x), "", as.character(x))
  x <- gsub(" ", " ", x)                 # non-breaking space -> normal space
  toks <- strsplit(x, "\\s*[,;]\\s*")
  lapply(toks, function(t) {
    t <- tolower(trimws(t))
    t[nzchar(t)]
  })
}
stage_tokens <- split_stage_tokens(triplets$higher_risk_stages)
invalid_stage <- vapply(stage_tokens, function(tokens) any(!tokens %in% niveles_nichd), logical(1))
if (any(invalid_stage)) {
  offending <- vapply(stage_tokens, function(tokens) paste(setdiff(tokens, niveles_nichd), collapse = ", "), character(1))
  print(triplets[invalid_stage, .(triplet_id, higher_risk_stages)][, invalid_tokens := offending[invalid_stage]][])
  stop("higher_risk_stages contiene etapas fuera de niveles_nichd.")
}

stages_present <- vapply(stage_tokens, function(tokens) any(nzchar(tokens)), logical(1))
if (any(stages_present & triplets$ontogenic_modulation != "yes")) {
  print(triplets[stages_present & ontogenic_modulation != "yes", .(triplet_id, ontogenic_modulation, higher_risk_stages)])
  stop("higher_risk_stages solo debe completarse cuando ontogenic_modulation = 'yes'.")
}

################################################################################
# Dictionary mapping
################################################################################

# 1. Drug coding . exact join on ATC 5th concept_name selected from the dropdown.
concept_dt <- fread(
  ruta_concept,
  quote = "",
  select = c("concept_id", "concept_name", "vocabulary_id", "concept_class_id", "concept_code", "invalid_reason")
)
atc5_dict <- unique(concept_dt[
  tolower(trimws(vocabulary_id)) == "atc" &
    toupper(trimws(concept_class_id)) == "ATC 5TH" &
    (is.na(invalid_reason) | trimws(invalid_reason) == ""),
  .(
    atc_concept_name = concept_name,
    atc_concept_code = concept_code,
    atc_concept_id = as.character(concept_id)
  )
], by = "atc_concept_name")

resolve_drug <- function(drug_names) {
  idx <- match(drug_names, atc5_dict$atc_concept_name)
  data.table(
    atc_concept_name = atc5_dict$atc_concept_name[idx],
    atc_concept_code = atc5_dict$atc_concept_code[idx],
    atc_concept_id = atc5_dict$atc_concept_id[idx]
  )
}

drug1_res <- resolve_drug(triplets$drug1)
drug2_res <- resolve_drug(triplets$drug2)

# Substance label kept for readability (drops the ATC route suffix).
triplets[, `:=`(
  drug1_name = trimws(sub(";.*", "", drug1)),
  drug1_atc = drug1_res$atc_concept_code,
  drug1_atc_original = drug1_res$atc_concept_code,
  drug1_atc_concept_id = drug1_res$atc_concept_id,
  drug1_atc_concept_name = drug1_res$atc_concept_name,
  drug1_mapping_source = "vocabulary_dropdown",
  drug2_name = trimws(sub(";.*", "", drug2)),
  drug2_atc = drug2_res$atc_concept_code,
  drug2_atc_original = drug2_res$atc_concept_code,
  drug2_atc_concept_id = drug2_res$atc_concept_id,
  drug2_atc_concept_name = drug2_res$atc_concept_name,
  drug2_mapping_source = "vocabulary_dropdown"
)]

unmapped_drugs_dt <- unique(rbindlist(list(
  triplets[is.na(drug1_atc_concept_id), .(triplet_id, drug_role = "drug1", drug = drug1)],
  triplets[is.na(drug2_atc_concept_id), .(triplet_id, drug_role = "drug2", drug = drug2)]
)))
if (nrow(unmapped_drugs_dt) > 0) {
  print(unmapped_drugs_dt)
  stop("Hay farmacos cuyo nombre no coincide con un ATC 5th del vocabulario.")
}

# 2. MedDRA event coding : helper resolves the finest available level upward to PT/HLT/HLGT as needed.
event_levels <- resolve_meddra_event_levels(
  triplets[, .(triplet_id, event_llt, event_pt, event_hlt, event_hlgt)]
)
triplets <- merge(
  triplets,
  event_levels[, .(triplet_id, meddra_pt, meddra_concept_id, meddra_concept_id_2, meddra_concept_id_3)],
  by = "triplet_id",
  all.x = TRUE
)

# 3. Event MedDRA code . annotate each PT with its MedDRA code for traceability.
meddra_code_map <- unique(concept_dt[
  tolower(trimws(vocabulary_id)) == "meddra",
  .(meddra_concept_id = as.character(concept_id), meddra_code = as.character(concept_code))
])
triplets <- merge(triplets, meddra_code_map, by = "meddra_concept_id", all.x = TRUE)

setcolorder(triplets, c("triplet_id", "control_type",
                        "drug1_name", "drug1_atc_original", "drug1_atc", "drug1_atc_concept_id",
                        "drug1_atc_concept_name", "drug1_mapping_source",
                        "drug2_name", "drug2_atc_original", "drug2_atc", "drug2_atc_concept_id",
                        "drug2_atc_concept_name", "drug2_mapping_source",
                        "meddra_pt", "meddra_code", "meddra_concept_id", "meddra_concept_id_2", "meddra_concept_id_3",
                        "interaction_type", "mechanism", "evidence_type", "evidence_level",
                        "pediatric_population", "age_range",
                        "ontogenic_modulation", "higher_risk_stages", "ontogeny_evidence",
                        "source_title",
                        "source_year", "source_type", "confidence_level",
                        "rationale", "comments"))
# Drop the raw dropdown columns now that drugs and events are resolved.
triplets[, c("drug1", "drug2", "event_llt", "event_pt", "event_hlt", "event_hlgt") := NULL]

################################################################################
# Quality control
################################################################################

# 1. Semantic duplicates. no two triplets may share the same resolved drug pair and event.
dup_triplets <- triplets[duplicated(triplets[, .(drug1_atc, drug2_atc, meddra_pt, meddra_concept_id_2, meddra_concept_id_3)])]
if (nrow(dup_triplets) > 0) {
  print(dup_triplets[, .(triplet_id, drug1_name, drug2_name, meddra_pt)])
  stop("Existen tripletes semanticos duplicados. Revisar la planilla de entrada.")
}

# 2. Source integrity. every triplet must have at least one source row.
source_counts <- sources[, .N, by = triplet_id]
triplets <- merge(triplets, source_counts, by = "triplet_id", all.x = TRUE)
stopifnot(all(!is.na(triplets$N) & triplets$N >= 1L))
triplets[, N := NULL]

# 3. Categorical constraints. all controlled-vocabulary fields must use valid levels.
stopifnot(all(triplets$control_type %in% niveles_control))
stopifnot(all(triplets$confidence_level %in% c("high", "moderate")))
stopifnot(all(triplets$interaction_type %in% c("pharmacokinetic", "pharmacodynamic", "mixed", "unknown", "pharmaceutical", "none")))

# Coherence: negative controls must have interaction_type = "none"; positive controls must declare one.
if (any(triplets$control_type == "negative" & triplets$interaction_type != "none")) {
  print(triplets[control_type == "negative" & interaction_type != "none", .(triplet_id, control_type, interaction_type)])
  stop("Los controles negativos deben tener interaction_type = 'none'.")
}
if (any(triplets$control_type == "positive" & triplets$interaction_type == "none")) {
  print(triplets[control_type == "positive" & interaction_type == "none", .(triplet_id, control_type, interaction_type)])
  stop("Los controles positivos no pueden tener interaction_type = 'none'.")
}
# Negative controls cannot claim age-modulated interaction risk.
if (any(triplets$control_type == "negative" & triplets$ontogenic_modulation == "yes")) {
  print(triplets[control_type == "negative" & ontogenic_modulation == "yes", .(triplet_id)])
  stop("Los controles negativos no pueden tener ontogenic_modulation = 'yes'.")
}

# 4. Referential integrity. every source row must link to a known triplet_id with a non-empty URL.
stopifnot(all(sources$triplet_id %in% triplets$triplet_id))
stopifnot(all(nchar(sources$URL) > 0))

# 5. Stable sort. deterministic row order for reproducible diffs.
setorder(triplets, triplet_id)
setorder(sources, triplet_id, citation)

################################################################################
# Results
################################################################################

fwrite(triplets, triplets_file)
fwrite(sources, sources_file)

summary_evidence <- triplets[, .N, by = evidence_type][order(-N, evidence_type)]
summary_mechanism <- triplets[, .(N = .N), by = .(interaction_type)][order(-N, interaction_type)]

cat("\nTriplets file:", triplets_file, "\n")
cat("Sources file:", sources_file, "\n")
cat("Curated triplets:", nrow(triplets), "\n")
cat("Positive controls:", sum(triplets$control_type == "positive"), "\n")
cat("Negative controls:", sum(triplets$control_type == "negative"), "\n")
cat("Unique pairs:", uniqueN(triplets[, .(drug1_name, drug2_name)]), "\n")
