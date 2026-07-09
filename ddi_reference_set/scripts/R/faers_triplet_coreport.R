################################################################################
# FAERS triplet co-report lookup (curation detectability gate)
# Script faers_triplet_coreport
################################################################################
#
# On-demand check for the positive/negative curation skills: 
# given a drug-drug-event triplet, report how many triplets reports are in FAERS 
# 
# mirrors the three-set intersection gam_benchmark uses to build the contingency cell

# The event is counted at MedDRA PT (the finest level). 
#
# Run from ddi_reference_set/. By id (no vocabulary read):
#   Rscript scripts/R/faers_triplet_coreport.R \
#     --drug1-id <atc_concept_id> --drug2-id <atc_concept_id> --event-id <meddra_pt_concept_id>
# By name (resolves through the OMOP vocabulary exactly like script 02, so ids
# match the ones the pipeline assigns):
#   Rscript scripts/R/faers_triplet_coreport.R \
#     --drug1 "methotrexate; systemic" --drug2 "trimethoprim; systemic" \
#     --event-pt "Bone marrow depression"
# The event flag is --event-id (a PT concept id) or one of
# --event-llt/--event-pt/--event-hlt/--event-hlgt (any level; resolved to its PT).

source("00_functions.R", local = TRUE)

# FAERS case-level table shared with gam_benchmark (relative path, not versioned).
ade_raw_file <- "../gam_benchmark/data/processed/ade_raw.csv"

################################################################################
# 1. Parse --key value arguments
################################################################################

args <- commandArgs(trailingOnly = TRUE)
opt <- list()
i <- 1L
while (i <= length(args)) {
  key <- args[[i]]
  if (startsWith(key, "--")) {
    val <- if (i + 1L <= length(args)) args[[i + 1L]] else ""
    opt[[sub("^--", "", key)]] <- val
    i <- i + 2L
  } else {
    i <- i + 1L
  }
}
getopt <- function(key) if (!is.null(opt[[key]]) && nzchar(trimws(opt[[key]]))) trimws(opt[[key]]) else NA_character_

################################################################################
# 2. Resolve drug and event ids 
################################################################################

# Drugs: build the ATC 5th dictionary only if a drug was given by name.
need_atc <- is.na(getopt("drug1-id")) || is.na(getopt("drug2-id"))
atc5_dict <- NULL
if (need_atc) {
  concept_dt <- fread(
    ruta_concept, quote = "",
    select = c("concept_id", "concept_name", "vocabulary_id", "concept_class_id", "invalid_reason")
  )
  atc5_dict <- unique(concept_dt[
    tolower(trimws(vocabulary_id)) == "atc" &
      toupper(trimws(concept_class_id)) == "ATC 5TH" &
      (is.na(invalid_reason) | trimws(invalid_reason) == ""),
    .(atc_concept_name = concept_name, atc_concept_id = as.character(concept_id))
  ], by = "atc_concept_name")
}

resolve_drug_id <- function(id_key, name_key) {
  if (!is.na(getopt(id_key))) return(as.integer(getopt(id_key)))
  drug_name <- getopt(name_key)
  if (is.na(drug_name)) stop(sprintf("Provide --%s or --%s.", id_key, name_key))
  idx <- match(drug_name, atc5_dict$atc_concept_name)
  if (is.na(idx)) stop(sprintf("Drug name not found as ATC 5th concept: %s", drug_name))
  as.integer(atc5_dict$atc_concept_id[idx])
}

# Returns the event's PT concept id: the level at which detectability is gated.
resolve_event_id <- function() {
  if (!is.na(getopt("event-id"))) return(as.integer(getopt("event-id")))
  ev <- data.table(
    triplet_id = "Q",
    event_llt = getopt("event-llt"), event_pt = getopt("event-pt"),
    event_hlt = getopt("event-hlt"), event_hlgt = getopt("event-hlgt")
  )
  if (ev[, all(is.na(c(event_llt, event_pt, event_hlt, event_hlgt)))]) {
    stop("Provide --event-id or one of --event-llt/--event-pt/--event-hlt/--event-hlgt.")
  }
  # Same resolver script 02 uses: rolls the entered level up to its PT concept id.
  as.integer(resolve_meddra_event_levels(ev)$meddra_concept_id)
}

drug1_id <- resolve_drug_id("drug1-id", "drug1")
drug2_id <- resolve_drug_id("drug2-id", "drug2")
event_id <- resolve_event_id()

################################################################################
# 3. Count triplet co-reports over the pediatric FAERS table (PT level)
################################################################################

ade <- fread(ade_raw_file,
             select = c("safetyreportid", "atc_concept_id", "meddra_concept_id", "nichd"))
ade <- ade[nichd %in% niveles_nichd]

res <- count_triplet_coreports(ade, drug1_id, drug2_id, event_id)

cat(sprintf("triplet: drug1_id=%d  drug2_id=%d  event_id=%d  (MedDRA level: PT)\n",
            drug1_id, drug2_id, event_id))
cat(sprintf("total_coreports: %d\n", res$total))
for (r in seq_len(nrow(res$per_stage))) {
  cat(sprintf("  %-18s %d\n", res$per_stage$nichd[r], res$per_stage$coreports[r]))
}
if (res$total == 0L) {
  cat("VERDICT: UNDETECTABLE - 0 co-reports, empty benchmark cell. Reject this triplet.\n")
} else {
  cat("VERDICT: detectable (>= 1 co-report).\n")
}
