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
#   Rscript agent/tools/faers_triplet_coreport.R \
#     --drug1-id <atc_concept_id> --drug2-id <atc_concept_id> --event-id <meddra_pt_concept_id>
# By name (resolves through the OMOP vocabulary exactly like the curation script, so ids
# match the ones the pipeline assigns):
#   Rscript agent/tools/faers_triplet_coreport.R \
#     --drug1 "methotrexate; systemic" --drug2 "trimethoprim; systemic" \
#     --event-pt "Bone marrow depression"
# The event flag is --event-id (a PT concept id) or one of
# --event-llt/--event-pt/--event-hlt/--event-hlgt (any level; resolved to its PT).
#
# --mode rank (only the two drugs, no event): ranked table of the PT events co-reported with
# the pair in pediatric FAERS, with triplet_coreports (pair + that event in the same case) and
# single_drug_event_max per event. Options:
# --top <n> (default 25), --min-triplet-coreports <n> (default 1).
#   Rscript agent/tools/faers_triplet_coreport.R --mode rank \
#     --drug1 "vincristine; parenteral" --drug2 "itraconazole; systemic"

source("00_functions.R", local = TRUE)

# FAERS case-level table produced by faers_parsing, read live (relative path, not versioned).
ade_raw_file <- "../faers_parsing/data/processed/ade_raw.csv"

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

# Output mode: "triplet" (one drug-drug-event query) or "rank" (drug-drug events table).
mode <- if (!is.na(getopt("mode"))) tolower(getopt("mode")) else "triplet"

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
  # Same resolver the curation script uses: rolls the entered level up to its PT concept id.
  as.integer(resolve_meddra_event_levels(ev)$meddra_concept_id)
}

drug1_id <- resolve_drug_id("drug1-id", "drug1")
drug2_id <- resolve_drug_id("drug2-id", "drug2")
event_id <- if (mode == "triplet") resolve_event_id() else NA_integer_

################################################################################
# 3. Count triplet co-reports over the pediatric FAERS table (PT level)
################################################################################

ade <- fread(ade_raw_file,
             select = c("safetyreportid", "atc_concept_id", "meddra_concept_id", "nichd"))
ade <- ade[nichd %in% niveles_nichd]

# --- rank mode: PT events co-reported with the pair, ranked by detectability -------------
if (mode == "rank") {
  # PT name lookup for the events co-reported with the pair.
  if (!exists("concept_dt")) {
    concept_dt <- fread(
      ruta_concept, quote = "",
      select = c("concept_id", "concept_name", "vocabulary_id", "concept_class_id", "invalid_reason")
    )
  }
  meddra_pt_dict <- unique(concept_dt[
    tolower(trimws(vocabulary_id)) == "meddra" &
      toupper(trimws(concept_class_id)) == "PT" &
      (is.na(invalid_reason) | trimws(invalid_reason) == ""),
    .(meddra_concept_id = as.integer(concept_id), event_pt = concept_name)
  ], by = "meddra_concept_id")

  top_n <- if (!is.na(getopt("top"))) as.integer(getopt("top")) else 25L
  min_co <- if (!is.na(getopt("min-triplet-coreports"))) as.integer(getopt("min-triplet-coreports")) else 1L

  # 1. Triplet co-reports per PT event: within the pair's co-reports, distinct reports
  #    where the pair and that event co-occur (drug1 + drug2 + event).
  reports_drug1 <- ade[atc_concept_id == drug1_id, unique(safetyreportid)]
  reports_drug2 <- ade[atc_concept_id == drug2_id, unique(safetyreportid)]
  pair_reports <- intersect(reports_drug1, reports_drug2)
  pair_events <- unique(ade[safetyreportid %in% pair_reports & !is.na(meddra_concept_id),
                            .(safetyreportid, meddra_concept_id)])
  ranked <- pair_events[, .(triplet_coreports = uniqueN(safetyreportid)), by = meddra_concept_id]

  # 2. single_drug_event_max: distinct pediatric reports where one drug alone co-occurs with
  #    the event, taken as the max over the two drugs (Kontsioti mono-drug stratifier).
  ev_ids <- ranked$meddra_concept_id
  supp1 <- ade[atc_concept_id == drug1_id & meddra_concept_id %in% ev_ids,
               .(s1 = uniqueN(safetyreportid)), by = meddra_concept_id]
  supp2 <- ade[atc_concept_id == drug2_id & meddra_concept_id %in% ev_ids,
               .(s2 = uniqueN(safetyreportid)), by = meddra_concept_id]
  ranked[supp1, on = "meddra_concept_id", s1 := i.s1]
  ranked[supp2, on = "meddra_concept_id", s2 := i.s2]
  ranked[is.na(s1), s1 := 0L]
  ranked[is.na(s2), s2 := 0L]
  ranked[, single_drug_event_max := pmax(s1, s2)]

  # 3. Name, filter and order by detectability.
  ranked[meddra_pt_dict, on = "meddra_concept_id", event_pt := i.event_pt]
  ranked[is.na(event_pt), event_pt := "(PT name unresolved)"]
  ranked <- ranked[triplet_coreports >= min_co][order(-triplet_coreports, event_pt)]

  cat(sprintf("pair: drug1_id=%d  drug2_id=%d  pair_coreports=%d  detectable_event_pts=%d  (PT, triplet_coreports>=%d)\n",
              drug1_id, drug2_id, length(pair_reports), nrow(ranked), min_co))
  cat(sprintf("%-45s %18s %22s\n", "event_pt", "triplet_coreports", "single_drug_event_max"))
  n_show <- min(top_n, nrow(ranked))
  for (r in seq_len(n_show)) {
    cat(sprintf("%-45s %18d %22d\n",
                substr(ranked$event_pt[r], 1L, 45L),
                ranked$triplet_coreports[r], ranked$single_drug_event_max[r]))
  }
  if (nrow(ranked) > n_show) {
    cat(sprintf("... (%d more event(s) with triplet_coreports >= %d)\n", nrow(ranked) - n_show, min_co))
  }
  quit(save = "no", status = 0L)
}

res <- count_triplet_coreports(ade, drug1_id, drug2_id, event_id)

cat(sprintf("triplet: drug1_id=%d  drug2_id=%d  event_id=%d  (MedDRA level: PT)\n",
            drug1_id, drug2_id, event_id))
cat(sprintf("total_triplet_coreports: %d\n", res$total))
for (r in seq_len(nrow(res$per_stage))) {
  cat(sprintf("  %-18s %d\n", res$per_stage$nichd[r], res$per_stage$triplet_coreports[r]))
}
if (res$total == 0L) {
  cat("VERDICT: UNDETECTABLE - 0 triplet co-reports, empty benchmark cell. Reject this triplet.\n")
} else {
  cat("VERDICT: detectable (>= 1 triplet co-report).\n")
}
