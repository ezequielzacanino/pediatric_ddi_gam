################################################################################
# Shared functions 
#
# Script 00_functions
################################################################################

################################################################################
# Configuration
################################################################################

library(pacman)
pacman::p_load(data.table, openxlsx)

# Shared OMOP vocabulary root at the workspace level; faers_parsing and gam_benchmark reference the same copy.
vocabulary_dir <- "../data/vocabulary/vocabulary_SNOMED_MEDDRA_RxNorm_ATC"
ruta_concept <- file.path(vocabulary_dir, "CONCEPT.csv")
ruta_concept_ancestor <- file.path(vocabulary_dir, "CONCEPT_ANCESTOR.csv")
# CONCEPT_RELATIONSHIP is only needed to translate MedDRA LLT entries to their PT
ruta_concept_relationship <- file.path(vocabulary_dir, "CONCEPT_RELATIONSHIP.csv")

niveles_nichd <- c(
  "term_neonatal", "infancy", "toddler", "early_childhood",
  "middle_childhood", "early_adolescence", "late_adolescence"
)

# Controlled vocabulary for the control type
niveles_control <- c("positive", "negative")

# Controlled evidence levels, ordered from strongest to weakest
# For a negative control the same tier describes the strength of the source documenting the *absence* of an interaction. 
# The tiers mirror the BNF / Micromedex evidence taxonomy.
evidence_levels <- c(
  "regulatory_label",          # FDA/EMA label or SmPC recognises (or omits) the interaction
  "controlled_study_or_meta",  # RCT, meta-analysis, systematic review, population/PK-modelling study
  "observational_study",       # retrospective / clinical / PK cohort study
  "case_series",               # case series
  "single_case_report",        # single case report
  "theoretical"                # mechanism / class extrapolation only, no observed case
)

################################################################################
# Workbook I/O
################################################################################

# Read a sheet from the curator workbook as a data.table.
read_workbook_sheet <- function(path, sheet) {
  dt <- as.data.table(read.xlsx(path, sheet = sheet))
  char_cols <- names(dt)[vapply(dt, is.character, logical(1))]
  for (col in char_cols) {
    set(dt, j = col, value = sub('^xml:space="preserve">', "", dt[[col]]))
  }
  dt[]
}

################################################################################
# Vocabulary-based MedDRA mapping
################################################################################

# Dictionary layer used by the curation script to resolve MedDRA PTs and roll them up to the HLT/HLGT identifiers

normalize_vocabulary_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)
  x[x %in% c("", "na")] <- NA_character_
  x
}

build_meddra_hierarchy_map <- function(
  rollup_level = "HLT",
  concept_path = ruta_concept,
  ancestor_path = ruta_concept_ancestor
) {

  rollup_class <- switch(
    rollup_level,
    "PT" = "pt",
    "HLT" = "hlt",
    "HLGT" = "hlgt",
    stop(sprintf("Unknown MedDRA roll-up level: %s", rollup_level))
  )

  concept_dt <- fread(
    concept_path,
    quote = "",
    select = c("concept_id", "concept_name", "vocabulary_id", "concept_class_id", "invalid_reason")
  )
  concept_dt[, `:=`(
    concept_id = as.character(concept_id),
    vocabulary_id_key = tolower(trimws(vocabulary_id)),
    concept_class_id_key = tolower(trimws(concept_class_id)),
    invalid_reason_key = fifelse(is.na(invalid_reason), "", trimws(invalid_reason)),
    concept_name_key = normalize_vocabulary_key(concept_name)
  )]

  meddra_concepts <- concept_dt[
    vocabulary_id_key == "meddra" & invalid_reason_key == ""
  ]

  pt_dt <- meddra_concepts[
    concept_class_id_key == "pt",
    .(
      meddra_pt = concept_name,
      meddra_pt_key = concept_name_key,
      meddra_concept_id = concept_id
    )
  ]

  if (rollup_level == "PT") {
    pt_dt[, `:=`(
      meddra_concept_id_2 = NA_character_,
      meddra_concept_id_3 = NA_character_,
      rollup_id = meddra_concept_id
    )]
    return(pt_dt[])
  }

  rollup_dt <- meddra_concepts[
    concept_class_id_key == rollup_class,
    .(
      rollup_id = concept_id,
      rollup_name = concept_name
    )
  ]

  ancestor_dt <- fread(
    ancestor_path,
    select = c(
      "ancestor_concept_id", "descendant_concept_id",
      "min_levels_of_separation", "max_levels_of_separation"
    )
  )
  ancestor_dt[, `:=`(
    ancestor_concept_id = as.character(ancestor_concept_id),
    descendant_concept_id = as.character(descendant_concept_id)
  )]

  path_dt <- merge(
    ancestor_dt,
    rollup_dt,
    by.x = "ancestor_concept_id",
    by.y = "rollup_id"
  )
  path_dt <- merge(
    path_dt,
    pt_dt,
    by.x = "descendant_concept_id",
    by.y = "meddra_concept_id"
  )

  path_dt[, ancestor_concept_id_num := as.numeric(ancestor_concept_id)]
  setorder(path_dt, meddra_pt_key, min_levels_of_separation, max_levels_of_separation, ancestor_concept_id_num)
  path_dt <- path_dt[, .SD[1L], by = meddra_pt_key]

  out <- path_dt[, .(
    meddra_pt,
    meddra_pt_key,
    meddra_concept_id = descendant_concept_id,
    rollup_id = ancestor_concept_id
  )]

  if (rollup_level == "HLT") {
    out[, `:=`(meddra_concept_id_2 = rollup_id, meddra_concept_id_3 = NA_character_)]
  } else {
    out[, `:=`(meddra_concept_id_2 = NA_character_, meddra_concept_id_3 = rollup_id)]
  }

  return(out[])
}

build_benchmark_meddra_map <- function(
  concept_path = ruta_concept,
  ancestor_path = ruta_concept_ancestor
) {
  pt_map <- build_meddra_hierarchy_map("PT", concept_path, ancestor_path)
  hlt_map <- build_meddra_hierarchy_map("HLT", concept_path, ancestor_path)[
    , .(meddra_pt_key, meddra_concept_id_2 = rollup_id)
  ]
  hlgt_map <- build_meddra_hierarchy_map("HLGT", concept_path, ancestor_path)[
    , .(meddra_pt_key, meddra_concept_id_3 = rollup_id)
  ]

  out <- merge(pt_map[, .(meddra_pt, meddra_pt_key, meddra_concept_id)], hlt_map, by = "meddra_pt_key", all.x = TRUE)
  out <- merge(out, hlgt_map, by = "meddra_pt_key", all.x = TRUE)
  out[]
}

################################################################################
# Spreadsheet-driven curation helpers
################################################################################

# script 00 uses build_vocabulary_picklists() to fill the dropdown lists
# curate_pediatric_ddi_reference_set uses resolve_meddra_event_levels() to turn the event at its finest MedDRA level into the PT/HLT/HLGT

# Build the vocabulary-backed dropdown lists used in the input template. 
# Only concepts that participate in the hierarchy are offered
# every term can be rolled up
build_vocabulary_picklists <- function(concept_path = ruta_concept) {
  concept_dt <- fread(
    concept_path,
    quote = "",
    select = c("concept_name", "vocabulary_id", "concept_class_id", "standard_concept", "invalid_reason")
  )
  concept_dt[, `:=`(
    voc = tolower(trimws(vocabulary_id)),
    cls = toupper(trimws(concept_class_id)),
    std = trimws(fifelse(is.na(standard_concept), "", standard_concept)),
    inv = trimws(fifelse(is.na(invalid_reason), "", invalid_reason))
  )]
  valid <- concept_dt[inv == ""]

  # Classification PT/HLT/HLGT are the nodes present in CONCEPT_ANCESTOR; 
  # LLT is taken in full because each LLT reaches a PT through the "Is a" relationship.
  meddra_terms <- function(class_id, classification_only = TRUE) {
    sub <- valid[voc == "meddra" & cls == class_id]
    if (classification_only) sub <- sub[std == "C"]
    sort(unique(sub$concept_name))
  }

  list(
    atc  = sort(unique(valid[voc == "atc" & cls == "ATC 5TH", concept_name])),
    llt  = meddra_terms("LLT", classification_only = FALSE),
    pt   = meddra_terms("PT"),
    hlt  = meddra_terms("HLT"),
    hlgt = meddra_terms("HLGT")
  )
}

# Map a free-text source_type descriptor to a controlled evidence_level
derive_evidence_level <- function(source_type) {
  s <- tolower(trimws(as.character(source_type)))
  out <- rep(NA_character_, length(s))
  out[is.na(out) & grepl("label", s)] <- "regulatory_label"
  out[is.na(out) & grepl("meta-analysis|systematic review|modelling|population pharmacokinetic|randomi", s)] <- "controlled_study_or_meta"
  out[is.na(out) & grepl("retrospective|clinical study|cohort|observational|pharmacokinetic study", s)] <- "observational_study"
  out[is.na(out) & grepl("case series", s)] <- "case_series"
  out[is.na(out) & grepl("case report", s)] <- "single_case_report"
  out[is.na(out) & grepl("theoretical|predicted", s)] <- "theoretical"
  out
}

# Roll the per-citation source_type up to a single controlled evidence_level per triplet 
# keeps the strongest tier across all of the triplet's sources. 
# When the workbook carries no explicit evidence_level
derive_triplet_evidence_level <- function(sources_dt) {
  src <- sources_dt[, .(triplet_id, source_type)]
  src[, level := derive_evidence_level(source_type)]
  src[, rank := match(level, evidence_levels)]
  src <- src[!is.na(rank)]
  src[, .(evidence_level = evidence_levels[min(rank)]), by = triplet_id]
}

# Translate MedDRA LLT concept_ids to their PT concept_id through the "Is a" relationship. 
map_llt_to_pt <- function(llt_ids, meddra_concepts, relationship_path = ruta_concept_relationship) {
  rel <- fread(
    relationship_path,
    quote = "",
    select = c("concept_id_1", "concept_id_2", "relationship_id", "invalid_reason")
  )
  rel[, `:=`(
    concept_id_1 = as.character(concept_id_1),
    concept_id_2 = as.character(concept_id_2)
  )]
  pt_ids <- meddra_concepts[cls == "PT", concept_id]
  out <- rel[
    relationship_id == "Is a" &
      (is.na(invalid_reason) | trimws(invalid_reason) == "") &
      concept_id_1 %in% llt_ids &
      concept_id_2 %in% pt_ids,
    .(llt_id = concept_id_1, pt_id = concept_id_2)
  ]
  # Deterministic tie-break: if an LLT ever resolved to more than one PT 
  # keep the lowest pt_id so the result does not depend on CONCEPT_RELATIONSHIP row order.
  setorder(out, llt_id, pt_id)
  unique(out, by = "llt_id")
}

# Resolve the adverse event of each triplet into PT/HLT/HLGT identifiers.
# For every triplet the finest non-empty level wins
resolve_meddra_event_levels <- function(
  events_dt,
  concept_path = ruta_concept,
  ancestor_path = ruta_concept_ancestor,
  relationship_path = ruta_concept_relationship
) {
  # 1. Pick the finest MedDRA level filled for each triplet.
  ev <- copy(events_dt)
  ev[, `:=`(entered_level = NA_character_, entered_name = NA_character_)]
  ev[!is.na(event_hlgt) & nzchar(trimws(event_hlgt)), `:=`(entered_level = "HLGT", entered_name = trimws(event_hlgt))]
  ev[!is.na(event_hlt)  & nzchar(trimws(event_hlt)),  `:=`(entered_level = "HLT",  entered_name = trimws(event_hlt))]
  ev[!is.na(event_pt)   & nzchar(trimws(event_pt)),   `:=`(entered_level = "PT",   entered_name = trimws(event_pt))]
  ev[!is.na(event_llt)  & nzchar(trimws(event_llt)),  `:=`(entered_level = "LLT",  entered_name = trimws(event_llt))]

  if (anyNA(ev$entered_level)) {
    print(ev[is.na(entered_level), .(triplet_id)])
    stop("There are triplets without any MedDRA level filled in the workbook.")
  }

  # 2. Resolve each entered term to its concept_id within its own class.
  concept_dt <- fread(
    concept_path,
    quote = "",
    select = c("concept_id", "concept_name", "vocabulary_id", "concept_class_id", "invalid_reason")
  )
  concept_dt[, `:=`(
    concept_id = as.character(concept_id),
    voc = tolower(trimws(vocabulary_id)),
    cls = toupper(trimws(concept_class_id)),
    inv = trimws(fifelse(is.na(invalid_reason), "", invalid_reason))
  )]
  meddra <- concept_dt[voc == "meddra" & inv == ""]
  meddra[, name_key := normalize_vocabulary_key(concept_name)]

  name_map <- unique(meddra[, .(cls, name_key, concept_id, concept_name)])
  ev[, name_key := normalize_vocabulary_key(entered_name)]
  ev <- merge(
    ev, name_map,
    by.x = c("entered_level", "name_key"), by.y = c("cls", "name_key"),
    all.x = TRUE, sort = FALSE
  )
  setnames(ev, "concept_id", "entered_id")

  # The name-based merge must stay 1:1 per triplet
  # guards against two MedDRA concepts collapsing to the same normalized name_key and multiplying a row.
  stopifnot(uniqueN(ev$triplet_id) == nrow(ev))

  unresolved <- ev[is.na(entered_id), .(triplet_id, entered_level, entered_name)]
  if (nrow(unresolved) > 0) {
    print(unresolved)
    stop("Some MedDRA events were not found in the vocabulary at the stated level.")
  }

  # 3. LLT entries are non-standard and absent from CONCEPT_ANCESTOR; 
  # translates them to their PT so the same ancestor-based roll-up applies to everything.
  ev[, `:=`(effective_id = entered_id, effective_level = entered_level)]
  llt_ids <- ev[entered_level == "LLT", unique(entered_id)]
  if (length(llt_ids) > 0) {
    llt_pt <- map_llt_to_pt(llt_ids, meddra, relationship_path)
    ev <- merge(ev, llt_pt, by.x = "entered_id", by.y = "llt_id", all.x = TRUE, sort = FALSE)
    missing_pt <- ev[entered_level == "LLT" & is.na(pt_id), .(triplet_id, entered_name)]
    if (nrow(missing_pt) > 0) {
      print(missing_pt)
      stop("Some LLT terms have no associated PT through the 'Is a' relationship.")
    }
    ev[entered_level == "LLT", `:=`(effective_id = pt_id, effective_level = "PT")]
    ev[, pt_id := NULL]
  }

  # 4. Roll the effective node up to PT/HLT/HLGT
  effective_ids <- ev[, unique(effective_id)]
  ancestor_dt <- fread(
    ancestor_path,
    select = c("ancestor_concept_id", "descendant_concept_id", "min_levels_of_separation", "max_levels_of_separation")
  )
  ancestor_dt[, `:=`(
    ancestor_concept_id = as.character(ancestor_concept_id),
    descendant_concept_id = as.character(descendant_concept_id)
  )]
  anc <- ancestor_dt[descendant_concept_id %in% effective_ids]
  anc <- merge(anc, meddra[, .(concept_id, cls, concept_name)], by.x = "ancestor_concept_id", by.y = "concept_id")
  anc <- anc[cls %in% c("PT", "HLT", "HLGT")]

  # Deterministic primary path per (effective node, target class).
  anc[, ancestor_concept_id_num := as.numeric(ancestor_concept_id)]
  setorder(anc, descendant_concept_id, cls, min_levels_of_separation, max_levels_of_separation, ancestor_concept_id_num)
  primary <- anc[, .SD[1L], by = .(descendant_concept_id, cls)]

  pt_lvl   <- primary[cls == "PT",   .(effective_id = descendant_concept_id, meddra_pt = concept_name, meddra_concept_id = ancestor_concept_id)]
  hlt_lvl  <- primary[cls == "HLT",  .(effective_id = descendant_concept_id, meddra_concept_id_2 = ancestor_concept_id)]
  hlgt_lvl <- primary[cls == "HLGT", .(effective_id = descendant_concept_id, meddra_concept_id_3 = ancestor_concept_id)]

  out <- ev[, .(triplet_id, entered_level, entered_name, effective_id)]
  out <- merge(out, pt_lvl,   by = "effective_id", all.x = TRUE, sort = FALSE)
  out <- merge(out, hlt_lvl,  by = "effective_id", all.x = TRUE, sort = FALSE)
  out <- merge(out, hlgt_lvl, by = "effective_id", all.x = TRUE, sort = FALSE)
  out[, effective_id := NULL]
  setorder(out, triplet_id)
  out[]
}

################################################################################
# FAERS triplet co-reporting 
################################################################################

# Distinct pediatric FAERS cases that co-report a drug-drug-event triplet
#
# The event is matched at MedDRA PT (`meddra_concept_id`), the finest level. 
# >= 1 at PT guarantees the triplet is detectable at every coarser benchmark roll-up too. 
count_triplet_coreports <- function(ade, drug1_id, drug2_id, event_id) {
  drug1_id <- as.integer(drug1_id)
  drug2_id <- as.integer(drug2_id)
  event_id <- as.integer(event_id)
  reports_drug1 <- ade[atc_concept_id == drug1_id, unique(safetyreportid)]
  reports_drug2 <- ade[atc_concept_id == drug2_id, unique(safetyreportid)]
  reports_event <- ade[meddra_concept_id == event_id, unique(safetyreportid)]
  triplet_reports <- intersect(intersect(reports_drug1, reports_drug2), reports_event)
  # Per-stage counts, guaranteeing one row per canonical NICHD stage in order.
  stage_reports <- unique(ade[safetyreportid %in% triplet_reports, .(safetyreportid, nichd)])
  per_stage <- data.table(nichd = niveles_nichd)
  per_stage[stage_reports[, .(coreports = uniqueN(safetyreportid)), by = nichd],
            on = "nichd", coreports := i.coreports]
  per_stage[is.na(coreports), coreports := 0L]
  list(total = length(triplet_reports), per_stage = per_stage)
}
