################################################################################
# Functions script
# Script: 00_functions.R
# To use: source("00_functions.R", local = TRUE)
################################################################################

################################################################################
# General pipeline configuration
################################################################################

set.seed(9427)

library(pacman)
pacman::p_load(data.table, pbapply, parallel, doParallel, foreach, doRNG,
  mgcv, MASS, akima, pROC, ggplot2, scales, svglite)

# NICHD stage level ordering used throughout all scripts
niveles_nichd <- c(
  "term_neonatal", "infancy", "toddler", "early_childhood",
  "middle_childhood", "early_adolescence", "late_adolescence"
)

# Reserve 25% of cores for the OS / other processes
n_cores <- max(1, floor(detectCores() * 0.25))

# Cross-project paths: ade_raw is produced by faers_parsing; vocabulary is shared at the workspace root.
ruta_ade_raw <- "../faers_parsing/data/processed/ade_raw.csv"
vocabulary_dir <- "../data/vocabulary/vocabulary_SNOMED_MEDDRA_RxNorm_ATC"
ruta_concept <- file.path(vocabulary_dir, "CONCEPT.csv")

# GAM formula parameters — defaults used throughout the pipeline
spline_individuales <- TRUE
include_sex <- FALSE
include_stage_sex <- FALSE
k_spline <- 7
include_nichd <- FALSE
nichd_spline <- FALSE
bs_type <- "cs"
select <- FALSE
method <- "fREML"

# Compact encoding of the active formula options; appended to output file names
suffix <- paste0(
  if (spline_individuales) "si" else "",
  if (include_sex) "s" else "",
  if (include_stage_sex) "ss" else "",
  if (include_nichd) "n" else "",
  if (nichd_spline) "ns" else "",
  bs_type
)

include_sex <- include_sex || include_stage_sex

Z90 <- qnorm(0.95)  # one-sided 95th percentile; used for 90% CIs (two-tailed)

# Percentile of the null distribution used as the signal-detection threshold
percentil <- "p95"

# Continuity correction for the stratified (classic) 2x2 estimators.
# When TRUE, a Haldane-Anscombe +0.5 is added to every cell before computing
classic_continuity_correction <- TRUE
continuity_correction_value <- 0.5

# Parametric-bootstrap draws for the GAM-AC CI
ac_bootstrap_n <- 1000

source("01_theme.R", local = TRUE)
theme_set(theme_base())

################################################################################
# Vocabulary-based drug mapping
################################################################################

# Builds the ATC canonical translation table from the shared OMOP vocabulary (CONCEPT.csv).
# Drugs sharing the same active compound (base name) collapse to a single canonical concept_id 
normalize_vocabulary_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)
  x[x %in% c("", "na")] <- NA_character_
  x
}

build_atc_vocabulary_map <- function(concept_path = ruta_concept) {
  if (!file.exists(concept_path)) {
    stop(sprintf("CONCEPT.csv was not found at %s", concept_path))
  }

  concept_dt <- fread(
    concept_path,
    quote = "",
    select = c(
      "concept_id", "concept_name", "vocabulary_id",
      "concept_class_id", "standard_concept", "concept_code", "invalid_reason"
    )
  )

  concept_dt[, `:=`(
    vocabulary_id_key = tolower(trimws(vocabulary_id)),
    concept_class_id_key = tolower(trimws(concept_class_id)),
    invalid_reason_key = fifelse(is.na(invalid_reason), "", trimws(invalid_reason)),
    atc_name_key = normalize_vocabulary_key(concept_name)
  )]

  atc_dt <- concept_dt[
    vocabulary_id_key == "atc" &
      invalid_reason_key == "" &
      grepl("atc", concept_class_id_key),
    .(
      atc_concept_code = as.character(concept_code),
      atc_concept_id = as.character(concept_id),
      atc_concept_name = concept_name,
      atc_name_key,
      standard_concept
    )
  ]

  atc_dt[, atc_signature := vapply(
    strsplit(gsub("\\s*-\\s*", " and ", atc_name_key), "\\s+and\\s+"),
    function(parts) paste(sort(trimws(parts)), collapse = " | "),
    character(1)
  )]

  # Deduplicate: prefers the standard concept; break ties by smallest concept_id.
  atc_dt[, `:=`(
    standard_priority = fifelse(standard_concept == "S", 1L, 0L),
    atc_concept_id_num = as.numeric(atc_concept_id)
  )]
  setorder(atc_dt, atc_concept_code, -standard_priority, atc_concept_id_num)
  atc_dt <- atc_dt[, .SD[1L], by = atc_concept_code]
  atc_dt[, c("standard_priority", "atc_concept_id_num") := NULL]

  return(atc_dt[])
}

build_drug_translation_table <- function() {
  atc_dt <- build_atc_vocabulary_map()
  atc_dt[, base_name := normalize_vocabulary_key(sub("[;,].*", "", atc_concept_name))]

  # One canonical concept_id per active compound (base name); smallest id is the tiebreaker.
  canonical_map <- atc_dt[
    ,
    .(canonical_id = atc_concept_id[which.min(as.numeric(atc_concept_id))]),
    by = base_name
  ]

  merge(
    atc_dt[, .(atc_concept_id, atc_concept_code, atc_concept_name, base_name, atc_signature)],
    canonical_map,
    by = "base_name",
    all.x = TRUE
  )
}

# MedDRA name lookup (meddra_concept_id -> concept_name) from CONCEPT.csv.
# Used to attach readable event names to meddra_concept_id in descriptive outputs.
build_meddra_name_map <- function(concept_path = ruta_concept) {
  if (!file.exists(concept_path)) {
    stop(sprintf("CONCEPT.csv was not found at %s", concept_path))
  }

  concept_dt <- fread(
    concept_path,
    quote = "",
    select = c("concept_id", "concept_name", "vocabulary_id", "invalid_reason")
  )

  meddra_dt <- concept_dt[
    tolower(trimws(vocabulary_id)) == "meddra" &
      (is.na(invalid_reason) | trimws(invalid_reason) == ""),
    .(meddra_concept_id = as.character(concept_id), meddra_name = concept_name)
  ]
  unique(meddra_dt, by = "meddra_concept_id")
}

################################################################################
# Function to build triplets
################################################################################

# Generates all (drugA, drugB, event) triplets for a single report.
# Returns NULL when the report has fewer than 2 drugs or no events.
# drugA <= drugB ordering is enforced to prevent duplicate pairs from reversed order.

make_triplets <- function(drug, event, report_id, nichd_stage) {
  
  if (length(drug) < 2 || length(event) < 1) return(NULL)
  
  drug <- unique(drug)
  event <- unique(event)
  
  if (length(drug) == 2) {
    combination <- matrix(c(min(drug), max(drug)), nrow = 1)
  } else {
    combination <- t(combn(drug, 2))
    combination <- t(apply(combination, 1, function(x) c(min(x), max(x))))
  }
  
  n_combination <- nrow(combination)
  n_events <- length(event)
  
  data.table(
    safetyreportid = report_id,
    drugA = rep(combination[,1], times = n_events),
    drugB = rep(combination[,2], times = n_events),
    meddra = rep(event, each = n_combination),
    nichd_num = nichd_stage
  )
}

################################################################################
# Effect size function (fold-change)
################################################################################

# Samples fold-changes from a shifted exponential: FC ~ 1 + Exp(lambda).
# Default lambda = 0.75 gives a typical range of [1, 10] with right skew.

fold_change <- function(n, lambda = 0.75) {
  1 + rexp(n, rate = lambda)
}

################################################################################
# Batch co-administration counts by NICHD stage
################################################################################

# Computes A+B co-administration counts by NICHD stage for a batch of triplets.
# Uses a single set of joins keyed on unique drug pairs

compute_coadmin_batch <- function(pairs_dt, ade_dt) {
  
  pair_meta <- unique(pairs_dt[, .(triplet_id, drugA, drugB, meddra)])
  
  if (nrow(pair_meta) == 0) {
    return(data.table(
      triplet_id = integer(),
      drugA = character(),
      drugB = character(),
      meddra = character(),
      nichd_num = integer(),
      nichd = character(),
      n_coadmin_stage = integer()
    ))
  }
  
  # Joins each unique drug once to avoid row explosion when many triplets share drugA or drugB.
  unique_pairs <- unique(pair_meta[, .(drugA, drugB)])
  unique_a <- unique(pair_meta[, .(drugA)])
  unique_b <- unique(pair_meta[, .(drugB)])
  
  reports_a <- ade_dt[
    unique_a,
    on = .(atc_concept_id = drugA),
    nomatch = 0L,
    .(drugA = i.drugA, safetyreportid, nichd_num)
  ]
  
  reports_b <- ade_dt[
    unique_b,
    on = .(atc_concept_id = drugB),
    nomatch = 0L,
    .(drugB = i.drugB, safetyreportid)
  ]
  
  # Inner-join A and B reports on safetyreportid to get report-level co-administrations.
  coadmin_drug <- reports_a[
    reports_b,
    on = .(safetyreportid),
    nomatch = 0L,
    allow.cartesian = TRUE,
    .(drugA, drugB = i.drugB, safetyreportid, nichd_num)
  ]
  
  # Restrict to requested pairs (filters out off-target combinations from the join).
  coadmin_drug <- coadmin_drug[
    unique_pairs,
    on = .(drugA, drugB),
    nomatch = 0L
  ]
  
  stage_counts <- coadmin_drug[
    , .(n_coadmin_stage = .N),
    by = .(drugA, drugB, nichd_num)
  ]
  
  full_grid <- CJ(
    row_id = seq_len(nrow(unique_pairs)),
    nichd_num = seq_along(niveles_nichd),
    unique = TRUE
  )
  full_grid[, `:=`(
    drugA = unique_pairs$drugA[row_id],
    drugB = unique_pairs$drugB[row_id]
  )]
  full_grid[, row_id := NULL]
  
  result_drug <- merge(
    full_grid,
    stage_counts,
    by = c("drugA", "drugB", "nichd_num"),
    all.x = TRUE
  )
  
  result_drug[is.na(n_coadmin_stage), n_coadmin_stage := 0L]
  result_drug[, nichd := niveles_nichd[nichd_num]]
  
  result <- merge(
    pair_meta,
    result_drug,
    by = c("drugA", "drugB"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )
  
  setcolorder(
    result,
    c("triplet_id", "drugA", "drugB", "meddra", "nichd_num", "nichd", "n_coadmin_stage")
  )
  
  return(result[order(triplet_id, nichd_num)])
}

################################################################################
# Dynamic pattern generation function
################################################################################

# Returns a stage pattern vector in [-1, 1] for signal injection.
# Shapes: uniform (constant 0), increase (tanh ramp up), decrease (tanh ramp down),
# plateau (bell, peak at central stages), inverse_plateau (U-shape, trough at central stages).
# N is always 7 (one value per NICHD stage).

generate_dynamic <- function(type, N = 7) {
  type <- as.character(type)
  if (type == "uniform") {
    return(rep(0, N))
  }
  if (type == "increase") {
    return(tanh(seq(-pi, pi, length.out = N)))
  }
  if (type == "decrease") {
    return(-tanh(seq(-pi, pi, length.out = N)))
  }
  if (type == "plateau") {
    return(c(
      tanh(seq(-pi, pi, length.out = floor(N/2))),
      tanh(seq(pi, -pi, length.out = ceiling(N/2)))
    ))
  }
  if (type == "inverse_plateau") {
    return(c(
      tanh(seq(pi, -pi, length.out = floor(N/2))),
      tanh(seq(-pi, pi, length.out = ceiling(N/2)))
    ))
  }
}

################################################################################
# Signal injection function
################################################################################

# Injects a synthetic DDI signal into co-administration reports and returns the
# report ids that should carry the event (no dataset copy is made).
#
# Injection model:
#   e_j  = additive base rate: P(event|A alone) + P(event|B alone) - product
#          (proxy; does not assume true independence)
#   t_ij = fold_change * e_j
#   p_dynamic(j) = t_ij + generate_dynamic(type, N=7) * t_ij   (clipped to [0.001, 0.999])
#   Y_new ~ Bernoulli(p_dynamic(stage_j)) per co-administration report
# Pre-existing events are preserved (e_final = max(e_old, e_new)).
#
# Returns: list with success, injection_success, n_injected, n_coadmin,
#   reports_to_mark, message, diagnostics.

inject_signal <- function(drugA_id, drugB_id, event_id, 
                          dynamic_type, fold_change, 
                          ade_raw_dt) {
  
  if (drugA_id > drugB_id) {   # canonical ordering: callers may not guarantee it
    temp <- drugA_id
    drugA_id <- drugB_id
    drugB_id <- temp
  }
  
  # ade_raw_dt is read-only: injection returns a report-id vector, so the shared
  # dataset is never mutated and triplets cannot contaminate each other.

  # 1. Co-administration reports
  reports_A <- unique(ade_raw_dt[atc_concept_id == drugA_id, safetyreportid])
  reports_B <- unique(ade_raw_dt[atc_concept_id == drugB_id, safetyreportid])
  reports_AB <- intersect(reports_A, reports_B)

  if (length(reports_AB) <= 0) {
    return(list(
      success = FALSE,
      injection_success = FALSE,
      n_injected = 0,
      n_coadmin = length(reports_AB),
      reports_to_mark = integer(0),
      message = sprintf(
        "insufficient co-administration: %d reports",
        length(reports_AB)
      ),
      diagnostics = list(
        reason = "insufficient_coadmin",
        n_coadmin = length(reports_AB),
        drugA = drugA_id,
        drugB = drugB_id,
        event = event_id
      )
    ))
  }
  
  # 2. Build the per-report injection table (co-administration reports only)
  target_reports <- unique(ade_raw_dt[
    safetyreportid %in% reports_AB,
    .(safetyreportid, nichd, nichd_num)
  ])

  # Flag pre-existing occurrences of the event so they are preserved, not double-counted.
  event_in_report <- unique(ade_raw_dt[
    meddra_concept_id == event_id, 
    safetyreportid
  ])
  
  target_reports[, e_old := as.integer(safetyreportid %in% event_in_report)]
  
  # 3. Base rate e_j via the additive independence proxy:
  #    e_j = P(event|A) + P(event|B) - P(event|A)*P(event|B)
  # Reports seen only in co-administration have no solo evidence; base rate defaults to 0
  # to avoid NaN from mean() over an empty set propagating to t_ij.
  reports_A_clean <- setdiff(reports_A, reports_AB)
  reports_B_clean <- setdiff(reports_B, reports_AB)
  p_baseA <- if (length(reports_A_clean) > 0) mean(reports_A_clean %in% event_in_report) else 0
  p_baseB <- if (length(reports_B_clean) > 0) mean(reports_B_clean %in% event_in_report) else 0
  p_base0 <- length(event_in_report) / length(unique(ade_raw_dt$safetyreportid))  # global rate (unused; kept for diagnostics)

  e_j <- p_baseA + p_baseB - (p_baseA * p_baseB)

  # t_ij: scaled effect size applied to the stage probabilities
  t_ij <- fold_change * e_j

  # 4. Stage-specific reporting probabilities:
  #    p_dynamic(j) = t_ij + generate_dynamic(type) * t_ij, clipped to [0.001, 0.999]
  N <- 7
  bprobs <- rep(t_ij, N)
  dy <- generate_dynamic(dynamic_type, N) * t_ij
  rprobs <- pmax(pmin(bprobs + dy, 0.999), 0.001)

  # 5. Stage probability table
  stage_probs <- data.table(
    nichd_num = 1:N,
    bprobs = bprobs,
    dy = dy,
    p_dynamic = rprobs
  )
  
  # 6. Draw Bernoulli outcomes and combine with pre-existing events
  target_reports <- merge(
    target_reports,
    stage_probs[, .(nichd_num, p_dynamic)],
    by = "nichd_num",
    all.x = TRUE
  )
  target_reports[, e_new := rbinom(.N, 1, p_dynamic)]
  target_reports[, e_final := pmax(e_old, e_new)]

  # 7. Reports to mark: newly injected only (pre-existing events already present)
  reports_to_mark <- target_reports[e_old == 0 & e_final == 1, safetyreportid]

  # At least 1 new injection required for success
  if (length(reports_to_mark) == 0) {
    return(list(
      success = FALSE,
      injection_success = FALSE,
      n_injected = 0,
      n_coadmin = length(reports_AB),
      reports_to_mark = integer(0),
      message = sprintf(
        "injection failed: 0 events generated (mean prob = %.4f, max = %.4f)",
        mean(target_reports$p_dynamic),
        max(target_reports$p_dynamic)
      ),
      diagnostics = list(
        reason = "zero_events_injected",
        low_probability_injection = TRUE,
        e_j = e_j,
        t_ij = t_ij,
        fold_change = fold_change,
        dynamic_type = dynamic_type,
        mean_p_dynamic = mean(target_reports$p_dynamic),
        max_p_dynamic = max(target_reports$p_dynamic),
        min_p_dynamic = min(target_reports$p_dynamic),
        n_eligible = nrow(target_reports[e_old == 0]),
        n_already_with_event = sum(target_reports$e_old),
        stage_probs = stage_probs
      )
    ))
  }
  
  # 8a. Requires at least 1 injection in the high-reporting stages of the intended dynamic.
  high_stages_by_dynamic <- list(
    "uniform" = 1:7,      # all relevant stages
    "increase" = c(6L, 7L),
    "decrease" = c(1L, 2L),
    "plateau" = c(3L, 4L, 5L),
    "inverse_plateau"  = c(1L, 7L)
  )
  
  high_stages <- high_stages_by_dynamic[[dynamic_type]]
  
  injection_by_stage_temp <- target_reports[
    safetyreportid %in% reports_to_mark, .N, by = nichd_num]
  
  n_injected_high <- injection_by_stage_temp[nichd_num %in% high_stages, sum(N, na.rm = TRUE)]
  
  # sum() on empty table returns NA, not 0
  if (length(n_injected_high) == 0 || is.na(n_injected_high)) n_injected_high <- 0L
  
  if (n_injected_high == 0L) {
    return(list(
      success = FALSE,
      injection_success = FALSE,
      n_injected = length(reports_to_mark),
      n_coadmin = length(reports_AB),
      reports_to_mark = integer(0),
      message = sprintf(
        "injection without signal in key stages: 0 events in stages %s (dynamic: %s, total injected: %d)",
        paste(high_stages, collapse = ","),
        dynamic_type,
        length(reports_to_mark)
      ),
      diagnostics = list(
      reason = "zero_events_in_high_stages",
      n_injected_high = n_injected_high,
      high_stages = high_stages,
      n_injected_total = length(reports_to_mark),
      e_j = e_j,
      t_ij = t_ij,
      fold_change = fold_change,
      dynamic_type = dynamic_type,
      stage_probs = stage_probs
    )
    ))
  }

  # 8b. reports_to_mark is returned to the caller and applied on-the-fly inbuild_eval_table
  injection_rate <- length(reports_to_mark) / nrow(target_reports[e_old == 0])
  
  # 9. Diagnostics bundle returned with success result
  diagnostics <- list(
    e_j = e_j,
    t_ij = t_ij,
    fold_change = fold_change,
    dynamic_type = dynamic_type,
    stage_probs = stage_probs,
    mean_p_dynamic = mean(target_reports$p_dynamic),
    max_p_dynamic = max(target_reports$p_dynamic),
    min_p_dynamic = min(target_reports$p_dynamic),
    n_eligible = nrow(target_reports),
    n_already_with_event = sum(target_reports$e_old),
    n_without_event = nrow(target_reports[e_old == 0]),
    n_new_events = length(reports_to_mark),
    injection_rate = length(reports_to_mark) / nrow(target_reports[e_old == 0]),
    injection_by_stage = injection_by_stage_temp
  )
  
  diagnostics$n_injected_high <- n_injected_high
  diagnostics$high_stages <- high_stages
  diagnostics$n_injected_total <- length(reports_to_mark)

  return(list(
    success = TRUE,
    injection_success = TRUE,
    n_injected = length(reports_to_mark),
    n_coadmin = length(reports_AB),
    reports_to_mark = reports_to_mark,
    message = sprintf(
      "injection successful: %d events injected across %d co-administration reports (rate: %.2f%%)",
      length(reports_to_mark),
      length(reports_AB),
      injection_rate * 100
    ),
    diagnostics = diagnostics
  ))
}

################################################################################
# Helper function: compute basic counts
################################################################################

# Returns n_events, n_events_coadmin, and n_coadmin for a (drugA, drugB, event) triplet.
# Injected events arrive as an explicit report-id vector

calc_basic_counts <- function(ade_data, drugA, drugB, meddra, simulated_reports = integer(0)) {
  r_a <- unique(ade_data[atc_concept_id == drugA, safetyreportid])
  r_b <- unique(ade_data[atc_concept_id == drugB, safetyreportid])
  r_coadmin <- intersect(r_a, r_b)
  r_ea <- unique(ade_data[meddra_concept_id == meddra, safetyreportid])
  # Merges injected event reports into the event set.
  if (length(simulated_reports) > 0) {
    r_ea <- union(r_ea, simulated_reports)
  } else if ("simulated_event" %in% names(ade_data)) {
    r_ea_sim <- unique(ade_data[simulated_event == TRUE & simulated_meddra == meddra, safetyreportid])
    r_ea <- union(r_ea, r_ea_sim)
  }
  list(
    n_events = length(r_ea),
    n_events_coadmin = length(intersect(r_coadmin, r_ea)),
    n_coadmin = length(r_coadmin)
  )
}

################################################################################
# Report-level evaluation table (shared by all per-triplet estimators)
################################################################################

# Builds the report-level modeling table once per (triplet, dataset).
# Shared by fit_gam, calculate_classic_ior, and calculate_classic_ac 
# avoids rebuilding the same membership sets three times. 
#
# Returns: one row per report with columns
#   safetyreportid, nichd, nichd_num, [sex], droga_a, droga_b, droga_ab, ea_ocurrio

build_eval_table <- function(ade_data, drugA_id, drugB_id, event_id,
                             simulated_reports = integer(0), include_sex = FALSE) {
  cols <- c("safetyreportid", "nichd", "nichd_num")
  if (include_sex) cols <- c(cols, "sex")

  # One row per report (the unit of analysis)
  dt <- unique(ade_data[, ..cols])

  # Membership sets for each exposure and the event
  reports_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reports_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])
  reports_ea <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])
  if (length(simulated_reports) > 0) {
    reports_ea <- union(reports_ea, simulated_reports)
  }

  # Binary exposure/outcome indicators
  dt[, droga_a := as.integer(safetyreportid %in% reports_a)]
  dt[, droga_b := as.integer(safetyreportid %in% reports_b)]
  dt[, droga_ab := as.integer(droga_a == 1L & droga_b == 1L)]
  dt[, ea_ocurrio := as.integer(safetyreportid %in% reports_ea)]

  if (include_sex) {
    dt[, sex := toupper(trimws(sex))]
    dt[sex == "M", sex := "MALE"]
    dt[sex == "F", sex := "FEMALE"]
    dt[, sex := factor(sex, levels = c("MALE", "FEMALE"))]
  }

  return(dt[])
}

################################################################################
# GAM fitting function
################################################################################

# Fits a logistic GAM for the drug-drug interaction signal, 
# Computes per-stage log-IOR (delta-method SE from the covariance matrix)
# Computes GAM-AC (parametric bootstrap). 
#
# Key parameters:
#   spline_individuales: splines for individual drug baseline risks (vs. linear)
#   nichd_spline: NICHD main effect as spline vs. linear term
#   bs_type: basis type ("cs" default); "tp" or "cr" also accepted
#   select: penalise-to-zero shrinkage for optional terms
#   k_spline: knot count — should match the number of NICHD stages (7)
#   method: GAM estimation method; keep as "fREML"

fit_gam <- function(drugA_id, drugB_id, event_id, ade_data,
                                 nichd_spline = TRUE,
                                 include_nichd = TRUE,
                                 spline_individuales = FALSE,
                                 bs_type = "cs",
                                 select = FALSE,
                                 include_sex = FALSE,
                                 include_stage_sex = FALSE,
                                 k_spline = 7,
                                 method = "fREML",
                                 simulated_reports = integer(0),
                                 eval_dt = NULL) {
  ###########
  # 1. Builds the report-level modeling dataset
  ###########

  # Reuse a precomputed evaluation table when provided; otherwise build it here.
  if (is.null(eval_dt)) {
    eval_dt <- build_eval_table(ade_data, drugA_id, drugB_id, event_id,
                                simulated_reports, include_sex)
  }
  datos_modelo <- eval_dt

  # Counts derived directly from the report-level table
  n_coadmin <- sum(datos_modelo$droga_ab)
  n_events_total <- sum(datos_modelo$ea_ocurrio)
  n_events_coadmin <- sum(datos_modelo$droga_ab == 1L & datos_modelo$ea_ocurrio == 1L)

  ###########
  # 2. Build formula from parameters
  ###########

  formula_parts <- "ea_ocurrio ~ "

  if (!spline_individuales) {
    formula_parts <- paste0(formula_parts, "droga_a + droga_b + ")
  } else {
    # Spline-smoothed individual drug baseline risks
    formula_parts <- paste0(
      formula_parts,
      sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_a) + ", 
              k_spline, bs_type),
      sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_b) + ", 
              k_spline, bs_type)
    )
  }

  if (include_nichd) {
    if (nichd_spline) {
      formula_parts <- paste0(
        formula_parts,
        sprintf("s(nichd_num, k = %d, bs = '%s') + ", k_spline, bs_type)
      )
    } else {
      formula_parts <- paste0(formula_parts, "nichd_num + ")
    }
  }

  # Interaction spline — key term; do not remove or modify
  formula_parts <- paste0(
    formula_parts,
    sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_ab)", k_spline, bs_type)
  )
  
  if (include_sex) {
    if (include_stage_sex) {
      # Stage-by-sex interaction spline
      formula_parts <- paste0(
        formula_parts,
        sprintf(" + s(nichd_num, k = %d, bs = '%s', by = sex)", 
                k_spline, bs_type)
      )
    } else {
      formula_parts <- paste0(formula_parts, " + sex")
    }
  }

  formula_final <- as.formula(formula_parts)

  ###########
  # 3. Model fitting
  ###########
  
  tryCatch({
    
    modelo <- bam(
      formula = formula_final,
      data = datos_modelo,
      family = binomial(link = "logit"),
      method = method,
      select = select,    
      discrete = TRUE,
      nthreads = 1  # avoid conflicts with the outer doParallel cluster in 10_augmentation
    )
    
    ###########
    # 3a. Computes log-IOR per NICHD stage
    ###########

    # All 4 exposure combinations x 7 stages for contrast computation
    grid_dif <- CJ(
      nichd_num = 1:7, 
      droga_a = c(0, 1), 
      droga_b = c(0, 1)
    )
    grid_dif[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]
    
    # Male as reference for sex (if in formula)
    if (include_sex) {
      grid_dif[, sex := factor("MALE", levels = c("MALE", "FEMALE"))]
    }

    # Point predictions on the link scale; SE comes from the covariance matrix below.
    grid_dif[, lp := predict(modelo, newdata = grid_dif, type = "link")]

    # Pivot to wide so contrast arithmetic is element-wise per stage.
    # With single value.var, dcast names columns "droga_a_droga_b" (e.g. "1_1").
    w_lp <- dcast(grid_dif, nichd_num ~ droga_a + droga_b, value.var = "lp")
    
    # log(IOR) = log(OR_11) - log(OR_10) - log(OR_01) + log(OR_00)
    log_ior <- w_lp[["1_1"]] - w_lp[["1_0"]] - w_lp[["0_1"]] + w_lp[["0_0"]]

    ###########
    # 3b. Log-IOR SE via the linear predictor covariance matrix
    ###########

    Xp <- predict(modelo, newdata = grid_dif, type = "lpmatrix")
    Vb <- vcov(modelo, unconditional = TRUE)  # unconditional = TRUE includes smoothing-parameter uncertainty
    cov_link <- Xp %*% Vb %*% t(Xp)
    
    log_ior_se <- numeric(7)
    for (stage in 1:7) {
      idx_00 <- which(grid_dif$nichd_num == stage & 
                        grid_dif$droga_a == 0 & grid_dif$droga_b == 0)
      idx_01 <- which(grid_dif$nichd_num == stage & 
                        grid_dif$droga_a == 0 & grid_dif$droga_b == 1)
      idx_10 <- which(grid_dif$nichd_num == stage & 
                        grid_dif$droga_a == 1 & grid_dif$droga_b == 0)
      idx_11 <- which(grid_dif$nichd_num == stage & 
                        grid_dif$droga_a == 1 & grid_dif$droga_b == 1)
      
      # Contrast vector
      cvec <- rep(0, nrow(grid_dif))
      cvec[c(idx_11, idx_10, idx_01, idx_00)] <- c(1, -1, -1, 1)
      
      # SE = sqrt(c' Sigma c), where c = contrast vector for IOR
      log_ior_se[stage] <- sqrt(max(
        as.numeric(t(cvec) %*% cov_link %*% cvec), 
        0
      ))
    }
    
    ###########
    # 3c. Confidence intervals and summary metrics
    ###########

    if (!exists("Z90")) {
      Z90 <- qnorm(0.95)
    }
    
    log_ior_lower90 <- log_ior - Z90 * log_ior_se
    log_ior_upper90 <- log_ior + Z90 * log_ior_se
    ior_values <- exp(log_ior)
    
    n_stages_significant <- sum(log_ior_lower90 > 0)
    max_ior <- max(ior_values)
    mean_ior <- mean(ior_values)

    ###########
    # 3d. Additive interaction contrast (Thakrar) per stage with 90% CI
    ###########
    
    # Evaluate on all 7 NICHD stages: spline extrapolation fills stages absent from data.
    stages <- 1:7
    
    # Prediction grid: 4 exposure combinations per stage
    nd_ac <- rbindlist(lapply(stages, function(s) {
      data.table(
        nichd_num = s,
        droga_a   = c(0, 1, 0, 1),
        droga_b   = c(0, 0, 1, 1),
        droga_ab  = c(0, 0, 0, 1)
      )
    }), use.names = TRUE)
    
    # Covariates for sex and NICHD main effect if included in formula
    if (include_sex) {
      nd_ac[, sex := factor(levels(datos_modelo$sex)[1], 
                             levels = levels(datos_modelo$sex))]
    }
    if (include_nichd && !nichd_spline) {
      nd_ac[, nichd := factor(niveles_nichd[nichd_num],
                                levels = niveles_nichd,
                                ordered = TRUE)]
    }
    
    # Parametric bootstrap for the contrast CI and its standard error.
    # Link predictions formerly computed here were never used downstream.
    X_ac <- predict(modelo, newdata = nd_ac, type = "lpmatrix")
    beta_hat <- coef(modelo)
    V_beta <- vcov(modelo, unconditional = TRUE)

    B <- ac_bootstrap_n

    # Simulate coefficients from their joint distribution
    # beta_sim ~ MVN(beta_hat, V_beta)
    beta_sims <- mvrnorm(n = B, mu = beta_hat, Sigma = V_beta)

    # Predicted reporting proportion per exposure cell for each simulated vector
    p_sims <- plogis(X_ac %*% t(beta_sims))

    # Additive interaction contrast (Thakrar) on predicted reporting proportions:
    # p11 - p10 - p01 + p00
    # detection uses the raw contrast and its 90% CI lower bound
    calc_add <- function(p) {
      # p: 4 predicted reporting proportions per stage: [p00, p10, p01, p11]
      p11 <- p[4]; p10 <- p[2]; p01 <- p[3]; p00 <- p[1]
      p11 - p10 - p01 + p00
    }

    # Per-stage summary from the bootstrap 
    ac_dt <- nd_ac[, {
      idx <- .I
      p_mat <- p_sims[idx, , drop = FALSE]
      add_sim <- apply(p_mat, 2, calc_add)

      data.table(
        AC = mean(add_sim),
        AC_lower90 = quantile(add_sim, 0.05),
        AC_upper90 = quantile(add_sim, 0.95)
      )
    }, by = nichd_num]

    # Safety join: ensures result is always length 7 and stage-ordered.
    ac_dt <- merge(data.table(nichd_num = 1:7), ac_dt, by = "nichd_num", all.x = TRUE)
    setorder(ac_dt, nichd_num)

    ac_values <- ac_dt$AC
    ac_lower90 <- ac_dt$AC_lower90
    ac_upper90 <- ac_dt$AC_upper90

    ###########
    # 4- Results
    ###########

    return(list(
      success = TRUE,
      n_events = n_events_total,
      n_coadmin = n_coadmin,
      n_events_coadmin = n_events_coadmin,
      log_ior = log_ior,
      log_ior_lower90 = log_ior_lower90,
      log_ior_upper90 = log_ior_upper90,
      log_ior_se = log_ior_se,
      ior_values = ior_values,
      n_stages_significant = n_stages_significant,
      max_ior = max_ior,
      mean_ior = mean_ior,
      ac_values = ac_values,
      ac_lower90 = ac_lower90,
      ac_upper90 = ac_upper90,
      n_stages_ac_significant = sum(ac_lower90 > 0, na.rm = TRUE),
      model_aic = AIC(modelo),
      model_deviance = deviance(modelo),
      formula_used = formula_parts,
      nichd_spline = nichd_spline,
      include_nichd = include_nichd,
      spline_individuales = spline_individuales,
      bs_type = bs_type,
      select = select,
      include_sex = include_sex,
      include_stage_sex = include_stage_sex,
      k_spline = k_spline
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE, 
      n_events = n_events_total, 
      n_coadmin = n_coadmin,
      error_msg = e$message,
      formula_attempted = formula_parts
    ))
  })
}

################################################################################
# Classical IOR calculation
################################################################################

# Computes stratified IOR using per-stage 2x2 contingency tables (Woolf method).
# IOR = (OR_11 * OR_00) / (OR_10 * OR_01); OR_00 = 1 by definition.
# Haldane-Anscombe continuity correction applied when classic_continuity_correction is TRUE.
# Returns a list with success and results_by_stage.

calculate_classic_ior <- function(drugA_id, drugB_id, event_id, ade_data,
                                   simulated_reports = integer(0), eval_dt = NULL) {

  # Reuses precomputed evaluation table when provided; otherwise build from ADE data.
  if (is.null(eval_dt)) {
    eval_dt <- build_eval_table(ade_data, drugA_id, drugB_id, event_id, simulated_reports)
  }
  datos_unicos <- eval_dt

  # Haldane-Anscombe +0.5 per cell when correction is enabled; cc=0 gives raw textbook form.
  cc <- if (classic_continuity_correction) continuity_correction_value else 0

  stage_results <- datos_unicos[, {

    # Exposure-group cell counts (raw, kept for diagnostics)
    # Group 11: A + B co-administration
    n_11_evento <- sum(droga_a == 1 & droga_b == 1 & ea_ocurrio == 1)
    n_11_no_evento <- sum(droga_a == 1 & droga_b == 1 & ea_ocurrio == 0)
    
    # Group 10: A only
    n_10_evento <- sum(droga_a == 1 & droga_b == 0 & ea_ocurrio == 1)
    n_10_no_evento <- sum(droga_a == 1 & droga_b == 0 & ea_ocurrio == 0)

    # Group 01: B only
    n_01_evento <- sum(droga_a == 0 & droga_b == 1 & ea_ocurrio == 1)
    n_01_no_evento <- sum(droga_a == 0 & droga_b == 1 & ea_ocurrio == 0)

    # Group 00: neither A nor B (reference)
    n_00_evento <- sum(droga_a == 0 & droga_b == 0 & ea_ocurrio == 1)
    n_00_no_evento <- sum(droga_a == 0 & droga_b == 0 & ea_ocurrio == 0)

    # Corrected counts for estimation; raw counts kept for diagnostic columns.
    a11 <- n_11_evento + cc; b11 <- n_11_no_evento + cc
    a10 <- n_10_evento + cc; b10 <- n_10_no_evento + cc
    a01 <- n_01_evento + cc; b01 <- n_01_no_evento + cc
    a00 <- n_00_evento + cc; b00 <- n_00_no_evento + cc

    or_11 <- (a11 / b11) / (a00 / b00)
    or_10 <- (a10 / b10) / (a00 / b00)
    or_01 <- (a01 / b01) / (a00 / b00)
    or_00 <- 1  # reference; OR_00 = 1 by definition
    # IOR = OR_11 / (OR_10 * OR_01) since OR_00 = 1
    ior_val <- (or_11 * or_00) / (or_10 * or_01)
    log_ior <- log(ior_val)
    
    # Variance on the log scale (Woolf method) on the corrected cells
    var_log_ior <- (1/a11 + 1/b11 +
                1/a10 + 1/b10 +
                1/a01 + 1/b01 +
                1/a00 + 1/b00)
    se_log_ior <- sqrt(var_log_ior)
    
    # 90% CI on the log scale, then exponentiated
    z90 <- qnorm(0.95)
    log_ior_lower90 <- log_ior - z90 * se_log_ior
    log_ior_upper90 <- log_ior + z90 * se_log_ior
    ior_lower90 <- exp(log_ior_lower90)
    ior_upper90 <- exp(log_ior_upper90)
    
    data.table(
      stage = nichd_num[1],
      ior_classic = ior_val,
      log_ior_classic = log_ior,
      ior_classic_lower90 = ior_lower90,
      ior_classic_upper90 = ior_upper90,
      log_ior_classic_lower90 = log_ior_lower90,
      log_ior_classic_upper90 = log_ior_upper90,
      se_log_ior_classic = se_log_ior,
      # Diagnostics
      n_11_evento = n_11_evento,
      n_11_total = n_11_evento + n_11_no_evento
    )
    
  }, by = nichd_num]
  
  setorder(stage_results, nichd_num)

  return(list(
    success = TRUE,
    results_by_stage = stage_results
  ))
}

################################################################################
# Stratified AC calculation
################################################################################

# Computes the stratified additive interaction contrast (Thakrar) per stage 
# R11 - R10 - R01 + R00 
# Detection uses the raw contrast and its 90% CI lower bound.
# Haldane-Anscombe correction applied when classic_continuity_correction is TRUE.
# Returns a list with success and results_by_stage.

calculate_classic_ac <- function(drugA_id, drugB_id, event_id, ade_data,
                                   simulated_reports = integer(0), eval_dt = NULL) {

  # Reuses the shared report-level table when provided 
  # otherwise build it from the read-only ADE table and the injected report-id vector.
  if (is.null(eval_dt)) {
    eval_dt <- build_eval_table(ade_data, drugA_id, drugB_id, event_id, simulated_reports)
  }
  datos_unicos <- eval_dt

  # Haldane-Anscombe +0.5 per group when correction is enabled; cc=0 gives raw form.
  cc <- if (classic_continuity_correction) continuity_correction_value else 0

  stage_results <- datos_unicos[, {

    # Exposure-group counts (raw, kept for diagnostics)
    # Group 11: A+B co-administration
    n_11_evento <- sum(droga_a == 1 & droga_b == 1 & ea_ocurrio == 1)
    n_11_total <- sum(droga_a == 1 & droga_b == 1)

    # Group 10: A only
    n_10_evento <- sum(droga_a == 1 & droga_b == 0 & ea_ocurrio == 1)
    n_10_total <- sum(droga_a == 1 & droga_b == 0)

    # Group 01: B only
    n_01_evento <- sum(droga_a == 0 & droga_b == 1 & ea_ocurrio == 1)
    n_01_total <- sum(droga_a == 0 & droga_b == 1)

    # Group 00: neither (reference)
    n_00_evento <- sum(droga_a == 0 & droga_b == 0 & ea_ocurrio == 1)
    n_00_total <- sum(droga_a == 0 & droga_b == 0)
    
    # Corrected denominators (cc = 0 reproduces the raw proportions).
    d11 <- n_11_total + 2 * cc
    d10 <- n_10_total + 2 * cc
    d01 <- n_01_total + 2 * cc
    d00 <- n_00_total + 2 * cc

    # Without correction (cc = 0) an empty exposure group leaves the risk undefined 
    if (d11 == 0 || d10 == 0 || d01 == 0 || d00 == 0) {
      data.table(
        stage = nichd_num[1],
        AC_classic = NA_real_,
        AC_classic_lower90 = NA_real_,
        AC_classic_upper90 = NA_real_,
        AC_classic_se = NA_real_,
        n_11_evento = n_11_evento,
        n_11_total = n_11_total,
        n_10_evento = n_10_evento,
        n_10_total = n_10_total,
        n_01_evento = n_01_evento,
        n_01_total = n_01_total,
        n_00_evento = n_00_evento,
        n_00_total = n_00_total,
        insufficient_data = TRUE
      )
    } else {  # cc > 0 guarantees all denominators are positive
      R11 <- (n_11_evento + cc) / d11
      R10 <- (n_10_evento + cc) / d10
      R01 <- (n_01_evento + cc) / d01
      R00 <- (n_00_evento + cc) / d00

      # Additive interaction contrast 
      ac_val <- R11 - R10 - R01 + R00

      # Binomial variance per proportion; boundary correction (0.25/n) at r = 0 or 1.
      var_r <- function(r, n) ifelse(r > 0 & r < 1, r*(1-r)/n, 0.25/n)
      var_R11 <- var_r(R11, d11)
      var_R10 <- var_r(R10, d10)
      var_R01 <- var_r(R01, d01)
      var_R00 <- var_r(R00, d00)

      # Var(contrast) is the sum of the four binomial variances 
      se_ac <- sqrt(var_R11 + var_R10 + var_R01 + var_R00)

      # 90% CI for the additive contrast (binomial SE); detection uses the lower bound
      z90 <- qnorm(0.95)
      ac_lower90 <- ac_val - z90 * se_ac
      ac_upper90 <- ac_val + z90 * se_ac

      data.table(
        stage = nichd_num[1],
        AC_classic = ac_val,
        AC_classic_lower90 = ac_lower90,
        AC_classic_upper90 = ac_upper90,
        AC_classic_se = se_ac,
        # Individual proportions for diagnostics
        R11 = R11, R10 = R10, R01 = R01, R00 = R00,
        # Diagnostic counts
        n_11_evento = n_11_evento, n_11_total = n_11_total,
        n_10_evento = n_10_evento, n_10_total = n_10_total,
        n_01_evento = n_01_evento, n_01_total = n_01_total,
        n_00_evento = n_00_evento, n_00_total = n_00_total,
        insufficient_data = FALSE
      )
    }
  }, by = nichd_num]
  setorder(stage_results, nichd_num)

  # All-NA means no stage had enough data for a valid estimate
  if (all(is.na(stage_results$AC_classic))) {
    return(list(
      success = FALSE,
      message = "Insufficient data in all stages",
      results_by_stage = stage_results
    ))
  }
  
  return(list(
    success = TRUE,
    results_by_stage = stage_results
  ))
}

################################################################################
# Bootstrap function by dynamic type and stage
################################################################################

# Bootstraps the mean difference in log-IOR between a dynamic type and the uniform baseline
# Returns NA when fewer than 3 observations are available.

bootstrap_dynamic_diff <- function(data, dynamic_type, stage_num, n_boot = 100) {

  target_data <- data[dynamic == dynamic_type & stage == stage_num, log_ior]
  uniform_data <- data[dynamic == "uniform" & stage == stage_num, log_ior]
  
  if (length(target_data) < 3 || length(uniform_data) < 3) {
    return(data.table(
      mean_diff = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    ))
  }
  
  boot_diffs <- replicate(n_boot, {
    target_sample <- sample(target_data, replace = TRUE)
    uniform_sample <- sample(uniform_data, replace = TRUE)
    mean(target_sample) - mean(uniform_sample)
  })
  data.table(
    mean_diff = mean(boot_diffs, na.rm = TRUE),
    ci_lower = quantile(boot_diffs, 0.025, na.rm = TRUE),
    ci_upper = quantile(boot_diffs, 0.975, na.rm = TRUE)
  )
}

################################################################################
# Bootstrap function by dynamic type and stage (AC)
################################################################################

# Additive-contrast variant of bootstrap_dynamic_diff.
# Resamples the median difference 

bootstrap_dynamic_diff_ac <- function(data, dynamic_type, stage_num, n_boot = 100) {

  target_data <- data[dynamic == dynamic_type & stage == stage_num, ac]
  uniform_data <- data[dynamic == "uniform" & stage == stage_num, ac]

  if (length(target_data) < 3 || length(uniform_data) < 3) {
    return(data.table(mean_diff = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
  }

  boot_diffs <- replicate(n_boot, {
    median(sample(target_data,  replace = TRUE)) -
    median(sample(uniform_data, replace = TRUE))
  })
    
  data.table(
    mean_diff = mean(boot_diffs, na.rm = TRUE),
    ci_lower  = quantile(boot_diffs, 0.025, na.rm = TRUE),
    ci_upper  = quantile(boot_diffs, 0.975, na.rm = TRUE)
  )
}

################################################################################
# Function to process a single positive triplet with sensitivity analysis
################################################################################

# Processes a single positive triplet across all downsampling reduction levels.
# Designed for parallel batch execution: injects the signal once, then fits
# GAM + classical IOR/AC on the full and each reduced dataset.
# Returns a combined data.table with one row per reduction level.

process_single_positive <- function(idx, pos_meta, ade_raw_dt, reduction_levels,
                                    reduced_idx_list,
                                    spline_individuales, include_sex, include_stage_sex,
                                    k_spline, bs_type, select, nichd_spline, base_seed = 9427) {
  # Unique seed per triplet (same scheme as 10_augmentation)
  set.seed(base_seed + idx)
  
  rowt <- pos_meta[idx]
  rowt$type <- "positive"

  # inject_signal returns a report-id vector; the dataset is never copied.
  inj_result <- tryCatch({
    inject_signal(
      drugA_id = rowt$drugA,
      drugB_id = rowt$drugB,
      event_id = rowt$meddra,
      dynamic_type = rowt$dynamic,
      fold_change = rowt$fold_change,
      ade_raw_dt = ade_raw_dt
    )
  }, error = function(e) {
    list(
      success = FALSE,
      injection_success = FALSE,
      n_injected = 0,
      n_coadmin = 0,
      reports_to_mark = integer(0),
      message = paste("injection error:", e$message),
      diagnostics = list(reason = "exception", error = e$message)
    )
  })

  inj_success <- inj_result$success
  n_injected_val <- inj_result$n_injected
  n_coadmin_val <- inj_result$n_coadmin
  diag_data <- list(inj_result$diagnostics)
  inj_message <- if(!is.null(inj_result$message)) inj_result$message else NA_character_
  
  t_ij_val <- if(inj_success && !is.null(inj_result$diagnostics$t_ij)) {
    inj_result$diagnostics$t_ij
  } else {
    NA_real_
  }

  rowt$t_ij <- t_ij_val

  if (!inj_success) {
    base_result <- data.table(
      triplet_id = idx,
      drugA = rowt$drugA,
      drugB = rowt$drugB,
      meddra = rowt$meddra,
      type = "positive",
      reduction_pct = 0,
      N = rowt$N,
      dynamic = rowt$dynamic,
      fold_change = rowt$fold_change,
      t_ij = t_ij_val,
      model_success = FALSE,
      injection_success = FALSE,
      n_injected = n_injected_val,
      n_coadmin = n_coadmin_val,
      n_events = NA_integer_,
      n_stages_significant = NA_integer_,
      max_ior = NA_real_,
      mean_ior = NA_real_,
      model_aic = NA_real_,
      stage = list(1:7),
      log_ior = list(rep(NA_real_, 7)),
      log_ior_lower90 = list(rep(NA_real_, 7)),
      ior_values = list(rep(NA_real_, 7)),
      classic_success = FALSE,
      log_ior_classic = list(rep(NA_real_, 7)),
      log_ior_classic_lower90 = list(rep(NA_real_, 7)),
      ior_classic = list(rep(NA_real_, 7)),
      ac_classic_success = FALSE,
      ac_values = list(rep(NA_real_, 7)),
      ac_lower90 = list(rep(NA_real_, 7)),
      ac_upper90 = list(rep(NA_real_, 7)),
      n_stages_ac_significant = NA_integer_,
      AC_classic = list(rep(NA_real_, 7)),
      AC_classic_lower90 = list(rep(NA_real_, 7)),
      AC_classic_upper90 = list(rep(NA_real_, 7)),
      AC_classic_se = list(rep(NA_real_, 7)),
      diagnostics = diag_data,
      spline_individuales = spline_individuales,
      nichd_spline = nichd_spline,
      include_sex = include_sex,
      include_stage_sex = include_stage_sex,
      k_spline = k_spline,
      bs_type = bs_type,
      select = select,
      formula_used = NA_character_,
      error_msg = inj_message
    )
    rm(inj_result); gc(verbose = FALSE)
    return(base_result)
  }

  reports_to_mark <- inj_result$reports_to_mark

  all_results <- list()

  # Base result (0% reduction): build the evaluation table once and reuse for all estimators.
  eval_base <- build_eval_table(ade_raw_dt, rowt$drugA, rowt$drugB, rowt$meddra,
                                reports_to_mark, include_sex)
  base_result <- fit_reduced_model(eval_base, rowt, 0)
  base_result$n_injected <- n_injected_val
  base_result$injection_success <- TRUE
  base_result$diagnostics <- diag_data

  all_results[[1]] <- base_result

  # Precomputed indices ensure identical row removal across all triplets
  for (red_pct in reduction_levels) {
    ade_reduced <- ade_raw_dt[reduced_idx_list[[as.character(red_pct)]]]

    eval_reduced <- build_eval_table(ade_reduced, rowt$drugA, rowt$drugB, rowt$meddra,
                                     reports_to_mark, include_sex)
    reduced_result <- fit_reduced_model(eval_reduced, rowt, red_pct)
    reduced_result$n_injected <- n_injected_val
    reduced_result$injection_success <- TRUE
    reduced_result$diagnostics <- diag_data

    all_results[[length(all_results) + 1]] <- reduced_result

    rm(ade_reduced, eval_reduced); gc(verbose = FALSE)
  }

  rm(inj_result); gc(verbose = FALSE)
  
  combined_results <- rbindlist(all_results, fill = TRUE)
  
  return(combined_results)
}

################################################################################
# Helper function: Reduce dataset by stage
################################################################################

# Reduces an augmented dataset by randomly removing a percentage of rows per stage
# 
# Parameters:
# ade_aug: augmented data.table
# reduction_pct: percentage to remove (e.g. 10 for 10%)
# nichd_col: column containing the NICHD stage
# seed: provided for reproducibility
#
# Return:
# Reduced data.table

# Returns the row indices to keep when downsampling by stage. 
# Exposed separately so a caller can precompute the indices once and reuse them across many triplets
reduce_indices_by_stage <- function(ade_aug, reduction_pct, nichd_col = "nichd", seed = NULL) {

  # Nothing to remove at 0% reduction: keeps every row
  if (reduction_pct <= 0) return(seq_len(nrow(ade_aug)))

  # Single deterministic seed for the whole reduction (reproducible)
  if (!is.null(seed)) set.seed(seed)

  keep_frac <- 1 - reduction_pct / 100

  # Per-stage row sampling in a single grouped pass; keep at least one row per stage that has data
  # .I yields the global row indices within each group.
  idx_keep <- ade_aug[, {
    n_keep <- max(1L, ceiling(.N * keep_frac))
    .I[sample(.N, min(n_keep, .N))]
  }, by = nichd_col]$V1

  sort(idx_keep)
}

reduce_dataset_by_stage <- function(ade_aug, reduction_pct, nichd_col = "nichd", seed = NULL) {
  if (reduction_pct <= 0) return(ade_aug)
  ade_aug[reduce_indices_by_stage(ade_aug, reduction_pct, nichd_col, seed)]
}

################################################################################
# Model fitting function on a reduced dataset
################################################################################

# Wrapper that fits GAM + classical IOR/AC on a precomputed evaluation table
# Assembles the result row 
# eval_dt must come from build_eval_table.

fit_reduced_model <- function(eval_dt, rowt, reduction_pct) {

  # Counts from the precomputed evaluation table (no extra ADE scans)
  counts_reduced <- list(
    n_events = sum(eval_dt$ea_ocurrio),
    n_coadmin = sum(eval_dt$droga_ab),
    n_events_coadmin = sum(eval_dt$droga_ab == 1L & eval_dt$ea_ocurrio == 1L)
  )

  model_res <- tryCatch({
    fit_gam(
      drugA_id = rowt$drugA,
      drugB_id = rowt$drugB,
      event_id = rowt$meddra,
      ade_data = NULL,
      spline_individuales = spline_individuales,
      include_sex = include_sex,
      include_stage_sex = include_stage_sex,
      k_spline = k_spline,
      bs_type = bs_type,
      select = select,
      nichd_spline = nichd_spline,
      include_nichd = include_nichd,
      eval_dt = eval_dt
    )
  }, error = function(e) {
    list(
      success = FALSE,
      n_events_total = counts_reduced$n_events,
      n_coadmin = counts_reduced$n_coadmin,
      error_msg = paste("reduced model error:", e$message)
    )
  })

  # Classical estimators reuse the same evaluation table
  classic_res <- tryCatch({
    calculate_classic_ior(rowt$drugA, rowt$drugB, rowt$meddra,
                          ade_data = NULL, eval_dt = eval_dt)
  }, error = function(e) {
    list(success = FALSE)
  })

  classic_ac <- tryCatch({
    calculate_classic_ac(rowt$drugA, rowt$drugB, rowt$meddra,
                           ade_data = NULL, eval_dt = eval_dt)
  }, error = function(e) {
    list(success = FALSE)
  })
  
  if (!model_res$success) {
    result <- data.table(
      triplet_id = rowt$triplet_id,
      drugA = rowt$drugA,
      drugB = rowt$drugB,
      meddra = rowt$meddra,
      type = rowt$type,
      reduction_pct = reduction_pct,
      N = counts_reduced$n_events_coadmin,
      dynamic = if(!is.null(rowt$dynamic)) rowt$dynamic else NA_character_,
      fold_change = if(!is.null(rowt$fold_change)) rowt$fold_change else NA_real_,
      t_ij = if(!is.null(rowt$t_ij)) rowt$t_ij else NA_real_,
      model_success = FALSE,
      injection_success = if(!is.null(rowt$injection_success)) rowt$injection_success else NA,
      n_injected = if(!is.null(rowt$n_injected)) rowt$n_injected else NA_integer_,
      n_coadmin = counts_reduced$n_coadmin,
      n_events = counts_reduced$n_events,
      n_stages_significant = NA_integer_,
      max_ior = NA_real_,
      mean_ior = NA_real_,
      model_aic = NA_real_,
      stage = list(1:7),
      log_ior = list(rep(NA_real_, 7)),
      log_ior_lower90 = list(rep(NA_real_, 7)),
      ior_values = list(rep(NA_real_, 7)),
      classic_success = classic_res$success,
      log_ior_classic = if(classic_res$success) list(classic_res$results_by_stage$log_ior_classic) else list(rep(NA_real_, 7)),
      log_ior_classic_lower90 = if(classic_res$success) list(classic_res$results_by_stage$log_ior_classic_lower90) else list(rep(NA_real_, 7)),
      ior_classic = if(classic_res$success) list(classic_res$results_by_stage$ior_classic) else list(rep(NA_real_, 7)),
      ac_classic_success = classic_ac$success,
      ac_values = list(rep(NA_real_, 7)),
      ac_lower90 = list(rep(NA_real_, 7)),
      ac_upper90 = list(rep(NA_real_, 7)),
      n_stages_ac_significant = NA_integer_,
      AC_classic = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic) else list(rep(NA_real_, 7)),
      AC_classic_lower90 = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic_lower90) else list(rep(NA_real_, 7)),
      AC_classic_upper90 = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic_upper90) else list(rep(NA_real_, 7)),
      AC_classic_se = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic_se) else list(rep(NA_real_, 7)),
      diagnostics = list(list(error = model_res$error_msg)),
      spline_individuales = spline_individuales,
      nichd_spline = nichd_spline,
      include_sex = include_sex,
      include_stage_sex = include_stage_sex,
      k_spline = k_spline,
      bs_type = bs_type,
      select = select,
      formula_used = if(!is.null(model_res$formula_attempted)) model_res$formula_attempted else NA_character_,
      error_msg = if(!is.null(model_res$error_msg)) model_res$error_msg else NA_character_
    )
    return(result)
  }
  
  result <- data.table(
    triplet_id = rowt$triplet_id,
    drugA = rowt$drugA,
    drugB = rowt$drugB,
    meddra = rowt$meddra,
    type = rowt$type,
    reduction_pct = reduction_pct,
    N = model_res$n_events_coadmin,
    dynamic = if(!is.null(rowt$dynamic)) rowt$dynamic else NA_character_,
    fold_change = if(!is.null(rowt$fold_change)) rowt$fold_change else NA_real_,
    t_ij = if(!is.null(rowt$t_ij)) rowt$t_ij else NA_real_,
    model_success = TRUE,
    injection_success = if(!is.null(rowt$injection_success)) rowt$injection_success else NA,
    n_injected = if(!is.null(rowt$n_injected)) rowt$n_injected else NA_integer_,
    n_coadmin = model_res$n_coadmin,
    n_events = model_res$n_events,
    n_stages_significant = model_res$n_stages_significant,
    max_ior = model_res$max_ior,
    mean_ior = model_res$mean_ior,
    model_aic = model_res$model_aic,
    stage = list(1:7),
    log_ior = list(model_res$log_ior),
    log_ior_lower90 = list(model_res$log_ior_lower90),
    ior_values = list(model_res$ior_values),
    classic_success = classic_res$success,
    log_ior_classic = if(classic_res$success) list(classic_res$results_by_stage$log_ior_classic) else list(rep(NA_real_, 7)),
    log_ior_classic_lower90 = if(classic_res$success) list(classic_res$results_by_stage$log_ior_classic_lower90) else list(rep(NA_real_, 7)),
    ior_classic = if(classic_res$success) list(classic_res$results_by_stage$ior_classic) else list(rep(NA_real_, 7)),
    ac_classic_success = classic_ac$success,
    ac_values = list(model_res$ac_values),
    ac_lower90 = list(model_res$ac_lower90),
    ac_upper90 = list(model_res$ac_upper90),
    n_stages_ac_significant = model_res$n_stages_ac_significant,
    AC_classic = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic) else list(rep(NA_real_, 7)),
    AC_classic_lower90 = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic_lower90) else list(rep(NA_real_, 7)),
    AC_classic_upper90 = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic_upper90) else list(rep(NA_real_, 7)),
    AC_classic_se = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic_se) else list(rep(NA_real_, 7)),
    diagnostics = list(list()),
    spline_individuales = spline_individuales,
    nichd_spline = nichd_spline,
    include_sex = include_sex,
    include_stage_sex = include_stage_sex,
    k_spline = k_spline,
    bs_type = bs_type,
    select = select,
    formula_used = model_res$formula_used
  )
  
  return(result)
}

################################################################################
# Permutation function for the null distribution
################################################################################

# Permutes drug and/or event labels within each NICHD stage to break the drug-event association 
# Uses pool_meta in 10_augmentation; permutation within-stage to preserve stage-level marginal distributions.

permute_pool <- function(pool_meta, niveles_nichd, 
                         perm_events = TRUE, 
                         perm_drugs = FALSE, 
                         seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  pool_copy <- copy(pool_meta)
  pool_copy[, drugs_perm := vector("list", .N)]
  pool_copy[, events_perm := vector("list", .N)]
  
  for (stage in niveles_nichd) {
    idx <- which(pool_copy$nichd == stage)
    
    if (length(idx) <= 1) {
      pool_copy$drugs_perm[idx] <- pool_copy$drugs[idx]
      pool_copy$events_perm[idx] <- pool_copy$events[idx]
      next
    }
    
    if (perm_drugs) {
      perm_idx_drugs <- sample(seq_along(idx), length(idx), replace = FALSE)
      pool_copy$drugs_perm[idx] <- pool_copy$drugs[idx[perm_idx_drugs]]
    } else {
      pool_copy$drugs_perm[idx] <- pool_copy$drugs[idx]
    }
    
    if (perm_events) {
      perm_idx_events <- sample(seq_along(idx), length(idx), replace = FALSE)
      pool_copy$events_perm[idx] <- pool_copy$events[idx[perm_idx_events]]
    } else {
      pool_copy$events_perm[idx] <- pool_copy$events[idx]
    }
  }
  
  pool_copy[, .(safetyreportid, nichd, nichd_num, drugs_perm, events_perm)]
}

# Replaces the permuted pool reports in ade_original with their permuted versions
# Preserves the per-report sex covariate (which the permutation does not carry)

reintroduce_permuted_reports <- function(ade_original, permuted_pool) {
  
  pool_report_ids <- unique(permuted_pool$safetyreportid)
  ade_without_pool <- ade_original[!safetyreportid %in% pool_report_ids]
  
  permuted_rows <- permuted_pool[, {
    drugs_vec <- drugs_perm[[1]]
    events_vec <- events_perm[[1]]
    
    if (length(drugs_vec) > 0 && length(events_vec) > 0) {
      # safetyreportid is supplied by `by`
      CJ(atc_concept_id = drugs_vec,
         meddra_concept_id = events_vec,
         nichd = nichd,
         nichd_num = nichd_num)
    } else {
      data.table()
    }
  }, by = safetyreportid]

  # Restores the per-report sex covariate, which the permutation metadata drops.
  # Without this the permuted rows carry NA sex and a sex-aware GAM formula loses them.
  if ("sex" %in% names(ade_original) && nrow(permuted_rows) > 0) {
    sex_by_report <- unique(ade_original[, .(safetyreportid, sex)])
    permuted_rows <- merge(permuted_rows, sex_by_report, by = "safetyreportid", all.x = TRUE)
  }

  rbindlist(list(ade_without_pool, permuted_rows), use.names = TRUE, fill = TRUE)
}

################################################################################
# Expansion function
################################################################################

# Unpacks list-type metric columns (log_ior, ac, etc.) into one row per triplet-stage
# merges the null-distribution thresholds, and adds the binary label.

expand_clean_all_metrics <- function(dt, label_val, null_thresholds_dt,
                                     use_threshold_ior = TRUE, use_threshold_ac = TRUE) {

  has_dynamic <- "dynamic" %in% names(dt)
  by_cols <- "triplet_id"
  if ("dynamic" %in% names(dt)) {
    by_cols <- c(by_cols, "dynamic", "t_ij")
  }
  if ("n_coadmin" %in% names(dt)) {
    by_cols <- c(by_cols, "n_coadmin")
  }
  if ("N" %in% names(dt)) {
    by_cols <- c(by_cols, "N")
  }

  expanded <- dt[, {
    stages <- unlist(stage)
    
    # GAM
    gam_log_ior <- unlist(log_ior)
    gam_log_ior_lower90 <- unlist(log_ior_lower90)
    gam_ac <- unlist(ac_values)
    gam_ac_lower90 <- unlist(ac_lower90)

    # Classical stratified metrics
    cls_log_ior <- unlist(log_ior_classic)
    cls_log_ior_lower90 <- unlist(log_ior_classic_lower90)
    cls_ac <- unlist(AC_classic)
    cls_ac_lower90 <- unlist(AC_classic_lower90)

    n <- min(length(stages), length(gam_log_ior), length(gam_log_ior_lower90),
             length(gam_ac), length(gam_ac_lower90),
             length(cls_log_ior), length(cls_log_ior_lower90),
             length(cls_ac), length(cls_ac_lower90))

    if (n > 0) {
      data.table(
        stage_num = stages[1:n],
        # GAM
        gam_log_ior = gam_log_ior[1:n],
        gam_log_ior_lower90 = gam_log_ior_lower90[1:n],
        gam_ac = gam_ac[1:n],
        gam_ac_lower90 = gam_ac_lower90[1:n],
        # Stratified
        classic_log_ior = cls_log_ior[1:n],
        classic_log_ior_lower90 = cls_log_ior_lower90[1:n],
        classic_ac = cls_ac[1:n],
        classic_ac_lower90 = cls_ac_lower90[1:n]
      )
    }
  }, by = by_cols]
  
  if (!has_dynamic) {
    expanded[, `:=`(dynamic = "control", t_ij = 0)]
  }
  
  expanded[, nichd := niveles_nichd[stage_num]]
  expanded[, label := label_val]
  
  # Merge with null distribution thresholds
  expanded <- merge(expanded, null_thresholds_dt, 
                   by.x = "stage_num", by.y = "stage", all.x = TRUE)
  
  # Store threshold usage parameters
  expanded[, `:=`(
    use_threshold_ior = use_threshold_ior,
    use_threshold_ac = use_threshold_ac
  )]
  return(expanded)
}

################################################################################
# Statistical power calculation function — GAM method
################################################################################

# Computes GAM-based statistical power over a 2D grid of (t_ij, n_coadmin) thresholds.
# At each grid point: power = TP / (TP + FN) at the triplet level (detected if ANY stage detects). 
# Returns the optimal threshold pair achieving target_power with maximum retention
# Returns the full power surface for visualization.

calculate_power_gam <- function(
  data_pos,
  target_power = 0.80,
  null_thresholds = NULL,
  metric_n = "n_coadmin",
  grid_resolution = 30,  
  use_threshold_ior = TRUE,   
  use_threshold_ac = TRUE,   
  detection = "double"  ) {    # "ior", "ac", or "double"
  
    library(data.table)

    ###########  
    # 1. Prepare positive data
    ###########
    pos_clean <- copy(data_pos)
  
    # Basic cleaning
    pos_clean <- pos_clean[is.finite(t_ij) & is.finite(get(metric_n))]
  
    # Merge thresholds if provided and specific columns are missing
    if (!is.null(null_thresholds)) {
      # Checks whether the threshold columns we need are missing
      need_merge <- FALSE
      if (use_threshold_ior && !"threshold_ior" %in% names(pos_clean)) need_merge <- TRUE
      if (use_threshold_ac && !"threshold_ac" %in% names(pos_clean)) need_merge <- TRUE
      if (need_merge) {
        # Ensures null_thresholds has a 'stage' column
        if (!"stage" %in% names(null_thresholds)) {
          stop("null_thresholds must have a column named 'stage'")
        }
        pos_clean <- merge(
          pos_clean,
          null_thresholds,
          by.x = "stage_num",
          by.y = "stage",
          all.x = TRUE
        )
      }
    }
  
  ###########
  # 2. 2D grid (t_ij x n_coadmin)
  ###########
  
  probs_grid <- seq(0, 0.95, length.out = grid_resolution)
  
  t_vals <- unique(quantile(pos_clean$t_ij, probs = probs_grid, na.rm = TRUE))
  n_vals <- unique(quantile(pos_clean[[metric_n]], probs = probs_grid, na.rm = TRUE))
  
  t_vals <- sort(unique(c(min(pos_clean$t_ij), t_vals)))
  n_vals <- sort(unique(c(min(pos_clean[[metric_n]]), n_vals)))
  
  search_grid <- CJ(t_threshold = t_vals, n_threshold = n_vals)
  
  ###########
  # 3. Power calculation
  ###########
  
  # Pre-computes detection flags according to the "detection" parameter
  
  # IOR criterion
  if (detection %in% c("ior", "double")) {
    if (use_threshold_ior) {
      if (!"threshold_ior" %in% names(pos_clean)) {
        stop("column threshold_ior not present in the data or in null_thresholds")
      }
      pos_clean[, ior_detected := (
        !is.na(gam_log_ior_lower90) & 
          gam_log_ior_lower90 > 0 & 
          gam_log_ior_lower90 > threshold_ior
      )]
    } else {
      pos_clean[, ior_detected := (
        !is.na(gam_log_ior_lower90) & gam_log_ior_lower90 > 0
      )]
    }
  }
  
  # AC criterion
  if (detection %in% c("ac", "double")) {
    if (use_threshold_ac) {
      if (!"threshold_ac" %in% names(pos_clean)) {
        stop("column threshold_ac not present in the data or in null_thresholds")
      }
      pos_clean[, ac_detected := (
        !is.na(gam_ac_lower90) &
          gam_ac_lower90 > 0 &
          gam_ac_lower90 > threshold_ac
      )]
    } else {
      pos_clean[, ac_detected := (
        !is.na(gam_ac_lower90) & gam_ac_lower90 > 0
      )]
    }
  }
  
  # Final detection flag based on the chosen mode
  if (detection == "ior") {
    pos_clean[, is_detected := ior_detected]
  } else if (detection == "ac") {
    pos_clean[, is_detected := ac_detected]
  } else {  # "double"
    pos_clean[, is_detected := ior_detected | ac_detected]
  }

  # Aggregates to triplet level: detected if ANY stage detects
  triplet_detection_gam <- pos_clean[, .(
    triplet_detected = any(is_detected, na.rm = TRUE),
    t_ij_triplet = unique(t_ij)[1],
    n_coadmin_triplet = unique(get(metric_n))[1]
  ), by = triplet_id]
  
  # Iterates over the search grid
  power_surface <- search_grid[, {
    
    # Triplet-level filter: t_ij >= t_thresh AND n >= n_thresh
    idx_subset <- which(
      triplet_detection_gam$t_ij_triplet >= t_threshold & 
      triplet_detection_gam$n_coadmin_triplet >= n_threshold
    )
    
    n_total <- length(idx_subset)
    
    if (n_total < 5) {
      list(tp = 0L, len = 0L, power = NA_real_)
    } else {
      n_tp <- sum(triplet_detection_gam$triplet_detected[idx_subset], na.rm = TRUE)
      list(
        tp = n_tp,
        len = n_total,
        power = n_tp / n_total
      )
    }
  }, by = .(t_threshold, n_threshold)]
  
  power_surface <- power_surface[!is.na(power)]
  
  ###########
  # 4. Identify optimal point
  ###########
  
  valid_configs <- power_surface[power >= target_power]
  
  if (nrow(valid_configs) == 0) {
    message(sprintf("Target power not reached (%.0f%%). Max: %.1f%%",
                    target_power*100, max(power_surface$power, na.rm=TRUE)*100))
    best_config <- power_surface[which.max(power)]
  } else {
    setorder(valid_configs, -len, power)
    best_config <- valid_configs[1]
  }
  
  t_star <- best_config$t_threshold
  n_star <- best_config$n_threshold
  achieved_power <- best_config$power
  n_retained <- best_config$len
  
  # Informational message based on detection mode
  detection_label <- switch(detection,
    "ior" = "IOR only",
    "ac" = "AC only",
    "double" = "IOR OR AC"
  )
  
  message(sprintf("\nOPTIMAL THRESHOLDS (GAM - %s):", detection_label))
  message(sprintf("  t_ij >= %.4f", t_star))
  message(sprintf("  %s >= %.1f", metric_n, n_star))
  message(sprintf("  Power achieved: %.1f%%", achieved_power * 100))
  message(sprintf("  Triplets retained: %d / %d (%.1f%%)",
                  n_retained,
                  uniqueN(triplet_detection_gam$triplet_id),
                  100 * n_retained / uniqueN(triplet_detection_gam$triplet_id)))
  
  ###########
  # 5. Build supersets
  ###########
  
  # Positive superset: ALL observations from triplets passing the filters
  triplets_passed <- triplet_detection_gam[
    t_ij_triplet >= t_star & n_coadmin_triplet >= n_star,
    triplet_id
  ]
  
  superset_pos <- pos_clean[triplet_id %in% triplets_passed]
  
  power_surface[, method := paste0("GAM-", toupper(detection))]
  
  return(list(
    power_surface = power_surface,
    t_star = t_star,
    n_star = n_star,
    superset_pos = superset_pos,
    achieved_power = achieved_power,
    criterion_type = "gam",
    detection_mode = detection,
    metric_n_used = metric_n
  ))
}

################################################################################
# Statistical power calculation function — Classical stratified method
################################################################################

# Computes statistical power by filtering on effect size and co-administration report count
#
# Classical methods are more prone to producing NAs or infinite CIs, so NA handling is parameterized
#
# Parameters:
# data_pos: data.table with expanded positive triplets
# data_neg: data.table with expanded negative triplets
# target_power: target power level (0.80)
# metric_n: column name to use as count filter ("n_coadmin" or "n_events")
# grid_resolution: number of steps for the search grid (default 30x30)
# detection: criterion for the superset — IOR, AC, or both
# na_remove: controls NA handling.
#  TRUE: superset includes ALL triplets meeting minimum detection criteria
#  FALSE: superset excludes NAs even if they meet the threshold criteria
#
# Return:
# power_surface: data.table with t_threshold, n_threshold, power, len
# t_star: optimal t_ij threshold
# n_star: optimal n_coadmin threshold
# superset_pos: positive triplets passing the optimal filters
# achieved_power: power achieved at the optimal point
#
# Implementation:
# See calculate_power_gam
# Key difference: power is computed at the stage level

calculate_power_classic <- function(
  data_pos,
  target_power = 0.80,
  null_thresholds = NULL,
  metric_n = "n_coadmin",
  grid_resolution = 30,
  detection = "double",  # "ior", "ac", or "double"
  na_remove = TRUE) {
    
  ###########
  # 1. Prepare positive data
  ###########

  detection <- match.arg(detection, choices = c("ior", "ac", "double"))

  pos_clean <- copy(data_pos)
  
  n_triplets_original_total <- uniqueN(pos_clean$triplet_id)
  n_obs_original_total <- nrow(pos_clean)

  # Counts before NA removal
  n_obs_before_na <- nrow(pos_clean)
  n_triplets_before_na <- uniqueN(pos_clean$triplet_id)

  # Removes NAs according to the detection criterion
  if (na_remove) {
    if (detection == "ior") {
      pos_clean <- pos_clean[!is.na(classic_log_ior_lower90)]
    } else if (detection == "ac") {
      pos_clean <- pos_clean[!is.na(classic_ac_lower90)]
    } else {  # double
      pos_clean <- pos_clean[!is.na(classic_log_ior_lower90) & !is.na(classic_ac_lower90)]
    }
  }
  # Logs how many were removed
  n_obs_after_na <- nrow(pos_clean)
  n_triplets_after_na <- uniqueN(pos_clean$triplet_id)
  n_triplets_lost_na <- n_triplets_before_na - n_triplets_after_na
 
  # Basic cleaning
  pos_clean <- pos_clean[is.finite(t_ij) & is.finite(get(metric_n))]
  
  if (!is.null(null_thresholds) && !"threshold" %in% names(pos_clean)) {
    pos_clean <- merge(
      pos_clean,
      null_thresholds,
      by.x = "stage_num",
      by.y = "stage",
      all.x = TRUE
    )
  }
  
  ###########
  # 2. 2D search grid
  ###########
  
  probs_grid <- seq(0, 0.95, length.out = grid_resolution)
  
  t_vals <- unique(quantile(pos_clean$t_ij, probs = probs_grid, na.rm = TRUE))
  n_vals <- unique(quantile(pos_clean[[metric_n]], probs = probs_grid, na.rm = TRUE))
  
  t_vals <- sort(unique(c(min(pos_clean$t_ij), t_vals)))
  n_vals <- sort(unique(c(min(pos_clean[[metric_n]]), n_vals)))
  
  search_grid <- CJ(t_threshold = t_vals, n_threshold = n_vals)
  
  ###########
  # 3. Power calculation
  ###########
  
  # Each row (triplet_id + stage) is an independent observation
  
  # Pre-compute detection flags according to the 'detection' parameter
  
  if (detection == "ior") {
    # IOR only: 90% CI lower bound > 0
    pos_clean[, is_detected := (
      !is.na(classic_log_ior_lower90) & classic_log_ior_lower90 > 0
    )]
    
  } else if (detection == "ac") {
    # AC (additive contrast) only: 90% CI lower bound > 0
    pos_clean[, is_detected := (
      !is.na(classic_ac_lower90) & classic_ac_lower90 > 0
    )]

  } else {  # detection == "double"
    # Double criterion: IOR (log-IOR CI > 0) OR AC (additive contrast CI > 0)
    pos_clean[, is_detected := (
      (!is.na(classic_log_ior_lower90) & classic_log_ior_lower90 > 0) |
      (!is.na(classic_ac_lower90) & classic_ac_lower90 > 0)
    )]
  }
    
  # Saves a copy with detection flags
  pos_all_with_detection <- copy(pos_clean)

  # Computed at the observation level (row): t_ij, n_coadmin, detection
  # Not aggregated to triplet level
  
  # Iterates over the search grid
  power_surface <- search_grid[, {
    
    # Observation-level filter: t_ij >= t_thresh AND n >= n_thresh
    idx_subset <- which(
      pos_clean$t_ij >= t_threshold & 
        pos_clean[[metric_n]] >= n_threshold
    )
    
    n_total <- length(idx_subset)
    
    if (n_total < 5) {
      list(tp = 0L, len = 0L, power = NA_real_)
    } else {
      n_tp <- sum(pos_clean$is_detected[idx_subset], na.rm = TRUE)
      list(
        tp = n_tp,
        len = n_total,
        power = n_tp / n_total
      )
    }
  }, by = .(t_threshold, n_threshold)]
  
  power_surface <- power_surface[!is.na(power)]
  
  ###########
  # 4. Identify optimal point
  ###########
  
  valid_configs <- power_surface[power >= target_power]
  
  if (nrow(valid_configs) == 0) {
    message(sprintf("Target power not reached (%.0f%%). Max: %.1f%%",
                    target_power*100, max(power_surface$power, na.rm=TRUE)*100))
    best_config <- power_surface[which.max(power)]
  } else {
    setorder(valid_configs, -len, power)
    best_config <- valid_configs[1]
  }
  
  t_star <- best_config$t_threshold
  n_star <- best_config$n_threshold
  achieved_power <- best_config$power
  n_retained <- best_config$len
  
  # Detection mode for the message
  detection_label <- switch(detection,
    "ior" = "IOR only",
    "ac" = "AC only",
    "double" = "IOR OR AC"
  )

  ###########
  # 5. Build supersets
  ###########
  
  # Positive superset: observations passing the filters
  superset_pos <- pos_all_with_detection[t_ij >= t_star & get(metric_n) >= n_star]
  
  # Superset metrics
  n_retained_total <- nrow(superset_pos)
  n_detected_total <- sum(superset_pos$is_detected, na.rm = TRUE)
  power_total <- ifelse(n_retained_total > 0, n_detected_total / n_retained_total, 0)

  # Count of unique triplets meeting the retention criterion
  n_triplets_retained <- uniqueN(superset_pos$triplet_id)
  n_triplets_total <- uniqueN(pos_all_with_detection$triplet_id)

  n_triplets_retained_vs_original <- n_triplets_retained
  pct_retained_vs_original <- 100 * n_triplets_retained_vs_original / n_triplets_original_total

  message(sprintf("\nThresholds (Stratified method, %s):", detection_label))
  message(sprintf("  t_ij >= %.4f", t_star))
  message(sprintf("  %s >= %.1f", metric_n, n_star))
  message(sprintf("  Power achieved: %.1f%%", achieved_power * 100))
  message(sprintf("  Triplets retained: %d / %d (%.1f%%)",
                  n_triplets_retained, n_triplets_total,
                  100 * n_triplets_retained / n_triplets_total))
  message(sprintf("  Power on full sample: %.1f%%", power_total * 100))

  # Lines to assess true retention, accounting for NAs excluded at the start of the function
  message(sprintf("  Triplets retained (vs original total): %d / %d (%.1f%%) [INCLUDES NA]",
                  n_triplets_retained_vs_original, n_triplets_original_total,
                  pct_retained_vs_original))
  
  power_surface[, method := paste0("Estratificado-", toupper(detection))]
  
  return(list(
    power_surface = power_surface,
    t_star = t_star,
    n_star = n_star,
    superset_pos = superset_pos,
    achieved_power = achieved_power,
    criterion_type = "classic",
    metric_n_used = metric_n,
    detection_mode = detection,
    na_remove = na_remove,
    achieved_power_total = power_total,  # power in full sample
    n_retained_grid = n_retained,  # retained in grid
    n_retained_total = n_retained_total,  # total retained
    n_detected_total = n_detected_total,  # total detected
    n_triplets_retained = n_triplets_retained,
    n_triplets_total = n_triplets_total
  ))
}

################################################################################
# Function to visualize the power surface over the full grid
################################################################################

# Generates a power surface heatmap with interpolation to a regular grid
# Faceted by method
#
# Parameters:
# power_result: list returned by calculate_power_gam() or calculate_power_classic()
# target_power: target power level (for visual reference)
# detection: detection type ("IOR", "AC", or "double")
# t_range: X-axis range (effect size)
# n_range: Y-axis range (co-administration report count)
# grid_size: interpolated grid size (default: 50x50 for smoothness)
# 
# Return:
# ggplot object with the full interpolated surface
#
# Implementation:
# Removes NAs and Inf values
# Interpolates to a regular grid using akima::interp (for missing data)
# Generates heatmap using geom_raster
# Returns optimal point statistics

plot_power_surface <- function(
  power_results_list,
  facet_by = "method", 
  target_power = 0.80, 
  detection = detection,
  t_range = c(0, 0.5),      
  n_range = c(0, 300),
  grid_size = 30) 
  {
  
  library(ggplot2)
  library(data.table)
  library(scales)
  library(akima)
  
  ###########
  # 1- Process multiple surfaces
  ###########
  
  all_surfaces <- rbindlist(lapply(names(power_results_list), function(met_name) {
    surface <- as.data.table(power_results_list[[met_name]]$power_surface)
    surface[, method_label := met_name]
    return(surface)
  }), fill = TRUE)
  
  # Extract optimal parameters for the subtitle
  opt_params <- lapply(names(power_results_list), function(met_name) {
    res <- power_results_list[[met_name]]
    data.table(
      method_label = met_name,
      t_star = res$t_star,
      n_star = res$n_star,
      achieved_power = res$achieved_power
    )
  })
  opt_params_dt <- rbindlist(opt_params)
  
  ###########
  # 2- Clean original data
  ###########
  
  # Remove NAs and non-finite values
  surface_clean <- all_surfaces[is.finite(t_threshold) & is.finite(n_threshold) & is.finite(power)]
  
  ###########
  # 3- Interpolate to a full regular grid
  ###########
  
  # Regular sequences with the same range for X and Y (to keep a square aspect)
  # Use the wider range across both axes to maintain proportion
  common_range <- c(0, max(t_range[2], n_range[2]))

  # Create regular grid sequences
  t_seq <- seq(t_range[1], t_range[2], length.out = grid_size)
  n_seq <- seq(n_range[1], n_range[2], length.out = grid_size)
  
  surfaces_interp <- lapply(unique(surface_clean$method_label), function(met) {
    surface_met <- surface_clean[method_label == met]

    # Normalizes each axis to [0, 1] before interpolating
    t_rng <- range(surface_met$t_threshold)
    n_rng <- range(surface_met$n_threshold)
    norm   <- function(v, r) if (diff(r) > 0) (v - r[1]) / diff(r) else v - r[1]
    denorm <- function(v, r) v * diff(r) + r[1]

    interp_result <- tryCatch({
      interp(
        x = norm(surface_met$t_threshold, t_rng),
        y = norm(surface_met$n_threshold, n_rng),
        z = surface_met$power,
        xo = norm(t_seq, t_rng),
        yo = norm(n_seq, n_rng),
        linear = TRUE,
        extrap = FALSE
      )
    }, error = function(e) NULL)

    if (is.null(interp_result)) return(NULL)

    dt <- data.table(expand.grid(
      t_threshold = denorm(interp_result$x, t_rng),
      n_threshold = denorm(interp_result$y, n_rng)
    ))
    dt[, power := as.vector(interp_result$z)]
    dt[, method_label := met]
    return(dt)
  })
  
  surface_plot <- rbindlist(surfaces_interp, fill = TRUE)
  surface_plot <- surface_plot[!is.na(power)]
  surface_plot[, power := pmin(pmax(power, 0), 1)]
  
  ###########
  # 5- Method labels
  ###########

  # Subtitle with per-method info
  subtitulo <- paste(opt_params_dt[, sprintf("%s: t>=%.4f, N>=%.0f (%.1f%%)", method_label, t_star, n_star, achieved_power*100)], collapse = " | ")

  method_label_clean <- switch(detection,
    "IOR" = "IOR detection",
    "AC" = "AC detection",
    detection
  )
  
  ###########
  # 6- Build the plot
  ###########
  
  p <- ggplot(surface_plot, aes(x = t_threshold, y = n_threshold, fill = power)) +
    
    # geom_raster() for regular grids
    geom_raster(interpolate = FALSE) +  # interpolate=TRUE smooths visually

    # Target-power boundary: 
    # the contour where power == target_power marks the >=80%-power region each method reaches
    geom_contour(aes(z = power), breaks = target_power, color = "white", linewidth = 0.6) +

    # Faceting
    facet_wrap(~ method_label, ncol = 2, scales = "fixed") +

    # Color scale
    scale_fill_viridis_c(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = percent_format(accuracy = 1), 
      name = "Power",
      option = "plasma",
      na.value = "white",
     guide = guide_colorbar(theme = theme( legend.text = element_text(angle = 45, hjust = 1, vjust = 1)))
    ) +
    
    # Axis scales
    scale_x_continuous(
      limits = t_range,
      expand = c(0, 0),
      breaks = pretty_breaks(n = 5)
    ) +
    
    scale_y_continuous(
      limits = n_range,
      expand = c(0, 0),
      breaks = pretty_breaks(n = 6)
    ) +

    # Square coordinate system
    coord_fixed(ratio = diff(t_range)/diff(n_range)) +

    # Labels
    labs(
      title = sprintf("Power Surface - %s", method_label_clean),
      subtitle = sprintf("Target: %.0f%% | %s", target_power * 100, subtitulo),
      x = expression("Effect Size " * (t[italic(ij)])),
      y = "Number of A-B reports",
      caption = sprintf("Grid: %dx%d cells", grid_size, grid_size)
    ) 
  return(p)
}

################################################################################
# Bootstrap-based metrics calculation function
################################################################################

# Computes classification metrics with 95% bootstrap confidence intervals
# 
# Parameters:
# dt: data.table with columns: triplet_id, detected, label
# n_boot: number of bootstrap replications
# aggregate_triplet: if TRUE, aggregates to triplet level using any(); if FALSE, keeps current granularity (stage or dynamic)
# 
# Return:
# data.table with metrics and 95% CIs

calculate_metrics <- function(dt, n_boot = 1000, aggregate_triplet = TRUE, score_type, score_type_auc) {
  
  # Aggregate to triplet level if requested
  if (aggregate_triplet) {
    dt_eval <- dt[, .(
      detected = any(detected, na.rm = TRUE),
      label = unique(label),
      score = max(get(score_type), na.rm = TRUE),
      score_auc = max(get(score_type_auc), na.rm = TRUE)
    ), by = triplet_id]
  } else {
    dt_eval <- dt[, .(triplet_id, detected, label, score = get(score_type), score_auc = get(score_type_auc))]
  }
  
  n_pos <- sum(dt_eval$label == 1)
  n_neg <- sum(dt_eval$label == 0)
  n_total <- nrow(dt_eval)
  
  auc_result <- tryCatch({
    # Filter NA and infinite values before computing ROC
    dt_roc <- dt_eval[is.finite(score_auc) & !is.na(score_auc) & !is.na(label)]
    
    if (nrow(dt_roc) < 2 || length(unique(dt_roc$label)) < 2) {
    stop("Insufficient data")
    }
    
    roc_data <- pROC::roc(
      response = dt_roc$label,
      predictor = dt_roc$score_auc,
      quiet = TRUE,
      direction = "<"
    )
    auc <- as.numeric(pROC::auc(roc_data))
    
    # Bootstrap for AUC
    b_idx_lab <- dt_roc$label
    b_idx_sc <- dt_roc$score_auc
    n_roc <- nrow(dt_roc)
    b_auc <- replicate(n_boot, {
      b_idx <- sample(n_roc, replace = TRUE)
      b_lab <- b_idx_lab[b_idx]
      b_sc <- b_idx_sc[b_idx]

      # Checks that both classes are present in the bootstrap sample
      n_p <- sum(b_lab == 1)
      n_n <- n_roc - n_p
      if (n_p == 0 || n_n == 0) {
        return(NA_real_)
      }

      r <- rank(b_sc)  # average ranks -> ties handled exactly as in pROC
      (sum(r[b_lab == 1]) - n_p * (n_p + 1) / 2) / (n_p * n_n)
    })
    
    list(
      auc = auc,
      auc_lower = unname(quantile(b_auc, 0.025, na.rm = TRUE)),
      auc_upper = unname(quantile(b_auc, 0.975, na.rm = TRUE))
    )
  }, error = function(e) {
    warning(sprintf("Error computing AUC: %s", e$message))
    list(auc = NA_real_, auc_lower = NA_real_, auc_upper = NA_real_)
  })
  
  # Bootstrap. 
  # Both triplet-level and row-level evaluation reduce to the same stratified resampling
  # Operating directly on the precomputed `detected` vectors avoids a per-replicate data.table join
  det_pos <- dt_eval[label == 1, detected]
  det_neg <- dt_eval[label == 0, detected]

  boot_stats <- replicate(n_boot, {
    b_tp <- sum(det_pos[sample.int(n_pos, n_pos, replace = TRUE)])
    b_fn <- n_pos - b_tp
    b_fp <- sum(det_neg[sample.int(n_neg, n_neg, replace = TRUE)])
    b_tn <- n_neg - b_fp

    b_sens <- ifelse((b_tp + b_fn) > 0, b_tp / (b_tp + b_fn), 0)
    b_spec <- ifelse((b_tn + b_fp) > 0, b_tn / (b_tn + b_fp), 0)
    b_ppv <- ifelse((b_tp + b_fp) > 0, b_tp / (b_tp + b_fp), 0)
    b_npv <- ifelse((b_tn + b_fn) > 0, b_tn / (b_tn + b_fn), 0)
    b_acc <- (b_tp + b_tn) / (n_pos + n_neg)
    b_f1 <- ifelse((b_tp + b_fp) > 0 && (b_tp + b_fn) > 0, 2 * b_tp / (2 * b_tp + b_fp + b_fn), 0)

    c(b_sens, b_spec, b_ppv, b_npv, b_acc, b_f1, b_tp, b_fn, b_fp, b_tn)
  })
  
  # Point estimates as bootstrap means (to avoid misalignment with CIs)
  sens_boot <- boot_stats[1, ]
  spec_boot <- boot_stats[2, ]
  ppv_boot <- boot_stats[3, ]
  npv_boot <- boot_stats[4, ]
  acc_boot <- boot_stats[5, ]
  f1_boot <- boot_stats[6, ]
  
  sens <- mean(sens_boot, na.rm = TRUE)
  spec <- mean(spec_boot, na.rm = TRUE)
  ppv <- mean(ppv_boot, na.rm = TRUE)
  npv <- mean(npv_boot, na.rm = TRUE)
  acc <- mean(acc_boot, na.rm = TRUE)
  f1 <- mean(f1_boot, na.rm = TRUE)
  
  # CI calculation
  calc_ci <- function(valores_boot) {
    c(
      point = mean(valores_boot, na.rm = TRUE),
      lower = unname(quantile(valores_boot, 0.025, na.rm = TRUE)),
      upper = unname(quantile(valores_boot, 0.975, na.rm = TRUE))
    )
  }
  
  sens_ci <- calc_ci(sens_boot)
  spec_ci <- calc_ci(spec_boot)
  ppv_ci <- calc_ci(ppv_boot)
  npv_ci <- calc_ci(npv_boot)
  acc_ci <- calc_ci(acc_boot)
  f1_ci <- calc_ci(f1_boot)
  
  # Counts from the original (non-bootstrap) dataset for reference
  tp_orig <- sum(dt_eval$detected & dt_eval$label == 1)
  fn_orig <- sum(!dt_eval$detected & dt_eval$label == 1)
  fp_orig <- sum(dt_eval$detected & dt_eval$label == 0)
  tn_orig <- sum(!dt_eval$detected & dt_eval$label == 0)
  
  data.table(
    sensitivity = sens_ci["point"], sensitivity_lower = sens_ci["lower"], sensitivity_upper = sens_ci["upper"],
    specificity = spec_ci["point"], specificity_lower = spec_ci["lower"], specificity_upper = spec_ci["upper"],
    PPV = ppv_ci["point"], PPV_lower = ppv_ci["lower"], PPV_upper = ppv_ci["upper"],
    NPV = npv_ci["point"], NPV_lower = npv_ci["lower"], NPV_upper = npv_ci["upper"],
    Accuracy = acc_ci["point"],
    F1 = f1_ci["point"], F1_lower = f1_ci["lower"], F1_upper = f1_ci["upper"],
    TP = tp_orig, FN = fn_orig, FP = fp_orig, TN = tn_orig,
    n_positives = n_pos, n_negatives = n_neg, n_total = n_total,
    AUC = auc_result$auc,
    AUC_lower = auc_result$auc_lower,
    AUC_upper = auc_result$auc_upper
  )
}

################################################################################
# ROC curve computation
################################################################################

# Computes ROC curve coordinates for a continuous interaction predictor.
#
# Signal classification elsewhere in the pipeline (detect_signal) is binary and driven by a fixed threshold 
# The ROC curve sweeps the full range of a continuous predictor, so it characterizes discrimination independently
# 
# Two predictors are of interest: the interaction point estimate and its 90% CI lower bound (the variable the fixed threshold acts on).
#
# Parameters:
#  dt: long-format data.table with a binary `label` column and the score column
#  score_col: name of the continuous predictor column
#  aggregate_triplet: if TRUE, collapse to triplet level using the max score, matching the triplet-level evaluation used for global metrics
#
# Return: 
# data.table with columns fpr, tpr, threshold and auc, or NULL when the data are insufficient

compute_roc_curve <- function(dt, score_col, aggregate_triplet = TRUE) {
  # 1. Reduces to one score-label pair per evaluated unit
  if (aggregate_triplet) {
    dt_eval <- dt[, .(
      label = unique(label),
      score = max(get(score_col), na.rm = TRUE)
    ), by = triplet_id]
  } else {
    dt_eval <- dt[, .(label, score = get(score_col))]
  }

  # 2. Drop non-finite scores; ROC requires both classes present
  dt_eval <- dt_eval[is.finite(score) & !is.na(label)]
  if (nrow(dt_eval) < 2 || length(unique(dt_eval$label)) < 2) {
    return(NULL)
  }

  # 3. Fit ROC; direction "<" because a higher score means stronger interaction
  roc_obj <- pROC::roc(
    response = dt_eval$label,
    predictor = dt_eval$score,
    quiet = TRUE,
    direction = "<"
  )

  # 4. Extract every threshold coordinate as plotting-ready long format
  co <- pROC::coords(
    roc_obj,
    x = "all",
    ret = c("threshold", "specificity", "sensitivity"),
    transpose = FALSE
  )
  setDT(co)

  data.table(
    threshold = co$threshold,
    fpr = 1 - co$specificity,  # false positive rate
    tpr = co$sensitivity,      # true positive rate
    auc = as.numeric(pROC::auc(roc_obj))
  )
}

################################################################################
# Data expansion function with reduction
################################################################################

# Expands data to long format
#
# Implementation:
# Loads data and creates objects according to the reduction level
# Filters out failed injections
# Expands to long format
# Filters to stages with high reporting according to the dynamic

expand <- function(red_pct) {
  suffix_file <- if(red_pct == 0) "" else paste0("_", red_pct)  # suffix based on reduction level
  
  ruta_pos <- paste0(ruta_base_sensitivity, "positive_triplets_results", suffix_file, ".rds")
  ruta_neg <- paste0(ruta_base_sensitivity, "negative_triplets_results", suffix_file, ".rds")
  
  pos_raw <- readRDS(ruta_pos)
  neg_raw <- readRDS(ruta_neg)
  
  # Successful injections only
  pos_valid <- pos_raw[injection_success == TRUE]
  
  # Expand to long format
  pos_exp <- expand_clean_all_metrics(pos_valid, 1, null_thresholds, use_threshold_ior, use_threshold_ac)
  neg_exp <- expand_clean_all_metrics(neg_raw, 0, null_thresholds, use_threshold_ior, use_threshold_ac)
  
  # Merge with stage classification
  pos_exp <- merge(pos_exp, stage_class, by = c("nichd", "dynamic"), all.x = TRUE)
  neg_exp[, class := 0]
  
  # Merge with co-administration data
  pos_exp <- merge(pos_exp, coadmin_stage_pos[, .(triplet_id, stage_num, n_coadmin_stage)], by = c("triplet_id", "stage_num"), all.x = TRUE)
  neg_exp <- merge(neg_exp, coadmin_stage_neg[, .(triplet_id, stage_num, n_coadmin_stage)], by = c("triplet_id", "stage_num"), all.x = TRUE)
  
  # High-reporting dataset (excluding uniform dynamic)
  pos_high <- pos_exp[class == 1]
  neg_high <- neg_exp[class == 0]
  
  list(
    pos_all = pos_exp,
    pos_high = pos_high,
    neg_high = neg_high,
    reduction_pct = red_pct
  )
}

################################################################################
# Signal detection function
################################################################################

# Detects signals according to the specified method
#
# Parameters:
# data: expanded triplet results (long format)
# thresholds: per-stage thresholds
# method: detection method ("gam" or "classic")
# criterion: detection criterion ("ior", "ac", or "double")
# use_null: if TRUE, applies null distribution thresholds
#

detect_signal <- function(dt, method_name, detection_type, use_null) {
  
  is_gam <- grepl("GAM", method_name) # checks for "GAM" in the method string
  
  # Determine columns based on the method
  if (is_gam) {
    ior_col <- "gam_log_ior_lower90"
    ac_col <- "gam_ac_lower90"     # additive contrast, 90% CI lower bound
    thresh_ior_col <- "threshold_ior"  # null distribution thresholds
    thresh_ac_col <- "threshold_ac"
  } else {
    ior_col <- "classic_log_ior_lower90"
    ac_col <- "classic_ac_lower90"  # additive contrast, 90% CI lower bound
    # Classic null thresholds. populated when 20_null runs classic methods
    thresh_ior_col  <- "threshold_classic_ior"
    thresh_ac_col <- "threshold_classic_ac"
  }

  # Compute detection flags by criterion type
  # use_null applies to both GAM (threshold_ior/ac) and classic (threshold_classic_ior/ac)
  if (detection_type == "IOR") {
    if (use_null && thresh_ior_col %in% names(dt)) {
      dt[, detected := !is.na(get(ior_col)) & get(ior_col) > 0 & get(ior_col) > get(thresh_ior_col)]
    } else {
      dt[, detected := !is.na(get(ior_col)) & get(ior_col) > 0]}
  } else if (detection_type == "AC") {
    if (use_null && thresh_ac_col %in% names(dt)) {
      dt[, detected := !is.na(get(ac_col)) & get(ac_col) > 0 & get(ac_col) > get(thresh_ac_col)]
    } else {
      dt[, detected := !is.na(get(ac_col)) & get(ac_col) > 0]}
  } else {
    # Double criterion: IOR OR AC
    ior_det <- if (use_null && thresh_ior_col %in% names(dt)) {
      !is.na(dt[[ior_col]]) & dt[[ior_col]] > 0 & dt[[ior_col]] > dt[[thresh_ior_col]]
    } else {
      !is.na(dt[[ior_col]]) & dt[[ior_col]] > 0}
    ac_det <- if (use_null && thresh_ac_col %in% names(dt)) {
      !is.na(dt[[ac_col]]) & dt[[ac_col]] > 0 & dt[[ac_col]] > dt[[thresh_ac_col]]
    } else {
      !is.na(dt[[ac_col]]) & dt[[ac_col]] > 0}
    dt[, detected := ior_det | ac_det]
  }
  
  # Replace NA with FALSE (i.e., not detected)
  dt[is.na(detected), detected := FALSE]
  
  return(dt)
}
