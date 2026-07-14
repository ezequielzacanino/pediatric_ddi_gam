################################################################################
# Semi-synthetic data generation script with sensitivity analysis
# Script 10_augmentation
################################################################################

source("00_functions.R", local = TRUE)

################################################################################
# Configuration
################################################################################

# Cohort sizes and signal-injection parameters
n_pos <- 1000
n_neg <- 10000
lambda_fc <- 0.75
dinamicas <- c("uniform","increase","decrease","plateau","inverse_plateau")

# Minimum data-density requirements for triplet inclusion
min_reports_triplet <- 2
min_nichd_with_rep <- 2
all_nichd_rep <- FALSE
max_events_per_pair <- 10000

# Sensitivity analysis: progressively reduce available reports to stress-test detection power
reduction_levels <- seq(10, 90, by = 10)  # 10%, 20% ... 90%

# Parallel batch parameters for positive triplets
batch_size_pos <- 25
save_interval <- 1   # checkpoint frequency (every N batches)

n_null_reports <- 100000  

output_dir <- paste0("./results/", suffix, "/augmentation_results/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Run-control flags (resume helpers)
################################################################################

# Reuse the cached triplet construction (triplets_dt, trip_summary) from a prior run
reuse_triplet_cache <- TRUE
triplet_cache_file  <- paste0(output_dir, "triplet_cache.rds")

# When FALSE, skip the positive candidate construction 
# reload pos_meta and positive scores from disk 
# Lets the negative section run without repeating the positive one.
run_positives <- FALSE

################################################################################
# Data loading
################################################################################

ade_raw_dt <- fread(ruta_ade_raw)

# Sex covariate is optional; normalize and factor only when used in the GAM formula
if (include_sex) {
  ade_raw_dt[, sex := toupper(trimws(sex))]
  ade_raw_dt[sex == "M", sex := "MALE"]
  ade_raw_dt[sex == "F", sex := "FEMALE"]
  ade_raw_dt[, sex := factor(sex, levels = c("MALE", "FEMALE"))]
  
  sex_summary <- ade_raw_dt[, .(n = .N), by = sex]
  message("\nSex distribution:")
  print(sex_summary)
}

message(sprintf("Dataset %s rows", format(nrow(ade_raw_dt), big.mark = ",")))

################################################################################
# Preprocessing 
################################################################################

# Map every ATC concept to its canonical ID 
# drugs with different formulation codes but the same active compound as the same entity
translation_table <- build_drug_translation_table()

cat(sprintf("Original IDs: %d\n", uniqueN(translation_table$atc_concept_id)))
cat(sprintf("Unique IDs: %d\n", uniqueN(translation_table$canonical_id)))

ade_raw_dt[, atc_concept_id := as.character(atc_concept_id)]

ade_raw_dt <- merge(
  ade_raw_dt, 
  translation_table[, .(atc_concept_id, canonical_id)], 
  by = "atc_concept_id", 
  all.x = TRUE
)

ade_raw_dt[!is.na(canonical_id), atc_concept_id := canonical_id]
ade_raw_dt[, canonical_id := NULL]

nrow_before <- nrow(ade_raw_dt)
ade_raw_dt <- unique(ade_raw_dt, by = c("safetyreportid", "atc_concept_id", "meddra_concept_id"))

ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
ade_raw_dt[, nichd_num := as.integer(nichd)]

# Keep only the columns used downstream. PSOCK workers receive a serialized copy of this table
keep_cols <- c("safetyreportid", "atc_concept_id", "meddra_concept_id", "nichd", "nichd_num")
if (include_sex) keep_cols <- c(keep_cols, "sex")
ade_raw_dt <- ade_raw_dt[, ..keep_cols]

# Index the columns repeatedly used for == subsetting so per-triplet lookups are binary searches
setindex(ade_raw_dt, atc_concept_id)
setindex(ade_raw_dt, meddra_concept_id)

################################################################################
# Candidate triplet construction
################################################################################

# Reuse the cached triplet construction when available
if (reuse_triplet_cache && file.exists(triplet_cache_file)) {
  message("Loading triplet construction from cache: ", triplet_cache_file)
  triplet_cache <- readRDS(triplet_cache_file)
  triplets_dt  <- triplet_cache$triplets_dt
  trip_summary <- triplet_cache$trip_summary
  rm(triplet_cache); gc()
} else {

# Extract unique drugs and events per report
drugs_by_report <- unique(ade_raw_dt[, .(safetyreportid, atc_concept_id)])
events_by_report <- unique(ade_raw_dt[, .(safetyreportid, meddra_concept_id)])

reports <- unique(ade_raw_dt[, .(safetyreportid, nichd, nichd_num)])
drugs_list <- drugs_by_report[, .(drugs = list(unique(atc_concept_id))), by = safetyreportid]
events_list <- events_by_report[, .(events = list(unique(meddra_concept_id))), by = safetyreportid]

reports_meta <- merge(reports, drugs_list, by = "safetyreportid", all.x = TRUE)
reports_meta <- merge(reports_meta, events_list, by = "safetyreportid", all.x = TRUE)

report_combo <- copy(reports_meta)  # copy to avoid modifying reports_meta in the lapply below

triplets_list <- pblapply(seq_len(nrow(report_combo)), function(i) {
  rowi <- report_combo[i]
  make_triplets(
    drug = rowi$drugs[[1]], 
    event = rowi$events[[1]], 
    report_id = rowi$safetyreportid, 
    nichd_stage = rowi$nichd
  )
})

triplets_dt <- rbindlist(triplets_list, use.names = TRUE)
rm(triplets_list); gc()

trip_counts_by_stage <- unique(triplets_dt[, .(drugA, drugB, meddra, nichd_num, safetyreportid)])[
  , .N, by = .(drugA, drugB, meddra, nichd_num)
]

trip_summary <- trip_counts_by_stage[, .(
  N = sum(N),
  n_stages = uniqueN(nichd_num),
  stages_with_data = list(nichd_num)
), by = .(drugA, drugB, meddra)]

# Persist the expensive triplet construction for reuse in later runs
saveRDS(list(triplets_dt = triplets_dt, trip_summary = trip_summary),
        triplet_cache_file)
message("Triplet construction cached to: ", triplet_cache_file)
}

# Pre-compute unique report-drug-stage triples once
coadmin_precomp <- unique(
  ade_raw_dt[, .(safetyreportid, atc_concept_id, nichd_num)]
)
setkey(coadmin_precomp, atc_concept_id)

################################################################################
# Positive cohort: candidate filtering, selection and batch modelling.
# Skipped when run_positives is FALSE
################################################################################

if (run_positives) {

if (all_nichd_rep) {
  # filter: all seven NICHD stages must have at least one report
  candidatos_pos <- trip_summary[
    N >= min_reports_triplet &
    n_stages == 7
  ]
  message(sprintf("Strict filter: all stages must have reports"))
} else {
  candidatos_pos <- trip_summary[
    N >= min_reports_triplet &
    n_stages >= min_nichd_with_rep
  ]
  message(sprintf("Flexible filter: at least %d stages with reports", min_nichd_with_rep))
}

# Capture funnel counts before filtering so they can be persisted for the manuscript
message("\nFiltering statistics:")
message(sprintf("  Total initial candidates: %s",
                format(nrow(trip_summary), big.mark = ",")))
message(sprintf("  Candidates meeting report filter (N >= %d): %s",
                min_reports_triplet,
                format(sum(trip_summary$N >= min_reports_triplet), big.mark = ",")))
message(sprintf("  Candidates meeting stage filter: %s",
                format(nrow(candidatos_pos), big.mark = ",")))

n_candidates_total <- nrow(trip_summary)
n_candidates_reports <- sum(trip_summary$N >= min_reports_triplet)
n_candidates_stages <- nrow(candidatos_pos)

stage_distribution <- candidatos_pos[, .N, by = n_stages][order(n_stages)]
message("\nCandidate distribution by number of stages:")
print(stage_distribution)

fwrite(stage_distribution, paste0(output_dir, "positive_candidates_by_n_stages.csv"))

# Cap events per drug pair to prevent any one pair from dominating the positive set
candidatos_pos <- candidatos_pos[,
  .SD[sample(.N, min(.N, max_events_per_pair))],
  by = .(drugA, drugB)
]
message(sprintf("Triplets after diversification: %d", nrow(candidatos_pos)))

n_pos_final <- min(n_pos, nrow(candidatos_pos))
positivos_sel <- candidatos_pos[sample(.N, n_pos_final)]

# Persist the positive-cohort construction funnel for manuscript reproducibility
positive_construction_funnel <- data.table(
  step = c("candidate_triplets_total", "meet_min_reports", "meet_min_stages",
           "post_diversification", "selected_positives_base_pairs"),
  n = c(n_candidates_total, n_candidates_reports, n_candidates_stages,
        nrow(candidatos_pos), n_pos_final)
)
fwrite(positive_construction_funnel, paste0(output_dir, "positive_construction_funnel.csv"))

message(sprintf("\nSelected positive triplets: %d", nrow(positivos_sel)))
message(sprintf("  Mean reports: %.1f", mean(positivos_sel$N)))
message(sprintf("  Mean stages: %.1f", mean(positivos_sel$n_stages)))

################################################################################
# Positive triplet selection
################################################################################

pos_meta_base <- positivos_sel[, .(drugA, drugB, meddra, N)]
pos_meta_base[, base_triplet_id := 1:.N]

pos_meta <- pos_meta_base[, {
  data.table(
    drugA = drugA,
    drugB = drugB,
    meddra = meddra,
    N = N,
    base_triplet_id = base_triplet_id,
    dynamic = dinamicas
  )
}, by = base_triplet_id]

pos_meta[, fold_change := fold_change(.N, lambda = lambda_fc)]
pos_meta[, triplet_id := 1:.N]

# Flag the 30 highest-report triplets for diagnostic plots
top30_ids <- pos_meta[order(-N)][1:min(.N, 30), triplet_id]
pos_meta[, is_top30 := triplet_id %in% top30_ids]

fwrite(pos_meta, paste0(output_dir, "positive_triplets_metadata.csv"))

# Also persist as RDS so a negatives-only run can reload pos_meta
saveRDS(pos_meta, paste0(output_dir, "pos_meta.rds"))

################################################################################
# Co-administration counts by NICHD stage for positive triplets
################################################################################

coadmin_stage_pos <- compute_coadmin_batch(
  pairs_dt = pos_meta[, .(triplet_id, drugA, drugB, meddra)],
  ade_dt = coadmin_precomp
)
setcolorder(coadmin_stage_pos, c("triplet_id", "drugA", "drugB", "meddra", "nichd_num", "nichd", "n_coadmin_stage"))

fwrite(coadmin_stage_pos, paste0(output_dir, "positive_coadmin_by_stage.csv"))

rm(coadmin_stage_pos)
gc()

################################################################################
# Batch processing of positive triplets with sensitivity analysis
################################################################################

# Precompute the per-level downsampling keep-indices once 
# The reduction is deterministic given the fixed seed 
# Every positive triplet reuses the same reduced datasets 
# workers then subset ade_raw_dt by these indices instead of re-sampling the full table per triplet and per level (C1)
  
rng_state_before <- .Random.seed
reduced_idx_list <- setNames(
  lapply(reduction_levels, function(rp) reduce_indices_by_stage(ade_raw_dt, rp, seed = 7113)),
  as.character(reduction_levels)
)
.Random.seed <- rng_state_before

# Set up parallel cluster for batch processing.

start_pos_cluster <- function() {
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

    clusterExport(cl, c("process_single_positive", "fit_reduced_model",
    "reduced_idx_list", "inject_signal", "fit_gam", "build_eval_table",
    "calculate_classic_ior", "calculate_classic_ac", "calc_basic_counts",
    "pos_meta", "ade_raw_dt", "reduction_levels", "niveles_nichd",
    "spline_individuales", "include_sex", "include_stage_sex",
    "k_spline", "bs_type", "select", "nichd_spline", "include_nichd", "Z90",
    "classic_continuity_correction", "continuity_correction_value", "ac_bootstrap_n",
    "generate_dynamic", "fold_change"),
    envir = parent.frame())

  clusterEvalQ(cl, {
    library(data.table)
    library(mgcv)
    library(MASS)
  })
  cl
}

# Recreate the worker pool every N batches to release accumulated mgcv memory.
cluster_recycle_interval <- 2

cl <- start_pos_cluster()

n_batches <- ceiling(nrow(pos_meta) / batch_size_pos)

# Resume support: detect existing positive checkpoints
existing_cp_pos <- list.files(
  output_dir,
  pattern = "checkpoint_positives_batch_\\d+\\.rds",
  full.names = TRUE
)

# Re-run the last checkpoint batch in case it was truncated by a crash
if (length(existing_cp_pos) > 0) {
  cp_batches_pos <- sort(as.integer(
    gsub(".*checkpoint_positives_batch_(\\d+)\\.rds", "\\1", existing_cp_pos)
  ))
  last_cp_pos <- max(cp_batches_pos)
  redo_batch <- last_cp_pos
  load_batches <- cp_batches_pos[cp_batches_pos < redo_batch]

  message(sprintf(
    "Last checkpoint detected: batch %d",
    last_cp_pos
  ))

  if (length(load_batches) > 0) {
    positives_scores_list <- setNames(
      lapply(load_batches, function(b)
        readRDS(paste0(output_dir, "checkpoint_positives_batch_", b, ".rds"))
      ),
      as.character(load_batches)
    )
  } else {
    positives_scores_list <- list()
  }

  start_batch_pos <- redo_batch
} else {
  positives_scores_list <- list()
  start_batch_pos <- 1
}


for (batch in start_batch_pos:n_batches) {

  # Recycle the worker pool before processing this batch (except the first one)
  if (batch > start_batch_pos &&
      (batch - start_batch_pos) %% cluster_recycle_interval == 0) {
    stopCluster(cl)
    gc(verbose = FALSE)
    cl <- start_pos_cluster()
  }

  start_idx <- (batch - 1) * batch_size_pos + 1
  end_idx <- min(batch * batch_size_pos, nrow(pos_meta))
  batch_indices <- start_idx:end_idx

  message(sprintf("Batch %d / %d (triplets %d-%d)", batch, n_batches, start_idx, end_idx))
  
  # Process batch in parallel. 
  # Reproducibility comes from the explicit per-triplet set.seed inside process_single_positive
  batch_results <- foreach(
    idx = batch_indices,
    .packages = c("data.table", "mgcv"),
    .errorhandling = "pass",
    .verbose = FALSE
  ) %dopar% {
    process_single_positive(
      idx, pos_meta, ade_raw_dt, reduction_levels,
      reduced_idx_list,
      spline_individuales, include_sex, include_stage_sex,
      k_spline, bs_type, select, nichd_spline,
      base_seed = 7113)
  }
  
  batch_results_clean <- Filter(function(x) !inherits(x, "error"), batch_results)

  if (length(batch_results_clean) > 0) {
    batch_dt <- rbindlist(batch_results_clean, fill = TRUE)
    positives_scores_list[[batch]] <- batch_dt
    message(sprintf("Batch %d done: %d successful triplets (base)", batch, sum(batch_dt$model_success & batch_dt$reduction_pct == 0, na.rm = TRUE)))

    if (batch %% save_interval == 0 || batch == n_batches) {
      checkpoint_file <- paste0(output_dir, "checkpoint_positives_batch_", batch, ".rds")
      saveRDS(batch_dt, checkpoint_file)
      message(sprintf("Checkpoint saved: %s", checkpoint_file))
    }
  }
  rm(batch_results, batch_results_clean, batch_dt)
  gc(verbose = FALSE)
}
stopCluster(cl)

positives_scores <- rbindlist(positives_scores_list, fill = TRUE)
rm(positives_scores_list)
gc()

message(sprintf("Total triplets processed: %d", nrow(positives_scores[reduction_pct == 0])))
message(sprintf("Successful (base): %d", sum(positives_scores$model_success & positives_scores$reduction_pct == 0, na.rm = TRUE)))

# Write one RDS and one flat CSV per reduction level for downstream scripts
for (red_pct in c(0, reduction_levels)) {
  suffix_file <- if (red_pct == 0) "" else paste0("_", red_pct)

  subset_data <- positives_scores[reduction_pct == red_pct]

  if (nrow(subset_data) > 0) {
    saveRDS(subset_data, paste0(output_dir, "positive_triplets_results", suffix_file, ".rds"))

    subset_csv <- copy(subset_data)

    if ("diagnostics" %in% names(subset_csv)) {
      subset_csv[, diagnostics := NULL]
    }

    # Flatten list columns to comma-separated strings for CSV compatibility
    list_cols <- names(subset_csv)[sapply(subset_csv, is.list)]
    
    for (col in list_cols) {
      subset_csv[, (col) := sapply(get(col), function(x) {
        if (is.null(x) || length(x) == 0) return(NA_character_)
        paste(x, collapse = ",")
      })]
    }
    fwrite(subset_csv, paste0(output_dir, "positive_triplets_results", suffix_file, ".csv"))
  }
}

gc()

message("\nSuccessful positives (base): ", sum(positives_scores$model_success & positives_scores$injection_success & positives_scores$reduction_pct == 0, na.rm = TRUE))

} else {
  # run_positives is FALSE: reload pos_meta and positives_scores from disk.
  # pos_meta drives negative pool construction (drug/event pools and positive exclusion);
  # positives_scores is needed for the final comparison, summary and sensitivity blocks.
  message("Skipping positive section: loading pos_meta and positive scores")
  pos_meta <- readRDS(paste0(output_dir, "pos_meta.rds"))

  # Reassemble positives_scores across all reduction levels from the per-level RDS
  pos_result_files <- list.files(
    output_dir,
    pattern = "positive_triplets_results(_\\d+)?\\.rds$",
    full.names = TRUE
  )
  positives_scores <- rbindlist(lapply(pos_result_files, readRDS), fill = TRUE)
}

################################################################################
# Negative triplet selection
################################################################################

# Restrict the negative pool to the same drug and event vocabulary
drugs_from_pos <- unique(c(pos_meta$drugA, pos_meta$drugB))
events_from_pos <- unique(pos_meta$meddra)

message(sprintf("Drug pool: %d", length(drugs_from_pos)))
message(sprintf("Event pool: %d", length(events_from_pos)))

# Build a lookup set of positive triplet identifiers to exclude from the negative pool 
# prevents contamination between the two sets
pos_triplet_ids <- paste(
  pmin(pos_meta$drugA, pos_meta$drugB),
  pmax(pos_meta$drugA, pos_meta$drugB),
  pos_meta$meddra,
  sep = "_"
)
pos_triplet_set <- unique(pos_triplet_ids)

# Enumerate all drug-pair x event combinations in chunks
# avoids a single cross-join that would exhaust memory for large pools
chunk_size <- 50
n_drugs <- length(drugs_from_pos)
n_chunks <- ceiling(n_drugs / chunk_size)

candidatos_neg_list <- list()
pb <- txtProgressBar(max = n_chunks, style = 3)

for (chunk_idx in 1:n_chunks) {
  start_idx <- (chunk_idx - 1) * chunk_size + 1
  end_idx <- min(chunk_idx * chunk_size, n_drugs)
  
  drugs_chunk <- drugs_from_pos[start_idx:end_idx]
  
  chunk_combinations <- CJ(
    drugA = drugs_chunk,
    drugB = drugs_from_pos,
    meddra = events_from_pos
  )
  # Enforce canonical ordering so (A,B) and (B,A) are the same pair
  chunk_combinations[, `:=`(
    drugA_ord = pmin(drugA, drugB),
    drugB_ord = pmax(drugA, drugB)
  )]
  chunk_combinations <- chunk_combinations[
    drugA == drugA_ord
  ]
  chunk_combinations[, `:=`(
    drugA = drugA_ord,
    drugB = drugB_ord,
    drugA_ord = NULL,
    drugB_ord = NULL
  )]

  chunk_combinations <- chunk_combinations[drugA != drugB]
  chunk_combinations[, triplet_id := paste(drugA, drugB, meddra, sep = "_")]
  # Exclude triplets already used as positives to prevent set overlap
  chunk_combinations <- chunk_combinations[!triplet_id %in% pos_triplet_set]
  # Keep only triplets observed in the data with sufficient reports
  chunk_candidates <- merge(
    chunk_combinations[, .(drugA, drugB, meddra, triplet_id)],
    trip_summary,
    by = c("drugA", "drugB", "meddra"),
    all.x = FALSE
  )

  chunk_candidates <- chunk_candidates[N >= min_reports_triplet]
  
  if (nrow(chunk_candidates) > 0) {
    candidatos_neg_list[[chunk_idx]] <- chunk_candidates
  }
  
  rm(chunk_combinations, chunk_candidates)
  gc(verbose = FALSE)
  
  setTxtProgressBar(pb, chunk_idx)
}
close(pb)

candidatos_neg <- rbindlist(candidatos_neg_list, use.names = TRUE)
rm(candidatos_neg_list)
gc()

candidatos_neg <- unique(candidatos_neg, by = "triplet_id")

# Count reports per NICHD stage for each negative candidate
# then apply the same minimum-stages filter used for the positive candidates
neg_ids <- paste(candidatos_neg$drugA, candidatos_neg$drugB, candidatos_neg$meddra, sep = "_")

neg_counts_by_stage <- unique(
  triplets_dt[paste(drugA, drugB, meddra, sep = "_") %in% neg_ids, 
              .(drugA, drugB, meddra, nichd_num, safetyreportid)]
)[, .(n_reports = .N), by = .(drugA, drugB, meddra, nichd_num)]

# triplets_dt holds every triplet for all 684k reports (several GB)
# Frees it now so it does not compete with the worker copies during the memory-heavy negative GAM passes below.
rm(triplets_dt)
gc(verbose = FALSE)

setnames(neg_counts_by_stage, "nichd_num", "nichd")

# nichd column holds stage names here; convert to integer index for the filter join
nichd_to_num <- setNames(1:7, niveles_nichd)
neg_counts_by_stage[, nichd_num := nichd_to_num[nichd]]

stage_counts_summary <- neg_counts_by_stage[, .(
  n_stages_with_data = uniqueN(nichd_num[n_reports > 0])
), by = .(drugA, drugB, meddra)]

setkey(candidatos_neg, drugA, drugB, meddra)
setkey(stage_counts_summary, drugA, drugB, meddra)

candidatos_neg_full <- merge(
  candidatos_neg,
  stage_counts_summary,
  by = c("drugA", "drugB", "meddra"),
  all.x = TRUE
)
candidatos_neg_filtered <- candidatos_neg_full[n_stages_with_data >= min_nichd_with_rep]

message(sprintf("\nFiltered negative candidates: %d (of %d)", nrow(candidatos_neg_filtered), nrow(candidatos_neg)))

n_neg_final <- min(n_neg, nrow(candidatos_neg_filtered))
selected_negatives <- candidatos_neg_filtered[sample(.N, n_neg_final)]

selected_negatives[, triplet_id := 1:.N]

# Persist the negative-cohort construction funnel for manuscript reproducibility
negative_construction_funnel <- data.table(
  step = c("drug_pool", "event_pool", "candidate_triplets_in_data",
           "meet_min_stages", "selected_negatives"),
  n = c(length(drugs_from_pos), length(events_from_pos), nrow(candidatos_neg),
        nrow(candidatos_neg_filtered), n_neg_final)
)
fwrite(negative_construction_funnel, paste0(output_dir, "negative_construction_funnel.csv"))

message(sprintf("\nSelected negatives: %d", nrow(selected_negatives)))
message(sprintf("  Mean reports: %.1f", mean(selected_negatives$N)))
message(sprintf("  Mean stages: %.1f", mean(selected_negatives$n_stages_with_data)))

# The candidate pools are no longer needed once selected_negatives 
rm(candidatos_neg, candidatos_neg_full, candidatos_neg_filtered,
   neg_counts_by_stage, stage_counts_summary, neg_ids)
gc(verbose = FALSE)

################################################################################
# Co-administration counts by NICHD stage for negative triplets
################################################################################

coadmin_stage_neg <- compute_coadmin_batch(
  pairs_dt = selected_negatives[, .(triplet_id, drugA, drugB, meddra)],
  ade_dt = coadmin_precomp
)
setcolorder(coadmin_stage_neg, c("triplet_id", "drugA", "drugB", "meddra", "nichd_num", "n_coadmin_stage"))

fwrite(coadmin_stage_neg, paste0(output_dir, "negative_coadmin_by_stage.csv"))

# Frees the co-admin table and precomputed index
rm(coadmin_stage_neg, coadmin_precomp)
gc()

# Reduced datasets are built one level at a time 

################################################################################
# Batch processing of negative triplets
################################################################################

batch_size_neg <- 20
n_batches <- ceiling(nrow(selected_negatives) / batch_size_neg)

# Recreates the worker pool every N batches to release memory
# Each restart re-exports the full ade_pass table to every worker
cluster_recycle_interval_neg <- 25

# Builds a fresh cluster for a negative pass and exports everything the workers need
start_neg_cluster <- function(ade_pass) {
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  clusterExport(cl, c(
    "fit_gam", "build_eval_table", "Z90", "niveles_nichd", "selected_negatives",
    "spline_individuales", "include_sex", "include_stage_sex",
    "k_spline", "nichd_spline", "include_nichd", "bs_type", "select", "method",
    "classic_continuity_correction", "continuity_correction_value", "ac_bootstrap_n",
    "calculate_classic_ior", "calculate_classic_ac", "calc_basic_counts"
  ), envir = .GlobalEnv)
  clusterEvalQ(cl, {
    library(data.table)
    library(mgcv)
    library(MASS)
  })
  # ade_pass is a function argument, not a global; uses the local environment
  clusterExport(cl, "ade_pass", envir = environment())
  cl
}

# Returns the batch number to resume from (last saved checkpoint, or 1 if none)
detect_neg_checkpoints <- function(out_dir, pass_tag) {
  pattern <- sprintf("checkpoint_neg_%s_batch_\\d+\\.rds", pass_tag)
  cp_files <- list.files(out_dir, pattern = pattern, full.names = TRUE)

  if (length(cp_files) == 0) return(1L)

  cp_nums <- sort(as.integer(
    gsub(sprintf(".*checkpoint_neg_%s_batch_(\\d+)\\.rds", pass_tag), "\\1", cp_files)
  ))
  last_cp <- max(cp_nums)
  message(sprintf("  [%s] Checkpoint found — resuming from batch %d", pass_tag, last_cp))
  return(last_cp)
}

# Runs all batches for one reduction-level pass.
# pass_tag: checkpoint filename suffix ("base", "red10", "red20", ...)
# ade_pass: dataset for this pass (full or reduced copy of ade_raw_dt)
# red_pct: stored in the reduction_pct output column
# Owns and recycles its cluster every cluster_recycle_interval_neg batches
run_neg_batch_pass <- function(pass_tag, ade_pass, red_pct) {

  cl <- start_neg_cluster(ade_pass)
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)

  # Re-runs from the last saved checkpoint: earlier batches are already on disk,
  start_batch <- detect_neg_checkpoints(output_dir, pass_tag)

  for (batch in start_batch:n_batches) {

    # Recycles the worker pool periodically (PSOCK RSS never shrinks on its own).
    # Counting from start_batch keeps the cadence correct on a resumed run.
    if (batch > start_batch &&
        (batch - start_batch) %% cluster_recycle_interval_neg == 0) {
      stopCluster(cl)
      gc(verbose = FALSE)
      cl <- start_neg_cluster(ade_pass)
    }

    start_idx <- (batch - 1) * batch_size_neg + 1
    end_idx <- min(batch * batch_size_neg, nrow(selected_negatives))
    batch_indices <- start_idx:end_idx

    message(sprintf("\n[%s] Batch %d / %d  (triplets %d-%d)",
                    pass_tag, batch, n_batches, start_idx, end_idx))
    # Reproducibility is provided by the explicit per-triplet set.seed below
    batch_results <- foreach(
      idx = batch_indices,
      .packages = c("data.table", "mgcv"),
      .errorhandling = "pass",
      .verbose = FALSE
    ) %dopar% {

      set.seed(7113 + idx)

      rowt <- selected_negatives[idx]
      rowt$type <- "negative"

      # Builds the report-level table once
      eval_dt <- build_eval_table(ade_pass, rowt$drugA, rowt$drugB, rowt$meddra,
                                  integer(0), include_sex)
      counts <- list(
        n_events = sum(eval_dt$ea_ocurrio),
        n_coadmin = sum(eval_dt$droga_ab),
        n_events_coadmin = sum(eval_dt$droga_ab == 1L & eval_dt$ea_ocurrio == 1L)
      )

      iter_result <- tryCatch({

        model_res <- fit_gam(
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

        classic_res <- calculate_classic_ior(
          rowt$drugA, rowt$drugB, rowt$meddra,
          ade_data = NULL, eval_dt = eval_dt
        )

        classic_ac <- calculate_classic_ac(
          rowt$drugA, rowt$drugB, rowt$meddra,
          ade_data = NULL, eval_dt = eval_dt
        )

        if (!model_res$success) {
          # GAM failed: returns NA placeholders for all model-derived columns
          data.table(
            triplet_id = rowt$triplet_id,
            drugA = rowt$drugA,
            drugB = rowt$drugB,
            meddra = rowt$meddra,
            type = "negative",
            reduction_pct = red_pct,
            N = counts$n_events_coadmin,
            model_success = FALSE,
            n_events = counts$n_events,
            n_coadmin = counts$n_coadmin,
            n_stages_significant = NA_integer_,
            max_ior = NA_real_,
            mean_ior = NA_real_,
            model_aic = NA_real_,
            stage = list(1:7),
            log_ior = list(rep(NA_real_, 7)),
            log_ior_lower90 = list(rep(NA_real_, 7)),
            ior_values = list(rep(NA_real_, 7)),
            formula_used = if (!is.null(model_res$formula_attempted)) model_res$formula_attempted else NA_character_,
            message = if (!is.null(model_res$error_msg)) model_res$error_msg else NA_character_,
            classic_success = classic_res$success,
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
            AC_classic_se = list(rep(NA_real_, 7))
          )
        } else {
          # GAM succeeded: stores all per-stage estimates
          data.table(
            triplet_id = rowt$triplet_id,
            drugA = rowt$drugA,
            drugB = rowt$drugB,
            meddra = rowt$meddra,
            type = "negative",
            reduction_pct = red_pct,
            N = model_res$n_events_coadmin,
            model_success = TRUE,
            n_coadmin  = model_res$n_coadmin,
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
            log_ior_classic = if (classic_res$success) list(classic_res$results_by_stage$log_ior_classic) else list(rep(NA_real_, 7)),
            log_ior_classic_lower90 = if (classic_res$success) list(classic_res$results_by_stage$log_ior_classic_lower90) else list(rep(NA_real_, 7)),
            ior_classic = if (classic_res$success) list(classic_res$results_by_stage$ior_classic) else list(rep(NA_real_, 7)),
            ac_classic_success = classic_ac$success,
            ac_values = list(model_res$ac_values),
            ac_lower90 = list(model_res$ac_lower90),
            ac_upper90 = list(model_res$ac_upper90),
            n_stages_ac_significant = model_res$n_stages_ac_significant,
            AC_classic = if (classic_ac$success) list(classic_ac$results_by_stage$AC_classic) else list(rep(NA_real_, 7)),
            AC_classic_lower90 = if (classic_ac$success) list(classic_ac$results_by_stage$AC_classic_lower90) else list(rep(NA_real_, 7)),
            AC_classic_upper90 = if (classic_ac$success) list(classic_ac$results_by_stage$AC_classic_upper90) else list(rep(NA_real_, 7)),
            AC_classic_se = if (classic_ac$success) list(classic_ac$results_by_stage$AC_classic_se) else list(rep(NA_real_, 7))
          )
        }

      }, error = function(e) {
        data.table(
          triplet_id = rowt$triplet_id,
          drugA = rowt$drugA,
          drugB = rowt$drugB,
          meddra = rowt$meddra,
          type = "negative",
          reduction_pct = red_pct,
          N = counts$n_events_coadmin,
          model_success = FALSE,
          n_events = counts$n_events,
          n_coadmin = counts$n_coadmin,
          error_msg = paste("Unhandled error:", e$message)
        )
      })

      # prevents the progressive growth that crashes the negative passes after many batches.
      rm(eval_dt)
      gc(verbose = FALSE)

      iter_result

    } # end foreach

    # Drops foreach-level condition objects (per-triplet errors are already caught by the inner tryCatch)
    batch_results_clean <- Filter(function(x) !inherits(x, "error"), batch_results)

    if (length(batch_results_clean) > 0) {
      batch_dt <- rbindlist(batch_results_clean, fill = TRUE)

      # Writes to disk immediately; this is the only in-memory copy of the batch
      checkpoint_file <- sprintf(
        "%scheckpoint_neg_%s_batch_%d.rds", output_dir, pass_tag, batch
      )
      saveRDS(batch_dt, checkpoint_file)

      message(sprintf("  Checkpoint saved: %s  (%d/%d successful)",
                      basename(checkpoint_file),
                      sum(batch_dt$model_success, na.rm = TRUE),
                      nrow(batch_dt)))
      rm(batch_dt)
    }

    rm(batch_results, batch_results_clean)
    gc(verbose = FALSE)

  } # ends batch loop
}

# run_neg_batch_pass owns its cluster; no cluster is set up here.
message(sprintf("\nPass 1 / %d: base dataset (reduction = 0%%) ",
                length(reduction_levels) + 1))

gc()

run_neg_batch_pass(pass_tag = "base", ade_pass = ade_raw_dt, red_pct = 0)

gc()

# Pass 2+: one iteration per reduction level 
# Each builds a single reduced copy of the dataset, runs all batches, then frees the copy before the next level
for (red_pct in reduction_levels) {

  pass_tag <- sprintf("red%d", red_pct)
  n_pass <- which(reduction_levels == red_pct) + 1

  message(sprintf("\nPass %d / %d: reduced dataset (%d%% reduction)",
                  n_pass, length(reduction_levels) + 1, red_pct))

  ade_reduced <- reduce_dataset_by_stage(ade_raw_dt, red_pct, seed = 7113)

  run_neg_batch_pass(pass_tag = pass_tag, ade_pass = ade_reduced, red_pct = red_pct)

  rm(ade_reduced)
  gc()
}

################################################################################
# Assemble final results from checkpoints and save per reduction level
################################################################################
# All pass results live on disk as per-batch checkpoint files

all_pass_tags <- c("base", sprintf("red%d", reduction_levels))

negatives_scores <- rbindlist(
  lapply(all_pass_tags, function(tag) {
    cp_files <- sort(list.files(
      output_dir,
      pattern = sprintf("checkpoint_neg_%s_batch_\\d+\\.rds", tag),
      full.names = TRUE
    ))
    if (length(cp_files) == 0) {
      warning(sprintf("No checkpoint files found for pass '%s'", tag))
      return(NULL)
    }
    rbindlist(lapply(cp_files, readRDS), fill = TRUE)
  }),
  fill = TRUE
)

message(sprintf("Total negative rows assembled: %s",
                format(nrow(negatives_scores), big.mark = ",")))
message(sprintf("Successful (base): %d",
                sum(negatives_scores$model_success & negatives_scores$reduction_pct == 0,
                    na.rm = TRUE)))

for (red_pct in c(0, reduction_levels)) {
  suffix_file <- if (red_pct == 0) "" else paste0("_", red_pct)

  subset_data <- negatives_scores[reduction_pct == red_pct]

  if (nrow(subset_data) > 0) {

    saveRDS(subset_data,
            paste0(output_dir, "negative_triplets_results", suffix_file, ".rds"))

    # Flattens list columns to comma-separated strings for CSV compatibility
    subset_csv <- copy(subset_data)
    list_cols <- names(subset_csv)[sapply(subset_csv, is.list)]

    for (col in list_cols) {
      subset_csv[, (col) := sapply(get(col), function(x) {
        if (is.null(x) || length(x) == 0) return(NA_character_)
        paste(x, collapse = ",")
      })]
    }

    fwrite(subset_csv,
           paste0(output_dir, "negative_triplets_results", suffix_file, ".csv"))
    rm(subset_csv)
  }
}

################################################################################
# Null pool creation
################################################################################

all_reports <- unique(ade_raw_dt$safetyreportid)

selected_null_reports <- sample(all_reports, n_null_reports)

null_pool_meta <- unique(ade_raw_dt[
  safetyreportid %in% selected_null_reports,
  .(safetyreportid, nichd, nichd_num)
])

null_stage_distribution <- null_pool_meta[, .N, by = nichd][order(nichd)]
message("\nReport distribution by NICHD stage:")
print(null_stage_distribution)

fwrite(null_stage_distribution, paste0(output_dir, "null_pool_reports_by_stage.csv"))

fwrite(null_pool_meta, paste0(output_dir, "null_pool_reports_metadata.csv"))

message("Total reports: ", nrow(null_pool_meta))

################################################################################
# Positive vs. Negative comparison
################################################################################

pos_success <- positives_scores[model_success == TRUE & injection_success == TRUE & reduction_pct == 0]
neg_success <- negatives_scores[model_success == TRUE & reduction_pct == 0]

message("\nSuccessful triplets (base):")
message(" Positives: ", nrow(pos_success))
message(" Negatives: ", nrow(neg_success))

if (nrow(pos_success) > 0 && nrow(neg_success) > 0) {
  
  pos_expanded <- pos_success[, {
    stages <- unlist(stage)
    log_iors <- unlist(log_ior)
    log_ior_l90s <- unlist(log_ior_lower90)
    
    n <- min(length(stages), length(log_iors), length(log_ior_l90s))
    if (n == 0) {
      data.table()
    } else {
      data.table(
        stage = stages[1:n],
        log_ior = log_iors[1:n],
        log_ior_lower90 = log_ior_l90s[1:n]
      )
    }
  }, by = .(triplet_id, type, dynamic)]
  
  neg_expanded <- neg_success[, {
    stages <- unlist(stage)
    log_iors <- unlist(log_ior)
    log_ior_l90s <- unlist(log_ior_lower90)
    
    n <- min(length(stages), length(log_iors), length(log_ior_l90s))
    if (n == 0) {
      data.table()
    } else {
      data.table(
        stage = stages[1:n],
        log_ior = log_iors[1:n],
        log_ior_lower90 = log_ior_l90s[1:n]
      )
    }
  }, by = .(triplet_id, type)]
  
  message("\nlog-IOR statistics (base):")
  message("  Positives. Mean: ", round(mean(pos_expanded$log_ior, na.rm = TRUE), 4))
  message("  Negatives. Mean: ", round(mean(neg_expanded$log_ior, na.rm = TRUE), 4))

  # Persist the pos-vs-neg log-IOR contrast (means, medians, Wilcoxon) 
  logior_contrast <- data.table(
    set = c("positive", "negative"),
    n_stage_rows = c(nrow(pos_expanded), nrow(neg_expanded)),
    mean_log_ior = c(mean(pos_expanded$log_ior, na.rm = TRUE),
                     mean(neg_expanded$log_ior, na.rm = TRUE)),
    median_log_ior = c(median(pos_expanded$log_ior, na.rm = TRUE),
                       median(neg_expanded$log_ior, na.rm = TRUE)),
    sd_log_ior = c(sd(pos_expanded$log_ior, na.rm = TRUE),
                   sd(neg_expanded$log_ior, na.rm = TRUE))
  )

  if (nrow(pos_expanded) > 30 && nrow(neg_expanded) > 30) {
    test_result <- wilcox.test(
      pos_expanded$log_ior,
      neg_expanded$log_ior,
      alternative = "greater"
    )
    message("\nWilcoxon test (positives > negatives):")
    message("  p-value: ", format.pval(test_result$p.value, digits = 3))
    logior_contrast[, wilcoxon_greater_pvalue := test_result$p.value]
  }

  fwrite(logior_contrast, paste0(output_dir, "positive_vs_negative_logior_contrast.csv"))

  all_expanded <- rbind(
    pos_expanded[, .(triplet_id, type, log_ior, stage)],
    neg_expanded[, .(triplet_id, type, log_ior, stage)]
  )
  
  p1 <- ggplot(all_expanded, aes(x = type, y = log_ior, fill = type)) +
    geom_violin(alpha = 0.6) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = "Log-IOR distribution by triplet type (Base)",
      x = "Type",
      y = "Log-IOR"
    )
  
  ggsave(paste0(output_dir, "comparison_log_ior_distribution.png"),
         p1, width = 8, height = 6, dpi = 300)
  
}

################################################################################
# Final summary
################################################################################

summary_stats <- data.table(
  metric = c(
    "config_min_reports_triplet",
    "config_min_nichd_with_rep",
    "config_all_nichd_rep",
    "reduction_levels_tested",
    "n_positive_total",
    "n_positive_injected",
    "n_positive_modeled_base",
    "n_negative_total",
    "n_negative_modeled_base",
    "mean_ior_positive_base",
    "mean_ior_negative_base",
    "mean_stages_sig_positive_base",
    "mean_stages_sig_negative_base"
  ),
  value = c(
    min_reports_triplet,
    min_nichd_with_rep,
    as.numeric(all_nichd_rep),
    paste(reduction_levels, collapse = ","),
    nrow(pos_meta),
    sum(positives_scores$injection_success & positives_scores$reduction_pct == 0, na.rm = TRUE),
    sum(positives_scores$model_success & positives_scores$reduction_pct == 0, na.rm = TRUE),
    nrow(selected_negatives),
    sum(negatives_scores$model_success & negatives_scores$reduction_pct == 0, na.rm = TRUE),
    mean(positives_scores[reduction_pct == 0]$mean_ior, na.rm = TRUE),
    mean(negatives_scores[reduction_pct == 0]$mean_ior, na.rm = TRUE),
    mean(positives_scores[reduction_pct == 0]$n_stages_significant, na.rm = TRUE),
    mean(negatives_scores[reduction_pct == 0]$n_stages_significant, na.rm = TRUE)
  )
)

print(summary_stats)
fwrite(summary_stats, paste0(output_dir, "summary_statistics.csv"))

################################################################################
# Aggregated sensitivity analysis
################################################################################

sensitivity_summary <- rbind(
  positives_scores[, .(
    type = "positive",
    reduction_pct = reduction_pct[1],
    n_total = .N,
    n_success = sum(model_success, na.rm = TRUE),
    mean_ior = mean(mean_ior, na.rm = TRUE),
    mean_stages_sig = mean(n_stages_significant, na.rm = TRUE)
  ), by = reduction_pct],
  negatives_scores[, .(
    type = "negative",
    reduction_pct = reduction_pct[1],
    n_total = .N,
    n_success = sum(model_success, na.rm = TRUE),
    mean_ior = mean(mean_ior, na.rm = TRUE),
    mean_stages_sig = mean(n_stages_significant, na.rm = TRUE)
  ), by = reduction_pct]
)

print(sensitivity_summary)
fwrite(sensitivity_summary, paste0(output_dir, "sensitivity_analysis_summary.csv"))

################################################################################
# Dynamics detection
################################################################################

# Reloads the base-level positive results
ruta_pos_results <- paste0("./results/", suffix, "/augmentation_results/positive_triplets_results.rds")
positives_scores <- readRDS(ruta_pos_results)

# Filters to successfully processed positive triplets
pos_for_dynamics <- positives_scores[
  model_success == TRUE & 
  injection_success == TRUE &
  reduction_pct == 0 &   
  !is.na(dynamic)
]

# Expands per-stage list columns into one row per triplet-stage
pos_dynamics_expanded <- pos_for_dynamics[, {
  stages <- unlist(stage)
  log_iors <- unlist(log_ior)

  # Guards against mismatched list lengths before indexing
  n <- min(length(stages), length(log_iors))
  
  if (n > 0) {
    data.table(
      stage = stages[1:n],
      log_ior = log_iors[1:n]
    )
  } else {
    data.table()
  }
}, by = .(triplet_id, dynamic, fold_change)]

pos_dynamics_expanded[, stage_name := niveles_nichd[stage]]

# Mean log-IOR by injected dynamic and developmental stage
dynamics_summary <- pos_dynamics_expanded[, .(
  mean_log_ior = mean(log_ior, na.rm = TRUE),
  sd_log_ior = sd(log_ior, na.rm = TRUE),
  n_triplets = uniqueN(triplet_id)
), by = .(dynamic, stage)]

dynamics_summary[, stage_name := niveles_nichd[stage]]

# Uses the uniform dynamic as the baseline; other dynamics are contrasted against it
uniform_baseline <- dynamics_summary[dynamic == "uniform", .(
  stage,
  baseline_log_ior = mean_log_ior
)]

dynamics_diff <- merge(
  dynamics_summary[dynamic != "uniform"],
  uniform_baseline,
  by = "stage",
  all.x = TRUE
)

dynamics_diff[, log_ior_diff := mean_log_ior - baseline_log_ior]

print(dynamics_diff[order(dynamic, stage), .(
  dynamic, 
  stage_name, 
  mean_log_ior, 
  baseline = baseline_log_ior,
  difference = log_ior_diff
)])

################################################################################
# Bootstrap confidence intervals
################################################################################

n_boot <- 100

dynamics_nonuniform <- unique(pos_dynamics_expanded[dynamic != "uniform", dynamic])
stages <- 1:7

bootstrap_results <- rbindlist(pblapply(dynamics_nonuniform, function(dyn) {
  rbindlist(lapply(stages, function(s) {
    boot_res <- bootstrap_dynamic_diff(pos_dynamics_expanded, dyn, s, n_boot)
    cbind(data.table(dynamic = dyn, stage = s), boot_res)
  }))
}))

dynamics_with_ci <- merge(
  dynamics_diff,
  bootstrap_results,
  by = c("dynamic", "stage"),
  all.x = TRUE
)

dynamics_with_ci[, stage_name := factor(stage_name, levels = niveles_nichd)]

fwrite(dynamics_with_ci, paste0(output_dir, "dynamics_recovery_analysis.csv"))

recovery_stats <- dynamics_with_ci[, .(
  mean_difference = mean(log_ior_diff, na.rm = TRUE),
  max_difference = max(abs(log_ior_diff), na.rm = TRUE),
  stages_significant = sum(ci_lower > 0 | ci_upper < 0, na.rm = TRUE)
), by = dynamic]

print(recovery_stats)

################################################################################
# Visualization
################################################################################

dynamic_colors <- c(
  "increase" = "#E41A1C",
  "decrease" = "#377EB8", 
  "plateau" = "#4DAF4A",
  "inverse_plateau" = "#984EA3"
)

nichd_labels <- c(
  "term_neonatal" = "Term neonatal",
  "infancy" = "Infancy",
  "toddler" = "Toddler",
  "early_childhood" = "Early childhood",
  "middle_childhood" = "Middle childhood",
  "early_adolescence" = "Early adolescence",
  "late_adolescence" = "Late adolescence"
)

p_dynamics_diff <- ggplot(
  dynamics_with_ci[!is.na(mean_diff)],
  aes(x = stage_name, y = log_ior_diff, color = dynamic, group = dynamic)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper, fill = dynamic),
    alpha = 0.2,
    color = NA
  ) +
  scale_color_manual(
    values = dynamic_colors,
    labels = c(
      "increase" = "Increase",
      "decrease" = "Decrease",
      "plateau" = "Plateau",
      "inverse_plateau" = "Inverse plateau"
    )
  ) +
  scale_fill_manual(
    values = dynamic_colors,
    labels = c(
      "increase" = "Increase",
      "decrease" = "Decrease",
      "plateau" = "Plateau",
      "inverse_plateau" = "Inverse plateau"
    )
  ) +
  scale_x_discrete(labels = nichd_labels) +
  labs(
    x = "Stage",
    y = "Δ Log(IOR) (vs uniform)",
    color = "Injected dynamic",
    fill = "Injected dynamic"
  ) 

ggsave(
  paste0(output_dir, "fig_dynamics_recovery.png"),
  p_dynamics_diff,
  width = 12,
  height = 8,
  dpi = 300
)

print(p_dynamics_diff)

################################################################################
# Dynamics detection using AC
################################################################################

# Mirrors the log-IOR dynamics block above using GAM-derived AC

# Expands per-stage AC list columns into one row per triplet-stage
pos_dynamics_expanded_ac <- pos_for_dynamics[, {
  stages <- unlist(stage)
  ac_vals <- unlist(ac_values)

  # Guards against mismatched list lengths
  n <- min(length(stages), length(ac_vals))
  
  if (n > 0) {
    data.table(
      stage = stages[1:n],
      ac = ac_vals[1:n]
    )
  } else {
    data.table()
  }
}, by = .(triplet_id, dynamic, fold_change)]

pos_dynamics_expanded_ac[, stage_name := niveles_nichd[stage]]

# Uses median/MAD rather than mean/SD because AC is right-skewed across triplets
dynamics_summary_ac <- pos_dynamics_expanded_ac[, .(
  mean_ac = median(ac, na.rm = TRUE),
  sd_ac = mad(ac, na.rm = TRUE),
  n_triplets = uniqueN(triplet_id)
), by = .(dynamic, stage)]

dynamics_summary_ac[, stage_name := niveles_nichd[stage]]

# Uniform dynamic as AC reference baseline
uniform_baseline_ac <- dynamics_summary_ac[dynamic == "uniform", .(
  stage,
  baseline_ac = mean_ac   # already median; name kept for join consistency
)]

# Delta AC vs. uniform baseline
dynamics_diff_ac <- merge(
  dynamics_summary_ac[dynamic != "uniform"],
  uniform_baseline_ac,
  by = "stage",
  all.x = TRUE
)

dynamics_diff_ac[, ac_diff := mean_ac - baseline_ac]

print(dynamics_diff_ac[order(dynamic, stage), .(
  dynamic, 
  stage_name, 
  mean_ac, 
  baseline = baseline_ac,
  difference = ac_diff
)])

################################################################################
# Bootstrap confidence intervals for AC
################################################################################

dynamics_nonuniform_ac <- unique(pos_dynamics_expanded_ac[dynamic != "uniform", dynamic])
stages_ac <- 1:7

bootstrap_results_ac <- rbindlist(pblapply(dynamics_nonuniform_ac, function(dyn) {
  rbindlist(lapply(stages_ac, function(s) {
    boot_res <- bootstrap_dynamic_diff_ac(pos_dynamics_expanded_ac, dyn, s, n_boot)
    cbind(data.table(dynamic = dyn, stage = s), boot_res)
  }))
}))

dynamics_with_ci_ac <- merge(
  dynamics_diff_ac,
  bootstrap_results_ac,
  by = c("dynamic", "stage"),
  all.x = TRUE
)

dynamics_with_ci_ac[, stage_name := factor(stage_name, levels = niveles_nichd)]

fwrite(dynamics_with_ci_ac, paste0(output_dir, "dynamics_recovery_analysis_ac.csv"))

recovery_stats_ac <- dynamics_with_ci_ac[, .(
  mean_difference = mean(ac_diff, na.rm = TRUE),
  max_difference = max(abs(ac_diff), na.rm = TRUE),
  stages_significant = sum(ci_lower > 0 | ci_upper < 0, na.rm = TRUE)
), by = dynamic]

print(recovery_stats_ac)

################################################################################
# Delta AC plot (vs. uniform dynamic)
################################################################################

p_dynamics_diff_ac <- ggplot(
  dynamics_with_ci_ac[!is.na(mean_diff)],
  aes(x = stage_name, y = ac_diff, color = dynamic, group = dynamic)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper, fill = dynamic),
    alpha = 0.2,
    color = NA
  ) +
  scale_color_manual(
    values = dynamic_colors,
    labels = c(
      "increase" = "Increase",
      "decrease" = "Decrease",
      "plateau" = "Plateau",
      "inverse_plateau" = "Inverse plateau"
    )
  ) +
  scale_x_discrete(labels = nichd_labels) +
  scale_fill_manual(
    values = dynamic_colors,
    labels = c(
      "increase" = "Increase",
      "decrease" = "Decrease",
      "plateau" = "Plateau",
      "inverse_plateau" = "Inverse plateau"
    )
  ) +
  labs(
    x = "Stage",
    y = "Δ AC (vs uniform)",
    color = "Injected dynamic",
    fill = "Injected dynamic"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )
  
ggsave(
  paste0(output_dir, "fig_dynamics_recovery_ac.png"),
  p_dynamics_diff_ac,
  width = 12,
  height = 8,
  dpi = 300
)

print(p_dynamics_diff_ac)
