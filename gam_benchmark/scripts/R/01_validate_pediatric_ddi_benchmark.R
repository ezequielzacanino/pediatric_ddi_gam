################################################################################
# GAM benchmark script
# Script 01_validate_pediatric_ddi_benchmark
################################################################################

source("00_functions.R", local = TRUE)

################################################################################
# Configuration
################################################################################

# Bootstrap iterations for sensitivity confidence intervals.
n_boot <- 1000

# Set TRUE to refit even when a cached RDS already exists.
overwrite_cache <- TRUE

# MedDRA roll-up level shared by ade_raw and the curated benchmark.
meddra_rollup_level <- "HLT"

# Curated set produced by the upstream ddi_reference_set project (script 01).
curated_set_dir <- "../ddi_reference_set/results/curated_pediatric_ddi_reference_set"
benchmark_triplets_path <- file.path(curated_set_dir, "curated_pediatric_ddi_triplets.csv")
benchmark_sources_path <- file.path(curated_set_dir, "curated_pediatric_ddi_sources.csv")

output_dir <- file.path("./results", suffix, "benchmark_validation")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

benchmark_ready_path <- file.path(output_dir, "benchmark_triplets_ready_for_modeling.csv")
benchmark_summary_path <- file.path(output_dir, "benchmark_summary.csv")
benchmark_unmapped_path <- file.path(output_dir, "benchmark_triplets_unmapped.csv")
benchmark_fit_rds_path <- file.path(
  output_dir,
  paste0("benchmark_fit_results_", tolower(meddra_rollup_level), ".rds")
)
benchmark_fit_status_path <- file.path(output_dir, "benchmark_fit_status.csv")
benchmark_expanded_path <- file.path(output_dir, "benchmark_stage_level_results.csv")
benchmark_triplet_detail_path <- file.path(output_dir, "benchmark_triplet_detection_detail.csv")
benchmark_metrics_path <- file.path(output_dir, "benchmark_metrics.csv")
benchmark_ontogeny_contrast_path <- file.path(output_dir, "benchmark_ontogeny_stage_contrast.csv")
benchmark_ontogeny_summary_path <- file.path(output_dir, "benchmark_ontogeny_summary.csv")

methods_dt <- data.table(
  method = c(
    "GAM-logIOR_nom", "GAM-logIOR", "GAM-AC_nom", "GAM-AC",
    "Estratificado-IOR", "Estratificado-IOR_null",
    "Estratificado-AC", "Estratificado-AC_null"
  ),
  detection_type = c("IOR", "IOR", "AC", "AC", "IOR", "IOR", "AC", "AC"),
  use_null = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)
)
methods_dt[, threshold_mode := fifelse(use_null, paste0("null_", percentil), "nominal")]

################################################################################
# Load data
################################################################################

# Translate benchmark drug/event terms to the same ATC and MedDRA identifiers used by ade_raw.
benchmark_ready <- prepare_benchmark_reference_set(
  benchmark_triplets_path = benchmark_triplets_path,
  benchmark_sources_path = benchmark_sources_path,
  output_path = benchmark_ready_path,
  rollup_level = meddra_rollup_level
)

benchmark_ready[, pair_key_original := paste(
  pmin(drug1_name, drug2_name),
  pmax(drug1_name, drug2_name),
  sep = " | "
)]

ade_data <- load_ade_modeling_data(rollup_level = meddra_rollup_level)
null_thresholds <- load_modeling_null_thresholds()

################################################################################
# Mapping diagnosis
################################################################################

# Separate PT/HLT/HLGT coverage counts so mapping failures are visible before fitting.
benchmark_summary <- data.table(
  n_raw_benchmark_triplets = nrow(benchmark_ready),
  n_raw_unique_pairs = benchmark_ready[, uniqueN(pair_key_original)],
  n_mapped_PT = benchmark_ready[, sum(mapped_PT, na.rm = TRUE)],
  coverage_PT = benchmark_ready[, mean(mapped_PT, na.rm = TRUE)],
  n_mapped_HLT = benchmark_ready[, sum(mapped_HLT, na.rm = TRUE)],
  coverage_HLT = benchmark_ready[, mean(mapped_HLT, na.rm = TRUE)],
  n_mapped_HLGT = benchmark_ready[, sum(mapped_HLGT, na.rm = TRUE)],
  coverage_HLGT = benchmark_ready[, mean(mapped_HLGT, na.rm = TRUE)],
  n_mapped_benchmark_triplets = benchmark_ready[, sum(mapping_success, na.rm = TRUE)],
  n_unmapped_benchmark_triplets = benchmark_ready[, sum(!mapping_success, na.rm = TRUE)],
  n_mapped_unique_pairs = benchmark_ready[mapping_success == TRUE, uniqueN(paste(drugA, drugB, sep = "|"))],
  meddra_rollup_level_used = meddra_rollup_level
)

fwrite(benchmark_summary, benchmark_summary_path)
fwrite(benchmark_ready[mapping_success == FALSE], benchmark_unmapped_path)

cat(sprintf("mapping summary [active roll-up: %s]\n", meddra_rollup_level))
cat(sprintf("- PT mapping:   %d/%d (%.1f%%)\n",
  benchmark_summary$n_mapped_PT,
  benchmark_summary$n_raw_benchmark_triplets,
  benchmark_summary$coverage_PT * 100
))
cat(sprintf("- HLT mapping:  %d/%d (%.1f%%)%s\n",
  benchmark_summary$n_mapped_HLT,
  benchmark_summary$n_raw_benchmark_triplets,
  benchmark_summary$coverage_HLT * 100,
  if (meddra_rollup_level == "HLT") " <- ACTIVE" else ""
))
cat(sprintf("- HLGT mapping: %d/%d (%.1f%%)%s\n",
  benchmark_summary$n_mapped_HLGT,
  benchmark_summary$n_raw_benchmark_triplets,
  benchmark_summary$coverage_HLGT * 100,
  if (meddra_rollup_level == "HLGT") " <- ACTIVE" else ""
))

################################################################################
# Benchmark fitting
################################################################################

# Fit each mapped triplet against the same ade_raw representation used by the GAM and stratified estimators.
benchmark_fit <- fit_benchmark_triplets(
  benchmark_ready_dt = benchmark_ready,
  ade_data = ade_data,
  cache_file = benchmark_fit_rds_path,
  overwrite_cache = overwrite_cache
)

benchmark_fit_status <- benchmark_fit[, .(
  mapping_success = unique(mapping_success)[1L],
  model_success = any(model_success, na.rm = TRUE),
  classic_success = any(classic_success, na.rm = TRUE)
), by = triplet_id]
fwrite(benchmark_fit_status, benchmark_fit_status_path)

benchmark_expanded <- expand_benchmark_results(benchmark_fit, null_thresholds)
fwrite(benchmark_expanded, benchmark_expanded_path)

################################################################################
# Detection
################################################################################

methods_cfg <- lapply(seq_len(nrow(methods_dt)), function(i) {
  list(
    method = methods_dt$method[i],
    detection = methods_dt$detection_type[i],
    use_null = methods_dt$use_null[i]
  )
})

benchmark_eval <- evaluate_benchmark_methods(
  benchmark_expanded_dt = benchmark_expanded,
  n_boot = n_boot,
  methods_cfg = methods_cfg
)

# For triplets with a curated expected high-risk window, compare the expected NICHD stages against the stages each method detects.
benchmark_ontogeny_contrast <- build_ontogeny_stage_contrast(benchmark_expanded, methods_cfg)
benchmark_ontogeny_summary <- summarize_ontogeny_contrast(benchmark_ontogeny_contrast)
fwrite(benchmark_ontogeny_contrast, benchmark_ontogeny_contrast_path)
fwrite(benchmark_ontogeny_summary, benchmark_ontogeny_summary_path)

triplet_method_grid <- merge(
  CJ(triplet_id = unique(benchmark_ready$triplet_id), method = methods_dt$method, unique = TRUE),
  methods_dt[, .(method, detection_type, threshold_mode, use_null)],
  by = "method",
  all.x = TRUE
)

benchmark_triplet_detail <- merge(
  triplet_method_grid,
  benchmark_eval$triplet_detail,
  by = c("triplet_id", "method", "detection_type", "threshold_mode"),
  all.x = TRUE
)
benchmark_triplet_detail[is.na(evaluable), evaluable := FALSE]
benchmark_triplet_detail[is.na(detected), detected := FALSE]
benchmark_triplet_detail[is.na(n_stages_detected), n_stages_detected := 0L]

meta_cols_out <- intersect(
  c(
    "triplet_id", "control_type", "drug1_name", "drug2_name", "meddra_pt",
    "interaction_type", "mechanism", "evidence_type", "evidence_level",
    "is_ime", "is_dme", "confidence_level",
    "ontogenic_modulation", "higher_risk_stages",
    "mapping_success", "drugA", "drugB", "meddra"
  ),
  names(benchmark_ready)
)
benchmark_triplet_detail <- merge(
  benchmark_triplet_detail,
  benchmark_ready[, ..meta_cols_out],
  by = "triplet_id",
  all.x = TRUE
)
benchmark_triplet_detail <- merge(
  benchmark_triplet_detail,
  benchmark_fit_status,
  by = c("triplet_id", "mapping_success"),
  all.x = TRUE
)
benchmark_triplet_detail[is.na(model_success), model_success := FALSE]
benchmark_triplet_detail[is.na(classic_success), classic_success := FALSE]
fwrite(benchmark_triplet_detail, benchmark_triplet_detail_path)

################################################################################
# Results
################################################################################

benchmark_metrics_list <- vector("list", nrow(methods_dt))
for (i in seq_len(nrow(methods_dt))) {
  method_row <- methods_dt[i]
  detail_i <- benchmark_triplet_detail[
    method == method_row$method,
    .(triplet_id, control_type, evaluable, detected, n_stages_detected, signal_score)
  ]
  benchmark_metrics_list[[i]] <- calculate_benchmark_metrics_from_triplet_detail(
    triplet_detail_dt = detail_i,
    method_name = method_row$method,
    detection_type = method_row$detection_type,
    use_null = method_row$use_null,
    n_boot = n_boot
  )
}

benchmark_metrics <- rbindlist(benchmark_metrics_list, fill = TRUE)
# Recovery counts are a positive-control summary.
benchmark_detection_summary <- benchmark_triplet_detail[control_type == "positive", .(
  n_mapped_triplets = sum(mapping_success, na.rm = TRUE),
  n_detected_mapped_triplets = sum(detected & mapping_success, na.rm = TRUE),
  detection_over_raw_benchmark = mean(detected, na.rm = TRUE),
  detection_over_mapped_benchmark = fifelse(
    sum(mapping_success, na.rm = TRUE) > 0,
    sum(detected & mapping_success, na.rm = TRUE) / sum(mapping_success, na.rm = TRUE),
    NA_real_
  )
), by = .(method, detection_type, threshold_mode)]

benchmark_metrics <- merge(
  benchmark_metrics,
  benchmark_detection_summary,
  by = c("method", "detection_type", "threshold_mode"),
  all.x = TRUE
)
benchmark_metrics[, names(benchmark_summary) := as.list(benchmark_summary)]
setorder(benchmark_metrics, detection_type, threshold_mode, method)
fwrite(benchmark_metrics, benchmark_metrics_path)

cat(sprintf("Benchmark metrics saved to %s\n", benchmark_metrics_path))
