################################################################################
# Plot generation script for the metrics produced by 30_metrics.R
# Script: 40_graphs.R
################################################################################

source("00_functions.R", local = TRUE)

output_dir <- paste0("./results/", suffix, "/metrics_results/")
fig_output_dir <- paste0(output_dir, "facet_figures/")

# Create output directory if it does not exist
dir.create(fig_output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Preprocessing
################################################################################

# Display labels for NICHD developmental stages
nichd_labels <- c(
  "term_neonatal" = "Term neonatal",
  "infancy" = "Infancy",
  "toddler" = "Toddler",
  "early_childhood" = "Early childhood",
  "middle_childhood" = "Middle childhood",
  "early_adolescence" = "Early adolescence",
  "late_adolescence" = "Late adolescence"
)

# Metric display labels
metric_labels <- c(
  "AUC" = "AUC",
  "sensitivity" = "Sensitivity",
  "specificity" = "Specificity",
  "PPV" = "PPV",
  "NPV" = "NPV",
  "F1" = "F1-Score"
)

# Signal dynamic display labels
dynamic_labels <- c(
  "uniform" = "Uniform",
  "increase" = "Increase",
  "plateau" = "Plateau",
  "decrease" = "Decrease",
  "inverse_plateau" = "Inverse plateau"
)

# Method pairs to compare: one pair per plot page (both null used in manuscript)
# GAM (null) vs Stratified (nominal): primary comparison
# GAM (null) vs Stratified (null): apples-to-apples null threshold comparison
# GAM (nominal) vs Stratified (nominal): purely nominal comparison
method_pairs <- list(
  list(
    name = "IOR", gam = "GAM-logIOR", classic = "Estratificado-IOR",
    gam_label = "GAM-IOR (null dist.)", classic_label = "Stratified-IOR (nominal)"),
  list(
    name = "AC", gam = "GAM-AC", classic = "Estratificado-AC",
    gam_label = "GAM-AC (null dist.)", classic_label = "Stratified-AC (nominal)"),
  list(
    name = "IOR_null", gam = "GAM-logIOR", classic = "Estratificado-IOR_null",
    gam_label = "GAM-IOR", classic_label = "Stratified-IOR"),
  list(
    name = "AC_null", gam = "GAM-AC", classic = "Estratificado-AC_null",
    gam_label = "GAM-AC", classic_label = "Stratified-AC"),
  list(
    name = "IOR_nom", gam = "GAM-logIOR_nom", classic = "Estratificado-IOR",
    gam_label = "GAM-IOR (nominal)", classic_label = "Stratified-IOR (nominal)"
  ),
  list(
    name = "AC_nom", gam = "GAM-AC_nom", classic = "Estratificado-AC",
    gam_label = "GAM-AC (nominal)", classic_label = "Stratified-AC (nominal)"
 )
)

# Dataset versions to process
dataset_versions <- c("original", "filtered", "intersection")

# Helper: readable title for each dataset version
version_title_label <- function(version) {
  switch(version,
    "original" = "Original dataset",
    "filtered" = "Power-filtered dataset",
    "intersection" = "Intersection dataset",
    version
  )
}

# Build color and linetype scales for a given method pair
make_scales <- function(pair) {
  list(
    color = setNames(c("#16A085", "#C0392B"), c(pair$gam, pair$classic)),
    linetype = setNames(c("solid", "dashed"), c(pair$gam, pair$classic)),
    labels = setNames(c(pair$gam_label, pair$classic_label), c(pair$gam, pair$classic))
  )
}

################################################################################
# Dynamic pattern plot
################################################################################

# Plot of tangential weight functions by signal dynamic type
dt_dyn_plot <- rbindlist(lapply(
  c("increase", "decrease", "plateau", "inverse_plateau"),
  function(d) data.table(
    stage = 1:7,
    value = generate_dynamic(d),
    dynamic = factor(
      c(increase = "Increasing", decrease = "Decreasing",
        plateau = "Plateau", inverse_plateau = "Valley")[d],
      levels = c("Increasing", "Decreasing", "Plateau", "Valley")
    )
  )
))

dynamics <- ggplot(dt_dyn_plot, aes(x = stage, y = value, color = dynamic)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~dynamic, nrow = 2) +
  scale_color_manual(values = c(
    "Increasing" = "#16A085",
    "Decreasing" = "#C0392B",
    "Plateau" = "#2980B9",
    "Valley" = "#8E44AD"
  ), guide = "none") +
  scale_x_continuous(breaks = 1:7, labels = nichd_labels, name = NULL) +
  scale_y_continuous(name = "Relative weight in reporting probability",
                     limits = c(-1.1, 1.1), breaks = c(-1, 0, 1)) +
  labs(title = NULL) +
  coord_cartesian(clip = "off") +
  theme( axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin  = margin(t = 5, r = 5, b = 5, l = 25)
  )
ggsave(paste0(output_dir, "dynamics.png"), dynamics, width = 15, height = 15, dpi = 300)


################################################################################
# Data loading
################################################################################

metrics_global <- list(
  original = fread(paste0(output_dir, "metrics_global_original.csv")),
  filtered = fread(paste0(output_dir, "metrics_global_filtered.csv")),
  intersection = fread(paste0(output_dir, "metrics_global_intersection.csv"))
)

metrics_dynamic <- list(
  original = fread(paste0(output_dir, "metrics_dynamic_original.csv")),
  filtered = fread(paste0(output_dir, "metrics_dynamic_filtered.csv")),
  intersection = fread(paste0(output_dir, "metrics_dynamic_intersection.csv"))
)

metrics_stage <- list(
  original = fread(paste0(output_dir, "metrics_stage_original.csv")),
  filtered = fread(paste0(output_dir, "metrics_stage_filtered.csv")),
  intersection = fread(paste0(output_dir, "metrics_stage_intersection.csv"))
)

# Add method_type column and convert nichd/dynamic to factors across all versions
for (v in dataset_versions) {
  metrics_global[[v]][, method_type := ifelse(grepl("GAM", method), "GAM", "Stratified")]

  dt <- metrics_dynamic[[v]]
  if ("dinamica" %in% names(dt)) setnames(dt, "dinamica", "dynamic")
  dt[, dynamic := factor(dynamic)]
  dt[, method_type := ifelse(grepl("GAM", method), "GAM", "Stratified")]
  metrics_dynamic[[v]] <- dt

  metrics_stage[[v]][, nichd := factor(nichd, levels = niveles_nichd)]
  metrics_stage[[v]][, method_type := ifelse(grepl("GAM", method), "GAM", "Stratified")]
}

################################################################################
# Data preparation for faceted plots
################################################################################

# Reshape stage-level metrics from wide to long format for ggplot2 faceting.
#
# Filters the data to the two methods in the given pair, then pivots
# all six performance metrics into a single long-format data.table.
#
# Arguments:
#  dt : data.table with per-stage metrics in wide format
#  pair : list with method identifiers and display labels (gam, classic, etc.)
#
# Returns: long-format data.table ready for ggplot2

prepare_facet_data <- function(dt, pair) {
  dt_filtered <- dt[method %in% c(pair$gam, pair$classic)]
  metrics <- c("AUC", "sensitivity", "specificity", "PPV", "NPV", "F1")
  # Pivot to long format, one row per method x stage x metric
  dt_long <- rbindlist(lapply(metrics, function(m) {
    lower_col <- paste0(m, "_lower")
    upper_col <- paste0(m, "_upper")

    # Skip metric if confidence interval columns are missing
    if (!(lower_col %in% names(dt_filtered)) || !(upper_col %in% names(dt_filtered))) {
      return(NULL)
    }

    data.table(
      method = dt_filtered$method,
      method_type = dt_filtered$method_type,
      nichd = dt_filtered$nichd,
      reduction_pct = dt_filtered$reduction_pct,
      metric = m,
      metric_label = metric_labels[m],
      value = dt_filtered[[m]],
      lower = dt_filtered[[lower_col]],
      upper = dt_filtered[[upper_col]]
    )
  }), use.names = TRUE)
  # Assign display label for each method
  dt_long[, method_label := ifelse(
    method == pair$gam,
    pair$gam_label,
    pair$classic_label
  )]
  # Set metric factor order to match the desired panel sequence
  dt_long[, metric_label := factor(
    metric_label,
    levels = metric_labels
  )]
  return(dt_long)
}

# Reshape dynamic-level metrics from wide to long format for faceting
prepare_facet_data_dynamic <- function(dt, pair) {
  dt_filtered <- dt[method %in% c(pair$gam, pair$classic)]
  metrics <- c("AUC", "sensitivity", "specificity", "PPV", "NPV", "F1")

  dt_long <- rbindlist(lapply(metrics, function(m) {
    lower_col <- paste0(m, "_lower")
    upper_col <- paste0(m, "_upper")
    if (!(lower_col %in% names(dt_filtered)) || !(upper_col %in% names(dt_filtered))) return(NULL)
    data.table(
      method = dt_filtered$method,
      method_type = dt_filtered$method_type,
      dynamic = dt_filtered$dynamic,
      reduction_pct = dt_filtered$reduction_pct,
      metric = m,
      metric_label = metric_labels[m],
      value = dt_filtered[[m]],
      lower = dt_filtered[[lower_col]],
      upper = dt_filtered[[upper_col]]
    )
  }), use.names = TRUE)

  dt_long[, method_label := ifelse(method == pair$gam, pair$gam_label, pair$classic_label)]
  dt_long[, metric_label := factor(metric_label, levels = metric_labels)]
  return(dt_long)
}

# Reshape global (unstratified) metrics from wide to long format for faceting
prepare_facet_data_global <- function(dt, pair) {
  dt_filtered <- dt[method %in% c(pair$gam, pair$classic)]
  metrics <- c("AUC", "sensitivity", "specificity", "PPV", "NPV", "F1")

  dt_long <- rbindlist(lapply(metrics, function(m) {
    lower_col <- paste0(m, "_lower")
    upper_col <- paste0(m, "_upper")
    if (!(lower_col %in% names(dt_filtered)) || !(upper_col %in% names(dt_filtered))) return(NULL)
    data.table(
      method = dt_filtered$method,
      method_type = dt_filtered$method_type,
      reduction_pct = dt_filtered$reduction_pct,
      metric = m,
      metric_label = metric_labels[m],
      value = dt_filtered[[m]],
      lower = dt_filtered[[lower_col]],
      upper = dt_filtered[[upper_col]]
    )
  }), use.names = TRUE)

  dt_long[, method_label := ifelse(method == pair$gam, pair$gam_label, pair$classic_label)]
  dt_long[, metric_label := factor(metric_label, levels = metric_labels)]
  return(dt_long)
}

################################################################################
# Main plot generation functions
################################################################################

# Faceted plot: metrics x NICHD stage
plot_facet_metrics <- function(dt_long, pair, version) {
  if (is.null(dt_long) || nrow(dt_long) == 0) {
    message(sprintf("No data for %s - %s (stage)", pair$name, version))
    return(NULL)
  }

  sc <- make_scales(pair)

  ggplot(dt_long, aes(x = reduction_pct, y = value, color = method, group = method)) +
    geom_point(size = 2.5, position = position_dodge(width = 5)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 4, alpha = 0.8, position = position_dodge(width = 5)) +
    facet_grid(metric_label ~ nichd, labeller = labeller(nichd = nichd_labels), scales = "free_y") +
    scale_color_manual(name = "Method", values = sc$color, labels = sc$labels) +
    scale_linetype_manual(name = "Method", values = sc$linetype, labels = sc$labels) +
    scale_x_continuous(breaks = seq(0, 90, by = 10), name = "Dataset reduction (%)") +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 4),
      name = "Metric value",
      # small multiplicative padding so errorbars don't clip at panel edge
      expand = expansion(mult = 0.08)
    ) +
    labs( title = sprintf("Stage-level metrics - %s vs %s", pair$gam_label, pair$classic_label), subtitle = version_title_label(version))
}

# Faceted plot: metrics x signal dynamics
plot_facet_metrics_dynamic <- function(dt_long, pair, version) {
  if (is.null(dt_long) || nrow(dt_long) == 0) {
    message(sprintf("No data for %s - %s (dynamic)", pair$name, version))
    return(NULL)
  }

  sc <- make_scales(pair)

  ggplot(dt_long, aes(x = reduction_pct, y = value, color = method, group = method)) +
    geom_point(size = 2.5, position = position_dodge(width = 5)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 4, alpha = 0.8, position = position_dodge(width = 5)) +
    facet_grid(metric_label ~ dynamic, labeller = labeller(dynamic = dynamic_labels), scales = "free_y") +
    scale_color_manual(name = "Method", values = sc$color, labels = sc$labels) +
    scale_linetype_manual(name = "Method", values = sc$linetype, labels = sc$labels) +
    scale_x_continuous(breaks = seq(0, 90, by = 10), name = "Dataset reduction (%)") +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 4),
      name = "Metric value",
      # small multiplicative padding so errorbars don't clip at panel edge
      expand = expansion(mult = 0.08)
    ) +
    labs( title = sprintf("Metrics by dynamic - %s vs %s", pair$gam_label, pair$classic_label), subtitle = version_title_label(version))
}

# Faceted plot: global (unstratified) metrics via facet_wrap per metric
plot_facet_metrics_global <- function(dt_long, pair, version) {
  if (is.null(dt_long) || nrow(dt_long) == 0) {
    message(sprintf("No data for %s - %s (global)", pair$name, version))
    return(NULL)
  }

  sc <- make_scales(pair)

  ggplot(dt_long, aes(x = reduction_pct, y = value, color = method, group = method)) +
    geom_point(size = 2.5, position = position_dodge(width = 5)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 4, alpha = 0.8, position = position_dodge(width = 5)) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
    scale_color_manual(name = "Method", values = sc$color, labels = sc$labels) +
    scale_linetype_manual(name = "Method", values = sc$linetype, labels = sc$labels) +
    scale_x_continuous(breaks = seq(0, 90, by = 10), name = "Dataset reduction (%)") +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 4),
      name = "Metric value",
      # small multiplicative padding so errorbars don't clip at panel edge
      expand = expansion(mult = 0.08)
    ) +
    labs(
      title = sprintf("Global metrics - %s vs %s", pair$gam_label, pair$classic_label),
      subtitle = version_title_label(version)
    )
}

################################################################################
# Plot saving
################################################################################

save_plot <- function(p, file_suffix, width, height) {
  png_path <- paste0(fig_output_dir, "fig_facet_", file_suffix, ".png")
  svg_path <- paste0(fig_output_dir, "fig_facet_", file_suffix, ".svg")
  ggsave(png_path, p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(svg_path, p, width = width, height = height, device = svglite)
}

# Stage-stratified faceted plots — iterate over all versions and method pairs
generate_all_facet_plots <- function() {
  for (version in dataset_versions) {
    dt_stage <- metrics_stage[[version]]
    for (pair in method_pairs) {
      dt_long <- prepare_facet_data(dt_stage, pair)
      p <- plot_facet_metrics(dt_long, pair, version)
      if (!is.null(p)) {
        save_plot(p, sprintf("%s_%s", tolower(pair$name), version), width = 16, height = 12)
      }
    }
  }
}

# Dynamic-stratified faceted plots — iterate over all versions and method pairs
generate_all_facet_plots_dynamic <- function() {
  for (version in dataset_versions) {
    dt_dynamic <- metrics_dynamic[[version]]
    for (pair in method_pairs) {
      dt_long <- prepare_facet_data_dynamic(dt_dynamic, pair)
      p <- plot_facet_metrics_dynamic(dt_long, pair, version)
      if (!is.null(p)) {
        save_plot(p, sprintf("%s_%s_dynamic", tolower(pair$name), version), width = 16, height = 12)
      }
    }
  }
}

# Global (unstratified) faceted plots — iterate over all versions and method pairs
generate_all_facet_plots_global <- function() {
  for (version in dataset_versions) {
    dt_global <- metrics_global[[version]]
    for (pair in method_pairs) {
      dt_long <- prepare_facet_data_global(dt_global, pair)
      p <- plot_facet_metrics_global(dt_long, pair, version)
      if (!is.null(p)) {
        save_plot(p, sprintf("%s_%s_global", tolower(pair$name), version), width = 12, height = 10)
      }
    }
  }
}

# Execute all faceted plot generators
generate_all_facet_plots()
generate_all_facet_plots_dynamic()
generate_all_facet_plots_global()

################################################################################
# Global ROC curves
################################################################################

# ROC curves 
# The continuous predictor is the interaction 90% CI lower bound
roc_global <- fread(paste0(output_dir, "roc_global_original.csv"))

# Method pairs compared in each ROC panel (GAM vs stratified, per measure).
# Only the base methods are needed: the _null/_nom variants share the same continuous score 
roc_pairs <- list(
  list(name = "IOR", gam = "GAM-logIOR", classic = "Estratificado-IOR"),
  list(name = "AC", gam = "GAM-AC", classic = "Estratificado-AC")
)

# Builds a single ROC panel comparing GAM vs stratified at baseline (no reduction)
plot_roc_pair <- function(roc_dt, pair) {
  dt <- roc_dt[method %in% c(pair$gam, pair$classic) & reduction_pct == 0 & predictor == "lower90"]
  if (nrow(dt) == 0) {
    message(sprintf("No ROC data for %s", pair$name))
    return(NULL)
  }

  dt[, method_label := ifelse(method == pair$gam, "GAM", "Stratified")]

  # AUC annotation, one entry per method
  auc_lab <- unique(dt[, .(method_label, auc)])
  auc_lab[, lab := sprintf("%s: AUC = %.3f", method_label, auc)]
  setorder(auc_lab, method_label)
  auc_lab[, y_pos := 0.05 + 0.06 * (.I - 1)]

  ggplot(dt, aes(x = fpr, y = tpr, color = method_label)) +
    # Chance diagonal as reference
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray60") +
    geom_line(linewidth = 0.9) +
    geom_text(
      data = auc_lab,
      aes(x = 0.55, y = y_pos, label = lab, color = method_label),
      inherit.aes = FALSE, hjust = 0, size = 3, show.legend = FALSE
    ) +
    scale_color_manual(name = "Method", values = c("GAM" = "#16A085", "Stratified" = "#C0392B")) +
    coord_equal() +
    scale_x_continuous("False positive rate (1 - specificity)", limits = c(0, 1)) +
    scale_y_continuous("True positive rate (sensitivity)", limits = c(0, 1)) +
    labs(title = sprintf("Global ROC - %s", pair$name),
         subtitle = sprintf("%s | 90%% CI lower bound", version_title_label("original")))
}

# Generates one ROC figure per measure (IOR, AC)
generate_roc_plots <- function() {
  for (pair in roc_pairs) {
    p <- plot_roc_pair(roc_global, pair)
    if (!is.null(p)) {
      save_plot(p, sprintf("roc_global_%s", tolower(pair$name)), width = 8, height = 7)
    }
  }
}

generate_roc_plots()

################################################################################
# Positive triplet plots
################################################################################

# File paths for positive triplet data 
# Drug canonicalization comes from the shared OMOP vocabulary via build_drug_translation_table() (00_functions.R)
ruta_positive_results <- paste0("./results/", suffix, "/augmentation_results/positive_triplets_results.rds")
output_dir_positive <- paste0(fig_output_dir, "positive_triplets/")

################################################################################
# Loading and preprocessing functions
################################################################################

# Load and preprocess positive triplet data.
#
# Arguments:
#   ruta_results : path to the RDS file with augmentation results
#   ruta_ade: path to the CSV file with raw ADE reports
#
# Returns: list with two elements:
#   $results: data.table filtered to baseline cases (reduction_pct == 0, successful)
#   $ade_raw: preprocessed ADE data.table
#
# Only baseline, fully successful triplets are retained for plotting

load_positive_data <- function(ruta_results, ruta_ade) {
  dt <- readRDS(ruta_results)

  # Keeps only baseline (no reduction) and fully successful cases
  dt <- dt[reduction_pct == 0 & model_success == TRUE & injection_success == TRUE]

  ade_raw_dt <- fread(ruta_ade)

  # Process sex variable if sex-stratified analysis is enabled
  if (include_sex) {
    ade_raw_dt[, sex := toupper(trimws(sex))]
    ade_raw_dt[sex == "M", sex := "MALE"]
    ade_raw_dt[sex == "F", sex := "FEMALE"]
    ade_raw_dt[, sex := factor(sex, levels = c("MALE", "FEMALE"))]

    sex_summary <- ade_raw_dt[, .(n = .N), by = sex]
    message("  Sex distribution:")
    print(sex_summary)
  }

  message(sprintf("  ADE dataset: %s rows", format(nrow(ade_raw_dt), big.mark = ",")))

  return(list(results = dt, ade_raw = ade_raw_dt))
}

# Apply drug ID canonicalization to the ADE dataset and set keys for fast lookup.
#
# Arguments:
#  ade_raw_dt: raw ADE data.table
#  translation_table: data.table mapping atc_concept_id -> canonical_id
#    (from build_drug_translation_table())
#
# Returns: processed data.table with key set on (atc_concept_id, meddra_concept_id, nichd_num)

translate_ade <- function(ade_raw_dt, translation_table) {
  # Character cast needed for the join key; vocabulary IDs may arrive as integer
  ade_raw_dt[, atc_concept_id := as.character(atc_concept_id)]

  ade_raw_dt <- merge(
    ade_raw_dt,
    translation_table[, .(atc_concept_id, canonical_id)],
    by = "atc_concept_id",
    all.x = TRUE
  )

  # Replaces original IDs with their canonical counterparts
  ade_raw_dt[!is.na(canonical_id), atc_concept_id := canonical_id]
  ade_raw_dt[, canonical_id := NULL]

  # Canonicalization may collapse multiple original IDs to one; deduplicate
  nrow_before <- nrow(ade_raw_dt)
  unique_cols <- c("safetyreportid", "atc_concept_id", "meddra_concept_id")
  if (include_sex) unique_cols <- c(unique_cols, "sex")

  ade_raw_dt <- unique(ade_raw_dt, by = unique_cols)
  message(sprintf("  Duplicates removed: %d", nrow_before - nrow(ade_raw_dt)))

  # Ordered factor preserves stage sequence; numeric version used as join key
  ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
  ade_raw_dt[, nichd_num := as.integer(nichd)]

  # Sets key for efficient subset lookups downstream
  setkey(ade_raw_dt, atc_concept_id, meddra_concept_id, nichd_num)

  return(ade_raw_dt)
}

################################################################################
# Report count calculation
################################################################################

# Count spontaneous reports by NICHD stage for a specific drug-drug-event triplet.
#
# Arguments:
#  drug_a: ATC concept ID of drug A
#  drug_b: ATC concept ID of drug B
#  meddra_event: MedDRA concept ID of the adverse event
#  ade_dt: preprocessed ADE data.table (output of translate_ade)
#
# Returns: data.table with 7 rows (one per stage) and columns
#   n_a, n_b, n_ab, n_event, n_event_ab
#
# Counts are computed on the fly from the ADE table rather than precomputed upstream.

calculate_triplet_counts <- function(drug_a, drug_b, meddra_event, ade_dt) {
  # Unique report IDs for each drug and the event
  ids_a <- unique(ade_dt[atc_concept_id == drug_a, safetyreportid])
  ids_b <- unique(ade_dt[atc_concept_id == drug_b, safetyreportid])
  ids_event <- unique(ade_dt[meddra_concept_id == meddra_event, safetyreportid])

  # Co-exposure and exclusive-exposure report sets
  ids_ab <- intersect(ids_a, ids_b)
  ids_a_only <- setdiff(ids_a, ids_b)
  ids_b_only <- setdiff(ids_b, ids_a)

  # Helper: count unique reports per NICHD stage for a given ID subset
  count_per_stage <- function(ids_subset) {
    if (length(ids_subset) == 0) {
      return(data.table(nichd_num = 1:7, n = 0L))
    }
    counts <- ade_dt[safetyreportid %in% ids_subset, .(n = uniqueN(safetyreportid)), by = nichd_num]
    # Fill in any missing stages with zero counts
    all_stages <- data.table(nichd_num = 1:7)
    counts <- merge(all_stages, counts, by = "nichd_num", all.x = TRUE)
    counts[is.na(n), n := 0L]
    return(counts[order(nichd_num)])
  }

  n_a_dt <- count_per_stage(ids_a_only)
  n_b_dt <- count_per_stage(ids_b_only)
  n_ab_dt <- count_per_stage(ids_ab)
  n_event_dt <- count_per_stage(ids_event)
  n_event_ab_dt <- count_per_stage(intersect(ids_event, ids_ab))

  result <- data.table(
    nichd_num = 1:7,
    n_a = n_a_dt$n,
    n_b = n_b_dt$n,
    n_ab = n_ab_dt$n,
    n_event = n_event_dt$n,
    n_event_ab = n_event_ab_dt$n
  )

  return(result)
}

################################################################################
# Triplet data expansion
################################################################################

# Expands a single triplet result row into a stage-level data.table with metrics and report counts merged together.
#
# Arguments:
#  row: single-row data.table from the results object
#  ade_dt: preprocessed ADE data.table (output of translate_ade)
#
# Returns: long data.table with one row per NICHD stage, containing GAM and
#  classic metric estimates, confidence bounds, and report counts.
#
# Metrics are extracted from list columns in the RDS results.
# Injected counts (from the augmentation diagnostics) are added to the observed event counts before plotting

expand_triplets_counts <- function(row, ade_dt) {
  stages <- unlist(row$stage)

  # GAM metrics
  gam_log_ior <- unlist(row$log_ior)
  gam_log_ior_lower <- unlist(row$log_ior_lower90)
  gam_log_ior_upper <- gam_log_ior + (gam_log_ior - gam_log_ior_lower)

  gam_ac <- unlist(row$ac_values)
  gam_ac_lower <- unlist(row$ac_lower90)
  gam_ac_upper <- unlist(row$ac_upper90)

  # Classic (stratified) metrics
  cls_log_ior <- unlist(row$log_ior_classic)
  cls_log_ior_lower <- unlist(row$log_ior_classic_lower90)
  cls_log_ior_upper <- cls_log_ior + (cls_log_ior - cls_log_ior_lower)

  cls_ac <- unlist(row$AC_classic)
  cls_ac_lower <- unlist(row$AC_classic_lower90)
  cls_ac_upper <- unlist(row$AC_classic_upper90)

  # Report counts computed on the fly from the ADE table
  counts_dt <- calculate_triplet_counts(row$drugA, row$drugB, row$meddra, ade_dt)

  # Extract injected counts per stage from diagnostics$injection_by_stage
  inj_by_stage <- tryCatch({
    diag <- row$diagnostics[[1]]
    if (!is.null(diag) && !is.null(diag$injection_by_stage)) {
      ibs <- merge(data.table(nichd_num = 1:7), diag$injection_by_stage, by = "nichd_num", all.x = TRUE)
      ibs[is.na(N), N := 0L]
      ibs[order(nichd_num)]$N
    } else { rep(0L, 7) }
  }, error = function(e) rep(0L, 7))

  # Add injected reports to observed event counts for plotting purposes
  counts_dt[, n_event_ab := n_event_ab + inj_by_stage]
  counts_dt[, n_event := n_event + inj_by_stage]

  # Guard against mismatched vector lengths across sources
  n <- min(length(stages), length(gam_log_ior), nrow(counts_dt))
  if (n == 0) return(NULL)

  result <- data.table(
    triplet_id = row$triplet_id,
    dynamic = row$dynamic,
    stage_num = stages[1:n],
    # GAM metrics
    gam_log_ior = gam_log_ior[1:n],
    gam_log_ior_lower = gam_log_ior_lower[1:n],
    gam_log_ior_upper = gam_log_ior_upper[1:n],
    gam_ac = gam_ac[1:n],
    gam_ac_lower = gam_ac_lower[1:n],
    gam_ac_upper = gam_ac_upper[1:n],
    # Classic metrics
    cls_log_ior = cls_log_ior[1:n],
    cls_log_ior_lower = cls_log_ior_lower[1:n],
    cls_log_ior_upper = cls_log_ior_upper[1:n],
    cls_ac = cls_ac[1:n],
    cls_ac_lower = cls_ac_lower[1:n],
    cls_ac_upper = cls_ac_upper[1:n],
    # Report counts
    n_a = counts_dt$n_a[1:n],
    n_b = counts_dt$n_b[1:n],
    n_ab = counts_dt$n_ab[1:n],
    n_event = counts_dt$n_event[1:n],
    n_event_ab = counts_dt$n_event_ab[1:n]
  )

  return(result)
}

################################################################################
# Triplet-level visualization
################################################################################

# Generate a dual-axis plot for a single triplet 
# metric trajectory (primary Y) overlaid on report count bars (secondary log Y)
#
# Arguments:
#  plot_dt: expanded stage-level data.table for this triplet
#  metric_col: column name of the main metric values
#  lower_col: column name of the lower confidence bound
#  upper_col: column name of the upper confidence bound
#  y_label: label for the primary Y axis
#  file_suffix: identifier string used in the output filename
#  y_limit: symmetric Y range for the primary axis (default +/-10)
#  max_count: expected maximum report count, used to calibrate the log scale (default 5000)
#  plot_title: optional plot title (defaults to dynamic + triplet ID)
#  plot_subtitle: optional subtitle (defaults to "A + B")
#  stratified: 
#     If TRUE, render the metric as per-stage error bars 
#     If FALSE (GAM), use the ribbon + connecting line

graph_metrics_counts <- function(plot_dt, metric_col, lower_col, upper_col,
                                  y_label, file_suffix, y_limit = 10, max_count = 5000,
                                  plot_title = NULL, plot_subtitle = NULL,
                                  stratified = FALSE) {
  # Pivot count columns to long format for bar rendering
  counts_long <- melt(
    plot_dt,
    id.vars = "stage_num",
    measure.vars = c("n_a", "n_b", "n_ab", "n_event", "n_event_ab"),
    variable.name = "metric",
    value.name = "count"
  )

  # Legend labels for each count group
  count_labels <- c(
    n_a = "A", n_b = "B", n_ab = "A-B",
    n_event = "Event", n_event_ab = "A-B-Event"
  )
  counts_long[, metric_label := factor(count_labels[as.character(metric)], levels = unname(count_labels))]

  # Compute dodged bar X positions manually (5 groups, evenly spaced within each stage)
  bar_w <- 0.12
  counts_long[, metric_idx := as.integer(metric)]
  counts_long[, x_center := stage_num + (metric_idx - 3) * bar_w]
  counts_long[, xmin := x_center - bar_w/2]
  counts_long[, xmax := x_center + bar_w/2]

  # Log-linear transformation for the secondary axis:
  # maps count = 1 -> -y_limit and count = max_count -> +y_limit
  transform_to_log <- function(count) {
    count_adj <- pmax(count, 1)
    -y_limit + (log10(count_adj) / log10(max_count)) * (2 * y_limit)
  }

  transform_from_log <- function(y) {
    ratio <- (y + y_limit) / (2 * y_limit)
    10^(ratio * log10(max_count))
  }

  counts_long[, ymax := transform_to_log(count)]
  counts_long[, ymin := -y_limit]

  # Metric CI rendering depends on method nature:
  # - Stratified: per-stage error bars and no connecting line, since each NICHD stage is estimated independently.
  # - GAM: CI ribbon plus a connecting line, reflecting the smooth cross-stage function
  metric_ci_layers <- if (stratified) {
    list(
      geom_errorbar(
        data = plot_dt,
        aes(x = stage_num, ymin = .data[[lower_col]], ymax = .data[[upper_col]]),
        color = "#2c3e50", width = 0.2, linewidth = 0.9
      )
    )
  } else {
    list(
      geom_ribbon(
        data = plot_dt,
        aes(x = stage_num, ymin = .data[[lower_col]], ymax = .data[[upper_col]]),
        alpha = 0.3, fill = "#2c3e50"
      ),
      geom_line(
        data = plot_dt,
        aes(x = stage_num, y = .data[[metric_col]]),
        color = "#2c3e50", linewidth = 1.2
      )
    )
  }

  p <- ggplot() +
    # Report count bars (secondary axis, log-scaled)
    geom_rect(
      data = counts_long,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = metric_label),
      alpha = 0.6, color = NA
    ) +
    # Zero reference line
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    # Metric CI (error bars for stratified, ribbon + line for GAM)
    metric_ci_layers +
    # Metric value points
    geom_point(
      data = plot_dt,
      aes(x = stage_num, y = .data[[metric_col]]),
      color = "#2c3e50", size = 4, fill = "white", shape = 21, stroke = 1.5
    ) +
    # Axis scales
    scale_x_continuous(
      breaks = 1:7,
      labels = nichd_labels,
      name = NULL
    ) +
    scale_y_continuous(
      name = y_label,
      limits = c(-y_limit, y_limit),
      breaks = {
        # Compute readable breaks regardless of the axis range, anchored at 0
        step <- 10^floor(log10(y_limit / 4))
        step <- step * if (y_limit / step > 8) 2 else 1
        pos_breaks <- seq(0, y_limit, by = step)
        sort(unique(c(-pos_breaks, pos_breaks)))
      },
      sec.axis = sec_axis(
        trans = transform_from_log,
        name = "Reports (log)",
        breaks = c(1, 10, 100, 1000, 10000),
        labels = function(x) format(x, big.mark = ",", scientific = FALSE)
      )
    ) +
    scale_fill_brewer(palette = "Set2", name = "") +
    # Plot labels
    labs(
      title = if (!is.null(plot_title)) plot_title
      else paste0(plot_dt$dynamic[1], " | Triplet: ", plot_dt$triplet_id[1]),
    subtitle = if (!is.null(plot_subtitle)) plot_subtitle else "A + B"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y.right = element_text(size = 8),
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "gray30"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 8)
    )
  return(p)
}

################################################################################
# Positive triplet plot generation
################################################################################

# Save a single positive triplet plot to disk.
#
# Arguments:
#  p: ggplot object
#  triplet_id: numeric triplet identifier
#  dynamic: dynamic label string (used in the filename)
#  suffix: method/metric identifier string (e.g. "GAM_LogIOR")
#  dir: output directory (defaults to output_dir_positive)

save_positive_graph <- function(p, triplet_id, dynamic, suffix,
                                 dir = output_dir_positive) {
  safe_name <- sprintf("Triplet_%d_%s_%s", triplet_id, dynamic, suffix)
  safe_name <- gsub("[^a-zA-Z0-9._-]", "_", safe_name)
  ggsave(
    filename = file.path(dir, paste0(safe_name, ".png")),
    plot = p,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
}

# Generate plots for all positive triplets from the augmentation pipeline.
#
# Source: positive_triplets_results.rds (baseline only: reduction_pct == 0)
# Counts: computed on the fly from ade_raw via calculate_triplet_counts()
# Metrics: GAM-LogIOR, GAM-AC, Classic-LogIOR, Classic-AC

generate_positive_graphs_from_results <- function() {
  dir.create(output_dir_positive, showWarnings = FALSE, recursive = TRUE)

  data_pos <- load_positive_data(ruta_positive_results, ruta_ade_raw)
  dt <- data_pos$results
  trans_table <- build_drug_translation_table()
  ade_processed <- translate_ade(data_pos$ade_raw, trans_table)

  pb <- txtProgressBar(min = 0, max = nrow(dt), style = 3)
  generated_graphs <- 0

  for (i in seq_len(nrow(dt))) {
    tryCatch({
      row <- dt[i]

      ptitle <- paste0(row$dynamic, " | Triplet: ", row$triplet_id)
      psubtitle <- "A + B"

      plot_dt <- expand_triplets_counts(row, ade_dt = ade_processed)
      if (is.null(plot_dt) || nrow(plot_dt) == 0) {
        message(sprintf("\n  Skipping triplet %d: no expandable data", i))
        next
      }

      if (any(is.finite(plot_dt$gam_log_ior_lower))) {
        p <- graph_metrics_counts(plot_dt, "gam_log_ior", "gam_log_ior_lower", "gam_log_ior_upper",
                                   "Log-IOR (GAM, 90% CI)", "GAM_LogIOR",
                                   plot_title = ptitle, plot_subtitle = psubtitle)
        save_positive_graph(p, plot_dt$triplet_id[1], plot_dt$dynamic[1], "GAM_LogIOR")
        generated_graphs <- generated_graphs + 1
      }

      if (any(is.finite(plot_dt$gam_ac_lower))) {
        ac_vals_gam <- c(plot_dt$gam_ac, plot_dt$gam_ac_lower, plot_dt$gam_ac_upper)
        ac_max_abs_gam <- max(abs(ac_vals_gam[is.finite(ac_vals_gam)]), na.rm = TRUE)
        # Additive contrast is a reporting-proportion difference; autoscale (no integer floor).
        y_limit_gam_ac <- max(ac_max_abs_gam * 1.2, 0.01)
        p <- graph_metrics_counts(plot_dt, "gam_ac", "gam_ac_lower", "gam_ac_upper",
                             "Additive contrast (GAM, 90% CI)", "GAM_AC",
                             y_limit = y_limit_gam_ac,
                             plot_title = ptitle, plot_subtitle = psubtitle)
        save_positive_graph(p, plot_dt$triplet_id[1], plot_dt$dynamic[1], "GAM_AC")
        generated_graphs <- generated_graphs + 1
      }

      if (any(is.finite(plot_dt$cls_log_ior_lower))) {
        p <- graph_metrics_counts(plot_dt, "cls_log_ior", "cls_log_ior_lower", "cls_log_ior_upper",
                                   "Log-IOR (Stratified, 90% CI)", "Classic_LogIOR",
                                   plot_title = ptitle, plot_subtitle = psubtitle,
                                   stratified = TRUE)
        save_positive_graph(p, plot_dt$triplet_id[1], plot_dt$dynamic[1], "Classic_LogIOR")
        generated_graphs <- generated_graphs + 1
      }

      if (any(is.finite(plot_dt$cls_ac_lower))) {
        ac_vals_cls <- c(plot_dt$cls_ac, plot_dt$cls_ac_lower, plot_dt$cls_ac_upper)
        ac_max_abs_cls <- max(abs(ac_vals_cls[is.finite(ac_vals_cls)]), na.rm = TRUE)
        # Additive contrast is a reporting-proportion difference; autoscale (no integer floor).
        y_limit_cls_ac <- max(ac_max_abs_cls * 1.2, 0.01)
        p <- graph_metrics_counts(plot_dt, "cls_ac", "cls_ac_lower", "cls_ac_upper",
                             "Additive contrast (Stratified, 90% CI)", "Classic_AC",
                             y_limit = y_limit_cls_ac,
                             plot_title = ptitle, plot_subtitle = psubtitle,
                             stratified = TRUE)
        save_positive_graph(p, plot_dt$triplet_id[1], plot_dt$dynamic[1], "Classic_AC")
        generated_graphs <- generated_graphs + 1
      }

    }, error = function(e) {
      message(sprintf("\n  Error in triplet %d: %s", i, e$message))
    })

    if (i %% 50 == 0) gc()
    setTxtProgressBar(pb, i)
  }
  close(pb)
  message(sprintf("\n  Generated graphs: %d", generated_graphs))
  return(invisible(generated_graphs))
}

generate_positive_graphs_from_results()
