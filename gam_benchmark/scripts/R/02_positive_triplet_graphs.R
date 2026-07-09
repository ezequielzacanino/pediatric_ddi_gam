################################################################################
# Per-triplet visual comparison of the curated positive controls.
# Script: 02_positive_triplet_graphs.R
#
# Reuses the dual-axis triplet plotting system from gam_validation/40_graphs.R
################################################################################

source("00_functions.R", local = TRUE)

pacman::p_load(ggplot2)

################################################################################
# Paths and configuration
################################################################################

# Roll-up level must match the one used by script 01 to name the fit cache.
meddra_rollup_level <- "HLT"

benchmark_dir <- file.path("./results", suffix, "benchmark_validation")
benchmark_fit_rds_path <- file.path(
  benchmark_dir,
  paste0("benchmark_fit_results_", tolower(meddra_rollup_level), ".rds")
)
output_dir_positive <- file.path(benchmark_dir, "positive_triplet_figures")

# Display labels for NICHD developmental stages (shared with 40_graphs.R).
nichd_labels <- c(
  "term_neonatal" = "Term neonatal",
  "infancy" = "Infancy",
  "toddler" = "Toddler",
  "early_childhood" = "Early childhood",
  "middle_childhood" = "Middle childhood",
  "early_adolescence" = "Early adolescence",
  "late_adolescence" = "Late adolescence"
)

################################################################################
# Report count calculation
################################################################################

# Count spontaneous reports by NICHD stage for a drug-drug-event triplet.
#
# Arguments:
#  drug_a, drug_b: canonical ATC concept IDs of the two drugs
#  meddra_event: MedDRA concept ID of the adverse event (at the active roll-up)
#  ade_dt: preprocessed ADE data.table (output of load_ade_modeling_data)
#
# Returns: data.table with 7 rows (one per stage) and columns
#  n_a, n_b, n_ab, n_event, n_event_ab

calculate_triplet_counts <- function(drug_a, drug_b, meddra_event, ade_dt) {
  ids_a <- unique(ade_dt[atc_concept_id == drug_a, safetyreportid])
  ids_b <- unique(ade_dt[atc_concept_id == drug_b, safetyreportid])
  ids_event <- unique(ade_dt[meddra_concept_id == meddra_event, safetyreportid])

  # Co-exposure and exclusive-exposure sets (n_ab excludes monotherapy reports)
  ids_ab <- intersect(ids_a, ids_b)
  ids_a_only <- setdiff(ids_a, ids_b)
  ids_b_only <- setdiff(ids_b, ids_a)

  # Fill all 7 stages even when a stage has zero reports
  count_per_stage <- function(ids_subset) {
    if (length(ids_subset) == 0) {
      return(data.table(nichd_num = 1:7, n = 0L))
    }
    counts <- ade_dt[safetyreportid %in% ids_subset, .(n = uniqueN(safetyreportid)), by = nichd_num]
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

# Expand a single fit row into a stage-level data.table with metric estimates, confidence bounds, and report counts.
#
#
# Arguments:
#  row: single-row data.table from the benchmark fit results
#  ade_dt: preprocessed ADE data.table (output of load_ade_modeling_data)

expand_triplet_counts <- function(row, ade_dt) {
  stages <- unlist(row$stage)

  # GAM metrics. Only the lower-90 bound is stored
  gam_log_ior <- unlist(row$log_ior)
  gam_log_ior_lower <- unlist(row$log_ior_lower90)
  gam_log_ior_upper <- gam_log_ior + (gam_log_ior - gam_log_ior_lower)

  gam_ac <- unlist(row$ac_values)
  gam_ac_lower <- unlist(row$ac_lower90)
  gam_ac_upper <- unlist(row$ac_upper90)

  # Stratified (classic) metrics
  cls_log_ior <- unlist(row$log_ior_classic)
  cls_log_ior_lower <- unlist(row$log_ior_classic_lower90)
  cls_log_ior_upper <- cls_log_ior + (cls_log_ior - cls_log_ior_lower)

  cls_ac <- unlist(row$AC_classic)
  cls_ac_lower <- unlist(row$AC_classic_lower90)
  cls_ac_upper <- unlist(row$AC_classic_upper90)

  # Report counts computed on the fly from the ADE table
  counts_dt <- calculate_triplet_counts(row$drugA.x, row$drugB.x, row$meddra.x, ade_dt)

  # Guard against mismatched vector lengths across sources
  n <- min(length(stages), length(gam_log_ior), nrow(counts_dt))
  if (n == 0) return(NULL)

  result <- data.table(
    triplet_id = row$triplet_id,
    stage_num = stages[1:n],
    # GAM metrics
    gam_log_ior = gam_log_ior[1:n],
    gam_log_ior_lower = gam_log_ior_lower[1:n],
    gam_log_ior_upper = gam_log_ior_upper[1:n],
    gam_ac = gam_ac[1:n],
    gam_ac_lower = gam_ac_lower[1:n],
    gam_ac_upper = gam_ac_upper[1:n],
    # Stratified metrics
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

# Dual-axis plot for a single triplet: 
# metric trajectory (primary Y) overlaid on report-count bars (secondary log Y)
#
# Arguments:
#  plot_dt: expanded stage-level data.table for this triplet
#  metric_col: column name of the main metric values
#  lower_col: column name of the lower confidence bound
#  upper_col: column name of the upper confidence bound
#  y_label: label for the primary Y axis
#  y_limit: symmetric Y range for the primary axis (default +/-10)
#  max_count: expected maximum report count, calibrates the log scale (default 5000)
#  plot_title, plot_subtitle : optional title/subtitle

graph_metrics_counts <- function(plot_dt, metric_col, lower_col, upper_col,
                                 y_label, y_limit = 10, max_count = 5000,
                                 plot_title = NULL, plot_subtitle = NULL) {
  counts_long <- melt(
    plot_dt,
    id.vars = "stage_num",
    measure.vars = c("n_a", "n_b", "n_ab", "n_event", "n_event_ab"),
    variable.name = "metric",
    value.name = "count"
  )

  count_labels <- c(
    n_a = "A", n_b = "B", n_ab = "A-B",
    n_event = "Event", n_event_ab = "A-B-Event"
  )
  counts_long[, metric_label := factor(count_labels[as.character(metric)], levels = unname(count_labels))]

  # Manual dodging: 5 groups evenly spaced within each stage
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

  p <- ggplot() +
    # Report count bars (secondary log axis, rendered below the metric line)
    geom_rect(
      data = counts_long,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = metric_label),
      alpha = 0.6, color = NA
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    # 90% CI ribbon for the metric
    geom_ribbon(
      data = plot_dt,
      aes(x = stage_num, ymin = .data[[lower_col]], ymax = .data[[upper_col]]),
      alpha = 0.3, fill = "#2c3e50"
    ) +
    geom_line(
      data = plot_dt,
      aes(x = stage_num, y = .data[[metric_col]]),
      color = "#2c3e50", linewidth = 1.2
    ) +
    geom_point(
      data = plot_dt,
      aes(x = stage_num, y = .data[[metric_col]]),
      color = "#2c3e50", size = 4, fill = "white", shape = 21, stroke = 1.5
    ) +
    scale_x_continuous(
      breaks = 1:7,
      labels = nichd_labels,
      name = NULL
    ) +
    scale_y_continuous(
      name = y_label,
      limits = c(-y_limit, y_limit),
      breaks = {
        step <- 10^floor(log10(y_limit / 4))
        step <- step * if (y_limit / step > 8) 2 else 1
        seq(-y_limit, y_limit, by = step)
      },
      sec.axis = sec_axis(
        trans = transform_from_log,
        name = "Reports (log)",
        breaks = c(1, 10, 100, 1000, 10000),
        labels = function(x) format(x, big.mark = ",", scientific = FALSE)
      )
    ) +
    scale_fill_brewer(palette = "Set2", name = "") +
    labs(
      title = plot_title,
      subtitle = plot_subtitle
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

# Save a single positive-triplet plot. Filenames lead with the triplet id and drug pair
save_positive_graph <- function(p, triplet_id, drug1, drug2, method_suffix,
                                dir = output_dir_positive) {
  safe_name <- sprintf("Triplet_%s_%s_%s_%s", triplet_id, drug1, drug2, method_suffix)
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

################################################################################
# Positive triplet plot generation
################################################################################

# Generate the four method plots for every curated positive control.
#
# Source: benchmark fit RDS, filtered to positive, mapped, successfully modeled.
# A method panel is drawn only when its CI bounds are finite.

generate_positive_triplet_graphs <- function() {
  dir.create(output_dir_positive, showWarnings = FALSE, recursive = TRUE)

  fit_dt <- readRDS(benchmark_fit_rds_path)
  fit_dt <- fit_dt[control_type == "positive" & mapping_success == TRUE & model_success == TRUE]
  message(sprintf("  Positive triplets to plot: %d", nrow(fit_dt)))

  # Same ID space (canonical drugs, rolled-up MedDRA) the models were fitted on.
  ade_processed <- load_ade_modeling_data(rollup_level = meddra_rollup_level)
  setkey(ade_processed, atc_concept_id, meddra_concept_id, nichd_num)

  pb <- txtProgressBar(min = 0, max = nrow(fit_dt), style = 3)
  generated_graphs <- 0

  for (i in seq_len(nrow(fit_dt))) {
    tryCatch({
      row <- fit_dt[i]

      ptitle <- sprintf("%s + %s", row$drug1_name, row$drug2_name)
      psubtitle <- sprintf("Event: %s | Triplet %s", row$meddra_pt, row$triplet_id)

      plot_dt <- expand_triplet_counts(row, ade_dt = ade_processed)
      if (is.null(plot_dt) || nrow(plot_dt) == 0) {
        message(sprintf("\n  Skipping triplet %s: no expandable data", row$triplet_id))
        next
      }

      if (any(is.finite(plot_dt$gam_log_ior_lower))) {
        p <- graph_metrics_counts(plot_dt, "gam_log_ior", "gam_log_ior_lower", "gam_log_ior_upper",
                                  "Log-IOR (GAM, 90% CI)",
                                  plot_title = ptitle, plot_subtitle = psubtitle)
        save_positive_graph(p, row$triplet_id, row$drug1_name, row$drug2_name, "GAM_LogIOR")
        generated_graphs <- generated_graphs + 1
      }

      if (any(is.finite(plot_dt$gam_ac_lower))) {
        ac_vals_gam <- c(plot_dt$gam_ac, plot_dt$gam_ac_lower, plot_dt$gam_ac_upper)
        ac_max_abs_gam <- max(abs(ac_vals_gam[is.finite(ac_vals_gam)]), na.rm = TRUE)
        # Additive contrast is a reporting-proportion difference; autoscale (no integer floor).
        y_limit_gam_ac <- max(ac_max_abs_gam * 1.2, 0.01)
        p <- graph_metrics_counts(plot_dt, "gam_ac", "gam_ac_lower", "gam_ac_upper",
                                  "Additive contrast (GAM, 90% CI)", y_limit = y_limit_gam_ac,
                                  plot_title = ptitle, plot_subtitle = psubtitle)
        save_positive_graph(p, row$triplet_id, row$drug1_name, row$drug2_name, "GAM_AC")
        generated_graphs <- generated_graphs + 1
      }

      if (any(is.finite(plot_dt$cls_log_ior_lower))) {
        p <- graph_metrics_counts(plot_dt, "cls_log_ior", "cls_log_ior_lower", "cls_log_ior_upper",
                                  "Log-IOR (Stratified, 90% CI)",
                                  plot_title = ptitle, plot_subtitle = psubtitle)
        save_positive_graph(p, row$triplet_id, row$drug1_name, row$drug2_name, "Stratified_LogIOR")
        generated_graphs <- generated_graphs + 1
      }

      if (any(is.finite(plot_dt$cls_ac_lower))) {
        ac_vals_cls <- c(plot_dt$cls_ac, plot_dt$cls_ac_lower, plot_dt$cls_ac_upper)
        ac_max_abs_cls <- max(abs(ac_vals_cls[is.finite(ac_vals_cls)]), na.rm = TRUE)
        # Additive contrast is a reporting-proportion difference; autoscale (no integer floor).
        y_limit_cls_ac <- max(ac_max_abs_cls * 1.2, 0.01)
        p <- graph_metrics_counts(plot_dt, "cls_ac", "cls_ac_lower", "cls_ac_upper",
                                  "Additive contrast (Stratified, 90% CI)", y_limit = y_limit_cls_ac,
                                  plot_title = ptitle, plot_subtitle = psubtitle)
        save_positive_graph(p, row$triplet_id, row$drug1_name, row$drug2_name, "Stratified_AC")
        generated_graphs <- generated_graphs + 1
      }

    }, error = function(e) {
      message(sprintf("\n  Error in triplet %s: %s", fit_dt$triplet_id[i], e$message))
    })

    if (i %% 50 == 0) gc()
    setTxtProgressBar(pb, i)
  }
  close(pb)
  message(sprintf("\n  Generated graphs: %d in %s", generated_graphs, output_dir_positive))
  return(invisible(generated_graphs))
}

generate_positive_triplet_graphs()
