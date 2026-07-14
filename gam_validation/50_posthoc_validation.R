################################################################################
# Post-hoc methodological validation script
# Script: 50_posthoc_validation.R
#
# runs targeted tests on the validation design. 
#
# Tests:
#   1. Injection-success & shape-fidelity audit: 
#         how often and why the injection fails, and whether the realized per-stage rate 
#         reproduces the intended temporal shape.
#   2. Shape-fidelity per-stage sensitivity: 
#         per-stage operating characteristics restricted to triplets that demonstrably reproduced their shape.
#   3. Adversarial non-smooth injection benchmark: 
#         re-inject non-smooth shapes and compare GAM vs stratified, 
#         to test whether the smooth (tanh) injection design structurally favors the GAM.
#   4. Null adequacy: 
#         whether the permutation null is an adequate per-stage reference for the real negatives, 
#         for both measures (IOR, AC) and both methods (GAM, stratified).
#   5. Findings summary.
#
################################################################################

source("00_functions.R", local = TRUE)

################################################################################
# Configuration
################################################################################

# Refitting test parameters 
# n_adv_pairs base pairs x n_shapes fits at reduction 0 only.
run_adversarial <- TRUE          # set FALSE to skip the refitting test entirely
n_adv_pairs <- 500               # distinct base drug pairs to re-inject (refit cost driver).

ac_bootstrap_n_posthoc <- 500    # lighter than the main pipeline's (1000)
adv_seed <- 12345                # base seed for the adversarial experiment

# Shape-fidelity parameters
dynamic_realized_z_crit <- 1.645  # one-sided z (alpha 0.05) for a shape-faithful triplet
flatness_alpha <- 0.05            # uniform arm counts as "flat" if its chi-square p >= this
n_boot_posthoc <- 1000            # matches the pipeline's (1000) for the refit metric CIs

ruta_aug <- paste0("./results/", suffix, "/augmentation_results/")
ruta_null <- paste0("./results/", suffix, "/null_distribution_results/")

ruta_pos_rds <- paste0(ruta_aug, "positive_triplets_results.rds")
ruta_pos_meta <- paste0(ruta_aug, "positive_triplets_metadata.csv")
ruta_coadmin_pos <- paste0(ruta_aug, "positive_coadmin_by_stage.csv")
ruta_coadmin_neg <- paste0(ruta_aug, "negative_coadmin_by_stage.csv")
ruta_neg_rds <- paste0(ruta_aug, "negative_triplets_results.rds")
ruta_null_dist <- paste0(ruta_null, "null_distribution.csv")
# GAM and stratified per-stage permutation-null thresholds 
ruta_null_thr <- paste0(ruta_null, "null_thresholds.csv")
ruta_null_thr_ac <- paste0(ruta_null, "null_thresholds_ac.csv")
ruta_null_thr_cls_ior <- paste0(ruta_null, "null_thresholds_classic_ior.csv")
ruta_null_thr_cls_ac <- paste0(ruta_null, "null_thresholds_classic_ac.csv")

output_dir <- paste0("./results/", suffix, "/posthoc_validation_results/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

nichd_labels <- c(
  term_neonatal = "Term neonatal", infancy = "Infancy", toddler = "Toddler",
  early_childhood = "Early childhood", middle_childhood = "Middle childhood",
  early_adolescence = "Early adolescence", late_adolescence = "Late adolescence"
)

# Findings accumulator
findings <- list()
add_finding <- function(justification, test, metric, value, guide = "") {
  findings[[length(findings) + 1]] <<- data.table(
    justification = justification, test = test, metric = metric, value = value, guide = guide
  )
}

message(sprintf("Post-hoc validation. Outputs -> %s", output_dir))

################################################################################
# Shared helper: unwrap one diagnostics cell from the positive results
################################################################################

# diagnostics stores a named list per triplet (built in inject_signal)
unwrap_diag <- function(diag_cell) {
  if (is.null(diag_cell)) return(NULL)
  if (is.list(diag_cell) && length(diag_cell) >= 1 && is.list(diag_cell[[1]])) {
    return(diag_cell[[1]])
  }
  diag_cell
}

################################################################################
# Shared helpers: null thresholds and the core method set
################################################################################

# The four per-stage permutation-null thresholds (GAM + stratified, IOR + AC)
# merged on stage — the second leg of the double threshold. NULL if any is missing.
load_null_thresholds <- function() {
  files <- c(ruta_null_thr, ruta_null_thr_ac, ruta_null_thr_cls_ior, ruta_null_thr_cls_ac)
  if (!all(file.exists(files))) return(NULL)
  thr_col <- paste0("threshold_", percentil)
  Reduce(function(a, b) merge(a, b, by = "stage"), list(
    fread(ruta_null_thr)[, .(stage, threshold_ior = get(thr_col))],
    fread(ruta_null_thr_ac)[, .(stage, threshold_ac = get(thr_col))],
    fread(ruta_null_thr_cls_ior)[, .(stage, threshold_classic_ior = get(thr_col))],
    fread(ruta_null_thr_cls_ac)[, .(stage, threshold_classic_ac = get(thr_col))]
  ))
}

# GAM vs stratified, IOR and AC. Same score columns and names as 30_metrics.R
# detect_signal + calculate_metrics return the identical schema
metodos_core <- list(
  list(nombre = "GAM-logIOR",         tipo = "IOR",  measure = "IOR",
       score_type = "gam_log_ior_lower90",     score_type_auc = "gam_log_ior"),
  list(nombre = "GAM-AC",           tipo = "AC", measure = "AC",
       score_type = "gam_ac_lower90",          score_type_auc = "gam_ac"),
  list(nombre = "Estratificado-IOR",  tipo = "IOR",  measure = "IOR",
       score_type = "classic_log_ior_lower90", score_type_auc = "classic_log_ior"),
  list(nombre = "Estratificado-AC", tipo = "AC", measure = "AC",
       score_type = "classic_ac_lower90",      score_type_auc = "classic_ac")
)

# High-reporting stage classification by dynamic 
# read by expand() via the global stage_class.
stage_class_table <- function() rbind(
  data.table(nichd = niveles_nichd, dynamic = "uniform", class = 1),
  data.table(nichd = niveles_nichd, dynamic = "increase", class = c(0, 0, NA, NA, NA, 1, 1)),
  data.table(nichd = niveles_nichd, dynamic = "decrease", class = c(1, 1, NA, NA, NA, 0, 0)),
  data.table(nichd = niveles_nichd, dynamic = "plateau", class = c(0, NA, 1, 1, 1, NA, 0)),
  data.table(nichd = niveles_nichd, dynamic = "inverse_plateau", class = c(1, NA, 0, 0, 0, NA, 1))
)

################################################################################
# Injection-success & shape-fidelity audit
################################################################################
# Produce direct evidence about the injection process itself
#
# How often, and why, does the injection fail
# Does the realized per-stage injection reproduce the intended temporal shape

module1 <- function() {
  if (!file.exists(ruta_pos_rds)) {
    message(" positive results not found - skipping.")
    return(invisible(NULL))
  }
  message("\n Injection-success audit")

  pos <- readRDS(ruta_pos_rds)
  pos <- pos[reduction_pct == 0]  # base level carries the injection diagnostics

  # 1. Failure funnel: success rate and failure reasons per dynamic
  failure_reason <- vapply(seq_len(nrow(pos)), function(i) {
    if (isTRUE(pos$injection_success[i])) return("success")
    d <- unwrap_diag(pos$diagnostics[[i]])
    if (!is.null(d$reason)) as.character(d$reason) else "unknown"
  }, character(1))
  pos[, failure_reason := failure_reason]

  success_by_dyn <- pos[, .(
    n_total = .N,
    n_success = sum(injection_success == TRUE, na.rm = TRUE),
    success_rate = round(mean(injection_success == TRUE, na.rm = TRUE) * 100, 1)
  ), by = dynamic][order(-success_rate)]
  fwrite(success_by_dyn, paste0(output_dir, "m1_injection_success_by_dynamic.csv"))

  reason_tab <- pos[injection_success == FALSE, .N, by = .(dynamic, failure_reason)][order(dynamic, -N)]
  fwrite(reason_tab, paste0(output_dir, "m1_injection_failure_reasons.csv"))

  overall_success <- round(mean(pos$injection_success == TRUE, na.rm = TRUE) * 100, 1)
  add_finding("injection success rate", "pct_successful_injections",
              as.character(overall_success), ">=80 target")

  # 2. Shape fidelity: realized injected rate per stage vs intended p_dynamic
  # Realized rate p_hat_s = n_injected_s / n_coadmin_stage_s is the MLE of the per-stage Bernoulli injection probability the design intends 
  coadmin <- if (file.exists(ruta_coadmin_pos)) fread(ruta_coadmin_pos) else NULL

  pos_ok <- pos[injection_success == TRUE]
  shape_rows <- rbindlist(lapply(seq_len(nrow(pos_ok)), function(i) {
    d <- unwrap_diag(pos_ok$diagnostics[[i]])
    sp <- d$stage_probs
    ibs <- d$injection_by_stage
    if (is.null(sp) || !is.data.table(sp)) return(NULL)
    out <- data.table(nichd_num = 1:7)
    out <- merge(out, sp[, .(nichd_num, p_dynamic)], by = "nichd_num", all.x = TRUE)
    if (!is.null(ibs) && is.data.table(ibs) && nrow(ibs) > 0) {
      out <- merge(out, ibs[, .(nichd_num, n_injected = N)], by = "nichd_num", all.x = TRUE)
    } else {
      out[, n_injected := 0L]
    }
    out[is.na(n_injected), n_injected := 0L]
    out[, `:=`(triplet_id = pos_ok$triplet_id[i], dynamic = pos_ok$dynamic[i])]
    out
  }), fill = TRUE)

  if (!is.null(coadmin) && "n_coadmin_stage" %in% names(coadmin)) {
    coadmin_key <- coadmin[, .(triplet_id, nichd_num, n_coadmin_stage)]
    shape_rows <- merge(shape_rows, coadmin_key, by = c("triplet_id", "nichd_num"), all.x = TRUE)
    shape_rows[is.na(n_coadmin_stage) | n_coadmin_stage == 0, n_coadmin_stage := NA_real_]
    shape_rows[, injected_rate := n_injected / n_coadmin_stage]
  } else {
    shape_rows[, injected_rate := as.numeric(n_injected)]
  }
  has_coadmin <- "n_coadmin_stage" %in% names(shape_rows)

  # Uniform dynamic. there is no slope to reproduce, so the realized rate should be flat. 
  # Pearson goodness-of-fit of the per-stage injected counts to the pooled rate under the flat null X2 ~ chi-square on (usable stages - 1) df. 
  # A faithful uniform injection FAILS to reject (large p). 
  # is_flat marks the bases on which the other dynamics can be judged
  fidelity_uniform <- shape_rows[dynamic == "uniform", {
    ok <- if (has_coadmin) is.finite(n_injected) & is.finite(n_coadmin_stage) & n_coadmin_stage > 0 else rep(FALSE, .N)
    if (sum(ok) >= 3) {
      n_s <- n_coadmin_stage[ok]
      p_bar <- sum(n_injected[ok]) / sum(n_s)
      x2 <- if (p_bar > 0 && p_bar < 1) sum((n_injected[ok] - n_s * p_bar)^2 / (n_s * p_bar * (1 - p_bar))) else NA_real_
      df <- sum(ok) - 1L
      .(flatness_x2 = x2, flatness_df = df, flatness_p = if (is.finite(x2)) 1 - pchisq(x2, df) else NA_real_)
    } else {
      .(flatness_x2 = NA_real_, flatness_df = NA_integer_, flatness_p = NA_real_)
    }
  }, by = triplet_id]
  fidelity_uniform[, is_flat := is.finite(flatness_p) & flatness_p >= flatness_alpha]
  fwrite(fidelity_uniform, paste0(output_dir, "m1_uniform_flatness.csv"))
  pct_flat <- round(mean(fidelity_uniform$is_flat) * 100, 1)
  add_finding("uniform injection flatness", "pct_uniform_triplets_flat",
              as.character(pct_flat), "p>=alpha flat")

  # Map each triplet to its base triplet
  # keep the bases whose uniform sibling is flat
  meta_link <- unique(fread(ruta_pos_meta)[, .(triplet_id, base_triplet_id)])
  flat_base_ids <- merge(fidelity_uniform[is_flat == TRUE, .(triplet_id)], meta_link,
                         by = "triplet_id")$base_triplet_id

  # Non-uniform dynamics: signed directional contrast 
  # w_s is the intended shape centred to sum zero 
  # T = sum_s w_s * p_hat_s is positive when the realized rate slopes the intended way 
  # Under the flat null T has mean 0 and binomial variance
  # sum_s w_s^2 * p_bar(1 - p_bar) / n_coadmin_s, giving a one-sided z.
  fidelity_nonuniform <- shape_rows[dynamic != "uniform", {
    ok <- if (has_coadmin) is.finite(injected_rate) & is.finite(p_dynamic) & is.finite(n_coadmin_stage) & n_coadmin_stage > 0 else rep(FALSE, .N)
    if (sum(ok) >= 3) {
      w <- p_dynamic[ok] - mean(p_dynamic[ok])
      n_s <- n_coadmin_stage[ok]
      p_bar <- sum(n_injected[ok]) / sum(n_s)
      var_t <- sum(w^2 * p_bar * (1 - p_bar) / n_s)
      z <- if (var_t > 0 && sd(w) > 0) sum(w * injected_rate[ok]) / sqrt(var_t) else NA_real_
    } else {
      z <- NA_real_
    }
    .(shape_z = z, shape_p = if (is.finite(z)) 1 - pnorm(z) else NA_real_, n_stages = sum(ok))
  }, by = .(triplet_id, dynamic)]
  fidelity_nonuniform <- merge(fidelity_nonuniform, meta_link, by = "triplet_id", all.x = TRUE)
  fidelity_nonuniform[, base_uniform_flat := base_triplet_id %in% flat_base_ids]
  fidelity_nonuniform[, shape_pass := base_uniform_flat & is.finite(shape_z) & shape_z >= dynamic_realized_z_crit]
  fwrite(fidelity_nonuniform, paste0(output_dir, "m1_triplet_shape_fidelity.csv"))

  # Summary and figure use only eligible (flat-base) triplets.
  eligible <- fidelity_nonuniform[base_uniform_flat == TRUE & is.finite(shape_z)]
  fidelity_summary <- eligible[, .(
    n_triplets = .N,
    median_shape_z = round(median(shape_z), 3),
    pct_shape_faithful = round(mean(shape_pass) * 100, 1)
  ), by = dynamic][order(dynamic)]
  fwrite(fidelity_summary, paste0(output_dir, "m1_shape_fidelity_by_dynamic.csv"))

  pct_faithful <- round(mean(eligible$shape_pass) * 100, 1)
  add_finding("injection shape fidelity (flat-base only)",
              "pct_nonuniform_triplets_shape_faithful", as.character(pct_faithful),
              "z>=1.645 faithful")

  # Shape-fidelity subset  
  # non-uniform triplets that reproduced their slope + uniform triplets that stayed flat.
  subset_ids <- unique(c(fidelity_nonuniform[shape_pass == TRUE, triplet_id],
                         fidelity_uniform[is_flat == TRUE, triplet_id]))
  fwrite(data.table(triplet_id = subset_ids), paste0(output_dir, "m1_shape_fidelity_subset.csv"))

  # Figure: per-triplet shape-fidelity z by dynamic 
  # dashed line is the one-sided alpha threshold, boxes above it reproduced the intended slope
  p_fid <- ggplot(eligible, aes(x = dynamic, y = shape_z)) +
    geom_hline(yintercept = dynamic_realized_z_crit, linetype = "dashed", color = "#C0392B") +
    geom_boxplot(outlier.alpha = 0.25, fill = "#BDC3C7") +
    labs(title = "Per-triplet injection shape fidelity (Module 1)",
         subtitle = sprintf("Directional-contrast z vs flat null (flat-base triplets); dashed = z = %.3f",
                            dynamic_realized_z_crit),
         x = "Injected dynamic", y = "Shape-fidelity z") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(paste0(output_dir, "m1_fig_shape_fidelity.png"), p_fid, width = 10, height = 6, dpi = 300, bg = "white")

  message(sprintf("  [M1] success=%.1f uniform_flat=%.1f shape_faithful_flatbase=%.1f n_eligible=%d",
                  overall_success, pct_flat, pct_faithful, nrow(eligible)))
  invisible(NULL)
}

################################################################################
# Shape-fidelity sensitivity 
################################################################################
# Sensitivity analysis on the operating characteristics. 
# re-computes the per-stage metrics for each method on the subset of positives with intended shape
# Reuses expand / detect_signal / calculate_metrics

module_shape_fidelity <- function() {
  ruta_subset <- paste0(output_dir, "m1_shape_fidelity_subset.csv")
  needed <- c(ruta_pos_rds, ruta_neg_rds, ruta_coadmin_pos, ruta_coadmin_neg, ruta_subset)
  if (!all(file.exists(needed))) {
    message("\n required inputs not found (positives/negatives/coadmin/subset) - skipping.")
    return(invisible(NULL))
  }
  null_thresholds_dt <- load_null_thresholds()
  if (is.null(null_thresholds_dt)) {
    message("\n null thresholds incomplete - skipping (double threshold unavailable).")
    return(invisible(NULL))
  }
  message("\n Shape-fidelity sensitivity (per-stage metrics)")

  # expand() reads these from the global environment (mirrors 30_metrics.R).
  ruta_base_sensitivity <<- ruta_aug
  use_threshold_ior <<- TRUE
  use_threshold_ac <<- TRUE
  null_thresholds <<- null_thresholds_dt
  coadmin_stage_pos <<- fread(ruta_coadmin_pos)
  coadmin_stage_neg <<- fread(ruta_coadmin_neg)
  # Ensure stage_num exists for older outputs that used nichd_num only.
  if (!"stage_num" %in% names(coadmin_stage_pos)) coadmin_stage_pos[, stage_num := nichd_num]
  if (!"stage_num" %in% names(coadmin_stage_neg)) coadmin_stage_neg[, stage_num := nichd_num]
  stage_class <<- stage_class_table()

  data0 <- expand(0)
  pos_high <- data0$pos_high
  neg_high <- data0$neg_high
  subset_ids <- fread(ruta_subset)$triplet_id
  pos_sub <- pos_high[triplet_id %in% subset_ids]
  n_all <- uniqueN(pos_high$triplet_id); n_sub <- uniqueN(pos_sub$triplet_id)

  # Per-stage metrics (no triplet aggregation) for one positive set, double threshold.
  per_stage_metrics <- function(pos_set, set_label) {
    rbindlist(lapply(metodos_core, function(met) {
      rbindlist(lapply(1:7, function(s) {
        dt <- rbind(pos_set[stage_num == s & class == 1], neg_high[stage_num == s], fill = TRUE)
        if (nrow(dt) == 0 || uniqueN(dt$label) < 2) return(NULL)
        dt <- detect_signal(dt, met$nombre, met$tipo, use_null = TRUE)
        m <- calculate_metrics(dt, n_boot_posthoc, aggregate_triplet = FALSE,
                               score_type = met$score_type, score_type_auc = met$score_type_auc)
        m[, `:=`(set = set_label, method = met$nombre, measure = met$measure,
                 stage = s, nichd = niveles_nichd[s])]
        m
      }), fill = TRUE)
    }), fill = TRUE)
  }

  metrics <- rbindlist(list(
    per_stage_metrics(pos_high, "all_success"),
    per_stage_metrics(pos_sub,  "shape_fidelity")
  ), fill = TRUE)
  setcolorder(metrics, c("set", "method", "measure", "stage", "nichd"))
  fwrite(metrics, paste0(output_dir, "m2_shape_fidelity_perstage_metrics.csv"))

  # GAM vs stratified 
  nichd_labels_m2 <- c(
    term_neonatal = "Term neonatal", infancy = "Infancy", toddler = "Toddler",
    early_childhood = "Early childhood", middle_childhood = "Middle childhood",
    early_adolescence = "Early adolescence", late_adolescence = "Late adolescence"
  )
  metric_labels_m2 <- c(AUC = "AUC", sensitivity = "Sensitivity", specificity = "Specificity")

  plot_m2_measure <- function(meas, gam_name, cls_name, gam_label, cls_label) {
    d <- metrics[set == "shape_fidelity" & method %in% c(gam_name, cls_name)]
    if (nrow(d) == 0) return(invisible(NULL))
    d[, nichd := factor(nichd, levels = niveles_nichd)]
    method_levels <- c(gam_name, cls_name)
    method_disp <- c(gam_label, cls_label)
    # Long format: one row per method x stage x metric, with 90% CI bounds.
    dl <- rbindlist(lapply(names(metric_labels_m2), function(mm) {
      data.table(
        method = factor(d$method, levels = method_levels),
        nichd = d$nichd,
        metric_label = factor(metric_labels_m2[[mm]], levels = unname(metric_labels_m2)),
        value = d[[mm]],
        lower = if (paste0(mm, "_lower") %in% names(d)) d[[paste0(mm, "_lower")]] else NA_real_,
        upper = if (paste0(mm, "_upper") %in% names(d)) d[[paste0(mm, "_upper")]] else NA_real_
      )
    }))
    p <- ggplot(dl, aes(x = nichd, y = value, color = method, group = method)) +
      geom_point(size = 2.2, position = position_dodge(width = 0.4)) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, alpha = 0.8,
                    na.rm = TRUE, position = position_dodge(width = 0.4)) +
      facet_grid(metric_label ~ ., scales = "free_y") +
      scale_x_discrete(labels = nichd_labels_m2, name = NULL) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 4), name = "Metric value",
                         expand = expansion(mult = 0.08)) +
      scale_color_manual(values = setNames(c("#16A085", "#C0392B"), method_levels),
                         labels = setNames(method_disp, method_levels), name = "Method") +
      labs(title = sprintf("Shape-fidelity subset per-stage operating characteristics - %s", meas),
           subtitle = sprintf("%s vs %s (double threshold)", gam_label, cls_label)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0(output_dir, sprintf("m2_fig_shape_fidelity_perstage_%s.png", tolower(meas))),
           p, width = 10, height = 9, dpi = 300, bg = "white")
  }
  plot_m2_measure("IOR", "GAM-logIOR", "Estratificado-IOR", "GAM-IOR", "Stratified-IOR")
  plot_m2_measure("AC", "GAM-AC", "Estratificado-AC", "GAM-AC", "Stratified-AC")

  ga <- metrics[set == "all_success"    & method == "GAM-logIOR", mean(AUC, na.rm = TRUE)]
  gs <- metrics[set == "shape_fidelity" & method == "GAM-logIOR", mean(AUC, na.rm = TRUE)]
  add_finding("shape-fidelity subset size", "n_triplets_subset_vs_all",
              sprintf("%d/%d", n_sub, n_all), "subset/all")
  add_finding("GAM-logIOR per-stage AUC: subset vs all", "mean_AUC_subset_minus_all",
              as.character(round(gs - ga, 3)), "higher=cleaner")

  message(sprintf("  [M2] subset=%d/%d GAM-logIOR mean per-stage AUC all=%.3f subset=%.3f",
                  n_sub, n_all, ga, gs))
  invisible(NULL)
}

################################################################################
# Adversarial non-smooth injection benchmark 
################################################################################
# Re-injects NON-smooth shapes that a spline does not favor 

# Non-smooth shapes, normalized to [-1, 1] like generate_dynamic()
# they plug into the same p_dynamic = t_ij * (1 + shape) used by inject_signal()
get_shape <- function(name, spike_stage = 4L, step_rising = TRUE) {
  switch(name,
    none = rep(0, 7),                                                   # no injection (negative arm)
    spike = { v <- rep(-1, 7); v[spike_stage] <- 1; v },                 # single-stage impulse at spike_stage
    step = if (step_rising) c(-1, -1, -1, 1, 1, 1, 1) else c(1, 1, 1, 1, -1, -1, -1),  # jump up or down
    sawtooth = c( 1, -1,  1, -1,  1, -1,  1),                                        # high-frequency alternation
    stop(sprintf("unknown shape: %s", name))
  )
}

adv_shapes <- c("none", "spike", "step", "sawtooth")
adv_active_stages <- list(none = integer(0), spike = 4L,
                          step = c(4L, 5L, 6L, 7L), sawtooth = c(1L, 3L, 5L, 7L))

# Worker
process_adv_pair <- function(pair_row, ade_raw_dt, shapes, get_shape_fn, bootstrap_n) {
  # 1. report-level table, shared by all arms of this pair.
  eval_real <- build_eval_table(ade_raw_dt, pair_row$drugA, pair_row$drugB,
                                pair_row$meddra, integer(0), include_sex = FALSE)
  if (nrow(eval_real) == 0 || sum(eval_real$droga_ab) < 3) return(NULL)

  # 2. Baseline rates from the real data, matching inject_signal()'s additive e_j.
  pA <- mean(eval_real[droga_a == 1 & droga_b == 0, ea_ocurrio])
  pB <- mean(eval_real[droga_a == 0 & droga_b == 1, ea_ocurrio])
  if (!is.finite(pA)) pA <- 0; if (!is.finite(pB)) pB <- 0
  e_j <- pA + pB - pA * pB
  t_ij <- pair_row$fold_change * e_j   # fold_change fixed per pair -> only the shape varies
  clamp <- function(p) pmin(pmax(p, 1e-4), 1 - 1e-4)

  out <- list()
  for (sh in shapes) {
    shp <- get_shape_fn(sh, pair_row$spike_stage, pair_row$step_rising)
    active <- which(shp > 0)            # active stages = the elevated portion of the shape
    p_dyn <- if (sh == "none") rep(0, 7) else clamp(t_ij * (1 + shp))

    # 3. Builds the augmented table: flip eligible co-admin rows with prob p_dyn
    # Requires >= 1 injected event in the active stages (injection_success).
    eval_aug <- copy(eval_real)
    if (sh != "none") {
      pmap <- data.table(nichd_num = 1:7, p_dyn = p_dyn)
      eval_aug <- merge(eval_aug, pmap, by = "nichd_num", all.x = TRUE)
      flip <- eval_aug$droga_ab == 1 & eval_aug$ea_ocurrio == 0 &
        runif(nrow(eval_aug)) < eval_aug$p_dyn
      if (sum(flip & eval_aug$nichd_num %in% active) < 1) next  # failed injection -> drop arm
      eval_aug[flip == TRUE, ea_ocurrio := 1L]
      eval_aug[, p_dyn := NULL]
    }

    # 4. Fits GAM and stratified estimators on the augmented table (reduction 0).
    g <- tryCatch(fit_gam(pair_row$drugA, pair_row$drugB, pair_row$meddra,
                          ade_data = NULL, spline_individuales = spline_individuales,
                          include_sex = include_sex, include_stage_sex = include_stage_sex,
                          k_spline = k_spline, bs_type = bs_type, select = select,
                          nichd_spline = nichd_spline, eval_dt = eval_aug),
                  error = function(e) list(success = FALSE))
    if (!isTRUE(g$success)) next
    ci <- tryCatch(calculate_classic_ior(pair_row$drugA, pair_row$drugB, pair_row$meddra,
                                         ade_data = NULL, eval_dt = eval_aug),
                   error = function(e) list(success = FALSE))
    cr <- tryCatch(calculate_classic_ac(pair_row$drugA, pair_row$drugB, pair_row$meddra,
                                          ade_data = NULL, eval_dt = eval_aug),
                   error = function(e) list(success = FALSE))

    res <- data.table(
      base_pair_id = pair_row$base_pair_id, shape = sh, stage_num = 1:7,
      gam_log_ior = g$log_ior, gam_log_ior_lower90 = g$log_ior_lower90,
      gam_ac = g$ac_values, gam_ac_lower90 = g$ac_lower90,
      classic_log_ior = if (isTRUE(ci$success)) ci$results_by_stage$log_ior_classic else NA_real_,
      classic_log_ior_lower90 = if (isTRUE(ci$success)) ci$results_by_stage$log_ior_classic_lower90 else NA_real_,
      classic_ac = if (isTRUE(cr$success)) cr$results_by_stage$AC_classic else NA_real_,
      classic_ac_lower90 = if (isTRUE(cr$success)) cr$results_by_stage$AC_classic_lower90 else NA_real_
    )
    res[, is_active := stage_num %in% active]
    out[[length(out) + 1]] <- res
  }
  if (length(out) == 0) return(NULL)
  rbindlist(out, fill = TRUE)
}

# intended adversarial shapes 
plot_adv_shapes <- function() {
  sh_dt <- rbindlist(lapply(setdiff(adv_shapes, "none"), function(sh) {
    data.table(shape = sh, nichd_num = 1:7, value = 1 + get_shape(sh),
               active = (1:7) %in% adv_active_stages[[sh]])
  }))
  sh_dt[, nichd := factor(niveles_nichd[nichd_num], levels = niveles_nichd, ordered = TRUE)]
  p <- ggplot(sh_dt, aes(x = nichd, y = value, group = shape)) +
    geom_line(color = "#34495E", linewidth = 0.9) +
    geom_point(aes(color = active), size = 2.5) +
    facet_wrap(~ shape, ncol = 1) +
    scale_color_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#34495E"),
                       labels = c(`TRUE` = "Active stage", `FALSE` = "Inactive"), name = NULL) +
    scale_x_discrete(labels = nichd_labels) +
    labs(title = "Intended adversarial injection shapes (Module 3)",
         subtitle = "Relative injection weight (1 + shape) per NICHD stage; red = active stages",
         x = "NICHD stage", y = "Relative injection weight") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(output_dir, "m3_fig_adversarial_shapes.png"), p, width = 9, height = 8, dpi = 300, bg = "white")
}

module3 <- function() {
  if (!run_adversarial) { message("\n disabled (run_adversarial = FALSE)."); return(invisible(NULL)) }
  if (!file.exists(ruta_pos_meta)) { message("  metadata not found - skipping."); return(invisible(NULL)) }
  if (!file.exists(ruta_ade_raw)) { message("  ade_raw not found - skipping."); return(invisible(NULL)) }
  message("\n Adversarial non-smooth injection benchmark (refitting)")
  plot_adv_shapes()  # intended adversarial shapes (static reference figure)

  ac_bootstrap_n <- ac_bootstrap_n_posthoc

  # 1. Loads and canonicalize ade_raw 
  ade_raw_dt <- fread(ruta_ade_raw)
  ade_raw_dt[, atc_concept_id := as.character(atc_concept_id)]
  translation_table <- build_drug_translation_table()
  ade_raw_dt <- merge(ade_raw_dt, translation_table[, .(atc_concept_id, canonical_id)],
                      by = "atc_concept_id", all.x = TRUE)
  ade_raw_dt[!is.na(canonical_id), atc_concept_id := canonical_id]
  ade_raw_dt[, canonical_id := NULL]
  ade_raw_dt <- unique(ade_raw_dt, by = c("safetyreportid", "atc_concept_id", "meddra_concept_id"))
  ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
  ade_raw_dt[, nichd_num := as.integer(nichd)]
  keep_cols <- c("safetyreportid", "atc_concept_id", "meddra_concept_id", "nichd", "nichd_num")
  ade_raw_dt <- ade_raw_dt[, ..keep_cols]
  setindex(ade_raw_dt, atc_concept_id)
  setindex(ade_raw_dt, meddra_concept_id)

  # 2. Samples distinct base drug pairs from the positive metadata.
  meta <- fread(ruta_pos_meta)
  meta[, `:=`(drugA = as.character(drugA), drugB = as.character(drugB),
              meddra = as.character(meddra))]
  base_pairs <- unique(meta, by = "base_triplet_id")[, .(drugA, drugB, meddra, fold_change)]
  set.seed(adv_seed)
  n_take <- min(n_adv_pairs, nrow(base_pairs))
  base_pairs <- base_pairs[sample(.N, n_take)]
  base_pairs[, base_pair_id := seq_len(.N)]
  # Stratify the adversarial shapes across base pairs: spike position cycles through all NICHD stages, step is split half rising / half falling.
  base_pairs[, spike_stage := ((base_pair_id - 1L) %% 7L) + 1L]
  base_pairs[, step_rising := base_pair_id %% 2L == 1L]
  message(sprintf(" re-injecting %d base pairs x %d shapes = up to %d fits",
                  n_take, length(adv_shapes), n_take * length(adv_shapes)))

  # 3. Parallel fan-out over base pairs (each worker injects all shape arms).
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  clusterExport(cl, c(
    "process_adv_pair", "get_shape", "build_eval_table", "fit_gam",
    "calculate_classic_ior", "calculate_classic_ac",
    "niveles_nichd", "spline_individuales", "include_sex", "include_stage_sex",
    "k_spline", "bs_type", "select", "nichd_spline", "method",
    "classic_continuity_correction", "continuity_correction_value",
    "ac_bootstrap_n", "adv_shapes", "ade_raw_dt", "adv_seed"
  ), envir = environment())
  clusterEvalQ(cl, { library(data.table); library(mgcv); library(MASS) })

  adv_results <- foreach(
    bp = seq_len(nrow(base_pairs)),
    .packages = c("data.table", "mgcv"),
    .errorhandling = "pass"
  ) %dopar% {
    set.seed(adv_seed + bp)
    process_adv_pair(base_pairs[bp], ade_raw_dt, adv_shapes, get_shape, ac_bootstrap_n)
  }
  stopCluster(cl)

  adv_results <- adv_results[!sapply(adv_results, function(x) inherits(x, "error") || is.null(x))]
  if (length(adv_results) == 0) { message(" no successful fits - skipping analysis."); return(invisible(NULL)) }
  adv <- rbindlist(adv_results, fill = TRUE)
  fwrite(adv, paste0(output_dir, "m3_adversarial_fits.csv"))

  # 4. Injection retention per shape
  retention <- rbindlist(lapply(setdiff(adv_shapes, "none"), function(sh) {
    data.table(shape = sh, n_pairs = n_take, n_success = adv[shape == sh, uniqueN(base_pair_id)])
  }))
  retention[, pct_retained := round(100 * n_success / n_pairs, 1)]
  fwrite(retention, paste0(output_dir, "m3_injection_retention.csv"))
  add_finding("adversarial injection retention", "min_pct_pairs_with_>=1_injected",
              as.character(min(retention$pct_retained)), "higher=less loss")

  # 5. Operating-characteristic metrics (like 30_metrics) under the double threshold
  null_thresholds <- load_null_thresholds()
  double_threshold <- !is.null(null_thresholds)
  if (!double_threshold) {
    null_thresholds <- data.table(stage = 1:7, threshold_ior = NA_real_, threshold_ac = NA_real_,
                                  threshold_classic_ior = NA_real_, threshold_classic_ac = NA_real_)
    message(" null thresholds incomplete - nominal CI (single threshold)")
  }

  adv_thr <- merge(adv, null_thresholds, by.x = "stage_num", by.y = "stage", all.x = TRUE)
  neg_pool <- copy(adv_thr[shape == "none"])[, `:=`(triplet_id = paste0("neg_", base_pair_id), label = 0L)]

  # Metrics 
  compute_metrics_set <- function(pos_shapes, scope_label) {
    pos <- copy(adv_thr[shape %in% pos_shapes & is_active == TRUE])
    if (nrow(pos) == 0) return(NULL)
    pos[, `:=`(triplet_id = paste0("pos_", base_pair_id, "_", shape), label = 1L)]
    rbindlist(lapply(metodos_core, function(met) {
      dt <- detect_signal(rbind(pos, neg_pool, fill = TRUE), met$nombre, met$tipo,
                          use_null = double_threshold)
      m <- calculate_metrics(dt, n_boot_posthoc, aggregate_triplet = TRUE,
                             score_type = met$score_type, score_type_auc = met$score_type_auc)
      m[, `:=`(scope = scope_label, method = met$nombre, measure = met$measure)]
      m
    }), fill = TRUE)
  }

  adversarial_shapes <- setdiff(adv_shapes, "none")
  metrics <- rbindlist(c(
    list(compute_metrics_set(adversarial_shapes, "global")),
    lapply(adversarial_shapes, function(sh) compute_metrics_set(sh, sh))
  ), fill = TRUE)
  setcolorder(metrics, c("scope", "method", "measure"))
  fwrite(metrics, paste0(output_dir, "m3_adversarial_metrics.csv"))

  # 6. Findings
  g <- metrics[scope == "global"]
  for (meas in c("IOR", "AC")) {
    gam_name <- if (meas == "IOR") "GAM-logIOR" else "GAM-AC"
    cls_name <- if (meas == "IOR") "Estratificado-IOR" else "Estratificado-AC"
    gr <- g[method == gam_name]; cr <- g[method == cls_name]
    if (nrow(gr) && nrow(cr)) {
      add_finding(sprintf("adversarial sensitivity GAM vs stratified (%s)", meas),
                  "sensitivity_GAM_minus_stratified",
                  as.character(round(gr$sensitivity - cr$sensitivity, 3)), "+favors GAM")
      add_finding(sprintf("adversarial AUC GAM vs stratified (%s)", meas),
                  "AUC_GAM_minus_stratified",
                  as.character(round(gr$AUC - cr$AUC, 3)), "+favors GAM")
      add_finding(sprintf("adversarial specificity GAM vs stratified (%s)", meas),
                  "specificity_GAM_vs_stratified",
                  sprintf("%.3f vs %.3f", gr$specificity, cr$specificity), "higher=fewer FP")
    }
  }

  # 7. Figure (one per measure)
  metric_labels_adv <- c(AUC = "AUC", sensitivity = "Sensitivity", specificity = "Specificity",
                         PPV = "PPV", NPV = "NPV", F1 = "F1-Score")
  plot_measure <- function(meas, gam_name, cls_name) {
    d <- g[method %in% c(gam_name, cls_name)]
    if (nrow(d) == 0) return(invisible(NULL))
    dl <- rbindlist(lapply(names(metric_labels_adv), function(mm) {
      data.table(
        method_type = ifelse(grepl("GAM", d$method), "GAM", "Stratified"),
        metric = factor(metric_labels_adv[[mm]], levels = unname(metric_labels_adv)),
        value = d[[mm]],
        lower = if (paste0(mm, "_lower") %in% names(d)) d[[paste0(mm, "_lower")]] else NA_real_,
        upper = if (paste0(mm, "_upper") %in% names(d)) d[[paste0(mm, "_upper")]] else NA_real_
      )
    }))
    p <- ggplot(dl, aes(x = method_type, y = value, fill = method_type)) +
      geom_col(width = 0.6) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, na.rm = TRUE) +
      facet_wrap(~ metric, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = c(GAM = "#16A085", Stratified = "#C0392B"), name = "Method") +
      labs(title = sprintf("Adversarial-subset operating characteristics - %s", meas),
           subtitle = "Non-smooth injected signals vs uninjected 'none' arm (double threshold)",
           x = NULL, y = "Metric value")
    ggsave(paste0(output_dir, sprintf("m3_fig_adversarial_metrics_%s.png", tolower(meas))),
           p, width = 11, height = 7, dpi = 300, bg = "white")
  }
  plot_measure("IOR", "GAM-logIOR", "Estratificado-IOR")
  plot_measure("AC", "GAM-AC", "Estratificado-AC")

  gi <- g[method == "GAM-logIOR"]; si <- g[method == "Estratificado-IOR"]
  if (nrow(gi) && nrow(si)) {
    message(sprintf(" adversarial (IOR): GAM sens=%.2f spec=%.2f | stratified sens=%.2f spec=%.2f",
                    gi$sensitivity, gi$specificity, si$sensitivity, si$specificity))
  }
  invisible(NULL)
}

################################################################################
# Null adequacy 
################################################################################
# Tests if the empirical null distribution is an adequate per-stage reference
#
# If real negatives sit systematically above the permutation null, threshold is optimistic 


# Score-column map: each method x measure pairs its permutation-null column (from null_distribution.csv)
# nominal: per-method significance floor for the double criterion (0 on the
# lower-bound scale for both the log-IOR and the additive contrast).
adequacy_methods <- list(
  list(label = "GAM-IOR",         null_col = "log_lower90",          neg_col = "log_ior_lower90",         nominal = 0),
  list(label = "GAM-AC",        null_col = "ac_lower90",         neg_col = "ac_lower90",            nominal = 0),
  list(label = "Stratified-IOR",  null_col = "classic_ior_lower90",  neg_col = "log_ior_classic_lower90", nominal = 0),
  list(label = "Stratified-AC", null_col = "classic_ac_lower90", neg_col = "AC_classic_lower90",    nominal = 0)
)

module4 <- function() {
  if (!file.exists(ruta_null_dist)) { message("\n null distribution not found - skipping."); return(invisible(NULL)) }
  message("\n Null adequacy (per-stage, both measures and methods)")

  null_dist <- fread(ruta_null_dist)

  if (file.exists(ruta_neg_rds)) {
    neg <- readRDS(ruta_neg_rds)
    neg <- neg[reduction_pct == 0 & model_success == TRUE]

    # Unwraps one negative lower-bound list-column to (triplet_id, stage, lower90).
    expand_neg <- function(col) {
      if (!col %in% names(neg)) return(NULL)
      neg[, {
        st <- unlist(stage); vv <- unlist(get(col))
        k <- min(length(st), length(vv))
        if (k > 0) data.table(stage = st[1:k], lower90 = vv[1:k]) else data.table()
      }, by = triplet_id]
    }

    # Per-stage realized FPR on real negatives at the null p95 threshold
    # one row per method x measure x stage. 
    # Under an adequate null this sits near 5%.
    cmp <- rbindlist(lapply(adequacy_methods, function(am) {
      if (!am$null_col %in% names(null_dist)) return(NULL)
      ne <- expand_neg(am$neg_col)
      if (is.null(ne)) return(NULL)
      rbindlist(lapply(1:7, function(s) {
        nv <- null_dist[stage == s][[am$null_col]]; nv <- nv[is.finite(nv)]
        gv <- ne[stage == s & is.finite(lower90), lower90]
        if (length(nv) < 10 || length(gv) < 10) return(NULL)
        thr <- quantile(nv, 0.95, na.rm = TRUE)        # permutation-null threshold (p95)
        # Effective detection threshold = the double criterion lower90 > 0 AND > thr
        eff_thr <- max(am$nominal, thr)
        data.table(
          method = am$label, stage = s,
          null_p95 = round(thr, 3),
          neg_median = round(median(gv), 3),
          # Exceedance rate of real negatives over the double-criterion threshold per stage
          neg_exceedance_rate = round(mean(gv > eff_thr) * 100, 1)
        )
      }), fill = TRUE)
    }), fill = TRUE)

    if (!is.null(cmp) && nrow(cmp) > 0) {
      fwrite(cmp, paste0(output_dir, "m4_null_vs_negatives_by_stage.csv"))
      for (lab in unique(cmp$method)) {
        mean_exc <- round(mean(cmp[method == lab, neg_exceedance_rate], na.rm = TRUE), 1)
        add_finding(sprintf("null adequacy vs real negatives (%s)", lab),
                    "mean_per-stage_neg_exceedance_rate", as.character(mean_exc), ">5 upper-bound FPR")
        message(sprintf("  %-15s mean per-stage negative exceedance = %.1f%%", lab, mean_exc))
      }
    }
  } else {
    message(" negative results not found - reporting null self-consistency only.")
    # Null self-consistency: by construction the per-stage p95 leaves ~5% above it.
    self <- rbindlist(lapply(adequacy_methods, function(am) {
      if (!am$null_col %in% names(null_dist)) return(NULL)
      null_dist[, {
        v <- get(am$null_col)
        .(method = am$label,
          pct_above_own_p95 = round(mean(v > quantile(v, 0.95, na.rm = TRUE), na.rm = TRUE) * 100, 1))
      }, by = stage]
    }), fill = TRUE)
    fwrite(self, paste0(output_dir, "m4_null_self_consistency.csv"))
    for (lab in unique(self$method)) {
      add_finding(sprintf("null self-consistency (%s)", lab), "mean_pct_above_own_p95",
                  as.character(round(mean(self[method == lab, pct_above_own_p95]), 1)), "~5 expected")
    }
  }
  invisible(NULL)
}

################################################################################
# Findings summary
################################################################################

run_summary <- function() {
  summary_dt <- rbindlist(findings, fill = TRUE)
  setcolorder(summary_dt, c("justification", "test", "metric", "value", "guide"))
  fwrite(summary_dt, paste0(output_dir, "posthoc_validation_summary.csv"))

  message("\n post-hoc summary")
  for (i in seq_len(nrow(summary_dt))) {
    message(sprintf("  [%s] %s = %s  (%s)",
                    summary_dt$justification[i], summary_dt$metric[i],
                    summary_dt$value[i], summary_dt$guide[i]))
  }
  invisible(NULL)
}

################################################################################
# Run all tests
################################################################################

module1()
module_shape_fidelity()
module3()
module4()
run_summary()
