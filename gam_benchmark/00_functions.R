################################################################################
# Shared functions for the gam_benchmark
################################################################################

################################################################################
# Configuration
################################################################################

set.seed(12345)

library(pacman)
pacman::p_load(data.table, pbapply, mgcv, MASS)

# Canonical NICHD stage ordering
niveles_nichd <- c(
  "term_neonatal", "infancy", "toddler", "early_childhood",
  "middle_childhood", "early_adolescence", "late_adolescence"
)

# Modeling table. Read live from the upstream producer (faers_parsing)
# The OMOP vocabulary is resolved below in the mapping block.
ruta_ade_raw <- "../faers_parsing/data/processed/ade_raw.csv"

# GAM formula parameters . consumed by fit_gam() and fit_reduced_model().
spline_individuales <- TRUE
include_sex <- FALSE
include_stage_sex <- FALSE
k_spline <- 7
include_nichd <- FALSE
nichd_spline <- FALSE
bs_type <- "cs"
select <- FALSE
method <- "fREML"

# Output filename suffix encoding the active model parametrization.
suffix <- paste0(
  if (spline_individuales) "si" else "",
  if (include_sex) "s" else "",
  if (include_stage_sex) "ss" else "",
  if (include_nichd) "n" else "",
  if (nichd_spline) "ns" else "",
  bs_type
)

Z90 <- qnorm(0.95)  # 90% CI multiplier (one-sided 5th/95th).

# Null-distribution percentile used as the detection threshold.
percentil <- "p95"

# Continuity correction for the stratified (classic) 2x2 estimators.
# When TRUE, a Haldane-Anscombe +0.5 is added to every cell before computing
classic_continuity_correction <- TRUE
continuity_correction_value <- 0.5

################################################################################
# Helper function: compute basic counts
################################################################################

# Computes report-level counts for a drug-event triplet.
# Merges real and simulated events when the augmented flag is present.
# Returns: n_events, n_events_coadmin, n_coadmin.

calc_basic_counts <- function(ade_data, drugA, drugB, meddra) {
  r_a <- unique(ade_data[atc_concept_id == drugA, safetyreportid])
  r_b <- unique(ade_data[atc_concept_id == drugB, safetyreportid])
  r_coadmin <- intersect(r_a, r_b)
  r_ea <- unique(ade_data[meddra_concept_id == meddra, safetyreportid])
  if("simulated_event" %in% names(ade_data)) {
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
# GAM fitting function
################################################################################

# Fits GAM for a drug-drug-event triplet and returns per-stage log-IOR and AC with 90% CIs.
#
# drugA_id, drugB_id, event_id: OMOP concept IDs for the triplet.
# ade_data: report-level data.table (may be augmented with simulated events).
# include_nichd: add NICHD stage as a baseline covariate.
# nichd_spline: model the NICHD baseline as a spline (FALSE = linear).
# spline_individuales: smooth individual drug baseline risks by stage.
# bs_type: spline basis ("cs", "tp", or "cr").
# select: allow penalty shrinkage to zero for each smooth term.
# include_sex, include_stage_sex: sex main effect and stage-by-sex interaction.
# k_spline: basis dimension; 7 matches the number of NICHD stages.
# method: BAM fitting criterion; fREML is faster than ML for large n.
#
# Returns: list with success flag, counts, log_ior/ac vectors (length 7), CIs,
#          and model diagnostics.

fit_gam <- function(drugA_id, drugB_id, event_id, ade_data,
                                 nichd_spline = TRUE,
                                 include_nichd = TRUE,
                                 spline_individuales = FALSE,
                                 bs_type = "cs",
                                 select = FALSE,
                                 include_sex = FALSE,
                                 include_stage_sex = FALSE,
                                 k_spline = 7,
                                 method = "fREML") {
  ###########
  # 1- Identify reports
  ###########
  
  reportes_droga_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reportes_droga_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])
  reportes_ea_real <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])
  reportes_coadmin <- intersect(reportes_droga_a, reportes_droga_b)
  
  # Include reports whose event was injected during signal augmentation.
  reportes_ea_sim <- if("simulated_event" %in% names(ade_data)) {
    unique(ade_data[
      simulated_event == TRUE & simulated_meddra == event_id,
      safetyreportid
    ])
  } else {
    integer(0)
  }
  # Re-query real events after the simulated set is built (safe to be redundant).
  reportes_ea_real <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])

  # Merge real and simulated event reports for a unified outcome set.
  reportes_ea <- union(reportes_ea_real, reportes_ea_sim)
  
  n_events_total <- length(reportes_ea)                         # outcome count regardless of exposure
  n_coadmin <- length(reportes_coadmin)                         # A+B co-administration count
  n_events_coadmin  <- length(intersect(reportes_coadmin, reportes_ea)) # outcome + co-administration count

  ###########
  # 2- Build the modeling dataset
  ###########
  
  # Minimum columns needed for the modeling table.
  cols_necesarias <- c("safetyreportid", "nichd", "nichd_num")

  if (include_sex) {
    cols_necesarias <- c(cols_necesarias, "sex")
  }
  
  datos_modelo <- unique(ade_data[, ..cols_necesarias])
  
  # Binary exposure and outcome indicators per report.
  datos_modelo[, ea_ocurrio := as.integer(safetyreportid %in% reportes_ea)]
  datos_modelo[, droga_a := as.integer(safetyreportid %in% reportes_droga_a)]
  datos_modelo[, droga_b := as.integer(safetyreportid %in% reportes_droga_b)]
  datos_modelo[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]
  
  # Reduce spline complexity when one or more pediatric stages are absent.
  n_unique_stages <- uniqueN(datos_modelo[!is.na(nichd_num), nichd_num])
  k_spline_model <- max(3L, min(k_spline, n_unique_stages - 1L))
  
  if (include_sex) {
    # Normalize single-letter codes to full labels before factoring.
    datos_modelo[, sex := toupper(trimws(sex))]
    datos_modelo[sex == "M", sex := "MALE"]
    datos_modelo[sex == "F", sex := "FEMALE"]
    datos_modelo[, sex := factor(sex, levels = c("MALE", "FEMALE"))]
  }
  
  ###########
  # 4- Build formula from parameters
  ###########
  
  formula_parts <- "ea_ocurrio ~ "

  # Individual drug baseline risk: linear or stage-varying spline.
  if (!spline_individuales) {
    formula_parts <- paste0(formula_parts, "droga_a + droga_b + ")
  } else {
    formula_parts <- paste0(
      formula_parts,
      sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_a) + ", 
              k_spline_model, bs_type),
      sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_b) + ", 
              k_spline_model, bs_type)
    )
  }

  # NICHD baseline: spline or linear, depending on parametrization.
  if (include_nichd) {
    if (nichd_spline) {
      formula_parts <- paste0(
        formula_parts,
        sprintf("s(nichd_num, k = %d, bs = '%s') + ", k_spline_model, bs_type)
      )
    } else {
      formula_parts <- paste0(formula_parts, "nichd_num + ")
    }
  }

  # Interaction spline (do not modify this term)
  formula_parts <- paste0(
    formula_parts,
    sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_ab)", k_spline_model, bs_type)
  )
  
  if (include_sex) {
    if (include_stage_sex) {
      # Stage-varying sex effect via a by-sex spline.
      formula_parts <- paste0(
        formula_parts,
        sprintf(" + s(nichd_num, k = %d, bs = '%s', by = sex)",
                k_spline_model, bs_type)
      )
    } else {
      formula_parts <- paste0(formula_parts, " + sex")
    }
  }

  formula_final <- as.formula(formula_parts)
  
  ###########
  # 5- Model fitting
  ###########
  
  tryCatch({
    
    modelo <- bam(
      formula = formula_final,
      data = datos_modelo,
      family = binomial(link = "logit"),
      method = method,
      select = select,
      discrete = TRUE,
      nthreads = 1  # Single thread avoids conflicts with the outer pblapply in 10_augmentation.
    )
    
    ###########
    # 6- Compute log-IOR per NICHD stage
    ###########
    
    # Full factorial prediction grid over stages and exposure groups.
    grid_dif <- CJ(
      nichd_num = 1:7, 
      droga_a = c(0, 1), 
      droga_b = c(0, 1)
    )
    grid_dif[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]
    
    if (include_sex) {
      # Hold sex at reference level (MALE) to isolate the drug-interaction contrast.
      grid_dif[, sex := factor("MALE", levels = c("MALE", "FEMALE"))]
    }

    # Predictions on the link (log-odds) scale.
    pred_dif <- predict(modelo, newdata = grid_dif, type = "link", se.fit = TRUE)
    grid_dif[, `:=`(lp = pred_dif$fit, se = pred_dif$se.fit)]
    
    # Wide pivot so each row holds the four link-scale predictions for one stage.
    w_lp <- dcast(grid_dif, nichd_num ~ droga_a + droga_b,
                  value.var = c("lp", "se"))

    # log(IOR) = lp_11 - lp_10 - lp_01 + lp_00 (additive contrast on the log-odds scale).
    log_ior <- w_lp$lp_1_1 - w_lp$lp_1_0 - w_lp$lp_0_1 + w_lp$lp_0_0
    
    ###########
    # 7- Compute log-IOR standard error using the covariance matrix
    ###########
    
    Xp <- predict(modelo, newdata = grid_dif, type = "lpmatrix")
    Vb <- vcov(modelo, unconditional = TRUE) # unconditional = TRUE includes smoothing-parameter uncertainty in the SE
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
      
      # IOR contrast vector: c'Vc gives the variance of the log-IOR linear combination.
      cvec <- rep(0, nrow(grid_dif))
      cvec[c(idx_11, idx_10, idx_01, idx_00)] <- c(1, -1, -1, 1)

      # SE(log-IOR) = sqrt(c' V c); clamp at 0 to guard against tiny numerical negatives.
      log_ior_se[stage] <- sqrt(max(
        as.numeric(t(cvec) %*% cov_link %*% cvec), 
        0
      ))
    }
    
    ###########
    # 8- Compute confidence intervals and metrics
    ###########
    
    # Fallback if Z90 was not initialised at the top of this script.
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
    # 9- Compute AC per stage with 90% CI
    ###########
    
    stages <- sort(unique(datos_modelo$nichd_num))

    # Four exposure combinations per stage for the AC prediction grid.
    nd_ac <- rbindlist(lapply(stages, function(s) {
      data.table(
        nichd_num = s,
        droga_a = c(0, 1, 0, 1),
        droga_b = c(0, 0, 1, 1),
        droga_ab = c(0, 0, 0, 1)
      )
    }), use.names = TRUE)
    
    # Add any extra covariates required by the formula, held at reference levels.
    if (include_sex) {
      nd_ac[, sex := factor(levels(datos_modelo$sex)[1],
                              levels = levels(datos_modelo$sex))]
    }
    if (include_nichd && !nichd_spline) {
      nd_ac[, nichd := factor(niveles_nichd[nichd_num],
                                levels = niveles_nichd,
                                ordered = TRUE)]
    }

    pred_ac <- predict(modelo, newdata = nd_ac, type = "link", se.fit = TRUE)
    nd_ac[, `:=`(
      eta = pred_ac$fit,
      se  = pred_ac$se.fit
    )]
    
    # Parametric bootstrap for AC CIs
    X_ac <- predict(modelo, newdata = nd_ac, type = "lpmatrix")
    beta_hat <- coef(modelo)
    V_beta <- vcov(modelo, unconditional = TRUE)

    B <- 2000  # bootstrap draws; matches the n_boot default in calculate_benchmark_metrics.

    # Draw beta from its asymptotic joint MVN distribution.
    beta_sims <- mvrnorm(n = B, mu = beta_hat, Sigma = V_beta)

    p_sims <- plogis(X_ac %*% t(beta_sims))

    # Additive interaction contrast (Thakrar) on predicted reporting proportions:
    # p11 - p10 - p01 + p00. 
    # detection uses the studentized contrast ac_z = contrast / SE. scale-free.
    calc_add <- function(p) {
      p11 <- p[4]; p10 <- p[2]; p01 <- p[3]; p00 <- p[1]
      p11 - p10 - p01 + p00
    }

    # Per-stage summary from the bootstrap distribution.
    ac_dt <- nd_ac[, {
      idx <- .I
      p_mat <- p_sims[idx, , drop = FALSE]
      add_sim  <- apply(p_mat, 2, calc_add)
      add_mean <- mean(add_sim)
      add_se   <- sd(add_sim)

      data.table(
        AC = add_mean,
        AC_lower90 = quantile(add_sim, 0.05),
        AC_upper90 = quantile(add_sim, 0.95),
        ac_z = add_mean / add_se
      )
    }, by = nichd_num]

    ac_values <- ac_dt$AC
    ac_lower90 <- ac_dt$AC_lower90
    ac_upper90 <- ac_dt$AC_upper90
    ac_z       <- ac_dt$ac_z
    
    ###########
    # 10- Results
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
      ac_z = ac_z,
      n_stages_ac_significant = sum(ac_z > qnorm(0.95), na.rm = TRUE),
      model_aic = AIC(modelo),
      model_deviance = deviance(modelo),
      formula_used = formula_parts,  # Stored for traceability in the result table.
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

# Computes the classical (stratified) IOR via 2x2 contingency tables per NICHD stage.
# IOR = (OR_11 x OR_00) / (OR_10 x OR_01), where OR_00 = 1 by definition.
# 90% CI via the Woolf method on the log scale.
# Returns: list(success, results_by_stage) with one row per stage.

calculate_classic_ior <- function(drugA_id, drugB_id, event_id, ade_data) {
  
  reportes_droga_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reportes_droga_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])

  reportes_ea_real <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])

  # Include injected events when the augmented flag is present.
  reportes_ea_sim <- if("simulated_event" %in% names(ade_data)) {
    unique(ade_data[
      simulated_event == TRUE & simulated_meddra == event_id,
      safetyreportid
    ])
  } else {
    integer(0)
  }

  reportes_ea <- union(reportes_ea_real, reportes_ea_sim)

  # One row per report with binary exposure and outcome indicators.
  datos_unicos <- unique(ade_data[, .(safetyreportid, nichd, nichd_num)])
  datos_unicos[, ea_ocurrio := as.integer(safetyreportid %in% reportes_ea)]
  datos_unicos[, droga_a := as.integer(safetyreportid %in% reportes_droga_a)]
  datos_unicos[, droga_b := as.integer(safetyreportid %in% reportes_droga_b)]

  # Haldane-Anscombe +0.5 per cell when correction is enabled; cc=0 gives raw textbook form.
  cc <- if (classic_continuity_correction) continuity_correction_value else 0

  stage_results <- datos_unicos[, {

    # Cell counts for each exposure stratum (event / no-event); raw, kept for diagnostics.
    n_11_evento    <- sum(droga_a == 1 & droga_b == 1 & ea_ocurrio == 1)
    n_11_no_evento <- sum(droga_a == 1 & droga_b == 1 & ea_ocurrio == 0)
    n_10_evento    <- sum(droga_a == 1 & droga_b == 0 & ea_ocurrio == 1)
    n_10_no_evento <- sum(droga_a == 1 & droga_b == 0 & ea_ocurrio == 0)
    n_01_evento    <- sum(droga_a == 0 & droga_b == 1 & ea_ocurrio == 1)
    n_01_no_evento <- sum(droga_a == 0 & droga_b == 1 & ea_ocurrio == 0)
    n_00_evento    <- sum(droga_a == 0 & droga_b == 0 & ea_ocurrio == 1)
    n_00_no_evento <- sum(droga_a == 0 & droga_b == 0 & ea_ocurrio == 0)

    # Corrected counts for estimation; raw counts kept for diagnostic columns.
    a11 <- n_11_evento + cc; b11 <- n_11_no_evento + cc
    a10 <- n_10_evento + cc; b10 <- n_10_no_evento + cc
    a01 <- n_01_evento + cc; b01 <- n_01_no_evento + cc
    a00 <- n_00_evento + cc; b00 <- n_00_no_evento + cc

    # ORs relative to the referent group (00); OR_00 = 1 by definition.
    or_11 <- (a11 / b11) / (a00 / b00)
    or_10 <- (a10 / b10) / (a00 / b00)
    or_01 <- (a01 / b01) / (a00 / b00)
    or_00 <- 1  # reference

    # IOR = OR_11 / (OR_10 * OR_01) since OR_00 = 1.
    ior_val <- (or_11 * or_00) / (or_10 * or_01)
    log_ior <- log(ior_val)

    # Woolf variance: sum of reciprocal cell counts on the log scale, on the corrected cells.
    var_log_ior <- (1/a11 + 1/b11 +
                    1/a10 + 1/b10 +
                    1/a01 + 1/b01 +
                    1/a00 + 1/b00)
    se_log_ior <- sqrt(var_log_ior)

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
      # Diagnostic counts for quality checks.
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

# Computes the stratified additive interaction contrast (Thakrar) per NICHD stage
# from 2x2 tables: R11 - R10 - R01 + R00 on the reporting proportions. 
# Detection uses the studentized contrast AC_classic_z = contrast / SE, with SE 
# Returns NA for a stage when any exposure stratum is empty.

calculate_classic_ac <- function(drugA_id, drugB_id, event_id, ade_data) {
  
  reportes_droga_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reportes_droga_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])

  reportes_ea_real <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])

  # Include injected events when the augmented flag is present.
  reportes_ea_sim <- if("simulated_event" %in% names(ade_data)) {
    unique(ade_data[
      simulated_event == TRUE & simulated_meddra == event_id,
      safetyreportid
    ])
  } else {
    integer(0)
  }

  reportes_ea <- union(reportes_ea_real, reportes_ea_sim)

  # One row per report with binary exposure and outcome indicators.
  datos_unicos <- unique(ade_data[, .(safetyreportid, nichd, nichd_num)])
  datos_unicos[, ea_ocurrio := as.integer(safetyreportid %in% reportes_ea)]
  datos_unicos[, droga_a := as.integer(safetyreportid %in% reportes_droga_a)]
  datos_unicos[, droga_b := as.integer(safetyreportid %in% reportes_droga_b)]

  # Haldane-Anscombe +0.5 per group when correction is enabled; cc=0 gives raw form.
  cc <- if (classic_continuity_correction) continuity_correction_value else 0

  stage_results <- datos_unicos[, {

    # Event and total counts per exposure stratum (raw, kept for diagnostics).
    n_11_evento <- sum(droga_a == 1 & droga_b == 1 & ea_ocurrio == 1)
    n_11_total  <- sum(droga_a == 1 & droga_b == 1)
    n_10_evento <- sum(droga_a == 1 & droga_b == 0 & ea_ocurrio == 1)
    n_10_total  <- sum(droga_a == 1 & droga_b == 0)
    n_01_evento <- sum(droga_a == 0 & droga_b == 1 & ea_ocurrio == 1)
    n_01_total  <- sum(droga_a == 0 & droga_b == 1)
    n_00_evento <- sum(droga_a == 0 & droga_b == 0 & ea_ocurrio == 1)
    n_00_total  <- sum(droga_a == 0 & droga_b == 0)

    # Corrected denominators (cc = 0 reproduces the raw proportions).
    d11 <- n_11_total + 2 * cc
    d10 <- n_10_total + 2 * cc
    d01 <- n_01_total + 2 * cc
    d00 <- n_00_total + 2 * cc

    # Without correction (cc = 0) an empty exposure group leaves the risk undefined.
    if (d11 == 0 || d10 == 0 || d01 == 0 || d00 == 0) {
      data.table(
        stage = nichd_num[1],
        AC_classic = NA_real_,
        AC_classic_lower90 = NA_real_,
        AC_classic_upper90 = NA_real_,
        AC_classic_se = NA_real_,
        AC_classic_z = NA_real_,
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
  }
    else {  # cc > 0 guarantees all denominators are positive
      R11 <- (n_11_evento + cc) / d11
      R10 <- (n_10_evento + cc) / d10
      R01 <- (n_01_evento + cc) / d01
      R00 <- (n_00_evento + cc) / d00

      # Additive interaction contrast (Thakrar): no division by the reference R00.
      ac_val <- R11 - R10 - R01 + R00

      # Binomial variance per proportion; boundary correction (0.25/n) at r = 0 or 1.
      var_r <- function(r, n) ifelse(r > 0 & r < 1, r*(1-r)/n, 0.25/n)
      var_R11 <- var_r(R11, d11)
      var_R10 <- var_r(R10, d10)
      var_R01 <- var_r(R01, d01)
      var_R00 <- var_r(R00, d00)

      # Var(contrast) = sum of the four binomial variances (independent groups);
      # the fragile R00^4 term of the AC ratio is gone.
      se_ac <- sqrt(var_R11 + var_R10 + var_R01 + var_R00)

      # Studentized contrast is the scale-free detection score; 90% CI for display.
      z90 <- qnorm(0.95)
      ac_z <- ac_val / se_ac
      ac_lower90 <- ac_val - z90 * se_ac
      ac_upper90 <- ac_val + z90 * se_ac

      data.table(
        stage = nichd_num[1],
        AC_classic = ac_val,
        AC_classic_lower90 = ac_lower90,
        AC_classic_upper90 = ac_upper90,
        AC_classic_se = se_ac,
        AC_classic_z = ac_z,
        # Individual proportions and counts retained for quality checks.
        R11 = R11, R10 = R10, R01 = R01, R00 = R00,
        n_11_evento = n_11_evento, n_11_total = n_11_total,
        n_10_evento = n_10_evento, n_10_total = n_10_total,
        n_01_evento = n_01_evento, n_01_total = n_01_total,
        n_00_evento = n_00_evento, n_00_total = n_00_total,
        insufficient_data = FALSE
      )
    }
  }, by = nichd_num]
  setorder(stage_results, nichd_num)

  # Return failure when every stage is NA (no usable data across the full age range).
  if (all(is.na(stage_results$AC_classic))) {
    return(list(
      success = FALSE,
      message = "Datos insuficientes en todas las etapas",
      results_by_stage = stage_results
    ))
  }
  
  return(list(
    success = TRUE,
    results_by_stage = stage_results
  ))
}

################################################################################
# Model fitting function on a reduced dataset
################################################################################

# Fits GAM and classical IOR/AC on one (possibly subset) dataset and packs
# all results into a single data.table row per triplet. Used by the reduction
# experiment loop to test sensitivity to sample-size changes.

fit_reduced_model <- function(ade_reduced, rowt, reduction_pct) {
  
  counts_reduced <- calc_basic_counts(ade_reduced, rowt$drugA, rowt$drugB, rowt$meddra)
  
  model_res <- tryCatch({
    fit_gam(
      drugA_id = rowt$drugA,
      drugB_id = rowt$drugB,
      event_id = rowt$meddra,
      ade_data = ade_reduced,
      spline_individuales = spline_individuales,
      include_sex = include_sex,
      include_stage_sex = include_stage_sex,
      k_spline = k_spline,
      bs_type = bs_type,
      select = select,
      nichd_spline = nichd_spline
    )
  }, error = function(e) {
    list(
      success = FALSE,
      n_events_total = counts_reduced$n_events,
      n_coadmin = counts_reduced$n_coadmin,
      error_msg = paste("Reduced-model error:", e$message)
    )
  })
  
  classic_res <- tryCatch({
    calculate_classic_ior(
      drugA_id = rowt$drugA,
      drugB_id = rowt$drugB,
      event_id = rowt$meddra,
      ade_data = ade_reduced
    )
  }, error = function(e) {
    list(success = FALSE)
  })
  
  classic_ac <- tryCatch({
    calculate_classic_ac(
      drugA_id = rowt$drugA,
      drugB_id = rowt$drugB,
      event_id = rowt$meddra,
      ade_data = ade_reduced
    )
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
      ac_z = list(rep(NA_real_, 7)),
      n_stages_ac_significant = NA_integer_,
      AC_classic = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic) else list(rep(NA_real_, 7)),
      AC_classic_lower90 = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic_lower90) else list(rep(NA_real_, 7)),
      AC_classic_upper90 = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic_upper90) else list(rep(NA_real_, 7)),
      AC_classic_se = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic_se) else list(rep(NA_real_, 7)),
      AC_classic_z = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic_z) else list(rep(NA_real_, 7)),
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
    ac_z = list(model_res$ac_z),
    n_stages_ac_significant = model_res$n_stages_ac_significant,
    AC_classic = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic) else list(rep(NA_real_, 7)),
    AC_classic_lower90 = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic_lower90) else list(rep(NA_real_, 7)),
    AC_classic_upper90 = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic_upper90) else list(rep(NA_real_, 7)),
    AC_classic_se = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic_se) else list(rep(NA_real_, 7)),
    AC_classic_z = if(classic_ac$success) list(classic_ac$results_by_stage$AC_classic_z) else list(rep(NA_real_, 7)),
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
# Expansion function
################################################################################

# Expands triplet results to long format with all metrics
#
# Parameters:
# dt: data.table with triplet results (with list-type columns)
# label_val: classification label (1 = positive, 0 = negative)
# null_thresholds_dt: data.table with null distribution thresholds per stage
# use_threshold_ior: if TRUE, applies the IOR null distribution threshold
# use_threshold_ac: if TRUE, applies the AC null distribution threshold
# 
# Return:
# Expanded data.table with one row per triplet-stage
# 
# Implementation:
# Unpacks list-type columns (stage, log_ior, etc.) into individual rows
# Merges with null distribution thresholds
# Retains triplet-level metadata (dynamic, t_ij, n_coadmin)

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
    
    # Unpack list-columns to vectors (stored as lists to allow variable stage lengths).
    gam_log_ior         <- unlist(log_ior)
    gam_log_ior_lower90 <- unlist(log_ior_lower90)
    gam_ac            <- unlist(ac_values)
    gam_ac_lower90    <- unlist(ac_lower90)
    gam_ac_z          <- unlist(ac_z)
    cls_log_ior         <- unlist(log_ior_classic)
    cls_log_ior_lower90 <- unlist(log_ior_classic_lower90)
    cls_ac            <- unlist(AC_classic)
    cls_ac_lower90    <- unlist(AC_classic_lower90)
    cls_ac_z          <- unlist(AC_classic_z)

    n <- min(length(stages), length(gam_log_ior), length(gam_log_ior_lower90),
             length(gam_ac), length(gam_ac_lower90), length(gam_ac_z),
             length(cls_log_ior), length(cls_log_ior_lower90),
             length(cls_ac), length(cls_ac_lower90), length(cls_ac_z))

    if (n > 0) {
      data.table(
        stage_num = stages[1:n],
        # GAM
        gam_log_ior = gam_log_ior[1:n],
        gam_log_ior_lower90 = gam_log_ior_lower90[1:n],
        gam_ac = gam_ac[1:n],
        gam_ac_lower90 = gam_ac_lower90[1:n],
        gam_ac_z = gam_ac_z[1:n],
        # Stratified
        classic_log_ior = cls_log_ior[1:n],
        classic_log_ior_lower90 = cls_log_ior_lower90[1:n],
        classic_ac = cls_ac[1:n],
        classic_ac_lower90 = cls_ac_lower90[1:n],
        classic_ac_z = cls_ac_z[1:n]
      )
    }
  }, by = by_cols]
  
  if (!has_dynamic) {
    expanded[, `:=`(dynamic = "control", t_ij = 0)]
  }
  
  expanded[, nichd := niveles_nichd[stage_num]]
  expanded[, label := label_val]
  
  expanded <- merge(expanded, null_thresholds_dt,
                   by.x = "stage_num", by.y = "stage", all.x = TRUE)

  # Record which threshold mode was active (for downstream interpretation).
  expanded[, `:=`(
    use_threshold_ior = use_threshold_ior,
    use_threshold_ac = use_threshold_ac
  )]
  return(expanded)
}

################################################################################
# Signal detection function
################################################################################

# Adds a `detected` flag and a continuous `signal_score` to dt in place.
# method_name: string; "GAM" substring selects the GAM lower-90 columns.
# detection_type: "IOR", "AC", or any other string (treated as IOR-OR-AC).
# use_null: if TRUE, an additional null-distribution threshold must also be exceeded.

detect_signal <- function(dt, method_name, detection_type, use_null) {

  is_gam <- grepl("GAM", method_name)

  # Select the appropriate lower-90 and null-threshold column names.
  if (is_gam) {
    ior_col        <- "gam_log_ior_lower90"
    ac_col       <- "gam_ac_z"          # studentized additive contrast
    thresh_ior_col <- "threshold_ior"
    thresh_ac_col <- "threshold_ac"
  } else {
    ior_col        <- "classic_log_ior_lower90"
    ac_col       <- "classic_ac_z"       # studentized additive contrast
    # Classic null thresholds are populated when 20_null.R runs the classic methods.
    thresh_ior_col  <- "threshold_classic_ior"
    thresh_ac_col <- "threshold_classic_ac"
  }

  # Nominal cut for the studentized additive contrast: z > qnorm(0.95) is the exact
  # equivalent of the former "90% CI lower bound > 0" test.
  z90 <- qnorm(0.95)

  if (detection_type == "IOR") {
    if (use_null && thresh_ior_col %in% names(dt)) {
      dt[, detected := !is.na(get(ior_col)) & get(ior_col) > 0 & get(ior_col) > get(thresh_ior_col)]
    } else {
      dt[, detected := !is.na(get(ior_col)) & get(ior_col) > 0]}
  } else if (detection_type == "AC") {
    if (use_null && thresh_ac_col %in% names(dt)) {
      dt[, detected := !is.na(get(ac_col)) & get(ac_col) > z90 & get(ac_col) > get(thresh_ac_col)]
    } else {
      dt[, detected := !is.na(get(ac_col)) & get(ac_col) > z90]}
  } else {
    # Any other string triggers the double criterion (IOR OR AC).
    ior_det <- if (use_null && thresh_ior_col %in% names(dt)) {
      !is.na(dt[[ior_col]]) & dt[[ior_col]] > 0 & dt[[ior_col]] > dt[[thresh_ior_col]]
    } else {
      !is.na(dt[[ior_col]]) & dt[[ior_col]] > 0}
    ac_det <- if (use_null && thresh_ac_col %in% names(dt)) {
      !is.na(dt[[ac_col]]) & dt[[ac_col]] > z90 & dt[[ac_col]] > dt[[thresh_ac_col]]
    } else {
      !is.na(dt[[ac_col]]) & dt[[ac_col]] > z90}
    dt[, detected := ior_det | ac_det]
  }

  dt[is.na(detected), detected := FALSE]

  # Continuous score for AUC: the studentized contrast (AC) or the log-IOR lower bound (IOR)
  # ranking is independent of the operating threshold. NA ranks lowest.
  score_col <- if (detection_type == "IOR") ior_col else if (detection_type == "AC") ac_col else NA_character_
  if (!is.na(score_col) && score_col %in% names(dt)) {
    dt[, signal_score := as.numeric(get(score_col))]
  } else {
    dt[, signal_score := NA_real_]
  }

  return(dt)
}

################################################################################
# Null Thresholds
################################################################################

load_modeling_null_thresholds <- function(
  percentil_sel = percentil,
  suffix_sel = suffix,
  require_classic = TRUE
) {
  
  thresh_col <- paste0("threshold_", percentil_sel)
  # Null thresholds are produced by gam_validation/20_null.R
  base_dir <- paste0("../gam_validation/results/", suffix_sel, "/null_distribution_results/")
  
  ruta_ior <- paste0(base_dir, "null_thresholds.csv")
  ruta_ac <- paste0(base_dir, "null_thresholds_ac.csv")
  ruta_cls_ior <- paste0(base_dir, "null_thresholds_classic_ior.csv")
  ruta_cls_ac <- paste0(base_dir, "null_thresholds_classic_ac.csv")
  
  if (!file.exists(ruta_ior) || !file.exists(ruta_ac)) {
    stop("Required GAM null thresholds were not found")
  }
  
  if (require_classic && (!file.exists(ruta_cls_ior) || !file.exists(ruta_cls_ac))) {
    stop("Required stratified null thresholds were not found")
  }
  
  null_thresholds_ior <- fread(ruta_ior)[, .(stage, threshold_ior = get(thresh_col))]
  null_thresholds_ac <- fread(ruta_ac)[, .(stage, threshold_ac = get(thresh_col))]
  
  if (file.exists(ruta_cls_ior) && file.exists(ruta_cls_ac)) {
    null_thresholds_cls_ior <- fread(ruta_cls_ior)[, .(stage, threshold_classic_ior = get(thresh_col))]
    null_thresholds_cls_ac <- fread(ruta_cls_ac)[, .(stage, threshold_classic_ac = get(thresh_col))]
  } else {
    null_thresholds_cls_ior <- data.table(stage = 1:7, threshold_classic_ior = NA_real_)
    null_thresholds_cls_ac <- data.table(stage = 1:7, threshold_classic_ac = NA_real_)
  }
  
  null_thresholds_dt <- merge(null_thresholds_ior, null_thresholds_ac, by = "stage")
  null_thresholds_dt <- merge(null_thresholds_dt, null_thresholds_cls_ior, by = "stage")
  null_thresholds_dt <- merge(null_thresholds_dt, null_thresholds_cls_ac, by = "stage")
  
  return(null_thresholds_dt)
}

################################################################################
# Benchmark fit
################################################################################

# Fits all successfully mapped benchmark triplets

fit_benchmark_triplets <- function(
  benchmark_ready_dt,
  ade_data,
  cache_file = NULL,
  overwrite_cache = FALSE
) {
  
  if (!is.null(cache_file) && file.exists(cache_file) && !overwrite_cache) {
    return(readRDS(cache_file))
  }
  
  benchmark_fit_dt <- copy(benchmark_ready_dt[mapping_success == TRUE])
  
  fit_list <- pblapply(seq_len(nrow(benchmark_fit_dt)), function(i) {
    rowt <- benchmark_fit_dt[i, .(triplet_id, drugA, drugB, meddra, type)]
    fit_reduced_model(
      ade_reduced = ade_data,
      rowt = rowt,
      reduction_pct = 0
    )
  })
  
  fit_results <- rbindlist(fit_list, fill = TRUE)
  fit_results <- merge(
    benchmark_ready_dt,
    fit_results,
    by = c("triplet_id", "type"),
    all.x = TRUE
  )
  
  if (!is.null(cache_file)) {
    saveRDS(fit_results, cache_file)
  }
  
  return(fit_results)
}

################################################################################
# Metrics
################################################################################

expand_benchmark_results <- function(benchmark_results_dt, null_thresholds_dt) {
  
  expanded_dt <- expand_clean_all_metrics(
    dt = benchmark_results_dt,
    label_val = 1,
    null_thresholds_dt = null_thresholds_dt,
    use_threshold_ior = TRUE,
    use_threshold_ac = TRUE
  )
  
  meta_cols <- intersect(
    c(
      "triplet_id", "type", "control_type", "model_success", "classic_success",
      "drug1_original", "drug2_original", "event_original",
      "meddra_preferred_term", "confidence_level", "source_title",
      "source_year", "drugA", "drugB", "meddra",
      # Names and ontogeny fields from the curated reference set
      "drug1_name", "drug2_name", "meddra_pt",
      "ontogenic_modulation", "higher_risk_stages"
    ),
    names(benchmark_results_dt)
  )
  
  expanded_dt <- merge(
    expanded_dt,
    unique(benchmark_results_dt[, ..meta_cols]),
    by = "triplet_id",
    all.x = TRUE
  )
  
  return(expanded_dt)
}

# Rank-based AUC (Mann-Whitney U) of a continuous score separating positive from negative controls
compute_auc <- function(scores, is_positive) {
  y <- as.logical(is_positive)
  s <- as.numeric(scores)
  if (length(unique(y[!is.na(y)])) < 2L) return(NA_real_)
  s[is.na(s)] <- -Inf
  r <- rank(s, ties.method = "average")
  n_pos <- sum(y, na.rm = TRUE)
  n_neg <- sum(!y, na.rm = TRUE)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  (sum(r[y]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

# Summarises benchmark performance from triplet-level detection.
calculate_benchmark_metrics_from_triplet_detail <- function(
  triplet_detail_dt,
  method_name,
  detection_type,
  use_null,
  n_boot = 2000
) {

  triplet_eval <- copy(triplet_detail_dt)

  if (!all(c("evaluable", "detected", "n_stages_detected") %in% names(triplet_eval))) {
    stop("triplet_detail_dt debe contener evaluable, detected y n_stages_detected")
  }
  # Backward compatible: a detail table without control_type is all positives.
  if (!"control_type" %in% names(triplet_eval)) triplet_eval[, control_type := "positive"]
  if (!"signal_score" %in% names(triplet_eval)) triplet_eval[, signal_score := NA_real_]

  is_pos <- triplet_eval$control_type == "positive"
  is_neg <- triplet_eval$control_type == "negative"
  detected <- triplet_eval$detected %in% TRUE
  evaluable <- triplet_eval$evaluable %in% TRUE

  n_total <- nrow(triplet_eval)
  n_positive <- sum(is_pos)
  n_negative <- sum(is_neg)
  n_evaluable <- sum(evaluable)
  has_neg <- n_negative > 0

  safe_div <- function(num, den) if (den > 0) num / den else NA_real_

  # Confusion matrix at the configured operating point (binary detected).
  tp <- sum(is_pos & detected)
  fn <- sum(is_pos & !detected)
  fp <- sum(is_neg & detected)
  tn <- sum(is_neg & !detected)

  # Sensitivity over positive controls (denominator = positives)
  # it is unchanged when the set has no negatives.
  sensitivity_total <- safe_div(tp, n_positive)
  eval_pos <- is_pos & evaluable
  sensitivity_evaluable <- safe_div(sum(eval_pos & detected), sum(eval_pos))

  specificity <- if (has_neg) safe_div(tn, n_negative) else NA_real_
  ppv <- if (has_neg) safe_div(tp, tp + fp) else NA_real_
  npv <- if (has_neg) safe_div(tn, tn + fn) else NA_real_
  f1 <- if (has_neg) safe_div(2 * tp, 2 * tp + fp + fn) else NA_real_
  accuracy <- if (has_neg) safe_div(tp + tn, n_total) else NA_real_
  balanced_accuracy <- if (has_neg && !is.na(sensitivity_total) && !is.na(specificity)) (sensitivity_total + specificity) / 2 else NA_real_
  youden_j <- if (has_neg && !is.na(sensitivity_total) && !is.na(specificity)) sensitivity_total + specificity - 1 else NA_real_
  auc <- if (has_neg) compute_auc(triplet_eval$signal_score, is_pos) else NA_real_

  # Bootstrap CIs over triplets (sensitivity, specificity, PPV, AUC).
  boot_fun <- function() {
    idx <- sample.int(n_total, n_total, replace = TRUE)
    b_pos <- is_pos[idx]; b_neg <- is_neg[idx]
    b_det <- detected[idx]; b_eval <- evaluable[idx]
    b_tp <- sum(b_pos & b_det); b_fp <- sum(b_neg & b_det); b_tn <- sum(b_neg & !b_det)
    b_eval_pos <- b_pos & b_eval
    c(
      safe_div(b_tp, sum(b_pos)),
      safe_div(sum(b_eval_pos & b_det), sum(b_eval_pos)),
      if (sum(b_neg) > 0) safe_div(b_tn, sum(b_neg)) else NA_real_,
      if (sum(b_neg) > 0) safe_div(b_tp, b_tp + b_fp) else NA_real_,
      if (sum(b_neg) > 0) compute_auc(triplet_eval$signal_score[idx], b_pos) else NA_real_
    )
  }
  boot_mat <- replicate(n_boot, boot_fun())
  ci <- function(row) {
    vals <- boot_mat[row, ]
    if (all(is.na(vals))) c(NA_real_, NA_real_) else quantile(vals, c(0.025, 0.975), na.rm = TRUE)
  }
  sens_total_ci <- ci(1); sens_eval_ci <- ci(2); spec_ci <- ci(3); ppv_ci <- ci(4); auc_ci <- ci(5)

  pos_stages <- triplet_eval$n_stages_detected[is_pos]

  data.table(
    method = method_name,
    detection_type = detection_type,
    threshold_mode = if (use_null) paste0("null_", percentil) else "nominal",
    n_benchmark_triplets = n_total,
    n_positive_controls = n_positive,
    n_negative_controls = n_negative,
    n_evaluable_triplets = n_evaluable,
    n_detected_triplets = tp,
    n_true_positive = tp,
    n_false_negative = fn,
    n_false_positive = fp,
    n_true_negative = tn,
    sensitivity_total = sensitivity_total,
    sensitivity_total_lower = sens_total_ci[[1]],
    sensitivity_total_upper = sens_total_ci[[2]],
    sensitivity_evaluable = sensitivity_evaluable,
    sensitivity_evaluable_lower = sens_eval_ci[[1]],
    sensitivity_evaluable_upper = sens_eval_ci[[2]],
    specificity = specificity,
    specificity_lower = spec_ci[[1]],
    specificity_upper = spec_ci[[2]],
    ppv = ppv,
    ppv_lower = ppv_ci[[1]],
    ppv_upper = ppv_ci[[2]],
    npv = npv,
    f1 = f1,
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    youden_j = youden_j,
    auc = auc,
    auc_lower = auc_ci[[1]],
    auc_upper = auc_ci[[2]],
    median_stages_detected = if (n_positive > 0) median(pos_stages, na.rm = TRUE) else NA_real_,
    max_stages_detected = if (n_positive > 0) max(pos_stages, na.rm = TRUE) else NA_real_
  )
}

calculate_benchmark_metrics <- function(
  dt,
  success_col,
  method_name,
  detection_type,
  use_null,
  n_boot = 2000
) {

  triplet_eval <- dt[, .(
    control_type = control_type[1L],
    evaluable = any(get(success_col), na.rm = TRUE),
    detected = any(detected, na.rm = TRUE),
    n_stages_detected = sum(detected, na.rm = TRUE),
    signal_score = if (all(is.na(signal_score))) NA_real_ else max(signal_score, na.rm = TRUE)
  ), by = triplet_id]

  calculate_benchmark_metrics_from_triplet_detail(
    triplet_detail_dt = triplet_eval,
    method_name = method_name,
    detection_type = detection_type,
    use_null = use_null,
    n_boot = n_boot
  )
}

################################################################################
# Benchmark evaluation
################################################################################

evaluate_benchmark_methods <- function(
  benchmark_expanded_dt,
  n_boot = 2000,
  methods_cfg = NULL
) {
  
  if (is.null(methods_cfg)) {
    methods_cfg <- list(
      list(method = "GAM-logIOR_nom", detection = "IOR", use_null = FALSE),
      list(method = "GAM-logIOR", detection = "IOR", use_null = TRUE),
      list(method = "GAM-AC_nom", detection = "AC", use_null = FALSE),
      list(method = "GAM-AC", detection = "AC", use_null = TRUE),
      list(method = "Estratificado-IOR", detection = "IOR", use_null = FALSE),
      list(method = "Estratificado-IOR_null", detection = "IOR", use_null = TRUE),
      list(method = "Estratificado-AC", detection = "AC", use_null = FALSE),
      list(method = "Estratificado-AC_null", detection = "AC", use_null = TRUE)
    )
  }
  
  metrics_list <- vector("list", length(methods_cfg))
  detail_list <- vector("list", length(methods_cfg))
  
  for (i in seq_along(methods_cfg)) {
    met <- methods_cfg[[i]]
    dt_met <- copy(benchmark_expanded_dt)
    dt_met <- detect_signal(
      dt = dt_met,
      method_name = met$method,
      detection_type = met$detection,
      use_null = met$use_null
    )
    
    success_col <- if (grepl("GAM", met$method)) "model_success" else "classic_success"
    
    metrics_list[[i]] <- calculate_benchmark_metrics(
      dt = dt_met,
      success_col = success_col,
      method_name = met$method,
      detection_type = met$detection,
      use_null = met$use_null,
      n_boot = n_boot
    )

    detail_list[[i]] <- dt_met[, .(
      evaluable = any(get(success_col), na.rm = TRUE),
      detected = any(detected, na.rm = TRUE),
      n_stages_detected = sum(detected, na.rm = TRUE),
      signal_score = if (all(is.na(signal_score))) NA_real_ else max(signal_score, na.rm = TRUE)
    ), by = triplet_id][, `:=`(
      method = met$method,
      detection_type = met$detection,
      threshold_mode = if (met$use_null) paste0("null_", percentil) else "nominal"
    )]
  }
  
  return(list(
    metrics = rbindlist(metrics_list, fill = TRUE),
    triplet_detail = rbindlist(detail_list, fill = TRUE)
  ))
}

################################################################################
# Ontogeny contrast (expected vs detected developmental window)
################################################################################

# Contrasts the curated expected high-risk developmental window against the NICHD stages 
# Only triplets flagged ontogenic_modulation == "yes" are evaluated
build_ontogeny_stage_contrast <- function(benchmark_expanded_dt, methods_cfg) {
  # The ontogeny fields only exist when the curated set was prepared with them
  ontogeny_cols <- c("ontogenic_modulation", "higher_risk_stages")
  if (!all(ontogeny_cols %in% names(benchmark_expanded_dt))) {
    warning("Ontogeny fields absent from benchmark results; skipping the ontogeny contrast. Refit the benchmark (overwrite_cache = TRUE) after updating the curated set.")
    return(data.table())
  }

  ontogenic <- benchmark_expanded_dt[ontogenic_modulation == "yes"]
  if (nrow(ontogenic) == 0) {
    return(data.table())
  }

  # Expected stages per triplet, parsed from the comma-separated curated field.
  expected_map <- unique(ontogenic[, .(triplet_id, higher_risk_stages)])
  expected_map[, expected_stages := lapply(
    strsplit(higher_risk_stages, "\\s*,\\s*"),
    function(stages) stages[nzchar(stages)]
  )]

  # Per-stage detection for each method (same detect_signal used in evaluation).
  detail_list <- vector("list", length(methods_cfg))
  for (i in seq_along(methods_cfg)) {
    met <- methods_cfg[[i]]
    dt_met <- detect_signal(copy(ontogenic), met$method, met$detection, met$use_null)
    detail_list[[i]] <- dt_met[, .(
      triplet_id, drug1_name, drug2_name, meddra_pt, stage_num, nichd, detected
    )][, `:=`(
      method = met$method,
      detection_type = met$detection,
      threshold_mode = if (met$use_null) paste0("null_", percentil) else "nominal"
    )]
  }

  contrast <- rbindlist(detail_list, fill = TRUE)
  contrast <- merge(
    contrast,
    expected_map[, .(triplet_id, higher_risk_stages, expected_stages)],
    by = "triplet_id", all.x = TRUE
  )
  contrast[, expected := mapply(function(stage, exp) stage %in% exp, nichd, expected_stages)]
  contrast[, expected_stages := NULL]
  setorder(contrast, triplet_id, method, stage_num)
  contrast[]
}

# Collapses the stage-level contrast to one row per (triplet, method)
summarize_ontogeny_contrast <- function(contrast_dt) {
  if (nrow(contrast_dt) == 0) {
    return(data.table())
  }
  contrast_dt[, .(
    expected_stages = paste(sort(unique(nichd[expected])), collapse = ","),
    detected_stages = paste(sort(unique(nichd[detected])), collapse = ","),
    n_expected = sum(expected, na.rm = TRUE),
    n_detected = sum(detected, na.rm = TRUE),
    n_expected_detected = sum(expected & detected, na.rm = TRUE),
    expected_window_detected = any(expected & detected, na.rm = TRUE),
    detection_outside_expected = any(detected & !expected, na.rm = TRUE)
  ), by = .(triplet_id, drug1_name, drug2_name, meddra_pt, method, detection_type, threshold_mode, higher_risk_stages)]
}

################################################################################
# Vocabulary-based benchmark mapping
################################################################################

# Dictionary layer used by scripts 01 and 02 to map benchmark drugs and events

# Shared OMOP vocabulary root at the workspace level.
vocabulary_dir <- "../data/vocabulary/vocabulary_SNOMED_MEDDRA_RxNorm_ATC"
ruta_concept <- file.path(vocabulary_dir, "CONCEPT.csv")
ruta_concept_ancestor <- file.path(vocabulary_dir, "CONCEPT_ANCESTOR.csv")

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
  
  # Order-independent signature for fuzzy matching of combination drug names.
  atc_dt[, atc_signature := vapply(
    strsplit(gsub("\\s*-\\s*", " and ", atc_name_key), "\\s+and\\s+"),
    function(parts) paste(sort(trimws(parts)), collapse = " | "),
    character(1)
  )]

  # Keep one row per ATC code, preferring the standard concept (S) with the lowest ID.
  atc_dt[, `:=`(
    standard_priority = fifelse(standard_concept == "S", 1L, 0L),
    atc_concept_id_num = as.numeric(atc_concept_id)
  )]
  setorder(atc_dt, atc_concept_code, -standard_priority, atc_concept_id_num)
  atc_dt <- atc_dt[, .SD[1L], by = atc_concept_code]
  atc_dt[, c("standard_priority", "atc_concept_id_num") := NULL]
  
  return(atc_dt[])
}

build_drug_translation_table <- function(drug_info_dt = NULL) {
  atc_dt <- build_atc_vocabulary_map()
  # Strip qualifiers after semicolon/comma to derive the base drug name.
  atc_dt[, base_name := normalize_vocabulary_key(sub("[;,].*", "", atc_concept_name))]

  # One canonical ID per base name (lowest numeric ID = oldest OMOP concept).
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

build_meddra_hierarchy_map <- function(
  rollup_level = "HLT",
  concept_path = ruta_concept,
  ancestor_path = ruta_concept_ancestor
) {
  if (!file.exists(concept_path)) {
    stop(sprintf("CONCEPT.csv was not found at %s", concept_path))
  }
  if (!file.exists(ancestor_path)) {
    stop(sprintf("CONCEPT_ANCESTOR.csv was not found at %s", ancestor_path))
  }
  
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
  
  # Prefer the shortest (most specific) ancestor path; break ties by lowest concept ID.
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

load_ade_modeling_data <- function(rollup_level = "HLT") {
  ade_raw_dt <- fread(ruta_ade_raw)

  if (include_sex) {
    # Normalize single-letter sex codes, consistent with fit_gam().
    ade_raw_dt[, sex := toupper(trimws(sex))]
    ade_raw_dt[sex == "M", sex := "MALE"]
    ade_raw_dt[sex == "F", sex := "FEMALE"]
    ade_raw_dt[, sex := factor(sex, levels = c("MALE", "FEMALE"))]
  }

  # Collapse ATC synonyms to a single canonical concept ID
  translation_table <- build_drug_translation_table()
  ade_raw_dt[, atc_concept_id := as.character(atc_concept_id)]
  ade_raw_dt <- merge(
    ade_raw_dt,
    translation_table[, .(atc_concept_id, canonical_id)],
    by = "atc_concept_id",
    all.x = TRUE
  )
  ade_raw_dt[!is.na(canonical_id), atc_concept_id := canonical_id]
  ade_raw_dt[, canonical_id := NULL]

  # Roll PT-level events up to HLT/HLGT so event identifiers match the benchmark.
  if (rollup_level != "PT") {
    event_map <- build_meddra_hierarchy_map(rollup_level)[, .(meddra_concept_id, rollup_id)]
    ade_raw_dt[, meddra_concept_id := as.character(meddra_concept_id)]
    ade_raw_dt <- merge(ade_raw_dt, event_map, by = "meddra_concept_id", all.x = TRUE)
    ade_raw_dt[!is.na(rollup_id), meddra_concept_id := rollup_id]
    ade_raw_dt[, rollup_id := NULL]
  }

  # De-duplicate after roll-up: a report may now have the same drug-event pair via multiple PTs mapping to the same HLT.
  ade_raw_dt <- unique(
    ade_raw_dt,
    by = c("safetyreportid", "atc_concept_id", "meddra_concept_id")
  )
  ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
  ade_raw_dt[, nichd_num := as.integer(nichd)]

  return(ade_raw_dt)
}

prepare_benchmark_reference_set <- function(
  benchmark_triplets_path = "../ddi_reference_set/results/curated_pediatric_ddi_reference_set/curated_pediatric_ddi_triplets.csv",
  benchmark_sources_path = "../ddi_reference_set/results/curated_pediatric_ddi_reference_set/curated_pediatric_ddi_sources.csv",
  output_path = NULL,
  rollup_level = "HLT"
) {
  benchmark_dt <- fread(benchmark_triplets_path)
  # Reference sets curated before control_type existed are all positive controls.
  if (!"control_type" %in% names(benchmark_dt)) benchmark_dt[, control_type := "positive"]
  benchmark_dt[is.na(control_type) | !nzchar(trimws(control_type)), control_type := "positive"]
  translation_table <- build_drug_translation_table()
  # One canonical ID per ATC code, lowest numeric value.
  atc_map <- unique(translation_table[, .(atc_concept_code, canonical_id)])
  atc_map[, canonical_id_num := as.numeric(canonical_id)]
  setorder(atc_map, atc_concept_code, canonical_id_num)
  atc_map <- atc_map[, .SD[1L], by = atc_concept_code]
  atc_map[, canonical_id_num := NULL]

  benchmark_dt <- merge(benchmark_dt, atc_map, by.x = "drug1_atc", by.y = "atc_concept_code", all.x = TRUE)
  setnames(benchmark_dt, "canonical_id", "drug1_id")
  benchmark_dt <- merge(benchmark_dt, atc_map, by.x = "drug2_atc", by.y = "atc_concept_code", all.x = TRUE)
  setnames(benchmark_dt, "canonical_id", "drug2_id")

  rollup_col <- switch(
    rollup_level,
    "PT" = "meddra_concept_id",
    "HLT" = "meddra_concept_id_2",
    "HLGT" = "meddra_concept_id_3",
    stop(sprintf("Unknown MedDRA roll-up level: %s", rollup_level))
  )
  
  benchmark_dt[, `:=`(
    mapped_PT = !is.na(drug1_id) & !is.na(drug2_id) & !is.na(meddra_concept_id),
    mapped_HLT = !is.na(drug1_id) & !is.na(drug2_id) & !is.na(meddra_concept_id_2),
    mapped_HLGT = !is.na(drug1_id) & !is.na(drug2_id) & !is.na(meddra_concept_id_3)
  )]
  benchmark_dt[, mapping_success := !is.na(drug1_id) & !is.na(drug2_id) & !is.na(get(rollup_col))]
  
  # Canonical drug pair order: drugA <= drugB by numeric ID, matching ade_raw convention.
  benchmark_dt[mapping_success == TRUE, `:=`(
    drugA = as.character(pmin(as.numeric(drug1_id), as.numeric(drug2_id))),
    drugB = as.character(pmax(as.numeric(drug1_id), as.numeric(drug2_id))),
    meddra = as.character(get(rollup_col)),
    type = "benchmark"
  )]
  benchmark_dt[mapping_success == FALSE, `:=`(
    drugA = NA_character_,
    drugB = NA_character_,
    meddra = NA_character_,
    type = "benchmark"
  )]
  
  if (file.exists(benchmark_sources_path)) {
    sources_dt <- fread(benchmark_sources_path)
    source_counts <- sources_dt[, .(n_sources = .N), by = triplet_id]
    benchmark_dt <- merge(benchmark_dt, source_counts, by = "triplet_id", all.x = TRUE)
  }
  
  if (!is.null(output_path)) {
    fwrite(benchmark_dt, output_path)
  }
  
  return(benchmark_dt[])
}

