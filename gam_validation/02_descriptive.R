################################################################################
# Descriptive analysis of the original dataset and the control sets
# Script: 02_descriptive.R
################################################################################

source("00_functions.R", local = TRUE)

################################################################################
# Configuration
################################################################################

# Max points per group for ECDF plots (avoids overloading the render)
max_plot_points <- 50000

ruta_augmentation <- paste0("./results/", suffix, "/augmentation_results/")
output_dir <- paste0("./results/", suffix, "/descriptive_results/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# NICHD stage axis labels (consistent with 30_metrics / 40_graphs)
nichd_labels <- c(
  term_neonatal = "Term neonatal",
  infancy = "Infancy",
  toddler = "Toddler",
  early_childhood = "Early childhood",
  middle_childhood = "Middle childhood",
  early_adolescence = "Early adolescence",
  late_adolescence = "Late adolescence"
)

# Dynamic display labels and colors (used in Part B)
dynamic_labels <- c(
  uniform = "Uniform", increase = "Increase", decrease = "Decrease",
  plateau = "Plateau", inverse_plateau = "Inverse plateau"
)
dynamic_colors <- c(
  uniform = "#95A5A6", increase = "#E41A1C", decrease = "#377EB8",
  plateau = "#4DAF4A", inverse_plateau = "#984EA3"
)

# Helpers: persist tables and figures with a confirmation message
save_table <- function(dt, filename) {
  fwrite(dt, paste0(output_dir, filename))
  message(sprintf("  table saved: %s (%d rows)", filename, nrow(dt)))
}
save_fig <- function(p, filename, width = 10, height = 7) {
  ggsave(paste0(output_dir, filename), p, width = width, height = height, dpi = 300, bg = "white")
  message(sprintf("  figure saved: %s", filename))
}

################################################################################
# A- Original dataset
################################################################################

# A0. Load and canonicalize, mirroring 10_augmentation.R
ade_raw_dt <- fread(ruta_ade_raw)
ade_raw_dt[, atc_concept_id := as.character(atc_concept_id)]

# Canonicalize ATC ids via the shared OMOP vocabulary to match the pipeline's drug mapping
translation_table <- build_drug_translation_table()
ade_raw_dt <- merge(
  ade_raw_dt,
  translation_table[, .(atc_concept_id, canonical_id)],
  by = "atc_concept_id", all.x = TRUE
)
ade_raw_dt[!is.na(canonical_id), atc_concept_id := canonical_id]
ade_raw_dt[, canonical_id := NULL]

# De-duplicate report-drug-event rows and build the ordered NICHD factor
ade_raw_dt <- unique(ade_raw_dt, by = c("safetyreportid", "atc_concept_id", "meddra_concept_id"))
ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
ade_raw_dt[, nichd_num := as.integer(nichd)]

# Canonical id -> representative drug name (for top-drug tables)
drug_names <- unique(translation_table[atc_concept_id == canonical_id,
                                       .(atc_concept_id = canonical_id, drug_name = atc_concept_name)])

# A1. Overall dimensions
n_reports <- uniqueN(ade_raw_dt$safetyreportid)
drugs_per_report <- ade_raw_dt[, .(n_drugs = uniqueN(atc_concept_id)), by = safetyreportid]
events_per_report <- ade_raw_dt[, .(n_events = uniqueN(meddra_concept_id)), by = safetyreportid]

dataset_overview <- data.table(
  metric = c("rows", "reports", "unique_drugs", "unique_events",
             "mean_drugs_per_report", "median_drugs_per_report",
             "mean_events_per_report", "median_events_per_report",
             "pct_reports_polypharmacy_ge2"),
  value = c(nrow(ade_raw_dt), n_reports,
            uniqueN(ade_raw_dt$atc_concept_id), uniqueN(ade_raw_dt$meddra_concept_id),
            round(mean(drugs_per_report$n_drugs), 3), median(drugs_per_report$n_drugs),
            round(mean(events_per_report$n_events), 3), median(events_per_report$n_events),
            round(mean(drugs_per_report$n_drugs >= 2) * 100, 2))
)
message(sprintf("Dataset: %s rows | %s reports | %s drugs | %s events",
  format(nrow(ade_raw_dt), big.mark = ","), format(n_reports, big.mark = ","),
  format(uniqueN(ade_raw_dt$atc_concept_id), big.mark = ","),
  format(uniqueN(ade_raw_dt$meddra_concept_id), big.mark = ",")))
print(dataset_overview)
save_table(dataset_overview, "dataset_overview.csv")

# A2. Reports per NICHD stage
reports_by_stage <- unique(ade_raw_dt[, .(safetyreportid, nichd, nichd_num)])[
  , .(n_reports = .N), by = .(nichd_num, nichd)][order(nichd_num)]
reports_by_stage[, pct_reports := round(100 * n_reports / sum(n_reports), 2)]
save_table(reports_by_stage, "reports_by_nichd_stage.csv")

p_stage <- ggplot(reports_by_stage, aes(x = nichd, y = n_reports)) +
  geom_col(fill = "#2980B9", alpha = 0.85, color = "white", width = 0.7) +
  geom_text(aes(label = sprintf("%s\n(%.1f%%)", scales::comma(n_reports), pct_reports)),
            vjust = -0.3, size = 3) +
  scale_x_discrete(labels = nichd_labels) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Reports per NICHD developmental stage", x = "NICHD stage", y = "Reports") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_stage, "fig_reports_by_nichd_stage.png", width = 11, height = 7)

# A3-A4. Drugs- and events-per-report distributions
distribution_summary <- function(x, label) {
  data.table(
    quantity = label,
    n_reports = length(x),
    min = min(x), p25 = quantile(x, 0.25), median = median(x),
    mean = round(mean(x), 3), p75 = quantile(x, 0.75),
    p95 = quantile(x, 0.95), max = max(x)
  )
}
per_report_summary <- rbind(
  distribution_summary(drugs_per_report$n_drugs, "drugs_per_report"),
  distribution_summary(events_per_report$n_events, "events_per_report")
)
save_table(per_report_summary, "per_report_distributions.csv")

p_drugs <- ggplot(drugs_per_report, aes(x = n_drugs)) +
  geom_histogram(binwidth = 1, fill = "#16A085", alpha = 0.85, color = "white") +
  scale_x_continuous(limits = c(0, quantile(drugs_per_report$n_drugs, 0.99))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Distribution of drugs per report (polypharmacy)",
       subtitle = "Truncated at the 99th percentile for readability",
       x = "Distinct drugs per report", y = "Reports")
save_fig(p_drugs, "fig_drugs_per_report.png")

# A5. Reports per drug (with names) -> top 25
reports_per_drug <- unique(ade_raw_dt[, .(safetyreportid, atc_concept_id)])[
  , .(n_reports = .N), by = atc_concept_id][order(-n_reports)]
reports_per_drug <- merge(reports_per_drug, drug_names, by = "atc_concept_id", all.x = TRUE)
setorder(reports_per_drug, -n_reports)
save_table(reports_per_drug[1:min(25, .N)], "top25_drugs.csv")

p_top_drugs <- ggplot(reports_per_drug[1:min(25, .N)],
  aes(x = reorder(ifelse(is.na(drug_name), atc_concept_id, drug_name), n_reports), y = n_reports)) +
  geom_col(fill = "#C0392B", alpha = 0.85) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Top 25 drugs by number of reports", x = NULL, y = "Reports")
save_fig(p_top_drugs, "fig_top25_drugs.png", width = 10, height = 9)

# A6. Reports per event (with MedDRA PT names) -> top 25
event_names <- build_meddra_name_map()
reports_per_event <- unique(ade_raw_dt[, .(safetyreportid, meddra_concept_id)])[
  , .(n_reports = .N), by = meddra_concept_id][order(-n_reports)]
reports_per_event[, meddra_concept_id := as.character(meddra_concept_id)]
reports_per_event <- merge(reports_per_event, event_names, by = "meddra_concept_id", all.x = TRUE)
setcolorder(reports_per_event, c("meddra_concept_id", "meddra_name", "n_reports"))
setorder(reports_per_event, -n_reports)
save_table(reports_per_event[1:min(25, .N)], "top25_events.csv")

# A7. Co-administration structure: all observed drug pairs and reports-per-pair
drugs_by_report <- unique(ade_raw_dt[!is.na(atc_concept_id), .(safetyreportid, atc_concept_id)])

drug_pairs_dt <- drugs_by_report[, {
  drug_list <- unique(atc_concept_id)
  if (length(drug_list) >= 2) {
    pairs <- t(combn(drug_list, 2))
    data.table(drugA = pmin(pairs[, 1], pairs[, 2]), drugB = pmax(pairs[, 1], pairs[, 2]))
  } else data.table()
}, by = safetyreportid]

coadmin_counts <- drug_pairs_dt[, .(n_reports = .N), by = .(drugA, drugB)][order(-n_reports)]

coadmin_overview <- data.table(
  metric = c("unique_drug_pairs", "mean_reports_per_pair", "median_reports_per_pair",
             "max_reports_per_pair", "pairs_with_ge2_reports"),
  value = c(nrow(coadmin_counts), round(mean(coadmin_counts$n_reports), 3),
            median(coadmin_counts$n_reports), max(coadmin_counts$n_reports),
            sum(coadmin_counts$n_reports >= 2))
)
save_table(coadmin_overview, "coadministration_overview.csv")

# ECDF: single drug vs co-administration report counts
combined_ecdf <- rbind(
  reports_per_drug[, .(entity = "Single drug", n_reports)],
  coadmin_counts[, .(entity = "Co-administration", n_reports)]
)
set.seed(7113)  # reproducible downsampling for ECDF plot
plot_ecdf <- combined_ecdf[, if (.N > max_plot_points) .SD[sample(.N, max_plot_points)] else .SD, by = entity]
p_ecdf <- ggplot(plot_ecdf, aes(x = n_reports, color = entity)) +
  stat_ecdf(linewidth = 1.2) +
  scale_x_log10(labels = scales::comma, breaks = c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000)) +
  scale_color_manual(values = c("Single drug" = "#C0392B", "Co-administration" = "#16A085")) +
  annotation_logticks(sides = "b") +
  labs(title = "ECDF: single-drug vs co-administration report counts",
       subtitle = sprintf("Median single drug: %.0f | Median co-administration: %.0f",
                          median(reports_per_drug$n_reports), median(coadmin_counts$n_reports)),
       x = "Number of reports (log scale)", y = "Cumulative proportion", color = "Entity")
save_fig(p_ecdf, "fig_drug_vs_coadmin_ecdf.png")

message("\nPart A (original dataset) complete.\n")

################################################################################
# B- Control sets and injection (requires 10_augmentation.R outputs)
################################################################################

ruta_pos <- paste0(ruta_augmentation, "positive_triplets_results.rds")

pos_meta <- fread(paste0(ruta_augmentation, "positive_triplets_metadata.csv"))
pos_results <- readRDS(ruta_pos)
pos_base <- pos_results[reduction_pct == 0]

coadmin_pos <- fread(paste0(ruta_augmentation, "positive_coadmin_by_stage.csv"))
coadmin_neg <- fread(paste0(ruta_augmentation, "negative_coadmin_by_stage.csv"))

message(sprintf("Positive base: %d", nrow(pos_base)))

# B1. Injection success rate per dynamic
injection_summary <- pos_base[, .(
  n_total = .N,
  n_successful = sum(injection_success == TRUE, na.rm = TRUE),
  success_rate = round(mean(injection_success == TRUE, na.rm = TRUE) * 100, 1)
), by = dynamic]

p_inj_rate <- ggplot(injection_summary,
  aes(x = reorder(dynamic, -success_rate), y = success_rate, fill = dynamic)) +
  geom_col(alpha = 0.85, color = "white", width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", success_rate, n_successful)),
            vjust = -0.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = dynamic_colors, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.15)), limits = c(0, 100)) +
  scale_x_discrete(labels = dynamic_labels) +
  labs(title = "Injection success rate per dynamic type",
       subtitle = "Triplets with at least one event injected in key stages",
       x = "Dynamic", y = "Success rate (%)")
save_fig(p_inj_rate, "fig_injection_success_rate.png")

# B2. Co-administration coverage per NICHD stage (positive vs negative)
pos_unique_trips <- pos_meta[, .SD[1L], by = base_triplet_id][, triplet_id]
coadmin_pos_unique <- coadmin_pos[triplet_id %in% pos_unique_trips]

pos_stage_cov <- coadmin_pos_unique[, .(
  set = "Positive",
  n_triplets = uniqueN(triplet_id[n_coadmin_stage > 0]),
  n_total = uniqueN(triplet_id),
  pct_with_coadmin = round(uniqueN(triplet_id[n_coadmin_stage > 0]) / uniqueN(triplet_id) * 100, 1)
), by = .(nichd_num, nichd)]
neg_stage_cov <- coadmin_neg[, .(
  set = "Negative",
  n_triplets = uniqueN(triplet_id[n_coadmin_stage > 0]),
  n_total = uniqueN(triplet_id),
  pct_with_coadmin = round(uniqueN(triplet_id[n_coadmin_stage > 0]) / uniqueN(triplet_id) * 100, 1)
), by = .(nichd_num, nichd)]
stage_cov <- rbind(pos_stage_cov, neg_stage_cov)
stage_cov[, nichd := factor(niveles_nichd[nichd_num], levels = niveles_nichd, ordered = TRUE)]
stage_cov[, set := factor(set, levels = c("Positive", "Negative"))]

p_cov <- ggplot(stage_cov, aes(x = nichd, y = n_triplets, fill = set)) +
  geom_col(position = position_dodge(width = 0.85), alpha = 0.85, color = "white", width = 0.75) +
  scale_fill_manual(values = c("Positive" = "#4DAF4A", "Negative" = "#E41A1C"), name = "Control set") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  scale_x_discrete(labels = nichd_labels) +
  labs(title = "Triplets with co-administration per NICHD stage",
       subtitle = "Unique triplets with n_coadmin > 0 at each stage",
       x = "NICHD stage", y = "Number of triplets") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_cov, "fig_coadmin_coverage_by_stage.png", width = 12, height = 7)

message("\nPart B (control sets) complete.\n")

message(sprintf("Descriptive analysis complete. Results saved to: %s", output_dir))
