################################################################################
# Dual-curation validation - human vs agent comparison
# Script 02_compare_dual_curation
################################################################################
#
# Reads the two completed blind workbooks (human and agent), maps every found
# triplet to a stable key (assigned pair + resolved MedDRA PT) and contrasts the
# two sets. Reports the basic agreement only: the per-triplet detail and a small
# summary (counts, triplet-level Jaccard, pair-level observed agreement).
#
# Matching key: (pair_id, PT concept id). pair_id is shared across both workbooks
# (same assigned pairs), and the event is rolled up to PT so entries made at
# different MedDRA levels still compare. A pair recorded as no_triplet_found
# contributes no key; pairs where both curators found nothing count as agreed
# empties at the pair level.
#
# Run from the dual_curation_validation/ root, after both workbooks are filled:
#   & 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' scripts\R\02_compare_dual_curation.R

source("../00_functions.R")
library(openxlsx)

################################################################################
# Configuration
################################################################################

voc_dir <- "../../data/vocabulary/vocabulary_SNOMED_MEDDRA_RxNorm_ATC"
concept_path <- file.path(voc_dir, "CONCEPT.csv")
ancestor_path <- file.path(voc_dir, "CONCEPT_ANCESTOR.csv")
relationship_path <- file.path(voc_dir, "CONCEPT_RELATIONSHIP.csv")

input_dir <- "./input"
results_dir <- "./results"
human_xlsx <- file.path(input_dir, "dual_curation_human.xlsx")
agent_xlsx <- file.path(input_dir, "dual_curation_agent.xlsx")
sample_key_file <- file.path(results_dir, "sampled_pairs_key.csv")

triplet_compare_file <- file.path(results_dir, "comparison_triplets.csv")
summary_file <- file.path(results_dir, "comparison_summary.csv")

for (f in c(human_xlsx, agent_xlsx, sample_key_file)) {
  if (!file.exists(f)) stop(sprintf("%s not found. Run script 01 and fill both workbooks first.", f))
}

event_cols <- c("event_llt", "event_pt", "event_hlt", "event_hlgt")

# Each pair sheet (pair_NN) carries the triplet table from this row down; rows
# 1-2 hold the read-only drug context. Kept in sync with script 01 (header_row).
header_row <- 3L

################################################################################
# 1. Read a curator's found triplets and resolve them to PT
################################################################################
# A "found triplet" is a row that is not flagged no_triplet_found and carries an
# event at some MedDRA level. Rows that only mark a pair as empty contribute no
# key. The workbook holds one sheet per assigned pair (pair_NN); the sheet name
# carries pair_id and the triplet table starts at header_row. Returns one row per
# (pair_id, PT) with the resolved PT name/id.

is_yes <- function(x) tolower(trimws(ifelse(is.na(x), "", x))) %in% c("yes", "y", "true", "1")
has_value <- function(x) !is.na(x) & nzchar(trimws(x))

read_curator_triplets <- function(path, role) {
  empty_result <- data.table(pair_id = integer(), meddra_concept_id = character(),
                             meddra_pt = character(), role = character())

  pair_sheets <- grep("^pair_\\d+$", getSheetNames(path), value = TRUE)
  per_sheet <- lapply(pair_sheets, function(sheet) {
    tbl <- tryCatch(
      as.data.table(read.xlsx(path, sheet = sheet, startRow = header_row)),
      error = function(e) NULL
    )
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
    tbl[, pair_id := as.integer(sub("^pair_", "", sheet))]
    tbl
  })
  dt <- rbindlist(per_sheet, fill = TRUE)
  if (nrow(dt) == 0) return(empty_result)

  for (col in c("no_triplet_found", event_cols)) {
    if (!col %in% names(dt)) dt[, (col) := NA_character_]
    dt[, (col) := as.character(get(col))]
  }

  any_event <- has_value(dt$event_llt) | has_value(dt$event_pt) |
    has_value(dt$event_hlt) | has_value(dt$event_hlgt)
  found <- dt[!is_yes(no_triplet_found) & any_event & !is.na(pair_id)]
  if (nrow(found) == 0) return(empty_result)

  # resolve_meddra_event_levels keys on triplet_id; give each row a unique one.
  found[, row_uid := paste0(role, "_", seq_len(.N))]
  resolved <- resolve_meddra_event_levels(
    found[, .(triplet_id = row_uid, event_llt, event_pt, event_hlt, event_hlgt)],
    concept_path = concept_path,
    ancestor_path = ancestor_path,
    relationship_path = relationship_path
  )
  out <- merge(
    found[, .(row_uid, pair_id)],
    resolved[, .(triplet_id, meddra_pt, meddra_concept_id)],
    by.x = "row_uid", by.y = "triplet_id"
  )
  out[, role := role]
  # Collapse duplicates: the same (pair, PT) entered twice by a curator is one.
  unique(out[, .(pair_id, meddra_concept_id, meddra_pt, role)])
}

human <- read_curator_triplets(human_xlsx, "human")
agent <- read_curator_triplets(agent_xlsx, "agent")

################################################################################
# 2. Triplet-level comparison (keyed on pair_id + PT)
################################################################################

key_cols <- c("pair_id", "meddra_concept_id")
both <- merge(
  human[, c(key_cols, "meddra_pt"), with = FALSE],
  agent[, c(key_cols, "meddra_pt"), with = FALSE],
  by = key_cols, all = TRUE, suffixes = c("_human", "_agent")
)
both[, meddra_pt := fcoalesce(meddra_pt_human, meddra_pt_agent)]
both[, in_human := !is.na(meddra_pt_human)]
both[, in_agent := !is.na(meddra_pt_agent)]
both[, status := fcase(
  in_human & in_agent, "matched",
  in_human & !in_agent, "human_only",
  !in_human & in_agent, "agent_only"
)]

# Attach the pair labels from the audit key for a readable output.
sample_key <- fread(sample_key_file)
triplet_compare <- merge(
  both[, .(pair_id, meddra_pt, meddra_concept_id, in_human, in_agent, status)],
  sample_key[, .(pair_id, drug1, drug2, pair_coreport)],
  by = "pair_id", all.x = TRUE
)
setcolorder(triplet_compare,
            c("pair_id", "drug1", "drug2", "meddra_pt", "meddra_concept_id",
              "in_human", "in_agent", "status", "pair_coreport"))
setorder(triplet_compare, pair_id, meddra_pt)
fwrite(triplet_compare, triplet_compare_file)

################################################################################
# 3. Summary metrics (counts, triplet Jaccard, pair-level observed agreement)
################################################################################

n_matched <- both[status == "matched", .N]
n_human_only <- both[status == "human_only", .N]
n_agent_only <- both[status == "agent_only", .N]
n_union <- n_matched + n_human_only + n_agent_only
jaccard <- if (n_union == 0) NA_real_ else round(n_matched / n_union, 4)

# Pair-level observed agreement: per assigned pair, did each curator find >= 1
# triplet? Pairs where both found nothing agree (both empty); start from the full
# pair list so those are counted. percent_observed_agreement is truth-free.
pair_found <- both[, .(human_nonempty = any(in_human), agent_nonempty = any(in_agent)),
                   by = pair_id]
all_pairs <- unique(sample_key[, .(pair_id)])
pair_found <- pair_found[all_pairs, on = "pair_id"]
pair_found[is.na(human_nonempty), human_nonempty := FALSE]
pair_found[is.na(agent_nonempty), agent_nonempty := FALSE]
pair_found[, agree := human_nonempty == agent_nonempty]
n_pairs <- nrow(pair_found)
percent_observed_agreement <- round(mean(pair_found$agree), 4)

summary_dt <- data.table(
  metric = c(
    "n_triplets_human", "n_triplets_agent", "n_matched",
    "n_human_only", "n_agent_only", "n_union", "jaccard_triplets",
    "n_pairs", "percent_observed_agreement"
  ),
  value = c(
    nrow(human), nrow(agent), n_matched,
    n_human_only, n_agent_only, n_union, jaccard,
    n_pairs, percent_observed_agreement
  )
)
fwrite(summary_dt, summary_file)

################################################################################
# 4. Report
################################################################################

cat("\nDual-curation comparison (human vs agent)\n")
cat("-----------------------------------------\n")
cat(sprintf("  triplets: human %d | agent %d | matched %d\n",
            nrow(human), nrow(agent), n_matched))
cat(sprintf("  human-only %d | agent-only %d | union %d | Jaccard %s\n",
            n_human_only, n_agent_only, n_union, jaccard))
cat(sprintf("  pairs: %d | observed agreement %s\n",
            n_pairs, percent_observed_agreement))
cat("\nOutputs:\n")
cat("  ", triplet_compare_file, "\n")
cat("  ", summary_file, "\n")
