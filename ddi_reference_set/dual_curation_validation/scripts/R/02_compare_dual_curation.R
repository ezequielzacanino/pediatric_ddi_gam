################################################################################
# Dual-curation validation - human vs agent comparison
# Script 02_compare_dual_curation
################################################################################
#
# Reads the two completed blind workbooks (human and agent) 
# maps every found triplet and contrasts the two sets. 
# Reports agreement
#
#
# Run from the dual_curation_validation/ root, after both workbooks are filled

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

# Each pair sheet (pair_NN) carries the triplet table from this row down 
# rows 1-2 hold the read-only drug context. Kept in sync with script 01 (header_row).
header_row <- 3L

################################################################################
# 1. Read found triplets and resolve to PT
################################################################################
# A "found triplet" is a row that is not flagged no_triplet_found
# The workbook holds one sheet per assigned pair (pair_NN)
# the sheet name carries pair_id 

is_yes <- function(x) tolower(trimws(ifelse(is.na(x), "", x))) %in% c("yes", "y", "true", "1")
has_value <- function(x) !is.na(x) & nzchar(trimws(x))

read_curator_triplets <- function(path, role) {
  empty_result <- data.table(pair_id = integer(), meddra_concept_id = character(),
                             meddra_pt = character(), hlt_concept_id = character(),
                             role = character())

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
    resolved[, .(triplet_id, meddra_pt, meddra_concept_id, hlt_concept_id = meddra_concept_id_2)],
    by.x = "row_uid", by.y = "triplet_id"
  )
  out[, role := role]
  # Collapse duplicates: the same (pair, PT) entered twice by a curator is one.
  unique(out[, .(pair_id, meddra_concept_id, meddra_pt, hlt_concept_id, role)])
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
# 3. Summary metrics 
################################################################################

# Triplet-level counts (pair, PT roll-up)
# "matched" needs both curators to have recorded the same (pair, PT) 
# human_only/agent_only are exclusive to one.
n_matched <- both[status == "matched", .N]
n_human_only <- both[status == "human_only", .N]
n_agent_only <- both[status == "agent_only", .N]

# Same triplet concordance after rolling both curators' events up to the HLT
human_hlt <- unique(human[!is.na(hlt_concept_id), .(pair_id, hlt_concept_id)])
agent_hlt <- unique(agent[!is.na(hlt_concept_id), .(pair_id, hlt_concept_id)])
n_matched_hlt <- nrow(merge(human_hlt, agent_hlt, by = c("pair_id", "hlt_concept_id")))

# Pair-level view
pair_found <- both[, .(human_nonempty = any(in_human), agent_nonempty = any(in_agent)),
                   by = pair_id]
all_pairs <- unique(sample_key[, .(pair_id)])
pair_found <- pair_found[all_pairs, on = "pair_id"]
pair_found[is.na(human_nonempty), human_nonempty := FALSE]
pair_found[is.na(agent_nonempty), agent_nonempty := FALSE]
pair_found[, agree := human_nonempty == agent_nonempty]
n_pairs <- nrow(pair_found)
# Matched pairs: both curators independently found >= 1 triplet for the same pair.
n_matched_pairs <- pair_found[human_nonempty & agent_nonempty, .N]
percent_observed_agreement <- round(mean(pair_found$agree), 4)

summary_dt <- data.table(
  metric = c(
    "n_triplets_human", "n_triplets_agent",
    "n_human_only", "n_agent_only", "n_matched", "n_matched_hlt",
    "n_pairs", "n_matched_pairs", "percent_observed_agreement"
  ),
  value = c(
    nrow(human), nrow(agent),
    n_human_only, n_agent_only, n_matched, n_matched_hlt,
    n_pairs, n_matched_pairs, percent_observed_agreement
  )
)
fwrite(summary_dt, summary_file)

################################################################################
# 4. Report
################################################################################

cat("\nDual-curation comparison (human vs agent)\n")
cat("-----------------------------------------\n")
cat(sprintf("  triplets: human %d | agent %d | matched %d (PT) / %d (HLT)\n",
            nrow(human), nrow(agent), n_matched, n_matched_hlt))
cat(sprintf("  human-only %d | agent-only %d\n",
            n_human_only, n_agent_only))
cat(sprintf("  pairs: %d | matched pairs %d | observed agreement %s\n",
            n_pairs, n_matched_pairs, percent_observed_agreement))
cat("\nOutputs:\n")
cat("  ", triplet_compare_file, "\n")
cat("  ", summary_file, "\n")
