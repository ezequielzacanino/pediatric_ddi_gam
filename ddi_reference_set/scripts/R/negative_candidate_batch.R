# Feeds negative-control candidates in pair-diverse batches (one drug pair does not
# dominate a batch), paged by a cursor. Single-drug counts are withheld so selection
# cannot be ranked by them; mono-drug attribution is decided against the label.
#   Rscript scripts/R/negative_candidate_batch.R --start <cursor> [--n <count>]

library(data.table)

csv_path <- "results/negative_control_candidates/negative_control_candidates.csv"
wb_path  <- "input/ddi_reference_input.xlsx"

args <- commandArgs(trailingOnly = TRUE)
opt <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
start <- as.integer(opt("--start", "1"))
n     <- min(as.integer(opt("--n", "5")), 25L)
if (is.na(start) || start < 1) stop("--start must be a positive integer")

d <- fread(csv_path)
d[, csv_row := .I]
d[, pair := tolower(paste(pmin(drug1, drug2), pmax(drug1, drug2), sep = "|"))]
d[, key  := paste(pair, tolower(trimws(event_pt)), sep = "|")]

# Skip triplets already in the workbook (either drug order).
done <- character(0)
if (requireNamespace("openxlsx", quietly = TRUE) && file.exists(wb_path)) {
  tr <- tryCatch(openxlsx::read.xlsx(wb_path, sheet = "triplets"), error = function(e) NULL)
  if (!is.null(tr) && nrow(tr))
    done <- tolower(paste(paste(pmin(tr$drug1, tr$drug2), pmax(tr$drug1, tr$drug2), sep = "|"),
                          trimws(tr$event_pt), sep = "|"))
}
d <- d[!(key %in% done)]

# Round-robin across pairs (occurrence rank, then random-order row) so consecutive
# candidates come from different drug pairs.
d[, rank := rowid(pair)]
setorder(d, rank, csv_row)
d[, pos := .I]

batch <- d[pos >= start & pos < start + n]
cols <- c("pos", "drug1", "drug2", "event_pt", "match_strategy",
          "known_interacting_pair", "matched_positive", "pair_coreport", "triplet_coreport")
if (!nrow(batch)) cat("no candidates at or beyond cursor", start, "\n") else print(batch[, ..cols])
cat("\nnext cursor: --start", start + n, "\n")
