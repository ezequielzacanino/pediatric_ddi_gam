# DDInter drug-drug interaction lookup for the negative-control screen.
#
# DDInter is pair-level (not event-level): it answers "does DDInter list a
# documented interaction between these two drugs, and at what risk Level?". A
# clean negative control needs the pair absent (or, if present, present for a
# different mechanism/event that the curator then rules out against the label).
#
# Data: the public DDInter bulk CSVs live under input/ddinter/ (columns
# DDInterID_A, Drug_A, DDInterID_B, Drug_B, Level). Source:
# http://ddinter.scbdd.com/static/media/download/ddinter_downloads_code_<ATC>.csv
#
# Run from the subproject root:
#   Rscript agent/tools/ddinter_lookup.R "<drug1>" "<drug2>"   # pair check
#   Rscript agent/tools/ddinter_lookup.R --find "<name>"       # resolve DDInter name
#
# Accepts the workbook ATC format ("sustancia; via"): only the substance (text
# before ";") is used for matching.

library(data.table)

ddinter_dir <- Sys.getenv("DDINTER_DIR", "input/ddinter")

# substance token from an ATC 5th dropdown value ("valproic acid; systemic" -> "valproic acid")
substance <- function(x) trimws(sub(";.*$", "", x))
norm <- function(x) gsub("\\s+", " ", tolower(trimws(x)))

load_ddinter <- function() {
  files <- list.files(ddinter_dir, pattern = "\\.csv$", full.names = TRUE)
  if (!length(files)) {
    stop("no DDInter CSVs in ", ddinter_dir, " (download the bulk files first)")
  }
  dt <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)
  # A pair repeats across each drug's category file; keep one row per unordered pair.
  dt[, `:=`(a = norm(Drug_A), b = norm(Drug_B))]
  unique(dt, by = c("a", "b", "Level"))
}

args <- commandArgs(trailingOnly = TRUE)

# --find mode: list DDInter drug names matching a substring, to reconcile naming.
if (length(args) >= 1 && args[1] == "--find") {
  query <- norm(substance(paste(args[-1], collapse = " ")))
  dt <- load_ddinter()
  names_all <- sort(unique(c(dt$Drug_A, dt$Drug_B)))
  hits <- names_all[grepl(query, norm(names_all), fixed = TRUE)]
  if (!length(hits)) cat("no DDInter drug name matches '", query, "'\n", sep = "") else cat(hits, sep = "\n")
  cat("\n")
  quit(save = "no")
}

if (length(args) < 2) {
  stop('usage: ddinter_lookup.R "<drug1>" "<drug2>"   |   ddinter_lookup.R --find "<name>"')
}

d1 <- norm(substance(args[1]))
d2 <- norm(substance(args[2]))
dt <- load_ddinter()

# Unordered pair match: {d1,d2} in either column order.
hit <- dt[(a == d1 & b == d2) | (a == d2 & b == d1)]

cat("drug1:", d1, "| drug2:", d2, "\n")
if (nrow(hit)) {
  cat("INTERACTION LISTED in DDInter:\n")
  for (i in seq_len(nrow(hit))) {
    cat("  ", hit$Drug_A[i], " x ", hit$Drug_B[i], " -> Level: ", hit$Level[i], "\n", sep = "")
  }
  cat("=> NOT a clean negative for this pair unless ruled out for THIS event vs label.\n")
} else {
  # Guard against a false 'absent' caused by a name mismatch.
  d1_seen <- d1 %in% c(dt$a, dt$b)
  d2_seen <- d2 %in% c(dt$a, dt$b)
  cat("NO INTERACTION LISTED for the pair in DDInter.\n")
  if (!d1_seen) cat("  warning: '", d1, "' not found as a DDInter drug name; check with --find before trusting absence.\n", sep = "")
  if (!d2_seen) cat("  warning: '", d2, "' not found as a DDInter drug name; check with --find before trusting absence.\n", sep = "")
}
