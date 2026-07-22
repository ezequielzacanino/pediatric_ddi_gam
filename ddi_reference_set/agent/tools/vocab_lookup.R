# Token-efficient vocabulary lookup for the curation agent.
#
# The workbook holds the controlled vocabularies as reference sheets used for the
# dropdowns (ref_atc, ref_llt, ref_pt, ref_hlt, ref_hlgt, ref_nichd). They total
# ~78k terms, so they must not be loaded into the agent context wholesale. This
# helper searches one level and prints only the matching dropdown strings, so a
# candidate drug/event can be resolved to its exact workbook value cheaply.
#
# Run from the subproject root:
#   Rscript scripts/R/vocab_lookup.R <level> <query> [--exact] [--max N]
#
#   level : atc | llt | pt | hlt | hlgt | nichd
#   query : text to match (case/space/underscore-insensitive substring)
#   --exact : report whether <query> is an exact dropdown value (MATCH/NO MATCH),
#             to validate a term before writing it to the workbook
#   --max N : cap the number of printed matches (default 50)

source("00_functions.R", local = TRUE)

workbook_path <- Sys.getenv("DDI_WORKBOOK", "input/ddi_reference_input.xlsx")

level_to_sheet <- c(
  atc   = "ref_atc",
  llt   = "ref_llt",
  pt    = "ref_pt",
  hlt   = "ref_hlt",
  hlgt  = "ref_hlgt",
  nichd = "ref_nichd"
)

# 1. Parse arguments: positional <level> <query>, plus optional --exact / --max N.
args <- commandArgs(trailingOnly = TRUE)
exact <- "--exact" %in% args
max_hits <- 50L
if ("--max" %in% args) {
  max_hits <- as.integer(args[which(args == "--max") + 1L])
}
positional <- args[!grepl("^--", args)]
positional <- positional[!positional %in% as.character(max_hits)]
level <- positional[1]
query <- paste(positional[-1], collapse = " ")

if (is.na(level) || !level %in% names(level_to_sheet) || !nzchar(query)) {
  stop("usage: vocab_lookup.R <atc|llt|pt|hlt|hlgt|nichd> <query> [--exact] [--max N]")
}

# 2. Read the requested reference sheet and normalize terms the same way the
#    curation script does, so matches respect its case/space/underscore rules.
terms <- read_workbook_sheet(workbook_path, level_to_sheet[[level]])[["term"]]
terms <- terms[!is.na(terms)]
key <- normalize_vocabulary_key(terms)
query_key <- normalize_vocabulary_key(query)

# 3a. Exact mode: validate a single candidate value against the dropdown.
if (exact) {
  hit <- terms[!is.na(key) & key == query_key]
  if (length(hit)) {
    cat("MATCH:", hit[1], "\n")
  } else {
    cat("NO MATCH\n")
  }
  quit(save = "no")
}

# 3b. Search mode: print matching dropdown strings only, capped at --max.
matches <- terms[!is.na(key) & grepl(query_key, key, fixed = TRUE)]
if (!length(matches)) {
  cat("no matches for '", query, "' in ", level_to_sheet[[level]], "\n", sep = "")
  quit(save = "no")
}
cat(matches[seq_len(min(length(matches), max_hits))], sep = "\n")
cat("\n")
if (length(matches) > max_hits) {
  cat("... (", length(matches) - max_hits, " more; refine the query or raise --max)\n", sep = "")
}
