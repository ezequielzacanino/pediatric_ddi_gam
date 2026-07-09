################################################################################
# Dataset generation
# Script 04_build_ade_raw
################################################################################

# libraries
pacman::p_load(tidyverse, data.table, mgcv, Rcpp, doParallel, DBI, RSQLite)

################################################################################
# Configuration
################################################################################

seed = 0
set.seed(seed)

script_args <- commandArgs(trailingOnly = FALSE)
script_path_arg <- "--file="
script_file_arg <- script_args[grep(script_path_arg, script_args)]
script_path <- dirname(normalizePath(sub(script_path_arg, "", script_file_arg)))
has_script_file <- length(script_path) > 0 && !is.na(script_path)
if (!has_script_file) {
  script_path <- getwd()
}

# Script lives in scripts/R cuando when run by Rscript
project_dir <- if (has_script_file) {
  normalizePath(file.path(script_path, "..", ".."), mustWork = TRUE)
} else {
  normalizePath(script_path, mustWork = TRUE)
}

cores <- as.integer(
  Sys.getenv(
    "FAERS_SIM_CORES",
    unset = max(1, min(50, parallel::detectCores(logical = FALSE)))
  )
)
registerDoParallel(cores=cores)

raw_data_dir <- file.path(project_dir, "data", "raw")
processed_data_dir <- file.path(project_dir, "data", "processed")
dir.create(processed_data_dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
# ADE pre-processing
################################################################################

ade_input_path <- file.path(
  raw_data_dir,
  "pediatric_patients_report_serious_reporter_drugs_reactions.csv.gz"
)
ade_output_path <- file.path(processed_data_dir, "ade_raw.csv")
ade_sqlite_path <- tempfile(
  pattern = "ade_raw_build_",
  tmpdir = tempdir(),
  fileext = ".sqlite"
)
ade_chunk_lines <- as.integer(
  Sys.getenv("FAERS_ADE_CHUNK_LINES", unset = "100000")
)
ade_export_chunk_rows <- as.integer(
  Sys.getenv("FAERS_ADE_EXPORT_CHUNK_ROWS", unset = "1000000")
)

ade_con <- gzfile(ade_input_path, open = "rt")
ade_header <- readLines(ade_con, n = 1, warn = FALSE)

sqlite_con <- DBI::dbConnect(RSQLite::SQLite(), ade_sqlite_path)

# 1- Project uses SQLite, main file is too big to be handled in the memory
invisible(DBI::dbExecute(sqlite_con, "PRAGMA journal_mode = WAL;"))
invisible(DBI::dbExecute(sqlite_con, "PRAGMA synchronous = OFF;"))
invisible(DBI::dbExecute(sqlite_con, "PRAGMA temp_store = MEMORY;"))
invisible(DBI::dbExecute(sqlite_con, "PRAGMA cache_size = -200000;"))

invisible(DBI::dbExecute(
  sqlite_con,
  paste(
    "CREATE TABLE ade_raw (",
    "row_key TEXT PRIMARY KEY,",
    "safetyreportid TEXT NOT NULL,",
    "ade TEXT NOT NULL,",
    "atc_concept_id INTEGER NOT NULL,",
    "meddra_concept_id INTEGER NOT NULL,",
    "nichd TEXT,",
    "sex TEXT,",
    "reporter_qualification TEXT,",
    "receive_date TEXT",
    ");"
  )
))

ade_chunk_index <- 0L
ade_rows_read <- 0L

repeat {
  ade_lines <- readLines(ade_con, n = ade_chunk_lines, warn = FALSE)

  if (length(ade_lines) == 0) {
    break
  }

  ade_chunk_index <- ade_chunk_index + 1L
  ade_rows_read <- ade_rows_read + length(ade_lines)

  # 2- Only reads neccesary columns for building ade_raw.
  ade_chunk <- fread(
    text = paste(c(ade_header, ade_lines), collapse = "\n"),
    select = c(
      "safetyreportid",
      "ATC_concept_id",
      "MedDRA_concept_id",
      "nichd",
      "patient_sex",
      "reporter_qualification",
      "receive_date"
    ),
    showProgress = FALSE
  )

  # 3- Normalizes types and deletes raws without standarized id.
  ade_chunk <- unique(
    ade_chunk[
      ,
      .(
        # 3a- deletes "-x" from id.
        safetyreportid = sub("-[^-]+$", "", as.character(safetyreportid)),
        ade = paste0(
          as.integer(ATC_concept_id)
        ),
        atc_concept_id = as.integer(ATC_concept_id),
        meddra_concept_id = as.integer(MedDRA_concept_id),
        nichd = fifelse(is.na(nichd), "", as.character(nichd)),
        sex = fifelse(is.na(patient_sex), "", as.character(patient_sex)),
        reporter_qualification = fifelse(
          is.na(reporter_qualification),
          "",
          as.character(reporter_qualification)
        ),
        receive_date = fifelse(
          is.na(receive_date),
          "",
          as.character(receive_date)
        )
      )
    ][
      !is.na(atc_concept_id) & !is.na(meddra_concept_id)
    ]
  )

  # 4- The key combines every output column to de-duplicate between blocks with INSERT OR IGNORE.
  # meddra_concept_id must be part of the key; otherwise distinct reactions for the same report-drug collapse into a single row.
  ade_chunk[
    ,
    row_key := paste(
      safetyreportid,
      ade,
      meddra_concept_id,
      nichd,
      sex,
      reporter_qualification,
      receive_date,
      sep = "|"
    )
  ]

  setcolorder(
    ade_chunk,
    c(
      "row_key",
      "safetyreportid",
      "ade",
      "atc_concept_id",
      "meddra_concept_id",
      "nichd",
      "sex",
      "reporter_qualification",
      "receive_date"
    )
  )

  DBI::dbWriteTable(
    sqlite_con,
    "ade_raw_stage",
    as.data.frame(ade_chunk),
    overwrite = TRUE
  )

  invisible(DBI::dbExecute(
    sqlite_con,
    paste(
      "INSERT OR IGNORE INTO ade_raw",
      "(row_key, safetyreportid, ade, atc_concept_id, meddra_concept_id,",
      "nichd, sex, reporter_qualification, receive_date)",
      "SELECT row_key, safetyreportid, ade, atc_concept_id, meddra_concept_id,",
      "nichd, sex, reporter_qualification, receive_date",
      "FROM ade_raw_stage;"
    )
  ))

  invisible(DBI::dbExecute(sqlite_con, "DROP TABLE ade_raw_stage;"))

  cat(
    sprintf(
      "Chunk %d procesado | filas leÃ­das acumuladas: %s\n",
      ade_chunk_index,
      format(ade_rows_read, big.mark = ",", scientific = FALSE)
    )
  )

  rm(ade_chunk, ade_lines)
  gc(verbose = FALSE)
}

close(ade_con)

ade_total_rows <- DBI::dbGetQuery(
  sqlite_con,
  "SELECT COUNT(*) AS n FROM ade_raw;"
)$n[[1]]

ade_result <- DBI::dbSendQuery(
  sqlite_con,
  paste(
    "SELECT",
    "safetyreportid,",
    "ade,",
    "atc_concept_id,",
    "meddra_concept_id,",
    "NULLIF(nichd, '') AS nichd,",
    "NULLIF(sex, '') AS sex,",
    "NULLIF(reporter_qualification, '') AS reporter_qualification,",
    "NULLIF(receive_date, '') AS receive_date",
    "FROM ade_raw;"
  )
)

ade_first_write <- TRUE

repeat {
  ade_export_chunk <- DBI::dbFetch(ade_result, n = ade_export_chunk_rows)

  if (nrow(ade_export_chunk) == 0) {
    break
  }

  fwrite(
    as.data.table(ade_export_chunk),
    ade_output_path,
    append = !ade_first_write,
    col.names = ade_first_write
  )

  ade_first_write <- FALSE
  rm(ade_export_chunk)
  gc(verbose = FALSE)
}

invisible(DBI::dbClearResult(ade_result))
invisible(DBI::dbDisconnect(sqlite_con))

if (file.exists(ade_sqlite_path)) {
  invisible(file.remove(ade_sqlite_path))
}

cat("\nade_raw generado en:\n")
cat(ade_output_path, "\n")
cat("Filas:", format(ade_total_rows, big.mark = ",", scientific = FALSE), "\n")
cat("Columnas: 8\n")

quit(save = "no", status = 0)
