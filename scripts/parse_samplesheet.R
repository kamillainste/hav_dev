#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: <samplesheet.tsv> <tmp_dir>")
}

samplesheet_path <- args[1]
tmp_dir <- args[2]

dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

message("Reading samplesheet: ", samplesheet_path)

df <- read.delim(
  samplesheet_path,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

## ------------------------------------------------------------
## STRICT SCHEMA
## ------------------------------------------------------------

required_cols <- c(
  "Sample",
  "Virus",
  "MatchRef",
  "Segmented",
  "Primers",
  "Reference",
  "Features",
  "MinCov",
  "Mismatch",
  "InputDir"
)

missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

message("Schema validation: OK")

## ------------------------------------------------------------
## CLEAN + VALIDATE
## ------------------------------------------------------------

df$Sample   <- trimws(df$Sample)
df$InputDir <- trimws(df$InputDir)

if (any(is.na(df$Sample) | df$Sample == "")) {
  stop("Invalid Sample values detected")
}

if (any(is.na(df$InputDir) | df$InputDir == "")) {
  stop("Invalid InputDir values detected")
}

if (any(!dir.exists(df$InputDir))) {
  missing <- df$InputDir[!dir.exists(df$InputDir)]
  stop("Missing InputDir paths (max 5): ",
       paste(head(missing, 5), collapse = ", "))
}

if (any(duplicated(df$Sample))) {
  dup <- df$Sample[duplicated(df$Sample)]
  stop("Duplicate Sample IDs: ", paste(unique(dup), collapse = ", "))
}

message("Samples: ", nrow(df))
message("Unique samples: ", length(unique(df$Sample)))

## ------------------------------------------------------------
## OUTPUT 1: LOOP FILE (NanoPlot + Chopper ONLY)
## ------------------------------------------------------------

loop_file <- file.path(tmp_dir, "tmp_samplesheet.tsv")

write.table(
  df[, c("Sample", "InputDir")],
  file = loop_file,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

message("Loop file written: ", loop_file)

## ------------------------------------------------------------
## OUTPUT 2: VALIDATED SOURCE SHEET (FOR VC PIPELINE)
## ------------------------------------------------------------

vc_source_file <- file.path(tmp_dir, "vc_source.tsv")

write.table(
  df,
  file = vc_source_file,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

message("VC source sheet written: ", vc_source_file)

## ------------------------------------------------------------
## CONTRACT OUTPUT (machine-friendly)
## ------------------------------------------------------------

cat(
  paste0(
    "LOOP_FILE=", loop_file, "\n",
    "VC_SOURCE=", vc_source_file, "\n"
  )
)