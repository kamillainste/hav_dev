#!/usr/bin/env Rscript
#
# Reads a samplesheet TSV and collects one consensus sequence per sample
# into a single multi-FASTA output file.
#
# Samplesheet format (TSV, header required):
#   Sample    — sequence ID used in the output FASTA header
#   InputDir  — directory containing the consensus file for that sample
#               (.fa, .fasta, or .TXT; first match is used)
#
# For .TXT files (raw Sanger contig text): all lines are concatenated
# as sequence.
# For .fa / .fasta files: the first sequence in the file is taken; the Sample
#   column overrides whatever header is in the file.
#
# Usage:
#   Rscript scripts/prepare_input_fasta.R <samplesheet.tsv> <output.fa>

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(seqinr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript scripts/prepare_input_fasta.R <samplesheet.tsv> <output.fa>",
    call. = FALSE
  )
}

samplesheet_path <- args[1]
output_fa        <- args[2]

# ── Read & validate samplesheet ───────────────────────────────────────────────
if (!file.exists(samplesheet_path)) {
  stop("Samplesheet not found: ", samplesheet_path, call. = FALSE)
}

df <- read_tsv(samplesheet_path, show_col_types = FALSE)

required_cols <- c("Sample", "InputDir")
missing_cols  <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Samplesheet missing required columns: ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}

df$Sample   <- trimws(df$Sample)
df$InputDir <- trimws(df$InputDir)

if (any(is.na(df$Sample) | df$Sample == "")) {
  stop("Empty Sample values detected in samplesheet", call. = FALSE)
}
if (any(duplicated(df$Sample))) {
  dups <- unique(df$Sample[duplicated(df$Sample)])
  stop("Duplicate Sample IDs: ", paste(dups, collapse = ", "), call. = FALSE)
}

message("Samples: ", nrow(df))

# ── Read one sequence per sample ──────────────────────────────────────────────
# Returns list(name, seq) or stops with an informative error.
read_sample_seq <- function(sample_name, input_dir) {
  if (!dir.exists(input_dir)) {
    stop("InputDir not found for sample '", sample_name, "': ", input_dir,
         call. = FALSE)
  }

  # First token only in FASTA headers (tools truncate at whitespace)
  clean_name <- str_replace(sample_name, "\\s+.*$", "")

  txt_files   <- list.files(input_dir, pattern = "\\.TXT$",
                            full.names = TRUE, ignore.case = FALSE)
  fasta_files <- list.files(input_dir, pattern = "\\.(fa|fasta)$",
                            full.names = TRUE, ignore.case = TRUE)

  if (length(txt_files) > 0) {
    raw <- readLines(txt_files[[1]])
    seq <- paste(str_remove_all(raw, "\\s+"), collapse = "")
    if (nchar(seq) == 0) {
      stop("Empty sequence in TXT file for sample '", sample_name, "'",
           call. = FALSE)
    }
    return(list(name = clean_name, seq = seq))
  }

  if (length(fasta_files) > 0) {
    seqs <- read.fasta(fasta_files[[1]], seqtype = "DNA",
                       as.string = TRUE, forceDNAtolower = FALSE)
    if (length(seqs) == 0) {
      stop("No sequences in FASTA for sample '", sample_name, "': ",
           fasta_files[[1]], call. = FALSE)
    }
    # Use Sample name from samplesheet, not the embedded FASTA header
    return(list(name = clean_name, seq = as.character(seqs[[1]])))
  }

  stop("No .fa, .fasta, or .TXT file found in InputDir for sample '",
       sample_name, "': ", input_dir, call. = FALSE)
}

all_seqs <- vector("list", nrow(df))
for (i in seq_len(nrow(df))) {
  all_seqs[[i]] <- read_sample_seq(df$Sample[i], df$InputDir[i])
}

# ── Write multi-FASTA ─────────────────────────────────────────────────────────
out_dir <- dirname(output_fa)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

fasta_lines <- unlist(lapply(
  all_seqs, function(s) c(paste0(">", s$name), s$seq)
))
writeLines(fasta_lines, output_fa)

message("Written ", length(all_seqs), " sequence(s) to: ", output_fa)
