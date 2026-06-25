#!/usr/bin/env Rscript
# Debug script to verify file paths in batch_report.Rmd
# 
# Usage: Rscript scripts/check_paths.R <batch_dir> <dataset_date> <batch_fa> <dataset_dir> <out_base>
#
# These are the SAME arguments that run_all_analyses.sh passes

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  cat("\n")
  cat("Usage: Rscript scripts/check_paths.R <batch_dir> <dataset_date> <batch_fa> <dataset_dir> <out_base>\n\n")
  cat("Example:\n")
  cat("  Rscript scripts/check_paths.R data/Batch-1 2026-04-10 data/Batch-1/batch.fa data/local_datasets/2026-04-10 data/Batch-1/output\n\n")
  quit(status = 1)
}

batch_dir_param <- args[1]
dataset_date_param <- args[2]
batch_fa_param <- args[3]
dataset_dir_param <- args[4]
out_base_param <- args[5]

cat("\n")
cat(paste0("=", strrep("=", 78), "\n"))
cat("FILE PATH DEBUG CHECK\n")
cat(paste0("=", strrep("=", 78), "\n\n"))

# Print input parameters
cat("INPUT PARAMETERS (from run_all_analyses.sh):\n")
cat("  batch_dir   :", batch_dir_param, "\n")
cat("  dataset_date:", dataset_date_param, "\n")
cat("  batch_fa    :", batch_fa_param, "\n")
cat("  dataset_dir :", dataset_dir_param, "\n")
cat("  out_base    :", out_base_param, "\n\n")

# Construct paths (same as batch_report.Rmd)
batch_dir    <- file.path(batch_dir_param)
out_dir      <- file.path(out_base_param)
trees_dir    <- file.path(out_dir, "trees")
dataset_base <- file.path(dataset_dir_param)
metadata_tsv <- file.path(dataset_base, "metadata_corrected.tsv")
batch_name   <- basename(batch_dir_param)
input_fasta  <- file.path(batch_dir, paste0(batch_name, ".fa"))

# Additional derived paths that will be used
lineages_dir <- file.path(out_dir, "lineages")
nextclade_tsv <- file.path(lineages_dir, "nextclade.tsv")
blast_results_tsv <- file.path(out_dir, "blast_results.tsv")

# Create a path check function
check_path <- function(name, path, should_exist = NA) {
  exists_flag <- "     "
  color <- ""
  
  if (!is.na(should_exist)) {
    path_exists <- file.exists(path)
    if (should_exist && !path_exists) {
      exists_flag <- "[✗] "
    } else if (should_exist && path_exists) {
      exists_flag <- "[✓] "
    } else if (!should_exist && path_exists) {
      exists_flag <- "[⚠] "
    } else {
      exists_flag <- "[?] "
    }
  }
  
  cat(sprintf("%s%-25s %s\n", exists_flag, name, path))
}

# Print all paths
cat("CONSTRUCTED PATHS:\n")
cat(paste0("─", strrep("─", 77), "\n"))

check_path("batch_dir", batch_dir, TRUE)
check_path("batch_fa", batch_fa_param, TRUE)
check_path("batch_name", batch_name)

cat("\n")

check_path("dataset_base", dataset_base, TRUE)
check_path("metadata_tsv", metadata_tsv, TRUE)

cat("\n")

check_path("out_dir", out_dir, TRUE)
check_path("trees_dir", trees_dir, TRUE)
check_path("lineages_dir", lineages_dir, TRUE)
check_path("blast_results_tsv", blast_results_tsv, TRUE)
check_path("nextclade_tsv", nextclade_tsv, TRUE)

cat("\n")
cat(paste0("─", strrep("─", 77), "\n"))

# Legend
cat("\nLEGEND:\n")
cat("  [✓]  File/directory exists (as expected)\n")
cat("  [✗]  File/directory MISSING (should exist for report)\n")
cat("  [⚠]  File/directory exists but shouldn't yet\n")
cat("  [?]  Status uncertain\n")
cat("  [  ]  No check performed\n\n")

# Test normalizePath
cat("PATH NORMALIZATION TEST:\n")
cat(paste0("─", strrep("─", 77), "\n"))
cat(sprintf("  batch_dir normalized: %s\n", normalizePath(batch_dir, mustWork = FALSE)))
cat(sprintf("  dataset_base normalized: %s\n", normalizePath(dataset_base, mustWork = FALSE)))
cat(sprintf("  out_dir normalized: %s\n", normalizePath(out_dir, mustWork = FALSE)))

cat(paste0("\n", strrep("=", 80), "\n\n"))
