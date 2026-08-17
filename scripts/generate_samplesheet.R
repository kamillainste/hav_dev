#!/usr/bin/env Rscript
#
# scripts/generate_samplesheet.R
#
# Generate a samplesheet TSV from FASTA filenames in a directory.
#
# FASTA filenames must follow the format: LabWareID_SekvensID-HAV.fa(sta)
#
# Usage:
#   Rscript scripts/generate_samplesheet.R <fasta_dir> <output_samplesheet>
#
# Arguments:
#   <fasta_dir>           Directory containing .fa/.fasta files
#   <output_samplesheet>  Path to output TSV file (created with Sample, LabWareID columns)
#
# Output:
#   TSV file with columns:
#   - Sample:     LabWareID_SekvensID-HAV (complete filename without extension)
#   - LabWareID:  LabWareID (extracted from filename)
#   - SekvensID:  SekvensID-HAV (extracted from filename)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript scripts/generate_samplesheet.R <fasta_dir> <output_samplesheet>\n", file = stderr())
  quit(status = 1)
}

fasta_dir <- args[1]
output_file <- args[2]

# Validate input directory
if (!dir.exists(fasta_dir)) {
  cat(sprintf("ERROR: FASTA directory not found: %s\n", fasta_dir), file = stderr())
  quit(status = 1)
}

# Find all FASTA files
fasta_files <- list.files(
  fasta_dir,
  pattern = "\\.(fa|fasta)$",
  full.names = FALSE,
  ignore.case = TRUE
)

if (length(fasta_files) == 0) {
  cat(sprintf("ERROR: No .fa or .fasta files found in %s\n", fasta_dir), file = stderr())
  quit(status = 1)
}

# Parse filenames: LabWareID_SekvensID-HAV.fa(sta)
samples <- data.frame(
  Sample = character(length(fasta_files)),
  LabWareID = character(length(fasta_files)),
  SekvensID = character(length(fasta_files)),
  stringsAsFactors = FALSE
)

for (i in seq_along(fasta_files)) {
  basename <- fasta_files[i]
  
  # Remove extension
  basename_no_ext <- sub("\\.(fa|fasta)$", "", basename, ignore.case = TRUE)
  
  # Split on first underscore: LabWareID_SekvensID-HAV
  underscore_pos <- regexpr("_", basename_no_ext, fixed = TRUE)
  
  if (underscore_pos == -1) {
    # No underscore at all: whole name is LabWareID, SekvensID missing
    labware_id <- basename_no_ext
    sequence_id <- ""
  } else {
    labware_id <- substr(basename_no_ext, 1, underscore_pos - 1)
    sequence_id <- substr(basename_no_ext, underscore_pos + 1, nchar(basename_no_ext))
  }
  
  # Missing LabWareID or SekvensID (e.g. "Posktr_.fa") defaults to "0"
  if (labware_id == "") {
    cat(sprintf("WARNING: Missing LabWareID in filename: %s — defaulting to '0'\n",
                fasta_files[i]), file = stderr())
    labware_id <- "0"
  }
  if (sequence_id == "") {
    cat(sprintf("WARNING: Missing SekvensID in filename: %s — defaulting to '0'\n",
                fasta_files[i]), file = stderr())
    sequence_id <- "0"
  }
  
  # Sample should be the full filename (without extension) so prepare_input_fasta.R can find it
  samples$Sample[i] <- basename_no_ext
  samples$LabWareID[i] <- labware_id
  samples$SekvensID[i] <- sequence_id
}

# Remove empty rows (from parsing failures)
samples <- samples[samples$Sample != "", ]

if (nrow(samples) == 0) {
  cat(sprintf("ERROR: No valid FASTA files parsed from %s\n", fasta_dir), file = stderr())
  quit(status = 1)
}

# Write samplesheet
write.table(
  samples,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(sprintf("Generated samplesheet with %d samples: %s\n", nrow(samples), output_file))
quit(status = 0)
