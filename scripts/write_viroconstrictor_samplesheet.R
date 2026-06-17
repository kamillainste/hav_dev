#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: <valid_samples.tsv> <input_samplesheet.tsv> <output_vc_samplesheet.tsv>")
}

valid_samples_path <- args[1]
full_samplesheet_path <- args[2]
out_path <- args[3]

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

message("Reading validated sample list: ", valid_samples_path)
valid <- read.delim(valid_samples_path, sep = "\t", header = FALSE,
                    col.names = c("Sample", "InputDir"),
                    stringsAsFactors = FALSE)

message("Reading full validated input samplesheet: ", full_samplesheet_path)
df <- read.delim(full_samplesheet_path,
                  sep = "\t",
                  header = TRUE,
                  stringsAsFactors = FALSE,
                  check.names = FALSE)

required_cols <- c(
  "Sample","Virus","MatchRef","Segmented","Primers",
  "Reference","Features","MinCov","Mismatch","InputDir"
)

missing <- setdiff(required_cols, colnames(df))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

df$Sample <- trimws(df$Sample)

valid_samples <- unique(valid$Sample)

df <- df[df$Sample %in% valid_samples, ]

if (nrow(df) == 0) {
  stop("No overlap between validated samples and samplesheet")
}

# Rename columns for the ViroConstrictor samplesheet
vc <- df

names(vc)[names(vc) == "MatchRef"] <- "Match-ref"
names(vc)[names(vc) == "MinCov"] <- "min-coverage"
names(vc)[names(vc) == "Mismatch"] <- "primer-mismatch-rate"

vc$InputDir <- NULL

required_final_cols <- c(
  "Sample",
  "Virus",
  "Match-ref",
  "Segmented",
  "Primers",
  "Reference",
  "Features",
  "min-coverage",
  "primer-mismatch-rate"
)

missing <- setdiff(required_final_cols, colnames(vc))
if (length(missing) > 0) {
  stop("VC schema mismatch: ", paste(missing, collapse = ", "))
}

vc <- vc[, required_final_cols]

write.table(
  vc,
  file = out_path,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

message("ViroConstrictor samplesheet written: ", out_path)
message("Samples included: ", nrow(vc))
message("Output path: ", out_path)