library(tidyverse)
library(seqinr)

base <- file.path(getwd(), "data/local_datasets/2026-04-10")

# --- Process FASTA ---
fasta_in  <- file.path(base, "2PA.fa")
fasta_out <- file.path(base, "input.fa")

lines <- readLines(fasta_in, encoding = "latin1")

out_lines <- character(length(lines))
for (i in seq_along(lines)) {
  if (startsWith(lines[[i]], ">")) {
    parts <- str_split_fixed(sub("^>", "", lines[[i]]), "\\|", n = 3)
    out_lines[[i]] <- paste0(">", parts[1, 2])
  } else {
    out_lines[[i]] <- lines[[i]]
  }
}

writeLines(out_lines, fasta_out, useBytes = TRUE)
cat("FASTA written:", fasta_out, "\n")

# --- Process CSV to TSV ---
csv_in  <- file.path(base, "export.csv")
tsv_out <- file.path(base, "metadata.tsv")

read_delim(csv_in, delim = ";", locale = locale(encoding = "latin1"), show_col_types = FALSE) |>
  filter(!is.na(Key), Key != "") |>
  select(
    id       = Key,
    genotype = Genotype,
    lineage  = OUTBREAK_VARIANT,
    date     = `Sample date`,
    country  = Source
  ) |>
  mutate(
    date = suppressWarnings(
      format(as.Date(date, format = "%d.%m.%Y"), "%Y-%m-%d")
    )
  ) |>
  write_tsv(tsv_out)

cat("TSV written:", tsv_out, "\n")
