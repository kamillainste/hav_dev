#!/usr/bin/env bash
# scripts/combine_sequences.sh
#
# Combines internal HAV sequences with external sequences from the SQLite
# database (data/local_datasets/externe_sek/externe.db) into a single
# combined.fa and combined_metadata.tsv for use in preparing the local BLAST
# dataset.
#
# Run after prepare_dataset.R has produced input.fa and metadata.tsv.
#
# Usage (from project root in WSL with HAVDEV conda active):
#   bash scripts/combine_sequences.sh <dataset_date>
#
# Example:
#   bash scripts/combine_sequences.sh 2026-04-10
#
# Inputs:
#   data/local_datasets/<date>/input.fa       -- internal sequences (prepare_dataset.R)
#   data/local_datasets/<date>/metadata.tsv   -- internal metadata (prepare_dataset.R)
#   data/local_datasets/externe_sek/externe.db -- external sequences SQLite database
#
# Outputs:
#   data/local_datasets/<date>/combined.fa           -- all sequences merged
#   data/local_datasets/<date>/combined_metadata.tsv -- all metadata merged

set -euo pipefail

DATASET_DATE="${1:?Usage: bash scripts/combine_sequences.sh <dataset_date>}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BASE_DIR="$PROJECT_DIR/data/local_datasets/$DATASET_DATE"
EXT_DB="$PROJECT_DIR/data/local_datasets/externe_sek/externe.db"
COMBINED_FA="$BASE_DIR/combined.fa"
COMBINED_META="$BASE_DIR/combined_metadata.tsv"

cd "$PROJECT_DIR"

# ── Checks ────────────────────────────────────────────────────────────────────
if [[ ! -f "$BASE_DIR/input.fa" ]]; then
  echo "ERROR: $BASE_DIR/input.fa not found." >&2
  echo "  Run prepare_dataset.R first: Rscript scripts/prepare_dataset.R $DATASET_DATE" >&2
  exit 1
fi

if [[ ! -f "$BASE_DIR/metadata.tsv" ]]; then
  echo "ERROR: $BASE_DIR/metadata.tsv not found." >&2
  echo "  Run prepare_dataset.R first: Rscript scripts/prepare_dataset.R $DATASET_DATE" >&2
  exit 1
fi

if [[ ! -f "$EXT_DB" ]]; then
  echo "ERROR: External sequence database not found: $EXT_DB" >&2
  exit 1
fi

echo "════════════════════════════════════════════════════════"
echo " Combining internal + external sequences — $DATASET_DATE"
echo "════════════════════════════════════════════════════════"

# ── Export external sequences from SQLite ─────────────────────────────────────
echo ""
echo "── Exporting external sequences from database ───────────────────────────"

TODAY=$(date +%Y-%m-%d)
EXT_FA="$BASE_DIR/ext_seqs_${TODAY}.fasta"
EXT_META="$BASE_DIR/ext_seqs_metadata_${TODAY}.tsv"

"$HOME/.conda/R_shared/bin/Rscript" - "$EXT_DB" "$EXT_FA" "$EXT_META" << 'REOF'
args <- commandArgs(trailingOnly = TRUE)
db_file   <- args[1]
fa_out    <- args[2]
meta_out  <- args[3]

suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
  library(readr)
})

con <- dbConnect(SQLite(), dbname = db_file)
df  <- dbGetQuery(con, "SELECT * FROM ext_seqs ORDER BY sequence_id")
dbDisconnect(con)

if (nrow(df) == 0) {
  cat("WARNING: No sequences in external database\n")
  file.create(fa_out)
  write_tsv(data.frame(), meta_out)
  quit(status = 0)
}

# Write FASTA
fasta_conn <- file(fa_out, open = "w")
for (i in seq_len(nrow(df))) {
  writeLines(paste0(">", df$sequence_id[i]), fasta_conn)
  writeLines(df$sequence[i], fasta_conn)
}
close(fasta_conn)

# Write metadata (standardise column names to match internal metadata.tsv)
meta <- data.frame(
  id      = df$sequence_id,
  genotype = df$genotype,
  lineage  = df$variant,
  date     = df$collection_date,
  country  = NA_character_,
  stringsAsFactors = FALSE
)
write_tsv(meta, meta_out)

cat(sprintf("Exported %d external sequences\n", nrow(df)))
cat(sprintf("  FASTA:    %s\n", fa_out))
cat(sprintf("  Metadata: %s\n", meta_out))
REOF

N_INT=$(grep -c "^>" "$BASE_DIR/input.fa")
N_EXT=$(grep -c "^>" "$EXT_FA" 2>/dev/null || echo 0)
echo "  Internal sequences: $N_INT"
echo "  External sequences: $N_EXT"

# ── Merge FASTA ───────────────────────────────────────────────────────────────
echo ""
echo "── Merging FASTA files ──────────────────────────────────────────────────"
cat "$BASE_DIR/input.fa" "$EXT_FA" > "$COMBINED_FA"
N_TOTAL=$(grep -c "^>" "$COMBINED_FA")
echo "  Combined: $N_TOTAL sequences → $COMBINED_FA"

# ── Merge metadata ────────────────────────────────────────────────────────────
echo ""
echo "── Merging metadata files ───────────────────────────────────────────────"
# Combine TSVs: internal metadata + external metadata (drop duplicate header)
head -n 1 "$BASE_DIR/metadata.tsv" > "$COMBINED_META"
tail -n +2 "$BASE_DIR/metadata.tsv" >> "$COMBINED_META"
tail -n +2 "$EXT_META" >> "$COMBINED_META"
N_META=$(tail -n +2 "$COMBINED_META" | wc -l)
echo "  Combined: $N_META metadata rows → $COMBINED_META"

echo ""
echo "════════════════════════════════════════════════════════"
echo " Done. Combined dataset ready at:"
echo "   $COMBINED_FA"
echo "   $COMBINED_META"
echo ""
echo " Next step — build the BLAST database:"
echo "   bash scripts/make_blast_db.sh $DATASET_DATE"
echo "════════════════════════════════════════════════════════"
