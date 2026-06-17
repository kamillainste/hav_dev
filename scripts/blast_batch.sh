#!/usr/bin/env bash
# scripts/blast_batch.sh
#
# BLASTs a batch FASTA of query sequences against the local HAV BLAST database
# built by make_blast_db.sh.
#
# Run from the project root with HAVDEV conda active:
#   conda activate /path/to/hav_dev/.conda/HAVDEV
#   bash scripts/blast_batch.sh <batch_fasta> [dataset_date]
#
# Arguments:
#   batch_fasta    Path to the multi-FASTA file (e.g. data/Batch-1/Batch-1.fa)
#   dataset_date   Dataset date folder (default: 2026-04-10)
#
# Output columns (tab-separated):
#   qseqid    query sequence ID
#   sseqid    subject (database hit) sequence ID
#   pident    % identity
#   length    alignment length
#   mismatch  number of mismatches (SNPs)
#   gapopen   number of gap-opening events
#   qstart    query alignment start
#   qend      query alignment end
#   sstart    subject alignment start
#   send      subject alignment end
#   evalue    E-value
#   bitscore  bit score
#   stitle    full FASTA header of hit (includes ID)
#   qlen      query sequence length
#   slen      subject sequence length
#   nident    number of identical positions
#   gaps      total number of gap characters
#   qcovs     % of query sequence covered by alignment
#
# All hits with e-value <= 1e-5 are reported (one best HSP per subject).
# Join with metadata_corrected.tsv on sseqid == id for genotype/lineage etc.
#
# To retrieve the full sequence for a hit, use the deduplicated FASTA:
#   grep -A2 '^><sseqid>' data/local_datasets/<date>/blast_db/input_dedup.fa
#   seqkit grep -p <sseqid> data/local_datasets/<date>/blast_db/input_dedup.fa

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
BATCH_FA="${1:?Usage: bash scripts/blast_batch.sh <batch_fasta> [dataset_date]}"
DATASET_DATE="${2:-2026-04-10}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DB_PATH="data/local_datasets/$DATASET_DATE/blast_db/hav"
THREADS=4
EVALUE="1e-5"

# Normalise batch FASTA to a project-root-relative path
if [[ "$BATCH_FA" == /* ]]; then
    # strip absolute PROJECT_DIR prefix if present
    BATCH_FA="${BATCH_FA#"$PROJECT_DIR/"}"
fi

# Derive output path alongside the batch FASTA
BATCH_DIR="$(dirname "$BATCH_FA")"
BATCH_NAME="$(basename "$BATCH_FA" .fa)"
OUT_TSV="$BATCH_DIR/${BATCH_NAME}_blast_results.tsv"

cd "$PROJECT_DIR"

# ── Checks ────────────────────────────────────────────────────────────────────
if [[ ! -f "$BATCH_FA" ]]; then
    echo "ERROR: batch FASTA not found: $PROJECT_DIR/$BATCH_FA" >&2
    exit 1
fi

if [[ ! -f "${DB_PATH}.nhr" && ! -f "${DB_PATH}.ndb" ]]; then
    echo "ERROR: BLAST database not found at $PROJECT_DIR/$DB_PATH" >&2
    echo "  Build it first: bash scripts/make_blast_db.sh $DATASET_DATE" >&2
    exit 1
fi

if ! command -v blastn &> /dev/null; then
    echo "ERROR: blastn not found. Activate the HAVDEV conda environment:" >&2
    echo "  conda activate /path/to/hav_dev/.conda/HAVDEV" >&2
    exit 1
fi

# ── Run BLAST ─────────────────────────────────────────────────────────────────
echo "── Running blastn ────────────────────────────────────────────────────────"
echo "  Query : $PROJECT_DIR/$BATCH_FA"
echo "  DB    : $PROJECT_DIR/$DB_PATH"
echo "  Output: $PROJECT_DIR/$OUT_TSV"
echo ""

# Write header
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    qseqid sseqid pident length mismatch gapopen \
    qstart qend sstart send evalue bitscore \
    stitle qlen slen nident gaps qcovs \
    > "$OUT_TSV"

blastn \
    -query       "$BATCH_FA" \
    -db          "$DB_PATH" \
    -out         "$OUT_TSV.tmp" \
    -evalue      "$EVALUE" \
    -max_hsps    1 \
    -num_threads "$THREADS" \
    -outfmt      "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen nident gaps qcovs"

cat "$OUT_TSV.tmp" >> "$OUT_TSV"
rm  "$OUT_TSV.tmp"

# ── Summary ───────────────────────────────────────────────────────────────────
N_HITS=$(tail -n +2 "$OUT_TSV" | wc -l)
N_QUERIES=$(grep -c "^>" "$BATCH_FA")

echo "── Done ──────────────────────────────────────────────────────────────────"
echo "  Queries  : $N_QUERIES"
echo "  Hits     : $N_HITS (e-value <= $EVALUE)"
echo "  Output   : $PROJECT_DIR/$OUT_TSV"
echo ""
echo "  To retrieve a hit sequence (from deduplicated FASTA):"
echo "    grep -A2 '^><sseqid>' data/local_datasets/$DATASET_DATE/blast_db/input_dedup.fa"
  echo ""
  echo "  To join with metadata:"
  echo "    data/local_datasets/$DATASET_DATE/metadata_corrected.tsv  (join on sseqid == id)"
