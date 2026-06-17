#!/usr/bin/env bash
# scripts/update_datasets.sh
#
# Refresh the local HAV dataset after adding new sequences to the source
# database. The local dataset has a single purpose: it is the reference
# collection that backs the BLAST database used to identify batch query
# sequences. There are no clade-separated datasets and no local reference
# trees — lineage assignment is handled entirely by the public-derived
# hav-vp1-2b-lineages Nextclade dataset (see build_vp1b_lineage_dataset.sh).
#
# ┌─────────────────────────────────────────────────────────────────────────┐
# │  WORKFLOW                                                                │
# │  1. Export from the sequence database into                              │
# │       data/local_datasets/<YYYY-MM-DD>/                                 │
# │         • export.csv  — metadata (Key, Genotype, OUTBREAK_VARIANT, …)   │
# │         • 2PA.fa      — all sequences in FASTA format                   │
# │  2. Run:  bash scripts/update_datasets.sh <YYYY-MM-DD>                  │
# │     It will:                                                            │
# │       a. Prepare input.fa + metadata.tsv      (prepare_dataset.R)       │
# │       b. Ensure metadata_corrected.tsv exists (for BLAST annotation)    │
# │       c. Build the BLAST database             (make_blast_db.sh)        │
# └─────────────────────────────────────────────────────────────────────────┘
#
# Prerequisites:
#   Run from the project root with HAVDEV conda active:
#     conda activate /path/to/hav_dev/.conda/HAVDEV
#     bash scripts/update_datasets.sh <YYYY-MM-DD>

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
DATASET_DATE="${1:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

if [[ -z "$DATASET_DATE" ]]; then
  echo "Usage: bash scripts/update_datasets.sh <YYYY-MM-DD>"
  echo "Example: bash scripts/update_datasets.sh $(date +%Y-%m-%d)"
  exit 1
fi

BASE_DIR="$PROJECT_DIR/data/local_datasets/$DATASET_DATE"

cd "$PROJECT_DIR"

echo "════════════════════════════════════════════════════════"
echo " Updating local BLAST dataset — $DATASET_DATE"
echo "════════════════════════════════════════════════════════"

# ── Check source export files exist ───────────────────────────────────────────
if [[ ! -f "$BASE_DIR/export.csv" ]] || [[ ! -f "$BASE_DIR/2PA.fa" ]]; then
  echo ""
  echo "ERROR: Missing source files in $BASE_DIR"
  echo "  Expected: export.csv and 2PA.fa"
  echo ""
  echo "Export these from the sequence database and place them in that directory,"
  echo "then re-run this script."
  exit 1
fi
echo ""
echo "Source files found in $BASE_DIR"

# ── Step 1: Prepare input.fa and metadata.tsv ─────────────────────────────────
echo ""
echo "── Step 1/3: Preparing input.fa and metadata.tsv ────────────────────────"
"$HOME/.conda/R_shared/bin/Rscript" scripts/prepare_dataset.R "$DATASET_DATE"

# ── Step 2: Ensure metadata_corrected.tsv exists ──────────────────────────────
# The batch report joins BLAST hits to metadata_corrected.tsv on sseqid == id.
# Lineage curation no longer feeds an internal designation step, so by default
# metadata_corrected.tsv is just a copy of metadata.tsv. Hand-edit it afterwards
# if you want to override any per-reference annotations.
echo ""
echo "── Step 2/3: Ensuring metadata_corrected.tsv ────────────────────────────"
META="$BASE_DIR/metadata.tsv"
META_CORR="$BASE_DIR/metadata_corrected.tsv"
if [[ ! -f "$META_CORR" ]]; then
  cp "$META" "$META_CORR"
  echo "  Created $META_CORR (copy of metadata.tsv)"
else
  echo "  Keeping existing $META_CORR (hand-curated — not overwritten)"
fi

# ── Step 3: Build the BLAST database ──────────────────────────────────────────
echo ""
echo "── Step 3/3: Building the BLAST database ────────────────────────────────"
bash scripts/make_blast_db.sh "$DATASET_DATE"

echo ""
echo "════════════════════════════════════════════════════════"
echo " Done. Local BLAST dataset ready at: $BASE_DIR"
echo ""
echo " Next: analyse a batch (BLAST + lineage assignment + report):"
echo "   bash scripts/hav_analyse.sh --mode sanger data/Batch-1 \\"
echo "     --samplesheet data/Batch-1/samplesheet.tsv --dataset-date $DATASET_DATE"
echo "════════════════════════════════════════════════════════"
