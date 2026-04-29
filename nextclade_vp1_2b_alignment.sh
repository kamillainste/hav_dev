#!/bin/bash

# Script: nextclade_vp1_2b_alignment.sh
# Description: Align multiple .fa sequences to vp1-2b-junction dataset using nextclade
#              with high-diversity alignment preset, including example sequences for context.
#              Outputs both full results and filtered results (new sequences only) to result.tsv.
# Usage: ./nextclade_vp1_2b_alignment.sh /path/to/sequences/folder

set -e

# ============================================================================
# Input validation
# ============================================================================

if [ $# -eq 0 ]; then
    echo "Error: Missing input folder argument"
    echo "Usage: $0 /path/to/sequences/folder"
    exit 1
fi

INPUT_FOLDER="$1"

# Check if input folder exists
if [ ! -d "$INPUT_FOLDER" ]; then
    echo "Error: Input folder does not exist: $INPUT_FOLDER"
    exit 1
fi

# Check if there are any .fa files
if [ -z "$(find "$INPUT_FOLDER" -maxdepth 1 -name '*.fa' -o -name '*.fasta')" ]; then
    echo "Error: No .fa or .fasta files found in $INPUT_FOLDER"
    exit 1
fi

# ============================================================================
# Setup paths
# ============================================================================

# Get the parent directory of the input folder
PARENT_DIR=$(dirname "$INPUT_FOLDER")

# Create output folder structure
OUTPUT_FOLDER="${PARENT_DIR}/Results"
mkdir -p "$OUTPUT_FOLDER"

# Get the absolute paths
INPUT_FOLDER=$(cd "$INPUT_FOLDER" && pwd)
OUTPUT_FOLDER=$(cd "$OUTPUT_FOLDER" && pwd)

echo "============================================================================"
echo "Nextclade VP1-2B Junction Alignment"
echo "============================================================================"
echo "Input folder:    $INPUT_FOLDER"
echo "Output folder:   $OUTPUT_FOLDER"
echo "============================================================================"
echo ""

# ============================================================================
# Prepare combined input (new sequences + example sequences)
# ============================================================================

TEMP_COMBINED="${OUTPUT_FOLDER}/.temp_combined.fa"
EXAMPLE_SEQS="data/vp1-2b-junction/example_sequences.fasta"

# Combine new sequences with example sequences
# Note: Using printf to add newlines between files (handles files without trailing newlines)
echo "Preparing combined sequence file with example sequences..."
(for f in "$INPUT_FOLDER"/*.fa; do cat "$f"; printf '\n'; done; cat "$EXAMPLE_SEQS") > "$TEMP_COMBINED"

# Extract list of example sequence names for filtering later
EXAMPLE_NAMES_FILE="${OUTPUT_FOLDER}/.example_names.txt"
grep "^>" "$EXAMPLE_SEQS" | sed 's/^>//' > "$EXAMPLE_NAMES_FILE"

# ============================================================================
# Run nextclade with combined sequences
# ============================================================================

echo "Running nextclade alignment with example sequences for context..."
nextclade run "$TEMP_COMBINED" \
    --input-dataset data/vp1-2b-junction/ \
    --alignment-preset high-diversity \
    --output-all "$OUTPUT_FOLDER/"

# ============================================================================
# Add mine_samples column to nextclade.tsv for Auspice coloring
# ============================================================================

echo "Adding mine_samples column to nextclade.tsv..."

# Temporary file for the modified TSV
TEMP_TSV="${OUTPUT_FOLDER}/.temp_nextclade.tsv"

# Add header with new column
head -1 "$OUTPUT_FOLDER/nextclade.tsv" | sed 's/$/\tmine_samples/' > "$TEMP_TSV"

# Add data rows with mine_samples values
tail -n +2 "$OUTPUT_FOLDER/nextclade.tsv" | while IFS=$'\t' read -r index seqname rest; do
    # Check if this sequence name is in the example names file
    if grep -Fxq "$seqname" "$EXAMPLE_NAMES_FILE"; then
        echo "${index}	${seqname}	${rest}	no"
    else
        echo "${index}	${seqname}	${rest}	yes"
    fi
done >> "$TEMP_TSV"

# Replace original nextclade.tsv with the modified version
mv "$TEMP_TSV" "$OUTPUT_FOLDER/nextclade.tsv"

# ============================================================================
# Extract only new sequences to result.tsv
# ============================================================================

echo "Extracting results for new sequences to result.tsv..."

# Get header from nextclade.tsv (now includes mine_samples column)
HEADER=$(head -1 "$OUTPUT_FOLDER/nextclade.tsv")

# Create result.tsv with header
echo "$HEADER" > "$OUTPUT_FOLDER/result.tsv"

# Add only rows where seqName (column 2) is NOT in the example sequences list
tail -n +2 "$OUTPUT_FOLDER/nextclade.tsv" | while IFS=$'\t' read -r index seqname rest mine_samples; do
    # Check if this sequence name is in the example names file
    if ! grep -Fxq "$seqname" "$EXAMPLE_NAMES_FILE"; then
        echo "${index}	${seqname}	${rest}	${mine_samples}"
    fi
done >> "$OUTPUT_FOLDER/result.tsv"

# ============================================================================
# Cleanup temporary files
# ============================================================================

rm "$TEMP_COMBINED" "$EXAMPLE_NAMES_FILE"

echo ""
echo "============================================================================"
echo "Alignment complete!"
echo "Full results:      $OUTPUT_FOLDER/nextclade.tsv"
echo "New sequences:     $OUTPUT_FOLDER/result.tsv"
echo "============================================================================"
