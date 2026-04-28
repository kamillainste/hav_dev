#!/bin/bash

###############################################################################
# Genotype-Aware HAV Alignment Workflow
# 
# This script performs two-pass alignment:
# 1. Pass 1: Quick genotype identification using whole-genome reference
# 2. Pass 2: Re-align each genotype with its specific reference
#
# Usage: ./genotype_aware_alignment.sh
###############################################################################

set -e

# Configuration
INPUT_DIR="data"
OUTPUT_PASS1="output_pass1_genotyping"
OUTPUT_FINAL_3A="output_final_3a"
OUTPUT_FINAL_IB="output_final_1b"
OUTPUT_FINAL_OTHER="output_final_other"
TEMP_DIR="temp_genotyping"
FINAL_OUTPUT="final_aligned_results"

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== HAV Genotype-Aware Alignment Workflow ===${NC}\n"

# Create directories
mkdir -p "$TEMP_DIR" "$FINAL_OUTPUT"

# ============================================================================
# PASS 1: Genotype identification using whole-genome reference
# ============================================================================
echo -e "${YELLOW}[PASS 1] Running genotype identification...${NC}"
echo "Input: $INPUT_DIR/Sekvens*"

nextclade run \
  --include-reference \
  --input-dataset "$INPUT_DIR/whole-genome/" \
  --output-all "$OUTPUT_PASS1/" \
  "$INPUT_DIR"/Sekvens* 2>&1 | grep -E "Processing|Finished|Writing" || true

echo -e "${GREEN}✓ Genotype identification complete${NC}\n"

# ============================================================================
# Parse genotypes and create sequence lists
# ============================================================================
echo -e "${YELLOW}[PASS 1] Parsing genotypes from results...${NC}"

# Create empty files
> "$TEMP_DIR/seqs_3a.txt"
> "$TEMP_DIR/seqs_1b.txt"
> "$TEMP_DIR/seqs_other.txt"
> "$TEMP_DIR/genotype_map.txt"  # For logging

# Parse TSV and extract genotypes
# Skip header (line 1) and read each line
while IFS=$'\t' read -r index seqName clade rest_of_fields; do
  # Extract just the base sequence name (remove pipe and after)
  base_name="${seqName%%\|*}"
  
  # Try to extract genotype from sequence header label (e.g., |IIIa|, |Ib, |IIa, |Ia)
  header_genotype=""
  if [[ "$seqName" =~ \|IIIa ]]; then
    header_genotype="3a"
  elif [[ "$seqName" =~ \|Ib ]]; then
    header_genotype="1b"
  elif [[ "$seqName" =~ \|IIa ]]; then
    header_genotype="2a"
  elif [[ "$seqName" =~ \|Ia ]]; then
    header_genotype="1a"
  fi
  
  # Determine final classification (prefer header label if present, otherwise use clade)
  final_genotype="$clade"
  source="clade"
  
  if [[ "$header_genotype" != "" && "$header_genotype" != "$clade" ]]; then
    # Header disagrees with clade - log this discrepancy
    echo "$seqName: header=$header_genotype, pass1=$clade" >> "$TEMP_DIR/genotype_map.txt"
    final_genotype="$header_genotype"
    source="header (override)"
  fi
  
  # Classify by final genotype
  if [[ "$final_genotype" == "3a" ]]; then
    echo "$base_name" >> "$TEMP_DIR/seqs_3a.txt"
    echo "  $seqName → 3a (source: $source)"
  elif [[ "$final_genotype" == "1b" || "$final_genotype" == "IB" ]]; then
    echo "$base_name" >> "$TEMP_DIR/seqs_1b.txt"
    echo "  $seqName → 1b (source: $source)"
  else
    echo "$base_name" >> "$TEMP_DIR/seqs_other.txt"
    echo "  $seqName → other: $clade (source: $source)"
  fi
done < <(tail -n +2 "$OUTPUT_PASS1/nextclade.tsv")

# Show discrepancies if any
if [ -s "$TEMP_DIR/genotype_map.txt" ]; then
  echo -e "${YELLOW}⚠ Discrepancies detected (header vs Pass 1):${NC}"
  cat "$TEMP_DIR/genotype_map.txt"
  echo ""
fi

echo ""

# ============================================================================
# PASS 2: Re-align each genotype with its specific reference
# ============================================================================
echo -e "${YELLOW}[PASS 2] Re-aligning with genotype-specific references...${NC}\n"

# Process 3a samples
if [ -s "$TEMP_DIR/seqs_3a.txt" ]; then
  echo -e "${BLUE}Processing genotype 3a samples:${NC}"
  seqs_3a=$(cat "$TEMP_DIR/seqs_3a.txt" | xargs -I {} find "$INPUT_DIR"/Sekvens* -maxdepth 0 -name "{}*" 2>/dev/null || echo "$INPUT_DIR/Sekvens${seqs_3a_list}")
  
  # Better approach: use the sequence names to find files
  seqs_3a_files=""
  while read -r seqname; do
    for file in "$INPUT_DIR"/Sekvens*; do
      if [[ -f "$file" && "$file" == *"$seqname"* ]]; then
        seqs_3a_files="$seqs_3a_files $file"
      fi
    done
  done < "$TEMP_DIR/seqs_3a.txt"
  
  # Also check in data directory
  for seqname in $(cat "$TEMP_DIR/seqs_3a.txt"); do
    if [[ -f "$INPUT_DIR/$seqname.fa" ]]; then
      seqs_3a_files="$seqs_3a_files $INPUT_DIR/$seqname.fa"
    fi
  done
  
  if [[ -n "$seqs_3a_files" ]]; then
    nextclade run \
      --include-reference \
      --input-dataset "$INPUT_DIR/3a_AJ299464.3/" \
      --output-all "$OUTPUT_FINAL_3A/" \
      $seqs_3a_files 2>&1 | grep -E "Processing|Finished|Writing" || true
    echo -e "${GREEN}✓ 3a re-alignment complete${NC}\n"
  else
    echo -e "${YELLOW}! No 3a sequences found${NC}\n"
  fi
else
  echo -e "${YELLOW}! No genotype 3a sequences detected${NC}\n"
fi

# Process 1b samples
if [ -s "$TEMP_DIR/seqs_1b.txt" ]; then
  echo -e "${BLUE}Processing genotype 1b samples:${NC}"
  
  seqs_1b_files=""
  while read -r seqname; do
    if [[ -f "$INPUT_DIR/$seqname.fa" ]]; then
      seqs_1b_files="$seqs_1b_files $INPUT_DIR/$seqname.fa"
    fi
  done < "$TEMP_DIR/seqs_1b.txt"
  
  if [[ -n "$seqs_1b_files" ]]; then
    nextclade run \
      --include-reference \
      --input-dataset "$INPUT_DIR/vp1-2b-junction/" \
      --output-all "$OUTPUT_FINAL_IB/" \
      $seqs_1b_files 2>&1 | grep -E "Processing|Finished|Writing" || true
    echo -e "${GREEN}✓ 1b re-alignment complete${NC}\n"
  else
    echo -e "${YELLOW}! No 1b sequences found${NC}\n"
  fi
else
  echo -e "${YELLOW}! No genotype 1b sequences detected${NC}\n"
fi

# Process other genotypes (if any)
if [ -s "$TEMP_DIR/seqs_other.txt" ]; then
  echo -e "${BLUE}Processing other genotypes:${NC}"
  
  seqs_other_files=""
  while read -r seqname; do
    if [[ -f "$INPUT_DIR/$seqname.fa" ]]; then
      seqs_other_files="$seqs_other_files $INPUT_DIR/$seqname.fa"
    fi
  done < "$TEMP_DIR/seqs_other.txt"
  
  if [[ -n "$seqs_other_files" ]]; then
    nextclade run \
      --include-reference \
      --input-dataset "$INPUT_DIR/whole-genome/" \
      --output-all "$OUTPUT_FINAL_OTHER/" \
      $seqs_other_files 2>&1 | grep -E "Processing|Finished|Writing" || true
    echo -e "${GREEN}✓ Other genotypes re-alignment complete${NC}\n"
  fi
else
  echo -e "${YELLOW}! No other genotypes detected${NC}\n"
fi

# ============================================================================
# Merge results
# ============================================================================
echo -e "${YELLOW}[MERGE] Combining all results...${NC}"

# Merge TSV files (keep header from first file)
header_printed=false
> "$FINAL_OUTPUT/nextclade_combined.tsv"

for tsv_file in "$OUTPUT_FINAL_3A"/nextclade.tsv "$OUTPUT_FINAL_IB"/nextclade.tsv "$OUTPUT_FINAL_OTHER"/nextclade.tsv; do
  if [[ -f "$tsv_file" ]]; then
    if [[ "$header_printed" == false ]]; then
      cat "$tsv_file" >> "$FINAL_OUTPUT/nextclade_combined.tsv"
      header_printed=true
    else
      tail -n +2 "$tsv_file" >> "$FINAL_OUTPUT/nextclade_combined.tsv"
    fi
  fi
done

# Merge FASTA files
> "$FINAL_OUTPUT/sequences_aligned.fasta"
for fasta_file in "$OUTPUT_FINAL_3A"/nextclade.aligned.fasta "$OUTPUT_FINAL_IB"/nextclade.aligned.fasta "$OUTPUT_FINAL_OTHER"/nextclade.aligned.fasta; do
  if [[ -f "$fasta_file" ]]; then
    cat "$fasta_file" >> "$FINAL_OUTPUT/sequences_aligned.fasta"
  fi
done

echo -e "${GREEN}✓ Results merged${NC}\n"

# ============================================================================
# Summary
# ============================================================================
echo -e "${BLUE}=== SUMMARY ===${NC}"
echo -e "Pass 1 results (genotyping): ${YELLOW}$OUTPUT_PASS1${NC}"
echo ""
echo -e "Pass 2 results:"
[[ -d "$OUTPUT_FINAL_3A" ]] && echo "  - Genotype 3a: ${YELLOW}$OUTPUT_FINAL_3A${NC}"
[[ -d "$OUTPUT_FINAL_IB" ]] && echo "  - Genotype 1b: ${YELLOW}$OUTPUT_FINAL_IB${NC}"
[[ -d "$OUTPUT_FINAL_OTHER" ]] && echo "  - Other: ${YELLOW}$OUTPUT_FINAL_OTHER${NC}"
echo ""
echo -e "Combined results: ${YELLOW}$FINAL_OUTPUT${NC}"
echo "  - nextclade_combined.tsv (all results merged)"
echo "  - sequences_aligned.fasta (all aligned sequences)"
echo ""
echo -e "${GREEN}✓ Workflow complete!${NC}"
