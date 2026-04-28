#!/bin/bash

###############################################################################
# Multi-Reference HAV Alignment Workflow
# 
# Runs all sequences against all genotype-specific references and selects
# the best match based on alignment quality metrics.
#
# Available genotypes: 1a, 1b, 3a
# Decision metric: alignmentScore (higher is better)
#
# Usage: ./multi_reference_alignment.sh
###############################################################################

set -e

# Configuration
INPUT_DIR="data"
GENOTYPES=("1a" "1b" "3a")
FINAL_OUTPUT="final_multi_reference_results"
TEMP_DIR="temp_multi_ref"
RESULTS_DIR="multi_ref_outputs"

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== Multi-Reference HAV Alignment Workflow ===${NC}\n"

# Create directories
mkdir -p "$TEMP_DIR" "$FINAL_OUTPUT" "$RESULTS_DIR"

# ============================================================================
# STEP 1: Run all sequences against each genotype reference
# ============================================================================
echo -e "${YELLOW}[STEP 1] Running all sequences against all references...${NC}\n"

for genotype in "${GENOTYPES[@]}"; do
  dataset_path="$INPUT_DIR"
  
  # Map genotype to dataset directory
  case "$genotype" in
    1a)
      # For 1a, check if dataset exists; if not use whole-genome as proxy
      if [ -d "$INPUT_DIR/hav_1a" ]; then
        dataset_path="$INPUT_DIR/hav_1a"
      else
        echo -e "${YELLOW}⚠ Note: 1a reference not found. Skipping genotype 1a.${NC}"
        continue
      fi
      ;;
    1b)
      dataset_path="$INPUT_DIR/vp1-2b-junction"
      ;;
    3a)
      dataset_path="$INPUT_DIR/3a_AJ299464.3"
      ;;
  esac
  
  if [ ! -d "$dataset_path" ]; then
    echo -e "${YELLOW}⚠ Skipping $genotype: dataset not found at $dataset_path${NC}"
    continue
  fi
  
  output_dir="$RESULTS_DIR/vs_${genotype}"
  echo -e "${BLUE}Processing genotype ${genotype}...${NC}"
  
  nextclade run \
    --include-reference \
    --input-dataset "$dataset_path" \
    --output-all "$output_dir/" \
    "$INPUT_DIR"/Sekvens* 2>&1 | grep -E "Processing|Finished" || true
  
  echo -e "${GREEN}✓ ${genotype} complete${NC}\n"
done

# ============================================================================
# STEP 2: Extract alignment scores and create comparison matrix
# ============================================================================
echo -e "${YELLOW}[STEP 2] Comparing alignment metrics...${NC}\n"

# Create header for comparison file
> "$TEMP_DIR/alignment_scores.txt"
echo "seqName,1a_score,1a_coverage,1a_subs,1b_score,1b_coverage,1b_subs,3a_score,3a_coverage,3a_subs,best_ref,best_score" >> "$TEMP_DIR/alignment_scores.txt"

# Get unique sequence names from any of the output files (they all have the same sequences)
first_result=$(ls "$RESULTS_DIR"/vs_*/nextclade.tsv | head -1)

if [ ! -f "$first_result" ]; then
  echo -e "${RED}Error: No alignment results found!${NC}"
  exit 1
fi

# Extract unique sequence names
unique_seqs=$(tail -n +2 "$first_result" | cut -f2 | sort -u)

# For each sequence, extract metrics from each genotype's output
while IFS= read -r seqname; do
  declare -A scores
  declare -A coverage
  declare -A subs
  
  best_ref="NONE"
  best_score=-999
  
  # Check each genotype
  for genotype in "${GENOTYPES[@]}"; do
    tsv_file="$RESULTS_DIR/vs_${genotype}/nextclade.tsv"
    
    if [ -f "$tsv_file" ]; then
      # Extract alignment metrics for this sequence
      metrics=$(grep -F "$seqname" "$tsv_file" | head -1)
      
      if [ -n "$metrics" ]; then
        # Extract columns: alignmentScore, coverage, totalSubstitutions
        # alignmentScore is column 15, coverage is column 19, totalSubstitutions is column 6
        alignment_score=$(echo "$metrics" | cut -f15)
        coverage=$(echo "$metrics" | cut -f19)
        total_subs=$(echo "$metrics" | cut -f6)
        
        scores[$genotype]=$alignment_score
        coverage[$genotype]=$coverage
        subs[$genotype]=$total_subs
        
        # Determine best match (highest alignment score) - using bash arithmetic
        if [ -z "$alignment_score" ] || [ "$alignment_score" = "" ]; then
          alignment_score=0
        fi
        # Convert to integer for comparison
        alignment_score_int=${alignment_score%.*}
        best_score_int=${best_score%.*}
        
        if (( alignment_score_int > best_score_int )); then
          best_score=$alignment_score
          best_ref=$genotype
        fi
      else
        scores[$genotype]="N/A"
        coverage[$genotype]="N/A"
        subs[$genotype]="N/A"
      fi
    else
      scores[$genotype]="SKIP"
      coverage[$genotype]="SKIP"
      subs[$genotype]="SKIP"
    fi
  done
  
  # Write comparison row
  # Handle non-numeric values in printf
  printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
    "$seqname" \
    "${scores[1a]:-N/A}" "${coverage[1a]:-N/A}" "${subs[1a]:-N/A}" \
    "${scores[1b]:-N/A}" "${coverage[1b]:-N/A}" "${subs[1b]:-N/A}" \
    "${scores[3a]:-N/A}" "${coverage[3a]:-N/A}" "${subs[3a]:-N/A}" \
    "$best_ref" \
    "$best_score" >> "$TEMP_DIR/alignment_scores.txt"
    
done < <(echo "$unique_seqs")

echo -e "${GREEN}✓ Comparison complete${NC}\n"

# ============================================================================
# STEP 3: Select best alignment for each sequence and create final output
# ============================================================================
echo -e "${YELLOW}[STEP 3] Creating final merged results...${NC}\n"

# Create final TSV (header + best match for each sequence)
> "$FINAL_OUTPUT/nextclade_best_matches.tsv"

# Get header from first available TSV
header_file=$(ls "$RESULTS_DIR"/vs_*/nextclade.tsv | head -1)
head -1 "$header_file" >> "$FINAL_OUTPUT/nextclade_best_matches.tsv"

# Add selected sequences from each genotype
for genotype in "${GENOTYPES[@]}"; do
  tsv_file="$RESULTS_DIR/vs_${genotype}/nextclade.tsv"
  
  if [ -f "$tsv_file" ]; then
    # Get sequences where this genotype was best match
    best_seqs=$(grep -f <(tail -n +2 "$TEMP_DIR/alignment_scores.txt" | awk -F',' -v g="$genotype" '$NF==g {print $1}') "$tsv_file" 2>/dev/null || true)
    
    if [ -n "$best_seqs" ]; then
      echo "$best_seqs" >> "$FINAL_OUTPUT/nextclade_best_matches.tsv"
    fi
  fi
done

# Create aligned sequences FASTA
> "$FINAL_OUTPUT/sequences_best_aligned.fasta"
for genotype in "${GENOTYPES[@]}"; do
  fasta_file="$RESULTS_DIR/vs_${genotype}/nextclade.aligned.fasta"
  
  if [ -f "$fasta_file" ]; then
    # Extract sequences that were selected as best match for this genotype
    best_seqs=$(tail -n +2 "$TEMP_DIR/alignment_scores.txt" | awk -F',' -v g="$genotype" '$NF==g {print $1}')
    
    while read -r seqname; do
      # Extract this sequence from the FASTA
      sed -n "/^>${seqname}/,/^>/p" "$fasta_file" | head -n -1 >> "$FINAL_OUTPUT/sequences_best_aligned.fasta" 2>/dev/null || true
    done < <(echo "$best_seqs")
  fi
done

echo -e "${GREEN}✓ Final results created${NC}\n"

# ============================================================================
# STEP 4: Create summary report
# ============================================================================
echo -e "${YELLOW}[STEP 4] Generating summary report...${NC}\n"

> "$FINAL_OUTPUT/genotype_selection_summary.txt"

echo "==== GENOTYPE SELECTION SUMMARY ====" >> "$FINAL_OUTPUT/genotype_selection_summary.txt"
echo "" >> "$FINAL_OUTPUT/genotype_selection_summary.txt"

# Count sequences per genotype
for genotype in "${GENOTYPES[@]}"; do
  count=$(tail -n +2 "$TEMP_DIR/alignment_scores.txt" | awk -F',' -v g="$genotype" '$NF==g' | wc -l)
  if [ $count -gt 0 ]; then
    echo "Genotype $genotype: $count sequences" >> "$FINAL_OUTPUT/genotype_selection_summary.txt"
  fi
done

echo "" >> "$FINAL_OUTPUT/genotype_selection_summary.txt"
echo "==== DETAILED SELECTION METRICS ====" >> "$FINAL_OUTPUT/genotype_selection_summary.txt"
echo "" >> "$FINAL_OUTPUT/genotype_selection_summary.txt"

# Show detailed metrics for each sequence
cat "$TEMP_DIR/alignment_scores.txt" | column -t -s',' >> "$FINAL_OUTPUT/genotype_selection_summary.txt"

# Display summary
echo ""
echo -e "${GREEN}=== SUMMARY ===${NC}"
echo ""
tail -n +2 "$TEMP_DIR/alignment_scores.txt" | awk -F',' '{print $11}' | sort | uniq -c | sort -rn | \
  while read count genotype; do
    echo "  $(printf "%-15s" "$genotype") → $count sequence(s)"
  done

echo ""
echo -e "Results directory: ${YELLOW}$RESULTS_DIR${NC}"
echo "  - vs_1a/nextclade.* (1a alignments)"
echo "  - vs_1b/nextclade.* (1b alignments)"
echo "  - vs_3a/nextclade.* (3a alignments)"
echo ""
echo -e "Final selected results: ${YELLOW}$FINAL_OUTPUT${NC}"
echo "  - nextclade_best_matches.tsv (best alignment per sequence)"
echo "  - sequences_best_aligned.fasta (best aligned sequences)"
echo "  - genotype_selection_summary.txt (detailed metrics and selection logic)"
echo ""
echo -e "${GREEN}✓ Workflow complete!${NC}"
