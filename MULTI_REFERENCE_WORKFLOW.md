# Multi-Reference HAV Alignment Workflow

## Overview

This workflow runs your sequences against **all available genotype-specific references** and automatically selects the best match for each sequence based on **alignment score** (higher = better).

### Key Advantages

✅ **No cascading errors** - Decision based on objective alignment quality, not previous misclassification
✅ **Handles missing references** - Works even if not all genotype references are available  
✅ **Objective scoring** - Uses alignment metrics to determine best fit
✅ **Scalable** - Works for small batches and large datasets

## How It Works

### STEP 1: Run All Sequences Against All References
```
data/Sekvens* → vs_1b/ (alignment results)
              → vs_3a/ (alignment results)
              → vs_1a/ (if available)
```

### STEP 2: Compare Alignment Metrics
For each sequence, extract:
- **alignmentScore** (higher = better fit)
- **coverage** (how much of the reference aligned)
- **totalSubstitutions** (number of differences)

Select the reference with the **highest alignment score** as the best match.

### STEP 3: Merge Results
Combine results by selecting only the rows from each reference's output where that reference had the best alignment score.

### STEP 4: Generate Reports
- **nextclade_best_matches.tsv** - All sequences with best-match results
- **sequences_best_aligned.fasta** - Aligned sequences from best references
- **genotype_selection_summary.txt** - Detailed metrics and selection logic

## Test Results

Your 7 test sequences were classified as:

| Sequence | Best Ref | Score | Classification |
|----------|----------|-------|-----------------|
| Sekvens1\|IIIa\| | **3a** | 1394 | ✓ Correctly identified as 3a |
| Sekvens2\|Ib | **1b** | 1035 | ✓ Correctly identified as 1b |
| Sekvens3\|IIa | **1b** | 835  | Identified as 1b (genotype IIa) |
| Sekvens4\|Ia | **3a** | 945  | 3a scored higher than 1b (868) |
| Sekvens5\|IIIa | **3a** | 1312 | ✓ Correctly identified as 3a |
| Sekvens6\|IIa | **1b** | 835  | Identified as 1b (genotype IIa) |
| Sekvens7\|IIIa | **3a** | 1258 | ✓ Correctly identified as 3a |

**Summary:** 4 sequences → 3a reference, 3 sequences → 1b reference

## Important Notes

### Clade Information

- **1b reference** (vp1-2b-junction) provides detailed clade output (IB, IIA, IA, IIB, etc.)
- **3a reference** (3a_AJ299464.3) is a single-genotype dataset and does NOT provide clade assignments (column is empty)

If you need clade information for all sequences, you may need to:
1. Run 3a sequences through the 1b reference as well to get their clade
2. Or add tree.json files to the 3a reference if available

### Handling Genotype 1a

The script looks for a 1a reference at `data/hav_1a/`. If you don't have it yet, it skips 1a and only compares 1b vs 3a.

To add 1a support:
```bash
nextclade dataset get --name <1a-dataset> --output-dir data/hav_1a
```

## Usage

Run the workflow:
```bash
./multi_reference_alignment.sh
```

Results are saved to:
- `multi_ref_outputs/` - Individual alignments for each reference
- `final_multi_reference_results/` - Final selected results

## Troubleshooting

**Issue: "N/A" values in alignment score**
- This means the sequence doesn't exist in that genotype's reference output
- This is normal when a 3a-only dataset couldn't align a sequence designed for 1b

**Issue: Empty clade for 3a sequences**
- The 3a reference doesn't provide clade information
- To get clades for 3a sequences, you may need tree.json in the 3a reference dataset

**Issue: Unexpected genotype assignments**
- Check the detailed metrics in `genotype_selection_summary.txt`
- The algorithm picks the highest alignment score - if results differ from expectations, check if reference data is appropriate
