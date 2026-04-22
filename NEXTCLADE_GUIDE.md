# NextClade Quick Start Guide

## What is NextClade?

NextClade is a tool for viral genome alignment, mutation calling, clade assignment, quality checks and phylogenetic placement.

## Basic Commands

### 1. Download a Dataset

```bash
# List available datasets
nextclade dataset list

# Download a specific dataset (e.g., SARS-CoV-2)
nextclade dataset get --name sars-cov-2 --output-dir data/sars-cov-2
```

### 2. Run Analysis

```bash
# Basic analysis
nextclade run \
  --input-dataset data/sars-cov-2 \
  --output-all output/ \
  sequences.fasta

# This will generate:
# - nextclade.tsv (main results)
# - nextclade.aligned.fasta (aligned sequences)
# - nextclade.json (detailed results)
# - nextclade.tree.json (phylogenetic tree)
```

### 3. Common Options

```bash
# Specify output format
nextclade run --output-tsv results.tsv input.fasta

# Run with quality checks
nextclade run --include-reference --input-dataset data input.fasta

# Parallel processing
nextclade run --jobs 4 input.fasta
```

## Typical Workflow

```bash
# 1. Create a data directory
mkdir -p data output

# 2. Download your dataset
nextclade dataset get --name <your-virus> --output-dir data/<your-virus>

# 3. Place your FASTA files in the project directory

# 4. Run analysis
nextclade run \
  --input-dataset data/<your-virus> \
  --output-all output/ \
  your_sequences.fasta

# 5. Review results in output/ directory
```

## Output Files Explained

- **nextclade.tsv** - Tab-separated values with mutations, clades, QC metrics
- **nextclade.aligned.fasta** - Your sequences aligned to reference
- **nextclade.json** - Detailed JSON output for programmatic access
- **nextclade.tree.json** - Phylogenetic tree placement

## Tips

- Use `--help` on any command for full options
- Results are reproducible - same input always gives same output
- Large datasets: use `--jobs` for parallel processing
- Quality control: check the QC columns in TSV output

## Resources

- Documentation: https://docs.nextstrain.org/projects/nextclade/
- GitHub: https://github.com/nextstrain/nextclade
- Datasets: https://github.com/nextstrain/nextclade_data

## Example: Hepatitis A Analysis

```bash
# Download HAV dataset (if available)
nextclade dataset list | grep -i hepatitis

# Or use custom reference genome
nextclade run \
  --input-ref reference.fasta \
  --output-all output/ \
  samples.fasta
```
