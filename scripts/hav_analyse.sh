#!/usr/bin/env bash
# scripts/hav_analyse.sh
#
# Master wrapper for the HAV batch analysis pipeline.
#
# Usage (from project root, HAVDEV conda active):
#   bash scripts/hav_analyse.sh --mode <sanger|wgs> <batch_dir> [OPTIONS]
#
# Modes:
#   sanger   Short-fragment Sanger sequences (VP1/2A or VP1/2B junction)
#   wgs      Whole-genome Nanopore tiling-PCR sequences (requires ViroConstrictor)
#
# Positional argument:
#   <batch_dir>   Output directory for this batch, relative to project root
#                 Sanger example: data/Batch-1
#                 WGS example:    data/HAV-2026_01
#
# Required for both modes:
#   --samplesheet <path>    TSV with columns Sample and InputDir.
#                           Sanger: InputDir = dir containing .fa/.fasta/.TXT
#                           WGS:    InputDir = dir containing raw FASTQ files
#                                   (also needs ViroConstrictor columns; see
#                                    data/HAV-2026_01/samplesheet.tsv as example)
#
# Common options:
#   --dataset-date <date>   Local dataset version to use (default: 2026-04-10)
#   --threads <n>           CPU threads (default: 4)
#   --n-neighbors <n>       Phylo-tree neighbors per query sequence (default: 30)
#   --skip-trees            Skip build_per_seq_trees.sh
#   --skip-report           Skip batch_report.Rmd rendering
#   -h, --help
#
# Sanger-specific options:
#   --primer-names <csv>    Comma-separated primer names — enables primer trimming
#   --primers-file <path>   Custom primer FASTA/CSV/TSV (optional supplement)
#
# WGS-specific options:
#   --skip-assembly         Skip ViroConstrictor; batch FASTA must already exist
#                           (--samplesheet still required unless FASTA exists)
#
# Workflow steps by mode:
#   sanger:
#     1. prepare_input_fasta.R   — samplesheet → <batch>.fa
#     2. run_all_analyses.sh     — BLAST + NextClade (with optional primer trim)
#     3. build_per_seq_trees.sh  — per-sequence ML trees
#     4. batch_report.Rmd        — HTML report
#
#   wgs:
#     1. run_pipeline.sh         — NanoPlot QC + Chopper + ViroConstrictor assembly
#     2. collect consensus        — ViroConstrictor all_consensus.fasta → <batch>.fa
#     3. run_all_analyses.sh     — BLAST + NextClade
#     4. build_per_seq_trees.sh  — per-sequence ML trees
#     5. batch_report.Rmd        — HTML report
#
# Prerequisites:
#   - HAVDEV conda environment active (nextclade, blastn, mafft, iqtree3)
#   - $HOME/.conda/R_shared  (for R analyses; built by setup_havdev_env.sh)
#   - Local dataset built: bash scripts/make_blast_db.sh <dataset_date>
#   - Community dataset present: data/nextclade_datasets/hav-vp1-2b-lineages/

set -euo pipefail

# ── Constants ─────────────────────────────────────────────────────────────────
RSCRIPT="$HOME/.conda/R_shared/bin/Rscript"
DEFAULT_DATASET_DATE="2026-04-10"
DEFAULT_THREADS=4
DEFAULT_N_NEIGHBORS=30

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  cat <<'USAGE'
Usage:
  bash scripts/hav_analyse.sh --mode <sanger|wgs> <batch_dir> [OPTIONS]

Modes:
  sanger   Short-fragment Sanger sequences (VP1/2A or VP1/2B junction)
  wgs      Whole-genome Nanopore tiling-PCR sequences

Positional:
  <batch_dir>              Batch directory relative to project root

Common options:
  --dataset-date <date>    Dataset version (default: 2026-04-10)
  --threads <n>            CPU threads (default: 4)
  --n-neighbors <n>        Tree neighbors per query (default: 30)
  --skip-trees             Skip phylogenetic tree building
  --skip-report            Skip HTML report generation
  -h, --help

Required (both modes):
  --samplesheet <path>     TSV with Sample + InputDir columns
                           Sanger: InputDir = dir with .fa/.fasta/.TXT file
                           WGS:    InputDir = dir with raw FASTQ files

Sanger options:
  --primer-names <csv>     Enable primer trimming; comma-separated primer names
  --primers-file <path>    Custom primer file (FASTA/CSV/TSV)

WGS options:
  --skip-assembly          Skip ViroConstrictor; batch FASTA must already exist

Examples:
  # Sanger batch without primer trimming
  bash scripts/hav_analyse.sh --mode sanger data/Batch-1 \
    --samplesheet data/Batch-1/samplesheet.tsv

  # Sanger batch with primer trimming
  bash scripts/hav_analyse.sh --mode sanger data/Batch-1 \
    --samplesheet data/Batch-1/samplesheet.tsv \
    --primer-names "HAV_6.1_codehop,HAV_10_codehop,HAV_8.2_codehop,HAV_11_codehop"

  # WGS batch — full pipeline from raw Nanopore reads
  bash scripts/hav_analyse.sh --mode wgs data/HAV-2026_01 \
    --samplesheet data/HAV-2026_01/samplesheet.tsv

  # WGS batch — assembly already done, run analyses only
  bash scripts/hav_analyse.sh --mode wgs data/HAV-2026_01 \
    --samplesheet data/HAV-2026_01/samplesheet.tsv --skip-assembly
USAGE
}

# ── Argument parsing ──────────────────────────────────────────────────────────
MODE=""
BATCH_DIR=""
DATASET_DATE="$DEFAULT_DATASET_DATE"
THREADS="$DEFAULT_THREADS"
N_NEIGHBORS="$DEFAULT_N_NEIGHBORS"
SKIP_TREES=0
SKIP_REPORT=0
PRIMER_NAMES=""
PRIMERS_FILE=""
SAMPLESHEET=""
SKIP_ASSEMBLY=0

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      MODE="${2:?--mode requires a value (sanger or wgs)}"
      shift 2 ;;
    --dataset-date)
      DATASET_DATE="${2:?--dataset-date requires a value}"
      shift 2 ;;
    --threads)
      THREADS="${2:?--threads requires a value}"
      shift 2 ;;
    --n-neighbors)
      N_NEIGHBORS="${2:?--n-neighbors requires a value}"
      shift 2 ;;
    --skip-trees)
      SKIP_TREES=1; shift ;;
    --skip-report)
      SKIP_REPORT=1; shift ;;
    --primer-names)
      PRIMER_NAMES="${2:?--primer-names requires a value}"
      shift 2 ;;
    --primers-file)
      PRIMERS_FILE="${2:?--primers-file requires a value}"
      shift 2 ;;
    --samplesheet)
      SAMPLESHEET="${2:?--samplesheet requires a value}"
      shift 2 ;;
    --skip-assembly)
      SKIP_ASSEMBLY=1; shift ;;
    -h|--help)
      usage; exit 0 ;;
    -*)
      echo "ERROR: Unknown option: $1" >&2; usage; exit 1 ;;
    *)
      POSITIONAL+=("$1"); shift ;;
  esac
done

[[ ${#POSITIONAL[@]} -ge 1 ]] && BATCH_DIR="${POSITIONAL[0]}"

# ── Resolve paths ─────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Normalise BATCH_DIR: strip absolute prefix and leading ./
BATCH_DIR="${BATCH_DIR#"$PROJECT_DIR/"}"
BATCH_DIR="${BATCH_DIR#./}"
BATCH_DIR="${BATCH_DIR%/}"

BATCH_NAME=$(basename "$BATCH_DIR")
BATCH_ABS="$PROJECT_DIR/$BATCH_DIR"
BATCH_FA="$BATCH_ABS/${BATCH_NAME}.fa"
OUT_BASE="$BATCH_ABS/output"

cd "$PROJECT_DIR"

# ── Validation ────────────────────────────────────────────────────────────────
if [[ -z "$MODE" ]]; then
  echo "ERROR: --mode <sanger|wgs> is required" >&2; usage; exit 1
fi
if [[ "$MODE" != "sanger" && "$MODE" != "wgs" ]]; then
  echo "ERROR: --mode must be 'sanger' or 'wgs', got: '$MODE'" >&2; exit 1
fi
if [[ -z "$BATCH_DIR" ]]; then
  echo "ERROR: <batch_dir> is required" >&2; usage; exit 1
fi
# Samplesheet is required for both modes.
# The only exception: WGS + --skip-assembly when the batch FASTA already exists.
SAMPLESHEET_OPTIONAL=0
if [[ "$MODE" == "wgs" && "$SKIP_ASSEMBLY" -eq 1 && -f "$BATCH_FA" ]]; then
  SAMPLESHEET_OPTIONAL=1
fi
if [[ -z "$SAMPLESHEET" && "$SAMPLESHEET_OPTIONAL" -eq 0 ]]; then
  echo "ERROR: --samplesheet <path> is required" >&2
  if [[ "$MODE" == "wgs" ]]; then
    echo "       Use --skip-assembly to reuse an existing batch FASTA without a sheet." >&2
  fi
  exit 1
fi
if [[ -n "$SAMPLESHEET" && ! -f "$SAMPLESHEET" ]]; then
  echo "ERROR: Samplesheet not found: $SAMPLESHEET" >&2; exit 1
fi

# HAVDEV tools must be on PATH
for tool in nextclade blastn mafft; do
  if ! command -v "$tool" &>/dev/null; then
    echo "ERROR: '$tool' not found — activate the HAVDEV conda environment first." >&2
    echo "  conda activate $PROJECT_DIR/.conda/HAVDEV" >&2
    exit 1
  fi
done

if [[ ! -x "$RSCRIPT" ]]; then
  echo "ERROR: Rscript not found at $RSCRIPT" >&2; exit 1
fi

DATASET_DIR="$PROJECT_DIR/data/local_datasets/$DATASET_DATE"
if [[ ! -d "$DATASET_DIR" ]]; then
  echo "ERROR: Local dataset not found: $DATASET_DIR" >&2
  echo "  Build with: bash scripts/make_blast_db.sh $DATASET_DATE" >&2
  exit 1
fi

mkdir -p "$OUT_BASE"

# ── Logging ───────────────────────────────────────────────────────────────────
LOGFILE="$OUT_BASE/hav_analyse.log"
exec > >(tee -a "$LOGFILE") 2>&1

step() {
  echo ""
  echo "▶ $*"
  echo "─────────────────────────────────────────────────────────────────"
}

echo "════════════════════════════════════════════════════════════════"
echo "HAV Analysis Pipeline"
echo "  Batch     : $BATCH_NAME"
echo "  Mode      : $MODE"
echo "  Dataset   : $DATASET_DATE"
echo "  Threads   : $THREADS"
echo "  Neighbors : $N_NEIGHBORS"
echo "  Started   : $(date)"
echo "  Log       : $LOGFILE"
echo "════════════════════════════════════════════════════════════════"

# ════════════════════════════════════════════════════════════════════════════
# WGS BRANCH
# ════════════════════════════════════════════════════════════════════════════
if [[ "$MODE" == "wgs" ]]; then

  if [[ "$SKIP_ASSEMBLY" -eq 0 ]]; then
    step "WGS 1/2: Nanopore QC + filtering + ViroConstrictor assembly"
    bash scripts/run_pipeline.sh "$SAMPLESHEET" "$BATCH_ABS" "$THREADS"
    echo "  Samplesheet : $SAMPLESHEET"
  else
    echo ""
    echo "  Skipping assembly (--skip-assembly)"
  fi

  step "WGS 2/2: Collecting ViroConstrictor consensus sequences → batch FASTA"
  VC_CONSENSUS="$BATCH_ABS/viroconstrictor/results/combined/all_samples/all_consensus.fasta"
  if [[ ! -f "$VC_CONSENSUS" ]]; then
    echo "ERROR: ViroConstrictor consensus not found:" >&2
    echo "  $VC_CONSENSUS" >&2
    echo "  Check that assembly completed, or run without --skip-assembly." >&2
    exit 1
  fi
  # ViroConstrictor writes headers like ">barcode84 HAV LC128713.1_IB mincov=20".
  # Keep only the first token (the sample name) so the sample is referred to
  # identically by every downstream tool: BLAST truncates qseqid at the first
  # space anyway, while Nextclade and the MSA/heatmap keep the full header —
  # stripping here makes all of them agree on "barcode84".
  awk '/^>/{print $1; next} {print}' "$VC_CONSENSUS" > "$BATCH_FA"
  N_SEQS=$(grep -c "^>" "$BATCH_FA")
  echo "  Collected $N_SEQS consensus sequences → $BATCH_FA (headers trimmed to sample name)"

fi

# ════════════════════════════════════════════════════════════════════════════
# SANGER BRANCH
# ════════════════════════════════════════════════════════════════════════════
if [[ "$MODE" == "sanger" ]]; then

  step "Sanger 1/1: Prepare batch FASTA from samplesheet"
  "$RSCRIPT" scripts/prepare_input_fasta.R "$SAMPLESHEET" "$BATCH_FA"
  if [[ ! -f "$BATCH_FA" ]]; then
    echo "ERROR: prepare_input_fasta.R completed but $BATCH_FA was not created." >&2
    exit 1
  fi

fi

# ════════════════════════════════════════════════════════════════════════════
# COMMON STEPS — BLAST + NextClade → trees → report
# ════════════════════════════════════════════════════════════════════════════
if [[ ! -f "$BATCH_FA" ]]; then
  echo "ERROR: Batch FASTA not found: $BATCH_FA" >&2
  exit 1
fi

# ── BLAST + NextClade ─────────────────────────────────────────────────────────
step "BLAST + NextClade"
ANALYSES_ARGS=("$BATCH_DIR" "$DATASET_DATE")

if [[ "$MODE" == "sanger" && -n "$PRIMER_NAMES" ]]; then
  ANALYSES_ARGS+=(--sanger --primer-names "$PRIMER_NAMES")
  [[ -n "$PRIMERS_FILE" ]] && ANALYSES_ARGS+=(--primers-file "$PRIMERS_FILE")
fi

bash scripts/run_all_analyses.sh "${ANALYSES_ARGS[@]}"

# ── Per-sequence phylogenetic trees ───────────────────────────────────────────
if [[ "$SKIP_TREES" -eq 0 ]]; then
  step "Per-sequence phylogenetic trees"
  bash scripts/build_per_seq_trees.sh "$BATCH_DIR" "$DATASET_DATE" "$N_NEIGHBORS"
else
  echo ""
  echo "  Skipping tree building (--skip-trees)"
fi

# ── HTML report ───────────────────────────────────────────────────────────────
if [[ "$SKIP_REPORT" -eq 0 ]]; then
  step "Generating HTML report"
  REPORT_OUT="$OUT_BASE/batch_report.html"
  # rmarkdown needs pandoc. It is bundled in the R_shared bin dir alongside
  # Rscript, but is not otherwise on PATH, so point rmarkdown at it explicitly.
  R_BIN_DIR="$(dirname "$RSCRIPT")"
  export RSTUDIO_PANDOC="$R_BIN_DIR"
  export PATH="$R_BIN_DIR:$PATH"
  "$RSCRIPT" -e "
    rmarkdown::render(
      '$PROJECT_DIR/scripts/batch_report.Rmd',
      params = list(
        batch_dir    = '$BATCH_DIR',
        dataset_date = '$DATASET_DATE',
        n_neighbors  = $N_NEIGHBORS
      ),
      output_file = '$REPORT_OUT'
    )
  "
  echo "  Report → $REPORT_OUT"
else
  echo ""
  echo "  Skipping report generation (--skip-report)"
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════════"
echo "Pipeline complete: $BATCH_NAME"
echo "  Finished    : $(date)"
echo ""
echo "  Batch FASTA       : $BATCH_FA"
echo "  BLAST results     : $OUT_BASE/blast_results.tsv"
echo "  NextClade results : $OUT_BASE/lineages/nextclade.tsv"
[[ "$SKIP_TREES"  -eq 0 ]] && echo "  Trees             : $OUT_BASE/trees/"
[[ "$SKIP_REPORT" -eq 0 ]] && echo "  Report            : $OUT_BASE/batch_report.html"
echo "════════════════════════════════════════════════════════════════"
