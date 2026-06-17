#!/usr/bin/env bash
set -euo pipefail

SAMPLESHEET="$1"
OUTDIR="$2"
THREADS="${3:-6}"
TMPDIR="$2/tmp"
CHOPPER_QUALITY=10
CHOPPER_MINLENGTH=200
RSCRIPT="$HOME/.conda/R_shared/bin/Rscript"

LOGFILE="${OUTDIR}/pipeline.log"
mkdir -p "$OUTDIR"
mkdir -p "$TMPDIR"

exec > >(tee -a "$LOGFILE") 2>&1

echo "Pipeline start: $(date)"

FILTERED_ROOT="${OUTDIR}/filtered_reads"
mkdir -p "$FILTERED_ROOT"

# Resolve tool binaries from their conda environments
NANOPLOT_BIN="$(conda run -n nanoplot which NanoPlot 2>/dev/null)" || { echo "ERROR: NanoPlot not found in conda env 'nanoplot'" >&2; exit 1; }
CHOPPER_BIN="$(conda run -n chopper which chopper 2>/dev/null)" || { echo "ERROR: chopper not found in conda env 'chopper'" >&2; exit 1; }


# Read and verify the samplesheet
R_OUT=$(
  "$RSCRIPT" scripts/parse_samplesheet.R \
    "$SAMPLESHEET" \
    "$TMPDIR"
)

eval "$R_OUT"

TMPFILE="$LOOP_FILE"
VC_SOURCE="$VC_SOURCE"

TMPFILE="${TMPDIR}/tmp_samplesheet.tsv"

# If tmp_samplesheet.tsv has been written to TMPDIR, continue
if [[ ! -f "$TMPFILE" ]]; then
  echo "Error: Temporary samplesheet for Nanoplot and Chopper not found. Exiting."
  exit 1
fi

########################################################
# Run NanoPlot and Chopper for each sample
########################################################

echo ""
echo "First running NanoPlot and Chopper on each sample"
echo ""

# Read the tsv file in TMPDIR and loop over
while IFS=$'\t' read -r Sample InputDir; do
  echo ""
  echo "Processing sample: $Sample"
  echo "Input directory: $InputDir"

  NANOPLOT_OUT="${OUTDIR}/Nanoplot/${Sample}"
  QC_DIR_RAW="${NANOPLOT_OUT}/raw"
  QC_DIR_FILTERED="${NANOPLOT_OUT}/filtered"
  mkdir -p "$QC_DIR_RAW"
  mkdir -p "$QC_DIR_FILTERED"

  FILTERED_FASTQ="${FILTERED_ROOT}/${Sample}.fastq.gz"

  echo "Running NanoPlot on the raw reads"
  "$NANOPLOT_BIN" \
    --fastq "$InputDir"/*.fastq.gz \
    --outdir "$QC_DIR_RAW" \
    --threads "$THREADS"

  echo "Running Chopper on the raw reads: Q${CHOPPER_QUALITY}, minlength ${CHOPPER_MINLENGTH}"
  zcat "$InputDir"/*.fastq.gz | \
    "$CHOPPER_BIN" \
        --quality "$CHOPPER_QUALITY" \
        --minlength "$CHOPPER_MINLENGTH" | \
    gzip > "$FILTERED_FASTQ"

  echo "Running NanoPlot on the filtered reads"
  "$NANOPLOT_BIN" \
    --fastq "$FILTERED_FASTQ" \
    --outdir "$QC_DIR_FILTERED" \
    --threads "$THREADS"

done < "$TMPFILE"

########################################################
# Run ViroConstrictor on all samples
########################################################

VC_DIR="${OUTDIR}/viroconstrictor"
mkdir -p "$VC_DIR"

# First validate which filtered fastqs are non-empty and can be included in the ViroConstrictor samplesheet. This is a safeguard against empty or missing files that would cause ViroConstrictor to fail. We will write a new samplesheet with only valid samples for ViroConstrictor.
VC_SAMPLESHEET="${VC_DIR}/vc_samplesheet.tsv"
VALID_SAMPLES_FILE="${TMPDIR}/valid_samples.tsv"

echo ""
echo "Validating filtered FASTQs before running ViroConstrictor"

# detect non-empty gzip FASTQs
> "$VALID_SAMPLES_FILE"

while IFS=$'\t' read -r Sample InputDir; do
  FILTERED_FASTQ="${FILTERED_ROOT}/${Sample}.fastq.gz"

  if [[ -s "$FILTERED_FASTQ" ]]; then
    n_reads=$(zcat "$FILTERED_FASTQ" | awk 'END{print NR+0}')
    if [[ "$n_reads" -gt 0 ]]; then
      echo -e "${Sample}\t${InputDir}" >> "$VALID_SAMPLES_FILE"
    else
      echo "EXCLUDE (0 reads): $Sample"
    fi
  else
    echo "EXCLUDE (missing): $Sample"
  fi
done < "$TMPFILE"

valid_count=$(wc -l < "$VALID_SAMPLES_FILE" || true)

echo "Valid samples: $valid_count"

"$RSCRIPT" scripts/write_viroconstrictor_samplesheet.R \
  "$VALID_SAMPLES_FILE" \
  "$SAMPLESHEET" \
  "$VC_SAMPLESHEET"

# Run the ViroConstrictor pipeline if we have valid samples
echo ""
echo "Starting the ViroConstrictor pipeline on all samples"
echo "Input : $FILTERED_ROOT"
echo "Sheet : $VC_SAMPLESHEET"

conda run -n viroconstrictor ViroConstrictor \
  --input "$FILTERED_ROOT" \
  --output "${OUTDIR}/viroconstrictor" \
  --samplesheet "$VC_SAMPLESHEET" \
  --platform nanopore \
  --amplicon-type end-to-end \
  --threads "$THREADS" \
  --skip-updates

echo ""
echo "Pipeline finished: $(date)"
