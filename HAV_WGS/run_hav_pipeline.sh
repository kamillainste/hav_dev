#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# run_hav_pipeline.sh — HAV whole-genome consensus from ONT amplicon reads
#                       (HAVNET WGS scheme, genotypes IA / IB / IIIA)
#
#   raw reads -> NanoPlot -> chopper (Q10, >= 1000 bp) -> NanoPlot
#             -> ViroConstrictor 1.6.6 (--amplicon-type fragmented,
#                                       fastp --cut_right removed)
#             -> consensus FASTA, coverage plots, summary table
#
# Run setup.sh once first. See README.md for installation and docs/METHOD.md
# for why these settings were chosen.
#
# Usage:
#   bash run_hav_pipeline.sh <output_dir> <samplesheet.tsv>
#         every samplesheet row carries a path_to_run_folder
#   bash run_hav_pipeline.sh <runs_dir> <output_dir> <samplesheet.tsv>
#         rows with an empty path_to_run_folder fall back to <runs_dir>
#
#   samplesheet — tab-separated. Columns are matched BY HEADER NAME, so order
#                 does not matter and extra columns are ignored:
#                   sample   (required)  sample identifier
#                   barcode  (required)  barcodeNN
#                   run      (required)  run folder name; also the output folder
#                   ct, gt, kommentar    (optional) carried into the summary
#                   path_to_run_folder   (optional) parent directory holding the
#                                        run folders; falls back to <runs_dir>
#
# Environment variables (all optional):
#   THREADS          threads per tool                     (default 6)
#   HAV_TMPDIR       apptainer scratch space; needs ~1-2 GB free per running
#                    container                            (default $TMPDIR or /tmp)
#   HAV_ENV_PREFIX   prefix of the conda env names; must match setup.sh (default hav_)
#
# RAW DATA LAYOUT. Nanopore writes the barcode folders two unpredictable levels
# below the run folder, so reads are looked for at:
#
#     <path_to_run_folder>/<run>/*/*/fastq_pass/<barcode>/*.fastq.gz
#
# The two "*" levels are the instrument sample-id and run-id directories. If that
# glob matches more than one barcode directory the sample is an ERROR and is
# skipped: sample+barcode+run is required to be unique, so an ambiguous match
# means the raw data does not agree with the samplesheet. The flat layout
# <path>/<run>/<barcode>/ is accepted as a fallback and reported as such.
#
# The samplesheet drives the run: only listed rows are processed. Rows whose
# directory is missing are reported as errors rather than silently analysed.
#
# Output naming: all files use SAMPLENAME_barcode_run as prefix, under
# <output_dir>/<run>/. Each run gets its own summary_table.tsv, consensuses/ and
# coverage_plots/; a combined summary_table_all_runs.tsv is written at the top
# level. ViroConstrictor keeps its intermediate files under
# <output_dir>/<run>/<prefix>/data/ (alignment, clipper, fastp and AmpliGone
# output, BAM, VCF, per-position depth).
#
# summary_table.tsv holds exactly one row per samplesheet row, whether the
# sample was analysed, skipped or failed. The "status" column says which;
# numeric columns are NA for anything that did not produce a result.
#
# No single sample can abort the run. Every external tool call is guarded, and
# samples with fewer than MIN_READS reads are skipped before NanoPlot sees them
# (NanoPlot exits 1 on an empty FASTQ and declares <5 reads unsupported).
#
# ViroConstrictor: "--target HAV" selects its built-in HEPATOVIRUS preset:
#     minimap2  --frag=no, --splice on the raw alignment
#     Clipper   --exclude-spliced --spliced-length-threshold 200
#     AmpliGone --virtual-primers
# plus the two changes from stock ViroConstrictor (docs/METHOD.md):
#   * --amplicon-type fragmented: AmpliGone's lookaround is 10 instead of the
#     10000 that end-to-end uses, which cuts read ends at primers far away.
#   * fastp --cut_right removed, via the patched package copy that setup.sh
#     builds in .build/. fastp still runs, but only drops reads < 100 bp.
#
# Chopper: Q10, minlength 1000. 1000 rather than 1500 because IIIA amplicon 4
# is only 1261 bp: at 1500 that amplicon was filtered out entirely, leaving a
# ~480 bp hole at 4440-4915. Do not raise it above ~1200 without checking the
# IIIA amplicon-4 span first.
#
# NEGATIVE CONTROLS: nothing here removes low-level carry-over. An NTC can get
# a partial (mosaic) consensus from reads of other samples on the flowcell.
# Check NTC read counts and consensuses by hand.
# ---------------------------------------------------------------------------

PIPELINE_VERSION="1.0.0"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REFERENCE="${SCRIPT_DIR}/assets/HAV_refs/IA_IB_IIIA_seqs_fixed.fasta"
PRIMERS="${SCRIPT_DIR}/assets/HAV_primers/HAVNET_wgs_primers_fixed.fasta"
PLOT_SCRIPT="${SCRIPT_DIR}/scripts/plot_coverage.R"
RELABEL_SCRIPT="${SCRIPT_DIR}/scripts/relabel_amplicons.py"

HAV_ENV_PREFIX="${HAV_ENV_PREFIX:-hav_}"
CONDA_ENV="${HAV_ENV_PREFIX}viroconstrictor"
THREADS="${THREADS:-6}"
MIN_COVERAGE=20
PRIMER_MISMATCH_RATE=0.1
CHOPPER_QUALITY=10
CHOPPER_MINLENGTH=1000
AMPLICON_TYPE="fragmented"

# Samples with fewer than this many reads (raw, or surviving chopper) are
# skipped rather than analysed. NanoPlot exits 1 on an empty FASTQ and prints
# "WARNING: this code is not intended for such small datasets" below 5 reads,
# so 5 keeps us out of the range the tool itself disclaims.
MIN_READS=5

# ---------------------------------------------------------------------------
# Argument parsing & validation
# ---------------------------------------------------------------------------

usage() {
    cat >&2 <<USAGE
Usage:
  $(basename "$0") <output_dir> <samplesheet.tsv>
        when every samplesheet row carries a path_to_run_folder

  $(basename "$0") <runs_dir> <output_dir> <samplesheet.tsv>
        when rows have no path_to_run_folder; <runs_dir> is the fallback

See README.md. Run setup.sh once before the first run.
USAGE
    exit 1
}

case "$#" in
    2)  RUNS_DIR=""
        OUTPUT_DIR="${1%/}"
        METADATA="$2" ;;
    3)  RUNS_DIR="${1%/}"
        OUTPUT_DIR="${2%/}"
        METADATA="$3"
        [[ -d "$RUNS_DIR" ]] || { echo "ERROR: runs_dir does not exist: $RUNS_DIR" >&2; exit 1; } ;;
    *)  usage ;;
esac
[[ -f "$METADATA" ]]  || { echo "ERROR: samplesheet not found: $METADATA" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "ERROR: reference file not found: $REFERENCE" >&2; exit 1; }
[[ -f "$PRIMERS" ]]   || { echo "ERROR: primers file not found: $PRIMERS" >&2; exit 1; }
# Snakemake shell commands inside ViroConstrictor do not quote paths.
case "$OUTPUT_DIR" in *" "*)
    echo "ERROR: output_dir must not contain spaces: $OUTPUT_DIR" >&2; exit 1 ;; esac

# ---------------------------------------------------------------------------
# Container scratch space. Every `apptainer exec` unpacks the .sif into a
# temporary rootfs (~1.1 GB) and deletes it when the container exits. Point
# HAV_TMPDIR at a disk with room if /tmp is small.
# ---------------------------------------------------------------------------
export APPTAINER_TMPDIR="${HAV_TMPDIR:-${TMPDIR:-/tmp}}"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
mkdir -p "$APPTAINER_TMPDIR"

# ---------------------------------------------------------------------------
# Environment checks: fail early with a pointer to setup.sh.
# ---------------------------------------------------------------------------
setup_hint() { echo "       Run setup.sh first (see README.md)." >&2; exit 1; }
command -v conda >/dev/null 2>&1 || { echo "ERROR: conda not found." >&2; setup_hint; }

# ViroConstrictor stops and asks questions when its profile is missing.
VC_PROFILE="${HOME}/.ViroConstrictor_defaultprofile.ini"
grep -qs '^repro_method *= *containers' "$VC_PROFILE" \
    || { echo "ERROR: $VC_PROFILE is missing or does not use containers." >&2; setup_hint; }

# ---------------------------------------------------------------------------
# Patched ViroConstrictor. setup.sh copies the installed 1.6.6 package to
# .build/ViroConstrictor and blanks the three Fastp_QC_SlidingWindow_Settings_*
# values (--cut_right ...) in preset_params.json
# (patches/viroconstrictor-1.6.6_no_cut_right.diff). PYTHONPATH makes
# `import ViroConstrictor` resolve to it, and the Snakefile and preset JSON are
# then found relative to its __file__; `conda run` passes PYTHONPATH through.
# The installed environment is untouched. Version, import path and the setting
# itself are checked before anything runs.
# ---------------------------------------------------------------------------
VC_PATCH_DIR="${SCRIPT_DIR}/.build"
VC_VERSION_EXPECTED="1.6.6"
[[ -f "${VC_PATCH_DIR}/ViroConstrictor/workflow/helpers/preset_params.json" ]] \
    || { echo "ERROR: patched ViroConstrictor missing under ${VC_PATCH_DIR}" >&2; setup_hint; }
export PYTHONPATH="${VC_PATCH_DIR}${PYTHONPATH:+:$PYTHONPATH}"

VC_CHECK="$(conda run -n "$CONDA_ENV" python -c '
import importlib.metadata as m, json, os, ViroConstrictor as v
d = os.path.dirname(v.__file__)
p = json.load(open(os.path.join(d, "workflow/helpers/preset_params.json")))
print(m.version("ViroConstrictor"), d,
      repr(p["DEFAULT"]["STAGE_MAIN"]["Fastp_QC_SlidingWindow_Settings_nanopore"]))
' 2>/dev/null)" || { echo "ERROR: could not import ViroConstrictor in conda env '$CONDA_ENV'" >&2; setup_hint; }
read -r VC_VERSION VC_IMPORTED VC_CUTRIGHT <<< "$VC_CHECK"
if [[ "$VC_VERSION" != "$VC_VERSION_EXPECTED" ]]; then
    echo "ERROR: env '$CONDA_ENV' has ViroConstrictor $VC_VERSION; the patch is for $VC_VERSION_EXPECTED" >&2
    exit 1
fi
if [[ "$VC_IMPORTED" != "${VC_PATCH_DIR}/ViroConstrictor" || "$VC_CUTRIGHT" != "''" ]]; then
    echo "ERROR: patched ViroConstrictor not in effect" >&2
    echo "       imported from: $VC_IMPORTED   fastp sliding-window setting: $VC_CUTRIGHT" >&2
    setup_hint
fi

# Tools in the other envs are called by absolute path.
env_bin() { conda run -n "${HAV_ENV_PREFIX}$1" which "$2" 2>/dev/null; }
NANOPLOT_BIN="$(env_bin nanoplot NanoPlot)" \
    || { echo "ERROR: NanoPlot not found in conda env '${HAV_ENV_PREFIX}nanoplot'" >&2; setup_hint; }
CHOPPER_BIN="$(env_bin chopper chopper)" \
    || { echo "ERROR: chopper not found in conda env '${HAV_ENV_PREFIX}chopper'" >&2; setup_hint; }
VC_PYTHON="$(env_bin viroconstrictor python)" \
    || { echo "ERROR: python not found in conda env '$CONDA_ENV'" >&2; setup_hint; }

# Coverage plotting is optional: a missing Rscript or plot script disables the
# plots with a warning rather than failing a run that is otherwise fine.
PLOTS_ENABLED=1
RSCRIPT_BIN="$(env_bin plot Rscript || true)"
if [[ -z "$RSCRIPT_BIN" ]]; then
    echo "WARNING: Rscript not found in conda env '${HAV_ENV_PREFIX}plot' — coverage plots will be skipped" >&2
    PLOTS_ENABLED=0
elif [[ ! -f "$PLOT_SCRIPT" ]]; then
    echo "WARNING: plot script not found: $PLOT_SCRIPT — coverage plots skipped" >&2
    PLOTS_ENABLED=0
fi

# ---------------------------------------------------------------------------
# Load sample rows from the samplesheet
#
# Columns are located BY HEADER NAME, so their order is free and unknown extra
# columns are ignored. Tabs are translated to \037 (unit separator) because tab
# is IFS whitespace: `IFS=$'\t' read` collapses consecutive tabs and would shift
# every field left on rows with an empty ct or gt.
# ---------------------------------------------------------------------------

HEADER_LINE="$(head -1 "$METADATA" | tr -d '\r')"
[[ -n "$HEADER_LINE" ]] || { echo "ERROR: samplesheet is empty: $METADATA" >&2; exit 1; }

declare -A COL=()
idx=0
while IFS= read -r name; do
    idx=$(( idx + 1 ))
    name="$(printf '%s' "$name" | tr -d '\r' | sed 's/^[[:space:]]*//; s/[[:space:]]*$//' | tr 'A-Z' 'a-z')"
    [[ -n "$name" ]] && COL["$name"]=$idx
done < <(printf '%s' "$HEADER_LINE" | tr '\t' '\n')

for required in sample barcode run; do
    [[ -n "${COL[$required]:-}" ]] || {
        echo "ERROR: samplesheet is missing a '$required' column" >&2
        echo "       found: $(printf '%s' "$HEADER_LINE" | tr '\t' ' ')" >&2
        exit 1; }
done

# field <row> <column-name>  -> value, or empty when the column is absent
field() {
    local row="$1" name="$2" i="${COL[$2]:-}"
    [[ -n "$i" ]] || { printf ''; return; }
    printf '%s' "$row" | cut -d$'\037' -f"$i"
}

mapfile -t META_ROWS < <(
    tail -n +2 "$METADATA" | tr -d '\r' | grep -v '^[[:space:]]*$' | tr '\t' '\037'
)

if [[ ${#META_ROWS[@]} -eq 0 ]]; then
    echo "ERROR: no sample rows found in $METADATA" >&2; exit 1
fi

# ---------------------------------------------------------------------------
# Samplesheet validation
#
# Uniqueness is a property of the SAMPLESHEET, so it is checked over every row
# BEFORE any processing. Checking it inside the processing loop made the verdict
# depend on whether the earlier row happened to be processable: a duplicate
# slipped through whenever the row it duplicated was itself skipped for an
# unrelated reason (a missing path, a missing run folder).
# ---------------------------------------------------------------------------

declare -A CLAIMED_BY=()     # "barcode|run" → sample that claimed it first
declare -A ROW_DUP=()        # row index → the sample this row duplicates
_ri=0
for _row in "${META_ROWS[@]}"; do
    _ri=$(( _ri + 1 ))
    _b="$(field "$_row" barcode)"; _r="$(field "$_row" run)"; _s="$(field "$_row" sample)"
    [[ -n "$_b" && -n "$_r" ]] || continue
    if [[ -n "${CLAIMED_BY["${_b}|${_r}"]:-}" ]]; then
        ROW_DUP[$_ri]="${CLAIMED_BY["${_b}|${_r}"]}"
    else
        CLAIMED_BY["${_b}|${_r}"]="$_s"
    fi
done
unset _ri _row _b _r _s

echo "Found ${#META_ROWS[@]} sample row(s) in $(basename "$METADATA")"
echo "Columns     : $(printf '%s' "$HEADER_LINE" | tr '\t' ' ')"
if [[ -n "$RUNS_DIR" ]]; then
    echo "Runs base   : $RUNS_DIR  (fallback when path_to_run_folder is empty)"
else
    echo "Runs base   : none given — every row must carry path_to_run_folder"
fi
echo "Output base : $OUTPUT_DIR"
echo "Metadata    : $METADATA"
echo "Min reads   : $MIN_READS"
GIT_REV="$(git -C "$SCRIPT_DIR" describe --always --dirty 2>/dev/null || echo "not a git checkout")"
echo "Pipeline    : hav-consensus ${PIPELINE_VERSION} (${GIT_REV})"
SETTINGS="chopper -q ${CHOPPER_QUALITY} --minlength ${CHOPPER_MINLENGTH}; ViroConstrictor ${VC_VERSION} --amplicon-type ${AMPLICON_TYPE} --target HAV --match-ref --min-coverage ${MIN_COVERAGE} --primer-mismatch-rate ${PRIMER_MISMATCH_RATE}; fastp --cut_right removed (${VC_IMPORTED})"
echo "Settings    : $SETTINGS"
if [[ ${#ROW_DUP[@]} -gt 0 ]]; then
    echo "WARNING     : ${#ROW_DUP[@]} row(s) reuse a barcode+run already claimed — they will be skipped" >&2
fi
echo ""

mkdir -p "$OUTPUT_DIR"

# Provenance of this output folder: what made it, from what, with which settings.
{
    echo "pipeline  : hav-consensus ${PIPELINE_VERSION} (${GIT_REV})"
    echo "made by   : ${SCRIPT_DIR}/$(basename "$0")"
    echo "started   : $(date -Is)  on $(hostname)"
    echo "command   : $0 $*"
    echo "samplesheet: $METADATA  (md5 $(md5sum < "$METADATA" | cut -d' ' -f1))"
    echo "runs base : ${RUNS_DIR:-none (path_to_run_folder per row)}"
    echo "reference : $REFERENCE"
    echo "primers   : $PRIMERS"
    echo "settings  : $SETTINGS"
} > "${OUTPUT_DIR}/pipeline_run_info.txt"

# ---------------------------------------------------------------------------
# Summary tables
#
# One per run, under <output_dir>/<run>/, identical in shape to the previous
# single-run output; plus one combined table at the top level.
# ---------------------------------------------------------------------------

SUMMARY_ALL="${OUTPUT_DIR}/summary_table_all_runs.tsv"
SUMMARY_FMT='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'
SUMMARY_HDR=( "prefix" "sample" "ct" "gt" "barcode" "run"
              "raw_reads" "trimmed_reads"
              "avg_depth" "breadth_of_coverage_pct_at_mincov${MIN_COVERAGE}" "reference_id"
              "status" "kommentar" )

declare -A RUN_SUMMARY=()   # run -> its summary_table.tsv path

# summary_for <run>  -> sets CUR_SUMMARY to that run's table, creating it on
# first use. It sets a global rather than printing a path on purpose: capturing
# the path with $( ) would run this in a subshell, and the RUN_SUMMARY
# registration would be silently discarded along with it.
CUR_SUMMARY=""
summary_for() {
    local run="$1"
    CUR_SUMMARY="${OUTPUT_DIR}/${run}/summary_table.tsv"
    if [[ -z "${RUN_SUMMARY[$run]:-}" ]]; then
        mkdir -p "${OUTPUT_DIR}/${run}" "${OUTPUT_DIR}/${run}/consensuses"
        # shellcheck disable=SC2059
        printf "$SUMMARY_FMT" "${SUMMARY_HDR[@]}" > "$CUR_SUMMARY"
        RUN_SUMMARY[$run]="$CUR_SUMMARY"
    fi
}

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# count_reads <fastq.gz> [<fastq.gz> ...]  →  total reads across all files
# Counted straight from the FASTQ rather than from NanoPlot, so the number is
# available even when NanoPlot is skipped or fails — and so the too-few-reads
# check below can run *before* NanoPlot is invoked.
count_reads() {
    local total
    total="$( { zcat -f -- "$@" 2>/dev/null || true; } | awk 'END { print int(NR/4) }' )"
    printf '%s\n' "${total:-0}"
}

# calc_depth_breadth <coverage_tsv> <min_coverage>
# Prints "avg_depth\tbreadth_pct" from a two-column (pos, depth) TSV.
calc_depth_breadth() {
    local tsv="$1" mincov="$2"
    awk -v mc="$mincov" '
        { sum += $2; n++; if ($2 >= mc) cov++ }
        END {
            if (n > 0) printf "%.2f\t%.2f\n", sum/n, (cov/n)*100
            else       printf "NA\tNA\n"
        }' "$tsv"
}

# make_coverage_plots <results_dir> <output_dir> <label>
# Runs plot_coverage.R over results_dir, writing one PNG per genotype found.
# Never fatal: plotting is a reporting step and must not take down a run whose
# analysis already succeeded. --vanilla keeps R's site-init chatter out of the
# log without affecting which libraries are visible.
make_coverage_plots() {
    local results_dir="$1" out_dir="$2" label="$3"
    [[ "$PLOTS_ENABLED" -eq 1 ]] || return 1
    if "$RSCRIPT_BIN" --vanilla "$PLOT_SCRIPT" \
            "$results_dir" "$out_dir" "$MIN_COVERAGE"; then
        return 0
    fi
    echo "  WARNING: coverage plotting failed for ${label}" >&2
    return 1
}

# write_row <prefix> <raw> <trimmed> <avg_depth> <breadth> <ref_id> <status>
# Emits one summary line into the current row's run table. The per-sample
# metadata fields come from the row being processed, so every samplesheet row
# yields exactly one line.
write_row() {
    summary_for "$RUN_NAME"
    # shellcheck disable=SC2059
    printf "$SUMMARY_FMT" \
        "$1" "$SAMPLE" "$CT" "$GT" "$BARCODE" "$RUN_NAME" \
        "$2" "$3" "$4" "$5" "$6" "$7" "$KOMMENTAR" \
        >> "$CUR_SUMMARY"
}

# ---------------------------------------------------------------------------
# Main loop: one pass over the metadata rows
# ---------------------------------------------------------------------------

TOTAL=0
PROCESSED=0
SKIPPED=0
FAILED=0

declare -A RUN_WARNED=()     # run name → already warned that its folder is missing
declare -A CLAIMED=()        # "run|barcode" → listed in the samplesheet
declare -A RUN_BASE=()       # run → the path_to_run_folder it was found under

for row in "${META_ROWS[@]}"; do
    SAMPLE="$(field "$row" sample)"
    CT="$(field "$row" ct)"
    GT="$(field "$row" gt)"
    BARCODE="$(field "$row" barcode)"
    RUN_NAME="$(field "$row" run)"
    KOMMENTAR="$(field "$row" kommentar)"
    ROW_BASE="$(field "$row" path_to_run_folder)"
    ROW_BASE="${ROW_BASE:-$RUNS_DIR}"
    ROW_BASE="${ROW_BASE%/}"
    NO_BASE=0
    [[ -z "$ROW_BASE" ]] && NO_BASE=1

    TOTAL=$(( TOTAL + 1 ))

    if [[ -z "$SAMPLE" || -z "$BARCODE" || -z "$RUN_NAME" ]]; then
        echo "ERROR: row ${TOTAL} has an empty sample/barcode/run field — skipping" >&2
        write_row "NA" "NA" "NA" "NA" "NA" "NA" "SKIPPED: incomplete samplesheet row"
        SKIPPED=$(( SKIPPED + 1 ))
        continue
    fi


    PREFIX="${SAMPLE}_${BARCODE}_${RUN_NAME}"
    CLAIMED["${RUN_NAME}|${BARCODE}"]=1
    RUN_BASE["$RUN_NAME"]="$ROW_BASE"

    # Duplicate barcode+run, decided in the validation pass above so the verdict
    # does not depend on whether the row it duplicates was processable.
    if [[ -n "${ROW_DUP[$TOTAL]:-}" ]]; then
        echo "ERROR: barcode $BARCODE in run $RUN_NAME is already claimed by ${ROW_DUP[$TOTAL]} — skipping $SAMPLE" >&2
        write_row "$PREFIX" "NA" "NA" "NA" "NA" "NA" "SKIPPED: duplicate barcode+run in samplesheet"
        SKIPPED=$(( SKIPPED + 1 ))
        continue
    fi

    if [[ "$NO_BASE" -eq 1 ]]; then
        echo "ERROR: $SAMPLE has an empty path_to_run_folder and no <runs_dir> fallback — skipping" >&2
        write_row "$PREFIX" "NA" "NA" "NA" "NA" "NA" \
                  "SKIPPED: no path_to_run_folder and no runs_dir fallback"
        SKIPPED=$(( SKIPPED + 1 ))
        continue
    fi

    run_dir="${ROW_BASE}/${RUN_NAME}"
    if [[ ! -d "$run_dir" ]]; then
        if [[ -z "${RUN_WARNED["$RUN_NAME"]:-}" ]]; then
            echo "ERROR: run folder not found: $run_dir" >&2
            RUN_WARNED["$RUN_NAME"]=1
        fi
        echo "  skipping $SAMPLE ($BARCODE, $RUN_NAME): run folder missing" >&2
        write_row "$PREFIX" "NA" "NA" "NA" "NA" "NA" "SKIPPED: run folder not found"
        SKIPPED=$(( SKIPPED + 1 ))
        continue
    fi

    # Nanopore layout: <run>/<sample-id>/<run-id>/fastq_pass/<barcode>
    # The two middle levels are instrument-assigned, so they are globbed.
    mapfile -t BC_HITS < <(
        find "$run_dir" -mindepth 4 -maxdepth 4 -type d \
             -path "*/fastq_pass/${BARCODE}" 2>/dev/null | sort
    )
    LAYOUT="fastq_pass"
    if [[ ${#BC_HITS[@]} -eq 0 ]]; then
        # fall back to the older flat layout <run>/<barcode>
        if [[ -d "${run_dir}/${BARCODE}" ]]; then
            BC_HITS=( "${run_dir}/${BARCODE}" )
            LAYOUT="flat (legacy)"
        fi
    fi

    if [[ ${#BC_HITS[@]} -eq 0 ]]; then
        echo "ERROR: no ${BARCODE} directory under $run_dir" >&2
        echo "       looked for */*/fastq_pass/${BARCODE} and ${BARCODE}" >&2
        write_row "$PREFIX" "NA" "NA" "NA" "NA" "NA" "SKIPPED: barcode folder not found"
        SKIPPED=$(( SKIPPED + 1 ))
        continue
    fi
    if [[ ${#BC_HITS[@]} -gt 1 ]]; then
        echo "ERROR: $BARCODE is ambiguous in run $RUN_NAME — ${#BC_HITS[@]} directories match:" >&2
        printf '         %s\n' "${BC_HITS[@]}" >&2
        echo "       sample+barcode+run must resolve to exactly one directory; skipping $SAMPLE" >&2
        write_row "$PREFIX" "NA" "NA" "NA" "NA" "NA" "SKIPPED: ambiguous barcode folder (${#BC_HITS[@]} matches)"
        SKIPPED=$(( SKIPPED + 1 ))
        continue
    fi
    barcode_dir="${BC_HITS[0]}"

    if ! compgen -G "${barcode_dir}/*.fastq.gz" > /dev/null; then
        echo "ERROR: no *.fastq.gz files in $barcode_dir" >&2
        echo "  skipping $SAMPLE ($BARCODE, $RUN_NAME)" >&2
        write_row "$PREFIX" "NA" "NA" "NA" "NA" "NA" "SKIPPED: no fastq.gz files"
        SKIPPED=$(( SKIPPED + 1 ))
        continue
    fi

    RUN_OUT="${OUTPUT_DIR}/${RUN_NAME}"
    summary_for "$RUN_NAME"                # creates <run>/ and <run>/consensuses/
    SAMPLE_OUT="${RUN_OUT}/${PREFIX}"

    echo "=========================================="
    echo "  Processing : $PREFIX  (${TOTAL}/${#META_ROWS[@]})"
    echo "    Input    : $barcode_dir  [${LAYOUT}]"
    echo "    Output   : $SAMPLE_OUT"
    echo "=========================================="

    QC_RAW_DIR="${SAMPLE_OUT}/qc_raw"
    FILTERED_DIR="${SAMPLE_OUT}/filtered"
    # Name the filtered FASTQ with PREFIX so ViroConstrictor uses it as the
    # sample identifier — this propagates to all output filenames and the
    # consensus FASTA header automatically.
    FILTERED_FASTQ="${FILTERED_DIR}/${PREFIX}.fastq.gz"
    QC_FILTERED_DIR="${SAMPLE_OUT}/qc_filtered"

    mkdir -p "$QC_RAW_DIR" "$FILTERED_DIR" "$QC_FILTERED_DIR"

    # --- Gate 1: enough raw reads to be worth processing? ---
    RAW_READS="$(count_reads "${barcode_dir}"/*.fastq.gz)"

    if [[ "$RAW_READS" -lt "$MIN_READS" ]]; then
        echo "  SKIP: only ${RAW_READS} raw read(s), need ${MIN_READS}" >&2
        write_row "$PREFIX" "$RAW_READS" "NA" "NA" "NA" "NA" \
            "SKIPPED: too few raw reads (${RAW_READS} < ${MIN_READS})"
        SKIPPED=$(( SKIPPED + 1 ))
        continue
    fi

    # --- Step 1: NanoPlot on raw reads (non-fatal) ---
    echo "  [1/4] NanoPlot QC (pre-filter)"
    if ! "$NANOPLOT_BIN" \
            --fastq "${barcode_dir}"/*.fastq.gz \
            --outdir "$QC_RAW_DIR" \
            --threads "$THREADS"; then
        echo "  WARNING: NanoPlot (pre-filter) failed for $PREFIX — continuing" >&2
    fi

    # --- Step 2: Chopper filtering ---
    echo "  [2/4] Chopper  Q${CHOPPER_QUALITY}  minlen ${CHOPPER_MINLENGTH}"
    if ! cat "${barcode_dir}"/*.fastq.gz \
            | gunzip -c \
            | "$CHOPPER_BIN" \
                --quality "$CHOPPER_QUALITY" \
                --minlength "$CHOPPER_MINLENGTH" \
            | gzip > "$FILTERED_FASTQ"; then
        echo "  ERROR: chopper failed for $PREFIX — skipping sample" >&2
        write_row "$PREFIX" "$RAW_READS" "NA" "NA" "NA" "NA" "FAILED: chopper error"
        FAILED=$(( FAILED + 1 ))
        continue
    fi

    # --- Gate 2: enough reads left after trimming? ---
    TRIMMED_READS="$(count_reads "$FILTERED_FASTQ")"
    echo "  Reads    : ${RAW_READS} raw → ${TRIMMED_READS} after trimming"

    if [[ "$TRIMMED_READS" -lt "$MIN_READS" ]]; then
        echo "  SKIP: only ${TRIMMED_READS} read(s) survived trimming, need ${MIN_READS}" >&2
        echo "        (skipping NanoPlot post-filter and ViroConstrictor)" >&2
        write_row "$PREFIX" "$RAW_READS" "$TRIMMED_READS" "NA" "NA" "NA" \
            "SKIPPED: too few reads after trimming (${TRIMMED_READS} < ${MIN_READS})"
        SKIPPED=$(( SKIPPED + 1 ))
        continue
    fi

    # --- Step 3: NanoPlot on filtered reads (non-fatal) ---
    echo "  [3/4] NanoPlot QC (post-filter)"
    if ! "$NANOPLOT_BIN" \
            --fastq "$FILTERED_FASTQ" \
            --outdir "$QC_FILTERED_DIR" \
            --threads "$THREADS"; then
        echo "  WARNING: NanoPlot (post-filter) failed for $PREFIX — continuing" >&2
    fi

    # --- Step 4: ViroConstrictor ---
    echo "  [4/4] ViroConstrictor"
    if ! conda run -n "$CONDA_ENV" ViroConstrictor \
            --input               "$FILTERED_DIR" \
            --output              "$SAMPLE_OUT" \
            --reference           "$REFERENCE" \
            --primers             "$PRIMERS" \
            --features            NONE \
            --platform            nanopore \
            --amplicon-type       "$AMPLICON_TYPE" \
            --target              HAV \
            --match-ref \
            --min-coverage        "$MIN_COVERAGE" \
            --primer-mismatch-rate "$PRIMER_MISMATCH_RATE" \
            --threads             "$THREADS" \
            --skip-updates; then
        echo "  WARNING: ViroConstrictor failed for $PREFIX — recording NA stats" >&2
        write_row "$PREFIX" "$RAW_READS" "$TRIMMED_READS" "NA" "NA" "NA" \
            "FAILED: ViroConstrictor error"
        FAILED=$(( FAILED + 1 ))
        continue
    fi

    # --- Post-processing: collect stats ---

    # Reset per-sample state explicitly: these persist across loop iterations,
    # so a sample that skips a branch below would otherwise inherit the
    # previous sample's reference, depth or coverage path.
    STATUS="OK"
    REF_ID="NA"
    AVG_DEPTH="NA"
    BREADTH="NA"
    COVERAGE_TSV=""

    # Find the chosen reference directory (RefID~<id>).
    # The -d test is load-bearing: find exits 1 on a missing directory and
    # pipefail propagates that through `sort | head`, so without it a
    # ViroConstrictor run that produced no data/Virus~HAV would abort the whole
    # script instead of falling through to the warning below.
    REF_DIR=""
    if [[ -d "${SAMPLE_OUT}/data/Virus~HAV" ]]; then
        REF_DIR="$(find "${SAMPLE_OUT}/data/Virus~HAV" \
                        -maxdepth 1 -type d -name "RefID~*" 2>/dev/null \
                   | sort | head -1 || true)"
    fi

    if [[ -z "$REF_DIR" ]]; then
        echo "  WARNING: no RefID~ directory found for $PREFIX" >&2
        STATUS="OK: no RefID directory"
    else
        REF_ID="$(basename "$REF_DIR" | sed 's/^RefID~//')"
        COVERAGE_TSV="${REF_DIR}/consensus/coverages/${PREFIX}_coverage.tsv"

        if [[ -f "$COVERAGE_TSV" ]]; then
            read -r AVG_DEPTH BREADTH < <(
                calc_depth_breadth "$COVERAGE_TSV" "$MIN_COVERAGE"
            )
        else
            echo "  WARNING: coverage TSV not found: $COVERAGE_TSV" >&2
            COVERAGE_TSV=""
            STATUS="OK: no coverage file"
        fi

        # --- Relabel amplicon coverage CSV columns by amplicon NUMBER ---
        # ViroConstrictor names every amplicon after the lowest-coordinate primer
        # in the BED, so an IIIA sample's columns can read HAVIB_*. See
        # scripts/relabel_amplicons.py; it is non-fatal and notes any mismatch in
        # <prefix>_ampliconcoverage_rename_notes.tsv.
        AMPLICON_COV="${REF_DIR}/primers/${PREFIX}_ampliconcoverage.csv"
        if [[ -f "$AMPLICON_COV" ]]; then
            GENOTYPE="${REF_ID##*_}"
            if ! "$VC_PYTHON" "$RELABEL_SCRIPT" "$AMPLICON_COV" "$GENOTYPE"; then
                echo "  WARNING: amplicon coverage rename failed for $PREFIX" >&2
            fi
        fi
    fi

    # --- Rename consensus: consensus.fasta → <prefix>.fasta ---
    CONSENSUS_SRC="${SAMPLE_OUT}/results/combined/by_sample/${PREFIX}/consensus.fasta"
    CONSENSUS_DST="${SAMPLE_OUT}/results/combined/by_sample/${PREFIX}/${PREFIX}.fasta"
    if [[ -f "$CONSENSUS_SRC" ]]; then
        cp "$CONSENSUS_SRC" "$CONSENSUS_DST"
        cp "$CONSENSUS_DST" "${RUN_OUT}/consensuses/${PREFIX}.fasta"
        echo "  Consensus : $(basename "$CONSENSUS_DST")"
        PROCESSED=$(( PROCESSED + 1 ))
    else
        echo "  WARNING: consensus not found: $CONSENSUS_SRC" >&2
        STATUS="OK: no consensus produced"
        FAILED=$(( FAILED + 1 ))
    fi

    # --- Coverage plot for this sample alone ---
    # Scoped to SAMPLE_OUT, so the PNG contains only this sample's depth trace.
    if [[ -n "$COVERAGE_TSV" ]]; then
        echo "  [plot] per-sample coverage"
        make_coverage_plots "$SAMPLE_OUT" "${SAMPLE_OUT}/coverage_plots" "$PREFIX" || true
    fi

    # ct / gt / kommentar come straight from the metadata row being processed
    write_row "$PREFIX" "$RAW_READS" "$TRIMMED_READS" \
        "$AVG_DEPTH" "$BREADTH" "$REF_ID" "$STATUS"

    echo "  Done: avg_depth=${AVG_DEPTH}x  breadth=${BREADTH}%  ref=${REF_ID}"
    echo ""
done

# ---------------------------------------------------------------------------
# Per-run coverage plots
# ---------------------------------------------------------------------------

PLOTS_MADE=0
ALL_PLOTS_MADE=0
if [[ "$PLOTS_ENABLED" -eq 1 && ${#RUN_SUMMARY[@]} -gt 0 ]]; then
    echo "=========================================="
    echo "Coverage plots"
    echo "=========================================="
    for run_name in $(printf '%s\n' "${!RUN_SUMMARY[@]}" | sort); do
        # A run every one of whose rows was skipped has a summary table but no
        # sample directories; plotting it would only create an empty folder and
        # emit a failure warning, so skip it.
        if ! compgen -G "${OUTPUT_DIR}/${run_name}/*/data" > /dev/null; then
            echo "  (no analysed samples in ${run_name} — no coverage plot)"
            continue
        fi
        # Scoped to one run, so each genotype PNG overlays only that run's samples.
        if make_coverage_plots "${OUTPUT_DIR}/${run_name}" \
                               "${OUTPUT_DIR}/${run_name}/coverage_plots" "$run_name"; then
            PLOTS_MADE=$(( PLOTS_MADE + 1 ))
        fi
    done

    # And one across every run together. Scoped to OUTPUT_DIR, which the plot
    # script walks recursively, so a genotype PNG here overlays every sample of
    # that genotype from all runs. The per-run coverage_plots/ folders hold only
    # PNGs and are ignored, since the scan matches *_coverage.tsv.
    if make_coverage_plots "$OUTPUT_DIR" "${OUTPUT_DIR}/coverage_plots" "all runs combined"; then
        ALL_PLOTS_MADE=1
    fi
    echo ""
fi

# ---------------------------------------------------------------------------
# Combined summary across every run
# ---------------------------------------------------------------------------

if [[ ${#RUN_SUMMARY[@]} -gt 0 ]]; then
    # shellcheck disable=SC2059
    printf "$SUMMARY_FMT" "${SUMMARY_HDR[@]}" > "$SUMMARY_ALL"
    for run_name in $(printf '%s\n' "${!RUN_SUMMARY[@]}" | sort); do
        tail -n +2 "${RUN_SUMMARY[$run_name]}" >> "$SUMMARY_ALL"
    done
fi

echo "=========================================="
echo "Finished."
echo "  Rows in metadata : ${TOTAL}"
echo "  Consensus built  : ${PROCESSED}"
echo "  Failed           : ${FAILED}"
echo "  Skipped          : ${SKIPPED}"
echo "  Runs in sheet    : ${#RUN_SUMMARY[@]}"
echo "  Results          : $OUTPUT_DIR"
for run_name in $(printf '%s\n' "${!RUN_SUMMARY[@]}" | sort); do
    echo "      ${run_name}/  ->  $(( $(wc -l < "${RUN_SUMMARY[$run_name]}") - 1 )) row(s)"
done
[[ -f "$SUMMARY_ALL" ]] && echo "  Combined summary : $SUMMARY_ALL"
if [[ "$PLOTS_MADE" -gt 0 ]]; then
    echo "  Coverage plots   : ${PLOTS_MADE} run(s), in <output>/<run>/coverage_plots"
fi
if [[ "$ALL_PLOTS_MADE" -eq 1 ]]; then
    echo "  Combined plots   : ${OUTPUT_DIR}/coverage_plots  (all runs together)"
fi

# Report barcode directories that exist on disk but are absent from the
# samplesheet. Only runs named in the samplesheet are walked - a run nobody
# listed is not this script's business - but every barcode inside them is.
UNLISTED=()
for run_name in "${!RUN_BASE[@]}"; do
    run_dir="${RUN_BASE[$run_name]}/${run_name}"
    [[ -d "$run_dir" ]] || continue
    while IFS= read -r bc_dir; do
        bc="$(basename "$bc_dir")"
        if [[ -z "${CLAIMED["${run_name}|${bc}"]:-}" ]]; then
            UNLISTED+=( "${run_name}/${bc}" )
        fi
    done < <(
        { find "$run_dir" -mindepth 4 -maxdepth 4 -type d -path "*/fastq_pass/barcode*" 2>/dev/null
          find "$run_dir" -mindepth 1 -maxdepth 1 -type d -name "barcode*" 2>/dev/null ; } | sort -u
    )
done

if [[ ${#UNLISTED[@]} -gt 0 ]]; then
    echo ""
    echo "  NOTE: ${#UNLISTED[@]} barcode folder(s) on disk are not in the metadata"
    echo "        table and were NOT analysed:"
    for entry in "${UNLISTED[@]}"; do
        echo "          $entry"
    done
fi
echo "=========================================="
