#!/usr/bin/env bash
set -euo pipefail

# Set up every conda environment the HAV pipeline needs:
#
#   1. HAVDEV   — project-local prefix env (./.conda/HAVDEV)
#                 Core bioinformatics tools + Python (environment.yml)
#   2. R_shared — shared prefix env ($HOME/.conda/R_shared)
#                 R + report/analysis packages (r_environment.yml)
#   3. WGS envs — named envs: nanoplot, chopper, viroconstrictor
#                 Nanopore assembly tools (envs/*.yml). Needed for --mode wgs.
#
# Usage:
#   bash scripts/setup_havdev_env.sh           # set up everything
#   bash scripts/setup_havdev_env.sh --core    # HAVDEV + R_shared only (skip WGS)
#
# Re-running updates existing environments in place (idempotent).

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

HAVDEV_PREFIX="$ROOT_DIR/.conda/HAVDEV"
HAVDEV_FILE="$ROOT_DIR/environment.yml"

RSHARED_PREFIX="$HOME/.conda/R_shared"
RSHARED_FILE="$ROOT_DIR/r_environment.yml"

WGS_ENV_DIR="$ROOT_DIR/envs"
WGS_ENVS=(nanoplot chopper viroconstrictor)

SETUP_WGS=1
for arg in "$@"; do
  case "$arg" in
    --core) SETUP_WGS=0 ;;
    -h|--help)
      sed -n '3,18p' "$0" | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "ERROR: unknown option: $arg" >&2; exit 1 ;;
  esac
done

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found in PATH. Install Miniconda/Miniforge first." >&2
  exit 1
fi

# Create or update a prefix-based environment from a yml file.
make_prefix_env() {
  local prefix="$1" file="$2" label="$3"
  if [[ ! -f "$file" ]]; then
    echo "ERROR: environment file not found: $file" >&2
    exit 1
  fi
  if [[ -d "$prefix" ]]; then
    echo "Updating $label at: $prefix"
    conda env update --prefix "$prefix" --file "$file" --prune
  else
    echo "Creating $label at: $prefix"
    conda env create --prefix "$prefix" --file "$file"
  fi
}

# Create or update a named environment from a yml file.
make_named_env() {
  local name="$1" file="$2"
  if [[ ! -f "$file" ]]; then
    echo "WARNING: env file not found, skipping '$name': $file" >&2
    return 0
  fi
  if conda env list | awk '{print $1}' | grep -qx "$name"; then
    echo "Updating env '$name'"
    conda env update -n "$name" --file "$file" --prune
  else
    echo "Creating env '$name'"
    conda env create -n "$name" --file "$file"
  fi
}

# ── 1. HAVDEV (core bioinformatics tools + Python) ────────────────────────────
echo "════════════════════════════════════════════════════════════════"
echo " 1/3  HAVDEV — core tools"
echo "════════════════════════════════════════════════════════════════"
mkdir -p "$ROOT_DIR/.conda"
make_prefix_env "$HAVDEV_PREFIX" "$HAVDEV_FILE" "HAVDEV (core tools)"

# ── 2. R_shared (R + report packages) ─────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════════"
echo " 2/3  R_shared — R + report packages"
echo "════════════════════════════════════════════════════════════════"
mkdir -p "$(dirname "$RSHARED_PREFIX")"
make_prefix_env "$RSHARED_PREFIX" "$RSHARED_FILE" "R_shared (R packages)"

# Safety net: install any R packages absent from conda channels into R_shared.
conda run --prefix "$RSHARED_PREFIX" Rscript - <<'RS'
req_cran <- c(
  "tidyverse", "seqinr", "rmarkdown", "knitr", "jsonlite", "ape",
  "phangorn", "kableExtra", "scales", "patchwork", "pheatmap",
  "DBI", "RSQLite"
)

req_bioc <- c("ggtree", "ggmsa", "msa", "Biostrings")

missing_cran <- req_cran[!vapply(req_cran, requireNamespace, logical(1), quietly = TRUE)]
missing_bioc <- req_bioc[!vapply(req_bioc, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (length(missing_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

cat("R package check complete.\n")
cat("  Missing CRAN before install:", paste(missing_cran, collapse = ", "), "\n")
cat("  Missing Bioc before install:", paste(missing_bioc, collapse = ", "), "\n")
RS

# ── 3. WGS assembly envs (named) ──────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════════"
echo " 3/3  WGS (Nanopore) environments"
echo "════════════════════════════════════════════════════════════════"
if [[ "$SETUP_WGS" -eq 1 ]]; then
  for env in "${WGS_ENVS[@]}"; do
    make_named_env "$env" "$WGS_ENV_DIR/${env}.yml"
  done
else
  echo "Skipping WGS environments (--core). Sanger mode only."
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════════"
echo "Setup complete."
echo "  HAVDEV   : $HAVDEV_PREFIX"
echo "             conda activate $HAVDEV_PREFIX"
echo "  R_shared : $RSHARED_PREFIX"
echo "             (used via \$HOME/.conda/R_shared/bin/Rscript)"
if [[ "$SETUP_WGS" -eq 1 ]]; then
  echo "  WGS envs : ${WGS_ENVS[*]}"
fi
echo "════════════════════════════════════════════════════════════════"
