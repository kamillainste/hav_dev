#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# setup.sh — one-time installation for hav-consensus. Safe to re-run.
#
# Usage:
#   bash setup.sh [--skip-containers]
#
# Steps:
#   1. create the conda environments in envs/ (existing ones are left as they are)
#   2. build the patched ViroConstrictor copy in .build/ from the installed
#      package + patches/viroconstrictor-1.6.6_no_cut_right.diff, checked by md5
#   3. make sure ~/.ViroConstrictor_defaultprofile.ini exists and uses containers
#      (without it ViroConstrictor stops and asks questions on its first run)
#   4. download ViroConstrictor's containers (~2.3 GB) and test that apptainer
#      can run them (skip with --skip-containers; the first pipeline run then
#      downloads them instead)
#
# Environment variables (all optional):
#   HAV_ENV_PREFIX     prefix of the conda env names          (default: hav_)
#   HAV_CONTAINER_DIR  container cache for a NEW profile      (default: ~/.viroconstrictor/containers)
#   HAV_TMPDIR         apptainer scratch space                (default: $TMPDIR or /tmp)
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HAV_ENV_PREFIX="${HAV_ENV_PREFIX:-hav_}"
VC_ENV="${HAV_ENV_PREFIX}viroconstrictor"
BUILD_DIR="${REPO_DIR}/.build"
PATCH_FILE="${REPO_DIR}/patches/viroconstrictor-1.6.6_no_cut_right.diff"
PROFILE="${HOME}/.ViroConstrictor_defaultprofile.ini"
NEW_CONTAINER_DIR="${HAV_CONTAINER_DIR:-${HOME}/.viroconstrictor/containers}"

# The patch is only valid for this exact file. md5 of preset_params.json in
# ViroConstrictor 1.6.6 (bioconda), before and after the patch.
VC_VERSION_EXPECTED="1.6.6"
MD5_ORIGINAL="dd7e4ac69dcd8724b7e3d5467ed24033"
MD5_PATCHED="7888fcbcc6dc8d03885bece08f77ce7f"

SKIP_CONTAINERS=0
for arg in "$@"; do
    case "$arg" in
        --skip-containers) SKIP_CONTAINERS=1 ;;
        -h|--help) sed -n '2,23p' "$0"; exit 0 ;;
        *) echo "ERROR: unknown option: $arg" >&2; exit 1 ;;
    esac
done

step() { printf '\n== %s\n' "$*"; }
die()  { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------
step "0/4 Checks"
[[ "$(uname -s)" == "Linux" ]] || die "run this inside Linux or WSL, not $(uname -s)"
case "$REPO_DIR" in
    /mnt/[a-z]/*) echo "WARNING: the repository is on a Windows drive ($REPO_DIR)."
                  echo "         Clone it inside WSL instead (e.g. ~/hav-consensus): Windows drives"
                  echo "         are slow from WSL and do not keep Linux file permissions." ;;
esac
if grep -lq $'\r' "$REPO_DIR"/*.sh "$REPO_DIR"/scripts/* 2>/dev/null; then
    die "scripts have Windows (CRLF) line endings. Clone the repository from inside WSL with git, not from Windows."
fi
command -v conda >/dev/null 2>&1 \
    || die "conda not found. Install Miniforge (https://conda-forge.org/download/), open a new terminal, re-run."
SOLVER="conda"
command -v mamba >/dev/null 2>&1 && SOLVER="mamba"
echo "repository : $REPO_DIR"
echo "conda      : $(conda --version)  (creating envs with $SOLVER)"
echo "env prefix : $HAV_ENV_PREFIX"

# ---------------------------------------------------------------------------
step "1/4 Conda environments"
# On linux-64 (incl. WSL) the envs are made from envs/lock/*.linux-64.txt: the
# exact package set the pipeline was tested with. The .yml files only pin the
# main tools; their dependencies drift over time (NanoPlot 1.32.1 broke on a
# newer matplotlib). Other platforms fall back to the .yml files.
env_exists() { conda env list | awk '{print $1}' | grep -qx "$1"; }
PLATFORM="$(uname -s)-$(uname -m)"
for yml in "$REPO_DIR"/envs/*.yml; do
    base="$(basename "$yml" .yml)"
    name="${HAV_ENV_PREFIX}${base}"
    lock="$REPO_DIR/envs/lock/${base}.linux-64.txt"
    if env_exists "$name"; then
        echo "$name: exists, left unchanged  (to rebuild: conda env remove -n $name, then re-run)"
    elif [[ "$PLATFORM" == "Linux-x86_64" && -f "$lock" ]]; then
        echo "$name: creating from envs/lock/$(basename "$lock") (exact tested versions) ..."
        conda create -q -y -n "$name" --file "$lock" \
            || die "could not create $name (see the conda output above)"
    else
        echo "$name: creating from envs/$(basename "$yml") ..."
        "$SOLVER" env create -q -y -n "$name" -f "$yml" \
            || die "could not create $name (see the conda output above)"
    fi
done

# ---------------------------------------------------------------------------
step "2/4 Patched ViroConstrictor (fastp --cut_right removed)"
# PYTHONPATH is cleared so this finds the INSTALLED package, not an older copy.
read -r VC_VERSION VC_INSTALLED < <(env -u PYTHONPATH conda run -n "$VC_ENV" python -c '
import importlib.metadata as m, os, ViroConstrictor as v
print(m.version("ViroConstrictor"), os.path.dirname(v.__file__))') \
    || die "cannot import ViroConstrictor in env $VC_ENV"
[[ "$VC_VERSION" == "$VC_VERSION_EXPECTED" ]] \
    || die "env $VC_ENV has ViroConstrictor $VC_VERSION; the patch is for $VC_VERSION_EXPECTED"
PRESET_REL="workflow/helpers/preset_params.json"
md5_of() { md5sum "$1" | cut -d' ' -f1; }
[[ "$(md5_of "$VC_INSTALLED/$PRESET_REL")" == "$MD5_ORIGINAL" ]] \
    || die "installed $PRESET_REL differs from the ViroConstrictor $VC_VERSION_EXPECTED file the patch was made for"

rm -rf "${BUILD_DIR}/ViroConstrictor"
mkdir -p "$BUILD_DIR"
cp -r "$VC_INSTALLED" "${BUILD_DIR}/ViroConstrictor"
find "${BUILD_DIR}/ViroConstrictor" -name __pycache__ -type d -prune -exec rm -rf {} +

conda run -n "$VC_ENV" python "${REPO_DIR}/scripts/apply_patch.py" \
    "$PATCH_FILE" "${BUILD_DIR}/ViroConstrictor/$PRESET_REL" \
    || die "could not apply $PATCH_FILE"
[[ "$(md5_of "${BUILD_DIR}/ViroConstrictor/$PRESET_REL")" == "$MD5_PATCHED" ]] \
    || die "patched $PRESET_REL does not have the expected md5"
echo "built      : ${BUILD_DIR}/ViroConstrictor  (from $VC_INSTALLED)"
echo "checked    : $PRESET_REL md5 $MD5_ORIGINAL -> $MD5_PATCHED"

# ---------------------------------------------------------------------------
step "3/4 ViroConstrictor profile"
# ViroConstrictor reads this fixed path. An existing profile belongs to the
# user and is never overwritten; it only has to use containers.
if [[ ! -s "$PROFILE" ]]; then
    cat > "$PROFILE" <<EOF
[COMPUTING]
compmode = local

[GENERAL]
auto_update = no
ask_for_update = no

[REPRODUCTION]
repro_method = containers
container_cache_path = ${NEW_CONTAINER_DIR}
EOF
    echo "written    : $PROFILE"
else
    echo "exists     : $PROFILE (left unchanged)"
fi
read -r PROFILE_OK REPRO CONTAINER_DIR < <(PYTHONPATH="$BUILD_DIR" conda run -n "$VC_ENV" python -c '
import configparser, sys
from ViroConstrictor.userprofile import AllOptionsGiven
c = configparser.ConfigParser(); c.read(sys.argv[1])
r = c.get("REPRODUCTION", "repro_method", fallback="-")
print(AllOptionsGiven(c), r, c.get("REPRODUCTION", "container_cache_path", fallback="-"))' "$PROFILE")
[[ "$PROFILE_OK" == "True" ]] \
    || die "$PROFILE is incomplete; ViroConstrictor would stop to ask questions. Move it aside and re-run setup.sh."
[[ "$REPRO" == "containers" ]] \
    || die "$PROFILE uses repro_method=$REPRO; this pipeline was validated with containers. Change it to 'containers' (and set container_cache_path) or move the file aside and re-run."
echo "containers : $CONTAINER_DIR"

# ---------------------------------------------------------------------------
step "4/4 Containers"
APPTAINER_BIN="$(conda run -n "$VC_ENV" which apptainer 2>/dev/null)" \
    || die "apptainer not found in env $VC_ENV"
echo "apptainer  : $APPTAINER_BIN ($(conda run -n "$VC_ENV" apptainer --version))"
export APPTAINER_TMPDIR="${HAV_TMPDIR:-${TMPDIR:-/tmp}}"
if [[ "$SKIP_CONTAINERS" -eq 1 ]]; then
    echo "skipped (--skip-containers): the first pipeline run downloads them"
else
    mkdir -p "$CONTAINER_DIR"
    PYTHONPATH="$BUILD_DIR" conda run --no-capture-output -n "$VC_ENV" python -c '
import sys, pathlib
from ViroConstrictor.workflow.helpers.containers import download_containers
sys.exit(download_containers(pathlib.Path(sys.argv[1]), verbose=True))' "$CONTAINER_DIR" \
        || die "container download failed (network? disk space in $CONTAINER_DIR?)"
    CLEAN_SIF="$(ls "$CONTAINER_DIR"/viroconstrictor_clean_*.sif 2>/dev/null | head -1)"
    [[ -n "$CLEAN_SIF" ]] || die "no viroconstrictor_clean_*.sif in $CONTAINER_DIR after download"
    conda run -n "$VC_ENV" apptainer exec "$CLEAN_SIF" fastp --version >/dev/null 2>&1 \
        || die "apptainer cannot run $CLEAN_SIF. See README.md, Troubleshooting."
    echo "tested     : apptainer exec $(basename "$CLEAN_SIF") fastp --version  OK"
fi

{
    echo "setup.sh finished $(date -Is) on $(hostname)"
    echo "env prefix       : $HAV_ENV_PREFIX"
    echo "ViroConstrictor  : $VC_VERSION from $VC_INSTALLED"
    echo "profile          : $PROFILE (repro_method=$REPRO)"
    echo "container cache  : $CONTAINER_DIR"
    echo "apptainer        : $APPTAINER_BIN"
} > "${BUILD_DIR}/setup_info.txt"

step "Setup complete"
echo "Run the pipeline with:"
echo "  bash run_hav_pipeline.sh <output_dir> <samplesheet.tsv>"
