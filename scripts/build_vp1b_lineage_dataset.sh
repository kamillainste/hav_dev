#!/usr/bin/env bash
# scripts/build_vp1b_lineage_dataset.sh
#
# Builds an enhanced local copy of the community hav-vp1-2b NextClade dataset
# that adds the FHI HAV lineage system (IA.1, IB.1.1, etc.) as a second
# clade-like attribute ("lineage_phylo"), without modifying the original broad
# genotype assignments ("clade_membership": IA, IB, ...).
#
# PUBLIC-ONLY DERIVATION
# ----------------------
# Lineages — both *which lineages exist* and *their defining mutations* — are
# designated directly on the public community hav-vp1-2b reference tree.  This
# script has NO dependency on data/local_datasets/ and no private FHI sequences
# feed into the lineage definitions.  The community dataset is fetched fresh from
# the Nextclade dataset server so the whole pipeline is reproducible from public
# sources alone.
#
# Because the public tree carries no outbreak labels, designation is purely
# criteria-based (monophyly + minimum size + >=1 defining mutation); the
# outbreak-seeding path of designate_lineages.py is intentionally NOT used here.
#
# Output dataset:  data/nextclade_datasets/hav-vp1-2b-lineages/
# Lineage register: data/lineage_designations.tsv  (derived from the public tree)
#
# Prerequisites:
#   conda activate <project>/.conda/HAVDEV   (nextclade, augur available)
#
# Usage (from project root):
#   bash scripts/build_vp1b_lineage_dataset.sh [options]
#
# Options:
#   --refresh                Force re-fetch of the community dataset even if present.
#   --no-fetch               Do not fetch; reuse an existing data/nextclade_datasets/hav-vp1-2b.
#   --genotypes "IA IB ..."  Genotypes to designate lineages for (default: "IA IB IIA IIIA").
#   --min-sequences <int>    Min leaf count per lineage (default: 3).
#   --max-depth <int>        Max lineage nesting depth (default: 3).

set -euo pipefail

PROJ_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SRC_DATASET="$PROJ_DIR/data/nextclade_datasets/hav-vp1-2b"
OUT_DATASET="$PROJ_DIR/data/nextclade_datasets/hav-vp1-2b-lineages"
WORK_DIR="$PROJ_DIR/data/nextclade_datasets/_work_vp1b_lineages"
LINEAGE_REGISTER="$PROJ_DIR/data/lineage_designations.tsv"
LINEAGE_PROPOSALS="$PROJ_DIR/data/lineage_proposals.tsv"
SCRIPTS="$PROJ_DIR/scripts"

COMMUNITY_NAME="community/masphl-bioinformatics/hav/vp1-2b-junction"

# ── Options ───────────────────────────────────────────────────────────────────
REFRESH=0
NO_FETCH=0
GENOTYPES="IA IB IIA IIIA"
MIN_SEQUENCES=3
MAX_DEPTH=3
while [[ $# -gt 0 ]]; do
  case "$1" in
    --refresh)       REFRESH=1; shift ;;
    --no-fetch)      NO_FETCH=1; shift ;;
    --genotypes)     GENOTYPES="${2:?--genotypes requires a value}"; shift 2 ;;
    --min-sequences) MIN_SEQUENCES="${2:?--min-sequences requires a value}"; shift 2 ;;
    --max-depth)     MAX_DEPTH="${2:?--max-depth requires a value}"; shift 2 ;;
    -h|--help)       sed -n '2,40p' "$0"; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; exit 1 ;;
  esac
done

echo "════════════════════════════════════════════════════"
echo " Build: hav-vp1-2b-lineages dataset (public-only)"
echo "════════════════════════════════════════════════════"
echo " Community: $COMMUNITY_NAME"
echo " Source:    $SRC_DATASET"
echo " Output:    $OUT_DATASET"
echo " Genotypes: $GENOTYPES"
echo ""

mkdir -p "$WORK_DIR" "$OUT_DATASET"

# ── 0. Fetch the community dataset (fresh, public source) ─────────────────────
echo "[0/6] Obtaining community dataset..."
if [[ "$NO_FETCH" -eq 1 ]]; then
  if [[ ! -f "$SRC_DATASET/tree.json" ]]; then
    echo "  ERROR: --no-fetch given but $SRC_DATASET/tree.json is missing." >&2
    exit 1
  fi
  echo "  --no-fetch: reusing existing $SRC_DATASET"
elif [[ "$REFRESH" -eq 1 || ! -f "$SRC_DATASET/tree.json" ]]; then
  echo "  Fetching '$COMMUNITY_NAME' → $SRC_DATASET"
  nextclade dataset get --name "$COMMUNITY_NAME" --output-dir "$SRC_DATASET"
else
  echo "  Reusing existing $SRC_DATASET (use --refresh to re-fetch)"
fi

# ── 1. Extract newick + nt_muts + metadata from the public tree.json ──────────
echo "[1/6] Extracting newick, nt_muts and metadata from public tree.json..."
NWK="$WORK_DIR/tree.nwk"
NT_MUTS="$WORK_DIR/nt_muts.json"
META="$WORK_DIR/metadata.tsv"
python3 "$SCRIPTS/extract_auspice_tree.py" \
    "$SRC_DATASET/tree.json" \
    "$NWK" \
    "$NT_MUTS" \
    "$META"

# ── 2. Designate lineages on the PUBLIC tree (per genotype, criteria-based) ────
echo "[2/6] Designating lineages on the public tree (per genotype)..."

# Back up an existing register, then start clean so the register is a faithful
# reflection of the current public tree (no stale, locally-derived rows).
if [[ -f "$LINEAGE_REGISTER" ]]; then
  _BACKUP="${LINEAGE_REGISTER%.tsv}_backup_$(date +%Y%m%d_%H%M%S).tsv"
  cp "$LINEAGE_REGISTER" "$_BACKUP"
  echo "  Existing register backed up → $(basename "$_BACKUP")"
fi
rm -f "$LINEAGE_REGISTER" "$LINEAGE_PROPOSALS"

CLADES_TSV="$WORK_DIR/lineage_clades.tsv"
printf 'clade\tgene\tsite\talt\n' > "$CLADES_TSV"

for GT in $GENOTYPES; do
  echo "  ── Genotype $GT ──"
  # NOTE: no --seed-from-outbreaks — the public tree has no outbreak labels.
  python3 "$SCRIPTS/designate_lineages.py" \
    "$NWK" \
    "$NT_MUTS" \
    "$META" \
    "$GT" \
    --register         "$LINEAGE_REGISTER" \
    --proposals        "$LINEAGE_PROPOSALS" \
    --filter-genotype  "$GT" \
    --min-sequences    "$MIN_SEQUENCES" \
    --max-depth        "$MAX_DEPTH" \
    --augur-clades-tsv "/tmp/lp_clades_${GT}.tsv"
  [[ -f "/tmp/lp_clades_${GT}.tsv" ]] && tail -n +2 "/tmp/lp_clades_${GT}.tsv" >> "$CLADES_TSV"
  rm -f "/tmp/lp_clades_${GT}.tsv"
done

N_CLADE_ROWS=$(tail -n +2 "$CLADES_TSV" | wc -l)
echo "  Lineage-defining mutation rows: $N_CLADE_ROWS"

# ── 3. Assign lineage_phylo labels via augur clades ───────────────────────────
echo "[3/6] Running augur clades (lineage_phylo)..."
LINEAGE_CLADES_JSON="$WORK_DIR/lineage_phylo_clades.json"
augur clades \
    --tree             "$NWK" \
    --mutations        "$NT_MUTS" \
    --clades           "$CLADES_TSV" \
    --membership-name  lineage_phylo \
    --label-name       lineage_phylo \
    --output-node-data "$LINEAGE_CLADES_JSON"

python3 - "$LINEAGE_CLADES_JSON" << 'PYEOF'
import json, sys, collections
with open(sys.argv[1]) as f:
    d = json.load(f)
counts = collections.Counter()
for node, attrs in d.get("nodes", {}).items():
    counts[attrs.get("lineage_phylo", "unassigned")] += 1
print(f"  lineage_phylo assignment summary ({len(counts)} distinct values):")
for val, n in sorted(counts.items(), key=lambda x: -x[1])[:20]:
    print(f"    {val}: {n} nodes")
PYEOF

# ── 4. Copy source dataset files and inject lineage_phylo into tree.json ───────
echo "[4/6] Patching tree.json with lineage_phylo annotations..."
# README.md is NOT copied from source — the lineages dataset has its own README.
for f in "$SRC_DATASET"/{CHANGELOG.md,example_sequences.fasta,genome_annotation.gff3,reference.fasta}; do
    cp "$f" "$OUT_DATASET/"
done

python3 "$SCRIPTS/patch_vp1b_tree.py" \
    "$SRC_DATASET/tree.json" \
    "$LINEAGE_CLADES_JSON" \
    "$OUT_DATASET/tree.json"

# ── 5. Write updated pathogen.json ────────────────────────────────────────────
echo "[5/6] Writing pathogen.json..."
python3 - "$SRC_DATASET/pathogen.json" "$OUT_DATASET/pathogen.json" << 'PYEOF'
import json, sys
from datetime import datetime, timezone

with open(sys.argv[1]) as f:
    p = json.load(f)

p["attributes"]["name"] = "Hepatitis A virus (with lineage designations)"
p["attributes"]["reference name"] = (
    "Hepatitis A virus, VP1-2B junction fragment of NCBI Reference Sequence (+ HAV lineage system)"
)
now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
p["version"]["updatedAt"] = now
p["version"]["tag"] = now.replace(":", "-")

with open(sys.argv[2], "w") as f:
    json.dump(p, f, indent=2)
print(f"  Wrote {sys.argv[2]}", file=__import__("sys").stderr)
PYEOF

# ── 6. Generate SILO lineage hierarchy YAML ───────────────────────────────────
echo "[6/6] Generating SILO lineage hierarchy (lineages.yaml)..."
LINEAGE_DEFS_DIR="$PROJ_DIR/data/lineage_definitions"
mkdir -p "$LINEAGE_DEFS_DIR"
SILO_YAML="$LINEAGE_DEFS_DIR/lineages.yaml"
python3 "$SCRIPTS/generate_silo_hierarchy.py" \
    "$LINEAGE_REGISTER" \
    "$SILO_YAML"
echo "  SILO hierarchy written to $SILO_YAML"

echo ""
echo "════════════════════════════════════════════════════"
echo " Done!  Dataset ready at:"
echo "   $OUT_DATASET"
echo ""
echo " Lineage register (derived from the public tree):"
echo "   $LINEAGE_REGISTER"
if [[ -f "$LINEAGE_PROPOSALS" ]]; then
echo " New lineage proposals to review:"
echo "   $LINEAGE_PROPOSALS"
echo "   python3 scripts/review_lineage_proposals.py $LINEAGE_PROPOSALS"
fi
echo ""
echo " SILO hierarchy (for LAPIS/Pathoplexus integration):"
echo "   $SILO_YAML"
echo ""
echo " Test with:"
echo "   nextclade run \\"
echo "     --input-dataset $OUT_DATASET \\"
echo "     --output-tsv /tmp/lineage_test.tsv \\"
echo "     <your_sequences.fasta>"
echo "════════════════════════════════════════════════════"
