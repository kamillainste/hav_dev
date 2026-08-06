#!/usr/bin/env bash
# scripts/build_per_seq_trees.sh
#
# Builds a per-query-sequence IQ-TREE phylogenetic tree containing the query
# and its N nearest neighbors drawn from:
#   - Top BLAST hits (ranked by bitscore)
#   - Nearest-node sequences from all NextClade runs (nearestNodeName from NDJSON)
#
# Called by run_all_analyses.sh with:
#   bash scripts/build_per_seq_trees.sh <batch_dir> <dataset_date> <n_neighbors> <out_base> <dataset_dir> <batch_fa>
#
# Arguments:
#   <batch_dir>    Batch directory path
#   <dataset_date> Dataset date (e.g., 2026-04-10)
#   <n_neighbors>  Number of nearest neighbors to include (default: 30)
#   <out_base>     Output base directory where results go
#   <dataset_dir>  Path to local_datasets/<dataset_date>
#   <batch_fa>     Path to batch FASTA file
#
# Writes per sequence to <out_base>/trees/<seqName>/:
#   neighbors.txt         — list of reference IDs used
#   combined.fa           — query + neighbor sequences (unaligned)
#   aligned.fa            — MAFFT alignment
#   aligned_trimmed.fa    — MAFFT alignment trimmed to the junction region
#   tree.treefile         — IQ-TREE best ML tree with bootstrap values
#   tree.iqtree           — IQ-TREE run summary
#   tree.log              — IQ-TREE log
#   tree.contree          — consensus tree
#
# Prerequisites:
#   - run_all_analyses.sh must have been run first (creates BLAST and NextClade results)
#   - HAVDEV conda: mafft, iqtree (iqtree3), seqkit, python3, jq

set -euo pipefail

BATCH_DIR="${1:?Usage: bash scripts/build_per_seq_trees.sh <batch_dir> [dataset_date] [n_neighbors] [output_base]}"
DATASET_DATE="${2:-2026-04-10}"
N_NEIGHBORS="${3:-30}"
OUT_BASE="${4:?out_base is required}"
DATASET_DIR="${5:?dataset_dir is required}"
BATCH_FA="${6:?batch_fasta is required}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DEDUP_FA="$DATASET_DIR/blast_db/input_dedup.fa"
BLAST_DB="$DATASET_DIR/blast_db/hav"

BATCH_NAME=$(basename "$BATCH_DIR")
TREES_DIR="$OUT_BASE/trees"

BLAST_TSV="$OUT_BASE/blast_results.tsv"
NC_LINEAGES_NDJSON="$OUT_BASE/lineages/nextclade.ndjson"

LINEAGES_DATASET="$PROJECT_DIR/data/nextclade_datasets/hav-vp1-2b-lineages"
LINEAGES_TREE_JSON="$LINEAGES_DATASET/tree.json"
COMMUNITY_SEQS_FA="$LINEAGES_DATASET/community_sequences.fasta"

THREADS=4

cd "$PROJECT_DIR"


# ── Validation ────────────────────────────────────────────────────────────────
for f in "$BATCH_FA" "$BLAST_TSV" "$DEDUP_FA"; do
  [[ -f "$f" ]] || { echo "ERROR: Required file not found: $f" >&2; exit 1; }
done

if ! command -v mafft &>/dev/null; then
  echo "ERROR: mafft not found. Ensure HAVDEV conda is active." >&2; exit 1
fi
# IQ-TREE v3 installs as 'iqtree' in the HAVDEV env
IQTREE_BIN="iqtree3"
if ! command -v "$IQTREE_BIN" &>/dev/null; then
  IQTREE_BIN="iqtree"
  command -v "$IQTREE_BIN" &>/dev/null || { echo "ERROR: iqtree not found." >&2; exit 1; }
fi

mkdir -p "$TREES_DIR"

# ── Python helper: collect neighbor IDs for a query ──────────────────────────
COLLECT_PY="$(mktemp /tmp/collect_neighbors.XXXXXX.py)"
cat > "$COLLECT_PY" << 'PYEOF'
"""
Collect up to N neighbor sequence IDs for a given query from:
  1. BLAST TSV  — top hits ranked by bitscore
  2. Community tree (hav-vp1-2b-lineages/tree.json) — leaf descendants of the
     Nextclade nearest node, ranked by accumulated branch mutations

Strategy: reserve up to N*2//3 slots for BLAST, fill remaining from community.

Usage: python3 script.py <query_id> <n> <blast_tsv> <nc_lineages_ndjson> <unused> [tree_json]
"""
import sys, json

query_id    = sys.argv[1]
n           = int(sys.argv[2])
blast_tsv   = sys.argv[3]
ndjson_path = sys.argv[4]   # nextclade.ndjson for nearestNodeName
# sys.argv[5] unused (legacy NONE)
tree_json   = sys.argv[6] if len(sys.argv) > 6 else "NONE"

# BLAST truncates FASTA headers at the first whitespace, so qseqid is always
# the first word of the full sequence name.  Match on that prefix.
blast_query_id = query_id.split()[0]

# ── 1. BLAST hits ─────────────────────────────────────────────────────────────
blast_hits = []   # [(id, bitscore)]
blast_seen = set()
try:
    with open(blast_tsv) as f:
        next(f)  # skip header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12 or parts[0] != blast_query_id:
                continue
            subject = parts[1]
            try:
                bitscore = float(parts[11])
            except ValueError:
                bitscore = 0.0
            if subject and not subject.startswith("NODE_") and subject not in blast_seen:
                blast_hits.append((subject, bitscore))
                blast_seen.add(subject)
except FileNotFoundError:
    pass

blast_hits.sort(key=lambda x: -x[1])
# Reserve up to 2/3 for BLAST; the rest for community.
# If community yields nothing, fall back to all-BLAST up to n.
n_blast_max    = n * 2 // 3
selected_blast = [h[0] for h in blast_hits[:n_blast_max]]
blast_fallback = [h[0] for h in blast_hits]  # full list for fallback

# ── 2. Community tree neighbors ────────────────────────────────────────────────
community_hits = []   # [(leaf_name, dist_from_nearest_node)]

if tree_json != "NONE" and ndjson_path != "NONE":
    # Get nearestNodeName from Nextclade NDJSON
    nearest_node_name = None
    try:
        with open(ndjson_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                if rec.get("seqName") == query_id:
                    nearest_node_name = rec.get("nearestNodeName")
                    break
    except FileNotFoundError:
        pass

    if nearest_node_name:
        try:
            with open(tree_json) as f:
                tree = json.load(f)

            def count_muts(node):
                return len(node.get("branch_attrs", {}).get("mutations", {}).get("nuc", []))

            found_node   = [None]
            found_parent = [None]

            def find_node(node, parent):
                if node.get("name") == nearest_node_name:
                    found_node[0]   = node
                    found_parent[0] = parent
                    return True
                for c in node.get("children", []):
                    if find_node(c, node):
                        return True
                return False

            find_node(tree["tree"], None)

            def collect_leaves(node, dist, results, exclude=None):
                if not node.get("children"):
                    name = node.get("name", "")
                    if name and not name.startswith("NODE_"):
                        if exclude is None or name not in exclude:
                            results.append((name, dist))
                else:
                    for c in node["children"]:
                        collect_leaves(c, dist + count_muts(c), results, exclude)

            if found_node[0]:
                collect_leaves(found_node[0], 0, community_hits)

                # If few descendants, expand to sibling subtrees from parent
                n_needed = n - len(selected_blast)
                if len(community_hits) < n_needed and found_parent[0]:
                    existing = {h[0] for h in community_hits}
                    for sibling in found_parent[0].get("children", []):
                        if sibling.get("name") != nearest_node_name:
                            collect_leaves(sibling, 0, community_hits, exclude=existing)
                            existing = {h[0] for h in community_hits}
        except Exception:
            pass

community_hits.sort(key=lambda x: x[1])
n_community_slots = n - len(selected_blast)
blast_set         = set(selected_blast)
selected_community = [
    h[0] for h in community_hits
    if h[0] not in blast_set
][:n_community_slots]

# If no community sequences found, fill remaining slots from BLAST
if not selected_community:
    selected_blast = blast_fallback[:n]
    selected_community = []

# ── Output ─────────────────────────────────────────────────────────────────────
for sid in selected_blast + selected_community:
    print(sid)
PYEOF

# ── Process each query sequence ───────────────────────────────────────────────
# Extract all query IDs from the batch FASTA
QUERY_IDS=()                          # Initialiserer tom array
while IFS= read -r line; do           # Les hver linje fra FASTA-filen
  [[ "$line" == ">"* ]] && {
    # Fjern ">" and strip ALL whitespace including newlines/carriage returns
    id="${line#>}"
    # Strip carriage returns first, then whitespace
    id="${id%$'\r'}"
    id=$(echo "$id" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    QUERY_IDS+=("$id")
  }
done < "$BATCH_FA"                    # Les fra BATCH_FA-filen

echo "════════════════════════════════════════════════════════════════"
echo "Building per-sequence trees: $BATCH_NAME"
echo "Sequences : ${#QUERY_IDS[@]}"
echo "Neighbors : up to $N_NEIGHBORS per sequence"
echo "Dataset   : $DATASET_DATE"
echo "════════════════════════════════════════════════════════════════"
echo ""

for QUERY in "${QUERY_IDS[@]}"; do
  # Normalize: trim whitespace and take only first token (before space)
  # This ensures folder names match what NextClade produces
  QUERY_NORMALIZED=$(echo "$QUERY" | awk '{print $1}' | xargs)
  SEQ_DIR="$TREES_DIR/$QUERY_NORMALIZED"
  mkdir -p "$SEQ_DIR"

  echo "── $QUERY ─────────────────────────────────────────────────────"

  # Collect neighbor IDs
  NC_LINEAGES_NDJSON_ARG="$NC_LINEAGES_NDJSON"
  [[ -f "$NC_LINEAGES_NDJSON_ARG" ]] || NC_LINEAGES_NDJSON_ARG="NONE"

  LINEAGES_TREE_JSON_ARG="$LINEAGES_TREE_JSON"
  [[ -f "$LINEAGES_TREE_JSON_ARG" ]] || LINEAGES_TREE_JSON_ARG="NONE"
#Disse filene er valgfrie. 
#Hvis NextClade-kjøringen ikke fullfører eller treet mangler, skal skriptet likevel kjøre 
#— bare uten den ekstra informasjonen. 
#Python-scriptet collect_neighbors.py sjekker deretter om verdien er "NONE" og 
#hopper over den delen av analysen.
 
  python3 "$COLLECT_PY" \
    "$QUERY_NORMALIZED" "$N_NEIGHBORS" \
    "$BLAST_TSV" \
    "$NC_LINEAGES_NDJSON_ARG" \
    "NONE" \
    "$LINEAGES_TREE_JSON_ARG" \
    > "$SEQ_DIR/neighbors.txt"

  N_IDS=$(wc -l < "$SEQ_DIR/neighbors.txt")
  echo "  Neighbors found: $N_IDS"

  if [[ "$N_IDS" -eq 0 ]]; then
    echo "  WARNING: No neighbors found — skipping $QUERY"
    continue
  fi

  # Extract query sequence by exact FASTA header text from batch FASTA.  
  # BLAST truncates headers at whitespace, so we must match on the full header line, not just the query ID.
  # seqkit name matching truncates at whitespace for many FASTA styles.
  # Strip carriage returns to handle Windows line endings.
  QUERY_FA="$SEQ_DIR/query.fa"
  awk -v q="$QUERY" '
    BEGIN {keep=0}
    /^>/ {
      header = substr($0,2)
      # Strip trailing carriage returns and whitespace to handle Windows line endings
      gsub(/\r$/, "", header)
      gsub(/[[:space:]]*$/, "", header)
      keep = (header == q)
    }
    keep {print}
  ' "$BATCH_FA" > "$QUERY_FA"

  if [[ ! -s "$QUERY_FA" ]]; then
    echo "  WARNING: Could not extract query sequence '$QUERY' from FASTA — skipping"
    continue
  fi

  # Extract neighbor sequences from local reference FASTA
  NEIGHBORS_FA="$SEQ_DIR/reference_seqs.fa"
  seqkit grep -f "$SEQ_DIR/neighbors.txt" "$DEDUP_FA" > "$NEIGHBORS_FA"

  # Also extract any community accessions from the downloaded community FASTA
  if [[ -f "$COMMUNITY_SEQS_FA" ]]; then
    seqkit grep -f "$SEQ_DIR/neighbors.txt" "$COMMUNITY_SEQS_FA" >> "$NEIGHBORS_FA"
  fi

  N_REF=$(grep -c "^>" "$NEIGHBORS_FA" 2>/dev/null || echo 0)
  echo "  Reference sequences extracted: $N_REF"

  if [[ "$N_REF" -eq 0 ]]; then
    echo "  WARNING: No reference sequences could be extracted — skipping $QUERY"
    continue
  fi

  # Combine query + references
  COMBINED_FA="$SEQ_DIR/combined.fa"
  cat "$QUERY_FA" "$NEIGHBORS_FA" > "$COMBINED_FA"
  TOTAL=$(grep -c "^>" "$COMBINED_FA")
  echo "  Total sequences for tree: $TOTAL"

  if [[ "$TOTAL" -lt 3 ]]; then
    echo "  WARNING: Need at least 3 sequences for a tree — skipping $QUERY"
    continue
  fi

  # Align with MAFFT
  ALIGNED_FA="$SEQ_DIR/aligned.fa"
  echo "  Aligning with MAFFT..."
  mafft --auto --quiet --thread "$THREADS" "$COMBINED_FA" > "$ALIGNED_FA"

  # Remove all-gap columns (trim alignment)
  TRIMMED_FA="$SEQ_DIR/aligned_trimmed.fa"
  python3 - "$ALIGNED_FA" "$TRIMMED_FA" << 'TRIM_PY'
import sys
seqs = {}
order = []
with open(sys.argv[1]) as f:
    cur_id = None
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            cur_id = line[1:]
            seqs[cur_id] = []
            order.append(cur_id)
        elif cur_id:
            seqs[cur_id].append(line)
seqs = {k: "".join(v) for k, v in seqs.items()}
# find non-all-gap columns
lens = set(len(v) for v in seqs.values())
if len(lens) != 1:
    # sequences differ in length — write as-is
    with open(sys.argv[2], "w") as out:
        for sid in order:
            out.write(f">{sid}\n{seqs[sid]}\n")
else:
    L = next(iter(lens))
    keep = [i for i in range(L) if any(seqs[s][i] not in "-N" for s in order)]
    with open(sys.argv[2], "w") as out:
        for sid in order:
            trimmed = "".join(seqs[sid][i] for i in keep)
            out.write(f">{sid}\n{trimmed}\n")
TRIM_PY

  ALN_LEN=$(python3 -c "
import sys
with open('$TRIMMED_FA') as f:
    for line in f:
        if not line.startswith('>'):
            print(len(line.strip()))
            break
" 2>/dev/null || echo 0)
  echo "  Alignment length: $ALN_LEN bp"

  # Build ML tree with IQ-TREE.
  # UFBoot requires at least 4 taxa, so fall back to no-bootstrap for n=3.
  TREE_PREFIX="$SEQ_DIR/tree"
  if [[ "$TOTAL" -lt 4 ]]; then
    echo "  Running IQ-TREE (GTR+G, no bootstrap; <4 taxa)..."
    IQTREE_ARGS=(
      -s   "$TRIMMED_FA"
      -m   GTR+G
      -T   AUTO
      --prefix "$TREE_PREFIX"
      -redo
    )
  else
    echo "  Running IQ-TREE (GTR+G, 1000 UFBoot)..."
    IQTREE_ARGS=(
      -s   "$TRIMMED_FA"
      -m   GTR+G
      -B   1000
      -T   AUTO
      --prefix "$TREE_PREFIX"
      -redo
    )
  fi

  # Remove old IQ-TREE output to allow --redo-free run
  rm -f "$TREE_PREFIX".* 2>/dev/null || true

  # Run IQ-TREE and capture errors
  if "$IQTREE_BIN" "${IQTREE_ARGS[@]}" 2>&1 | tee "$SEQ_DIR/iqtree.log"; then
    IQTREE_STATUS=0
  else
    IQTREE_STATUS=$?
  fi

  if [[ -f "$TREE_PREFIX.treefile" ]]; then
    echo "  ✓ Tree written: $TREE_PREFIX.treefile"
  else
    echo "  ✗ WARNING: IQ-TREE did not produce a treefile for $QUERY"
    if [[ $IQTREE_STATUS -ne 0 ]]; then
      echo "    IQ-TREE exit code: $IQTREE_STATUS"
    fi
    if [[ -f "$SEQ_DIR/iqtree.log" ]]; then
      echo "    Last 10 lines of IQ-TREE log:"
      tail -10 "$SEQ_DIR/iqtree.log" | sed 's/^/      /'
    fi
  fi
  echo ""
done

rm -f "$COLLECT_PY"

echo "════════════════════════════════════════════════════════════════"
echo "Per-sequence trees complete."
echo "  Output: $TREES_DIR/<seqName>/tree.treefile"
echo ""
echo "  Next step — generate report:"
echo "    \"\$HOME/.conda/R_shared/bin/Rscript\" -e \"rmarkdown::render('scripts/batch_report.Rmd',"
echo "      params=list(batch_dir='$BATCH_DIR', dataset_date='$DATASET_DATE'))\""
echo "════════════════════════════════════════════════════════════════"
