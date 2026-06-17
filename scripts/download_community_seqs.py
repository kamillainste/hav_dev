#!/usr/bin/env python3
"""
One-time download of community sequences for the hav-vp1-2b-lineages dataset.
Extracts all leaf accessions from tree.json and fetches them from NCBI.

Usage: python3 scripts/download_community_seqs.py [project_root]
Output: data/nextclade_datasets/hav-vp1-2b-lineages/community_sequences.fasta
"""
import json, time, sys, os, re
from urllib import request, parse

PROJECT_ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
TREE_JSON  = os.path.join(PROJECT_ROOT, "data/nextclade_datasets/hav-vp1-2b-lineages/tree.json")
OUT_FASTA  = os.path.join(PROJECT_ROOT, "data/nextclade_datasets/hav-vp1-2b-lineages/community_sequences.fasta")

NCBI_EMAIL  = "hav.analysis@fhi.no"
BATCH_SIZE  = 200
SLEEP_SEC   = 0.4   # stay comfortably below 3 req/sec without API key
EUTILS_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# ── Extract leaf names ────────────────────────────────────────────────────────
print(f"Loading {TREE_JSON}...")
with open(TREE_JSON) as f:
    tree = json.load(f)

leaves = []
def find_leaves(node):
    if not node.get("children"):
        name = node.get("name", "")
        if name:
            leaves.append(name)
    else:
        for c in node["children"]:
            find_leaves(c)

find_leaves(tree["tree"])
print(f"Total leaves: {len(leaves)}")

# Filter to NCBI-like accessions (start with letters followed by digits)
accessions = [l for l in leaves if re.match(r'^[A-Za-z]{1,2}\d', l) and not l.startswith("NODE_")]
skipped    = [l for l in leaves if l not in accessions]
print(f"NCBI accessions to download: {len(accessions)}")
if skipped:
    print(f"Skipped (non-accession names): {skipped[:5]}")

# ── Download in batches ───────────────────────────────────────────────────────
def fetch_batch(ids, attempt=1):
    data = parse.urlencode({
        "db":      "nucleotide",
        "id":      ",".join(ids),
        "rettype": "fasta",
        "retmode": "text",
        "email":   NCBI_EMAIL,
        "tool":    "hav_pipeline",
    }).encode()
    try:
        req = request.Request(EUTILS_URL, data=data, method="POST")
        with request.urlopen(req, timeout=60) as resp:
            return resp.read().decode("utf-8", errors="replace")
    except Exception as e:
        if attempt < 4:
            wait = 2 ** attempt
            print(f"  Retry {attempt} after {wait}s (error: {e})", file=sys.stderr)
            time.sleep(wait)
            return fetch_batch(ids, attempt + 1)
        print(f"  ERROR: giving up on batch after {attempt} attempts: {e}", file=sys.stderr)
        return ""

n_batches = (len(accessions) + BATCH_SIZE - 1) // BATCH_SIZE
total_downloaded = 0

with open(OUT_FASTA, "w") as out:
    for i in range(0, len(accessions), BATCH_SIZE):
        batch = accessions[i : i + BATCH_SIZE]
        batch_num = i // BATCH_SIZE + 1
        print(f"  Batch {batch_num}/{n_batches}: downloading {len(batch)} sequences...")
        fasta = fetch_batch(batch)
        if fasta.strip():
            out.write(fasta)
            n = fasta.count(">")
            total_downloaded += n
            print(f"    -> {n} sequences received")
        time.sleep(SLEEP_SEC)

print(f"\nDone. Total sequences downloaded: {total_downloaded}/{len(accessions)}")
print(f"Output: {OUT_FASTA}")
