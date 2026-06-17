#!/usr/bin/env python3
"""
Create a combined metadata TSV for Microreact by merging:
  - Reference metadata (metadata_corrected.tsv)
  - Batch query sequences from a Nextclade results TSV

Usage:
  python3 scripts/make_microreact_metadata.py \
    data/local_datasets/2026-04-10/metadata_corrected.tsv \
    data/Batch-1/local/all/nextclade.tsv \
    data/Batch-1/local/all/metadata_combined.tsv
"""
import csv
import sys

ref_file, nextclade_file, out_file = sys.argv[1], sys.argv[2], sys.argv[3]

# Read reference metadata
with open(ref_file, encoding="utf-8") as f:
    ref_rows = list(csv.DictReader(f, delimiter="\t"))

# Read batch nextclade results — only seqName + clade, skip alignment failures
batch_rows = []
with open(nextclade_file, encoding="utf-8") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        seq = row.get("seqName", "").strip()
        clade = row.get("clade", "").strip()
        if not seq:
            continue
        batch_rows.append({
            "id":       seq,
            "genotype": "batch",
            "lineage":  clade if clade else "unclassified",
            "date":     "NA",
            "country":  "NA",
        })

ref_ids = {r["id"] for r in ref_rows}
new_batch = [r for r in batch_rows if r["id"] not in ref_ids]

all_rows = ref_rows + new_batch
fieldnames = ["id", "genotype", "lineage", "date", "country"]

with open(out_file, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                            extrasaction="ignore")
    writer.writeheader()
    writer.writerows(all_rows)

print(f"Written {len(all_rows)} rows ({len(ref_rows)} reference + {len(new_batch)} batch) to {out_file}")
