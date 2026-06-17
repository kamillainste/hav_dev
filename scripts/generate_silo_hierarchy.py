#!/usr/bin/env python3
"""
generate_silo_hierarchy.py

Generates a SILO/LAPIS lineage hierarchy file (lineages.yaml) from the
canonical HAV lineage register (lineage_designations.tsv).

The output YAML encodes the parent-child DAG of all active lineages in the
SILO lineage_definitions format:
  https://github.com/GenSpectrum/LAPIS-SILO/blob/main/documentation/lineage_definitions.md

This enables monophyletic (wildcard) search in LAPIS/SILO/Pathoplexus:
  e.g. "all sequences in lineage IA.1 and all descendants"

Format:
  <lineage_label>:
    aliases: []
    parents:
    - <parent_label>

Special root nodes:
  None      -- ultimate root; unassigned sequences belong here
  unassigned -- explicit handle for unclassified sequences

Ordering: topological (parents before children), then alphabetical within depth.

Usage:
  python3 generate_silo_hierarchy.py <lineage_designations.tsv> <out_lineages.yaml>
"""

import csv
import sys
from collections import OrderedDict
from pathlib import Path


def infer_parent(lineage: str) -> str:
    """
    Infer parent lineage from hierarchical Pango-style name.
    Examples:
      IA.1.2.1 -> IA.1.2
      IA.1     -> IA
      IA       -> None   (genotype root, no dot)
    """
    if "." not in lineage:
        # Genotype root (IA, IB, IIA, IIIA)
        return "None"
    return lineage.rsplit(".", 1)[0]


def depth(lineage: str) -> int:
    """Number of levels below root. 'IA'=0, 'IA.1'=1, 'IA.1.2'=2 ..."""
    if "." not in lineage:
        return 0
    return lineage.count(".")


def build_hierarchy(tsv_path: str) -> OrderedDict:
    """
    Read lineage_designations.tsv and build the SILO hierarchy dict.
    Only non-deprecated rows are included.
    Genotype roots are inferred from the 'genotype' column.
    """
    lineages: dict[str, dict] = {}
    genotypes: set[str] = set()

    with open(tsv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("deprecated", "").strip().upper() == "TRUE":
                continue
            name = row["lineage"].strip()
            geno = row["genotype"].strip()
            genotypes.add(geno)
            lineages[name] = {"aliases": [], "parents": [infer_parent(name)]}

    # Add genotype root nodes (IA, IB, IIA, IIIA) if not already present
    for geno in genotypes:
        if geno not in lineages:
            lineages[geno] = {"aliases": [], "parents": ["None"]}

    # Special root nodes
    result = OrderedDict()
    result["None"] = {"aliases": [], "parents": []}
    result["unassigned"] = {"aliases": [], "parents": []}

    # Sort remaining by depth then alphabetically
    def sort_key(name: str):
        return (depth(name), name)

    for name in sorted(lineages.keys(), key=sort_key):
        result[name] = lineages[name]

    return result


def write_yaml(hierarchy: OrderedDict, out_path: str) -> None:
    """Write the hierarchy in SILO lineages.yaml format (hand-crafted, no pyyaml dependency)."""
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        for name, attrs in hierarchy.items():
            f.write(f"{name}:\n")
            f.write(f"  aliases: {attrs['aliases']}\n")
            parents = attrs["parents"]
            if parents:
                f.write(f"  parents:\n")
                for p in parents:
                    f.write(f"  - {p}\n")
            else:
                f.write(f"  parents: []\n")


def main() -> None:
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <lineage_designations.tsv> <out_lineages.yaml>",
              file=sys.stderr)
        sys.exit(1)

    tsv_path = sys.argv[1]
    out_path = sys.argv[2]

    hierarchy = build_hierarchy(tsv_path)
    write_yaml(hierarchy, out_path)

    n_lineages = sum(1 for k in hierarchy if k not in ("None", "unassigned"))
    print(f"  Wrote {n_lineages} lineage entries + None/unassigned roots to {out_path}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
