"""Extract newick tree, augur-format nt_muts.json, and (optionally) a metadata
TSV from an Auspice v2 tree.json.

This reverses `augur export v2`, pulling out:
  1. tree.nwk  — newick topology with branch lengths (uses `branch_attrs.div` or falls back
                  to node_attrs.div differences)
  2. nt_muts.json — augur ancestral node-data format:
                    {"nodes": {"NODE_X": {"muts": ["C15T", ...]}, ...}}
  3. metadata.tsv (optional) — one row per tip with columns
                    id, genotype, lineage, date, country
                    derived purely from tip node_attrs (clade_membership → genotype,
                    geo_loc_name → country, date → date).  The 'lineage' (outbreak)
                    column is always empty because the public community tree carries
                    no outbreak labels — this is what makes designation public-only.

Usage:
    python3 extract_auspice_tree.py <tree.json> <out_tree.nwk> <out_nt_muts.json> [out_metadata.tsv]
"""
import csv
import json
import sys
import re
from pathlib import Path


def _div(node):
    """Return the divergence value for a node (float or 0)."""
    nd = node.get("node_attrs", {}).get("div")
    if isinstance(nd, dict):
        return float(nd.get("value", 0))
    if nd is not None:
        return float(nd)
    return 0.0


def _newick(node, parent_div=0.0):
    """Recursively build a newick string from an Auspice v2 tree node."""
    name = node.get("name", "")
    # Escape name characters that break newick
    safe_name = re.sub(r"[(),;:\s]", "_", name)
    div = _div(node)
    bl = max(div - parent_div, 0.0)

    children = node.get("children", [])
    if children:
        child_parts = ",".join(_newick(c, div) for c in children)
        return f"({child_parts}){safe_name}:{bl:.8f}"
    else:
        return f"{safe_name}:{bl:.8f}"


def _extract_muts(node, nodes_dict):
    """Walk tree, collect branch mutations into augur nt_muts format."""
    name = node.get("name", "")
    nuc_muts = (
        node.get("branch_attrs", {})
            .get("mutations", {})
            .get("nuc", [])
    )
    nodes_dict[name] = {"muts": list(nuc_muts)}
    for child in node.get("children", []):
        _extract_muts(child, nodes_dict)


def _attr(node, key):
    """Return a node_attrs value, unwrapping the {'value': ...} envelope."""
    v = node.get("node_attrs", {}).get(key)
    if isinstance(v, dict):
        return v.get("value", "") or ""
    return v if v is not None else ""


def _extract_metadata(node, rows):
    """Walk tree, collect one metadata row per tip.

    genotype  ← clade_membership   (IA, IB, IIA, IIIA, …)
    country   ← geo_loc_name       (text before the first ':' — 'Japan: …' → 'Japan')
    date      ← date               (may be partial, e.g. '1985-XX-XX')
    lineage   ← always ''          (no outbreak labels in the public community tree)
    """
    if not node.get("children"):
        name = node.get("name", "")
        if name:
            geo = _attr(node, "geo_loc_name")
            country = geo.split(":")[0].strip() if geo else ""
            rows.append({
                "id":       name,
                "genotype": _attr(node, "clade_membership"),
                "lineage":  "",
                "date":     _attr(node, "date"),
                "country":  country,
            })
    else:
        for child in node.get("children", []):
            _extract_metadata(child, rows)


def main():
    if len(sys.argv) not in (4, 5):
        print("Usage: extract_auspice_tree.py <tree.json> <out.nwk> <out_nt_muts.json> "
              "[out_metadata.tsv]", file=sys.stderr)
        sys.exit(1)

    tree_json_path, nwk_path, nt_muts_path = sys.argv[1], sys.argv[2], sys.argv[3]
    metadata_path = sys.argv[4] if len(sys.argv) == 5 else None

    with open(tree_json_path) as f:
        tree_json = json.load(f)

    root = tree_json["tree"]

    # --- Write newick ---
    nwk = _newick(root) + ";"
    Path(nwk_path).write_text(nwk)
    print(f"  Wrote newick: {nwk_path}", file=sys.stderr)

    # --- Write nt_muts.json ---
    nodes = {}
    _extract_muts(root, nodes)
    nt_muts = {"nodes": nodes}
    with open(nt_muts_path, "w") as f:
        json.dump(nt_muts, f)
    print(f"  Wrote nt_muts.json: {nt_muts_path} ({len(nodes)} nodes)", file=sys.stderr)

    # --- Write metadata TSV (optional) ---
    if metadata_path:
        rows = []
        _extract_metadata(root, rows)
        with open(metadata_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f, fieldnames=["id", "genotype", "lineage", "date", "country"],
                delimiter="\t", lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        n_gt = sum(1 for r in rows if r["genotype"])
        print(f"  Wrote metadata: {metadata_path} "
              f"({len(rows)} tips, {n_gt} with genotype)", file=sys.stderr)


if __name__ == "__main__":
    main()
