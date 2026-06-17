"""Inject lineage_phylo annotations from augur clades output into an Auspice v2 tree.json.

This script:
  1. Reads the original hav-vp1-2b tree.json (source).
  2. Reads lineage_phylo_clades.json from augur clades.
  3. Injects lineage_phylo into node_attrs for every node.
  4. Adds the required meta declarations for NextClade (extensions.nextclade.clade_node_attrs)
     and Auspice (meta.clade_node_attrs) so the new column appears in the UI and TSV output.
  5. Writes the patched tree.json to the output path.

Usage:
    python3 patch_vp1b_tree.py <src_tree.json> <lineage_phylo_clades.json> <out_tree.json>
"""
import json
import sys
from pathlib import Path

# Map of broad clade_membership values to the lineage prefix they are allowed to carry.
# Any lineage_phylo value that does not start with the expected prefix is a
# cross-genotype convergent misassignment.
_GENO_PREFIX = {
    "IA":   "IA",
    "IB":   "IB",
    "IIA":  "IIA",
    "IIIA": "IIIA",
}


def _inject(node, lp_nodes, parent_final_lp="unassigned", parent_augur_lp="unassigned"):
    """Recursively inject lineage_phylo into every node's node_attrs.

    Two-stage correction is applied after taking the augur-assigned label:

    1. Cascade clearance: augur propagates clade labels top-down, so if the
       "clade founder" node receives a wrong label (e.g. IA.1.2.1 on an
       unassigned-clade internal node due to a convergent T21C), all descendants
       that merely *inherited* the same label would also be wrong.  We detect
       inheritance by comparing a node's augur label to its parent's augur label.
       When they match (inherited rather than newly founded), we use the parent's
       *final* (post-correction) label instead of the augur label, propagating
       any correction made higher up.

    2. Genotype-lineage consistency: if the (possibly cascade-corrected) label
       still has a genotype prefix that does not match the node's broad
       clade_membership (e.g. IIIA.1.2 on an IA-broad node due to convergent
       T27C, or IA.1.2.1 on an IC/unassigned node), it is replaced with the
       parent's final label.  Nodes with no recognised genotype (IC, IIB,
       unassigned) are set to "unassigned".
    """
    name = node.get("name", "")
    attrs = node.setdefault("node_attrs", {})

    broad = attrs.get("clade_membership", {}).get("value", "")
    augur_lp = lp_nodes.get(name, {}).get("lineage_phylo", "unassigned")

    lp = augur_lp

    # Stage 1: cascade clearance
    # If this node's augur label == parent's augur label, the label was inherited
    # (not founded here).  Use the parent's final label so that any correction
    # made to an ancestor propagates downward.
    if lp != "unassigned" and lp == parent_augur_lp:
        lp = parent_final_lp

    # Stage 2: genotype-lineage consistency
    expected_prefix = _GENO_PREFIX.get(broad)  # None for IC, IIB, unassigned, etc.
    if lp != "unassigned":
        if expected_prefix is None:
            # Node has no recognised genotype (IC, IIB, unassigned …) → clear
            lp = "unassigned"
        elif not lp.startswith(expected_prefix):
            # Cross-genotype misassignment → fall back to parent's final label
            lp = parent_final_lp

    attrs["lineage_phylo"] = {"value": lp}
    for child in node.get("children", []):
        # Pass our own augur_lp as the child's parent_augur_lp so that
        # children can detect whether they inherited from us.
        _inject(child, lp_nodes, lp, augur_lp)


def main():
    if len(sys.argv) != 4:
        print("Usage: patch_vp1b_tree.py <src_tree.json> <lineage_phylo_clades.json> <out_tree.json>",
              file=sys.stderr)
        sys.exit(1)

    src_path, clades_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

    with open(src_path) as f:
        tree = json.load(f)

    with open(clades_path) as f:
        clades = json.load(f)

    lp_nodes = clades.get("nodes", {})
    print(f"  lineage_phylo assignments loaded: {len(lp_nodes)} nodes", file=sys.stderr)

    # --- Inject lineage_phylo into all tree nodes ---
    _inject(tree["tree"], lp_nodes)

    # --- Update meta for NextClade and Auspice ---
    meta = tree.setdefault("meta", {})

    # NextClade requires extensions.nextclade.clade_node_attrs for extra TSV columns
    ext = meta.setdefault("extensions", {})
    nc_ext = ext.setdefault("nextclade", {})
    cna = nc_ext.setdefault("clade_node_attrs", [])
    existing_nc = [c.get("name") for c in cna]
    if "lineage_phylo" not in existing_nc:
        cna.insert(0, {"name": "lineage_phylo", "displayName": "Lineage (phylo)"})

    # Auspice branch labels and legend
    top_cna = meta.setdefault("clade_node_attrs", [])
    existing_top = [c.get("name") for c in top_cna]
    if "lineage_phylo" not in existing_top:
        top_cna.insert(0, {
            "name": "lineage_phylo",
            "displayName": "Lineage (phylo)",
            "hideInLegend": True
        })
    if "clade_membership" not in existing_top:
        top_cna.append({"name": "clade_membership", "displayName": "Genotype"})

    # Display defaults
    dd = meta.setdefault("display_defaults", {})
    dd.setdefault("branch_label", "lineage_phylo")
    dd.setdefault("color_by", "clade_membership")
    dd.setdefault("distance_measure", "div")

    # Add lineage_phylo coloring if not present
    colorings = meta.setdefault("colorings", [])
    coloring_keys = [c.get("key") for c in colorings]
    if "lineage_phylo" not in coloring_keys:
        colorings.insert(0, {
            "key": "lineage_phylo",
            "title": "Lineage (phylo)",
            "type": "categorical"
        })

    # --- Write output (minified) ---
    with open(out_path, "w") as f:
        json.dump(tree, f, separators=(",", ":"))

    # Summary stats
    assigned = sum(
        1 for v in lp_nodes.values()
        if v.get("lineage_phylo", "unassigned") != "unassigned"
    )
    print(f"  Injected lineage_phylo: {assigned}/{len(lp_nodes)} nodes assigned (non-unassigned)",
          file=sys.stderr)
    print(f"  Wrote patched tree.json: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
