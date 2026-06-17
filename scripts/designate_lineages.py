#!/usr/bin/env python3
"""
designate_lineages.py — Propose phylogenetic lineage designations for HAV.

Traverses a reference tree (e.g. the public hav-vp1-2b tree extracted by
extract_auspice_tree.py) and identifies monophyletic groups that qualify as
named lineages under a Pango-inspired hierarchical naming scheme:

    IA.1        — first top-level lineage within genotype IA
    IA.1.1      — first sub-lineage of IA.1
    IA.2        — second top-level lineage (separate clade from IA.1)

Qualification criteria (configurable via CLI or DEFAULT_* constants):
    • ≥ MIN_SEQUENCES leaf sequences in the clade
    • ≥ MIN_DEFINING_MUTS branch mutations at this node (from augur ancestral)
    • bootstrap support ≥ MIN_SUPPORT (skipped when tree has no support values)

The script compares proposed designations against an existing register
(lineage_designations.tsv) by first matching on defining_nt_muts (stable across
minor topology changes), then falling back to founder_node name. Known lineages
keep their names; new candidates are written to lineage_proposals.tsv for review.

The --seed-from-outbreaks flag bypasses qualification criteria for existing
outbreak clusters: the MRCA of each outbreak group is designated automatically.
It requires outbreak labels in the metadata 'lineage' column and is therefore
only usable on local/private trees that carry such labels — NOT on the public
community tree, where designation is purely criteria-based.

Usage (public-only build, via build_vp1b_lineage_dataset.sh — inputs extracted
from the community hav-vp1-2b tree.json):
    python3 scripts/designate_lineages.py \\
        <work>/tree.nwk \\
        <work>/nt_muts.json \\
        <work>/metadata.tsv \\
        IA \\
        --register        data/lineage_designations.tsv \\
        --proposals       data/lineage_proposals.tsv \\
        --filter-genotype IA
"""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Optional

from Bio import Phylo

# ── Default thresholds (change here to adjust project-wide defaults) ──────────
DEFAULT_MIN_SEQUENCES     = 3
DEFAULT_MIN_SUPPORT       = 70.0   # bootstrap 0–100; set 0 to disable check
DEFAULT_MIN_DEFINING_MUTS = 1
DEFAULT_MAX_DEPTH         = 3      # e.g. IA → IA.1 → IA.1.1 (depth 3)

REGISTER_COLS = [
    "lineage", "genotype", "founder_node", "defining_nt_muts",
    "num_sequences", "earliest_date", "latest_date", "countries",
    "outbreaks", "notes", "deprecated",
]


# ── Data class ─────────────────────────────────────────────────────────────────

@dataclass
class Candidate:
    node_name: str
    leaves: list
    defining_muts: list       # branch mutations at this node
    support: Optional[float]
    earliest_date: str
    latest_date: str
    countries: str
    outbreaks: list           # outbreak names (lineage col) present in this clade
    seeded: bool = False      # True if designated via --seed-from-outbreaks


# ── I/O helpers ───────────────────────────────────────────────────────────────

def load_tree(path: str):
    return Phylo.read(path, "newick")


def load_nt_muts(path: str) -> dict:
    """Return {node_name: [mut_str, ...]} from augur ancestral JSON."""
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    return {
        name: info.get("muts", [])
        for name, info in data.get("nodes", {}).items()
    }


def load_metadata(path: str) -> dict:
    """Return {seq_id: row_dict} from metadata TSV."""
    meta: dict = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[row["id"]] = row
    return meta


def load_register(path: str) -> list:
    """Load existing lineage register TSV. Returns [] if file absent."""
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ── Tree helpers ───────────────────────────────────────────────────────────────

def build_parent_map(tree) -> dict:
    """Return {child_node_name: parent_node_name | None} for every node."""
    parent_map: dict = {tree.root.name: None}
    stack = [tree.root]
    while stack:
        node = stack.pop()
        for child in node.clades:
            parent_map[child.name] = node.name
            stack.append(child)
    return parent_map


def get_leaves(node) -> list:
    return sorted(t.name for t in node.get_terminals())


def get_support(node) -> Optional[float]:
    conf = getattr(node, "confidence", None)
    if conf is None:
        return None
    try:
        return float(conf)
    except (TypeError, ValueError):
        return None


def find_mrca(tree, tip_names: list):
    """Return the MRCA clade for a list of tip names, or None."""
    clades = [c for c in tree.get_terminals() if c.name in tip_names]
    if not clades:
        return None
    if len(clades) == 1:
        return clades[0]
    return tree.common_ancestor(clades)


# ── Date helpers ───────────────────────────────────────────────────────────────

def _parse_date_tuple(s: str) -> Optional[tuple]:
    """Parse 'YYYY-MM-DD' → (year, month, day) or None."""
    if not s or s in ("NA", "XXXX-XX-XX", ""):
        return None
    try:
        parts = s.split("-")
        if len(parts) != 3:
            return None
        return tuple(int(p) if (p.isdigit() and int(p) != 0) else 0 for p in parts)
    except Exception:
        return None


def _date_str(dt: tuple) -> str:
    return f"{dt[0]:04d}-{dt[1]:02d}-{dt[2]:02d}"


def _leaf_date_tuples(leaf_ids: list, meta: dict) -> list:
    dates = []
    for lid in leaf_ids:
        d = _parse_date_tuple(meta.get(lid, {}).get("date", ""))
        if d and d[0] > 0:
            dates.append(d)
    return dates


def _leaf_countries(leaf_ids: list, meta: dict) -> str:
    cs = set()
    for lid in leaf_ids:
        c = meta.get(lid, {}).get("country", "")
        if c and c not in ("NA", "Human", ""):
            cs.add(c)
    return ",".join(sorted(cs))


def _leaf_outbreaks(leaf_ids: list, meta: dict) -> list:
    """Return sorted unique outbreak names (metadata 'lineage' column) in clade."""
    names: set = set()
    for lid in leaf_ids:
        lin = meta.get(lid, {}).get("lineage", "")
        if lin and lin not in ("NA", ""):
            names.add(lin)
    return sorted(names)


def _make_candidate(node, leaves: list, nt_muts: dict, meta: dict,
                    seeded: bool = False) -> "Candidate":
    muts = sorted(nt_muts.get(node.name, []))
    dates = _leaf_date_tuples(leaves, meta)
    return Candidate(
        node_name=node.name,
        leaves=leaves,
        defining_muts=muts,
        support=get_support(node),
        earliest_date=_date_str(min(dates)) if dates else "NA",
        latest_date=_date_str(max(dates)) if dates else "NA",
        countries=_leaf_countries(leaves, meta),
        outbreaks=_leaf_outbreaks(leaves, meta),
        seeded=seeded,
    )


# ── Candidate discovery ────────────────────────────────────────────────────────

def find_qualifying_nodes(
    tree,
    nt_muts: dict,
    meta: dict,
    min_seqs: int,
    min_support: float,
    min_muts: int,
    seed_from_outbreaks: bool,
    filter_genotype: Optional[str] = None,
) -> list:
    """Return list of Candidates (qualifying internal nodes, not the root)."""
    candidates: list = []
    seen_nodes: set = set()

    filter_gt = filter_genotype.upper() if filter_genotype else None

    def _passes_genotype_filter(leaves: list) -> bool:
        """True if all leaves with a known genotype have the target genotype."""
        if filter_gt is None:
            return True
        for leaf in leaves:
            gt = meta.get(leaf, {}).get("genotype", "").upper()
            if gt and gt not in ("NA", "") and gt != filter_gt:
                return False
        return True

    # ── Seed from existing outbreak clusters ──────────────────────────────────
    if seed_from_outbreaks:
        outbreak_tips: dict = defaultdict(list)
        for leaf in tree.get_terminals():
            lin = meta.get(leaf.name, {}).get("lineage", "")
            if lin and lin not in ("NA", ""):
                outbreak_tips[lin].append(leaf.name)

        for outbreak, tips in sorted(outbreak_tips.items()):
            if len(tips) < 2:
                continue
            mrca = find_mrca(tree, tips)
            if mrca is None or mrca.is_terminal() or not mrca.name:
                continue
            if mrca.name == tree.root.name:
                continue  # don't designate the root
            if mrca.name in seen_nodes:
                continue
            leaves = get_leaves(mrca)
            if not _passes_genotype_filter(leaves):
                continue
            seen_nodes.add(mrca.name)
            candidates.append(_make_candidate(mrca, leaves, nt_muts, meta, seeded=True))

    # ── Criteria-based discovery (BFS / level-order) ──────────────────────────
    for node in tree.find_clades(order="level"):
        if node.is_terminal() or not node.name:
            continue
        if node.name == tree.root.name:
            continue
        if node.name in seen_nodes:
            continue

        leaves = get_leaves(node)
        if not _passes_genotype_filter(leaves):
            continue
        if len(leaves) < min_seqs:
            continue
        muts = nt_muts.get(node.name, [])
        if len(muts) < min_muts:
            continue
        support = get_support(node)
        if support is not None and min_support > 0 and support < min_support:
            continue

        seen_nodes.add(node.name)
        candidates.append(_make_candidate(node, leaves, nt_muts, meta))

    return candidates


# ── Hierarchy ──────────────────────────────────────────────────────────────────

def build_lineage_hierarchy(candidates: list, parent_map: dict) -> dict:
    """
    For each candidate, find its closest qualifying ancestor (lineage parent).
    Returns {node_name: parent_node_name | None}.
    None means top-level lineage (no qualifying ancestor).
    """
    cand_names = {c.node_name for c in candidates}
    hierarchy: dict = {}
    for cand in candidates:
        current = parent_map.get(cand.node_name)  # immediate tree parent
        found_parent = None
        while current is not None:
            if current in cand_names:
                found_parent = current
                break
            current = parent_map.get(current)
        hierarchy[cand.node_name] = found_parent
    return hierarchy


# ── Name assignment ────────────────────────────────────────────────────────────

def _muts_key(muts: list) -> frozenset:
    return frozenset(m for m in muts if m)


def assign_lineage_names(
    candidates: list,
    hierarchy: dict,
    genotype: str,
    existing_register: list,
    max_depth: int,
) -> dict:
    """
    Assign hierarchical lineage names (IA.1, IA.1.1, …) to candidates.

    Existing names are preserved by matching on defining_nt_muts first, then
    founder_node. New candidates receive fresh numbers. Numbers at each level
    are assigned in chronological order (earliest collection date).

    Returns {node_name: lineage_name}.
    """
    gt_prefix = genotype.upper()

    # Build lookup tables from existing register (same genotype, not deprecated)
    existing_by_muts: dict = {}
    existing_by_node: dict = {}
    for row in existing_register:
        if row.get("deprecated", "FALSE").upper() == "TRUE":
            continue
        if row.get("genotype", "").upper() != gt_prefix:
            continue
        muts_str = row.get("defining_nt_muts", "")
        mkey = _muts_key(muts_str.split(",") if muts_str else [])
        if mkey:
            existing_by_muts[mkey] = row["lineage"]
        fn = row.get("founder_node", "")
        if fn:
            existing_by_node[fn] = row["lineage"]

    # Build children lists within the lineage hierarchy
    children: dict = defaultdict(list)
    roots: list = []
    for cand in candidates:
        parent = hierarchy.get(cand.node_name)
        if parent is None:
            roots.append(cand)
        else:
            children[parent].append(cand)

    def _sort_key(c: Candidate) -> str:
        return c.earliest_date if c.earliest_date != "NA" else "9999-99-99"

    def _get_existing_name(cand: Candidate) -> Optional[str]:
        mkey = _muts_key(cand.defining_muts)
        if mkey and mkey in existing_by_muts:
            return existing_by_muts[mkey]
        return existing_by_node.get(cand.node_name)

    name_map: dict = {}

    def _assign_level(node_list: list, prefix: str, depth: int) -> None:
        if depth > max_depth or not node_list:
            return

        node_list = sorted(node_list, key=_sort_key)

        # Collect numbers already reserved by pre-existing names at this level
        taken: set = set()
        for cand in node_list:
            existing = _get_existing_name(cand)
            if existing:
                parts = existing.split(".")
                prefix_parts = prefix.split(".")
                # Only grab the next-level number (e.g. "1" from "IA.1" when prefix="IA")
                if len(parts) == len(prefix_parts) + 1 and existing.startswith(prefix + "."):
                    try:
                        taken.add(int(parts[-1]))
                    except ValueError:
                        pass

        counter = 1
        for cand in node_list:
            existing = _get_existing_name(cand)
            if existing:
                name_map[cand.node_name] = existing
            else:
                while counter in taken:
                    counter += 1
                new_name = f"{prefix}.{counter}"
                taken.add(counter)
                counter += 1
                name_map[cand.node_name] = new_name

            assigned = name_map[cand.node_name]
            child_list = sorted(children.get(cand.node_name, []), key=_sort_key)
            _assign_level(child_list, assigned, depth + 1)

    roots.sort(key=_sort_key)
    _assign_level(roots, gt_prefix, 1)
    return name_map


# ── Register I/O ───────────────────────────────────────────────────────────────

def bridge_source_to_clades_tsv(
    tree,
    nt_muts: dict,
    source_assignments: dict,  # {seq_id: lineage_name}
    path: str,
) -> None:
    """
    Bridge mode for the all-genotype tree.

    For each unique lineage name in source_assignments (e.g. 'IA.1', 'IB.1.1'),
    find the MRCA of all tips with that name in this tree, then write the MRCA's
    branch mutations as an augur clades-format TSV.  The original per-genotype
    lineage names are preserved verbatim — no 'ALL.*' renaming occurs.

    This re-expresses per-genotype lineage definitions in the all-tree's reference
    (NC_001489.1) coordinates so that 'augur clades' can annotate the all-tree
    correctly without any coordinate bridging.
    """
    import re
    MUT_RE = re.compile(r'^([A-Za-z-])([0-9]+)([A-Za-z-])$')

    # Group tip IDs by lineage name
    tips_by_lineage: dict = defaultdict(list)
    for seq_id, lin in source_assignments.items():
        if lin:
            tips_by_lineage[lin].append(seq_id)

    # Only keep tips that are actually in this tree
    tree_tips = {t.name for t in tree.get_terminals()}

    rows: list = []
    n_placed = 0
    n_skipped = 0
    for lineage_name, tip_ids in sorted(tips_by_lineage.items()):
        present = [t for t in tip_ids if t in tree_tips]
        if len(present) < 2:
            n_skipped += 1
            continue
        mrca = find_mrca(tree, present)
        if mrca is None or not mrca.name:
            n_skipped += 1
            continue
        muts = nt_muts.get(mrca.name, [])
        for mut in muts:
            m = MUT_RE.match(mut.strip())
            if not m:
                continue
            _ref, site, alt = m.group(1), m.group(2), m.group(3)
            rows.append({"clade": lineage_name, "gene": "nuc",
                         "site": site, "alt": alt})
        n_placed += 1

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["clade", "gene", "site", "alt"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(
        f"  Bridge clades TSV  : {path}  "
        f"({len(rows)} mutation rows, {n_placed} lineages placed, "
        f"{n_skipped} skipped <2 tips in tree)",
        file=sys.stderr,
    )


def write_augur_clades_tsv(name_map: dict, candidates: list, path: str) -> None:
    """
    Write a clade-definitions TSV in augur clades format from the designation results.

    augur clades format (tab-separated):
        clade   gene    site    alt
        IA.1    nuc     45      G

    Mutations come from the 'defining_muts' field of each Candidate (branch mutations
    at the founding node, expressed relative to the tree's reference sequence in the
    same 1-based coordinate system that augur ancestral uses).
    """
    import re
    MUT_RE = re.compile(r'^([A-Za-z-])([0-9]+)([A-Za-z-])$')

    cand_by_node = {c.node_name: c for c in candidates}

    rows: list = []
    for node_name, lineage_name in sorted(name_map.items(), key=lambda x: x[1]):
        cand = cand_by_node.get(node_name)
        if cand is None:
            continue
        for mut in cand.defining_muts:
            m = MUT_RE.match(mut.strip())
            if not m:
                continue
            _ref, site, alt = m.group(1), m.group(2), m.group(3)
            rows.append({"clade": lineage_name, "gene": "nuc",
                         "site": site, "alt": alt})

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["clade", "gene", "site", "alt"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def get_tip_assignments(tree, name_map: dict) -> dict:
    """
    Return {leaf_name: lineage_phylo} by walking the tree pre-order from the root.
    Each leaf inherits the lineage of its nearest designated ancestor.
    """
    assignments: dict = {}

    def _walk(node, lineage: Optional[str]) -> None:
        assigned = name_map.get(node.name) if node.name else None
        current = assigned if assigned else lineage
        if node.is_terminal():
            if current and node.name:
                assignments[node.name] = current
        else:
            for child in node.clades:
                _walk(child, current)

    _walk(tree.root, None)
    return assignments


def build_register_rows(
    name_map: dict,
    candidates: list,
    genotype: str,
    existing_register: list,
) -> tuple:
    """
    Build the full updated register and the list of newly proposed lineages.

    Returns (all_rows, new_rows):
        all_rows  — complete register to write back to lineage_designations.tsv
        new_rows  — only newly proposed lineages (for lineage_proposals.tsv)
    """
    gt_prefix = genotype.upper()
    cand_by_node = {c.node_name: c for c in candidates}
    assigned_names = set(name_map.values())

    # Index existing register for this genotype by lineage name
    existing_by_lineage: dict = {}
    other_genotype_rows: list = []
    for row in existing_register:
        if row.get("genotype", "").upper() == gt_prefix:
            existing_by_lineage[row["lineage"]] = row
        else:
            other_genotype_rows.append(row)

    # Build updated rows for this genotype
    current_rows: list = []
    new_rows: list = []

    for node_name, lineage_name in sorted(name_map.items(), key=lambda x: x[1]):
        cand = cand_by_node[node_name]
        is_new = lineage_name not in existing_by_lineage
        existing_notes = existing_by_lineage.get(lineage_name, {}).get("notes", "")
        row: dict = {
            "lineage":          lineage_name,
            "genotype":         gt_prefix,
            "founder_node":     node_name,
            "defining_nt_muts": ",".join(cand.defining_muts),
            "num_sequences":    str(len(cand.leaves)),
            "earliest_date":    cand.earliest_date,
            "latest_date":      cand.latest_date,
            "countries":        cand.countries,
            "outbreaks":        ",".join(cand.outbreaks),
            "notes":            existing_notes or ("seeded" if cand.seeded else ""),
            "deprecated":       "FALSE",
        }
        current_rows.append(row)
        if is_new:
            new_rows.append(row)

    # Mark lineages no longer found in the tree as deprecated
    deprecated_rows: list = []
    for lineage_name, row in existing_by_lineage.items():
        if lineage_name not in assigned_names:
            deprecated_row = dict(row)
            deprecated_row["deprecated"] = "TRUE"
            deprecated_rows.append(deprecated_row)

    all_rows = other_genotype_rows + current_rows + deprecated_rows
    return all_rows, new_rows


def write_tsv(path: str, rows: list, cols: list) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=cols, delimiter="\t",
            extrasaction="ignore", lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Designate phylogenetic lineage names for HAV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("tree_file",     help="Newick tree from augur refine")
    parser.add_argument("nt_muts_file",  help="Node mutations JSON from augur ancestral")
    parser.add_argument("metadata_file", help="TSV with id, genotype, lineage, date, country")
    parser.add_argument("genotype",      help="Genotype label (e.g. Ia, Ib, IIa)")
    parser.add_argument(
        "--register", default="data/lineage_designations.tsv",
        help="Path to lineage register TSV (created/updated). "
             "Default: data/lineage_designations.tsv",
    )
    parser.add_argument(
        "--proposals", default="data/lineage_proposals.tsv",
        help="Path to write new candidate proposals. "
             "Default: data/lineage_proposals.tsv",
    )
    parser.add_argument(
        "--seed-from-outbreaks", action="store_true", default=False,
        help="Designate the MRCA of each existing outbreak cluster "
             "(metadata 'lineage' column) as a founding node, bypassing "
             "normal qualification criteria. Recommended for the initial run.",
    )
    parser.add_argument("--min-sequences", type=int,   default=DEFAULT_MIN_SEQUENCES,
                        help=f"Min leaf count per clade (default: {DEFAULT_MIN_SEQUENCES})")
    parser.add_argument("--min-support",   type=float, default=DEFAULT_MIN_SUPPORT,
                        help=f"Min bootstrap support 0–100 (default: {DEFAULT_MIN_SUPPORT}; 0=disable)")
    parser.add_argument("--min-muts",      type=int,   default=DEFAULT_MIN_DEFINING_MUTS,
                        help=f"Min branch mutations at node (default: {DEFAULT_MIN_DEFINING_MUTS})")
    parser.add_argument("--max-depth",     type=int,   default=DEFAULT_MAX_DEPTH,
                        help=f"Max nesting depth (default: {DEFAULT_MAX_DEPTH})")
    parser.add_argument(
        "--assignments", default=None, metavar="PATH",
        help="Path to write a per-sequence lineage_phylo TSV (columns: id, lineage_phylo).",
    )
    parser.add_argument(
        "--augur-clades-tsv", default=None, metavar="PATH",
        help="Path to write an augur clades-format TSV (clade / gene / site / alt) derived "
             "from the defining mutations of each designated lineage.  Pass this file to "
             "'augur clades' to assign lineage_phylo via the official augur toolchain.",
    )
    parser.add_argument(
        "--filter-genotype", default=None, metavar="GT",
        help="Restrict lineage designation to nodes where all tips with a known "
             "genotype belong to this genotype (e.g. Ia, Ib). Useful when running "
             "on the all-genotype tree with --seed-from-outbreaks to prevent "
             "cross-genotype designation.",
    )
    parser.add_argument(
        "--lineage-source-tsv", default=None, metavar="PATH",
        help="TSV with columns 'id' and 'lineage_phylo'.  When provided, the 'lineage' "
             "values in metadata are overridden with these assignments before running "
             "designation.  With --bridge-from-source, these assignments are used directly "
             "as lineage names (no new naming).",
    )
    parser.add_argument(
        "--bridge-from-source", action="store_true", default=False,
        help="Bridge mode for the all-genotype tree.  Requires --lineage-source-tsv and "
             "--augur-clades-tsv.  Instead of running the full designation algorithm, finds "
             "the MRCA of each lineage group in THIS tree and writes branch mutations to the "
             "augur clades TSV using the original lineage names (IA.1, IB.1 etc.).  No "
             "register or proposals are written.",
    )
    args = parser.parse_args()

    # ── Load inputs ───────────────────────────────────────────────────────────
    print(f"  Loading tree        : {args.tree_file}", file=sys.stderr)
    tree = load_tree(args.tree_file)

    print(f"  Loading mutations   : {args.nt_muts_file}", file=sys.stderr)
    nt_muts = load_nt_muts(args.nt_muts_file)

    print(f"  Loading metadata    : {args.metadata_file}", file=sys.stderr)
    meta = load_metadata(args.metadata_file)

    # Load lineage-source-tsv (used for bridge mode and/or seeding overrides)
    source_lineages: dict = {}
    if args.lineage_source_tsv:
        with open(args.lineage_source_tsv, encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                seq_id = (row.get("id") or "").strip()
                lin    = (row.get("lineage_phylo") or "").strip()
                if seq_id and lin:
                    source_lineages[seq_id] = lin
        print(f"  Lineage source      : {args.lineage_source_tsv} "
              f"({len(source_lineages)} assignments loaded)", file=sys.stderr)

    # ── Bridge mode: map existing lineage names onto this tree's mutations ────
    # Used for the all-genotype tree so that IA.1, IB.1 etc. are preserved as
    # names rather than creating ALL.1, ALL.2 etc.
    if args.bridge_from_source:
        if not args.lineage_source_tsv or not args.augur_clades_tsv:
            print("ERROR: --bridge-from-source requires both --lineage-source-tsv "
                  "and --augur-clades-tsv", file=sys.stderr)
            sys.exit(1)
        bridge_source_to_clades_tsv(tree, nt_muts, source_lineages, args.augur_clades_tsv)
        return

    # Override metadata lineage column from source assignments (seeding)
    if source_lineages:
        n_overrides = 0
        for seq_id, lin in source_lineages.items():
            if seq_id in meta:
                meta[seq_id]["lineage"] = lin
                n_overrides += 1
        print(f"  Lineage overrides   : {n_overrides} applied", file=sys.stderr)

    existing_register = load_register(args.register)
    n_existing = sum(
        1 for r in existing_register
        if r.get("genotype", "").upper() == args.genotype.upper()
        and r.get("deprecated", "FALSE").upper() != "TRUE"
    )

    n_tips     = sum(1 for _ in tree.get_terminals())
    n_internal = sum(1 for _ in tree.get_nonterminals())
    print(f"  Tree                : {n_tips} tips, {n_internal} internal nodes",
          file=sys.stderr)
    print(f"  Register            : {len(existing_register)} total rows, "
          f"{n_existing} active for {args.genotype}", file=sys.stderr)

    # ── Discover qualifying nodes ─────────────────────────────────────────────
    candidates = find_qualifying_nodes(
        tree, nt_muts, meta,
        min_seqs=args.min_sequences,
        min_support=args.min_support,
        min_muts=args.min_muts,
        seed_from_outbreaks=args.seed_from_outbreaks,
        filter_genotype=args.filter_genotype,
    )
    n_seeded = sum(1 for c in candidates if c.seeded)
    print(f"  Qualifying nodes    : {len(candidates)} "
          f"({n_seeded} seeded, {len(candidates) - n_seeded} criteria-based)",
          file=sys.stderr)

    if not candidates:
        print(
            "  No qualifying nodes found. Try:\n"
            "    --seed-from-outbreaks  to use existing outbreak clusters\n"
            "    --min-sequences 2      to lower the sequence threshold\n"
            "    --min-muts 0           to remove the mutation requirement",
            file=sys.stderr,
        )
        sys.exit(0)

    # ── Build hierarchy and assign names ─────────────────────────────────────
    parent_map = build_parent_map(tree)
    hierarchy  = build_lineage_hierarchy(candidates, parent_map)
    name_map   = assign_lineage_names(
        candidates, hierarchy, args.genotype, existing_register, args.max_depth,
    )

    # ── Print summary ─────────────────────────────────────────────────────────
    print(f"\n  {'Lineage':<16} {'Node':<18} {'N':>4}  Date range"
          f"  {'Outbreaks (first 3)'}",
          file=sys.stderr)
    print(f"  {'─'*80}", file=sys.stderr)
    for cand in sorted(candidates, key=lambda c: name_map.get(c.node_name, "")):
        name = name_map.get(cand.node_name, "?")
        tag  = "*" if cand.seeded else " "
        outbreaks_preview = ",".join(cand.outbreaks[:3])
        if len(cand.outbreaks) > 3:
            outbreaks_preview += f" (+{len(cand.outbreaks) - 3})"
        print(
            f"  {name:<16} {cand.node_name:<18} {len(cand.leaves):>4}"
            f"  {cand.earliest_date}–{cand.latest_date}"
            f"  {outbreaks_preview}{tag}",
            file=sys.stderr,
        )
    print(f"  (* = seeded from outbreak cluster)", file=sys.stderr)

    # ── Write outputs ─────────────────────────────────────────────────────────
    all_rows, new_rows = build_register_rows(
        name_map, candidates, args.genotype, existing_register,
    )

    write_tsv(args.register, all_rows, REGISTER_COLS)
    print(f"\n  Register written : {args.register}  ({len(all_rows)} rows)",
          file=sys.stderr)

    write_tsv(args.proposals, new_rows, REGISTER_COLS)
    if new_rows:
        print(f"  Proposals written: {args.proposals}  ({len(new_rows)} new lineages)",
              file=sys.stderr)
        print(f"\n  → Review proposals with: "
              f"python3 scripts/review_lineage_proposals.py {args.proposals}",
              file=sys.stderr)
    else:
        print(f"  No new proposals (all candidates already in register)",
              file=sys.stderr)

    # ── Per-sequence assignments (for use with the all-genotype tree) ─────────
    if args.assignments:
        tip_assignments = get_tip_assignments(tree, name_map)
        with open(args.assignments, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
            writer.writerow(["id", "lineage_phylo"])
            for seq_id, lin in sorted(tip_assignments.items()):
                writer.writerow([seq_id, lin])
        print(
            f"  Assignments written: {args.assignments}  ({len(tip_assignments)} sequences)",
            file=sys.stderr,
        )

    # ── augur clades TSV (official tool input) ────────────────────────────────
    if args.augur_clades_tsv:
        write_augur_clades_tsv(name_map, candidates, args.augur_clades_tsv)
        n_rows = sum(len(c.defining_muts) for c in candidates if c.node_name in name_map)
        print(
            f"  augur clades TSV   : {args.augur_clades_tsv}  ({n_rows} mutation rows)",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
