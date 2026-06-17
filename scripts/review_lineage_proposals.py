#!/usr/bin/env python3
"""
review_lineage_proposals.py — Display new HAV lineage proposals for curator review.

Reads lineage_proposals.tsv written by designate_lineages.py and prints a
human-readable summary. Proposals are already pre-added to the register
(lineage_designations.tsv), so no action is needed to accept them.

To REJECT a proposal:
    1. Open data/lineage_designations.tsv
    2. Find the row for the lineage
    3. Set deprecated=TRUE

Usage:
    python3 scripts/review_lineage_proposals.py [proposals_file]

Default proposals file: data/lineage_proposals.tsv
"""

import csv
import os
import sys

PROPOSALS_DEFAULT = "data/lineage_proposals.tsv"
REGISTER_DEFAULT  = "data/lineage_designations.tsv"


def main() -> None:
    proposals_file = sys.argv[1] if len(sys.argv) > 1 else PROPOSALS_DEFAULT

    if not os.path.exists(proposals_file):
        print(f"No proposals file found at: {proposals_file}")
        print("Run designate_lineages.py first to generate proposals.")
        return

    rows: list = []
    with open(proposals_file, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    if not rows:
        print("No pending proposals — nothing to review.")
        return

    # Group by genotype for cleaner display
    by_genotype: dict = {}
    for row in rows:
        gt = row.get("genotype", "?")
        by_genotype.setdefault(gt, []).append(row)

    print()
    print("=" * 72)
    print(f"  HAV Lineage Proposals  ({len(rows)} new designation(s))")
    print("=" * 72)

    for gt in sorted(by_genotype):
        gt_rows = by_genotype[gt]
        print(f"\n  Genotype {gt}  ({len(gt_rows)} proposal(s))")
        print("  " + "─" * 68)

        for row in sorted(gt_rows, key=lambda r: r["lineage"]):
            seeded = "seeded" in row.get("notes", "").lower()
            origin = "seeded from outbreak cluster" if seeded else "criteria-based discovery"

            print(f"\n  Lineage  : {row['lineage']}  [{origin}]")
            print(f"  Node     : {row['founder_node']}")

            muts = row.get("defining_nt_muts", "")
            if muts:
                # Wrap long mutation lists
                muts_list = muts.split(",")
                if len(muts_list) <= 6:
                    print(f"  Mutations: {muts}")
                else:
                    print(f"  Mutations: {','.join(muts_list[:6])} (+{len(muts_list)-6} more)")
            else:
                print(f"  Mutations: (none)")

            print(f"  Sequences: {row['num_sequences']}")
            print(f"  Dates    : {row['earliest_date']}  →  {row['latest_date']}")

            countries = row.get("countries", "")
            if countries:
                print(f"  Countries: {countries}")

            outbreaks_str = row.get("outbreaks", "")
            if outbreaks_str:
                outbreaks = outbreaks_str.split(",")
                preview   = ", ".join(outbreaks[:5])
                extra     = f" (+{len(outbreaks) - 5} more)" if len(outbreaks) > 5 else ""
                print(f"  Outbreaks: {preview}{extra}")

    print()
    print("─" * 72)
    print("  These proposals have been pre-added to the register.")
    print(f"  Register : {REGISTER_DEFAULT}")
    print()
    print("  To REJECT a proposal:")
    print(f"    1. Open {REGISTER_DEFAULT}")
    print("    2. Find the lineage row and set  deprecated=TRUE")
    print()
    print("  To apply accepted proposals, rebuild the lineages dataset:")
    print("    bash scripts/build_vp1b_lineage_dataset.sh")
    print("─" * 72)
    print()


if __name__ == "__main__":
    main()
