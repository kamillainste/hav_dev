"""Relabel ViroConstrictor's amplicon-coverage columns by amplicon NUMBER.

Usage: python relabel_amplicons.py <ampliconcoverage.csv> <genotype>

ViroConstrictor (amplicon_covs.py, _create_amplicon_names_list) takes the prefix
of the FIRST row of the coordinate-sorted primer BED and stamps it on every
amplicon. The prefix therefore says nothing about which primer defined an
amplicon: an IIIA sample comes out labelled HAVIB_* because HAVIB_1_LEFT maps
~24 bp upstream of HAVIIIA_1_LEFT and sorts first. Column count is also not
guaranteed to be six, nor the numbers contiguous, so renaming by position would
mislabel every column after a missing one. This matches on the number inside
the original name, and records any column whose original genotype prefix
differs from the assigned reference in <csv>_rename_notes.tsv.

Never fatal for the pipeline: on any error it warns and leaves the file as is.
"""
import re
import sys

try:
    import pandas as pd

    path, gt = sys.argv[1], sys.argv[2]
    df = pd.read_csv(path, index_col=0)

    new_cols, notes = [], []
    for col in df.columns:
        m = re.match(r"^HAV([A-Za-z0-9]+)_0*(\d+)\s*$", str(col).strip())
        if not m:
            new_cols.append(col)
            notes.append(f"{col}\tunparsed\tleft unchanged")
            continue
        src_gt, num = m.group(1), int(m.group(2))
        new_cols.append(f"HAV{gt}_{num:03d}")
        if src_gt != gt:
            notes.append(
                f"{col}\tHAV{gt}_{num:03d}\tprefix was {src_gt}; ViroConstrictor names EVERY "
                "amplicon after the lowest-coordinate primer in the BED, so this is not a "
                "per-amplicon match"
            )

    if len(set(new_cols)) != len(new_cols):
        print("WARNING: amplicon rename would collide; leaving file unchanged", file=sys.stderr)
    else:
        df.columns = new_cols
        df.to_csv(path)

    expected = set(range(1, 7))
    present = {int(c.split("_")[-1]) for c in new_cols if re.match(r"^HAV.+_\d+$", str(c))}
    missing = sorted(expected - present)
    if missing:
        notes.append(f"-\t-\tamplicon(s) absent from ViroConstrictor output: {missing}")

    if notes:
        with open(path.replace(".csv", "_rename_notes.tsv"), "w") as fh:
            fh.write("original_column\tnew_column\tnote\n")
            fh.write("\n".join(notes) + "\n")
        for n in notes:
            print("  NOTE amplicon rename: " + n.replace("\t", " | "), file=sys.stderr)
except Exception as e:  # reporting step: never take the sample down
    print(f"WARNING: Failed to relabel amplicon coverage columns: {e}", file=sys.stderr)
