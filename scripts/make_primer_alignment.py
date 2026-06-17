#!/usr/bin/env python3
"""
Create a positional MSA FASTA showing primer binding sites on the reference.

Uses the best (highest match_fraction) binding site per primer from the GFF3 file,
and positions each primer sequence (RC for reverse-strand primers) at its binding
coordinates relative to the reference region spanned by the primers.

Output: data/Primers/PCR_primers/primer_alignment.fasta (gap-padded MSA, same length for all seqs)
"""

import argparse
import re
from pathlib import Path

# ── IUPAC complements (needed to reverse-complement degenerate primers) ────────
RC_MAP = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R",   # puRine  <-> pYrimidine
    "S": "S", "W": "W",   # Strong  <-> Weak
    "K": "M", "M": "K",   # Keto    <-> aMino
    "B": "V", "V": "B",   # not-A   <-> not-T
    "D": "H", "H": "D",   # not-C   <-> not-G
    "N": "N",              # aNy
    "5": "N",              # '5' codehop degenerate treated as N
    "I": "I",              # Inosine
    "-": "-", ".": ".",
}


def reverse_complement(seq: str) -> str:
    return "".join(RC_MAP.get(b.upper(), "N") for b in reversed(seq.upper()))


def read_fasta(path: Path) -> dict[str, str]:
    """Return {header_line: sequence} (header includes '>' prefix stripped)."""
    seqs: dict[str, str] = {}
    current_header = None
    current_seq: list[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_header is not None:
                    seqs[current_header] = "".join(current_seq)
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
    if current_header is not None:
        seqs[current_header] = "".join(current_seq)
    return seqs


def parse_primer_name(header: str) -> tuple[str, str, str]:
    """Return (name, direction, comment) from primer FASTA header."""
    parts = header.split()
    name = parts[0]
    meta = dict(kv.split("=", 1) for kv in parts[1:] if "=" in kv)
    return name, meta.get("direction", ""), meta.get("comment", "")


def best_matches_from_gff3(gff3_path: Path) -> dict[str, dict]:
    """
    Return, per primer name, the single best match (highest match_fraction)
    found anywhere in the GFF3.  Ties broken by longest match_length.
    """
    best: dict[str, dict] = {}
    with open(gff3_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue
            start, end = int(parts[3]), int(parts[4])
            score = float(parts[5])
            strand = parts[6]
            attrs = parts[8]

            name_m = re.search(r"Name=([^;]+)", attrs)
            if not name_m:
                continue
            name = name_m.group(1)
            mlen = int(re.search(r"match_length=(\d+)", attrs).group(1))

            if name not in best or (score, mlen) > (best[name]["score"], best[name]["mlen"]):
                best[name] = {
                    "start": start, "end": end,
                    "score": score, "mlen": mlen, "strand": strand,
                }
    return best


def build_alignment(
    ref_seq: str,
    primer_seqs: dict[str, str],
    best: dict[str, dict],
    primer_meta: dict[str, tuple[str, str]],
    margin: int = 20,
) -> list[tuple[str, str]]:
    """
    Build gap-padded sequences aligned to the reference.

    Returns list of (header, gapped_seq) pairs.  The window shown is:
      [min_primer_start - margin, max_primer_end + margin]
    (1-based, inclusive, clipped to reference bounds)
    """
    primer_names = list(best.keys())
    win_start = max(1, min(best[n]["start"] for n in primer_names) - margin)
    win_end   = min(len(ref_seq), max(best[n]["end"] for n in primer_names) + margin)
    win_len   = win_end - win_start + 1

    records: list[tuple[str, str]] = []

    # Reference
    ref_region = ref_seq[win_start - 1 : win_end]  # 0-based slice
    records.append((f"Reference NC_001489.1 [{win_start}-{win_end}]", ref_region))

    # Primers (order: outer-F, inner-F, inner-R, outer-R)
    primer_order = [
        ("HAV_6.1_codehop",  "First-outer",   "F"),
        ("HAV_8.2_codehop",  "nested-inner",  "F"),
        ("HAV_11_codehop",   "nested-inner",  "R"),
        ("HAV_10_codehop",   "First-outer",   "R"),
    ]

    for pname, comment, expected_dir in primer_order:
        if pname not in best:
            continue
        m = best[pname]
        p_start, p_end = m["start"], m["end"]   # 1-based GFF3 coords
        strand = m["strand"]

        # Primer sequence (from FASTA)
        raw_seq = primer_seqs.get(pname, "")
        if not raw_seq:
            continue

        # For minus-strand hits, the primer as written anneals to the
        # complementary strand → reverse-complement it to read 5'→3' along
        # the reference so it lines up visually.
        if strand == "-":
            display_seq = reverse_complement(raw_seq)
            rc_note = " [RC]"
        else:
            display_seq = raw_seq.upper()
            rc_note = ""

        # Pad with gaps to position within the window
        left_pad  = p_start - win_start       # gaps before primer
        right_pad = win_len - left_pad - len(display_seq)
        if right_pad < 0:
            # Primer extends past window (shouldn't happen with margin ≥ 0)
            display_seq = display_seq[:win_len - left_pad]
            right_pad = 0

        gapped = "-" * left_pad + display_seq + "-" * right_pad
        assert len(gapped) == win_len, f"{pname}: {len(gapped)} vs {win_len}"

        dir_tag = f" (direction={expected_dir})"
        score_tag = f" score={m['score']:.2f} match={m['mlen']}bp"
        header = (
            f"{pname}{rc_note} [{p_start}-{p_end}]"
            f"{dir_tag} {comment}{score_tag}"
        )
        records.append((header, gapped))

    return records, win_start, win_end


def write_fasta(records: list[tuple[str, str]], output: Path, line_width: int = 60):
    with open(output, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i : i + line_width] + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference-fasta", required=True)
    ap.add_argument("--primers-fasta", required=True)
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--output-fasta", required=True)
    ap.add_argument("--margin", type=int, default=20,
                    help="Flanking reference bases to include on each side (default: 20)")
    args = ap.parse_args()

    # ── Load data ──────────────────────────────────────────────────────────────
    ref_data = read_fasta(Path(args.reference_fasta))
    ref_header, ref_seq = next(iter(ref_data.items()))
    ref_seq = ref_seq.upper()
    print(f"Reference: {ref_header[:60]} ({len(ref_seq)} bp)")

    primer_data = read_fasta(Path(args.primers_fasta))
    primer_seqs: dict[str, str] = {}
    primer_meta: dict[str, tuple[str, str]] = {}
    for hdr, seq in primer_data.items():
        name, direction, comment = parse_primer_name(hdr)
        primer_seqs[name] = seq.upper()
        primer_meta[name] = (direction, comment)
        print(f"  Primer: {name:25s} {direction}  ({comment})  {len(seq)} bp")

    best = best_matches_from_gff3(Path(args.gff3))
    print("\nBest binding sites per primer:")
    for name, m in sorted(best.items()):
        print(f"  {name:25s}  {m['start']:5d}-{m['end']:5d}  "
              f"strand={m['strand']}  score={m['score']:.2f}  match_len={m['mlen']}")

    # ── Build and write alignment ───────────────────────────────────────────────
    records, win_start, win_end = build_alignment(
        ref_seq, primer_seqs, best, primer_meta, margin=args.margin
    )

    write_fasta(records, Path(args.output_fasta))
    print(f"\nAlignment window: {win_start}-{win_end}  ({win_end - win_start + 1} columns)")
    print(f"Wrote {len(records)} sequences to: {args.output_fasta}")

    # ── ASCII overview ─────────────────────────────────────────────────────────
    print("\nPositional overview (reference window):")
    print(f"{'Position':>8}  {'Sequence':<40}")
    print("-" * 60)
    for header, seq in records:
        # Find non-gap span
        stripped = seq.lstrip("-")
        left_gaps = len(seq) - len(stripped)
        right_gaps = len(seq.rstrip("-")) - (len(seq) - len(stripped))
        primer_len = len(seq) - left_gaps - (len(seq) - len(seq.rstrip("-")))
        abs_start = win_start + left_gaps
        abs_end   = win_start + len(seq) - len(seq.rstrip("-")) - 1
        label = header.split("[")[0].strip()[:35]
        print(f"{abs_start:>8}-{abs_end:<6}  {label}")


if __name__ == "__main__":
    main()
