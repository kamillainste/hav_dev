#!/usr/bin/env python3
"""Generate GFF3 file of primer binding coordinates on reference sequence.

Finds primer matches on the reference using IUPAC-compatible base matching,
and outputs a GFF3 file with primer binding locations.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


IUPAC: Dict[str, set[str]] = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "U": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
    "I": {"A", "C", "G", "T"},
    "5": {"A", "C", "G", "T"},
}

RC_MAP = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "U": "A",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "D": "H",
    "H": "D",
    "V": "B",
    "N": "N",
    "I": "N",
    "5": "N",
}


@dataclass
class Primer:
    name: str
    sequence: str
    direction: str
    region: str = ""
    comment: str = ""
    source: str = ""


def parse_header(header: str) -> Tuple[str, Dict[str, str]]:
    """Parse FASTA header into name and key=value metadata."""
    parts = header.split()
    if not parts:
        raise ValueError("Invalid FASTA header")
    name = parts[0]
    metadata: Dict[str, str] = {}
    for token in parts[1:]:
        if "=" in token:
            k, v = token.split("=", 1)
            metadata[k.strip().lower()] = v.strip()
    return name, metadata


def read_fasta(path: str) -> List[Tuple[str, str]]:
    """Read FASTA file and return list of (name, sequence) tuples."""
    records: List[Tuple[str, str]] = []
    name: Optional[str] = None
    seq_parts: List[str] = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_parts).upper()))
                name = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.replace(" ", "").upper())
    if name is not None:
        records.append((name, "".join(seq_parts).upper()))
    return records


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    seq = seq.upper()
    bases = list(seq)
    rc_bases = [RC_MAP.get(b, "N") for b in bases]
    return "".join(reversed(rc_bases))


def base_compatible(seq_base: str, primer_base: str) -> bool:
    """Check if sequence base is compatible with primer base using IUPAC."""
    seq_base = seq_base.upper()
    primer_base = primer_base.upper()
    seq_chars = IUPAC.get(seq_base, {seq_base})
    primer_chars = IUPAC.get(primer_base, {primer_base})
    return len(seq_chars & primer_chars) > 0


def find_primer_matches(
    reference: str, primer: Primer, min_match_len: int = 8
) -> List[Dict]:
    """Find all matches of primer on reference and return match records."""
    matches = []
    ref = reference.upper()
    primer_seq = primer.sequence.upper()

    # Determine which sequence to search
    if primer.direction == "F":
        search_seq = primer_seq
    else:
        search_seq = reverse_complement(primer_seq)

    # Search for all possible matches
    for pos in range(len(ref) - len(search_seq) + 1):
        ref_window = ref[pos : pos + len(search_seq)]

        # Check compatibility at this position
        matches_here = sum(
            1
            for i in range(len(search_seq))
            if base_compatible(ref_window[i], search_seq[i])
        )

        # If match is good enough, record it
        if matches_here >= min_match_len:
            match_frac = matches_here / len(search_seq)
            matches.append(
                {
                    "primer_name": primer.name,
                    "primer_direction": primer.direction,
                    "primer_region": primer.region,
                    "primer_comment": primer.comment,
                    "primer_length": len(search_seq),
                    "ref_start": pos + 1,  # 1-based GFF3 convention
                    "ref_end": pos + len(search_seq),
                    "match_length": matches_here,
                    "match_fraction": match_frac,
                    "strand": "+" if primer.direction == "F" else "-",
                }
            )

    return matches


def load_primers(path: str) -> List[Primer]:
    """Load primers from FASTA file."""
    primers = []
    for header, seq in read_fasta(path):
        name, meta = parse_header(header)
        primer = Primer(
            name=name,
            sequence=seq,
            direction=meta.get("direction", "").upper(),
            region=meta.get("region", ""),
            comment=meta.get("comment", ""),
            source=path,
        )
        primers.append(primer)
    return primers


def write_gff3(matches: List[Dict], output_path: str, reference_name: str) -> None:
    """Write GFF3 file of primer matches."""
    with open(output_path, "w", encoding="utf-8") as f:
        # GFF3 header
        f.write("##gff-version 3\n")
        f.write(f"##sequence-region {reference_name} 1 {reference_name}\n")

        # Sort matches by start position
        matches.sort(key=lambda x: x["ref_start"])

        # Write records
        for i, match in enumerate(matches, 1):
            seqname = reference_name
            source = "primer_alignment"
            feature = "primer_binding_site"
            start = match["ref_start"]
            end = match["ref_end"]
            score = match["match_fraction"]
            strand = match["strand"]
            phase = "."

            # Build attributes
            attributes = (
                f"ID=primer_{i};"
                f"Name={match['primer_name']};"
                f"match_length={match['match_length']};"
                f"primer_length={match['primer_length']};"
                f"match_fraction={match['match_fraction']:.2f}"
            )
            if match["primer_region"]:
                attributes += f";region={match['primer_region']}"
            if match["primer_comment"]:
                attributes += f";comment={match['primer_comment']}"

            # Write GFF3 line
            fields = [
                seqname,
                source,
                feature,
                str(start),
                str(end),
                f"{score:.3f}",
                strand,
                phase,
                attributes,
            ]
            f.write("\t".join(fields) + "\n")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate GFF3 file of primer binding coordinates on reference"
    )
    parser.add_argument(
        "--reference-fasta",
        required=True,
        help="Reference FASTA file",
    )
    parser.add_argument(
        "--primers-fasta",
        default="data/Primers/PCR_primers/primers.fa",
        help="Primers FASTA file (default: data/Primers/PCR_primers/primers.fa)",
    )
    parser.add_argument(
        "--output-gff3",
        required=True,
        help="Output GFF3 file",
    )
    parser.add_argument(
        "--min-match-length",
        type=int,
        default=8,
        help="Minimum match length in nucleotides (default: 8)",
    )
    return parser.parse_args()


def main():
    """Main entry point."""
    args = parse_args()

    # Load reference
    ref_recs = read_fasta(args.reference_fasta)
    if not ref_recs:
        raise ValueError("No sequences found in reference FASTA")
    ref_name, ref_seq = ref_recs[0]
    print(f"Loaded reference: {ref_name} ({len(ref_seq)} bp)")

    # Load primers
    primers = load_primers(args.primers_fasta)
    print(f"Loaded {len(primers)} primers")

    # Find all matches
    all_matches = []
    for primer in primers:
        matches = find_primer_matches(ref_seq, primer, min_match_len=args.min_match_length)
        all_matches.extend(matches)
        if matches:
            print(f"  {primer.name}: {len(matches)} match(es)")

    print(f"Total matches: {len(all_matches)}")

    # Write GFF3
    write_gff3(all_matches, args.output_gff3, ref_name)
    print(f"Wrote {args.output_gff3}")


if __name__ == "__main__":
    main()
