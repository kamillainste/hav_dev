#!/usr/bin/env python3
"""Trim primer-covered sequence ends using HAV reference coordinates.

This script is intentionally standalone and opt-in. It supports two primer sources:
1) Built-in primer FASTA (default: data/Primers/PCR_primers/primers.fa)
2) Optional custom primer file (FASTA or CSV/TSV) to add/override primers

Primer coordinates are resolved against the reference sequence if not explicitly
provided. Trimming is then performed by coordinate mapping from an alignment to
reference (MAFFT --add --keeplength), not by direct sequence-end matching.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


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
    "I": {"A", "C", "G", "T"},  # inosine-like wildcard
    "5": {"A", "C", "G", "T"},  # lab notation seen in some CODEHOP exports
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
    ref_start: Optional[int] = None
    ref_end: Optional[int] = None
    source: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Trim primer-covered ends by reference coordinates after aligning queries "
            "to reference with MAFFT."
        )
    )
    parser.add_argument("--input-fasta", required=True, help="Input multi-FASTA")
    parser.add_argument("--output-fasta", required=True, help="Output trimmed FASTA")
    parser.add_argument(
        "--reference-fasta",
        required=True,
        help="Reference FASTA used for coordinate system",
    )
    parser.add_argument(
        "--primer-names",
        required=True,
        help="Comma-separated primer names to apply (e.g. HAV_6.1_codehop,HAV_10_codehop)",
    )
    parser.add_argument(
        "--builtin-primers",
        default="data/Primers/PCR_primers/primers.fa",
        help="Built-in primer file (FASTA with metadata headers)",
    )
    parser.add_argument(
        "--custom-primers",
        default="",
        help="Optional custom primer file (FASTA or CSV/TSV). Overrides built-in by name.",
    )
    parser.add_argument(
        "--report-tsv",
        default="",
        help="Optional trimming report TSV output path",
    )
    parser.add_argument(
        "--resolved-primers-tsv",
        default="",
        help="Optional output TSV of resolved primers with coordinates",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=2,
        help="Threads for MAFFT (default: 2)",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=100,
        help="Warn when trimmed sequence is shorter than this (default: 100)",
    )
    return parser.parse_args()


def read_fasta(path: str) -> List[Tuple[str, str]]:
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


def write_fasta(path: str, records: Iterable[Tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def normalize_sequence(seq: str) -> str:
    seq = seq.upper().replace(" ", "")
    seq = seq.replace("-", "")
    seq = seq.replace(".", "")
    return seq


def parse_header(header: str) -> Tuple[str, Dict[str, str]]:
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


def detect_delimiter(path: str) -> str:
    ext = Path(path).suffix.lower()
    if ext == ".tsv":
        return "\t"
    if ext == ".csv":
        return ","
    # Best effort by first line
    with open(path, "r", encoding="utf-8") as handle:
        first = handle.readline()
    if "\t" in first:
        return "\t"
    return ","


def load_primers(path: str) -> Dict[str, Primer]:
    primers: Dict[str, Primer] = {}
    ext = Path(path).suffix.lower()

    if ext in {".fa", ".fasta", ".fna"}:
        for header, seq in read_fasta(path):
            name, meta = parse_header(header)
            primer = Primer(
                name=name,
                sequence=normalize_sequence(seq),
                direction=meta.get("direction", "").upper(),
                region=meta.get("region", ""),
                comment=meta.get("comment", ""),
                ref_start=(int(meta["ref_start"]) if "ref_start" in meta and meta["ref_start"].isdigit() else None),
                ref_end=(int(meta["ref_end"]) if "ref_end" in meta and meta["ref_end"].isdigit() else None),
                source=path,
            )
            primers[name] = primer
        return primers

    delimiter = detect_delimiter(path)
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            row_l = {k.lower().strip(): (v.strip() if v is not None else "") for k, v in row.items()}
            name = row_l.get("name", "")
            if not name:
                continue
            primer = Primer(
                name=name,
                sequence=normalize_sequence(row_l.get("sequence", "")),
                direction=row_l.get("direction", "").upper(),
                region=row_l.get("region", ""),
                comment=row_l.get("comment", ""),
                ref_start=(int(row_l["ref_start"]) if row_l.get("ref_start", "").isdigit() else None),
                ref_end=(int(row_l["ref_end"]) if row_l.get("ref_end", "").isdigit() else None),
                source=path,
            )
            primers[name] = primer
    return primers


def validate_primer(primer: Primer) -> None:
    if not primer.sequence:
        raise ValueError(f"Primer '{primer.name}' has empty sequence")
    invalid = sorted({ch for ch in primer.sequence if ch not in IUPAC})
    if invalid:
        raise ValueError(
            f"Primer '{primer.name}' contains unsupported nucleotide symbols: {''.join(invalid)}"
        )
    if primer.direction not in {"F", "R"}:
        raise ValueError(
            f"Primer '{primer.name}' has invalid direction '{primer.direction}'. Expected F or R."
        )


def reverse_complement_ambiguous(seq: str) -> str:
    return "".join(RC_MAP.get(base, "N") for base in reversed(seq))


def match_ambiguous(ref_window: str, primer_seq: str) -> bool:
    if len(ref_window) != len(primer_seq):
        return False
    for ref_base, pr_base in zip(ref_window, primer_seq):
        if ref_base not in {"A", "C", "G", "T", "N"}:
            return False
        allowed = IUPAC.get(pr_base, {"A", "C", "G", "T"})
        if ref_base != "N" and ref_base not in allowed:
            return False
    return True


def infer_coordinates(reference_seq: str, primer: Primer) -> Tuple[int, int]:
    seq = primer.sequence
    if primer.direction == "R":
        seq = reverse_complement_ambiguous(seq)

    matches: List[int] = []
    plen = len(seq)
    for i in range(0, len(reference_seq) - plen + 1):
        if match_ambiguous(reference_seq[i : i + plen], seq):
            matches.append(i)

    if len(matches) == 0:
        raise ValueError(
            f"Could not map primer '{primer.name}' to reference using its sequence"
        )
    if len(matches) > 1:
        raise ValueError(
            f"Primer '{primer.name}' maps to multiple reference locations ({len(matches)} matches)"
        )

    start0 = matches[0]
    start = start0 + 1
    end = start0 + plen
    return start, end


def resolve_primers(
    reference_seq: str,
    selected_names: List[str],
    builtin_path: str,
    custom_path: str,
) -> Dict[str, Primer]:
    all_primers: Dict[str, Primer] = {}

    if builtin_path:
        if not os.path.exists(builtin_path):
            raise FileNotFoundError(f"Built-in primer file not found: {builtin_path}")
        all_primers.update(load_primers(builtin_path))

    if custom_path:
        if not os.path.exists(custom_path):
            raise FileNotFoundError(f"Custom primer file not found: {custom_path}")
        all_primers.update(load_primers(custom_path))

    missing = [name for name in selected_names if name not in all_primers]
    if missing:
        raise ValueError(
            "Unknown primer name(s): "
            + ", ".join(missing)
            + ". Add them to data/Primers/PCR_primers/primers.fa or pass --custom-primers."
        )

    selected: Dict[str, Primer] = {name: all_primers[name] for name in selected_names}

    for primer in selected.values():
        validate_primer(primer)
        if primer.ref_start is None or primer.ref_end is None:
            start, end = infer_coordinates(reference_seq, primer)
            primer.ref_start = start
            primer.ref_end = end

        if primer.ref_start < 1 or primer.ref_end < 1 or primer.ref_start > primer.ref_end:
            raise ValueError(
                f"Invalid coordinates for primer '{primer.name}': {primer.ref_start}-{primer.ref_end}"
            )

    return selected


def mafft_align_to_reference(
    reference_name: str,
    reference_seq: str,
    query_records: List[Tuple[str, str]],
    threads: int,
) -> Dict[str, str]:
    with tempfile.TemporaryDirectory(prefix="hav_primer_trim_") as tmpdir:
        ref_path = os.path.join(tmpdir, "reference.fa")
        query_path = os.path.join(tmpdir, "queries.fa")
        aln_path = os.path.join(tmpdir, "aligned.fa")

        write_fasta(ref_path, [(reference_name, reference_seq)])
        write_fasta(query_path, query_records)

        cmd = [
            "mafft",
            "--quiet",
            "--thread",
            str(max(1, threads)),
            "--add",
            query_path,
            "--keeplength",
            ref_path,
        ]

        with open(aln_path, "w", encoding="utf-8") as out_handle:
            proc = subprocess.run(cmd, stdout=out_handle, stderr=subprocess.PIPE, text=True, check=False)

        if proc.returncode != 0:
            raise RuntimeError(
                "MAFFT failed while aligning queries to reference. "
                f"Exit code {proc.returncode}. stderr: {proc.stderr.strip()}"
            )

        aligned = read_fasta(aln_path)
        return {name: seq for name, seq in aligned}


def build_ref_column_map(ref_aligned: str) -> List[Optional[int]]:
    mapping: List[Optional[int]] = []
    pos = 0
    for ch in ref_aligned:
        if ch != "-":
            pos += 1
            mapping.append(pos)
        else:
            mapping.append(None)
    return mapping


def find_col_for_left_cut(ref_map: List[Optional[int]], left_cut_pos: int) -> Optional[int]:
    col = None
    for idx, ref_pos in enumerate(ref_map):
        if ref_pos is not None and ref_pos <= left_cut_pos:
            col = idx
    return col


def find_col_for_right_cut(ref_map: List[Optional[int]], right_cut_pos: int) -> Optional[int]:
    for idx, ref_pos in enumerate(ref_map):
        if ref_pos is not None and ref_pos >= right_cut_pos:
            return idx
    return None


def count_non_gap_prefix(seq_aln: str, col_inclusive: int) -> int:
    return sum(1 for ch in seq_aln[: col_inclusive + 1] if ch != "-")


def count_non_gap_suffix(seq_aln: str, col_inclusive: int) -> int:
    return sum(1 for ch in seq_aln[col_inclusive:] if ch != "-")


def write_resolved_primers(path: str, primers: Dict[str, Primer]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "name",
                "sequence",
                "direction",
                "region",
                "comment",
                "ref_start",
                "ref_end",
                "source",
            ]
        )
        for name in sorted(primers):
            p = primers[name]
            writer.writerow(
                [
                    p.name,
                    p.sequence,
                    p.direction,
                    p.region,
                    p.comment,
                    p.ref_start,
                    p.ref_end,
                    p.source,
                ]
            )


def write_report(path: str, rows: List[Dict[str, str]]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fields = [
        "seq_name",
        "input_length",
        "output_length",
        "trimmed_left",
        "trimmed_right",
        "left_ref_cut",
        "right_ref_cut",
        "left_primers",
        "right_primers",
        "warnings",
    ]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    args = parse_args()

    reference_records = read_fasta(args.reference_fasta)
    if len(reference_records) == 0:
        raise ValueError(f"No reference sequence found in {args.reference_fasta}")
    reference_name, reference_seq = reference_records[0]
    reference_seq = normalize_sequence(reference_seq).replace("U", "T")

    query_records = read_fasta(args.input_fasta)
    if len(query_records) == 0:
        raise ValueError(f"No query sequences found in {args.input_fasta}")

    selected_names = [x.strip() for x in args.primer_names.split(",") if x.strip()]
    if not selected_names:
        raise ValueError("No primer names were provided")

    selected = resolve_primers(
        reference_seq=reference_seq,
        selected_names=selected_names,
        builtin_path=args.builtin_primers,
        custom_path=args.custom_primers,
    )

    left_primers = [p for p in selected.values() if p.direction == "F"]
    right_primers = [p for p in selected.values() if p.direction == "R"]

    if not left_primers and not right_primers:
        raise ValueError("Selected primers contain neither forward nor reverse primers")

    left_cut = max((p.ref_end for p in left_primers), default=None)
    right_cut = min((p.ref_start for p in right_primers), default=None)

    aligned = mafft_align_to_reference(
        reference_name=reference_name,
        reference_seq=reference_seq,
        query_records=[(n, normalize_sequence(s).replace("U", "T")) for n, s in query_records],
        threads=args.threads,
    )

    if reference_name not in aligned:
        raise RuntimeError("Reference sequence was not found in MAFFT output")

    ref_aln = aligned[reference_name]
    ref_map = build_ref_column_map(ref_aln)

    left_col = find_col_for_left_cut(ref_map, left_cut) if left_cut is not None else None
    right_col = find_col_for_right_cut(ref_map, right_cut) if right_cut is not None else None

    trimmed_records: List[Tuple[str, str]] = []
    report_rows: List[Dict[str, str]] = []

    for q_name, q_raw_seq in query_records:
        q_seq = normalize_sequence(q_raw_seq).replace("U", "T")
        q_aln = aligned.get(q_name)
        if q_aln is None:
            raise RuntimeError(f"Query '{q_name}' missing from MAFFT output")

        warnings: List[str] = []

        trimmed_left = 0
        trimmed_right = 0

        if left_col is not None:
            trimmed_left = count_non_gap_prefix(q_aln, left_col)
        elif left_cut is not None:
            warnings.append("left_cut_not_mappable")

        if right_col is not None:
            trimmed_right = count_non_gap_suffix(q_aln, right_col)
        elif right_cut is not None:
            warnings.append("right_cut_not_mappable")

        if trimmed_left + trimmed_right > len(q_seq):
            warnings.append("trim_overlap")
            kept = ""
            trimmed_left = len(q_seq)
            trimmed_right = 0
        else:
            kept = q_seq[trimmed_left : len(q_seq) - trimmed_right if trimmed_right > 0 else len(q_seq)]

        if len(kept) < args.min_length:
            warnings.append("short_after_trim")

        trimmed_records.append((q_name, kept))

        report_rows.append(
            {
                "seq_name": q_name,
                "input_length": str(len(q_seq)),
                "output_length": str(len(kept)),
                "trimmed_left": str(trimmed_left),
                "trimmed_right": str(trimmed_right),
                "left_ref_cut": str(left_cut) if left_cut is not None else "",
                "right_ref_cut": str(right_cut) if right_cut is not None else "",
                "left_primers": ",".join(sorted(p.name for p in left_primers)),
                "right_primers": ",".join(sorted(p.name for p in right_primers)),
                "warnings": ",".join(warnings),
            }
        )

    os.makedirs(os.path.dirname(os.path.abspath(args.output_fasta)), exist_ok=True)
    write_fasta(args.output_fasta, trimmed_records)

    if args.report_tsv:
        write_report(args.report_tsv, report_rows)

    if args.resolved_primers_tsv:
        write_resolved_primers(args.resolved_primers_tsv, selected)

    print(f"Trimmed {len(trimmed_records)} sequence(s): {args.output_fasta}")
    if args.report_tsv:
        print(f"Report written: {args.report_tsv}")
    if args.resolved_primers_tsv:
        print(f"Resolved primer coordinates: {args.resolved_primers_tsv}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
