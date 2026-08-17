#!/usr/bin/env python3
"""Split designed p35 sequences by whether the rebuilt helix contains R or K.

The helix is defined as residues 147 through the residue immediately before the
DPYRVK motif.

Example:
  python filter_no_rk.py \\
      --fasta ../02_proteinmpnn/seqs_renamed \\
      --pass-fasta scores/helix_no_RK.fa \\
      --fail-fasta scores/helix_with_RK.fa
"""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            sequences[name] = "".join(chunks)
    return sequences


def helix_segment(sequence: str, start_pos: int = 147, min_length: int = 14, max_length: int = 25) -> str | None:
    motif = "DPYRVK"
    start_idx = start_pos - 1
    for end_idx in range(start_idx + min_length - 1, min(start_idx + max_length, len(sequence))):
        if sequence[end_idx + 1 : end_idx + 7] == motif:
            return sequence[start_idx : end_idx + 1]
    return None


def write_fasta(path: Path, records: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n{sequence}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, help="FASTA file or directory")
    parser.add_argument("--pass-fasta", default="helix_no_RK.fa")
    parser.add_argument("--fail-fasta", default="helix_with_RK.fa")
    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    files = sorted(fasta_path.glob("*.fa")) if fasta_path.is_dir() else [fasta_path]

    passed = {}
    failed = {}
    skipped = 0
    for fasta in files:
        for header, sequence in parse_fasta(fasta).items():
            if "sample=" not in header.lower() and "GGGGGGGG" in sequence[146:171]:
                skipped += 1
                continue
            segment = helix_segment(sequence)
            if segment is None:
                skipped += 1
                continue
            record_name = header if header.startswith("mIL12") else f"{fasta.stem} {header}"
            if "R" not in segment and "K" not in segment:
                passed[record_name] = sequence
                print(f"[PASS] {record_name}  helix={segment}")
            else:
                failed[record_name] = sequence
                print(f"[FAIL] {record_name}  helix={segment}")

    write_fasta(Path(args.pass_fasta), passed)
    write_fasta(Path(args.fail_fasta), failed)
    total = len(passed) + len(failed)
    print(f"No R/K: {len(passed)} / {total}; skipped {skipped}")


if __name__ == "__main__":
    main()
