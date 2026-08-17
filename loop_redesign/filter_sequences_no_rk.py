#!/usr/bin/env python3
"""Split loop-redesign sequences by whether p35 residues 147-160 contain R or K.

Positions are 1-based inclusive: 147-160 is the wild-type loop that ProteinMPNN
redesigned. The DPYRVK motif starts at residue 161 and is not part of the filter.
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


def write_fasta(path: Path, records: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n{sequence}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", default="seqs/mIL12_redesignloop_combined.fa")
    parser.add_argument("--pass-fasta", default="seqs/mIL12_redesignloop_no_RK_147-160.fa")
    parser.add_argument("--fail-fasta", default="seqs/mIL12_redesignloop_with_RK_147-160.fa")
    parser.add_argument("--start", type=int, default=147, help="1-based start residue")
    parser.add_argument("--end", type=int, default=160, help="1-based inclusive end residue")
    args = parser.parse_args()

    start_idx = args.start - 1
    sequences = parse_fasta(Path(args.fasta))
    passed = {}
    failed = {}
    for name, sequence in sequences.items():
        if len(sequence) < args.end:
            print(f"Skip short sequence {name} (len={len(sequence)})")
            continue
        fragment = sequence[start_idx:args.end]
        if "R" not in fragment and "K" not in fragment:
            passed[name] = sequence
            print(f"[PASS] {name}  147-160={fragment}")
        else:
            failed[name] = sequence
            print(f"[FAIL] {name}  147-160={fragment}")

    write_fasta(Path(args.pass_fasta), passed)
    write_fasta(Path(args.fail_fasta), failed)
    total = len(passed) + len(failed)
    print(f"No R/K: {len(passed)} / {total} ({100 * len(passed) / total:.1f}%)" if total else "No sequences")


if __name__ == "__main__":
    main()
