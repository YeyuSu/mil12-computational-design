#!/usr/bin/env python3
"""Merge ProteinMPNN loop-redesign FASTA files (T=0.2 and T=0.3).

Skips the first (wild-type / backbone) record in each file and renames samples
to mIL12_redesignloop_0_2_sampleN / mIL12_redesignloop_0_3_sampleN.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def collect_samples(path: Path, prefix: str) -> list[tuple[str, str]]:
    records = []
    header = None
    sequence = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence)))
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
        if header is not None:
            records.append((header, "".join(sequence)))

    designed = []
    for header, seq in records:
        match = re.search(r"sample=(\d+)", header)
        if not match:
            continue
        designed.append((f"{prefix}_sample{match.group(1)}", seq))
    return designed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t02", default="seqs/mIL12_B147-160_AF3predicted_mIL12_T=0.2.fa")
    parser.add_argument("--t03", default="seqs/mIL12_B147-160_AF3predicted_mIL12_T=0.3.fa")
    parser.add_argument("--output", default="seqs/mIL12_redesignloop_combined.fa")
    args = parser.parse_args()

    combined = collect_samples(Path(args.t02), "mIL12_redesignloop_0_2")
    combined += collect_samples(Path(args.t03), "mIL12_redesignloop_0_3")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        for name, sequence in combined:
            handle.write(f">{name}\n{sequence}\n")
    print(f"Wrote {len(combined)} sequences to {output}")


if __name__ == "__main__":
    main()
