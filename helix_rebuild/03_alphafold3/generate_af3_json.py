#!/usr/bin/env python3
"""Build AlphaFold 3 sequence-only JSON inputs from designed p35 FASTA files.

The template must contain protein chains A (p40, unchanged) and B (p35, replaced).
Helix designs in this repository were predicted without MSAs.

Example:
  python generate_af3_json.py \\
      --fasta ../02_proteinmpnn/seqs_renamed \\
      --template ../../../inputs/af3_sequence_template.json \\
      --outdir ./json_inputs
"""

from __future__ import annotations

import argparse
import json
import re
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


def is_backbone_record(header: str, sequence: str) -> bool:
    if "sample=" in header.lower():
        return False
    designed = sequence[146:171]
    return "G" * 10 in designed


def json_name(header: str, fasta_stem: str) -> str:
    sample = re.search(r"sample=(\d+)", header)
    if sample and fasta_stem.startswith("mIL12_alphahelix_14-25aa_"):
        return f"{fasta_stem}_{sample.group(1)}"
    if header.startswith("mIL12_") or header.startswith("mil12_"):
        return re.sub(r"[^A-Za-z0-9._-]+", "_", header)
    if sample:
        return f"{fasta_stem}_{sample.group(1)}"
    return re.sub(r"[^A-Za-z0-9._-]+", "_", header)[:80]


def replace_chain_b(template: dict, name: str, sequence: str) -> dict:
    data = json.loads(json.dumps(template))
    data["name"] = name
    replaced = False
    for entry in data.get("sequences", []):
        protein = entry.get("protein", {})
        if protein.get("id") != "B":
            continue
        protein["sequence"] = sequence
        protein.pop("unpairedMsa", None)
        protein.pop("pairedMsa", None)
        replaced = True
        break
    if not replaced:
        raise ValueError("Template has no protein chain with id 'B'")
    return data


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, help="Designed p35 FASTA file or directory")
    parser.add_argument(
        "--template",
        default="../../../inputs/af3_sequence_template.json",
        help="AF3 JSON template with chains A (p40) and B (p35)",
    )
    parser.add_argument("--outdir", default="./json_inputs", help="Output directory")
    parser.add_argument(
        "--keep-backbone",
        action="store_true",
        help="Also emit JSON for poly-G backbone records",
    )
    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    fasta_files = sorted(fasta_path.glob("*.fa")) if fasta_path.is_dir() else [fasta_path]
    template = json.loads(Path(args.template).read_text(encoding="utf-8"))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    n_written = 0
    for fasta in fasta_files:
        for header, sequence in parse_fasta(fasta).items():
            if not args.keep_backbone and is_backbone_record(header, sequence):
                continue
            name = json_name(header, fasta.stem)
            payload = replace_chain_b(template, name, sequence)
            (outdir / f"{name}.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
            n_written += 1

    print(f"Generated {n_written} JSON files in {outdir}")


if __name__ == "__main__":
    main()
